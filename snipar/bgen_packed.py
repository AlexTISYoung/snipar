"""Gather selected samples from common diploid, biallelic BGEN blocks."""

from functools import lru_cache
import platform
import struct
import sys
import zlib

import numpy as np

_ZSTD_ERROR = zlib.error
_zstd_decompress = None
if sys.version_info >= (3, 14):
    try:
        from compression import zstd
    except ImportError:
        zstd = None

    if zstd is not None:
        _ZSTD_ERROR = zstd.ZstdError

        def _zstd_decompress(data, size):
            decoder = zstd.ZstdDecompressor()
            result = decoder.decompress(data, max_length=size + 1)
            return result if decoder.eof else b''
else:
    try:
        import zstandard
    except ImportError:
        zstandard = None

    if zstandard is not None:
        _ZSTD_ERROR = zstandard.ZstdError

        def _zstd_decompress(data, size):
            # A known frame content size takes precedence over max_output_size.
            # Check it first to bound the output allocation in both cases.
            content_size = zstandard.frame_content_size(data)
            if content_size not in (-1, size):
                return b''
            return zstandard.ZstdDecompressor().decompress(data, max_output_size=size)


_BLOCK_HEADER = struct.Struct('<IHBB')
_UINT32 = struct.Struct('<I')
_NATIVE_PHASED_8 = platform.machine().lower() not in ('amd64', 'x86_64')


def subset_samples(samples, count, selected_count, *, minimum_count=256):
    """Choose the packed path only when it can avoid substantial decoding work."""
    if count < minimum_count or not 0 < selected_count <= count // 4:
        return None
    if isinstance(samples, slice):
        return np.arange(*samples.indices(count), dtype=np.intp)
    if samples.ndim != 1:
        return None
    if samples.dtype.kind == 'b':
        return np.flatnonzero(samples) if len(samples) == count else None
    if samples.dtype.kind not in 'iu' or np.any(samples >= count) or np.any(samples < -count):
        return None
    if np.any(samples < 0):
        samples = np.where(samples < 0, samples + count, samples)
    return samples.astype(np.intp, copy=False)


@lru_cache(maxsize=None)
def _probability_table(bits, native):
    probabilities = np.arange(1 << bits, dtype=np.float64) / ((1 << bits) - 1)
    if native and bits == 8:
        # bgen's scalar and unphased SIMD code use its seven-decimal LUT.
        np.round(probabilities, 7, out=probabilities)
    probabilities = probabilities.astype(np.float32)
    probabilities.flags.writeable = False
    return probabilities


class PackedBgenSubset:
    """Read fixed-width states for the requested diploid samples from each block."""

    @classmethod
    def from_stream(cls, stream, count, samples, *, native=False, native_table=False):
        # The layout and compression flags are the final word of the file header.
        stream.seek(4)
        length_bytes = stream.read(4)
        if len(length_bytes) != 4:
            return None
        header_size, = _UINT32.unpack(length_bytes)
        if header_size < 20:
            return None
        stream.seek(header_size)
        flag_bytes = stream.read(4)
        if len(flag_bytes) != 4:
            return None
        flags, = _UINT32.unpack(flag_bytes)
        compression = flags & 3
        if (flags >> 2) & 15 != 2 or compression == 3:
            return None
        if compression == 2 and _zstd_decompress is None:
            return None
        return cls(stream, count, samples, compression, native, native_table)

    def __init__(self, stream, count, samples, compression, native, native_table):
        self.stream = stream
        self.count = count
        self.samples = samples
        self.compression = compression
        self.native = native
        self.native_table = native_table
        self.active = True
        self.minimum_length = 10 + count + (count + 3) // 4
        self.maximum_length = 10 + 5 * count

    def _block(self, offset):
        self.stream.seek(offset)
        header_size = 8 if self.compression else 4
        header = self.stream.read(header_size)
        if len(header) != header_size:
            return None
        compressed_size, = _UINT32.unpack_from(header)
        if self.compression:
            length, = _UINT32.unpack_from(header, 4)
            size = compressed_size - 4
        else:
            length = size = compressed_size
        if not self.minimum_length <= length <= self.maximum_length or not 0 < size <= 2 * length + 1024:
            return None
        data = self.stream.read(size)
        if len(data) != size:
            return None
        if self.compression:
            try:
                if self.compression == 1:
                    decoder = zlib.decompressobj()
                    data = decoder.decompress(data, length + 1)
                    if not decoder.eof:
                        return None
                else:
                    data = _zstd_decompress(data, length)
            except (zlib.error, _ZSTD_ERROR):
                return None
        return data if len(data) == length else None

    def read(self, offset):
        """Return the selected probabilities and phase, or defer to the general reader."""
        if not self.active:
            return None
        result = self._read(offset)
        # Files commonly use one layout throughout. Keep an unsupported file
        # on the regular reader after one probe to avoid repeated inflation.
        self.active = result is not None
        return result

    def _read(self, offset):
        data = self._block(offset)
        if data is None or _BLOCK_HEADER.unpack_from(data) != (self.count, 2, 2, 2):
            return None
        phased, bits = data[self.count + 8:self.count + 10]
        if phased not in (0, 1) or bits not in (1, 2, 4, 8, 16):
            return None
        if self.native and bits != 1 and not (bits == 8 and self.native_table and
                                              (not phased or _NATIVE_PHASED_8)):
            return None
        if len(data) != self.count + 10 + (2 * self.count * bits + 7) // 8:
            return None

        ploidy = np.frombuffer(data, np.uint8, self.count, 8)[self.samples]
        if np.any((ploidy & 127) != 2):
            return None
        start = self.count + 10
        maximum = (1 << bits) - 1
        if bits < 8:
            sample_bits = 2 * bits * self.samples
            packed = np.frombuffer(data, np.uint8, offset=start)[sample_bits // 8]
            shifts = sample_bits % 8
            first = (packed >> shifts) & maximum
            second = (packed >> (shifts + bits)) & maximum
        else:
            dtype = np.uint8 if bits == 8 else '<u2'
            states = np.frombuffer(data, dtype, offset=start).reshape(self.count, 2)[self.samples]
            first, second = states.T
        table = _probability_table(bits, self.native)
        result = np.empty((len(self.samples), 4 if phased else 3), dtype=np.float32)
        result[:, 0] = table[first]
        if phased:
            result[:, 1] = table[maximum - first]
            result[:, 2] = table[second]
            result[:, 3] = table[maximum - second]
        else:
            result[:, 1] = table[second]
            remaining = maximum - first.astype(np.int32) - second
            if np.any(remaining < 0):
                return None
            result[:, 2] = table[remaining]
        result[ploidy == 130] = np.nan
        return result, bool(phased)
