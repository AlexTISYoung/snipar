import errno
import hashlib
from math import isqrt
import os
from tempfile import TemporaryDirectory
from typing import Any, Dict, Optional, Union
from pathlib import Path
import numpy as np
import pandas as pd
from os import path
from typing import Tuple

from snipar.docgen import get_parser_doc
from snipar.bgen_filter import BGEN_SUBSET_SUFFIX, open_bgen_subset

def make_id_dict(x,col=0):
    """
    Make a dictionary that maps from the values in the given column (col) to their row-index in the input array
    """
    if len(x.shape)>1:
        x = x[:,col]
    return {identifier: i for i, identifier in enumerate(x)}

def convert_str_array(x):
    """
    Decode an array of UTF-8 strings.
    """
    x = np.asarray(x)
    if x.dtype.kind == 'U' or not x.size:
        return x.astype(str, copy=False)
    if x.dtype.kind == 'O':
        return np.array([value.decode('utf-8') if isinstance(value, bytes) else str(value)
                         for value in x.flat]).reshape(x.shape)
    return np.char.decode(np.asarray(x, dtype='S'), 'utf-8')

def encode_str_array(x):
    """
    Encode an array of strings as UTF-8 bytes.
    """
    x = np.asarray(x)
    if x.dtype.kind == 'S' or not x.size:
        return x.astype('S', copy=False)
    return np.char.encode(x.astype(str, copy=False), 'utf-8')

def parse_obsfiles(obsfiles, obsformat='bed', append = True, wildcard = '@', chromosomes=None):
    obs_files = []
    chroms = []
    if wildcard in obsfiles:
        if chromosomes is None:
            chromosomes = np.arange(1,23)
        for i in chromosomes:
            obsfile = obsfiles.replace(wildcard, str(i))+'.'+obsformat
            if path.exists(obsfile):
                if append:
                    obs_files.append(obsfile)
                else:
                    obs_files.append(obsfile[:-(len(obsformat)+1)])
                chroms.append(i)
        if len(obs_files)==0:
            raise(ValueError('Observed genotype files not found'))
        else:
            print(str(len(obs_files))+' files found')
    else:
            obsfile = obsfiles+'.'+obsformat
            if path.exists(obsfile):
                obs_files = [obsfile]
                chroms = [0]
            else:
                raise(ValueError(obsfile+' not found'))
    return np.array(obs_files), np.array(chroms,dtype=int)

def parse_filelist(obsfiles, impfiles, obsformat, chromosomes=None):
    if impfiles is None:
        raise ValueError('Must provide a prefix for the imputed genotype files')
    obs_files = []
    imp_files = []
    chroms = []
    if '@' in obsfiles or '@' in impfiles:
        if chromosomes is None:
            chromosomes = np.arange(1,23)
        for i in chromosomes:
            obsfile = obsfiles.replace('@', str(i))+'.'+obsformat
            impfile = impfiles.replace('@', str(i))+'.hdf5'
            if path.exists(impfile) and path.exists(obsfile):
                obs_files.append(obsfile)
                imp_files.append(impfile)
                chroms.append(i)
        if len(imp_files)==0:
            raise(ValueError('Observed/imputed genotype files not found'))
        else:
            print(str(len(imp_files))+' matched observed and imputed genotype files found')
    else:
            obsfile = obsfiles+'.'+obsformat
            impfile = impfiles+'.hdf5'
            if path.exists(obsfile) and path.exists(impfile):
                obs_files = [obsfile]
                imp_files = [impfile]
                chroms = [0]
            else:
                if not path.exists(obsfile):
                    raise(ValueError(obsfile+' not found'))
                if not path.exists(impfile):
                    raise(ValueError(impfile+' not found'))
    return np.array(obs_files), np.array(imp_files), np.array(chroms,dtype=int)

def outfile_name(outprefix,outsuffix,chrom=None):
    if '@' in outprefix:
        if chrom is None:
            raise(ValueError('Must provide chromosome number with wildcard character'))
        outprefix = outprefix.replace('@', str(chrom))
        return outprefix+outsuffix
    elif chrom is not None:
        return outprefix+'_chr_'+str(chrom)+outsuffix
    else:
        return outprefix+outsuffix

def coord2linear(ind1: int, ind2: int) -> int:
    row_ind, col_ind = sorted((int(ind1), int(ind2)), reverse=True)
    if col_ind < 0:
        raise ValueError('Matrix coordinates must be nonnegative')
    return row_ind * (row_ind + 1) // 2 + col_ind


def linear2coord(ind: int) -> Tuple[int, int]:
    ind = int(ind)
    if ind < 0:
        raise ValueError('Packed matrix index must be nonnegative')
    row = (isqrt(8 * ind + 1) - 1) // 2
    return row, ind - row * (row + 1) // 2

def read_oxford_sample(filename):
    """Read individual IDs, preserving leading zeros and skipping the column-type row."""
    samples = pd.read_csv(filename, sep=r"\s+", dtype=str, skiprows=[1], keep_default_na=False)
    return samples["ID_2"].to_numpy()


def _legacy_bgen(filename, *, lazy=False, allow_complex=False, verbose=False):
    """Use bgen-reader's alternate index path when the adjacent index cannot be written."""
    import bgen_reader
    from bgen_reader._environment import BGEN_READER_CACHE_HOME

    source = Path(filename).resolve()
    suffix = '.complex.metadata2.mmm' if allow_complex else '.metadata2.mmm'
    if lazy:
        suffix = '.metafile'

    def read(index):
        if lazy:
            return bgen_reader.read_bgen(source, metafile_filepath=index, verbose=verbose)
        return bgen_reader.open_bgen(source, metadata_filepath=index, allow_complex=allow_complex, verbose=verbose)

    stat = source.stat()

    def current(index):
        try:
            return index.stat().st_mtime_ns >= stat.st_mtime_ns
        except FileNotFoundError:
            return False

    def metadata_error(error):
        # A current timestamp does not guarantee that the metadata mmap is
        # readable. Rebuild it only when the failure came from its parser.
        traceback = error.__traceback__
        while traceback is not None:
            module = traceback.tb_frame.f_globals.get('__name__')
            if module == 'bgen_reader._multimemmap':
                return True
            traceback = traceback.tb_next
        return False

    def open_at(index):
        if current(index):
            try:
                return read(index)
            except FileNotFoundError:
                pass
            except Exception as error:
                if lazy or not metadata_error(error):
                    raise
        # bgen-reader uses a fixed temporary filename. Publish a completed
        # index from a private directory when multiple workers build at once.
        with TemporaryDirectory(dir=index.parent) as directory:
            temporary = Path(directory) / 'index'
            if lazy:
                bgen_reader.create_metafile(source, temporary, verbose=verbose)
            else:
                read(temporary).close()
            if temporary.stat().st_mtime_ns < stat.st_mtime_ns:
                os.utime(temporary, ns=(stat.st_mtime_ns, stat.st_mtime_ns))
            try:
                temporary.replace(index)
            except PermissionError as error:
                # Windows can keep an index published by another worker mapped.
                if not current(index):
                    raise
                try:
                    return read(index)
                except Exception:
                    raise error
        if lazy:
            # Mirror bgen-reader's invalidation after replacing a .metafile.
            from bgen_reader import _genotype, _metafile
            for module in (_genotype, _metafile):
                with module.lock:
                    module.cache.clear()
        return read(index)

    try:
        return open_at(source.with_name(source.name + suffix))
    except OSError as error:
        if error.errno not in (errno.EACCES, errno.EPERM, errno.EROFS):
            raise

    # Version the key so an updated input never replaces another reader's
    # mapped index. Full paths distinguish identically named chromosomes.
    identity = os.fsencode(source) + f'\0{stat.st_size}:{stat.st_mtime_ns}'.encode()
    key = hashlib.sha256(identity).hexdigest()
    return open_at(BGEN_READER_CACHE_HOME.resolve() / 'metafile' / (key + suffix))


def _checked_bgen_sample(filename, nsamples):
    samples = read_oxford_sample(filename)
    if len(samples) != nsamples:
        raise ValueError('sample file length and bgen file length do not match.')
    return samples


class _bgen:
    """Simple wrapper class of bgen_reader.open_bgen that controls access to sample ids."""
    def __init__(self, 
                 filepath: Union[str, Path],
                 samples_filepath: Optional[Union[str, Path]] = None,
                 allow_complex: bool = False,
                 verbose: bool = True, reader=None):
        self._bgen = reader if reader is not None else _legacy_bgen(filepath, verbose=verbose, allow_complex=allow_complex)
        try:
            self._samples = _checked_bgen_sample(samples_filepath, self._bgen.nsamples)
        except BaseException:
            self._bgen.close()
            raise
    
    @property
    def samples(self) -> np.ndarray:
        return self._samples
    
    def __getattr__(self, name: str) -> Any:
        """Access to attributes other than sample ids."""
        return getattr(self._bgen, name)
    
    def __enter__(self):
        return self

    def __exit__(self, *args):
        self._bgen.close()


def open_bgen(filename: str, verbose: bool = False):
    """Open BGEN genotypes or a temporary SNIPar variant selection."""
    if os.fsdecode(filename).endswith(BGEN_SUBSET_SUFFIX):
        return open_bgen_subset(filename, _open_bgen, verbose=verbose)
    return _open_bgen(filename, verbose=verbose)


def _open_bgen(filename, verbose=False):
    try:
        import bgen_reader
    except ModuleNotFoundError as exc:
        if exc.name != "bgen_reader":
            raise
        from snipar.bgen_io import BgenFile
        return BgenFile(filename)
    bgen = _legacy_bgen(filename, verbose=verbose)
    if len(bgen.samples) and bgen.samples[0] == 'sample_0':
        print('WARNING: Sample ids in bgen file are generic. Trying to read the corresponding .sample file ...')
        samples_filepath = Path(filename).with_suffix('.sample')
        bgen = _bgen(filename, samples_filepath=samples_filepath, reader=bgen)
    return bgen


def read_bgen(filename: str, verbose: bool = False) -> Dict:
    """Wrapper of bgen_reader.read_bgen that checks if sample ids make sense."""
    try:
        import bgen_reader
    except ModuleNotFoundError as exc:
        if exc.name != "bgen_reader":
            raise
        raise ImportError("The lazy read_bgen API requires snipar[legacy-bgen]. "
                          "Use open_bgen for array access.") from exc
    bgen = _legacy_bgen(filename, lazy=True, verbose=verbose)
    if len(bgen['samples']) and bgen['samples'][0] == 'sample_0':
        print('WARNING: Sample ids in bgen file are generic. Trying to read the corresponding .sample file ...')
        samples_filepath = Path(filename).with_suffix('.sample')
        if not path.exists(samples_filepath):
            raise FileNotFoundError(f'{samples_filepath} does not exist.')
        bgen['samples'] = _checked_bgen_sample(samples_filepath, len(bgen['samples']))
    return bgen
