"""Write marker posteriors and summarize sibling sharing across chromosomes."""

from contextlib import contextmanager
import gzip
from pathlib import Path
from tempfile import TemporaryDirectory

import h5py
import numpy as np

from snipar.utilities import encode_str_array, outfile_name


def _interval_weights(positions, genetic_map):
    """Trapezoidal weights on the span of the retained markers."""
    steps = np.stack((np.diff(genetic_map), np.diff(positions) / 1e6)) * .5
    weights = np.zeros((2, len(positions)))
    weights[:, :-1] += steps
    weights[:, 1:] += steps
    return weights


class Relatedness:
    """Accumulate cM- and Mb-weighted state lengths in original pair order."""

    def __init__(self, pairs):
        self.pairs = pairs
        self.lengths = np.zeros((len(pairs), 2, 3))
        self.spans = np.zeros((len(pairs), 2))
        self.chromosomes = np.zeros(len(pairs), dtype=np.int32)

    @contextmanager
    def chromosome_batch(self, rows, posterior_prefix=None):
        """Return an inference callback and atomically finish a posterior file."""
        weights = span = output = scratch = target = None
        written = 0

        def store(chrom, pairs, snps, positions, genetic_map, probabilities, pair_indices):
            nonlocal weights, span, output, scratch, target, written
            if weights is None:
                weights = _interval_weights(positions, genetic_map)
                span = weights.sum(axis=1)
                if posterior_prefix is not None:
                    target = Path(outfile_name(posterior_prefix, '.ibd.posterior.hdf5', chrom))
                    scratch = TemporaryDirectory(prefix='.snipar-ibd-', dir=target.parent)
                    output = h5py.File(Path(scratch.name) / target.name, 'w')
                    output.attrs.update(format_version=1, chromosome=str(chrom))
                    output['positions'] = np.asarray(positions, dtype=np.int64)
                    output['variant_ids'] = encode_str_array(snps)
                    output['genetic_map_cm'] = genetic_map
                    pair_dtype = encode_str_array(self.pairs[rows]).dtype
                    output.create_dataset('pairs', (0, 2), maxshape=(len(rows), 2), dtype=pair_dtype)
                    output.create_dataset('probabilities', (0, len(positions), 3),
                                          maxshape=(len(rows), len(positions), 3), dtype=np.float32,
                                          chunks=(1, min(len(positions), 1024), 3),
                                          compression='lzf', shuffle=True)
            selected = rows[pair_indices]
            self.lengths[selected] += np.einsum('pmk,dm->pdk', probabilities, weights, dtype=np.float64)
            self.spans[selected] += span
            self.chromosomes[selected] += 1
            if output is not None:
                stop = written + len(pairs)
                output['pairs'].resize(stop, axis=0)
                output['pairs'][written:stop] = encode_str_array(pairs)
                output['probabilities'].resize(stop, axis=0)
                output['probabilities'][written:stop] = probabilities
                written = stop

        try:
            yield store
            if output is not None:
                output.close()
                (Path(scratch.name) / target.name).replace(target)
        finally:
            if output is not None:
                output.close()
            if scratch is not None:
                scratch.cleanup()

    def write(self, prefix):
        target = Path(str(prefix).replace('@', 'genome') + '.ibd.relatedness.tsv.gz')
        with TemporaryDirectory(prefix='.snipar-ibd-', dir=target.parent) as directory:
            pending = Path(directory) / target.name
            with gzip.open(pending, 'wt', compresslevel=3, encoding='utf-8') as output:
                columns = ('ID1', 'ID2', 'N_CHR', 'Span_cM', 'IBD0_cM', 'IBD1_cM', 'IBD2_cM',
                           'PropIBD_cM', 'Span_Mb', 'IBD0_Mb', 'IBD1_Mb', 'IBD2_Mb', 'PropIBD_Mb')
                output.write('\t'.join(columns) + '\n')
                for row in np.flatnonzero(self.chromosomes):
                    values = []
                    for lengths, span in zip(self.lengths[row], self.spans[row]):
                        share = (.5 * lengths[1] + lengths[2]) / span if span > 0 else np.nan
                        values.extend(f'{value:.10g}' for value in (span, *lengths, share))
                    output.write('\t'.join((*self.pairs[row], str(self.chromosomes[row]), *values)) + '\n')
            pending.replace(target)
        return target
