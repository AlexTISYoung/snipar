#!/usr/bin/env python
"""Write proband, parental, and cross-generation relatedness matrices."""

import argparse
from contextlib import ExitStack

import h5py
import numpy as np
from bed_reader import open_bed

from snipar.read.bed import _read_iids
from snipar.resources import available_cpus
from snipar.utilities import encode_str_array


def _observed_trios(bedfile, ped, min_maf, max_missing):
    bedfile = str(bedfile)
    if bedfile.lower().endswith('.bed'):
        bedfile = bedfile[:-4]
    bedfile += '.bed'
    pedigree_ids = ped[:, 1:4]
    requested = set(pedigree_ids.flat)
    with open_bed(bedfile, count_A1=True) as bed:
        index = {iid: row for row, iid in enumerate(_read_iids(bedfile, bed)) if iid in requested}
        rows = np.fromiter((index.get(iid, -1) for iid in pedigree_ids.flat), dtype=np.intp,
                           count=pedigree_ids.size).reshape(pedigree_ids.shape)
        keep = (rows >= 0).all(axis=1) & (ped[:, 2:4] != '0').all(axis=1)
        if not np.any(keep):
            raise ValueError('No individuals with both parents genotyped found')
        observed, mapping = np.unique(rows[keep], return_inverse=True)
        gts = bed.read(index=(observed, slice(None)), dtype='float32', order='C', num_threads=available_cpus())
    counts = np.isfinite(gts).sum(axis=0)
    freqs = np.nansum(gts, axis=0) / np.maximum(2 * counts, 1)
    variants = ((freqs >= min_maf) & (freqs <= 1 - min_maf)
                & (100 * ((len(observed) - counts) / len(observed)) <= max_missing))
    if not variants.all():
        gts = np.compress(variants, gts, axis=1)
    np.copyto(gts, 2 * freqs[variants], where=np.isnan(gts))
    proband, father, mother = mapping.reshape(-1, 3).T
    parent = gts[father]
    parent += gts[mother]
    return gts[proband], parent, ped[keep, :2]


def _with_imputation(bedfile, imputation, ids, min_maf, max_missing):
    from snipar.read import get_gts_matrix

    matrix = get_gts_matrix(bedfile=bedfile, par_gts_f=imputation, ids=ids, parsum=True)
    proband = matrix.gts[:, 0, :].filled(np.nan)
    parent = matrix.gts[:, 1, :].filled(np.nan)
    called = np.isfinite(proband)
    freqs = np.nansum(proband, axis=0) / np.maximum(2 * called.sum(axis=0), 1)
    variants = ((freqs >= min_maf) & (freqs <= 1 - min_maf)
                & (100 * (~called).mean(axis=0) <= max_missing))
    return (np.compress(variants, proband, axis=1), np.compress(variants, parent, axis=1),
            np.column_stack((matrix.fams, matrix.ids)))


def _standardize(dosages):
    called = np.isfinite(dosages)
    count = called.sum(axis=0)
    mean = np.nansum(dosages, axis=0, dtype=np.float64) / np.maximum(count, 1)
    centered = dosages.astype(np.float32, copy=False)
    centered -= mean
    np.copyto(centered, 0, where=~called)
    sd = np.sqrt(np.einsum('ij,ij->j', centered, centered, dtype=np.float64) / len(dosages))
    valid = sd > 0
    np.divide(centered, sd, out=centered, where=valid)
    return centered, valid


def _write_matrices(proband, parent, pedigree, prefix, hdf5):
    n, variants = proband.shape
    names = ('R', 'R_par', 'R_o_par')
    roots = [f'{prefix}_{name}' for name in names]
    # Use bounded work buffers even when the three output matrices exceed RAM.
    block_rows = max(1, min(n, (32 * 1024 ** 2) // (n * proband.itemsize)))
    with ExitStack() as stack:
        if hdf5:
            output = stack.enter_context(h5py.File(f'{prefix}_grms.hdf5', 'w'))
            output['N'] = variants
            output['ped'] = encode_str_array(pedigree)
            matrices = [output.create_dataset(name, (n, n), dtype='f4', compression='gzip',
                                              compression_opts=1, shuffle=True) for name in names]
        else:
            matrices = [stack.enter_context(open(f'{root}.grm.bin', 'wb')) for root in roots]
            ns = [stack.enter_context(open(f'{root}.grm.N.bin', 'wb')) for root in roots]
            for root in roots:
                np.savetxt(f'{root}.grm.id', pedigree, fmt='%s')
            np.savetxt('RDR_GRMs.txt', roots, fmt='%s')
        for start in range(0, n, block_rows):
            end = min(start + block_rows, n)
            p, m = proband[start:end], parent[start:end]
            cols = n if hdf5 else end
            for i, (left, right) in enumerate(((p, proband), (m, parent), (p, parent))):
                values = left @ right[:cols].T
                if i == 2:
                    values += m @ proband[:cols].T
                    values *= np.sqrt(0.5) / variants
                else:
                    values /= variants
                if hdf5:
                    matrices[i][start:end, :] = values
                else:
                    packed = np.concatenate([row[:j + 1] for j, row in enumerate(values, start)])
                    packed.astype('<f4', copy=False).tofile(matrices[i])
                    np.full(packed.size, variants, dtype='<f4').tofile(ns[i])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('gts', help='Path to BED file with observed genotypes')
    parser.add_argument('outprefix', help='Prefix for output files')
    parser.add_argument('--ped', help='Pedigree file (FID, IID, father, mother)')
    parser.add_argument('--pargts', help='HDF5 file from impute.py with parental genotypes and pedigree')
    parser.add_argument('--min_maf', type=float, default=.01, help='Minimum minor allele frequency (default 0.01)')
    parser.add_argument('--max_missing', type=float, default=1,
                        help='Maximum percentage of missing calls (default 1)')
    parser.add_argument('--hdf5_out', action='store_true', help='Write full matrices to HDF5')
    args = parser.parse_args(argv)
    if args.ped is None and args.pargts is None:
        parser.error('Provide --ped, --pargts, or both')
    if not 0 <= args.min_maf <= .5 or not 0 <= args.max_missing <= 100:
        parser.error('MAF must be within [0, 0.5] and missingness within [0, 100]')
    ped = np.loadtxt(args.ped, dtype=str, skiprows=1, ndmin=2) if args.ped else None
    if ped is not None and (ped.shape[0] == 0 or ped.shape[1] < 4):
        parser.error('Pedigree must contain FID, IID, father and mother columns')
    if args.pargts:
        proband, parent, pedigree = _with_imputation(
            args.gts, args.pargts, None if ped is None else ped[:, 1], args.min_maf, args.max_missing,
        )
    else:
        proband, parent, pedigree = _observed_trios(args.gts, ped, args.min_maf, args.max_missing)
    proband, valid_proband = _standardize(proband)
    parent, valid_parent = _standardize(parent)
    valid = valid_proband & valid_parent
    if not np.any(valid):
        parser.error('No variants with variation in both proband and parental dosages passed filters')
    print(f'Computing relatedness matrices for {len(pedigree)} individuals and {np.sum(valid)} SNPs')
    if not valid.all():
        proband = np.compress(valid, proband, axis=1)
        parent = np.compress(valid, parent, axis=1)
    _write_matrices(proband, parent, pedigree, args.outprefix, args.hdf5_out)


if __name__ == '__main__':
    main()
