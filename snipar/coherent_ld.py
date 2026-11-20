"""Genotype-derived LD scores for the coherent correlation command.

The existing SNIPar LD estimator is imported only when genotype input is used.
This adapter resolves files, checks chromosome identity, and aligns the returned
scores to summary-statistic rows without changing the LD calculation.
"""

from pathlib import Path

import numpy as np
import pandas as pd


def ldscores_from_bed(bedfile, chrom, ld_wind, ld_out=None):
    """Load the genotype dependencies only for the ``--bed`` input route."""
    try:
        from snipar.ld import ldscores_from_bed as compute
    except ImportError as exc:
        raise RuntimeError(
            "Computing LD scores from --bed requires SNIPar's genotype dependencies; "
            "install them or supply --ldscores instead."
        ) from exc
    return compute(bedfile, chrom, ld_wind, ld_out)


def _bed_files(pattern, chromosomes):
    pattern = str(pattern).replace('#', '@')
    requests = [(pattern.replace('@', str(chrom)), chrom) for chrom in chromosomes] if '@' in pattern else [(pattern, None)]
    files = []
    for value, requested_chrom in requests:
        bed = Path(value if value.endswith('.bed') else value + '.bed')
        if not bed.is_file():
            if '@' in pattern:
                continue
            raise ValueError('Observed genotype file not found: ' + str(bed))
        for suffix in ('.bim', '.fam'):
            companion = bed.with_suffix(suffix)
            if not companion.is_file():
                raise ValueError('Missing PLINK companion file: ' + str(companion))
        files.append((bed, requested_chrom))
    if not files:
        raise ValueError('No observed genotype files match the retained chromosomes: ' + pattern)
    return files


def _bim_metadata(bed, requested_chrom):
    frame = pd.read_csv(bed.with_suffix('.bim'), sep=r'\s+', header=None,
                        usecols=[0, 1, 3], names=['CHR', 'SNP', 'BP'],
                        dtype={'SNP': str})
    if frame.empty or frame['SNP'].isna().any() or frame['SNP'].duplicated().any():
        raise ValueError('BIM files must contain nonmissing, unique SNP IDs: ' + str(bed))
    for name in ('CHR', 'BP'):
        values = pd.to_numeric(frame[name], errors='coerce').to_numpy(dtype=float)
        if (not np.isfinite(values).all() or np.any(values != np.floor(values))
                or np.any(values < (1 if name == 'CHR' else 0))):
            raise ValueError('Invalid ' + name + ' values in ' + str(bed.with_suffix('.bim')))
        frame[name] = values.astype(int)
    chroms = frame['CHR'].unique()
    if len(chroms) != 1:
        raise ValueError('--bed requires one chromosome per file; split ' + str(bed) + ' by chromosome.')
    chrom = int(chroms[0])
    if requested_chrom is not None and chrom != requested_chrom:
        raise ValueError('BIM chromosome does not match the @ filename for ' + str(bed))
    return frame.set_index('SNP'), chrom


def compute_ldscores(bed_pattern, snp_ids, chromosomes, ld_wind=1.0, ld_out=None):
    """Return LD scores in input-row order and file provenance.

    ``bed_pattern`` accepts a PLINK prefix or a ``.bed`` filename, with ``@`` or ``#``
    expanded over the input chromosomes. Missing chromosome files and absent
    SNPs leave NaN scores for the caller's ordinary LD filtering. Each file
    must contain a single chromosome. SNPs are matched by both ID and chromosome.
    Set ``chromosomes=None`` for legacy single-file summary statistics without a
    chromosome column: SNP IDs alone are matched, and must be unique across BEDs.

    ``ld_out`` writes all computed scores in LDSC format, using
    ``PREFIX<chromosome>.l2.ldscore.gz`` or substituting a supplied ``@``/``#``.
    The returned metadata contains ``bed_files`` and ``ldscore_files`` lists.
    """
    snp_ids = np.asarray(snp_ids, dtype=str)
    if snp_ids.ndim != 1:
        raise ValueError('SNP IDs must be a one-dimensional array.')
    unknown_chromosomes = chromosomes is None
    chromosomes = np.full(len(snp_ids), np.nan) if unknown_chromosomes else np.asarray(chromosomes, dtype=float)
    if chromosomes.shape != snp_ids.shape:
        raise ValueError('SNP IDs and chromosomes must be matching one-dimensional arrays.')
    if not np.isfinite(ld_wind) or ld_wind <= 0:
        raise ValueError('ld_wind must be finite and positive.')
    valid_chrom = np.isfinite(chromosomes) & (chromosomes >= 1) & (chromosomes == np.floor(chromosomes))
    requested = list(range(1, 23)) if unknown_chromosomes else sorted(set(chromosomes[valid_chrom].astype(int)))
    if not requested:
        raise ValueError('No valid summary-statistic chromosomes for LD-score computation.')
    scores = np.full(len(snp_ids), np.nan)
    metadata = {'bed_files': [], 'ldscore_files': []}
    seen_ids = set()
    for bed, requested_chrom in _bed_files(bed_pattern, requested):
        bim, chrom = _bim_metadata(bed, requested_chrom)
        rows = np.arange(len(snp_ids)) if unknown_chromosomes else np.flatnonzero(valid_chrom & (chromosomes == chrom))
        if not len(rows):
            continue
        print('Computing LD scores for chromosome ' + str(chrom), flush=True)
        # Export separately: the legacy writer formats its chromosome array
        # into the filename rather than using the scalar chromosome number.
        values, identifiers = ldscores_from_bed(str(bed), chrom, ld_wind, None)
        values = np.asarray(values, dtype=float)
        identifiers = np.asarray(identifiers, dtype=str)
        if values.ndim != 1 or identifiers.shape != values.shape:
            raise ValueError('LD-score backend returned incompatible SNP and score arrays.')
        if pd.Index(identifiers).has_duplicates:
            raise ValueError('LD-score backend returned duplicate SNP IDs for ' + str(bed))
        if not np.isin(identifiers, bim.index.to_numpy()).all():
            raise ValueError('LD-score backend returned SNP IDs absent from ' + str(bed.with_suffix('.bim')))
        if unknown_chromosomes and seen_ids.intersection(identifiers):
            raise ValueError('LD-score SNP IDs must be unique across BED files when summary chromosomes are unknown.')
        seen_ids.update(identifiers)
        lookup = pd.Series(values, index=identifiers)
        matched = pd.Index(snp_ids[rows]).isin(lookup.index)
        if not matched.any():
            raise ValueError('No SNPs in common between LD scores and summary statistics for chromosome ' + str(chrom))
        scores[rows[matched]] = lookup.reindex(snp_ids[rows[matched]]).to_numpy(dtype=float)
        metadata['bed_files'].append(str(bed))
        if ld_out is not None:
            output_pattern = str(ld_out).replace('#', '@')
            prefix = output_pattern.replace('@', str(chrom)) if '@' in output_pattern else output_pattern + str(chrom)
            output = Path(prefix + '.l2.ldscore.gz')
            output.parent.mkdir(parents=True, exist_ok=True)
            table = bim.loc[identifiers].reset_index()[['CHR', 'SNP', 'BP']]
            table['L2'] = values
            table.to_csv(output, sep='\t', index=False, compression='gzip', na_rep='NA')
            metadata['ldscore_files'].append(str(output))
    if not metadata['bed_files']:
        raise ValueError('No genotype chromosomes overlap the summary statistics.')
    return scores, metadata
