from concurrent.futures import ThreadPoolExecutor

import numpy as np
from snipar.gtarray import gtarray
from snipar.utilities import make_id_dict, open_bgen

from snipar.read.variants import bgen_variant_ids, match_variants, observed_variant_indices, sibling_sample_indices


def _match_variants(gts_f, par_gts_f, snp_ids=None, start=0, end=None):
    sid, aliases = bgen_variant_ids(gts_f)
    alleles = np.array([x.split(',') for x in gts_f.allele_ids])
    return match_variants(sid, np.asarray(gts_f.chromosomes), np.asarray(gts_f.positions),
                          alleles, par_gts_f, snp_ids, start, end, aliases)


def match_observed_and_imputed_snps(gts_f, par_gts_f, snp_ids=None, start=0, end=None):
    matched = _match_variants(gts_f, par_gts_f, snp_ids, start, end)
    return matched.with_imputed_mask(par_gts_f['bim_values'].shape[0])


def get_snps(gts_f, snp_ids=None, match=False):
    """
    Used in get_gts_matrix_given_ped to match observed and imputed SNPs and return SNP information on shared SNPs.
    Removes SNPs that have duplicated SNP ids.
    in_obs_sid contains the SNPs in the imputed genotypes that are present in the observed SNPs
    obs_sid_index contains the index in the observed SNPs of the common SNPs
    """
    obs_sid, _ = bgen_variant_ids(gts_f)
    obs_sid_index = observed_variant_indices(obs_sid, snp_ids, match)
    if not obs_sid_index.size:
        raise ValueError('No SNPs found in bgen file')
    alleles = np.array([pair.split(',') for pair in np.asarray(gts_f.allele_ids)[obs_sid_index]])
    return (np.asarray(gts_f.chromosomes)[obs_sid_index], obs_sid[obs_sid_index],
            np.asarray(gts_f.positions)[obs_sid_index], alleles, obs_sid_index)


def _first_allele_dosage(reader, samples, variants, *, num_threads=None, order=None):
    if hasattr(reader, 'read_dosage'):
        return reader.read_dosage(samples, variants, order=order or 'F')
    if hasattr(reader, '_vaddr'):
        samples, n_samples = _read_index(samples, len(reader.samples))
        variants, n_variants = _read_index(variants, len(reader._vaddr))
        result = np.empty((n_samples, n_variants), dtype=np.float32, order=order or 'F')

        def store(index, genotype):
            target = result[:, index]
            probabilities = genotype.probability
            p_first = probabilities[samples, 0]
            if genotype.phased:
                np.add(p_first, probabilities[samples, 2], out=target)
            else:
                np.multiply(p_first, 2, out=target)
                np.add(target, probabilities[samples, 1], out=target)
            target[genotype.missing[samples]] = np.nan

        _legacy_diploid(reader, samples, variants, num_threads or 1, store)
        return result
    kwargs = {} if num_threads is None else {'num_threads': num_threads}
    probabilities = reader.read((samples, variants), np.float32, **kwargs)
    result = np.empty(probabilities.shape[:2], dtype=np.float32, order=order or 'F')
    if probabilities.shape[2] == 4:
        np.add(probabilities[:, :, 0], probabilities[:, :, 2], out=result)
    elif probabilities.shape[2] == 3:
        np.multiply(probabilities[:, :, 0], 2, out=result)
        np.add(result, probabilities[:, :, 1], out=result)
    else:
        raise ValueError('Expected diploid, biallelic BGEN genotypes')
    return result


def bgen_is_phased(reader, variant):
    if hasattr(reader, 'is_phased'):
        return bool(reader.is_phased(variant))
    return bool(reader.phased[variant])


def read_diploid_probabilities(reader, samples, variants, phased, *, num_threads=1):
    """Read a hard-call tile and validate the genotype layout without scanning the file."""
    if hasattr(reader, 'is_phased'):
        return reader.read((samples, variants), np.float32, expected_phased=phased)
    if not hasattr(reader, '_vaddr'):
        return reader.read((samples, variants), np.float32)

    samples, n_samples = _read_index(samples, len(reader.samples))
    variants, n_variants = _read_index(variants, len(reader._vaddr))
    result = np.empty((n_samples, n_variants, 4 if phased else 3), dtype=np.float32, order='F')

    def store(index, genotype):
        if genotype.phased != phased:
            raise ValueError('BGEN hard calls require a consistently phased or unphased file')
        # Excluded samples may have higher ploidy and additional probability states.
        result[:, index, :] = genotype.probability[samples, :result.shape[2]]

    _legacy_diploid(reader, samples, variants, num_threads, store)
    return result


def read_hard_calls(reader, samples, variants, phased, missing, *, num_threads=1):
    """Write BGEN hard calls directly to a compact sample-by-variant matrix."""
    samples, n_samples = _read_index(samples, len(reader.samples))
    variant_count = len(reader.rsids)
    variants, n_variants = _read_index(variants, variant_count)
    result = np.empty((n_samples, n_variants), dtype=np.int8, order='F')
    n_probabilities = 4 if phased else 3
    if hasattr(reader, '_vaddr'):
        def store(index, genotype):
            if genotype.phased != phased:
                raise ValueError('BGEN hard calls require a consistently phased or unphased file')
            values = genotype.probability[samples, :n_probabilities]
            hard_calls(values, phased, missing, out=result[:, index])

        _legacy_diploid(reader, samples, variants, num_threads, store)
        return result

    # Keep the native probability tiles and temporary masks small regardless
    # of how many variants the caller needs to retain in the hard-call matrix.
    if isinstance(variants, slice):
        variants = np.arange(*variants.indices(variant_count), dtype=np.intp)
    elif variants.dtype.kind == 'b':
        variants = np.flatnonzero(variants)
    tile = max(1, 16 * 1024**2 // (4 * n_probabilities * max(1, n_samples)))
    for start in range(0, n_variants, tile):
        stop = min(start + tile, n_variants)
        probabilities = read_diploid_probabilities(reader, samples, variants[start:stop], phased,
                                                   num_threads=num_threads)
        hard_calls(probabilities, phased, missing, out=result[:, start:stop])
        del probabilities
    return result


def _read_index(index, size):
    if index is None:
        return slice(None), size
    if isinstance(index, slice):
        return index, len(range(*index.indices(size)))
    index = np.atleast_1d(index)
    return index, np.count_nonzero(index) if index.dtype.kind == 'b' else len(index)


def _legacy_diploid(reader, samples, variants, num_threads, store):
    """Validate selected variants during decode instead of decompressing a chromosome up front."""

    from cbgen import bgen_file

    offsets = reader._vaddr[variants]
    alleles = reader.nalleles[variants]
    if np.any(alleles != 2):
        raise ValueError('Expected diploid, biallelic BGEN genotypes')

    def fill(start, stop):
        with bgen_file(reader._filepath) as source:
            for column in range(start, stop):
                genotype = source.read_genotype(int(offsets[column]), precision=32)
                if np.any(genotype.ploidy[samples] != 2):
                    raise ValueError('Expected diploid, biallelic BGEN genotypes')
                store(column, genotype)

    threads = min(num_threads, len(offsets))
    if threads == 1:
        fill(0, len(offsets))
    elif threads > 1:
        boundaries = np.linspace(0, len(offsets), threads + 1, dtype=int)
        with ThreadPoolExecutor(max_workers=threads) as pool:
            jobs = [pool.submit(fill, start, stop) for start, stop in zip(boundaries[:-1], boundaries[1:])]
            for job in jobs:
                job.result()


def hard_calls(probabilities, phased, missing, *, out=None):
    """Keep genotype calls only when their posterior exceeds 0.99."""
    if out is None:
        out = np.empty(probabilities.shape[:-1], dtype=np.int8, order='F')
    threshold = np.float32(.99)
    if phased:
        first = probabilities[..., 0] > threshold
        second = probabilities[..., 2] > threshold
        called = first | (probabilities[..., 1] > threshold)
        called &= second | (probabilities[..., 3] > threshold)
        out[:] = first
        out += second
        out[~called] = missing
    else:
        out.fill(missing)
        for state, dose in enumerate((2, 1, 0)):
            np.copyto(out, dose, where=probabilities[..., state] > threshold)
    return out


def get_gts_matrix_given_ped(ped, imp_fams, bgenfile, par_gts_f=None ,snp_ids=None, ids=None, sib=False, sib_diff=False, parsum=False, start=0, end=None, include_unrel=False, robust=False, trios_sibs=False, match_snp_ids=False, verbose=False, print_sample_info = False, matched_variants=None):
    """Read one aligned BED/BGEN family design."""
    from snipar.read import prepare_gts_matrix

    with prepare_gts_matrix(
        ped=ped, imp_fams=imp_fams, bgenfile=bgenfile, par_gts_f=par_gts_f, ids=ids,
        sib=sib, sib_diff=sib_diff, parsum=parsum, include_unrel=include_unrel, robust=robust,
        trios_sibs=trios_sibs, verbose=verbose or print_sample_info, _split_pedigree=False,
    ) as reader:
        variants = matched_variants
        if variants is None:
            variants = reader.match_variants(snp_ids, match_snp_ids=match_snp_ids, start=start, end=end)
        return reader.read(variants)


def read_sibs_from_bgen(bgenfile, sibpairs):
    with open_bgen(bgenfile, verbose=True) as bgen:
        ids = bgen.samples
        snp_ids, _ = bgen_variant_ids(bgen)
        sibindices = sibling_sample_indices(sibpairs, make_id_dict(ids), 'BGEN')
        selected_ids = ids[sibindices]
        snp_ids = snp_ids.copy()
        positions = np.array(bgen.positions)
        gts = _first_allele_dosage(bgen, sibindices, None)
    return gtarray(garray=gts, ids=selected_ids, sid=snp_ids, pos=positions)


def read_PO_pairs_from_bgen(ped,bgenfile):
    with open_bgen(bgenfile, verbose=False) as bgen:
        return _read_PO_pairs(ped, bgen, bgenfile)


def _read_PO_pairs(ped, bgen, bgenfile):
    ids = bgen.samples
    id_dict = make_id_dict(ids)
    # SNP IDs
    snp_ids, _ = bgen_variant_ids(bgen)
    ## Find parent-offspring pairs
    # genotyped individuals
    genotyped = np.array([x in id_dict for x in ped[:, 1]])
    ped = ped[genotyped, :]
    # with genotyped father
    father_genotyped = np.array([x in id_dict for x in ped[:, 2]])
    # with genotyped mother
    mother_genotyped = np.array([x in id_dict for x in ped[:, 3]])
    # either
    opg = np.logical_or(father_genotyped, mother_genotyped)
    opg_ped = ped[opg, :]
    # number of pairs
    npair = np.sum(father_genotyped) + np.sum(mother_genotyped)
    if npair == 0:
        raise(ValueError('No parent-offspring pairs in  '+str(bgenfile)+' for genotype error probability estimation'))
    print(f'{npair} parent-offspring pairs found in {bgenfile}')
    if npair*snp_ids.shape[0] < 10**5:
        print('Warning: limited information for estimation of genotyping error probability.')
    ## Read genotypes
    all_ids = np.unique(np.hstack((opg_ped[:, 1],
                                   ped[father_genotyped, 2],
                                   ped[mother_genotyped, 3])))
    all_ids_indices = np.sort(np.array([id_dict[x] for x in all_ids]))
    gts = _first_allele_dosage(bgen, all_ids_indices, None)
    return gtarray(gts, ids=ids[all_ids_indices], sid=snp_ids.copy()), opg_ped, npair
