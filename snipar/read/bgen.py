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
