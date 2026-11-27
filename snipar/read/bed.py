from pathlib import Path

import numpy as np
import pandas as pd
from pysnptools.snpreader import Bed
from snipar.gtarray import gtarray
from snipar.utilities import make_id_dict
from snipar.read.variants import match_variants, observed_variant_indices, sibling_sample_indices


def _read_iids(filename, decoder):
    """Load the IID column without retaining other local PLINK sample properties."""
    fam = Path(filename).with_suffix('.fam')
    if not fam.is_file() or fam.stat().st_size == 0:
        return decoder.iid
    column = pd.read_csv(fam, sep=r'\s+', header=None, usecols=[1], dtype=str, na_filter=False)
    return column[1].to_numpy(dtype=str)


def _match_variants(gts_f, par_gts_f, bim, snp_ids=None, start=0, end=None):
    alleles = np.loadtxt(bim, dtype=str, usecols=(4, 5), ndmin=2)
    return match_variants(gts_f.sid, gts_f.pos[:, 0].astype(int), gts_f.pos[:, 2].astype(int),
                          alleles, par_gts_f, snp_ids, start, end)


def match_observed_and_imputed_snps(gts_f, par_gts_f, bim, snp_ids=None, start=0, end=None):
    matched = _match_variants(gts_f, par_gts_f, bim, snp_ids, start, end)
    return matched.with_imputed_mask(par_gts_f['bim_values'].shape[0])


def get_snps(gts_f,bim,snp_ids=None, match=False):
    obs_sid = gts_f.sid
    obs_sid_index = observed_variant_indices(obs_sid, snp_ids, match)
    if not obs_sid_index.size:
        raise ValueError('No SNPs found in bed file')
    alleles = np.loadtxt(bim, dtype=str, usecols=(4, 5), ndmin=2)[obs_sid_index]
    positions = gts_f.pos[obs_sid_index]
    return (positions[:, 0].astype(int), obs_sid[obs_sid_index], positions[:, 2].astype(int),
            alleles, obs_sid_index)


def get_gts_matrix_given_ped(ped, imp_fams, bedfile, par_gts_f=None ,snp_ids=None, ids=None, sib=False, sib_diff=False, parsum=False, start=0, end=None, include_unrel=False, robust=False, trios_sibs=False, match_snp_ids=False, verbose=False, print_sample_info = False, matched_variants=None):
    """Read one aligned BED/BGEN family design."""
    from snipar.read import prepare_gts_matrix

    with prepare_gts_matrix(
        ped=ped, imp_fams=imp_fams, bedfile=bedfile, par_gts_f=par_gts_f, ids=ids,
        sib=sib, sib_diff=sib_diff, parsum=parsum, include_unrel=include_unrel, robust=robust,
        trios_sibs=trios_sibs, verbose=verbose or print_sample_info, _split_pedigree=False,
    ) as reader:
        variants = matched_variants
        if variants is None:
            variants = reader.match_variants(snp_ids, match_snp_ids=match_snp_ids, start=start, end=end)
        return reader.read(variants)


def read_sibs_from_bed(bedfile,sibpairs):
    bed = Bed(bedfile, count_A1=True)
    ids = bed.iid
    id_dict = make_id_dict(ids, 1)
    sibindices = sibling_sample_indices(sibpairs, id_dict, 'BED')
    gts = bed[sibindices, :].read(dtype=np.float32).val
    return gtarray(garray=gts, ids=ids[sibindices, 1], sid=bed.sid, pos=bed.pos[:, 2].astype(int))


def read_PO_pairs_from_bed(ped,bedfile):
    # Read bed
    bed = Bed(bedfile, count_A1=True)
    ids = bed.iid
    id_dict = make_id_dict(ids, 1)
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
        raise(ValueError('No parent-offspring pairs in  '+str(bedfile)+' for genotype error probability estimation'))
    print(str(npair)+' parent-offspring pairs found in '+bedfile)
    if npair*bed.sid.shape[0] < 10**5:
        print('Warning: limited information for estimation of genotyping error probability.')
    ## Read genotypes
    all_ids = np.unique(np.hstack((opg_ped[:, 1],
                                   ped[father_genotyped, 2],
                                   ped[mother_genotyped, 3])))
    all_ids_indices = np.sort(np.array([id_dict[x] for x in all_ids]))
    gts = bed[all_ids_indices, :].read().val
    return gtarray(gts,ids = ids[all_ids_indices, 1], sid=bed.sid), opg_ped, npair
