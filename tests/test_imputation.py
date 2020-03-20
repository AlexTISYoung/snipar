import h5py
import numpy as np
import pandas as pd
from pysnptools.snpreader import Bed
from scipy.stats import norm
#testing the imputation result for whole genome
def imputation_test(chromosomes,
                   imputed_prefix = 'test_data/parent_imputed_chr',
                   expected_prefix = "../UKBioRDE_revision/data/tmp/filtered_ukb_chr",
                   pedigree_address = "../UKBioRDE_revision/data/tmp/pedigree"
                   ):
    #Data files for chromosome i should be named in this fashion: "prefix{i}"
    chromosomes_expected_genes = []
    chromosomes_imputed_genes = []
    for chromosome in chromosomes:
        with h5py.File(imputed_prefix+str(chromosome),'r') as f:
            gts = np.array(f["imputed_par_gts"])
            fids = np.array(f["families"])
        expected = Bed(expected_prefix+str(chromosome)+".bed")
        expected_gts = expected.read().val
        expected_ids = expected.iid
        iid_to_bed_index = {i:index for index, i in enumerate(expected_ids[:,1])}
        ped = pd.read_csv(pedigree_address, sep = " ").astype(str)
        #fids of control families start with _
        index_of_families_in_imputation = {fid[1:]:index for index,fid in enumerate(fids)}
        control_families = [fid[1:] for fid in set(ped["FID"].values.tolist()) if fid.startswith("_")]
        #for each family select id of the parents
        parents_of_control_families = ped.groupby("FID").agg({
                                    'FATHER_ID':lambda x: ([a for a in list(x) if a in ped["IID"].tolist()]+[None])[0],
                                    'MOTHER_ID':lambda x: ([a for a in list(x) if a in ped["IID"].tolist()]+[None])[0],
                                    }).loc[control_families]
        mother_indexes = [iid_to_bed_index[parents_of_control_families.loc[i]["MOTHER_ID"]] for i in control_families]
        father_indexes = [iid_to_bed_index[parents_of_control_families.loc[i]["FATHER_ID"]] for i in control_families]
        expected_parent_gts = expected_gts[mother_indexes,:]+expected_gts[father_indexes,:]
        expected_genes = expected_parent_gts.reshape((1,-1))
        index_of_control_families_in_imputation = [index_of_families_in_imputation[i] for i in control_families]
        imputed_genes = gts[index_of_control_families_in_imputation,:].reshape((1,-1))
        mask = ~(np.isnan(expected_genes) | np.isnan(imputed_genes))
        expected_genes = expected_genes[mask]
        imputed_genes = imputed_genes[mask]
        chromosomes_expected_genes.append(expected_genes)
        chromosomes_imputed_genes.append(imputed_genes)
    whole_expected_genes = np.concatenate(chromosomes_expected_genes)
    whole_imputed_genes = np.concatenate(chromosomes_imputed_genes)
    covs = np.cov(whole_expected_genes, whole_imputed_genes)
    coef = covs[0,1]/covs[1,1]
    residual_var = np.var(whole_expected_genes - coef*whole_imputed_genes)
    s2 = residual_var/(len(control_families)*22*2*covs[1,1])
    z = (1-coef)/np.sqrt(s2)
    q = norm.cdf(z)
    p_value = min(q, 1-q)
    #TODO compute z correctly(find the correct sd)
    return coef, z, p_value
