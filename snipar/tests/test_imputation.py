import h5py, code
import numpy as np
import pandas as pd
from pysnptools.snpreader import Bed
from scipy.stats import norm
import logging
#testing the imputation result for whole genome
def imputation_test(chromosomes,
                   imputed_prefix = '',
                   expected_prefix = "../UKBioRDE_revision/data/tmp/filtered_ukb_chr",
                   start = None,
                   end = None
                   ):
    #Data files for chromosome i should be named in this fashion: "prefix{i}"
    chromosomes_expected_genes_o = []
    chromosomes_expected_genes_pm = []
    chromosomes_imputed_genes_o = []
    chromosomes_imputed_genes_pm = []
    for chromosome in chromosomes:
        with h5py.File(imputed_prefix+str(chromosome)+".hdf5",'r') as f:
            non_duplicates = np.array(f["non_duplicates"])
            if start is not None:
                gts = np.array(f["imputed_par_gts"])[:,start:end]
                non_duplicates = non_duplicates[start:end]+start
            else:
                gts = np.array(f["imputed_par_gts"])
            fids = np.array(f["families"]).astype(str)
            ped_array = np.array(f["pedigree"]).astype(str)
            ped = pd.DataFrame(ped_array[1:], columns = ped_array[0])            
        iids_in_ped = set(ped["IID"])
        expected = Bed(expected_prefix+str(chromosome)+".bed", count_A1 = True)
        expected_gts = expected[:, non_duplicates].read().val
        expected_ids = expected.iid
        iid_to_bed_index = {i:index for index, i in enumerate(expected_ids[:,1])}
        #fids of control families start with _
        #this has the predix _*_
        index_of_families_in_imputation = {fid:index for index,fid in enumerate(fids)}
        # no parent control starts with _o_
        # only has father control starts with _p_
        # only has father control starts with _m_
        control_o_families = set()
        control_p = set()
        control_m = set()
        for index, row in ped.iterrows():
            if row["FID"].startswith("_o_"):
                control_o_families.add(row["FID"][3:])
            if row["FID"].startswith("_p_"):
                control_p.add(row["FID"][3:])
            if row["FID"].startswith("_m_"):
                control_m.add(row["FID"][3:])        
        control_o_families = list(control_o_families)
        control_p = list(control_p)
        control_m = list(control_m)
        control_pm_families = control_p + control_m
        #for each family select id of the parents
        ped = ped[ped["FID"].isin(control_o_families) | ped["FID"].isin(control_pm_families)]
        parent_ids = ped.groupby("FID").agg({
                                    'FATHER_ID':lambda x: ([a for a in list(x) if a in iids_in_ped]+[None])[0],
                                    'MOTHER_ID':lambda x: ([a for a in list(x) if a in iids_in_ped]+[None])[0],
                                    })
        parents_of_control_o_families = parent_ids.loc[control_o_families]
        mother_indexes_control_o = [iid_to_bed_index[parents_of_control_o_families.loc[i, "MOTHER_ID"]] for i in control_o_families]
        father_indexes_control_o = [iid_to_bed_index[parents_of_control_o_families.loc[i, "FATHER_ID"]] for i in control_o_families]
        expected_parent_gts_control_o = (expected_gts[mother_indexes_control_o,:]+expected_gts[father_indexes_control_o,:])/2
        expected_genes_o = expected_parent_gts_control_o.reshape((1,-1))
        index_of_control_families_in_imputation_o = [index_of_families_in_imputation["_o_"+i] for i in control_o_families]
        imputed_genes_o = gts[index_of_control_families_in_imputation_o,:].reshape((1,-1))
        mask_o = ~(np.isnan(expected_genes_o) | np.isnan(imputed_genes_o))
        expected_genes_o = expected_genes_o[mask_o]
        imputed_genes_o = imputed_genes_o[mask_o]
        parent_of_control_m = parent_ids.loc[control_m]
        parent_of_control_p = parent_ids.loc[control_p]
        father_indexes_control_m = [iid_to_bed_index[parent_of_control_m.loc[i, "FATHER_ID"]] for i in control_m]
        mother_indexes_control_p = [iid_to_bed_index[parent_of_control_p.loc[i, "MOTHER_ID"]] for i in control_p]
        expected_parent_gts_control_pm = expected_gts[mother_indexes_control_p + father_indexes_control_m, :]
        expected_genes_pm = expected_parent_gts_control_pm.reshape((1,-1))
        index_of_control_families_in_imputation_pm = [index_of_families_in_imputation["_p_" + i] for i in control_p] + [index_of_families_in_imputation["_m_" + i] for i in control_m]
        imputed_genes_pm = gts[index_of_control_families_in_imputation_pm,:].reshape((1,-1))
        mask_pm = ~(np.isnan(expected_genes_pm) | np.isnan(imputed_genes_pm))
        expected_genes_pm = expected_genes_pm[mask_pm]
        imputed_genes_pm = imputed_genes_pm[mask_pm]
        chromosomes_expected_genes_o.append(expected_genes_o)
        chromosomes_expected_genes_pm.append(expected_genes_pm)
        chromosomes_imputed_genes_o.append(imputed_genes_o)
        chromosomes_imputed_genes_pm.append(imputed_genes_pm)
        
    whole_expected_genes_o = np.concatenate(chromosomes_expected_genes_o)
    whole_imputed_genes_o = np.concatenate(chromosomes_imputed_genes_o)
    whole_expected_genes_pm = np.concatenate(chromosomes_expected_genes_pm)
    whole_imputed_genes_pm = np.concatenate(chromosomes_imputed_genes_pm)

    covs_o = np.cov(whole_expected_genes_o, whole_imputed_genes_o)
    # print("whole_expected_genes_o")
    # print(whole_expected_genes_o)
    # print("whole_imputed_genes_o")
    # print(whole_imputed_genes_o)
    coef_o = covs_o[0,1]/covs_o[1,1]
    residual_var_o = np.var(whole_expected_genes_o - coef_o*whole_imputed_genes_o)
    s2_o = residual_var_o/(len(control_o_families)*22*2*covs_o[1,1])
    z_o = (1-coef_o)/np.sqrt(s2_o)
    q_o = norm.cdf(z_o)
    p_value_o = min(q_o, 1-q_o)
    # print("whole_expected_genes_pm")
    # print(whole_expected_genes_pm)
    # print("whole_imputed_genes_pm")
    # print(whole_imputed_genes_pm)
    # print("pm corr", np.corrcoef(whole_expected_genes_pm, whole_imputed_genes_pm))
    # print("o corr", np.corrcoef(whole_expected_genes_o, whole_imputed_genes_o))
    covs_pm = np.cov(whole_expected_genes_pm, whole_imputed_genes_pm)
    coef_pm = covs_pm[0,1]/covs_pm[1,1]
    residual_var_pm = np.var(whole_expected_genes_pm - coef_pm*whole_imputed_genes_pm)
    s2_pm = residual_var_pm/(len(control_pm_families)*22*2*covs_pm[1,1])
    z_pm = (1-coef_pm)/np.sqrt(s2_pm)
    q_pm = norm.cdf(z_pm)
    p_value_pm = min(q_pm, 1-q_pm)
    logging.info(f"returning {((coef_o, coef_pm), (z_o, z_pm), (p_value_o, p_value_pm))}")
    #TODO compute z correctly(find the correct sd)
    return (coef_o, coef_pm), (z_o, z_pm), (p_value_o, p_value_pm)

