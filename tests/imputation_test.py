import h5py
import numpy as np
import pandas as pd
from pysnptools.snpreader import Bed
#testing the imputation result for whole genome
chromosomes = list(range(22,23))
chromosomes_expected_genes = []
chromosomes_imputed_genes = []
#TODO handle parent imputed name
for chromosome in chromosomes:
    with h5py.File('test_data/parent_imputed_chr'+str(chromosome),'r') as f:
        gts = np.array(f["imputed_par_gts"])
        fids = np.array(f["families"])
    expected = Bed("../UKBioRDE_revision/data/tmp/filtered_ukb_chr"+str(chromosome)+".bed")
    expected_gts = expected.read().val
    expected_ids = expected.iid
    iid_to_bed_index = {i:index for index, i in enumerate(expected_ids[:,0])}
    ped = pd.read_csv("../UKBioRDE_revision/data/tmp/pedigree", sep = " ")
    index_of_control_families_in_imputation = [index for index,fid in enumerate(fids) if fid.startswith("_")] 
    control_families = [fid[1:] for fid in set(ped["FID"].values.tolist()) if fid.startswith("_")]
    parents_of_control_families = ped.groupby("FID").agg({'FATHER_ID':lambda x: ([a for a in list(x) if not a.endswith("_P")]+[None])[0],
                                'MOTHER_ID':lambda x: ([a for a in list(x) if not a.endswith("_M")]+[None])[0],
                                }).loc[control_families]
    mother_indexes = [iid_to_bed_index[parents_of_control_families.loc[i]["MOTHER_ID"]] for i in control_families]
    father_indexes = [iid_to_bed_index[parents_of_control_families.loc[i]["FATHER_ID"]] for i in control_families]
    expected_parent_gts = expected_gts[mother_indexes,:]+expected_gts[father_indexes,:]
    expected_genes = expected_parent_gts.reshape((1,-1))
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

