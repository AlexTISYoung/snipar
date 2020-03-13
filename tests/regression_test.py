from bin.cython_impute_from_sibs import impute, prepare_data
import pandas as pd
import numpy as np
from pysnptools.snpreader import Bed
import os
import numpy as np
import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s')
logging.root.setLevel(logging.NOTSET)
# python example/simulate_pop.py 1000 0.05 40000 0 0 0 "test_data/generated"
os.system('python example/simulate_pop.py 1000 0.05 40000 0 0 0 "test_data/generated"')
# echo -e "FID IID FATHER_ID MOTHER_ID\n$(cat test_data/generated_fams.ped)" > test_data/generated_fams.ped
os.system('echo -e "FID IID FATHER_ID MOTHER_ID\n$(cat test_data/generated_fams.ped)" > test_data/generated_fams.ped')
# plink1 --file test_data/generated --make-bed --out test_data/generated
os.system('plink1 --file test_data/generated --make-bed --out test_data/generated')


columns = ["FID", "IID", "FATHER_ID", "MOTHER_ID", "sex", "phenotype"] + ["genotype_"+str(i) for i in range(1000)]
ped = pd.read_csv("test_data/generated.ped", sep = " ", header = None, names = columns)
ped = ped[["FID", "IID", "FATHER_ID", "MOTHER_ID"]]
fathers = ped["IID"].str.endswith("_P")
mothers = ped["IID"].str.endswith("_M")
parents = ped[mothers | fathers]
sibs = ped[(~mothers) & (~fathers)]
sibs.to_csv("test_data/generated_sibs.ped", sep = " ")
parents.to_csv("test_data/generated_parents.ped", sep = " ")
with open("test_data/generated_sibs.txt", "w") as f:
	for i, j in sibs[["FID", "IID"]].values.tolist():
		f.write(str(i) + " " + str(j) + "\n")
# plink1 --bfile test_data/generated --keep test_data/generated_sibs.txt --make-bed --out test_data/generated_sibs
os.system('plink1 --bfile test_data/generated --keep test_data/generated_sibs.txt --make-bed --out test_data/generated_sibs')
# plink1 --bfile test_data/generated --remove test_data/generated_sibs.txt --make-bed --out test_data/generated_parents
os.system('plink1 --bfile test_data/generated --remove test_data/generated_sibs.txt --make-bed --out test_data/generated_parents')
sibships, iid_to_bed_index, gts, ibd, pos = prepare_data("test_data/generated_sibs.ped",
														"test_data/generated_sibs",
														"test_data/generated.segments.gz",
														1)
gts = gts.astype(float)
pos = pos.astype(int)
imputed_fids, imputed_par_gts = impute(sibships, iid_to_bed_index, gts, ibd, pos)
expected_parents = Bed("test_data/generated_parents.bed")
expected_parents_gts = expected_parents.read().val
expected_parents_ids = expected_parents.iid
expected_parents_sum = expected_parents_gts[[bool(i%2) for i in range(80000)]] + expected_parents_gts[[not bool(i%2) for i in range(80000)]]
expected_parents_sum = expected_parents_sum[[int(t) for t in imputed_fids],:]
expected_genotypes = expected_parents_sum.reshape((1,-1))
imputed_genotypes = imputed_par_gts.reshape((1,-1))
covs = np.cov(expected_genotypes, imputed_genotypes)
coef = covs[0,1]/covs[1,1]
residual_var = np.var(expected_genotypes - coef*imputed_genotypes)
s2 = residual_var/(1000*covs[1,1])
z = (1-coef)/np.sqrt(s2)


#TODO handle file addresses in config