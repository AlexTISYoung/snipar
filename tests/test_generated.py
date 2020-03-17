import sys
from sibreg.bin.impute_from_sibs import impute, prepare_data	
import pandas as pd
import numpy as np
from pysnptools.snpreader import Bed
import os
import numpy as np
import logging
import unittest
from scipy.stats import norm

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s')
logging.root.setLevel(logging.NOTSET)
class TestGenerated(unittest.TestCase):
	def test_generate_and_regress(self):
		#requies plink
		#TODO add plink to the package requirments
		#TODO test on windows
		number_of_snps = 1000
		min_f = 0.05
		number_of_families = 100
		filename = "test_data/generated"
		if not os.path.exists(os.path.dirname(filename)):
			try:
				os.makedirs(os.path.dirname(filename))
			except OSError as exc: # Guard against race condition
				if exc.errno != errno.EEXIST:
					raise
		#generating population
		os.system('python example/simulate_pop.py ' + str(number_of_snps) + ' ' + str(min_f) + ' ' + str(number_of_families) + ' 0 0 0 "test_data/generated"')
		# Adding header to the pedigree file
		os.system('echo -e "FID IID FATHER_ID MOTHER_ID\n$(cat test_data/generated_fams.ped)" > test_data/generated_fams.ped')
		#convert the generated data to a bed file
		os.system('plink/plink --noweb --file test_data/generated --make-bed --out test_data/generated')
		columns = ["FID", "IID", "FATHER_ID", "MOTHER_ID", "sex", "phenotype"] + ["genotype_"+str(i) for i in range(number_of_snps)]
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
		#writing sibs only
		#TODO handle plink path
		os.system('plink/plink --noweb --bfile test_data/generated --keep test_data/generated_sibs.txt --make-bed --out test_data/generated_sibs')
		#writing parents only
		os.system('plink/plink --noweb --bfile test_data/generated --remove test_data/generated_sibs.txt --make-bed --out test_data/generated_parents')
		sibships, iid_to_bed_index, gts, ibd, pos, sid = prepare_data("test_data/generated_sibs.ped",
																"test_data/generated_sibs",
																"test_data/generated.segments.gz",
																1)
		gts = gts.astype(float)
		pos = pos.astype(int)
		imputed_fids, imputed_par_gts = impute(sibships, iid_to_bed_index, gts, ibd, pos, sid)
		expected_parents = Bed("test_data/generated_parents.bed")
		expected_parents_gts = expected_parents.read().val
		expected_parents_ids = expected_parents.iid
		parent1 = expected_parents_gts[[bool(i%2) for i in range(2*number_of_families)]]
		parent2 = expected_parents_gts[[not bool(i%2) for i in range(2*number_of_families)]]
		expected_parents_sum = parent1+parent2							   
		expected_parents_sum = expected_parents_sum[[int(t) for t in imputed_fids],:]
		expected_genotypes = expected_parents_sum.reshape((1,-1))
		imputed_genotypes = imputed_par_gts.reshape((1,-1))
		covs = np.cov(expected_genotypes, imputed_genotypes)
		coef = covs[0,1]/covs[1,1]
		residual_var = np.var(expected_genotypes - coef*imputed_genotypes)
		s2 = residual_var/(number_of_snps*covs[1,1])
		#TODO should i divide by number_of_snps*covs[1,1]*number_of_families
		z = (1-coef)/np.sqrt(s2)
		q = norm.cdf(z)
		p_val = min(q, 1-q)
		self.assertGreaterEqual(p_val, 0.01)
