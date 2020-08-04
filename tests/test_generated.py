import sys
from sibreg.bin.impute_from_sibs import impute
from sibreg.bin.preprocess_data import prepare_data
import pandas as pd
import numpy as np
from pysnptools.snpreader import Bed
import os
import numpy as np
import logging
import unittest
from scipy.stats import norm
import time

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s')
logging.root.setLevel(logging.NOTSET)
class TestGenerated(unittest.TestCase):
	def test_generate_and_regress(self):
		#requies plink
		number_of_snps = 1000
		min_f = 0.05
		number_of_families = 100
		filename = "outputs/tmp/generated"
		if not os.path.exists(os.path.dirname(filename)):
			try:
				os.makedirs(os.path.dirname(filename))
			except OSError as exc: # Guard against race condition
				if exc.errno != errno.EEXIST:
					raise
		start = 0
		interval = 0
		#generating population
		os.system('python example/simulate_pop.py ' + str(number_of_snps) + ' ' + str(min_f) + ' ' + str(number_of_families) + ' 0 0 0 "outputs/tmp/generated"')
		# Adding header to the ped.igree file
		os.system('echo -e "FID IID FATHER_ID MOTHER_ID\n$(cat outputs/tmp/generated_fams.ped)" > outputs/tmp/generated_fams.ped')
		#convert the generated data to a bed file
		os.system('plink/plink --noweb --file outputs/tmp/generated --make-bed --out outputs/tmp/generated')
		columns = ["FID", "IID", "FATHER_ID", "MOTHER_ID", "sex", "phenotype"] + ["genotype_"+str(i) for i in range(number_of_snps)]
		ped = pd.read_csv("outputs/tmp/generated.ped", sep = " ", header = None, names = columns)
		ped = ped[["FID", "IID", "FATHER_ID", "MOTHER_ID"]].astype(str)
		only_remove_father_ids = [str(i)+"_P" for i in range(number_of_families//4)]
		only_remove_mother_ids = [str(i)+"_M" for i in range(number_of_families//4, number_of_families//2)]
		remove_both_parents_ids = [str(i)+"_M" for i in range(number_of_families//2, number_of_families)] + [str(i)+"_P" for i in range(number_of_families//2, number_of_families)]
		parents = ped[ped["IID"].str.endswith("_P") | ped["IID"].str.endswith("_M")]
		sibs = ped[~ped["IID"].isin(only_remove_father_ids+only_remove_mother_ids+remove_both_parents_ids)]
		sibs.to_csv("outputs/tmp/generated_sibs.ped", sep = " ")
		parents.to_csv("outputs/tmp/generated_parents.ped", sep = " ")
		with open("outputs/tmp/generated_sibs.txt", "w") as f:
			for i, j in sibs[["FID", "IID"]].values.tolist():
				f.write(str(i) + " " + str(j) + "\n")
		
		with open("outputs/tmp/generated_parents.txt", "w") as f:
			for i, j in ped[["FID", "IID"]].values.tolist():
				if j.endswith("_P") or j.endswith("_M"):
					f.write(str(i) + " " + str(j) + "\n")
		#writing sibs only
		os.system('plink/plink --noweb --bfile outputs/tmp/generated --keep outputs/tmp/generated_sibs.txt --make-bed --out outputs/tmp/generated_sibs')
		#writing parents only
		os.system('plink/plink --noweb --bfile outputs/tmp/generated --keep outputs/tmp/generated_parents.txt --make-bed --out outputs/tmp/generated_parents')
		ibd = pd.read_csv("outputs/tmp/generated.segments.gz", sep = "\t")
		sibships, iid_to_bed_index, phased_gts, unphased_gts, ibd, pos, chromosomes, hdf5_output_dict = prepare_data(sibs,
																													None,
																													"outputs/tmp/generated_sibs",
																													ibd)
		print("unphased_gts", str(unphased_gts.astype("i"))[:200])
		print("pos", str(pos.astype(int))[:200])
		print("sibships", str(sibships)[:200])
		print("iid_to_bed_index", str(iid_to_bed_index)[:200])
		print("phased_gts", str(phased_gts)[:200])
		print("unphased_gts", str(unphased_gts)[:200])
		print("ibd", str(ibd)[:200])
		print("pos", str(pos)[:200])
		print("hdf5_output_dict", str(hdf5_output_dict)[:200])
		print("chromosomes", str(chromosomes)[:200])
		imputed_fids, imputed_par_gts = impute(sibships, iid_to_bed_index, phased_gts, unphased_gts, ibd, pos, hdf5_output_dict, str(chromosomes), threads = 2)
		print("imputed_par_gts", imputed_par_gts)
		expected_parents = Bed("outputs/tmp/generated_parents.bed", count_A1 = True)
		expected_parents_gts = expected_parents.read().val
		expected_parents_ids = expected_parents.iid
		father = expected_parents_gts[[bool(i%2) for i in range(2*number_of_families)]]
		father = father[[int(t) for t in imputed_fids],:]
		mother = expected_parents_gts[[not bool(i%2) for i in range(2*number_of_families)]]
		mother = mother[[int(t) for t in imputed_fids],:]		
		expected_parents = np.zeros(imputed_par_gts.shape)
		no_parent = ~sibships["has_father"] & ~sibships["has_mother"]
		only_mother = ~sibships["has_father"] & sibships["has_mother"]
		only_father = ~sibships["has_mother"] & sibships["has_father"]
		expected_parents[no_parent] = (mother[no_parent] + father[no_parent])/2
		expected_parents[only_mother] = mother[only_mother]
		expected_parents[only_father] = father[only_father]
		expected_genotypes = expected_parents.reshape((1,-1))
		imputed_genotypes = imputed_par_gts.reshape((1,-1))
		covs = np.cov(expected_genotypes, imputed_genotypes)
		coef = covs[0,1]/covs[1,1]
		print(covs)
		residual_var = np.var(expected_genotypes - coef*imputed_genotypes)
		s2 = residual_var/(number_of_snps*covs[1,1])
		#TODO should i divide by number_of_snps*covs[1,1]*number_of_families
		z = (1-coef)/np.sqrt(s2)
		q = norm.cdf(z)
		p_val = min(q, 1-q)
		self.assertGreaterEqual(p_val, 0.01)

	# def test_empty(self):
	# 	#requies plink
	# 	number_of_snps = 10
	# 	min_f = 0.3
	# 	number_of_families = 1000
	# 	filename = "outputs/tmp/generated"
	# 	if not os.path.exists(os.path.dirname(filename)):
	# 		try:
	# 			os.makedirs(os.path.dirname(filename))
	# 		except OSError as exc: # Guard against race condition
	# 			if exc.errno != errno.EEXIST:
	# 				raise
	# 	start = 0
	# 	interval = 0
	# 	#generating population
	# 	os.system('python example/simulate_pop.py ' + str(number_of_snps) + ' ' + str(min_f) + ' ' + str(number_of_families) + ' 0 0 0 "outputs/tmp/generated"')
	# 	# Adding header to the ped.igree file
	# 	os.system('echo -e "FID IID FATHER_ID MOTHER_ID\n$(cat outputs/tmp/generated_fams.ped)" > outputs/tmp/generated_fams.ped')
	# 	#convert the generated data to a bed file
	# 	os.system('plink/plink --noweb --file outputs/tmp/generated --make-bed --out outputs/tmp/generated')
	# 	columns = ["FID", "IID", "FATHER_ID", "MOTHER_ID", "sex", "phenotype"] + ["genotype_"+str(i) for i in range(number_of_snps)]
	# 	ped = pd.read_csv("outputs/tmp/generated.ped", sep = " ", header = None, names = columns)
	# 	ped = ped[["FID", "IID", "FATHER_ID", "MOTHER_ID"]].astype(str)
	# 	only_remove_father_ids = [str(i)+"_P" for i in range(number_of_families//4)]
	# 	only_remove_mother_ids = [str(i)+"_M" for i in range(number_of_families//4, number_of_families//2)]
	# 	remove_both_parents_ids = [str(i)+"_M" for i in range(number_of_families//2, number_of_families)] + [str(i)+"_P" for i in range(number_of_families//2, number_of_families)]
	# 	# only_remove_father_ids = [str(i)+"_M" for i in range(number_of_families)]
	# 	# only_remove_mother_ids = [str(i)+"_P" for i in range(number_of_families)]
	# 	# remove_both_parents_ids = []
	# 	parents = ped[ped["IID"].str.endswith("_P") | ped["IID"].str.endswith("_M")]
	# 	sibs = ped[~ped["IID"].isin(only_remove_father_ids+only_remove_mother_ids+remove_both_parents_ids)]
	# 	sibs.to_csv("outputs/tmp/generated_sibs.ped", sep = " ")
	# 	parents.to_csv("outputs/tmp/generated_parents.ped", sep = " ")
	# 	with open("outputs/tmp/generated_sibs.txt", "w") as f:
	# 		for i, j in sibs[["FID", "IID"]].values.tolist():
	# 			f.write(str(i) + " " + str(j) + "\n")
		
	# 	with open("outputs/tmp/generated_parents.txt", "w") as f:
	# 		for i, j in ped[["FID", "IID"]].values.tolist():
	# 			if j.endswith("_P") or j.endswith("_M"):
	# 				f.write(str(i) + " " + str(j) + "\n")
	# 	#writing sibs only
	# 	os.system('plink/plink --noweb --bfile outputs/tmp/generated --keep outputs/tmp/generated_sibs.txt --make-bed --out outputs/tmp/generated_sibs')
	# 	#writing parents only
	# 	os.system('plink/plink --noweb --bfile outputs/tmp/generated --keep outputs/tmp/generated_parents.txt --make-bed --out outputs/tmp/generated_parents')
	# 	ibd = pd.read_csv("outputs/tmp/generated.segments.gz", sep = "\t")
	# 	print(ibd.columns)
	# 	print(ibd)

	# 	sibships, iid_to_bed_index, phased_gts, unphased_gts, ibd, pos, chromosomes, hdf5_output_dict = prepare_data(sibs,
	# 																												None,
	# 																												"outputs/tmp/generated_sibs",
	# 																												ibd)
	# 	print("unphased_gts", str(unphased_gts.astype("i"))[:200])
	# 	print("pos", str(pos.astype(int))[:200])
	# 	print("sibships", str(sibships)[:200])
	# 	print("iid_to_bed_index", str(iid_to_bed_index)[:200])
	# 	print("phased_gts", str(phased_gts)[:200])
	# 	print("unphased_gts", str(unphased_gts)[:200])
	# 	print("ibd", str(ibd)[:200])
	# 	print("pos", str(pos)[:200])
	# 	print("hdf5_output_dict", str(hdf5_output_dict)[:200])
	# 	print("chromosomes", str(chromosomes)[:200])
	# 	imputed_fids, imputed_par_gts = impute(sibships, iid_to_bed_index, phased_gts, unphased_gts, ibd, pos, hdf5_output_dict, str(chromosomes), threads = 2)
	# 	print("imputed_par_gts", imputed_par_gts)
	# 	print("imputed_par_gts sum", np.sum(imputed_par_gts))
	# 	expected_parents = Bed("outputs/tmp/generated_parents.bed", count_A1 = True)
	# 	expected_parents_gts = expected_parents.read().val
	# 	expected_parents_ids = expected_parents.iid
	# 	print("expected_gts")
	# 	print(expected_parents_gts)
	# 	print(np.sum(expected_parents_gts))
	# 	father = expected_parents_gts[[bool(i%2) for i in range(2*number_of_families)]]
	# 	father = father[[int(t) for t in imputed_fids],:]
	# 	mother = expected_parents_gts[[not bool(i%2) for i in range(2*number_of_families)]]
	# 	mother = mother[[int(t) for t in imputed_fids],:]		
	# 	expected_parents = np.zeros(imputed_par_gts.shape)
	# 	no_parent = ~sibships["has_father"] & ~sibships["has_mother"]
	# 	only_mother = ~sibships["has_father"] & sibships["has_mother"]
	# 	only_father = ~sibships["has_mother"] & sibships["has_father"]

	# 	expected_parents[no_parent] = (mother[no_parent] + father[no_parent])/2
	# 	expected_parents[only_mother] = mother[only_mother]
	# 	expected_parents[only_father] = father[only_father]

	# 	# expected_parents[no_parent] = (mother[no_parent] + father[no_parent])/2
	# 	# expected_parents = expected_parents[no_parent]
	# 	# imputed_par_gts = imputed_par_gts[no_parent]


	# 	# expected_parents[only_mother] = mother[only_mother]
	# 	# expected_parents = expected_parents[only_mother]
	# 	# imputed_par_gts = imputed_par_gts[only_mother]


	# 	# expected_parents[only_father] = mother[only_father]
	# 	# expected_parents = expected_parents[only_father]
	# 	# imputed_par_gts = imputed_par_gts[only_father]
	# 	# print("shape", imputed_par_gts.shape)

	# 	expected_genotypes = expected_parents.reshape((1,-1))
	# 	imputed_genotypes = imputed_par_gts.reshape((1,-1))
	# 	s = (~np.isnan(expected_genotypes)) & (~np.isnan(imputed_genotypes))
	# 	expected_genotypes = expected_genotypes[s]
	# 	imputed_genotypes = imputed_genotypes[ s]
	# 	covs = np.cov(expected_genotypes, imputed_genotypes)
	# 	coef = covs[0,1]/covs[1,1]
	# 	print(covs)
	# 	print("coef", coef)
	# 	residual_var = np.var(expected_genotypes - coef*imputed_genotypes)
	# 	s2 = residual_var/(number_of_snps*covs[1,1])
	# 	#TODO should i divide by number_of_snps*covs[1,1]*number_of_families
	# 	z = (1-coef)/np.sqrt(s2)
	# 	q = norm.cdf(z)
	# 	p_val = min(q, 1-q)
	# 	self.assertGreaterEqual(p_val, 0.01)
