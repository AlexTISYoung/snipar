python bin/impute_from_sibs_setup.py build_ext --inplace
for i in {1..22}
do
	echo $i
	../.conda_envs/sibreg/bin/python bin/impute_runner.py --out  "test_data/imputation_result_chr_$i" $i "test_data/IBD.segments.gz" "../UKBioRDE_revision/data/tmp/filtered_ukb_chr$i" "test_data/pedigree" "test_data/parent_imputed_chr$i"
done
