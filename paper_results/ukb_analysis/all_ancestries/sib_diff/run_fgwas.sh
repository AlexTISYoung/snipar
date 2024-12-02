source ${snipar_env}/bin/activate

set -e

for i in {1..10} {14..18}; do
	name=$(head -n 1 ${PHEN_PATH}/processed_traits_noadj_allancestry.txt | awk -v var=$i '{ print $((2+var)) }')
	if [ ! -d "$outdir/$name" ] 
	then
		mkdir $outdir/$name
	fi
	gwas.py \
		${PHEN_PATH}/processed_traits_noadj_allancestry.txt \
		--imp ${IMP_PATH}/all_ancestries/chr@ \
		--bgen ${BGEN_PATH}/all_ancestries/ukb_hap_chr@_v2_filtered \
		--covar ${PHEN_PATH}/covariates_allancestry.txt \
		--chr_range 1-22 \
		--phen_index $i \
		--cpus 4 \
		--ibdrel_path ${KING_IBDSEG_PATH}/ukb_v2_3rd \
		--sib_diff \
		 --out $outdir/$name/@
done