source ${snipar_env}/bin/activate

set -e


for i in {1..10} {14..17} 19; do
	name=$(head -n 1 ${PHEN_PATH}/processed_traits_noadj.txt | awk -v var=$i '{ print $((2+var)) }')
	if [ ! -d "$outdir/$name" ] 
	then
		mkdir $outdir/$name
	fi
	gwas.py \
		${PHEN_PATH}/processed_traits_noadj.txt \
		--imp ${IMP_PATH}/european/chr_@ \
		--bgen ${BGEN_PATH}/european/ukb_hap_chr@_v2_filtered \
		--covar ${PHEN_PATH}/covariates.txt \
		--ibdrel_path ${KING_IBDSEG_PATH}/ukb_v2_3rd \
		--chr_range 1-22 \
		--phen_index $i \
		 --out $outdir/$name/@

   done

