path2ldscreg="/homes/nber/harij/gitrepos/SNIPar/ldsc_reg"

dirs=( "from_chr1_to_chr23_start0_end50_run0_p0-0_ab_corr0-0_vb0-25_length2"
"from_chr1_to_chr23_start0_end50_run0_p0-0_ab_corr1-0_vb0-25_length2"
"from_chr1_to_chr23_start0_end50_run0_p0-5_ab_corr0-0_vb0-25_length2"
"from_chr1_to_chr23_start0_end50_run0_p0-5_ab_corr0-5_vb0-25_length2"
"from_chr1_to_chr23_start0_end50_run0_p0-5_ab_corr1-0_vb0-25_length2"
"from_chr1_to_chr3_start0_end20_run0_p0-0_ab_corr1-0_vb0-25_length2") 

for dir in "${dirs[@]}"
do
    echo ${dir}
    
    python ${path2ldscreg}/run_estimates.py /disk/genetics/ukb/alextisyoung/haplotypes/simulated_pops/${dir}/gen_0_gen_1_phenotype.hdf5 \
            -ldsc "/disk/genetics/ukb/alextisyoung/haplotypes/simulated_pops/${dir}/ldscores/*[0-9].l2.ldscore.gz" \
            -l "/homes/nber/harij/gitrepos/SNIPar/ldsc_reg/multigen_sim/${dir}.log" \
            --jkse \
            --jkse_blocksize 20 \
            --jkse_cores 4
done