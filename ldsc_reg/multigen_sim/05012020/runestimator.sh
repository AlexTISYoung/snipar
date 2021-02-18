path2ldscreg="/homes/nber/harij/gitrepos/SNIPar/ldsc_reg"

for dir in /disk/genetics/ukb/alextisyoung/haplotypes/simulated_pops_large/from_chr1_to_chr23_start0_endNone_run0_p0-0_ab_corr0-5_vb0-25_length2/phenotype_dir_par_corr_0.5/*/
do

    runno=$(echo ${dir} | cut -d'/' -f 10)

    echo "Run Number: $runno"

    echo "Population Effect"
    python ${path2ldscreg}/run_estimates.py "${dir}/chr_*.hdf5" \
        -ldsc "/disk/genetics/ukb/alextisyoung/haplotypes/simulated_pops_large/from_chr1_to_chr23_start0_endNone_run0_p0-0_ab_corr0-5_vb0-25_length2/ldscores/*[0-9].l2.ldscore.gz" \
        -l "${path2ldscreg}/multigen_sim/05012020/${runno}_pop.log" \
        --jkse \
        --jkse_blocksize 1000 \
        --jkse_cores 48 \
        -maf 1


    echo "Average Parental Effect"
    python ${path2ldscreg}/run_estimates.py "${dir}/chr_*.hdf5" \
        -ldsc "/disk/genetics/ukb/alextisyoung/haplotypes/simulated_pops_large/from_chr1_to_chr23_start0_endNone_run0_p0-0_ab_corr0-5_vb0-25_length2/ldscores/*[0-9].l2.ldscore.gz" \
        -l "${path2ldscreg}/multigen_sim/05012020/${runno}_avgparental.log" \
        -e "direct_plus_averageparental" \
        --jkse \
        --jkse_blocksize 1000 \
        --jkse_cores 48 \
        -maf 1

done