path2ldscreg="/homes/nber/harij/gitrepos/SNIPar/ldsc_reg"

for dir in /disk/genetics/ukb/alextisyoung/haplotypes/multi_sims/from_chr1_to_chr23_start0_end50_run*
do

    runno=${dir:89:2}

    echo "Run Number: $runno"
    python ${path2ldscreg}/run_estimates.py "${dir}/fgwas/chr_*.hdf5" \
        -ldsc "${dir}/ldscores/*[0-9].l2.ldscore.gz" \
        -l "${path2ldscreg}/multigen_sim/multisims/${runno}.log" \
        --jkse \
        --jkse_blocksize 20 \
        --jkse_cores 24 \
        -maf 1

done