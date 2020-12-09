path2ldscreg="/homes/nber/harij/gitrepos/SNIPar/ldsc_reg"


# for dir in /disk/genetics/ukb/alextisyoung/GS20k_sumstats/traits/[0-9]*/
# do
#     traitno=$(echo ${dir:54:2} | tr '/' ' ' )
#     echo $traitno 
#     python ${path2ldscreg}/run_estimates.py "$dir/chr_*.hdf5" \
#         -ldsc "/disk/genetics/ukb/alextisyoung/GS20k_sumstats/ldscores/*[0-9].l2.ldscore.gz" \
#         -l "${path2ldscreg}/generation_scotland/$traitno.log" \
#         --jkse \
#         --jkse_blocksize 1000 \
#         --jkse_cores 24

# done

for dir in /disk/genetics/ukb/alextisyoung/GS20k_sumstats/traits/[0-9]*/
do
    traitno=$(echo ${dir:54:2} | tr '/' ' ' )
    echo $traitno 
    python ${path2ldscreg}/run_estimates.py "$dir/chr_*.hdf5" \
        -ldsc "/disk/genetics/ukb/alextisyoung/GS20k_sumstats/ldscores/*[0-9].l2.ldscore.gz" \
        -l "${path2ldscreg}/generation_scotland/${traitno}_nobound.log" \
        --no-rbound \
        --jkse \
        --jkse_blocksize 1000 \
        --jkse_cores 24

done