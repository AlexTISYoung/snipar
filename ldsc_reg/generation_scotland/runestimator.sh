path2ldscreg="/homes/nber/harij/gitrepos/SNIPar/ldsc_reg"


# run for direct + population effects, bounds on r
# for dir in /disk/genetics/ukb/alextisyoung/GS20k_sumstats/traits/[0-9]*/
# do
#     traitno=$(echo ${dir:54:2} | tr '/' ' ' )
#     echo $traitno 
#     python ${path2ldscreg}/run_estimates.py "$dir/chr_*.hdf5" \
#         -ldsc "/disk/genetics/ukb/alextisyoung/GS20k_sumstats/ldscores/*[0-9].l2.ldscore.gz" \
#         -l "${path2ldscreg}/generation_scotland/$traitno.log" \
#         --jkse \
            # -maf 0 \
#         --jkse_blocksize 1000 \
#         --jkse_cores 24

# done

# run for direct + population effects, no bounds on r
# for dir in /disk/genetics/ukb/alextisyoung/GS20k_sumstats/traits/[0-9]*/
# do
#     traitno=$(echo ${dir:54:2} | tr '/' ' ' )
#     echo $traitno 
#     python ${path2ldscreg}/run_estimates.py "$dir/chr_*.hdf5" \
#         -ldsc "/disk/genetics/ukb/alextisyoung/GS20k_sumstats/ldscores/*[0-9].l2.ldscore.gz" \
#         -l "${path2ldscreg}/generation_scotland/${traitno}_nobound.log" \
#         --no-rbound \
#           -maf 0 \
#         --jkse \
#         --jkse_blocksize 1000 \
#         --jkse_cores 24

# done

# run for direct + average parental effects, no bounds on r
for dir in /disk/genetics/ukb/alextisyoung/GS20k_sumstats/traits/[0-9]*/
do
    traitno=$(echo ${dir:54:2} | tr '/' ' ' )
    echo $traitno 
    python ${path2ldscreg}/run_estimates.py "$dir/chr_*.hdf5" \
        -ldsc "/disk/genetics/ukb/alextisyoung/GS20k_sumstats/ldscores/*[0-9].l2.ldscore.gz" \
        -l "${path2ldscreg}/generation_scotland/${traitno}_nobound_avgparental.log" \
        -e "direct_plus_averageparental" \
        -maf 0 \
        --no-rbound \
        --jkse \
        --jkse_blocksize 1000 \
        --jkse_cores 24

done