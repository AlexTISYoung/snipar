path2ldscreg="/homes/nber/harij/gitrepos/SNIPar/ldsc_reg"


# for dir in /disk/genetics/ukb/alextisyoung/haplotypes/relatives/traits/[0-9]*/
# do
#     traitno=$(echo ${dir:60:2} | tr '/' ' ' )
#     echo $traitno 
#     python ${path2ldscreg}/run_estimates.py "$dir/chr_*.hdf5" \
#         -ldsc "/disk/genetics/ukb/alextisyoung/haplotypes/relatives/bedfiles/ldscores/*[0-9].l2.ldscore.gz" \
#         -l "${path2ldscreg}/ukb/$traitno.log" \
#         -maf 0 \
#         --jkse \
#         --jkse_blocksize 1000 \
#         --jkse_cores 24

# done


# running without a bound on r
for dir in /disk/genetics/ukb/alextisyoung/haplotypes/relatives/traits/[0-9]*/
do
    traitno=$(echo ${dir:60:2} | tr '/' ' ' )
    echo $traitno 
    echo "Estimating for population effect, no r bound."
    python ${path2ldscreg}/run_estimates.py "$dir/chr_*.hdf5" \
        -ldsc "/disk/genetics/ukb/alextisyoung/haplotypes/relatives/bedfiles/ldscores/*[0-9].l2.ldscore.gz" \
        -l "${path2ldscreg}/ukb/${traitno}_norbound.log" \
        -maf 0 \
        --no-rbound \
        --jkse \
        --jkse_blocksize 1000 \
        --jkse_cores 24

done



# running without a bound on r, summ stats are direct + average parental effects
# for dir in /disk/genetics/ukb/alextisyoung/haplotypes/relatives/traits/[0-9]*/
# do
#     traitno=$(echo ${dir:60:2} | tr '/' ' ' )
#     echo $traitno 
#     python ${path2ldscreg}/run_estimates.py "$dir/chr_*.hdf5" \
#         -ldsc "/disk/genetics/ukb/alextisyoung/haplotypes/relatives/bedfiles/ldscores/*[0-9].l2.ldscore.gz" \
#         -l "${path2ldscreg}/ukb/${traitno}_norbound_avgparental.log" \
#         --no-rbound \
#         -maf 0 \
#         -e "direct_plus_averageparental" \
#         --jkse \
#         --jkse_blocksize 1000 \
#         --jkse_cores 24

# done