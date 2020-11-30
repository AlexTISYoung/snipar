sniparpath="/homes/nber/harij/gitrepos/SNIPar"

python "${sniparpath}/ldsc_reg/run_estimates.py" \
        "/disk/genetics/ukb/alextisyoung/haplotypes/relatives/EA/chr_*.hdf5" \
        -ldsc "/disk/genetics/ukb/alextisyoung/haplotypes/relatives/bedfiles/ldscores/*[0-9].l2.ldscore.gz" \
        -l "${sniparpath}/ldsc_reg/estimating_EA.log" 