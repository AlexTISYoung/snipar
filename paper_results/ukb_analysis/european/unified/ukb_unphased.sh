gwas.py \
    /disk/genetics/ukb/alextisyoung/phenotypes/processed_traits_noadj.txt \
    --imp /disk/genetics/ukb/jguan/ukb_analysis/output/parent_imputed/robust/chr@ \
    --bed /disk/genetics/ukb/jguan/ukb_analysis/output/plink/v4/ukb_imp_chr@_filtered \
    --covar /disk/genetics/ukb/alextisyoung/phenotypes/covariates.txt \
    --impute_unrel \
    --chr_range 22 \
    --cpus 20 \
    --threads 1 \
    --cpus 30 \
    --phen_index 8 \
    --out /disk/genetics/ukb/alextisyoung/junming/unphased_gwas/8/chr_@ \
    --ibdrel_path /disk/genetics/ukb/jguan/preprocessed/king/ukb_v2_3rd

--ibdrel_path /disk/genetics/ukb/jguan/preprocessed/king/ukb_v2_3rd
