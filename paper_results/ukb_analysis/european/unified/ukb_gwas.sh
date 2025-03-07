conda activate junming

for phenotype in Non_HDL cigarettes.per.day HDL SBP DBP self.rated.health Cognitive.ability Neuroticism \
                 AAFB household.income subjective.well.being drinks.per.week height menarche FEV1 \
                 depressive.symptoms ever.smoked myopia morning.person NC
do
    if [ "$phenotype" == "AAFB" ] || [ "$phenotype" == "menarche" ]
    then
        covariates="/disk/genetics/ukb/jguan/ukb_analysis/output/covariates_nosex.txt"
    else
        covariates="/disk/genetics/ukb/jguan/ukb_analysis/output/covariates_allancestry.txt"
    fi
    
    gwas.py \
        /disk/genetics/ukb/alextisyoung/phenotypes/processed_traits_noadj.txt \
        --imp /disk/genetics/ukb/jguan/ukb_analysis/output/parent_imputed/v2/chr_@ \
        --bgen /disk/genetics/ukb/jguan/ukb_analysis/output/bgen/ukb_hap_chr@_v2_filtered \
        --covar "$covariates" \
        --cpus 20 \
        --phen "$phenotype" \
        --out "/disk/genetics/ukb/alextisyoung/junming/unified_ukb_gwas/${phenotype}_phased_chr_@" \
        --ibdrel_path /disk/genetics/ukb/jguan/preprocessed/king/ukb_v2_3rd \
        --impute_unrel \
        --vc_out "/disk/genetics/ukb/alextisyoung/junming/unified_ukb_gwas/${phenotype}_phased_chr_@_vc"
done

gwas.py \
    /disk/genetics/ukb/alextisyoung/phenotypes/processed_traits_noadj.txt \
    --imp /disk/genetics/ukb/jguan/ukb_analysis/output/parent_imputed/v2/chr_@ \
    --bgen /disk/genetics/ukb/jguan/ukb_analysis/output/bgen/ukb_hap_chr@_v2_filtered \
    --covar /disk/genetics/ukb/alextisyoung/phenotypes/covariates.txt \
    --cpus 20 \
    --phen_index 8 \
    --chr_range 22 \
    --out "/disk/genetics/ukb/alextisyoung/junming/unified_ukb_gwas/BMI_phased_chr_@" \
    --ibdrel_path /disk/genetics/ukb/jguan/preprocessed/king/ukb_v2_3rd \
    --impute_unrel \
    --vc_out "/disk/genetics/ukb/alextisyoung/junming/unified_ukb_gwas/BMI_phased_chr_@_vc"