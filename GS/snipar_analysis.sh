# VCF imputed from HRC: --
# Filter VCF based for MAF>1%, R2>0.99, bi-allelic 
# Convert to phased bgen
# Infer IBD with snipar
# Impute parental genotypes with snipar
# Analyse traits with snipar
# Compute PGI using EA4 weights with SBayesR
# /disk/genetics/sibling_consortium/GS20k/alextisyoung/grandpar
gpardir='/disk/genetics/sibling_consortium/GS20k/alextisyoung/grandpar'
hapdir=$gpardir'/haplotypes'
plink2='/homes/nber/alextisyoung/plink2'
plink='/homes/nber/alextisyoung/plink'
qctool='/disk/genetics/ukb/alextisyoung/qctool/build/release/qctool_v2.0.7'
king='/homes/nber/alextisyoung/king'
gcta64='/homes/nber/alextisyoung/gcta_1.93.2beta/gcta64'
conda activate snipar_env
cd $gpardir
### Filter VCF for phased haplotypes of SNPs with MAF>1%, Rsq>0.99, AvgCall>0.99, HWE<10^(-6), bi-alleleic
for i in {1..22}
do
$gpardir/vcftools-vcftools-581c231/bin/bin/vcftools --gzvcf /disk/genetics/sibling_consortium/GS20k/aokbay/imputed/HRC/vcf/chr$i.dose.vcf.gz \
     --snps $hapdir/chr_$i'_MAF_0.01_call_0.99_Rsq_0.99.txt' \
     --remove-indels --maf 0.01 --hwe 0.000001 --phased --recode \
     --stdout  --min-alleles 2 --max-alleles 2 | gzip -c > $hapdir/chr_$i.vcf.gz
done
### Convert VCF to phased BGEN file ###
for i in {1..22}
do
/disk/genetics/ukb/alextisyoung/qctool/build/release/qctool_v2.0.7 -g $hapdir/chr_$i.vcf.gz -og $hapdir/chr_$i.bgen -os $hapdir/chr_$i.sample
done
### Convert to bed
for i in {1..22}
do
$plink2 --vcf $hapdir/chr_$i.vcf.gz --make-bed --out $hapdir/bedfiles/chr_$i
echo $hapdir/bedfiles/chr_$i >> $hapdir/bedfiles/merge_list.txt
done
### Merge
$plink --merge-list $hapdir/bedfiles/merge_list.txt --make-bed --out $hapdir/bedfiles/autosome

### Infer relations with KING
$king -b $hapdir/bedfiles/autosome.bed --related --cpus 20 --prefix $gpardir/king
### Infer all pairwise IBD segments
$king -b $hapdir/bedfiles/autosome.bed --ibdseg --cpus 80 --prefix $gpardir/king

### Load snipar python virtualenv ###
source $gpardir/env/bin/activate
mkdir $gpardir/ibd $gpardir/imputed $gpardir/traits
ibd.py --bed $hapdir/bedfiles/chr_@ --king $gpardir/king.kin0 --agesex $gpardir/agesex.txt --ld_out --threads 20 --out $gpardir/ibd/chr_@
impute.py --ibd $gpardir/ibd/chr_@.ibd --bgen $hapdir/chr_@_haps --king $gpardir/king.kin0 --agesex $gpardir/agesex.txt --threads 40 --out $gpardir/imputed/chr_@ -c

### GWAS ###
for i in {1..5}
do
mkdir $gpardir/traits/$i
gwas.py $gpardir/processed_traits_noadj.txt --out $gpardir/traits/$i/chr_@ --bgen $hapdir/chr_@_haps --imp $gpardir/imputed/chr_@ --covar $gpardir/covariates.fam --phen_index $i --threads 40
done

### PGS ###
rm $gpardir/gctb_2.02_Linux/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_sparse_mldm_list.txt
for i in {1..22}
do
echo $gpardir/gctb_2.02_Linux/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_chr$i'_v3_50k.ldm.sparse' >> $gpardir/gctb_2.02_Linux/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_sparse_mldm_list.txt
done

### Depression sumstats ###
mkdir $gpardir/pgs/depression
wget 'https://datashare.ed.ac.uk/bitstream/handle/10283/3203/PGC_UKB_depression_genome-wide.txt?sequence=3&isAllowed=y' -P $gpardir/pgs/depression

### Ever-smoker sumstats ###
mkdir $gpardir/pgs/ever_smoke
wget 'https://conservancy.umn.edu/bitstream/handle/11299/241912/EUR_stratified.zip?sequence=20&isAllowed=y' -P $gpardir/pgs/ever_smoke

### UKB Height ###
mkdir $gpardir/pgs/height
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/variants.tsv.bgz -O variants.tsv.bgz -P $gpardir/pgs/height
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/50_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O 50_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -P $gpardir/pgs/height
### UKB BMI ###
mkdir $gpardir/pgs/bmi
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/21001_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O 21001_irnt.gwas.imputed_v3.both_sexes.tsv.bgz

# Compute EA4 PGS weights with gctb SBayesR
$gpardir/gctb_2.03beta_Linux/gctb --sbayes R \
     --mldm $gpardir/gctb_2.02_Linux/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_sparse_mldm_list.txt \
     --pi 0.95,0.02,0.02,0.01 \
     --gamma 0.0,0.01,0.1,1 \
     --gwas-summary $gpardir/pgs/depression/PGC_UKB_depression.ma \
     --chain-length 10000 \
     --burn-in 2000 \
     --out-freq 10 \
     --out $gpardir/pgs/depression/PGC_UKB_depression \
     --exclude-mhc \
     --unscale-genotype \
     --impute-n 

# Compute depression PGS weights with gctb SBayesR
$gpardir/gctb_2.03beta_Linux/gctb --sbayes R \
     --mldm $gpardir/gctb_2.02_Linux/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_sparse_mldm_list.txt \
     --pi 0.95,0.02,0.02,0.01 \
     --gamma 0.0,0.01,0.1,1 \
     --gwas-summary $gpardir/EA4_excl_UKBrel_STR_GS_2020_08_21.ma \
     --chain-length 10000 \
     --burn-in 2000 \
     --out-freq 10 \
     --out $gpardir/pgs/EA4_excl_UKBrel_STR_GS_2020_08_21_hm3 \
     --exclude-mhc \
     --unscale-genotype \
     --impute-n 

# Compute ever-smoke PGS weights with gctb SBayesR
$gpardir/gctb_2.03beta_Linux/gctb --sbayes R \
     --mldm $gpardir/gctb_2.02_Linux/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_sparse_mldm_list.txt \
     --pi 0.95,0.02,0.02,0.01 \
     --gamma 0.0,0.01,0.1,1 \
     --gwas-summary $gpardir/pgs/ever_smoke/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR.ma \
     --chain-length 10000 \
     --burn-in 2000 \
     --out-freq 10 \
     --out $gpardir/pgs/ever_smoke/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR \
     --exclude-mhc \
     --unscale-genotype \
     --p-value 0.4 \
     --robust --rsq 0.9

# Compute height PGS weights with gctb SBayesR
$gpardir/gctb_2.03beta_Linux/gctb --sbayes R \
     --mldm $gpardir/gctb_2.02_Linux/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_sparse_mldm_list.txt \
     --pi 0.95,0.02,0.02,0.01 \
     --gamma 0.0,0.01,0.1,1 \
     --gwas-summary $gpardir/pgs/height/height_UKB.ma \
     --chain-length 10000 \
     --burn-in 2000 \
     --out-freq 10 \
     --out $gpardir/pgs/height/height_UKB \
     --exclude-mhc \
     --unscale-genotype \
     --impute-n 

# Compute BMI weights with gctb SBayesR
$gpardir/gctb_2.03beta_Linux/gctb --sbayes R \
     --mldm $gpardir/gctb_2.02_Linux/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_sparse_mldm_list.txt \
     --pi 0.95,0.02,0.02,0.01 \
     --gamma 0.0,0.01,0.1,1 \
     --gwas-summary $gpardir/pgs/bmi/BMI_UKB.ma \
     --chain-length 10000 \
     --burn-in 2000 \
     --out-freq 10 \
     --out $gpardir/pgs/bmi/BMI_UKB \
     --exclude-mhc \
     --unscale-genotype \
     --impute-n 

## Compute PGS
Rscript sbayesr_to_snipar.R
pgs.py $gpardir/pgs/EA4_hm3 --weights $gpardir/pgs/EA4_excl_UKBrel_STR_GS_2020_08_21_hm3.txt --bgen $hapdir/chr_@_haps --imp $gpardir/imputed/chr_@ --beta_col beta --grandpar
#pgs.py $gpardir/pgs/EA4_hm3_sib --weights $gpardir/pgs/EA4_excl_UKBrel_STR_GS_2020_08_21_hm3.txt --bgen $hapdir/chr_@_haps --imp $gpardir/imputed/chr_@ --beta_col beta --fit_sib
pgs.py $gpardir/pgs/depression/PGC_UKB_depression --weights $gpardir/pgs/depression/PGC_UKB_depression.txt --bgen $hapdir/chr_@_haps --imp $gpardir/imputed/chr_@ --beta_col beta --grandpar --threads 20
pgs.py $gpardir/pgs/height/height_UKB --weights $gpardir/pgs/height/height_UKB.txt --bgen $hapdir/chr_@_haps --imp $gpardir/imputed/chr_@ --beta_col beta --grandpar --threads 30
pgs.py $gpardir/pgs/bmi/BMI_UKB --weights $gpardir/pgs/bmi/BMI_UKB.txt --bgen $hapdir/chr_@_haps --imp $gpardir/imputed/chr_@ --beta_col beta --grandpar --threads 30
pgs.py $gpardir/pgs/ever_smoke/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR --weights $gpardir/pgs/ever_smoke/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR.txt --bgen $hapdir/chr_@_haps --imp $gpardir/imputed/chr_@ --beta_col beta --grandpar --threads 30


# Estimated correlation between maternal and paternal EA4 PGSs: 0.1481
for pgs in EA4_hm3
do
for i in {1..16}
do
pgs.py $gpardir/pgs/$pgs/$i --pgs $gpardir/pgs/$pgs'.pgs.txt' --phenofile $gpardir/processed_traits_noadj_noukb.txt --covar $gpardir/covariates.fam  --gen_models 1-3 --phen_index $i --scale_pgs --scale_phen --sparse_thres 0.05 --ibdrel_path $gpardir/king
#pgs.py $gpardir/pgs/$pgs/$i'_sib' --pgs $gpardir/pgs/$pgs'_sib.pgs.txt' --phenofile $gpardir/processed_traits_noadj_noukb.txt --covar $gpardir/covariates.fam  --gen_models 2 --phen_index $i --scale_pgs --scale_phen --sparse_thres 0.05 --ibdrel_path $gpardir/king --fit_sib
done
done

# Depression: Estimated correlation between maternal and paternal PGSs: 0.035 S.E.=0.0112
for i in {1..16}
do
pgs.py $gpardir/pgs/depression/$i --pgs $gpardir/pgs/depression/PGC_UKB_depression.pgs.txt --phenofile $gpardir/processed_traits_noadj_noukb.txt --covar $gpardir/covariates.fam  --gen_models 1-3 --phen_index $i --scale_pgs --scale_phen --sparse_thres 0.05 --ibdrel_path $gpardir/king --threads 20
done

# Height: Estimated correlation between maternal and paternal PGSs: 0.0729 S.E.=0.0108
for i in {1..16}
do
pgs.py $gpardir/pgs/height/$i --pgs $gpardir/pgs/height/height_UKB.pgs.txt --phenofile $gpardir/processed_traits_noadj_noukb.txt --covar $gpardir/covariates.fam  --gen_models 1-3 --phen_index $i --scale_pgs --scale_phen --sparse_thres 0.05 --ibdrel_path $gpardir/king --threads 20
done

# BMI: Estimated correlation between maternal and paternal PGSs: 0.0488 S.E.=0.011
for i in {1..16}
do
pgs.py $gpardir/pgs/bmi/$i --pgs $gpardir/pgs/bmi/BMI_UKB.pgs.txt --phenofile $gpardir/processed_traits_noadj_noukb.txt --covar $gpardir/covariates.fam  --gen_models 1-3 --phen_index $i --scale_pgs --scale_phen --sparse_thres 0.05 --ibdrel_path $gpardir/king --threads 20
done

# Ever-smoke: Estimated correlation between maternal and paternal PGSs: 0.0395 S.E.=0.0111
for i in {1..16}
do
pgs.py $gpardir/pgs/ever_smoke/$i --pgs $gpardir/pgs/ever_smoke/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR.pgs.txt --phenofile $gpardir/processed_traits_noadj_noukb.txt --covar $gpardir/covariates.fam  --gen_models 1-3 --phen_index $i --scale_pgs --scale_phen --sparse_thres 0.05 --ibdrel_path $gpardir/king --threads 20
done

# Logistic LMM
Rscript logistic_lmm.R 
