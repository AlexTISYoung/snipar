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
# Compute PGS weights with gctb SBayesR
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
     
## Compute PGS
Rscript sbayesr_to_snipar.R
pgs.py $gpardir/pgs/EA4_hm3 --weights $gpardir/pgs/EA4_excl_UKBrel_STR_GS_2020_08_21_hm3.txt --bgen $hapdir/chr_@_haps --imp $gpardir/imputed/chr_@ --beta_col beta --grandpar
# Estimated correlation between maternal and paternal PGSs: 0.1481
for i in {1..16}
do
pgs.py $gpardir/pgs/results/$i --pgs $gpardir/pgs/EA4_hm3.pgs.txt --phenofile $gpardir/processed_traits_noadj.txt --covar $gpardir/covariates.fam  --gen_models 1-3 --phen_index $i --scale_pgs --scale_phen --bpg --ibdrel_path $gpardir/king
done
