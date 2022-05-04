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

### Filter VCF for phased haplotypes of SNPs with MAF>1%, Rsq>0.99, AvgCall>0.99, HWE<10^(-6), bi-alleleic
for i in {1..22}
do
$gpardir/vcftools-vcftools-581c231/bin/bin/vcftools --gzvcf /disk/genetics/sibling_consortium/GS20k/aokbay/imputed/HRC/vcf/chr$i.dose.vcf.gz --snps $hapdir/chr_$i'_MAF_0.01_call_0.99_Rsq_0.99.txt' --remove-indels --maf 0.01 --hwe 0.000001 --phased --recode --stdout | gzip -c > $hapdir/chr_$i.vcf.gz
done
### Convert VCF to phased BGEN file ###
for i in {1..22}
do
$qctool -g $hapdir/chr_$i.vcf.gz -og $hapdir/chr_$i.bgen -os $hapdir/chr_$i.sample
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
$king -b $hapdir/bedfiles/autosome.bed --related --degree 1 --cpus 20 --prefix $gpardir/king

### Load snipar python virtualenv ###
source $gpardir/env/bin/activate
mkdir $gpardir/ibd $gpardir/imputed $gpardir/traits
ibd.py --bed $hapdir/bedfiles/chr_@ --king $gpardir/king.kin0 --agesex $gpardir/agesex.txt --ld_out --threads 20 --out $gpardir/ibd/chr_@
impute.py --ibd $gpardir/ibd/chr_@.ibd --bgen $hapdir/chr_@_haps --king $gpardir/king.kin0 --agesex $gpardir/agesex.txt --threads 40 --out $gpardir/imputed/chr_@
for i in {14..16}
do
mkdir traits/$i
gwas.py ../../processed_traits_noadj.txt --out traits/$i/chr_@ --bgen ../genotypes/haplotypes/chr_@_haps --imp imputed/chr_@ --covar ../../covariates.fam --phen_index $i --threads 20
done