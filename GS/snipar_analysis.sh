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
$gpardir/vcftools-vcftools-581c231/bin/bin/vcftools --gzvcf /disk/genetics/sibling_consortium/GS20k/aokbay/imputed/HRC/vcf/chr$i.dose.vcf.gz --snps $hapdir/chr_$i'_MAF_0.01_call_0.99_Rsq_0.99.txt' --remove-indels --maf 0.01 --hwe 0.000001 --phased --recode --stdout | gzip -c > $hapdir/chr_$i.vcf.gz
done
### Convert VCF to phased BGEN file ###
for i in {1..22}
do
idone
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
impute.py --ibd $gpardir/ibd/chr_@.ibd --bgen $hapdir/chr_@_haps --king $gpardir/king.kin0 --agesex $gpardir/agesex.txt --threads 40 --out $gpardir/imputed/chr_@ -c

### GWAS ###
for i in {15..1}
do
mkdir $gpardir/traits/$i
gwas.py $gpardir/processed_traits_noadj.txt --out $gpardir/traits/$i/chr_@ --bgen $hapdir/chr_@_haps --imp $gpardir/imputed/chr_@ --covar $gpardir/covariates.fam --phen_index $i --threads 40
done

### PGS ###
# Compute PGS
pgs.py $gpardir/pgs/GS_EA_13_weights_LDpred_p1 --weights $gpardir/pgs/GS_EA_13_weights_LDpred_p1.0000e+00_matched.txt --bgen $hapdir/chr_@_haps --imp $gpardir/imputed/chr_@ --compute_controls
# Compute grandparental PGS
R3script $gpardir/impute_gpar_PGS_GS.R

### Compute GRM ###
$gcta64 --make-grm-bin --bfile /disk/genetics/sibling_consortium/GS20k/GS20k_TopStrand --maf 0.05 --thread-num 40 --out $gpardir/grms/R
python $gpardir/make_grms.py

### Compute variance components ###
mkdir $gpardir/grms/varcomps
$gcta64 --mgrm $gpardir/grms/mgrm.txt --reml --reml-no-lrt --pheno $gpardir/processed_traits_noadj.txt --mpheno 16 --qcovar $gpardir/pgs/GS_EA_13_weights_LDpred_p1.pgs.with_covariates.txt --out $gpardir/grms/varcomps/16 --thread-num 20

### Estimate grandparental PGS model ###

