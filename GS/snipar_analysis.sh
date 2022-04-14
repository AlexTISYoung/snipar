# VCF imputed from HRC: --
# Filter VCF based for MAF>1%, R2>0.99, bi-allelic 
# Convert to phased bgen
# Infer IBD with snipar
# Impute parental genotypes with snipar
# Analyse traits with snipar
# Compute PGI using EA4 weights with SBayesR
# /disk/genetics/sibling_consortium/GS20k/alextisyoung/grandpar

### Filter VCF for phased haplotypes of SNPs with MAF>1%, Rsq>0.99, AvgCall>0.99, HWE<10^(-6), bi-alleleic
for i in {1..21}
do
/disk/genetics/sibling_consortium/GS20k/alextisyoung/grandpar/vcftools-vcftools-581c231/bin/bin/vcftools --gzvcf /disk/genetics/sibling_consortium/GS20k/aokbay/imputed/HRC/vcf/chr$i.dose.vcf.gz --snps /disk/genetics/sibling_consortium/GS20k/alextisyoung/grandpar/haplotypes/chr_$i'_MAF_0.01_call_0.99_Rsq_0.99.txt' --remove-indels --maf 0.01 --hwe 0.000001 --phased --recode --stdout | gzip -c > haplotypes/chr_$i.vcf.gz
done
### Convert VCF to phased BGEN file ###
/disk/genetics/ukb/alextisyoung/qctool/build/release/qctool_v2.0.7 

# /disk/genetics/sibling_consortium/GS20k/alextisyoung/HM3/tidy
source env/bin/activate
ibd.py --bgen ../genotypes/haplotypes/chr_@_haps --pedigree ../../pedigree.txt --out ibd/chr_@ --ld_out --threads 20
impute.py --ibd ibd/chr_@.ibd --bgen ../genotypes/haplotypes/chr_@_haps --pedigree ../../pedigree.txt --threads 20 --out imputed/chr_@
for i in {14..16}
do
mkdir traits/$i
gwas.py ../../processed_traits_noadj.txt --out traits/$i/chr_@ --bgen ../genotypes/haplotypes/chr_@_haps --imp imputed/chr_@ --covar ../../covariates.fam --phen_index $i --threads 20
done