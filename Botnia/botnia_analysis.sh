# VCF imputed from HRC: --
# Filter VCF based for MAF>1%, R2>0.99, bi-allelic 
# Convert to phased bgen
# Infer IBD with snipar
# Impute parental genotypes with snipar
# Analyse traits with snipar
# Compute PGI using EA4 weights with SBayesR

# Define directories
gpardir='/ludc/Active_Projects/BOTNIA_AYoung_analysis/Private/'
hapdir=$gpardir'/haplotypes'
vcf_dir='/ludc/Raw_Data_Archive/Chip_Genotyping/BOTEX_GWAS/Imputation/IMPUTE3_Michigan_imputation/'
# Define tools
vcftools='/ludc/Tools/Software/VCFtools/0.1.16-13/bin/vcftools'
plink2='/ludc/Tools/Software/Plink/v2.00a2.3LM/plink2'
plink='/ludc/Tools/Software/Plink/v1.90p/plink'
qctool='/ludc/Tools/Software/qctool/2.0.6/qctool'
king=$gpardir'king/king'
#gcta64='/homes/nber/alextisyoung/gcta_1.93.2beta/gcta64'

### Filter VCF for phased haplotypes of SNPs with MAF>1%, Rsq>0.99, AvgCall>0.99, HWE<10^(-6), bi-alleleic
# Read info files
Rscript read_vcf_info.R
# filter vcfs
for i in {6..22}
do
$vcftools --gzvcf $vcf_dir/chr$i.dose.vcf.gz --snps $hapdir/chr_$i'_MAF_0.01_call_0.99_Rsq_0.99.txt' \
    --remove-indels --maf 0.01 --hwe 0.000001 --phased --recode --stdout | \
    gzip -c > $hapdir/chr_$i.vcf.gz&
done
### Convert to bed
mkdir $hapdir/bedfiles
for i in {1..22}
do
$plink2 --vcf $hapdir/chr_$i.vcf.gz --make-bed --out $hapdir/bedfiles/chr_$i&
echo $hapdir/bedfiles/chr_$i >> $hapdir/bedfiles/merge_list.txt
done
### Merge
# Remove multi-allelic variants
for i in {1..22}
do
$plink --bfile $hapdir/bedfiles/chr_$i --exclude $hapdir/bedfiles/autosome-merge.missnp --make-bed --out $hapdir/bedfiles/chr_$i
done
# 
$plink --merge-list $hapdir/bedfiles/merge_list.txt --make-bed --out $hapdir/bedfiles/autosome --exclude $hapdir/bedfiles/autosome-merge.missnp
### Infer relations with KING
$king -b $hapdir/bedfiles/autosome.bed --related --cpus 20 --prefix $gpardir/king
# 
#               MZ      PO      FS      2nd
#  =====================================================
#  Inference    28     4778    6805      5
### Find unrelated subsample
$king -b $hapdir/bedfiles/pruned/autosome.bed --unrelated --cpus 20 --prefix $gpardir/unrelated


### Perform PCA to identify European samples
mkdir $hapdir/bedfiles/pca
mkdir $hapdir/bedfiles/pca/data
# Get population labels
wget https://personal.broadinstitute.org/hhuang//public//GINGER/pop/igsr_samples.tsv > $hapdir/bedfiles/pca/data/
# perform PCA
Rscript botnia_PCA.R

## Filter non-EUR from .bed ##
for i in {1..22}
do 
$plink --bfile $hapdir/bedfiles/chr_$i --remove $hapdir/bedfiles/pca/non_european_samples.txt --make-bed --out $hapdir/bedfiles/chr_$i
done 

### Convert VCF to phased BGEN file ###
for i in {1..22}
do
$qctool -g $hapdir/chr_$i.vcf.gz -og $hapdir/chr_$i.bgen -os $hapdir/chr_$i.sample -excl-snpids $hapdir/bedfiles/autosome-merge.missnp -excl-samples  $hapdir/bedfiles/pca/non_european_samples.txt > chr_$i.log&
done
### Recode bgen to use proper sample IDs and remove multi-allelic/duplicated SNPs ###
for i in {1..22}
do
$qctool -g $hapdir/chr_$i.bgen -s $hapdir/chr_$i.sample -og $hapdir/chr_$i'_haps.bgen' \
    -os $hapdir/chr_$i'_haps.sample' -excl-snpids $hapdir/bedfiles/autosome-merge.missnp \
    -excl-samples  $hapdir/bedfiles/pca/non_european_samples.txt > chr_$i.log&
done

### Load snipar python virtualenv ###
## install snipar
python3.9 -m venv $gpardir/env
source $gpardir/env/bin/activate
pip install snipar
## make directories
mkdir $gpardir/ibd $gpardir/imputed $gpardir/traits
## Infer IBD
# LD prune variants to reduce computation time
for i in {1..22}
do
$plink --bfile $hapdir/bedfiles/chr_$i --maf 0.05 --indep 10 10 10 --out $hapdir/bedfiles/pruned/chr_$i
$plink --bfile chr_$i --extract $hapdir/bedfiles/pruned/chr_$i.prune.in --make-bed --out $hapdir/bedfiles/pruned/chr_$i
done
#ibd.py --bed $hapdir/bedfiles/chr_@ --king $gpardir/king.kin0 --agesex $gpardir/phenotypes/agesex.txt --ld_out --threads 20 --out $gpardir/ibd/chr_@ 
ibd.py --bed $hapdir/bedfiles/pruned/chr_@ --king $gpardir/king.kin0 --agesex $gpardir/phenotypes/agesex.txt --ld_out --threads 20 --out $gpardir/ibd/chr_@_pruned 
#Estimated mean genotyping error probability: 0.000151

### Infer pairwise IBD segments
$plink --merge-list $hapdir/bedfiles/pruned/merge_list.txt --make-bed --out $hapdir/bedfiles/pruned/autosome
$king -b $hapdir/bedfiles/pruned/autosome.bed --ibdseg --cpus 80 --prefix $gpardir/king
### 

## Impute
impute.py --ibd $gpardir/ibd/chr_@_pruned.ibd --bgen $hapdir/chr_@_haps --king $gpardir/king.kin0 --agesex $gpardir/phenotypes/agesex.txt --threads 40 --out $gpardir/imputed/chr_@ -c --chr_range 6

### GWAS ###
for i in {1..5}
do
mkdir $gpardir/traits/$i
gwas.py $gpardir/processed_traits_noadj.txt --out $gpardir/traits/$i/chr_@ --bgen $hapdir/chr_@_haps --imp $gpardir/imputed/chr_@ --covar $gpardir/covariates.fam --phen_index $i --threads 40
done

mkdir $gpardir/traits/EA
gwas.py $gpardir/phenotypes/EA.txt --out $gpardir/traits/EA/chr_@ --bgen $hapdir/chr_@_haps --imp $gpardir/imputed/chr_@ --threads 40



### PGS ###
rm $gpardir/gctb_2.03beta_Linux/ukb_50k_bigset_2.8M/ukb50k_2.8M_shrunk_sparse.mldmlist
for i in {1..22}
do
echo $gpardir/gctb_2.03beta_Linux/ukb_50k_bigset_2.8M/ukb50k_shrunk_chr$i'_mafpt01.ldm.sparse' >> $gpardir/gctb_2.03beta_Linux/ukb_50k_bigset_2.8M/ukb50k_2.8M_shrunk_sparse.mldmlist
done
# Compute PGS weights with gctb SBayesR
$gpardir/gctb_2.03beta_Linux/gctb --sbayes R \
     --mldm $gpardir/gctb_2.03beta_Linux/ukb_50k_bigset_2.8M/ukb50k_2.8M_shrunk_sparse.mldmlist \
     --pi 0.95,0.02,0.02,0.01 \
     --gamma 0.0,0.01,0.1,1 \
     --gwas-summary $gpardir/pgs/EA4.ma \
     --chain-length 10000 \
     --burn-in 2000 \
     --out-freq 10 \
     --out $gpardir/pgs/EA4_2.8m \
     --exclude-mhc \
     --unscale-genotype
     
# Compute PGS
Rscript sbayesr_to_snipar.R
pgs.py $gpardir/pgs/EA4_2.8m --weights $gpardir/pgs/EA4_2.8m.txt --bgen $hapdir/chr_@_haps --imp $gpardir/imputed/chr_@ --beta_col beta --grandpar
# Estimated correlation between maternal and paternal PGSs: 0.0732

### Estimate grandparental PGS model ###
mkdir $gpardir/pgs/results
for i in {1..8}
do
pgs.py $gpardir/pgs/results/$i --pgs $gpardir/pgs/EA4_2.8m.pgs.txt --phenofile $gpardir/phenotypes/processed_traits_noadj.txt --covar $gpardir/phenotypes/covariates.txt  --gen_models 1-3 --phen_index $i --scale_pgs --scale_phen --sparse_thres 0.025 --ibdrel_path $gpardir/king
done

### Perform family based GWAS ###
mkdir $gpardir/phenotypes/gwas
for i in {1..8}
do
mkdir $gpardir/phenotypes/gwas/$i/
gwas.py $gpardir/phenotypes/processed_traits_noadj.txt --out $gpardir/phenotypes/gwas/$i/chr_@ --bgen $hapdir/chr_@_haps --imp $gpardir/imputed/chr_@ --covar $gpardir/phenotypes/covariates.txt --phen_index $i 
done