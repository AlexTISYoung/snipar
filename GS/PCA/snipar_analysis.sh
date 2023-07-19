# /disk/genetics/sibling_consortium/GS20k/alextisyoung/HM3/tidy
source env/bin/activate
ibd.py --bgen ../genotypes/haplotypes/chr_@_haps --pedigree ../../pedigree.txt --out ibd/chr_@ --ld_out --threads 20
impute.py --ibd ibd/chr_@.ibd --bgen ../genotypes/haplotypes/chr_@_haps --pedigree ../../pedigree.txt --threads 20 --out imputed/chr_@
for i in {14..16}
do
mkdir traits/$i
gwas.py ../../processed_traits.fam --out traits/$i/chr_@ --bgen ../genotypes/haplotypes/chr_@_haps --imp imputed/chr_@ --covar ../../covariates.fam --phen_index $i --threads 20
done