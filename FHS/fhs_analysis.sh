### Grandparental analysis in FHS #33
gpardir='/var/genetics/data/fhs/private/v32/raw/extracted/v20'
hapdir=$gpardir'/bgen'
plink2='/homes/nber/alextisyoung/plink2'
plink='/homes/nber/alextisyoung/plink'
qctool='/disk/genetics/ukb/alextisyoung/qctool/build/release/qctool_v2.0.7'
king='/homes/nber/alextisyoung/king'
gcta64='/homes/nber/alextisyoung/gcta_1.93.2beta/gcta64'
gsdir='/disk/genetics/sibling_consortium/GS20k/alextisyoung/grandpar'
pgsdir='/var/genetics/data/fhs/private/v32/raw/pgs'

### Infer all pairwise IBD segments
$king -b $hapdir/bedfiles/autosome.bed --ibdseg --cpus 80 --prefix $gpardir/king

### Load snipar python virtualenv ###
source $gsdir/env/bin/activate

### PGS ###
# Compute PGS weights with gctb SBayesR
$gsdir/gctb_2.03beta_Linux/gctb --sbayes R \
     --mldm $gsdir/gctb_2.02_Linux/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_sparse_mldm_list.txt \
     --pi 0.95,0.02,0.02,0.01 \
     --gamma 0.0,0.01,0.1,1 \
     --gwas-summary $pgsdir/EA4_excl_UKBrel_STR_GS_2020_08_21.ma \
     --chain-length 10000 \
     --burn-in 2000 \
     --out-freq 10 \
     --out $pgsdir/EA4_excl_UKBrel_STR_GS_2020_08_21_hm3 \
     --exclude-mhc \
     --unscale-genotype \
     --impute-n 
     
## Compute PGS
pgs.py $pgsdir/EA4_excl_UKBrel_STR_GS_2020_08_21_hm3 --weights $pgsdir/EA4_excl_UKBrel_STR_GS_2020_08_21_hm3.snpRes --bed $gpardir/bedfiles/chr@_filtered --imp $gpardir/parent_imputed/chr_@ --beta_col A1Effect --SNP Name --grandpar
# Estimated correlation between maternal and paternal PGSs: 0.1481
for i in {1..16}
do
pgs.py $gpardir/pgs/results/$i --pgs $gpardir/pgs/EA4_hm3.pgs.txt --phenofile $gpardir/processed_traits_noadj.txt --covar $gpardir/covariates.fam  --gen_models 1-3 --phen_index $i --scale_pgs --scale_phen --bpg --sparse_thres 0.025 --ibdrel_path $gpardir/king
done
