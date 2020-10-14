ldsc_repo="/disk/genetics2/pub/repo/ssgac/ldsc_mod" 
snipar_repo="/homes/nber/harij/gitrepos/SNIPar"
genotypes="/disk/genetics/ukb/alextisyoung/genotypes"

python ${ldsc_repo}/ldsc.py --bfile ${genotypes}/chr_22 --l2 --ld-wind-cm 1 --out ${snipar_repo}/ldsc_reg/inferz/ldsc_reg/22 --yes-really