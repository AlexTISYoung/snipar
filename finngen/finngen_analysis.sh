
hm3_dir=/finngen/library-red/finngen_R10/genotype_plink_1.0/data/
pca_outlier=/finngen/library-red/finngen_R10/pca_1.0/data/finngen_R10_rejected.txt
kin0='/finngen/library-red/finngen_R10/kinship_2.0/data/finngen_R10.kin0'
kin0_filtered='~/king_1st_degree_no_mz_no_pca_outlier.kin0'
ids_filtered='~/ids_1st_degree_no_mz_no_pca_outlier.txt'
agesex_raw='/finngen/library-red/finngen_R10/phenotype_1.0/data/finngen_R10_baseline_1.0.txt.gz'
agesex='~/agesex.txt'
outdir='/finngen/red/alex_young/'
vcf_dir='/finngen/library-red/finngen_R10/genotype_1.0/data/'

########### PGS ############
source snipar_env/bin/activate

##### Externalizing PGS #######
for (i in 1:22){
    bim=read.table(paste('~/haplotypes/chr_',i,'.bim',sep=''))
    chr_i = d[d$chr==i,]
    p = rbind(p,data.frame(CHR=i,
                SNP=bim[match(chr_i$pos_hg38,bim[,4]),2],
                            A1=chr_i$a1,A2=chr_i$a2,
                            beta=chr_i$sbayesr_SBayesR))
    print(i)
}

p=p[!is.na(p$SNP),]
write.table(p,'~/pgs/externalizing/externalizing_sbayesr.txt',quote=F,row.names=F)


##### PGI Catalog ######
gmap = read.table('~/pgs/HM3_SNPs_ChrPosID_hg19tohg38.map',header=T)
#d = read.table('PGI_catalog/',header=T)
d$pos = gmap[match(d$sid,gmap$ChrPosID_hg19),1]
d = d[!is.na(d$pos),]
d$pos = sapply(d$pos,function(x) strsplit(x,':')[[1]][2])
p=c()
for (i in 1:22){
    bim=read.table(paste('~/haplotypes/chr_',i,'.bim',sep=''))
    chr_i = d[d$chr==paste('chrom',i,sep='_'),]
    p = rbind(p,data.frame(CHR=i,
                SNP=bim[match(chr_i$pos,bim[,4]),2],
                            A1=chr_i$nt1,A2=chr_i$nt2,
                            beta=chr_i$ldpred_beta))
    print(i)
}

p = p[!is.na(p$SNP),]
write.table(p,'~/pgs/PGI_catalog/EVERSMOKE2_hg38.txt',quote=F,row.names=F)

### EA SBayesR ###
gmap = read.table('~/pgs/HM3_SNPs_ChrPosID_hg19tohg38.map',header=T)
d = read.table('~/pgs/PGI_catalog/EA6-single_weights_SBayesR.txt',header=T)
d$pos = gmap[match(d$Chrom.Position,gmap$ChrPosID_hg19),1]
d = d[!is.na(d$pos),]
d$pos = sapply(d$pos,function(x) strsplit(x,':')[[1]][2])
p=c()
for (i in 1:22){
    bim=read.table(paste('~/haplotypes/chr_',i,'.bim',sep=''))
    chr_i = d[d$Chrom==i,]
    p = rbind(p,data.frame(CHR=i,
                SNP=bim[match(chr_i$posbim[,4]),2],
                            A1=chr_i$A1,A2=chr_i$A2,
                            beta=chr_i$A1Effect))
    print(i)
}
p = p[!is.na(p$SNP),]
write.table(p,'~/pgs/PGI_catalog/EA6_hg38.txt',quote=F,row.names=F)

######## Compute #########
pgs.py ~/pgs/EA6 --weights ~/pgs/PGI_catalog/EA6_hg38.txt --bgen ~/haplotypes/chr_@ --imp ~/haplotypes/imputed/chr_@ --beta_col beta --batch_size 1000 --grandpar
pgs.py ~/pgs/EA6_sib --weights ~/pgs/PGI_catalog/EA6_hg38.txt --bgen ~/haplotypes/chr_@ --imp ~/haplotypes/imputed/chr_@ --beta_col beta --batch_size 1000 --fit_sib
pgs.py ~/pgs/EA4/EA4_sib --weights ~/pgs/EA4/EA4_sbayesr.txt --bgen ~/haplotypes/chr_@ --imp ~/haplotypes/imputed/chr_@ --fit_sib --beta_col beta --batch_size 30000
pgs.py ~/pgs/EA4/EA4 --weights ~/pgs/EA4/EA4_sbayesr.txt --bgen ~/haplotypes/chr_@ --imp ~/haplotypes/imputed/chr_@ --beta_col beta --batch_size 30000 --grandpar
pgs.py ~/pgs/externalizing/externalizing --weights ~/pgs/externalizing/externalizing_sbayesr.txt --bgen ~/haplotypes/chr_@ --imp ~/haplotypes/imputed/chr_@ --beta_col beta --batch_size 20000 --grandpar
pgs.py ~/pgs/externalizing/externalizing_sib --weights ~/pgs/externalizing/externalizing_sbayesr.txt --bgen ~/haplotypes/chr_@ --imp ~/haplotypes/imputed/chr_@ --beta_col beta --batch_size 1000 --fit_sib
pgs.py ~/pgs/DEP1 --weights ~/pgs/PGI_catalog/DEP1_hg38.txt --bgen ~/haplotypes/chr_@ --imp ~/haplotypes/imputed/chr_@ --beta_col beta --batch_size 1000 --grandpar
pgs.py ~/pgs/DEP1_sib --weights ~/pgs/PGI_catalog/DEP1_hg38.txt --bgen ~/haplotypes/chr_@ --imp ~/haplotypes/imputed/chr_@ --beta_col beta --batch_size 1000 --fit_sib
pgs.py ~/pgs/ADHD1 --weights ~/pgs/PGI_catalog/ADHD1_hg38.txt --bgen ~/haplotypes/chr_@ --imp ~/haplotypes/imputed/chr_@ --beta_col beta --batch_size 1000 --grandpar
pgs.py ~/pgs/ADHD1_sib --weights ~/pgs/PGI_catalog/ADHD1_hg38.txt --bgen ~/haplotypes/chr_@ --imp ~/haplotypes/imputed/chr_@ --beta_col beta --batch_size 1000 --fit_sib
pgs.py ~/pgs/AFB2 --weights ~/pgs/PGI_catalog/AFB2_hg38.txt --bgen ~/haplotypes/chr_@ --imp ~/haplotypes/imputed/chr_@ --beta_col beta --batch_size 1000 --grandpar
pgs.py ~/pgs/AFB2_sib --weights ~/pgs/PGI_catalog/AFB2_hg38.txt --bgen ~/haplotypes/chr_@ --imp ~/haplotypes/imputed/chr_@ --beta_col beta --batch_size 1000 --fit_sib
pgs.py ~/pgs/EVERSMOKE2 --weights ~/pgs/PGI_catalog/EVERSMOKE2_hg38.txt --bgen ~/haplotypes/chr_@ --imp ~/haplotypes/imputed/chr_@ --beta_col beta --batch_size 1000 --grandpar
pgs.py ~/pgs/EVERSMOKE2_sib --weights ~/pgs/PGI_catalog/EVERSMOKE2_hg38.txt --bgen ~/haplotypes/chr_@ --imp ~/haplotypes/imputed/chr_@ --beta_col beta --batch_size 1000 --fit_sib

for pgs in EA4 externalizing ADHD1 AFB2 EVERSMOKE2 DEP1
do
for i in {1..18}
do
pgs.py ~/pgs/$pgs/$i --pheno ~/phenotypes/processed_traits.txt --pgs ~/pgs/$pgs'.pgs.txt' --gen_models 1-3 --scale_phen --scale_pgs --phen_index $i --ibdrel_path ~/king
pgs.py ~/pgs/$pgs/$i'_sib' --pheno ~/phenotypes/processed_traits.txt --pgs ~/pgs/$pgs'_sib.pgs.txt' --fit_sib --scale_phen --scale_pgs --phen_index $i --ibdrel_path ~/king
done
done