### Get SNP IDs from SBayesR to snipar ###
setwd('/var/genetics/data/fhs/private/v32/raw/pgs')
d=read.table('EA4_excl_UKBrel_STR_GS_2020_08_21_hm3.snpRes',header=T)
bim=read.table('/var/genetics/data/fhs/private/v32/raw/extracted/v20/bedfiles/autosome_filtered.bim',header=F)
bim=cbind(bim,paste(bim[,1],':',bim[,4],sep=''))
d$chrpos = paste(d$Chrom,d$Position,sep=':')
d_out = data.frame(SNP=d$chrpos,
                    A1=d$A1,A2=d$A2,beta=d$A1Effect)
write.table(d_out,'/var/genetics/data/fhs/private/v32/raw/pgs/EA4_excl_UKBrel_STR_GS_2020_08_21_hm3.weights.txt',quote=F,row.names=F)
