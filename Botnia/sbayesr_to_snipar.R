### Get SNP IDs from SBayesR to snipar ###
setwd('/ludc/Active_Projects/BOTNIA_AYoung_analysis/Private/pgs/')
d=read.table('EA4_hm3.snpRes',header=T)
bim=read.table('../haplotypes/bedfiles/autosome.bim',header=F)
bim=cbind(bim,paste(bim[,1],':',bim[,4],sep=''))
d$chrpos = paste(d$Chrom,d$Position,sep=':')
d_out = data.frame(SNP=bim[match(d$chrpos,bim[,7]),2],
                    A1=d$A1,A2=d$A2,beta=d$A1Effect)
write.table(d_out,'EA4_hm3.txt',quote=F,row.names=F)
