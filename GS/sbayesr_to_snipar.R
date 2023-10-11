### Get SNP IDs from SBayesR to snipar ###
convert_sbayes_r_to_snipar = function(sbayes_out,outfile,bimfile=bimfile){
    d=read.table(sbayes_out,header=T)
    bim=read.table(bimfile,header=F)
    bim=cbind(bim,paste(bim[,1],':',bim[,4],sep=''))
    d$chrpos = paste(d$Chrom,d$Position,sep=':')
    d_out = data.frame(SNP=bim[match(d$chrpos,bim[,7]),2],
                        A1=d$A1,A2=d$A2,beta=d$A1Effect)
    d_out = d_out[!is.na(d_out$SNP),]
    write.table(d_out,outfile,quote=F,row.names=F)
}

setwd('/disk/genetics/sibling_consortium/GS20k/alextisyoung/grandpar/pgs/')
bimfile = '/disk/genetics/sibling_consortium/GS20k/alextisyoung/grandpar/haplotypes/bedfiles/autosome.bim'

ea4_sbayesr = 'EA4_excl_UKBrel_STR_GS_2020_08_21_hm3.snpRes'
ea4_outfile = 'EA4_excl_UKBrel_STR_GS_2020_08_21_hm3.txt'

depression_sbayesr = 'depression/PGC_UKB_depression.snpRes'
depression_outfile = 'depression/PGC_UKB_depression.txt'

height_sbayesr = 'height/height_UKB.snpRes'
height_outfile = 'height/height_UKB.txt'

#convert_sbayes_r_to_snipar(ea4_sbayesr,ea4_outfile,bimfile)
#convert_sbayes_r_to_snipar(depression_sbayesr,depression_outfile,bimfile)
#convert_sbayes_r_to_snipar(height_sbayesr,height_outfile,bimfile)
#convert_sbayes_r_to_snipar('bmi/BMI_UKB.snpRes','bmi/BMI_UKB.txt',bimfile)
convert_sbayes_r_to_snipar('ever_smoke/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR.snpRes',
                            'ever_smoke/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR.txt',bimfile)