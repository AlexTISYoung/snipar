setwd('/ludc/Active_Projects/BOTNIA_AYoung_analysis/Private/')
vcf_dir='/ludc/Raw_Data_Archive/Chip_Genotyping/BOTEX_GWAS/Imputation/IMPUTE3_Michigan_imputation/'

for (i in 2:22){
info = read.table(paste(vcf_dir,'chr',i,'.info.gz',sep=''),header=T)
snp_out = info$SNP[info$MAF>0.01 & info$AvgCall>0.99 & info$Rsq>0.99]
write.table(snp_out,paste('haplotypes/chr',i,'MAF_0.01_call_0.99_Rsq_0.99.txt',sep='_'),
            quote=F,row.names=F,col.names=F)
print(i)
print(length(snp_out))}
