# "/disk/genetics/sibling_consortium/GS20k/alextisyoung/grandpar/traits/15
#python_path = '/ludc/Tools/Software/Python/versions/2.7.16/bin/python'
#ldsc_path = '/ludc/Active_Projects/BOTNIA_AYoung_analysis/Private//ldsc/ldsc.py'

python_path = '/usr/local/bin/python2.7'
ldsc_path = '~/ldsc/ldsc.py'
ldscores = '/disk/genetics/ukb/alextisyoung/hapmap3/haplotypes/imputed_phased/bedfiles/ldscores/'
phenotype =  'height'

d=c()
for (i in 1:22){
d = rbind(d,read.table(paste('chr_',i,'.sumstats.gz',sep=''),header=T))
}

direct = data.frame(CHR=d$chromosome,
                    SNP=d$SNP,
                    A1=d$A1,
                    A2=d$A2,
                    N=d$direct_N,
                    Z=d$direct_Z,
                    P=10^(-d$direct_log10_P))

write.table(direct,paste(phenotype,'direct.sumstats',sep='.'),quote=F,row.names=F)

population = data.frame(CHR=d$chromosome,
                    SNP=d$SNP,
                    A1=d$A1,
                    A2=d$A2,
                    N=d$population_N,
                    Z=d$population_Z,
                    P=10^(-d$population_log10_P))

giant_height = data.frame(CHR=giant$CHR,
                    SNP=giant$SNP,
                    A1=giant$Tested_Allele,
                    A2=giant$Other_Allele,
                    N=giant$N,
                    Z=giant$BETA/giant$SE,
                    P=giant$P)

write.table(giant_height,'Wood.height.sumstats',quote=F,row.names=F)
    

write.table(population,paste(phenotype,'population.sumstats',sep='.'),quote=F,row.names=F)

command = paste(python_path,ldsc_path,'--h2',paste(phenotype,'direct.sumstats',sep='.'),'--ref-ld-chr',
                    ldscores,'--w-ld-chr',ldscores,'--out',paste(phenotype,'direct',sep='.'))
system(command)

command = paste(python_path,ldsc_path,'--h2',paste(phenotype,'population.sumstats',sep='.'),'--ref-ld-chr',
                    ldscores,'--w-ld-chr',ldscores,'--out',paste(phenotype,'population',sep='.'))
system(command)
