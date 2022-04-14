# "/disk/genetics/sibling_consortium/GS20k/alextisyoung/HM3/tidy/traits/15
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

write.table(direct,'direct.sumstats',quote=F,row.names=F)

population = data.frame(CHR=d$chromosome,
                    SNP=d$SNP,
                    A1=d$A1,
                    A2=d$A2,
                    N=d$population_N,
                    Z=d$population_Z,
                    P=10^(-d$population_log10_P))

write.table(population,'population.sumstats',quote=F,row.names=F)

command = paste('/usr/local/bin/python ~/ldsc/ldsc.py --h2 direct.sumstats --ref-ld-chr ../../ibd/ --w-ld-chr ../../ibd/ --out direct')
system(command)

command = paste('/usr/local/bin/python ~/ldsc/ldsc.py --h2 population.sumstats --ref-ld-chr ../../ibd/ --w-ld-chr ../../ibd/ --out population')
system(command)