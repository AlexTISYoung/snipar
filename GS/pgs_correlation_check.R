library(rhdf5)
setwd('/disk/genetics/sibling_consortium/GS20k/alextisyoung/grandpar/pgs')

## Read control PGS
pgsfile_o = 'EA4_excl_UKBrel_STR_GS_2020_08_21_hm3.pgs.control_sibling.txt'

pgs_o = read.table(pgsfile_o,header=T,stringsAsFactors=F)
pgs_o = pgs_o[order(pgs_o[,1]),]

# Restrict to families with two genotyped siblings
pgs_o[,1] = sapply(pgs_o[,1],function(x) strsplit(x,'_o_')[[1]][2])
fams = unique(pgs_o[,1])
fam_counts = sapply(fams,function(x) sum(pgs_o[,1]==x))
f2 = fams[fam_counts==2]
pgs_o = pgs_o[pgs_o[,1]%in%f2,]

## Read true parental PGS
pgsfile = 'EA4_excl_UKBrel_STR_GS_2020_08_21_hm3.pgs.txt'
pgs = read.table(pgsfile,header=T,stringsAsFactors=F)

pgs_o = cbind(pgs_o, pgs[match(pgs_o$IID,pgs$IID),c('paternal','maternal')])
pgs_o$parental_obs = pgs_o$paternal+pgs_o$maternal

r = cor(pgs_o$paternal,pgs_o$maternal)

## Compute correlations
R_obs = cor(pgs_o[,c('proband','parental','parental_obs')])

# predicted
R_pred = rbind(c(1,sqrt((2+r)/3),sqrt((1+r)/2)),
               c(sqrt((2+r)/3),1,sqrt((3/2)*(1+r)/(2+r))),
               c(sqrt((1+r)/2),sqrt((3/2)*(1+r)/(2+r)),1)              
               )

corr_ratio = R_obs/R_pred

write.table(corr_ratio,'observed_to_predicted_correlation_ratio.txt',quote=F)

## Compute variance ratios
vars = apply(pgs_o[,c('proband','parental','parental_obs')],2,var)
var_ratios = c(vars[2]/vars[1],vars[3]/vars[2],vars[3]/vars[1])
predicted_ratios = c(1.5*(1+r/2),4*(1+r)/(3*(1+r/2)),2*(1+r))

observed_to_predicted_ratios = var_ratios/predicted_ratios