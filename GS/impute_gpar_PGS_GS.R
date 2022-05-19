library(rhdf5)
setwd('/disk/genetics/sibling_consortium/GS20k/alextisyoung/grandpar/pgs')

pgsfile = 'GS_EA_13_weights_LDpred_p1.pgs.txt'
impfile = '../imputed/chr_1.hdf5'

pgs = read.table(pgsfile,header=T,stringsAsFactors=F)

# combine pgs with covariates
covars = read.table('../covariates.fam',header=T,stringsAsFactors=F)
covars = covars[match(pgs[,2],covars[,2]),]
covars$pgs = scale(pgs[,3])
write.table(covars,'GS_EA_13_weights_LDpred_p1.pgs.with_covariates.txt',quote=F,row.names=F)

# Match with pedigree
ids = pgs[,2]

ped = t(h5read(impfile,'pedigree'))
dimnames(ped)[[2]] = ped[1,]
ped = ped[2:dim(ped)[1],]
ped = ped[sapply(ped[,1],function(x) substr(x,1,1)!='_'),]

pgs$gpaternal = NA
pgs$gmaternal = NA

dimnames(pgs)[[1]] = ids

# Find imputed grandparents
# estimate r 
bpg = ped[ped[,3]%in%ped[,2] & ped[,4]%in%ped[,2],]
pgs_bpg = pgs[ids%in%bpg[,2],]
r = cor(pgs_bpg$paternal,pgs_bpg$maternal)

a = (1+r)/(1+r/2)

gimpute = matrix(-1,nrow=dim(pgs)[1],ncol=4)
dimnames(gimpute)[[2]] = c('pp','pm','mp','mm')
dimnames(gimpute)[[1]] = dimnames(pgs)[[1]]

for (i in 1:dim(pgs)[1]){
  par_ids = ped[ped[,2]==ids[i],3:4]
  if (par_ids[1]%in%dimnames(pgs)[[1]]){
    pgs_index = match(par_ids[1],dimnames(pgs)[[1]])
    if (par_ids[1]%in%ped[,2]){
      gpar_genotyped = as.integer(ped[match(par_ids[1],ped[,2]),5:6]=='True')
    }    
    scale_factors = c(1,1)*gpar_genotyped+a*c(1,1)*(1-gpar_genotyped)
    pgs$gpaternal[i] = sum(a*pgs[pgs_index,c('paternal','maternal')])
    gimpute[i,1:2] = gpar_genotyped
  } else {pgs$gpaternal[i] = (1+r)*pgs$paternal[i]}
  if (par_ids[2]%in%dimnames(pgs)[[1]]){
    pgs_index = match(par_ids[2],dimnames(pgs)[[1]])
    if (par_ids[2]%in%ped[,2]){
      gpar_genotyped = as.integer(ped[match(par_ids[2],ped[,2]),5:6]=='True')
    }    
    scale_factors = c(1,1)*gpar_genotyped+a*c(1,1)*(1-gpar_genotyped)
    pgs$gmaternal[i] = sum(a*pgs[pgs_index,c('paternal','maternal')])
    gimpute[i,3:4] = gpar_genotyped
  } else {pgs$gmaternal[i] = (1+r)*pgs$maternal[i]}
}

in_bpg = dimnames(pgs)[[1]]%in%bpg[,2]
pgs = pgs[in_bpg,]
gimpute = gimpute[in_bpg,]

# scale
pgs[,-c(1:2)] = pgs[,-c(1:2)]/sd(pgs[,3],na.rm=T)

# Add covariates
covars = covars[match(pgs[,2],covars[,2]),]
pgs_out = cbind(covars[,1:(dim(covars)[2]-1)],pgs[,-c(1:2)])

write.table(pgs_out,'pgs_gpar.txt',quote=F,row.names=F)

