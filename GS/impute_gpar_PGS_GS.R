library(rhdf5)
setwd('/disk/genetics/sibling_consortium/GS20k/alextisyoung/grandpar/pgs')

opg_am_adj = function(pgi_imp,pgi_obs,r,n){
  rcoef = r/((2^n)*(1+r)-r)
  return(rcoef*pgi_obs+(1+rcoef)*pgi_imp)
}

npg_am_adj = function(pgi_imp,r,n){
  rcoef = (1+r)/(1+(1-(1/2)^(n-1))*r)
  return(rcoef*pgi_imp)
}

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

# Fam sizes
fams = unique(ped[,1])
fsizes = sapply(fams,function(x) sum(ped[,1]==x))
names(fsizes) = fams

pgs$gpp = NA
pgs$gpm = NA
pgs$gmp = NA 
pgs$gmm = NA

dimnames(pgs)[[1]] = ids

# Find imputed grandparents
# estimate r 
bpg = ped[ped[,3]%in%ped[,2] & ped[,4]%in%ped[,2],]
pgs_bpg = pgs[ids%in%bpg[,2],]
r = cor(pgs_bpg$paternal,pgs_bpg$maternal)

gimpute = matrix(NA,nrow=dim(pgs)[1],ncol=4)
dimnames(gimpute)[[2]] = c('pp','pm','mp','mm')
dimnames(gimpute)[[1]] = dimnames(pgs)[[1]]

for (i in 1:dim(pgs)[1]){
  par_ids = ped[ped[,2]==pgs[i,2],3:4]
  # If father genotyped and in pedigree, compute paternal grandparental PGIs
  if (par_ids[1]%in%pgs[,2] & par_ids[1]%in%ped[,2]){
      pgs_index = match(par_ids[1],dimnames(pgs)[[1]])
      father_ped = ped[match(par_ids[1],ped[,2]),]
      father_fsize = fsizes[father_ped[1]]
      gpar_genotyped = ped[match(par_ids[1],ped[,2]),5:6]=='True'
      gimpute[i,1:2] = gpar_genotyped
      # If father genotyped, set as paternal grandfather pgs
      if (gpar_genotyped[1]){
        pgs$gpp[i] = pgs[pgs_index,'paternal']
      } else {
        # If mother genotyped, use adjustment given imputation with one parent; 
        # otherwise, use adjustment for imputation without parents
        if (gpar_genotyped[2]){
          pgs$gpp[i] = opg_am_adj(pgs[pgs_index,'paternal'],pgs[pgs_index,'maternal'],r,father_fsize)
        } else {
          pgs$gpp[i] = npg_am_adj(sum(pgs[pgs_index,c('paternal','maternal')]),r,father_fsize)/2
        }  
      }
      # If mother genotyped, set as paternal grandmother pgs
      if (gpar_genotyped[2]){
        pgs$gpm[i] = pgs[pgs_index,'maternal']
      } else {
        # If father genotyped, use adjustment given imputation with one parent; 
        # otherwise, use adjustment for imputation without parents
        if (gpar_genotyped[1]){
          pgs$gpm[i] = opg_am_adj(pgs[pgs_index,'maternal'],pgs[pgs_index,'paternal'],r,father_fsize)
        } else {
          pgs$gpm[i] = npg_am_adj(sum(pgs[pgs_index,c('paternal','maternal')]),r,father_fsize)/2
        }  
      }
  } else {pgs[i,c('gpp','gpm')] = (1+r)*pgs[i,c('paternal')]/2
          gimpute[i,1:2] = FALSE
  }
  # If mother genotyped and in pedigree, compute maternal grandparental PGIs
  if (par_ids[2]%in%dimnames(pgs)[[1]] & par_ids[2]%in%ped[,2]){
      pgs_index = match(par_ids[2],dimnames(pgs)[[1]])
      mother_ped = ped[match(par_ids[2],ped[,2]),]
      mother_fsize = fsizes[mother_ped[1]]
      gpar_genotyped = ped[match(par_ids[2],ped[,2]),5:6]=='True'
      gimpute[i,3:4] = gpar_genotyped
      # If father genotyped, set as paternal grandfather pgs
      if (gpar_genotyped[1]){
        pgs$gmp[i] = pgs[pgs_index,'paternal']
      } else {
        # If mother genotyped, use adjustment given imputation with one parent; 
        # otherwise, use adjustment for imputation without parents
        if (gpar_genotyped[2]){
          pgs$gmp[i] = opg_am_adj(pgs[pgs_index,'paternal'],pgs[pgs_index,'maternal'],r,mother_fsize)
        } else {
          pgs$gmp[i] = npg_am_adj(sum(pgs[pgs_index,c('paternal','maternal')]),r,mother_fsize)/2
        }  
      }
      # If mother genotyped, set as paternal grandmother pgs
      if (gpar_genotyped[2]){
        pgs$gmm[i] = pgs[pgs_index,'maternal']
      } else {
        # If father genotyped, use adjustment given imputation with one parent; 
        # otherwise, use adjustment for imputation without parents
        if (gpar_genotyped[1]){
          pgs$gmm[i] = opg_am_adj(pgs[pgs_index,'maternal'],pgs[pgs_index,'paternal'],r,mother_fsize)
        } else {
          pgs$gmm[i] = npg_am_adj(sum(pgs[pgs_index,c('paternal','maternal')]),r,mother_fsize)/2
        }  
      }
  } else {
    pgs[i,c('gmp','gmm')] = (1+r)*pgs[i,c('maternal')]/2
    gimpute[i,3:4] = FALSE
  }
}

in_bpg = dimnames(pgs)[[1]]%in%bpg[,2]
pgs = pgs[in_bpg,]
gimpute = gimpute[in_bpg,]

# scale
pgs[,-c(1:2)] = pgs[,-c(1:2)]/sd(pgs[,3],na.rm=T)

# Add covariates
covars = covars[match(pgs[,2],covars[,2]),]
pgs_out = cbind(covars[,1:(dim(covars)[2]-1)],pgs[,-c(1:2)])

write.table(pgs_out,'pgs_gpar_full.txt',quote=F,row.names=F)

## Load effects
effects = read.table('fpgs_gpar.effects.txt',header=F,row.names=1)[29:35,]
v = as.matrix(read.table('fpgs_gpar.vcov.txt')[29:35,29:35])

A = rbind(diag(7),
          c(0,0.5,0.5,0,0,0,0),
          c(0,0,0,0.25,0.25,0.25,0.25))

all_effects = A%*%effects[,1]
all_effects = cbind(all_effects,sqrt(diag(A%*%v%*%t(A))))
dimnames(all_effects)[[1]] = c(dimnames(effects)[[1]],'parental','grandparental')
all_effects = cbind(all_effects,all_effects[,1]/all_effects[,2])
all_effects = cbind(all_effects,1-pchisq(all_effects[,3]^2,1))
dimnames(all_effects)[[2]] = c('est','S.E.','t','P')

