library(rhdf5)

args = commandArgs(trailingOnly = T)

sib_file = args[1]
#sib_file = 'h2_quad_0.8.hdf5'
effect_file = args[2]
#effect_file = 'h2_quad_0.8.effects.txt'
out = args[3]
#out = 'h2_trio_0.5.meta.hdf5'

if (is.null(args[8])){
  no_sib = FALSE
} else {no_sib = as.logical(args[8])}

beta_list = list()
beta_cov_list = list()

###### Sibs only ######
print('Reading estimates from pGWAS')
xtx = h5read(sib_file,'xtx')
xty = h5read(sib_file,'xty')
sigma2 = h5read(sib_file,'sigma2')

psize = dim(xtx)[1]

beta = matrix(NA,nrow=dim(xty)[2],ncol=psize)
beta_cov = array(NA,dim=c(dim(xty)[2],psize,psize))
beta_se = matrix(NA,nrow=dim(xty)[2],ncol=psize)

for (i in 1:dim(xty)[2]){
  xtx_i = xtx[,,i]
  xty_i = xty[,i]
  if (no_sib){
    xtx_i = xtx_i[-3,-3]
    xty_i = xty_i[-3]
  }
  beta[i,] = solve(xtx_i,xty_i)
  beta_cov[i,,] = sigma2*solve(xtx_i)
  beta_se[i,] = sqrt(diag(beta_cov[i,,]))
}

beta_list[['sibs']] = beta[,-1]
beta_cov_list[['sibs']] = beta_cov[,-1,-1]

###### PO_no_sib ######
print('Reading estimates from poGWAS')
xtx = h5read(po_file_no_sib,'xtx')
xty = h5read(po_file_no_sib,'xty')
sigma2 = h5read(po_file_no_sib,'sigma2')

psize = dim(xtx)[1]

beta = matrix(NA,nrow=dim(xty)[2],ncol=psize)
beta_cov = array(NA,dim=c(dim(xty)[2],psize,psize))
beta_se = matrix(NA,nrow=dim(xty)[2],ncol=psize)

for (i in 1:dim(xty)[2]){
  beta[i,] = solve(xtx[,,i],xty[,i])
  beta_cov[i,,] = sigma2*solve(xtx[,,i])
  beta_se[i,] = sqrt(diag(beta_cov[i,,]))
}

beta_list[['po_no_sib']] = beta[,-1]
beta_cov_list[['po_no_sib']] = beta_cov[,-1,-1]

###### PO_sib ######
print('Reading estimates from poGWAS')
xtx = h5read(po_file_sib,'xtx')
xty = h5read(po_file_sib,'xty')
sigma2 = h5read(po_file_sib,'sigma2')

psize = dim(xtx)[1]

beta = matrix(NA,nrow=dim(xty)[2],ncol=psize)
beta_cov = array(NA,dim=c(dim(xty)[2],psize,psize))
beta_se = matrix(NA,nrow=dim(xty)[2],ncol=psize)

for (i in 1:dim(xty)[2]){
  beta[i,] = solve(xtx[,,i],xty[,i])
  beta_cov[i,,] = sigma2*solve(xtx[,,i])
  beta_se[i,] = sqrt(diag(beta_cov[i,,]))
}

beta_list[['po_sib']] = beta[,-1]
beta_cov_list[['po_sib']] = beta_cov[,-1,-1]

###### both parents (no_sib) ######
print('Reading estimates from triGWAS')
xtx = h5read(bpg_file_no_sib,'xtx')
xty = h5read(bpg_file_no_sib,'xty')
sigma2 = h5read(bpg_file_no_sib,'sigma2')

psize = dim(xtx)[1]

beta = matrix(NA,nrow=dim(xty)[2],ncol=psize)
beta_cov = array(NA,dim=c(dim(xty)[2],psize,psize))
beta_se = matrix(NA,nrow=dim(xty)[2],ncol=psize)

for (i in 1:dim(xty)[2]){
  beta[i,] = solve(xtx[,,i],xty[,i])
  beta_cov[i,,] = sigma2*solve(xtx[,,i])
  beta_se[i,] = sqrt(diag(beta_cov[i,,]))
}

beta_list[['bpg_no_sib']] = beta[,-1]
beta_cov_list[['bpg_no_sib']] = beta_cov[,-1,-1]

###### both parents (sib) ######
print('Reading estimates from triGWAS')
xtx = h5read(bpg_file_sib,'xtx')
xty = h5read(bpg_file_sib,'xty')
sigma2 = h5read(bpg_file_sib,'sigma2')

psize = dim(xtx)[1]

beta = matrix(NA,nrow=dim(xty)[2],ncol=psize)
beta_cov = array(NA,dim=c(dim(xty)[2],psize,psize))
beta_se = matrix(NA,nrow=dim(xty)[2],ncol=psize)

for (i in 1:dim(xty)[2]){
  beta[i,] = solve(xtx[,,i],xty[,i])
  beta_cov[i,,] = sigma2*solve(xtx[,,i])
  beta_se[i,] = sqrt(diag(beta_cov[i,,]))
}

beta_list[['bpg_sib']] = beta[,-1]
beta_cov_list[['bpg_sib']] = beta_cov[,-1,-1]

#### Meta analysis
print('Meta analysing effects')
  meta = array(NA,dim=c(4,dim(xtx)[3],2))
  dimnames(meta)[[1]] = c('direct','sib','paternal','maternal')

# Direct
for (i in 1:dim(xty)[2]){
  delta_vec = c(beta_list[['sibs']][i,1],beta_list[['po_no_sib']][i,1],beta_list[['po_sib']][i,1],
                beta_list[['bpg_no_sib']][i,1],beta_list[['bpg_sib']][i,1])
  delta_se_vec = sqrt(c(beta_cov_list[['sibs']][i,1,1],beta_cov_list[['po_no_sib']][i,1,1],beta_cov_list[['po_sib']][i,1,1]
                      ,beta_cov_list[['bpg_no_sib']][i,1,1],beta_cov_list[['bpg_sib']][i,1,1]))
  weights = 1/delta_se_vec^2
  meta[1,i,1] = sum(weights*delta_vec)/sum(weights)
  meta[1,i,2] = sqrt(sum(weights)^(-1))
}
# Sib
for (i in 1:dim(xty)[2]){
  delta_vec = c(beta_list[['sibs']][i,2],beta_list[['po_sib']][i,2],
                beta_list[['bpg_sib']][i,2])
  delta_se_vec = sqrt(c(beta_cov_list[['sibs']][i,2,2],beta_cov_list[['po_sib']][i,2,2]
                      ,beta_cov_list[['bpg_sib']][i,2,2]))
  weights = 1/delta_se_vec^2
  meta[2,i,1] = sum(weights*delta_vec)/sum(weights)
  meta[2,i,2] = sqrt(sum(weights)^(-1))
}

##### Check #####
effects = read.table(effect_file)
# Direct
r = lm(meta[1,,1]~effects[,1])
s = summary(r)$coefficients
print(paste('bias for direct effects: ',
            round(s[2,1]-1,4),' (',round(s[2,2],4),' S.E.)',sep=''))
# Sib
r = lm(meta[2,,1]~effects[,2])
s = summary(r)$coefficients
print(paste('bias for sib effects: ',
            round(s[2,1]-1,4),' (',round(s[2,2],4),' S.E.)',sep=''))

h5out = paste(out,'hdf5',sep='.')
h5createFile(h5out)
h5write(meta,h5out,'meta')
H5close()

##### Meta analysis ######
# print('Meta analysing effects')
# A = rbind(c(1,0,0),c(0,0.5,0.5))
# beta_meta = matrix(NA,nrow=dim(xty)[2],ncol=3)
# beta_cov_meta = array(NA,dim=c(dim(xty)[2],3,3))
# beta_se_meta =  matrix(NA,nrow=dim(xty)[2],ncol=3)
# for (i in 1:dim(xty)[2]){
#   sigma_inv = list()
#   sigma_inv[['sibs']] = solve(beta_cov_list[['sibs']][i,,])
#   sigma_inv[['po']] = solve(beta_cov_list[['po']][i,,])
#   sigma_inv[['bpg']] = solve(beta_cov_list[['bpg']][i,,])
#   beta_cov_meta[i,,] = solve(t(A)%*%sigma_inv[['sibs']]%*%A+sigma_inv[['po']]+sigma_inv[['bpg']])
#   beta_meta[i,] = beta_cov_meta[i,,]%*%(t(A)%*%sigma_inv[['sibs']]%*%beta_list[['sibs']][i,]+
#                                           sigma_inv[['po']]%*%beta_list[['po']][i,]+
#                                           sigma_inv[['bpg']]%*%beta_list[['bpg']][i,])
#   beta_se_meta[i,] = sqrt(diag(beta_cov_meta[i,,]))
# }