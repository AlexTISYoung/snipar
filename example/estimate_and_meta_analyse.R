library(rhdf5)

args = commandArgs(trailingOnly = T)

sib_file = args[1]
po_file = args[2]
bpg_file = args[3]
effect_file = args[4]
out = args[5]

if (is.null(args[6])){
  no_sib = FALSE
} else {no_sib = args[6]}

beta_list = list()
beta_cov_list = list()

###### Sibs only ######
print('Reading estimates from pGWAS')
xtx = h5read(sib_file,'xtx')
xty = h5read(sib_file,'xty')
sigma2 = h5read(sib_file,'sigma2')

beta = matrix(NA,nrow=dim(xty)[2],ncol=3)
beta_cov = array(NA,dim=c(dim(xty)[2],3,3))
beta_se = matrix(NA,nrow=dim(xty)[2],ncol=3)

for (i in 1:dim(xty)[2]){
  xtx_i = xtx[,,i]
  xty_i = xty[,i]
  if (!no_sib){
    xtx_i = xtx_i[-3,-3]
    xty_i = xty_i[-3]
  }
  beta[i,] = solve(xtx_i,xty_i)
  beta_cov[i,,] = sigma2*solve(xtx_i)
  beta_se[i,] = sqrt(diag(beta_cov[i,,]))
}

beta_list[['sibs']] = beta[,-1]
beta_cov_list[['sibs']] = beta_cov[,-1,-1]

###### PO ######
print('Reading estimates from poGWAS')
xtx = h5read(po_file,'xtx')
xty = h5read(po_file,'xty')
sigma2 = h5read(po_file,'sigma2')

beta = matrix(NA,nrow=dim(xty)[2],ncol=4)
beta_cov = array(NA,dim=c(dim(xty)[2],4,4))
beta_se = matrix(NA,nrow=dim(xty)[2],ncol=4)

for (i in 1:dim(xty)[2]){
  beta[i,] = solve(xtx[,,i],xty[,i])
  beta_cov[i,,] = sigma2*solve(xtx[,,i])
  beta_se[i,] = sqrt(diag(beta_cov[i,,]))
}

beta_list[['po']] = beta[,-1]
beta_cov_list[['po']] = beta_cov[,-1,-1] 

###### both parents ######
print('Reading estimates from triGWAS')
xtx = h5read(bpg_file,'xtx')
xty = h5read(bpg_file,'xty')
sigma2 = h5read(bpg_file,'sigma2')

beta = matrix(NA,nrow=dim(xty)[2],ncol=4)
beta_cov = array(NA,dim=c(dim(xty)[2],4,4))
beta_se = matrix(NA,nrow=dim(xty)[2],ncol=4)

for (i in 1:dim(xty)[2]){
  beta[i,] = solve(xtx[,,i],xty[,i])
  beta_cov[i,,] = sigma2*solve(xtx[,,i])
  beta_se[i,] = sqrt(diag(beta_cov[i,,]))
}

beta_list[['bpg']] = beta[,-1]
beta_cov_list[['bpg']] = beta_cov[,-1,-1] 

##### Meta analysis ######
print('Meta analysing effects')
A = rbind(c(1,0,0),c(0,0.5,0.5))
beta_meta = matrix(NA,nrow=dim(xty)[2],ncol=3)
beta_cov_meta = array(NA,dim=c(dim(xty)[2],3,3))
beta_se_meta =  matrix(NA,nrow=dim(xty)[2],ncol=3)
for (i in 1:dim(xty)[2]){
  sigma_inv = list()
  sigma_inv[['sibs']] = solve(beta_cov_list[['sibs']][i,,])
  sigma_inv[['po']] = solve(beta_cov_list[['po']][i,,])
  sigma_inv[['bpg']] = solve(beta_cov_list[['bpg']][i,,])
  beta_cov_meta[i,,] = solve(t(A)%*%sigma_inv[['sibs']]%*%A+sigma_inv[['po']]+sigma_inv[['bpg']])
  beta_meta[i,] = beta_cov_meta[i,,]%*%(t(A)%*%sigma_inv[['sibs']]%*%beta_list[['sibs']][i,]+
                                          sigma_inv[['po']]%*%beta_list[['po']][i,]+
                                          sigma_inv[['bpg']]%*%beta_list[['bpg']][i,]) 
  beta_se_meta[i,] = sqrt(diag(beta_cov_meta[i,,]))
}

##### Check #####
effects = read.table(effect_file)
effect_names = c('direct','paternal','maternal')
for (i in 1:3){
  r = lm(beta_meta[,i]~effects[,i])
  s = summary(r)$coefficients
  print(paste('bias for ',effect_names[i],' effects: ',
              round(s[2,1]-1,3),' (',round(s[2,2],4),' S.E.)',sep=''))
}

h5out = paste(out,'hdf5',sep='.')
h5createFile(h5out)
h5write(beta_meta,h5out)
h5write(beta_se_meta,h5out)
h5write(beta_cov_meta,h5out)
H5close()