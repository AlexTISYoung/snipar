library(rhdf5)

beta_list = list()
beta_cov_list = list()

sib_file = 'sibs/sim_60k_parental.hdf5'
po_file = 'one_parent_genotyped/sim_60k_parental.hdf5'
bpg_file = 'both_parents_genotyped/sim_60k_parental.hdf5'

###### Sibs only ######
xtx = h5read(sib_file,'xtx')
xty = h5read(sib_file,'xty')
sigma2 = h5read(sib_file,'sigma2')

beta = matrix(NA,nrow=dim(xty)[2],ncol=3)
beta_cov = array(NA,dim=c(dim(xty)[2],3,3))
beta_se = matrix(NA,nrow=dim(xty)[2],ncol=3)

for (i in 1:dim(xty)[2]){
  beta[i,] = solve(xtx[-3,-3,i],xty[-3,i])
  beta_cov[i,,] = sigma2*solve(xtx[-3,-3,i])
  beta_se[i,] = sqrt(diag(beta_cov[i,,]))
}

beta_list[['sibs']] = beta[,-1]
beta_cov_list[['sibs']] = beta_cov[,-1,-1]

###### PO ######
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
effects = read.table('sim_60k_parental.effects.txt')
r = lm(beta_meta[,1]~effects[,2])
plot(r)
z = (beta_meta[,1]-effects[,2])/beta_se_meta[,1]
1-pchisq(sum(z^2),length(z))

