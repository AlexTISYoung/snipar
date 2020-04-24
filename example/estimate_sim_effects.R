library(rhdf5)

args = commandArgs(trailingOnly = T)

sib_file = args[1]
effect_file = args[2]
out = args[3]

sib_file = 'h2_quad_0.8.hdf5'
effect_file = 'h2_quad_0.8.effects.txt'
out = 'h2_quad_0.8.effects.hdf5'

###### Sibs only ######
print('Reading estimates')
xtx = h5read(sib_file,'xtx')
xty = h5read(sib_file,'xty')

psize = dim(xtx)[1]

beta = matrix(NA,nrow=dim(xty)[2],ncol=psize)
beta_cov = array(NA,dim=c(dim(xty)[2],psize,psize))
beta_se = matrix(NA,nrow=dim(xty)[2],ncol=psize)

for (i in 1:dim(xty)[2]){
  xtx_i = xtx[,,i]
  xty_i = xty[,i]
  beta[i,] = solve(xtx_i,xty_i)
  beta_cov[i,,] = solve(xtx_i)
  beta_se[i,] = sqrt(diag(beta_cov[i,,]))
}


##### Check #####
effects = read.table(effect_file)
# Direct
r = lm(beta[,2]~effects[,1])
s = summary(r)$coefficients
print(paste('bias for direct effects: ',
            round(s[2,1]-1,4),' (',round(s[2,2],4),' S.E.)',sep=''))
# paternal
r = lm(beta[,3]~effects[,2])
s = summary(r)$coefficients
print(paste('bias for paternal effects: ',
            round(s[2,1]-1,4),' (',round(s[2,2],4),' S.E.)',sep=''))

# maternal
r = lm(beta[,4]~effects[,3])
s = summary(r)$coefficients
print(paste('bias for maternal effects: ',
            round(s[2,1]-1,4),' (',round(s[2,2],4),' S.E.)',sep=''))

# Chi-square test
z = (beta[,-1]-effects)/beta_se[,-1]
p = 1-pchisq(sum(z^2),dim(z)[1]*dim(z)[2])
print(paste('Chi-square test p-value: ',round(p,4),sep=''))