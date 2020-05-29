library(rhdf5)

args = commandArgs(trailingOnly = T)

sib_file = args[1]
effect_file = args[2]
out = args[3]

sib_file = 'test_data/h2_quad_0.8.hdf5'
effect_file = 'test_data/h2_quad_0.8.effects.txt'

###### Sibs only ######
print('Reading estimates')
beta = t(h5read(sib_file,'estimate'))
beta_se = t(h5read(sib_file,'estimate_ses'))

##### Check #####
effects = read.table(effect_file)
# Direct
r = lm(beta[,1]~effects[,1])
s = summary(r)$coefficients
print(paste('bias for direct effects: ',
            round(s[2,1]-1,4),' (',round(s[2,2],4),' S.E.)',sep=''))
# paternal
r = lm(beta[,2]~effects[,2])
s = summary(r)$coefficients
print(paste('bias for paternal effects: ',
            round(s[2,1]-1,4),' (',round(s[2,2],4),' S.E.)',sep=''))

# maternal
r = lm(beta[,3]~effects[,3])
s = summary(r)$coefficients
print(paste('bias for maternal effects: ',
            round(s[2,1]-1,4),' (',round(s[2,2],4),' S.E.)',sep=''))

z = (beta-effects)/beta_se
p = 1-pchisq(sum(z^2),dim(z)[1]*dim(z)[2])
print(paste('Chi-square test p-value: ',round(p,digits=4),sep=''))