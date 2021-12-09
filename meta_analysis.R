###########################################################################################################################
# Script for meta-analysing PGI effect estimates (direct, paternal, maternal) from UK Biobank, Generation Scotland, and STR
###########################################################################################################################
####### Input files ########
# UK Biobank estimates (output from EA4_PGI_analysis.py) in directory UKB/ with paths UKB/gls_1_effects.txt, UKB/gls_1_effects.txt, ...
# UK Biobank variance-covariance matrices of effect estimates in direct UKB/ with paths UKB/gls_1_vcov.txt, UKB/gls_1_vcov.txt, ...
# Generation Scotland estimates (output from EA4_PGI_analysis.py) in directory GS/ with paths GS/gls_1_effects.txt, GS/gls_1_effects.txt, ...
# Generation Scotland variance-covariance matrices of effect estimates in direct GS/ with paths GS/gls_1_vcov.txt, GS/gls_1_vcov.txt, ...
# STR estimates (output from fPGS.py with --parsum option) in directory STR_EA4_allscores/PGS_EA_5_LDpred_p1/
############################

# Function to approximate variance of ratio
var_ratio_approx = function(x,y,vx,vy,cxy){
  x2 = x^2; y2 = y^2
  return((x2/y2)*(vx/x2-2*cxy/(x*y)+vy/y2))
}

traits = list()
files = list()
ests = list()
vcovs = list()
parcor = list()

#### UKB traits
traits[['UKB']] = c('SWB','glucose','non-HDL','HDL','SBP','DBP','SRH','FEV1','ever-smoked','CPD',
                    'DPW','height','BMI','cognition','neuroticism','AAFB','NC_M','NC_F','household income',
                    'depression','EA','hourly income','NC')

files[['UKB']] = paste('UKB/gls_',c(1:23),sep='')

ests[['UKB']] = matrix(NA,nrow=23,ncol=3)
dimnames(ests[['UKB']])[[1]] = traits[['UKB']]

vcovs[['UKB']] = array(NA,dim=c(23,3,3))
dimnames(vcovs[['UKB']])[[1]] = traits[['UKB']]

parcor[['UKB']] = 1.899243488793618817e-01

for (i in 1:23){
  ests[['UKB']][i,] = read.table(paste(files[['UKB']][i],'effects.txt',sep='_'))[1:3,1]
  vcovs[['UKB']][i,,] = as.matrix(read.table(paste(files[['UKB']][i],'vcov.txt',sep='_')))
}

##### GS traits

traits[['GS']] = c('glucose','non-HDL','HDL','FEV1','height','BMI','CPD','cognition','vocab','EA_old','SBP','DBP','neuroticism','EA')

files[['GS']] = paste('GS/gls_',c(1:14),sep='')

ests[['GS']] = matrix(NA,nrow=14,ncol=3)
dimnames(ests[['GS']])[[1]] = traits[['GS']]

vcovs[['GS']] = array(NA,dim=c(14,3,3))
dimnames(vcovs[['GS']])[[1]] = traits[['GS']]

parcor[['GS']] = 1.815641856948789767e-01

for (i in 1:14){
  ests[['GS']][i,] = read.table(paste(files[['GS']][i],'effects.txt',sep='_'))[1:3,1]
  vcovs[['GS']][i,,] = as.matrix(read.table(paste(files[['GS']][i],'vcov.txt',sep='_')))
}


##### STR traits

traits[['STR']] = c('BMI','CPD','depression','DPW','EA','ever-smoked','height','household income','NC','SRH','cognition','SWB')

files[['STR']] = paste('STR_EA4_allscores/PGS_EA_5_LDpred_p1/',
                       c('bmi','cigsperday','depressiveSymptoms','drinks','eduYears','eversmoker','height','hhInc','numberChildren','SRH','stdIQ','SWB'),
                       sep='')

ests[['STR']] = matrix(NA,nrow=12,ncol=2)
dimnames(ests[['STR']])[[1]] = traits[['STR']]

vcovs[['STR']] = array(NA,dim=c(12,2,2))
dimnames(vcovs[['STR']])[[1]] = traits[['STR']]

parcor[['STR']] = 2*0.5877-1

for (i in 1:12){
  est_i = read.table(paste(files[['STR']][i],'_EA_5_LDpred_p1.sibdiff.pgs_effects.txt',sep=''),row.names=1)
  v = matrix(0,nrow=2,ncol=2)
  diag(v) = est_i[,2]^2
  vcovs[['STR']][i,,] = v 
  ests[['STR']][i,] = est_i[,1]
}

###### Meta-analysis 
tlist = unique(c(traits[['UKB']],traits[['GS']],traits[['STR']]))
tlist = tlist[!tlist=='EA_old']
meta_effects = array(NA,dim=c(length(tlist),3))
dimnames(meta_effects)[[1]] = tlist
meta_vcov = array(NA,dim=c(length(tlist),3,3))
dimnames(meta_vcov)[[1]] = tlist
datasets = matrix(FALSE,nrow=length(tlist),ncol=3)
dimnames(datasets)[[1]] = tlist
dimnames(datasets)[[2]] = c('UKB','GS','STR')
c_sib_STR = (1+0.5*(1-parcor[['STR']])/(1+parcor[['STR']]))^(-1)
A = rbind(c(1,0,0),c(1,c_sib_STR,c_sib_STR))
for (i in 1:24){
  t = tlist[i]
  V = 0
  theta = 0
  if (t%in%traits[['UKB']]){
    v = solve(vcovs[['UKB']][t,,])
    V = V+v
    theta =theta+v%*%ests[['UKB']][t,]
    print(mean(abs(diag(v))))
    datasets[t,'UKB'] = T
  }
  if (t%in%traits[['GS']]){
    v = solve(vcovs[['GS']][t,,])
    V = V+v
    theta =theta+v%*%ests[['GS']][t,]
    datasets[t,'GS'] = T
  }
  if (t%in%traits[['STR']]){
    v = solve(vcovs[['STR']][t,,])
    V = V+t(A)%*%v%*%A
    theta = theta+t(A)%*%v%*%ests[['STR']][t,]
    datasets[t,'STR'] = T
  }
  meta_vcov[t,,] = solve(V)
  meta_effects[t,] = meta_vcov[t,,]%*%theta
}

r_meta = mean(c(parcor[['UKB']],parcor[['GS']],parcor[['STR']]))

B = rbind(c(1,0,0),c(1,(1+r_meta)/2,(1+r_meta)/2),c(0,(1+r_meta)/2,(1+r_meta)/2),c(0,-1,1))
results = array(NA,dim=c(dim(meta_effects)[1],10))
results[,c(1,3,5,7)] = meta_effects%*%t(B)
results[,9] = results[,1]/results[,3]
dimnames(results)[[1]] = dimnames(meta_effects)[[1]]
dimnames(results)[[2]] = c('direct','direct_SE','population','population_SE','population-direct','population-direct_SE','maternal-paternal','maternal-paternal_SE',
                           'direct.population.ratio','direct.population.ratio_SE')
for (i in 1:dim(results)[1]){
  vc = B%*%meta_vcov[i,,]%*%t(B)
  results[i,c(2,4,6,8)] = sqrt(diag(vc))
  results[i,10] = sqrt(var_ratio_approx(results[i,'direct'],results[i,'population'],
                                   results[i,'direct_SE']^2,results[i,'population_SE']^2,
                                   vc[1,2]))
}
results = data.frame(results)

results$population.direct.P = -log10(exp(1))*pchisq((results[,'population.direct']/results[,'population.direct_SE'])^2, 1, lower.tail = FALSE, log.p = TRUE)
results$maternal.paternal.P = -log10(exp(1))*pchisq((results[,'maternal.paternal']/results[,'maternal.paternal_SE'])^2, 1, lower.tail = FALSE, log.p = TRUE)

results$phenotype = dimnames(results)[[1]]
presults = data.frame()
for (i in 1:dim(results)[1]){
  results_i = data.frame(trait=rep(results$phenotype[i],2),
                         effect = c('direct','population'),
                         est = c(results[i,'direct'],results[i,'population']),
                         se = c(results[i,'direct_SE'],results[i,'population_SE']))
  presults = rbind(presults,results_i)
}

presults$lower = presults$est+qnorm(0.025)*presults$se
presults$upper = presults$est-qnorm(0.025)*presults$se
presults$effect = factor(presults$effect,levels = c('direct', 'population'))
ass_effects = presults[presults$effect=='population','est']^2
names(ass_effects) = presults[presults$effect=='population','trait']
ass_effects = ass_effects[order(-ass_effects)]
presults$trait = factor(presults$trait,levels=names(ass_effects))

library(ggplot2)
pgsplot=ggplot(presults, aes(fill=effect, y=est, x=trait)) + 
  geom_bar(position="dodge", stat="identity")+
  geom_errorbar( aes(ymin=lower, ymax=upper),position=position_dodge(.9),
                 width=0.5)+theme(axis.text.x = element_text(angle = 90))

ggsave('pgs_results.pdf',pgsplot,width=10,height=5)


################ Results Table ###################
library('sjPlot')
results=results[order(-results$population^2),]
results = results[,c(dim(results)[2],1:(dim(results)[2]-1))]
results$phenotype[match("household income",results$phenotype)] ='household_income'
results$phenotype[match("hourly income",results$phenotype)] ='hourly_income'
write.table(results,'pgs_results_table.txt',quote=F,row.names=F)
tab_df(results,title='PGS Decomposition',file='pgs_results_table.doc',digits=4)

tab_df(data.frame(phenotype=dimnames(datasets)[[1]],datasets),title='Datasets for meta-analysis',file='datasets.doc')




