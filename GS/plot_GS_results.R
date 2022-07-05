setwd('~/Google Drive/grandparental/GS/')

# Get trait names
traits = read.table('trait_names.txt')
traits[traits[,2]=='EA_quals',2] = 'EA'
traits = traits[-match(c('EA_years','EA_years_mid'),traits[,2]),]
n_phen = dim(traits)[1]
# Collect results
results = data.frame(direct=rep(NA,n_phen),
                     direct_SE=rep(NA,n_phen),
                     paternal=rep(NA,n_phen),
                     paternal_SE=rep(NA,n_phen),
                     maternal=rep(NA,n_phen),
                     maternal_SE=rep(NA,n_phen),
                     grand_paternal=rep(NA,n_phen),
                     grand_paternal_SE=rep(NA,n_phen),
                     grand_maternal=rep(NA,n_phen),
                     grand_maternal_SE=rep(NA,n_phen),
                     parental=rep(NA,n_phen),
                     parental_SE=rep(NA,n_phen),
                     grandparental=rep(NA,n_phen),
                     grandparental_SE=rep(NA,n_phen))

dimnames(results)[[1]] = traits[,2]

results=results[order(-results$direct^2),]

gpar_effect_names = c('proband','paternal','maternal','gpp','gpm','gmp','gmm')

A = matrix(0,nrow=7,ncol=7)
A[1:3,1:3] = diag(3)
A[4,] = c(0,0,0,0.5,0.5,0,0)
A[5,] = c(0,0,0,0,0,0.5,0.5)
A[6,] = c(0,0.5,0.5,0,0,0,0)
A[7,] = c(0,0,0,0.25,0.25,0.25,0.25)

plot_ests = data.frame()
plot_ses = data.frame()

for (i in 1:dim(traits)[1]){
  # Read estimates
  ests = read.table(paste(traits[i,1],3,'effects.txt',sep='.'),row.names=1)
  ests = ests[gpar_effect_names,]
  # Read variance-covariance matrix
  vcov = read.table(paste(traits[i,1],3,'vcov.txt',sep='.'),row.names=1,header=T)
  vcov = vcov[gpar_effect_names,gpar_effect_names]
  # Transform
  r_ests = A%*%ests[,1]
  r_ests = cbind(r_ests,sqrt(diag(A%*%as.matrix(vcov)%*%t(A))))
  # Save in table
  results[i,] = as.vector(t(r_ests))
  # save in format for ggplot2
  plot_ests = rbind(plot_ests,data.frame(phenotype=rep(dimnames(results)[[1]][i],3),
                                         effect=c('direct','parental','grandparental'),
                                         est=r_ests[c(1,6,7),1]))
  plot_ses = rbind(plot_ses,data.frame(phenotype=rep(dimnames(results)[[1]][i],3),
                                        effect=c('direct','parental','grandparental'),
                                        SE=r_ests[c(1,6,7),2]))
}

plot_ests$effect = factor(plot_ests$effect,
  levels=c('direct','parental','grandparental'))

plot_ests$phenotype = factor(plot_ests$phenotype,
                             levels=dimnames(results)[[1]])

plot_lower = plot_ests[,3]+qnorm(0.025)*plot_ses[,3]
plot_upper = plot_ests[,3]-qnorm(0.025)*plot_ses[,3]

write.csv(results,'GS_results.csv',quote=F)

## Plot
library(ggplot2)
gpar_plot = ggplot(plot_ests,aes(x = phenotype,y = est,fill = effect)) +
  geom_col(position = position_dodge(0.75),width=0.75)+
  geom_errorbar(aes(ymin = plot_lower,ymax=plot_upper),position=position_dodge(0.75),width=0.25, size=0.5)+
theme_minimal() + theme(axis.line = element_line(color="black"),
                        axis.ticks = element_line(color="black"),
                        panel.border = element_blank(),
                        axis.text.x = element_text(angle = 45, hjust=1))


ggsave('gpar_plot.pdf',plot=gpar_plot,width=10,height=5)