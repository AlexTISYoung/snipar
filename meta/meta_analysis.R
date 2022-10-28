setwd('~/Google Drive/grandparental/')

# Function for fixed effects meta-analysis
fe_meta = function(ests,ses){
  weights = 1/ses^2
  est = sum(weights*ests,na.rm=T)/sum(weights,na.rm=T)
  se = sqrt(1/sum(weights,na.rm=T))
  return(c(est,se))
}

gs_results = read.csv('GS/GS_results.csv')
dimnames(gs_results)[[1]] = gs_results[,1]
botnia_results = read.csv('botnia/botnia_results.csv')
#moba_results = read.csv('moba_maths_reading.csv')
fhs_results = read.csv('FHS_SNIPAR_Results.csv')

moba_results = data.frame(phenotype='maths_reading_grade_5',direct=0.189,direct_SE=0.006,
                          parental=0.042,parental_SE=0.017,grandparental=0.062,grandparental_SE=0.017)


# ea_results = rbind(data.frame(study='MoBa',phenotype='Maths and Reading Grade 5',
#                               estimate=0.042,estimate_SE=0.017),
#                    data.frame(study='Generation Scotland',phenotype='Educational attainment (years)',
#                               estimate=gs_results['EA','parental'],estimate_SE=gs_results['EA','parental_SE'])
#                    )
# 
# library(metafor)
# ea_meta = rma(yi=ea_results$estimate,sei=ea_results$estimate_SE,method='EE',slab=ea_results$study)
# forest(ea_meta,xlab='Standardized effect',digits=3)


phenotype_names = data.frame(cohort=c('botnia','fhs','gs'),
                        blood_glucose=c('glucose','BG','Glucose'),
                        HDL=c('HDL','HDL','HDL'),
                        non_HDL=c('non_HDL','NoHDL','Non_HDL'),
                        height=c('height','HGT','height'),
                        BMI=c('BMI','BMI','BMI'),
                        DBP=c('DBP','DBP','DBP'),
                        SBP=c('SBP','SBP','SBP'),
                        FEV1=c(NA,'FEV1','FEV1'),
                        ever_smoker=c(NA,'EVSMK','ever.smoked'),
                        cigarettes_per_day=c(NA,'CPD.LT','cigarettes.per.day'),
                        cognitive_ability=c(NA,NA,'cog'),
                        vocabulary=c(NA,NA,'vocab'),
                        neuroticism=c(NA,NA,'neuroticism'))
                  
phenotypes = dimnames(phenotype_names)[[2]][-1]      
meta_results = data.frame(phenotypes=phenotypes,direct=rep(NA,length(phenotypes)),
                          direct_SE=rep(NA,length(phenotypes)),
                          parental=rep(NA,length(phenotypes)),
                          parental_SE=rep(NA,length(phenotypes)),
                          grandparental=rep(NA,length(phenotypes)),
                          grandparental_SE=rep(NA,length(phenotypes)))

for (i in 1:length(phenotypes)){
  cohort_names = phenotype_names[,i+1]
  estimates = matrix(NA,nrow=3,ncol=3)
  estimate_ses = matrix(NA,nrow=3,ncol=3)
  if (!is.na(cohort_names[1])){
    estimates[1,]=unlist(botnia_results[botnia_results[,1]==cohort_names[1],c('direct','parental','grandparental')])
    estimate_ses[1,]=unlist(botnia_results[botnia_results[,1]==cohort_names[1],c('direct_SE','parental_SE','grandparental_SE')])
  }
  if (!is.na(cohort_names[2])){
    estimates[2,1] = fhs_results[fhs_results$PHENO==cohort_names[2] & fhs_results$EFFECT_TYPE=='proband','BETA']
    estimates[2,2] = fhs_results[fhs_results$PHENO==cohort_names[2] & fhs_results$EFFECT_TYPE=='parental','BETA']
    estimate_ses[2,1] = fhs_results[fhs_results$PHENO==cohort_names[2] & fhs_results$EFFECT_TYPE=='proband','SE']
    estimate_ses[2,2] = fhs_results[fhs_results$PHENO==cohort_names[2] & fhs_results$EFFECT_TYPE=='parental','SE']
  }
  if (!is.na(cohort_names[3])){
    estimates[3,]=unlist(gs_results[gs_results[,1]==cohort_names[3],c('direct','parental','grandparental')])
    estimate_ses[3,]=unlist(gs_results[gs_results[,1]==cohort_names[3],c('direct_SE','parental_SE','grandparental_SE')])
  }
  # Meta-analysis
  meta_results[i,c('direct','direct_SE')] = fe_meta(estimates[,1],estimate_ses[,1])
  meta_results[i,c('parental','parental_SE')] = fe_meta(estimates[,2],estimate_ses[,2])
  meta_results[i,c('grandparental','grandparental_SE')] = fe_meta(estimates[,3],estimate_ses[,3])
}

write.csv(meta_results,'meta_results.csv')


plot_ests = data.frame()
plot_ses = data.frame()

for (i in 1:length(phenotypes)){
  plot_ests = rbind(plot_ests,data.frame(phenotype=rep(phenotypes[i],3),
                                         effect=c('direct','parental','grandparental'),
                                         est=as.vector(t(meta_results[i,c('direct','parental','grandparental')]))))
  plot_ses = rbind(plot_ses,data.frame(phenotype=rep(phenotypes[i],3),
                                         effect=c('direct','parental','grandparental'),
                                         est=as.vector(t(meta_results[i,c('direct_SE','parental_SE','grandparental_SE')]))))
}

plot_ests$phenotype = factor(plot_ests$phenotype,
                             levels=meta_results$phenotypes[order(-meta_results$direct^2)])

plot_ests$effect = factor(plot_ests$effect,
                          levels=c('direct','parental','grandparental'))

plot_lower = plot_ests[,3]+qnorm(0.025)*plot_ses[,3]
plot_upper = plot_ests[,3]-qnorm(0.025)*plot_ses[,3]

## Plot
library(ggplot2)
gpar_plot = ggplot(plot_ests,aes(x = phenotype,y = est,fill = effect)) +
  geom_col(position = position_dodge(0.75),width=0.75)+
  geom_errorbar(aes(ymin = plot_lower,ymax=plot_upper),position=position_dodge(0.75),width=0.25, size=0.5)+
  theme_minimal() + theme(axis.line = element_line(color="black"),
                          axis.ticks = element_line(color="black"),
                          panel.border = element_blank(),
                          axis.text.x = element_text(angle = 45, hjust=1))


ggsave('gpar_plot_meta.pdf',plot=gpar_plot,width=10,height=5)




