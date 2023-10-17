setwd('~/snipar/meta')

# Read gen models 1-3 function
read_gen_models = function(gen1_effects,gen2_effects,gen2_vcov,gen3_effects,gen3_vcov,sign_flip=FALSE){
  ## Estimates to output 
  results = matrix(NA,nrow=11,ncol=2)
  dimnames(results)[[1]] = c('population','direct','paternal_NTC','maternal_NTC','average_NTC','paternal','maternal','parental','grandpaternal','grandmaternal','grandparental')
  dimnames(results)[[2]] = c('estimates','SE')
  ## Check files exist
  if (!file.exists(gen1_effects)){
    print(paste(gen1_effects,'does not exist'))
    return(results)
  } else if (!file.exists(gen2_effects)){
    print(paste(gen2_effects,'does not exist'))
    return(results)
  } else if (!file.exists(gen2_vcov)){
    print(paste(gen2_vcov,'does not exist'))
    return(results)
  } else if (!file.exists(gen3_effects)){
    print(paste(gen3_effects,'does not exist'))
    return(results)
  } else if (!file.exists(gen3_vcov)){
    print(paste(gen3_vcov,'does not exist'))
    return(results)
  } else {
  ## Get population effect
  results_1gen = read.table(gen1_effects,row.names=1)
  results['population',1:2] = as.matrix(results_1gen['proband',])
  ## Get non-transmitted coefficients
  results_2gen = read.table(gen2_effects,row.names=1)
  results_2gen = results_2gen[c('proband','paternal','maternal'),1:2]
  results_2gen_vcov = read.table(gen2_vcov,row.names=1)
  results_2gen_vcov = results_2gen_vcov[c('proband','paternal','maternal'),c('proband','paternal','maternal')]
  # Transform
  A_ntc = matrix(0,nrow=4,ncol=3)
  A_ntc[1:3,1:3] = diag(3)
  A_ntc[4,] = c(0,0.5,0.5)
  results_2gen_vcov = A_ntc%*%as.matrix(results_2gen_vcov)%*%t(A_ntc)
  results[c('direct','paternal_NTC','maternal_NTC','average_NTC'),1] = A_ntc%*%as.matrix(results_2gen[,1])
  results[c('direct','paternal_NTC','maternal_NTC','average_NTC'),2] = sqrt(diag(results_2gen_vcov))
  ## Get 3 generation results
  results_effects = read.table(gen3_effects,row.names=1)
  results_vcov = read.table(gen3_vcov,row.names=1,header=T)
  if ('gp'%in%dimnames(results_effects)[[1]]){gparsum=TRUE} else if ('gpp'%in%dimnames(results_effects)[[1]]){gparsum=FALSE} else {stop()}
  if (gparsum){
    results_effects = results_effects[c('paternal','maternal','gp','gm'),] 
    results_vcov = results_vcov[c('paternal','maternal','gp','gm'),
                                c('paternal','maternal','gp','gm')]
  } else{
    results_effects = results_effects[c('paternal','maternal','gpp','gpm','gmp','gmm'),]
    results_vcov = results_vcov[c('paternal','maternal','gpp','gpm','gmp','gmm'),
                                c('paternal','maternal','gpp','gpm','gmp','gmm')]
  }
  # Transform and save
  if (!gparsum){
    # Transformation matrix for full model
    A_3 = matrix(0,nrow=6,ncol=6)
    A_3[1:2,1:2] = diag(2)
    A_3[3,1:2] = c(0.5,0.5)
    A_3[4:6,3:6] = rbind(c(0.5,0.5,0,0),
                         c(0,0,0.5,0.5),
                         c(0.25,0.25,0.25,0.25))} else{
                           # Transformation matrix for grandparental sum model
                           A_3 = matrix(0,nrow=6,ncol=4)
                           A_3[1:2,1:2] = diag(2)
                           A_3[3,1:2] = c(0.5,0.5)
                           A_3[4,3] = 1
                           A_3[5,4] = 1
                           A_3[6,3:4] = c(0.5,0.5)}
  results_effects = A_3%*%as.matrix(results_effects[,1])
  results_vcov = A_3%*%as.matrix(results_vcov)%*%t(A_3)
  results[c('paternal','maternal','parental','grandpaternal','grandmaternal','grandparental'),1]=results_effects
  results[c('paternal','maternal','parental','grandpaternal','grandmaternal','grandparental'),2]=sqrt(diag(results_vcov))
  # Return
  if (sign_flip){results[,1] = -results[,1]}
  return(results)}
}

# Function for fixed effects meta-analysis
fe_meta = function(ests,ses){
  weights = 1/ses^2
  est = sum(weights*ests,na.rm=T)/sum(weights,na.rm=T)
  se = sqrt(1/sum(weights,na.rm=T))
  return(c(est,se))
}

# MoBa
moba_dir = '../moba/MoBa_3gen_results_050923'
moba_traits = c('Achievement','Height','BMI','ADHD','Depression')

# Finngen
finngen_dir = '../finngen/pgs_out/'
finngen_traits = c('height','NC','BMI','ever_smoker',
                'depression','AAFB','AAFB_WOMEN',
                'AAFB_MEN','NC_WOMEN','NC_MEN',
                'ADHD','asthma','eczema','hypertension',
                'alcohol_use_disorder','allergic_rhinitis',
                'migraine','copd')
# Botnia
#botnia_dir = '../Botnia/'
#botnia_traits = read.table('../Botnia/traits.txt',header=F)
# GS 
gs_dir = '../GS/pgs/'
#gs_traits = read.table('GS/trait_names.txt',header=F)           
gs_traits = c("Glucose","Non_HDL", "HDL","height","FEV1","BMI","ever.smoked",
                "cigarettes.per.day","cog","vocab","SBP","DBP","neuroticism",
                "EA_years","EA_years_mid","EA_quals")
#gs_traits = cbind(1:length(gs_traits),gs_traits)
# FSH
fhs_csv = read.csv('../FHS/FHS_Update/fhs_GPall_effects.csv')
fhs_csv_vcov = read.csv('../FHS/FHS_Update/fhs_GPall_vcov.csv')

phenotype_names = data.frame(cohort=c('botnia','fhs','gs','moba','finngen'),
                        educational_attainment=c('EA',NA,'EA_quals', NA, NA),
                        math_reading_attainment_age_10=c(NA,NA,NA,'Achievement',NA),
                        blood_glucose=c('glucose','BG.All','Glucose',NA,NA),
                        HDL=c('HDL','HDL','HDL',NA,NA),
                        non_HDL=c('non_HDL','NonHDL','Non_HDL',NA,NA),
                        height_adult=c('height','HGT','height',NA,'height'),
                        height_age_8=c(NA,NA,NA,'Height',NA),
                        BMI_adult=c('BMI','BMI','BMI',NA,'BMI'),
                        BMI_age_8=c(NA,NA,NA,'BMI',NA),
                        DBP=c('DBP','DBP','DBP',NA,NA),
                        SBP=c('SBP','SBP','SBP',NA,NA),
                        FEV1=c(NA,'FEV1','FEV1',NA,NA),
                        ever_smoker=c(NA,'EVSMK','ever.smoked',NA,'ever_smoker'),
                        cigarettes_per_day=c(NA,'CPD.Cur','cigarettes.per.day',NA,NA),
                        cognitive_ability=c(NA,NA,'cog',NA,NA),
                        vocabulary=c(NA,NA,'vocab',NA,NA),
                        number_of_children_women=c(NA,'NEB',NA,NA,'NC_WOMEN'),
                        number_of_children_men=c(NA,NA,NA,NA,'NC_MEN'),
                        age_at_first_birth_women=c(NA,"AFB",NA,NA,'AAFB_WOMEN'),
                        depression=c(NA,'MDD',NA,NA,'depression'),
                        depressive_symptoms=c(NA,'CESD','neuroticism','Depression',NA),
                        ADHD=c(NA,NA,NA,'ADHD','ADHD'),
                        hypertension=c(NA,NA,NA,NA,'hypertension'),
                        alcohol_use_disorder=c(NA,NA,NA,NA,'alcohol_use_disorder'))

phenotypes = dimnames(phenotype_names)[[2]][-1] 

# Record which phenotypes in which cohort
in_cohort = t(phenotype_names)
dimnames(in_cohort)[[2]] = in_cohort[1,]
in_cohort = in_cohort[-1,]
in_cohort = !is.na(in_cohort)

pgs_names = data.frame(cohort=c('botnia','fhs','gs','moba','finngen'),
                              ADHD=c(NA,"ADHD_Demontis",NA,'ADHD_Demontis_2023','ADHD1'),
                              AAFB=c(NA,'AFB_Repo',NA,NA,'AFB2'),
                              BMI=c('BMI','BMI_UKB','bmi','BMI_GIANT_2018','bmi'),
                              depression=c(NA,'MDD_Howard','depression','PGC_UKB_depresssion','DEP1'),
                              EA4= c('EA4','EA_Okbay','EA4_hm3','EA4','EA4'),
                              ever_smoker= c('ever_smoker','EVSMK_UKB','ever_smoke',NA,'EVERSMOKE2'),
                              externalizing = c(NA,NA,NA,NA,'externalizing'),
                              height=c('height','HGT_UKB','height','height_yengo_2022','height'),
                              number_of_children_women=c(NA,'NEB_UKB',NA,NA,'NEB2'))

pgss = dimnames(pgs_names)[[2]][-1] 

## Meta-analysis results table
# Effects
effect_names =  c('population','direct','paternal_NTC','maternal_NTC','average_NTC',
                  'paternal','maternal','parental','grandpaternal','grandmaternal','grandparental')

meta_results = array(NA,dim=c(length(pgss),length(phenotypes),length(effect_names)*3))
meta_colnames = c()
for (j in 1:length(effect_names)){
  meta_colnames=c(meta_colnames,
                  c(effect_names[j],
                    paste(effect_names[j],'SE',sep='_'),
                    paste(effect_names[j],'log10P',sep='_')))
}
dimnames(meta_results)= list(pgs=pgss,phenotype=phenotypes,statistic=meta_colnames)


## Cohort specific results tables
botnia_results = meta_results[,in_cohort[,'botnia'],]
fhs_results = meta_results[,in_cohort[,'fhs'],]
gs_results = meta_results[,in_cohort[,'gs'],]
moba_results = meta_results[,in_cohort[,'moba'],]

# Transformation matrix for FHS
A_full = matrix(0,nrow=7,ncol=7)
A_full[1:3,1:3] = diag(3)
A_full[4,1:3] = c(0,0.5,0.5)
A_full[5:7,4:7] = rbind(c(0.5,0.5,0,0),
                     c(0,0,0.5,0.5),
                     c(0.25,0.25,0.25,0.25))

for (j in 1:length(pgss)){
  pgs_names_j = pgs_names[,j+1]
  print(paste('PGS',pgss[j]))
  for (i in 1:length(phenotypes)){
    print(phenotypes[i])
    # Get PGS name for cohort
    # Get names for each cohort
    cohort_names = phenotype_names[,i+1]
    # Estimates
    estimates = matrix(NA,nrow=11,ncol=5)
    dimnames(estimates)[[1]] = effect_names
    dimnames(estimates)[[2]] = phenotype_names[,1]
    # Standard errors
    estimate_ses = matrix(NA,nrow=11,ncol=5)
    dimnames(estimate_ses) = dimnames(estimates)
    ## Read results from each cohort and transform
    # # Botnia
    # if (!is.na(cohort_names[1])){
    #   # Find trait index
    #   trait_index = botnia_traits[match(cohort_names[1],botnia_traits[,2]),1]
    #   # Read estimates and transform
    #   botnia_estimates = read_gen_models(paste(botnia_dir,trait_index,'.1.effects.txt',sep=''),
    #                                     paste(botnia_dir,trait_index,'.2.effects.txt',sep=''),
    #                                     paste(botnia_dir,trait_index,'.2.vcov.txt',sep=''),
    #                                     paste(botnia_dir,trait_index,'.3.effects.txt',sep=''),
    #                                     paste(botnia_dir,trait_index,'.3.vcov.txt',sep=''))
    #   # Store for meta-analysis
    #   if (dimnames(phenotype_names)[[2]][i+1]!="educational_attainment"){
    #   estimates[,1] = botnia_estimates[,1]
    #   estimate_ses[,1] = botnia_estimates[,2]}
    #   # Store in botnia results
    #   botnia_results[match(phenotypes[i],botnia_results[,1]),effect_names] = botnia_estimates[,1]
    #   botnia_results[match(phenotypes[i],botnia_results[,1]),paste(effect_names,'SE',sep='_')] = botnia_estimates[,2]
    # }
    # FHS
    if (!is.na(cohort_names[2]) & !is.na(pgs_names_j[2])){
      # Effects
      fhs_effects = fhs_csv[fhs_csv$PHENO==cohort_names[2] & fhs_csv$PGS==pgs_names_j[2],]
      fhs_vcov = fhs_csv_vcov[fhs_csv_vcov$PHENO==cohort_names[2] & fhs_csv_vcov$PGS==pgs_names_j[2],]
      # Get population effect
      estimates['population',2] = as.matrix(fhs_effects[fhs_effects$GEN_MOD==1 & fhs_effects$VAR=='proband',c('BETA')])
      estimate_ses['population',2] = as.matrix(fhs_effects[fhs_effects$GEN_MOD==1 & fhs_effects$VAR=='proband',c('SE')])
      # Get 2 generation effects
      fhs_2gen = fhs_effects[fhs_effects$GEN_MOD==2,]
      dimnames(fhs_2gen)[[1]] = fhs_2gen$VAR
      fhs_2gen_vcov = fhs_vcov[fhs_vcov$GEN_MOD==2,]
      dimnames(fhs_2gen_vcov)[[1]] = fhs_2gen_vcov$VAR
      # Transform
      A_ntc = matrix(0,nrow=3,ncol=2)
      A_ntc[1:2,1:2] = diag(2)
      A_ntc[3,] = c(0.5,0.5)
      fhs_2gen_vcov = A_ntc%*%as.matrix(fhs_2gen_vcov[c('paternal','maternal'),c('paternal','maternal')])%*%t(A_ntc)
      estimates[c('paternal_NTC','maternal_NTC','average_NTC'),2] = A_ntc%*%as.matrix(fhs_2gen[c('paternal','maternal'),'BETA'])
      estimate_ses[c('paternal_NTC','maternal_NTC','average_NTC'),2] = sqrt(diag(fhs_2gen_vcov))
      # Get 3 generation effects
      fhs_effects = fhs_effects[fhs_effects$GEN_MOD==3,]
      dimnames(fhs_effects)[[1]] = fhs_effects$VAR
      fhs_effects=fhs_effects[c('proband','paternal','maternal','gpp','gpm','gmp','gmm'),'BETA']
      # vcov
      fhs_vcov = fhs_vcov[fhs_vcov$GEN_MOD==3,]
      dimnames(fhs_vcov)[[1]] = fhs_vcov$VAR
      fhs_vcov=fhs_vcov[c('proband','paternal','maternal','gpp','gpm','gmp','gmm'),
                        c('proband','paternal','maternal','gpp','gpm','gmp','gmm')]
      # Transform
      fhs_effects = A_full%*%as.matrix(fhs_effects)
      fhs_vcov = A_full%*%as.matrix(fhs_vcov)%*%t(A_full)
      # Save
      estimates[c('direct','paternal','maternal','parental','grandpaternal','grandmaternal','grandparental'),2] = fhs_effects
      estimate_ses[c('direct','paternal','maternal','parental','grandpaternal','grandmaternal','grandparental'),2] = sqrt(diag(fhs_vcov))
      # Store in fhs results
      #fhs_results[match(phenotypes[i],fhs_results[,1]),effect_names] = estimates[,2]
      #fhs_results[match(phenotypes[i],fhs_results[,1]),paste(effect_names,'SE',sep='_')] =  estimate_ses[,2] 
    }
    # GS
    if (!is.na(cohort_names[3]) & !is.na(pgs_names_j[3])){
      # Find trait index
      trait_index = match(cohort_names[3],gs_traits)
      # Read estimates and transform
      gs_estimates = read_gen_models(paste0(gs_dir,pgs_names_j[3],'/',trait_index,'.1.effects.txt'),
                                        paste0(gs_dir,pgs_names_j[3],'/',trait_index,'.2.effects.txt'),
                                        paste0(gs_dir,pgs_names_j[3],'/',trait_index,'.2.vcov.txt'),
                                        paste0(gs_dir,pgs_names_j[3],'/',trait_index,'.3.effects.txt'),
                                        paste0(gs_dir,pgs_names_j[3],'/',trait_index,'.3.vcov.txt'))
      # Store for meta-analysis
      estimates[,3] = gs_estimates[,1]
      estimate_ses[,3] = gs_estimates[,2]
      # Store in GS results
      #gs_results[match(phenotypes[i],gs_results[,1]),effect_names] = gs_estimates[,1]
      #gs_results[match(phenotypes[i],gs_results[,1]),paste(effect_names,'SE',sep='_')] = gs_estimates[,2]
    }
    # MoBa
    if (!is.na(cohort_names[4]) & !is.na(pgs_names_j[4])){
      trait_index = match(cohort_names[4],moba_traits)
      moba_estimates = read_gen_models(paste0(moba_dir,'/',pgs_names_j[4],'_',trait_index,'.1.effects.txt',sep=''),
                                    paste0(moba_dir,'/',pgs_names_j[4],'_',trait_index,'.2.effects.txt',sep=''),
                                    paste0(moba_dir,'/',pgs_names_j[4],'_',trait_index,'.2.vcov.txt',sep=''),
                                    paste0(moba_dir,'/',pgs_names_j[4],'_',trait_index,'.3.effects.txt',sep=''),
                                    paste0(moba_dir,'/',pgs_names_j[4],'_',trait_index,'.3.vcov.txt',sep=''))
      # Store for meta-analysis
      estimates[,4] = moba_estimates[,1]
      estimate_ses[,4] = moba_estimates[,2]
      # Store in MoBa results
      #moba_results[match(phenotypes[i],moba_results[,1]),effect_names] = moba_estimates[,1]
      #moba_results[match(phenotypes[i],moba_results[,1]),paste(effect_names,'SE',sep='_')] = moba_estimates[,2]
    }
    # Finngen
    if (!is.na(cohort_names[5]) & !is.na(pgs_names_j[5])){
      # Find trait index
      trait_index = match(cohort_names[5],finngen_traits)
      # Read estimates and transform
      if (pgs_names_j[5]=='EA4'){sign_flip=TRUE} else {sign_flip=FALSE}
      finngen_estimates = read_gen_models(paste0(finngen_dir,pgs_names_j[5],'/',trait_index,'.1.effects.txt'),
                                        paste0(finngen_dir,pgs_names_j[5],'/',trait_index,'.2.effects.txt'),
                                        paste0(finngen_dir,pgs_names_j[5],'/',trait_index,'.2.vcov.txt'),
                                        paste0(finngen_dir,pgs_names_j[5],'/',trait_index,'.3.effects.txt'),
                                        paste0(finngen_dir,pgs_names_j[5],'/',trait_index,'.3.vcov.txt'),sign_flip=sign_flip)
      # Store for meta-analysis
      estimates[,5] = finngen_estimates[,1]
      estimate_ses[,5] = finngen_estimates[,2]
      # Store in GS results
      #gs_results[match(phenotypes[i],gs_results[,1]),effect_names] = gs_estimates[,1]
      #gs_results[match(phenotypes[i],gs_results[,1]),paste(effect_names,'SE',sep='_')] = gs_estimates[,2]
    }
    # Meta-analysis
    for (effect_j in effect_names){
      meta_results[j,i,c(effect_j,paste(effect_j,'SE',sep='_'))] = fe_meta(estimates[effect_j,],estimate_ses[effect_j,])
    }
  }
}

# Calculate P-values
for (effect_name in effect_names){
  meta_results[,,paste(effect_name,'log10P',sep='_')] = -log10(exp(1))*pchisq((meta_results[,,effect_name]/meta_results[,,paste(effect_name,'SE',sep='_')])^2,1,lower.tail=F,log.p=T)
}

dimnames(meta_results)[[1]] = meta_results$phenotype

meta_results = data.frame(in_cohort[match(dimnames(meta_results)[[1]],dimnames(in_cohort)[[1]]),],
                          meta_results[,-1])

write.csv(meta_results,'meta_results.csv',row.names=T)

# Write cohort specific results
fhs_results = fhs_results[,-c(grep('log10P',dimnames(fhs_results)[[2]]))]
write.csv(fhs_results,'fhs_results.csv',row.names=F)
moba_results = moba_results[,-c(grep('log10P',dimnames(moba_results)[[2]]))]
write.csv(moba_results,'moba_results.csv',row.names=F)
gs_results = gs_results[,-c(grep('log10P',dimnames(gs_results)[[2]]))]
write.csv(gs_results,'gs_results.csv',row.names=F)
botnia_results = botnia_results[,-c(grep('log10P',dimnames(botnia_results)[[2]]))]
write.csv(botnia_results,'botnia_results.csv',row.names=F)


plot_ests = data.frame()
plot_ses = data.frame()

for (i in 1:length(phenotypes)){
  plot_ests = rbind(plot_ests,data.frame(phenotype=rep(phenotypes[i],11),
                                         effect=effect_names,
                                         est=as.vector(t(meta_results[i,effect_names])),
                                         ses=as.vector(t(meta_results[i,paste(effect_names,'SE',sep='_')]))))
#  plot_ses = rbind(plot_ses,data.frame(phenotype=rep(phenotypes[i],11),
#                                         effect=effect_names,
#                                         est=as.vector(t(meta_results[i,paste(effect_names,'SE',sep='_')]))))
}

plot_ests$phenotype = factor(plot_ests$phenotype,
                             levels=dimnames(meta_results)[[1]][order(meta_results$population^2)])

plot_ests$effect = factor(plot_ests$effect,
                          levels=rev(effect_names))

plot_ests$lower = plot_ests$est+qnorm(0.025)*plot_ests$ses
plot_ests$upper = plot_ests$est-qnorm(0.025)*plot_ests$ses

## Plot
library(ggplot2)
require(gridExtra)
require(ggplot2)
require(reshape2)
require(ggrepel)
pop_decomp_effects = c('population','direct','parental','grandparental')
gpar_plot = ggplot(plot_ests[plot_ests$effect%in%pop_decomp_effects,],aes(x = est,y = phenotype,fill = effect)) +
  geom_col(position = position_dodge(0.75),width=0.75)+
  geom_errorbarh(aes(xmin = plot_ests$lower[plot_ests$effect%in%c(pop_decomp_effects)],
                     xmax=plot_ests$upper[plot_ests$effect%in%c(pop_decomp_effects)]),position=position_dodge(0.75),width=0.15, size=0.25)+
  theme_minimal() + theme(axis.line = element_line(color="black"),
                          axis.ticks = element_line(color="black"),
                          panel.border = element_blank(),
                          axis.text.x = element_text(angle = 45, hjust=1))


ggsave('gpar_pop_decomp.pdf',plot=gpar_plot,width=10,height=15)

## Plot indirect only
indirect_effects = c('paternal','maternal','parental','average_NTC','paternal_NTC','maternal_NTC')
gpar_indirect_plot = ggplot(plot_ests[plot_ests$effect%in%indirect_effects & plot_ests$phenotype=='math_reading_attainment_age_10',],
                            aes(x = phenotype,y = est,fill = effect)) +
  geom_col(position = position_dodge(0.75),width=0.75)+
  geom_errorbar(aes(ymin = plot_ests$lower[plot_ests$effect%in%indirect_effects & plot_ests$phenotype=='math_reading_attainment_age_10'],
                     ymax=plot_ests$upper[plot_ests$effect%in%indirect_effects & plot_ests$phenotype=='math_reading_attainment_age_10']),position=position_dodge(0.75),width=0.15, size=0.25)+
  theme_minimal() + theme(axis.line = element_line(color="black"),
                          axis.ticks = element_line(color="black"),
                          panel.border = element_blank(),
                          axis.text.x = element_text(angle = 45, hjust=1))


ggsave('gpar_indirect_plot_meta.pdf',plot=gpar_indirect_plot,width=7,height=5)


#### Detailed MoBa results ####
var_ratio_approx = function(x,y,vx,vy,cxy){
  x2 = x^2; y2 = y^2
  return((x2/y2)*(vx/x2-2*cxy/(x*y)+vy/y2))
}
edu_effects = read.table(paste(moba_dir,'/EA4_3gen_Edu.3.effects.txt',sep=''))
edu_vcov = read.table(paste(moba_dir,'/EA4_3gen_Edu.3.vcov.txt',sep=''))
edu_pop_effect = read.table(paste(moba_dir,'/EA4_3gen_Edu.1.effects.txt',sep=''),row.names=1)

r=0.125

A_3 = matrix(0,nrow=9,ncol=5)
A_3[1:3,1:3] = diag(3)
A_3[4,1:3] = c(0,0.5,0.5)
A_3[5,4] = 1
A_3[6,5] = 1
A_3[7,4:5] = c(0.5,0.5)
A_3[8,1:3] = c(0,-1,1)
A_3[9,1:3] = c(1,(1+r)/2,(1+r)/2)

# Transform
edu_effects = A_3%*%as.matrix(edu_effects[-1,2])
edu_vcov = A_3%*%as.matrix(edu_vcov[-1,-1])%*%t(A_3)
edu_effects = cbind(edu_effects,sqrt(diag(edu_vcov)))

dimnames(edu_effects)[[1]]=c('direct','paternal','maternal','parental','grandpaternal','grandmaternal','grandparental','maternal_minus_paternal','population')
dimnames(edu_effects)[[2]] = c('estimate','SE')
dimnames(edu_vcov)[[1]]=dimnames(edu_effects)[[1]]
dimnames(edu_vcov)[[2]]=dimnames(edu_effects)[[1]]

indirect_direct_ratio = edu_effects['parental',1]/edu_effects['direct',1]
indirect_direct_ratio_var = var_ratio_approx(edu_effects['parental',1],edu_effects['direct',1],
                                             edu_effects['parental',2]^2,edu_effects['direct',2]^2,
                                             edu_vcov['parental','direct'])

print(c(indirect_direct_ratio,sqrt(indirect_direct_ratio_var)))

