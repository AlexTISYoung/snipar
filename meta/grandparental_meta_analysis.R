setwd('~/snipar/meta')

# Approximate variance of ratio 
var_ratio_approx = function(x,y,vx,vy,cxy){
  x2 = x^2; y2 = y^2
  return((x2/y2)*(vx/x2-2*cxy/(x*y)+vy/y2))
}

# Read gen models 1-3 function
read_gen_models = function(gen1_effects,gen2_effects,gen2_vcov,gen3_effects,gen3_vcov,sign_flip=FALSE,lme4=FALSE){
  ## Estimates to output 
  parameter_names = c('population','direct','paternal_NTC','maternal_NTC','average_NTC','maternal_minus_paternal','maternal_minus_paternal_direct_ratio',
  'direct_3','paternal','maternal','parental','grandpaternal','grandmaternal','grandparental',
  'parental_direct_ratio','paternal_direct_ratio','maternal_direct_ratio')
  results = matrix(NA,nrow=length(parameter_names),ncol=2)
  dimnames(results)[[1]] = parameter_names
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
  if (lme4){
    results_1gen = read.table(gen1_effects,skip=1,header=F,row.names=1)
  } else {results_1gen = read.table(gen1_effects,row.names=1)}
  results['population',1:2] = as.matrix(results_1gen['proband',1:2])
  ## Get non-transmitted coefficients
  if (lme4){
   results_2gen = read.table(gen2_effects,row.names=1,skip=1,header=F) 
  } else {results_2gen = read.table(gen2_effects,row.names=1)}
  results_2gen = results_2gen[c('proband','paternal','maternal'),1:2]
  results_2gen_vcov = read.table(gen2_vcov,row.names=1)
  results_2gen_vcov = results_2gen_vcov[c('proband','paternal','maternal'),c('proband','paternal','maternal')]
  # Transform
  A_ntc = matrix(0,nrow=5,ncol=3)
  A_ntc[1:3,1:3] = diag(3)
  A_ntc[4,] = c(0,0.5,0.5)
  A_ntc[5,] = c(0,-1,1)
  results_2gen_vcov = A_ntc%*%as.matrix(results_2gen_vcov)%*%t(A_ntc)
  gen2_effects = c('direct','paternal_NTC','maternal_NTC','average_NTC','maternal_minus_paternal')
  dimnames(results_2gen_vcov)[[1]] = gen2_effects
  dimnames(results_2gen_vcov)[[2]] = gen2_effects
  results[gen2_effects,1] = A_ntc%*%as.matrix(results_2gen[,1])
  results[gen2_effects,2] = sqrt(diag(results_2gen_vcov))
  # Maternal minus paternal direct ratio
  results['maternal_minus_paternal_direct_ratio',1] = results['maternal_minus_paternal',1]/results['direct_3',1]
  results['maternal_minus_paternal_direct_ratio',2] = sqrt(var_ratio_approx(results['maternal_minus_paternal',1],results['direct',1],
                                                        results['maternal_minus_paternal',2]^2,results['direct',2]^2,
                                                        results_2gen_vcov['maternal_minus_paternal','direct']))
  ## Get 3 generation results
  if (lme4){
    results_effects = read.table(gen3_effects,skip=1,header=F,row.names=1)
  } else {
    results_effects = read.table(gen3_effects,row.names=1)}
  results_vcov = read.table(gen3_vcov,row.names=1,header=T)
  if ('gp'%in%dimnames(results_effects)[[1]]){gparsum=TRUE} else if ('gpp'%in%dimnames(results_effects)[[1]]){gparsum=FALSE} else {stop()}
  if (gparsum){
    results_effects = results_effects[c('proband','paternal','maternal','gp','gm'),] 
    results_vcov = results_vcov[c('proband','paternal','maternal','gp','gm'),
                                c('proband','paternal','maternal','gp','gm')]
  } else{
    results_effects = results_effects[c('proband','paternal','maternal','gpp','gpm','gmp','gmm'),1:2]
    results_vcov = results_vcov[c('proband','paternal','maternal','gpp','gpm','gmp','gmm'),
                                c('proband','paternal','maternal','gpp','gpm','gmp','gmm')]
  }
  # Transform and save
  if (!gparsum){
    # Transformation matrix for full model
    A_3 = matrix(0,nrow=7,ncol=7)
    A_3[1:3,1:3] = diag(3)
    A_3[4,2:3] = c(0.5,0.5)
    A_3[5:7,4:7] = rbind(c(0.5,0.5,0,0),
                         c(0,0,0.5,0.5),
                         c(0.25,0.25,0.25,0.25))} else{
                           # Transformation matrix for grandparental sum model
                           A_3 = matrix(0,nrow=7,ncol=5)
                           A_3[1:3,1:3] = diag(3)
                           A_3[4,2:3] = c(0.5,0.5)
                           A_3[5,4] = 1
                           A_3[6,5] = 1
                           A_3[7,4:5] = c(0.5,0.5)}
  results_effects = A_3%*%as.matrix(results_effects[,1])
  gen3_effects = c('direct_3','paternal','maternal','parental','grandpaternal','grandmaternal','grandparental')
  results_vcov = A_3%*%as.matrix(results_vcov)%*%t(A_3)
  dimnames(results_vcov)[[1]] = gen3_effects
  dimnames(results_vcov)[[2]] = gen3_effects
  # Compute indirect to direct ratios
  results[gen3_effects,1]=results_effects
  results[c('direct_3','paternal','maternal','parental','grandpaternal','grandmaternal','grandparental'),2]=sqrt(diag(results_vcov))
  # Parental direct ratio
  results['parental_direct_ratio',1] = results['parental',1]/results['direct_3',1]
  results['parental_direct_ratio',2] = sqrt(var_ratio_approx(results['parental',1],results['direct_3',1],
                                                        results['parental',2]^2,results['direct_3',2]^2,
                                                        results_vcov['parental','direct_3']))
  # Paternal direct ratio
  results['paternal_direct_ratio',1] = results['paternal',1]/results['direct_3',1]
  results['paternal_direct_ratio',2] = sqrt(var_ratio_approx(results['paternal',1],results['direct_3',1],
                                                        results['paternal',2]^2,results['direct_3',2]^2,
                                                        results_vcov['paternal','direct_3']))
  # Maternal direct ratio
  results['maternal_direct_ratio',1] = results['maternal',1]/results['direct_3',1]
  results['maternal_direct_ratio',2] = sqrt(var_ratio_approx(results['maternal',1],results['direct_3',1],
                                                        results['maternal',2]^2,results['direct_3',2]^2,
                                                        results_vcov['maternal','direct_3']))                               
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

fe_meta_Z = function(ests,ses){
  weights = 1/ses^2
  ests = ests/ses
  est = sum(weights*ests,na.rm=T)/sqrt(sum(weights^2,na.rm=T))
  se = 1 
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
botnia_dir = '../Botnia/results/'
botnia_traits = read.table('../Botnia/results/traits.txt',header=F)
botnia_traits=c('EA','glucose','HDL','non_HDL','height','BMI','DBP','SBP')
# GS 
gs_dir = '../GS/pgs/'
#gs_traits = read.table('GS/trait_names.txt',header=F)           
gs_traits = c("Glucose","Non_HDL", "HDL","height","FEV1","BMI","ever.smoked",
                "cigarettes.per.day","cog","vocab","SBP","DBP","neuroticism",
                "EA_years","EA_years_mid","EA_quals")
#gs_traits = cbind(1:length(gs_traits),gs_traits)
# FSH
fhs_dir = '../FHS/FHS_GPall_Mods/'
fhs_binary_dir = '../FHS/FHS_Binary_Phenos/'

phenotype_names = data.frame(cohort=c('botnia','fhs','gs','moba','finngen'),
                        educational_attainment=c('EA',NA,'EA_quals', 'Achievement', NA),
                        blood_glucose=c('glucose','BG.All','Glucose',NA,NA),
                        HDL=c('HDL','HDL','HDL',NA,NA),
                        non_HDL=c('non_HDL','NonHDL','Non_HDL',NA,NA),
                        height=c(NA,'HGT','height','Height','height'),
                        BMI=c(NA,'BMI','BMI','BMI','BMI'),
                        DBP=c('DBP','DBP','DBP',NA,NA),
                        SBP=c('SBP','SBP','SBP',NA,'hypertension'),
                        FEV1=c(NA,'FEV1','FEV1',NA,NA),
                        ever_smoker=c(NA,'EVSMK','ever.smoked',NA,'ever_smoker'),
                        cigarettes_per_day=c(NA,'CPD.Max','cigarettes.per.day',NA,NA),
                        cognitive_ability=c(NA,NA,'cog',NA,NA),
                        vocabulary=c(NA,NA,'vocab',NA,NA),
                        age_at_first_birth_women=c(NA,'AFB',NA,NA,'AAFB_WOMEN'),
                        number_of_children_women=c(NA,'NEB',NA,NA,NA),
                        depressive_symptoms=c(NA,'CESD','neuroticism','Depression','depression'),
                        ADHD=c(NA,NA,NA,'ADHD','ADHD'),
                        alcohol_use_disorder=c(NA,NA,NA,NA,'alcohol_use_disorder'))

phenotypes = dimnames(phenotype_names)[[2]][-1]
binary_phens = c('ADHD','ever_smoker','hypertension','alcohol_use_disorder','depression') 

# Record which phenotypes in which cohort
in_cohort = t(phenotype_names)
dimnames(in_cohort)[[2]] = in_cohort[1,]
in_cohort = in_cohort[-1,]
in_cohort = !is.na(in_cohort)

pgs_names = data.frame(cohort=c('botnia','fhs','gs','moba','finngen'),
                              ADHD=c(NA,NA,NA,'ADHD_Demontis_2023','ADHD1'),
                              AAFB=c(NA,'AFB_Repo_pgs',NA,NA,'AFB2'),
                              BMI=c('BMI_UKB','BMI_UKB_pgs','bmi','BMI_GIANT_2018','bmi'),
                              depression=c(NA,'MDD_pgs','depression','PGC_UKB_depresssion','DEP1'),
                              EA4= c('EA4_2.8m','EA4_pgs','EA4_hm3','EA4','EA4'),
                              ever_smoker= c(NA,'EVSMK_UKB_pgs','ever_smoke',NA,'EVERSMOKE2'),
                              externalizing = c(NA,NA,NA,NA,'externalizing'),
                              height=c('height_UKB','HGT_UKB_pgs','height','height_yengo_2022','height'),
                              number_of_children_women=c(NA,'NEB_Repo_pgs',NA,NA,'NEBwomen2'))

pgss = dimnames(pgs_names)[[2]][-1] 

## Meta-analysis results table
# Effects
effect_names =  c('population','direct','paternal_NTC','maternal_NTC','average_NTC','maternal_minus_paternal',
                  'direct_3','paternal','maternal','parental','grandpaternal','grandmaternal','grandparental',
                  'parental_direct_ratio','paternal_direct_ratio','maternal_direct_ratio','maternal_minus_paternal_direct_ratio')

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
finngen_results = meta_results[,in_cohort[,'finngen'],]

for (j in 1:length(pgss)){
  pgs_names_j = pgs_names[,j+1]
  print(paste('PGS',pgss[j]))
  for (i in 1:length(phenotypes)){
    print(phenotypes[i])
    # Get PGS name for cohort
    # Get names for each cohort
    cohort_names = phenotype_names[,i+1]
    # Estimates
    estimates = matrix(NA,nrow=length(effect_names),ncol=5)
    dimnames(estimates)[[1]] = effect_names
    dimnames(estimates)[[2]] = phenotype_names[,1]
    # Standard errors
    estimate_ses = matrix(NA,nrow=length(effect_names),ncol=5)
    dimnames(estimate_ses) = dimnames(estimates)
    ## Read results from each cohort and transform
    # # Botnia
    if (!is.na(cohort_names[1]) & !is.na(pgs_names_j[1])){
      # Find trait index
      trait_index = trait_index = match(cohort_names[1],botnia_traits)
      # Read estimates and transform
      botnia_estimates = read_gen_models(paste(botnia_dir,pgs_names_j[1],'/',trait_index,'.1.effects.txt',sep=''),
                                        paste(botnia_dir,pgs_names_j[1],'/',trait_index,'.2.effects.txt',sep=''),
                                        paste(botnia_dir,pgs_names_j[1],'/',trait_index,'.2.vcov.txt',sep=''),
                                        paste(botnia_dir,pgs_names_j[1],'/',trait_index,'.3.effects.txt',sep=''),
                                        paste(botnia_dir,pgs_names_j[1],'/',trait_index,'.3.vcov.txt',sep=''))
      # Store for meta-analysis
      estimates[,1] = botnia_estimates[,1]
      estimate_ses[,1] = botnia_estimates[,2]
      # Store in botnia results
      botnia_results[j,match(phenotypes[i],dimnames(botnia_results)$phenotype),effect_names] = botnia_estimates[,1]
      botnia_results[j,match(phenotypes[i],dimnames(botnia_results)$phenotype),paste(effect_names,'SE',sep='_')] = botnia_estimates[,2]
    }
    # FHS
    if (!is.na(cohort_names[2]) & !is.na(pgs_names_j[2])){
      # Read estimates and transform
      if (phenotypes[i]%in%binary_phens){
        fhs_estimates = read_gen_models(paste(fhs_binary_dir,pgs_names_j[2],'/',cohort_names[2],'.1.coefficients.txt',sep=''),
                                        paste(fhs_binary_dir,pgs_names_j[2],'/',cohort_names[2],'.2.coefficients.txt',sep=''),
                                        paste(fhs_binary_dir,pgs_names_j[2],'/',cohort_names[2],'.2.vcov.txt',sep=''),
                                        paste(fhs_binary_dir,pgs_names_j[2],'/',cohort_names[2],'.3.coefficients.txt',sep=''),
                                        paste(fhs_binary_dir,pgs_names_j[2],'/',cohort_names[2],'.3.vcov.txt',sep=''),lme4=TRUE)
      } else {
        fhs_estimates = read_gen_models(paste(fhs_dir,pgs_names_j[2],'/',cohort_names[2],'.1.effects.txt',sep=''),
                                        paste(fhs_dir,pgs_names_j[2],'/',cohort_names[2],'.2.effects.txt',sep=''),
                                        paste(fhs_dir,pgs_names_j[2],'/',cohort_names[2],'.2.vcov.txt',sep=''),
                                        paste(fhs_dir,pgs_names_j[2],'/',cohort_names[2],'.3.effects.txt',sep=''),
                                        paste(fhs_dir,pgs_names_j[2],'/',cohort_names[2],'.3.vcov.txt',sep=''))}
      # Store for meta-analysis
      estimates[,2] = fhs_estimates[,1]
      estimate_ses[,2] = fhs_estimates[,2]
      # Store in fhs results
      fhs_results[j,match(phenotypes[i],dimnames(fhs_results)$phenotype),effect_names] = fhs_estimates[,1]
      fhs_results[j,match(phenotypes[i],dimnames(fhs_results)$phenotype),paste(effect_names,'SE',sep='_')] = fhs_estimates[,2]
    }
    # GS
    if (!is.na(cohort_names[3]) & !is.na(pgs_names_j[3])){
      # Find trait index
      trait_index = match(cohort_names[3],gs_traits)
      # Read estimates and transform
      if (phenotypes[i]%in%binary_phens){
        gs_estimates = read_gen_models(paste0(gs_dir,pgs_names_j[3],'/',cohort_names[3],'.1.coefficients.txt'),
                                        paste0(gs_dir,pgs_names_j[3],'/',cohort_names[3],'.2.coefficients.txt'),
                                        paste0(gs_dir,pgs_names_j[3],'/',cohort_names[3],'.2.vcov.txt'),
                                        paste0(gs_dir,pgs_names_j[3],'/',cohort_names[3],'.3.coefficients.txt'),
                                        paste0(gs_dir,pgs_names_j[3],'/',cohort_names[3],'.3.vcov.txt'),lme4=TRUE)
      } else {
        gs_estimates = read_gen_models(paste0(gs_dir,pgs_names_j[3],'/',trait_index,'.1.effects.txt'),
                                        paste0(gs_dir,pgs_names_j[3],'/',trait_index,'.2.effects.txt'),
                                        paste0(gs_dir,pgs_names_j[3],'/',trait_index,'.2.vcov.txt'),
                                        paste0(gs_dir,pgs_names_j[3],'/',trait_index,'.3.effects.txt'),
                                        paste0(gs_dir,pgs_names_j[3],'/',trait_index,'.3.vcov.txt'))}
      # Store for meta-analysis
      estimates[,3] = gs_estimates[,1]
      estimate_ses[,3] = gs_estimates[,2]
      # Store in GS results
      gs_results[j,match(phenotypes[i],dimnames(gs_results)$phenotype),effect_names] = gs_estimates[,1]
      gs_results[j,match(phenotypes[i],dimnames(gs_results)$phenotype),paste(effect_names,'SE',sep='_')] = gs_estimates[,2]
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
      moba_results[j,match(phenotypes[i],dimnames(moba_results)$phenotype),effect_names] = moba_estimates[,1]
      moba_results[j,match(phenotypes[i],dimnames(moba_results)$phenotype),paste(effect_names,'SE',sep='_')] = moba_estimates[,2]
    }
    # Finngen
    if (!is.na(cohort_names[5]) & !is.na(pgs_names_j[5])){
      # Find trait index
      trait_index = match(cohort_names[5],finngen_traits)
      # Read estimates and transform
      if (pgs_names_j[5]=='EA4'){sign_flip=TRUE} else {sign_flip=FALSE}
      if (phenotypes[i]%in%binary_phens & !phenotypes[i]=='ADHD'){
        finngen_estimates = read_gen_models(paste0(finngen_dir,pgs_names_j[5],'/',cohort_names[5],'.1.coefficients.txt'),
                                        paste0(finngen_dir,pgs_names_j[5],'/',cohort_names[5],'.2.coefficients.txt'),
                                        paste0(finngen_dir,pgs_names_j[5],'/',cohort_names[5],'.2.vcov.txt'),
                                        paste0(finngen_dir,pgs_names_j[5],'/',cohort_names[5],'.3.coefficients.txt'),
                                        paste0(finngen_dir,pgs_names_j[5],'/',cohort_names[5],'.3.vcov.txt'),sign_flip=sign_flip,lme4=TRUE)
      } else {
        finngen_estimates = read_gen_models(paste0(finngen_dir,pgs_names_j[5],'/',trait_index,'.1.effects.txt'),
                                        paste0(finngen_dir,pgs_names_j[5],'/',trait_index,'.2.effects.txt'),
                                        paste0(finngen_dir,pgs_names_j[5],'/',trait_index,'.2.vcov.txt'),
                                        paste0(finngen_dir,pgs_names_j[5],'/',trait_index,'.3.effects.txt'),
                                        paste0(finngen_dir,pgs_names_j[5],'/',trait_index,'.3.vcov.txt'),sign_flip=sign_flip)}
      # Store for meta-analysis
      estimates[,5] = finngen_estimates[,1]
      estimate_ses[,5] = finngen_estimates[,2]
      # Store in GS results
      finngen_results[j,match(phenotypes[i],dimnames(finngen_results)$phenotype),effect_names] = finngen_estimates[,1]
      finngen_results[j,match(phenotypes[i],dimnames(finngen_results)$phenotype),paste(effect_names,'SE',sep='_')] = finngen_estimates[,2]
    }
    # Meta-analysis
    for (effect_j in effect_names){
      if (effect_j%in%c('parental_direct_ratio','paternal_direct_ratio','maternal_direct_ratio','maternal_minus_paternal_direct_ratio')){
        meta_results[j,i,c(effect_j,paste(effect_j,'SE',sep='_'))] = fe_meta(estimates[effect_j,],estimate_ses[effect_j,])
      } else {
        meta_results[j,i,c(effect_j,paste(effect_j,'SE',sep='_'))] = fe_meta_Z(estimates[effect_j,],estimate_ses[effect_j,])
      }
    }
  }
}

# Calculate P-values
for (effect_name in effect_names){
  meta_results[,,paste(effect_name,'log10P',sep='_')] = -log10(exp(1))*pchisq((meta_results[,,effect_name]/meta_results[,,paste(effect_name,'SE',sep='_')])^2,1,lower.tail=F,log.p=T)
  botnia_results[,,paste(effect_name,'log10P',sep='_')] = -log10(exp(1))*pchisq((botnia_results[,,effect_name]/botnia_results[,,paste(effect_name,'SE',sep='_')])^2,1,lower.tail=F,log.p=T)
  fhs_results[,,paste(effect_name,'log10P',sep='_')] = -log10(exp(1))*pchisq((fhs_results[,,effect_name]/fhs_results[,,paste(effect_name,'SE',sep='_')])^2,1,lower.tail=F,log.p=T)
  gs_results[,,paste(effect_name,'log10P',sep='_')] = -log10(exp(1))*pchisq((gs_results[,,effect_name]/gs_results[,,paste(effect_name,'SE',sep='_')])^2,1,lower.tail=F,log.p=T)
  moba_results[,,paste(effect_name,'log10P',sep='_')] = -log10(exp(1))*pchisq((moba_results[,,effect_name]/moba_results[,,paste(effect_name,'SE',sep='_')])^2,1,lower.tail=F,log.p=T)
  finngen_results[,,paste(effect_name,'log10P',sep='_')] = -log10(exp(1))*pchisq((finngen_results[,,effect_name]/finngen_results[,,paste(effect_name,'SE',sep='_')])^2,1,lower.tail=F,log.p=T)
}

meta_results[is.infinite(meta_results)] = NA

## Plot
# EA plot
require(ggplot2)
ea_plot = rbind(cbind(phenotype='EA',cohort='GS',
                value=gs_results['EA4','educational_attainment','parental'],
                SE=gs_results['EA4','educational_attainment','parental_SE']),
               cbind(phenotype='EA',cohort='Botnia',
                value=botnia_results['EA4','educational_attainment','parental'],
                SE=botnia_results['EA4','educational_attainment','parental_SE']),
              cbind(phenotype='math_reading_achievement_age_10',cohort='MoBa',
                value=moba_results['EA4','educational_attainment','parental'],
                SE=moba_results['EA4','educational_attainment','parental_SE'])
)
ea_plot = data.frame(ea_plot)
ea_plot$pgs = 'EA4'
# AFB plot
afb_plot = rbind(cbind(phenotype='AAFB (women)',cohort='FHS',
                value=fhs_results['AAFB','age_at_first_birth_women','parental'],
                SE=fhs_results['AAFB','age_at_first_birth_women','parental_SE']),
               cbind(phenotype='AAFB (women)',cohort='Finngen',
                value=finngen_results['AAFB','age_at_first_birth_women','parental'],
                SE=finngen_results['AAFB','age_at_first_birth_women','parental_SE'])
)
afb_plot = data.frame(afb_plot)
afb_plot$pgs = 'AAFB'

plot_df = rbind(ea_plot,afb_plot)
plot_df$value=as.numeric(plot_df$value)
plot_df$SE=as.numeric(plot_df$SE)
plot_df$lower = plot_df$value+qnorm(0.025)*plot_df$SE
plot_df$upper = plot_df$value-qnorm(0.025)*plot_df$SE

gplot = ggplot(plot_df, aes(fill=phenotype, y=value, x=cohort)) +
    geom_errorbar(aes(ymin = plot_df$lower, ymax = plot_df$upper))+
    geom_point(position="dodge", stat="identity")+coord_flip()+theme_minimal() + theme(axis.line = element_line(color="black"),
                          axis.ticks = element_line(color="black"),
                          panel.border = element_blank(),
                          axis.text.x = element_text(angle = 45, hjust=1))+
                          geom_hline(yintercept=0,linetype='dashed',width=2)+facet_wrap(~pgs)



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

