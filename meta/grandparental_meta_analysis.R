setwd('meta')
# Approximate variance of ratio 
var_ratio_approx = function(x,y,vx,vy,cxy){
  x2 = x^2; y2 = y^2
  return((x2/y2)*(vx/x2-2*cxy/(x*y)+vy/y2))
}

# gen1_effects = '../GS/pgs/ever_smoke/ever.smoked.1.coefficients.txt'
# gen2_effects = '../GS/pgs/ever_smoke/ever.smoked.2.coefficients.txt'
# gen3_paternal = '../GS/pgs/ever_smoke/ever.smoked.3.paternal.coefficients.txt' 
# gen3_maternal = '../GS/pgs/ever_smoke/ever.smoked.3.maternal.coefficients.txt' 
# gen2_vcov = '../GS/pgs/ever_smoke/ever.smoked.2.vcov.txt'
# gen3_paternal_vcov = '../GS/pgs/ever_smoke/ever.smoked.3.paternal.vcov.txt'  
# gen3_maternal_vcov = '../GS/pgs/ever_smoke/ever.smoked.3.maternal.vcov.txt'  

# Read gen models 1-3 function
read_gen_models = function(gen1_effects,gen2_effects,gen2_vcov,gen3_paternal,gen3_paternal_vcov,gen3_maternal,gen3_maternal_vcov,sign_flip=FALSE,lme4=FALSE){
  ## Estimates to output 
  parameter_names = c('population','direct','paternal_NTC','maternal_NTC','average_NTC','average_NTC_direct_ratio',
  'maternal_minus_paternal','maternal_minus_paternal_direct_ratio',
  'paternal','maternal','grandpaternal','grandmaternal')
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
  } else if (!file.exists(gen3_paternal)){
    print(paste(gen3_paternal,'does not exist'))
    return(results)
  } else if (!file.exists(gen3_maternal)){
    print(paste(gen3_maternal,'does not exist'))
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
  # Average NTC direct ratio
  results['average_NTC_direct_ratio',1] = results['average_NTC',1]/results['direct',1]
  results['average_NTC_direct_ratio',2] = sqrt(var_ratio_approx(results['average_NTC',1],results['direct',1],
                                                        results['average_NTC',2]^2,results['direct',2]^2,
                                                        results_2gen_vcov['average_NTC','direct']))
  # Maternal minus paternal direct ratio
  results['maternal_minus_paternal_direct_ratio',1] = results['maternal_minus_paternal',1]/results['direct',1]
  results['maternal_minus_paternal_direct_ratio',2] = sqrt(var_ratio_approx(results['maternal_minus_paternal',1],results['direct',1],
                                                        results['maternal_minus_paternal',2]^2,results['direct',2]^2,
                                                        results_2gen_vcov['maternal_minus_paternal','direct']))
  ## Get 3 generation results
  # paternal
  if (lme4){
    results_effects = read.table(gen3_paternal,skip=1,header=F,row.names=1)
  } else {
    results_effects = read.table(gen3_paternal,row.names=1)}
  results_vcov = read.table(gen3_paternal_vcov,row.names=1,header=T)
  if ('gp'%in%dimnames(results_effects)[[1]]){gparsum=TRUE} else if ('gpp'%in%dimnames(results_effects)[[1]]){gparsum=FALSE} else {stop()}
  if (gparsum){
    results_effects = results_effects[c('proband','paternal','gp'),] 
    results_vcov = results_vcov[c('proband','paternal','gp'),
                                c('proband','paternal','gp')]
  } else{
    results_effects = results_effects[c('proband','paternal','gpp','gmp'),1:2]
    results_vcov = results_vcov[c('proband','paternal','gpp','gmp'),
                                c('proband','paternal','gpp','gmp')]
  }
  # Transform and save
  if (!gparsum){
    # Transformation matrix for full model
    A_3 = matrix(0,nrow=2,ncol=4)
    A_3[1,2] = 1
    A_3[2,] = c(0,0,0.5,0.5)} else{
                           # Transformation matrix for grandparental sum model
                           A_3 = matrix(0,nrow=2,ncol=3)
                           A_3[1:2,2:3]=diag(2)}
  results_effects = A_3%*%as.matrix(results_effects[,1])
  gen3_effects = c('paternal','grandpaternal')
  results_vcov = A_3%*%as.matrix(results_vcov)%*%t(A_3)
  dimnames(results_vcov)[[1]] = gen3_effects
  dimnames(results_vcov)[[2]] = gen3_effects
  # Compute indirect to direct ratios
  results[gen3_effects,1]=results_effects
  results[gen3_effects,2]=sqrt(diag(results_vcov))           
  # maternal
  if (lme4){
    results_effects = read.table(gen3_maternal,skip=1,header=F,row.names=1)
  } else {
    results_effects = read.table(gen3_maternal,row.names=1)}
  results_vcov = read.table(gen3_maternal_vcov,row.names=1,header=T)
  if ('gm'%in%dimnames(results_effects)[[1]]){gparsum=TRUE} else if ('gmm'%in%dimnames(results_effects)[[1]]){gparsum=FALSE} else {stop()}
  if (gparsum){
    results_effects = results_effects[c('proband','maternal','gm'),] 
    results_vcov = results_vcov[c('proband','maternal','gm'),
                                c('proband','maternal','gm')]
  } else{
    results_effects = results_effects[c('proband','maternal','gpm','gmm'),1:2]
    results_vcov = results_vcov[c('proband','maternal','gpm','gmm'),
                                c('proband','maternal','gpm','gmm')]
  }
  # Transform and save
  if (!gparsum){
    # Transformation matrix for full model
    A_3 = matrix(0,nrow=2,ncol=4)
    A_3[1,2] = 1
    A_3[2,] = c(0,0,0.5,0.5)} else{
                           # Transformation matrix for grandparental sum model
                           A_3 = matrix(0,nrow=2,ncol=3)
                           A_3[1:2,2:3]=diag(2)}
  results_effects = A_3%*%as.matrix(results_effects[,1])
  gen3_effects = c('maternal','grandmaternal')
  results_vcov = A_3%*%as.matrix(results_vcov)%*%t(A_3)
  dimnames(results_vcov)[[1]] = gen3_effects
  dimnames(results_vcov)[[2]] = gen3_effects
  # Compute indirect to direct ratios
  results[gen3_effects,1]=results_effects
  results[gen3_effects,2]=sqrt(diag(results_vcov))                   
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
moba_dir = '../moba'
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
botnia_dir = '../Botnia/pgs_results/'
botnia_traits = read.table('../Botnia/pgs_results/traits.txt',header=F)
botnia_traits=c('EA','glucose','HDL','non_HDL','height','BMI','DBP','SBP')
# GS 
gs_dir = '../GS/pgs/'
#gs_traits = read.table('GS/trait_names.txt',header=F)           
gs_traits = c("Glucose","Non_HDL", "HDL","height","FEV1","BMI","ever.smoked",
                "cigarettes.per.day","cog","vocab","SBP","DBP","neuroticism",
                "EA_years","EA_years_mid","EA_quals")
#gs_traits = cbind(1:length(gs_traits),gs_traits)
# FSH
fhs_dir = '../FHS/FHS_Unipar_Results/'
fhs_binary_dir = '../FHS/FHS_Unipar_Results/'

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
                        number_of_children_women=c(NA,'NEB',NA,NA,'NC_WOMEN'),
                        depressive_symptoms=c(NA,'CESD','neuroticism','Depression','depression'),
                        ADHD=c(NA,NA,NA,'ADHD','ADHD'),
                        alcohol_use_disorder=c(NA,'AUD',NA,NA,'alcohol_use_disorder'))

phenotypes = dimnames(phenotype_names)[[2]][-1]
binary_phens = c('ADHD','ever_smoker','hypertension','alcohol_use_disorder','depression') 

# Record which phenotypes in which cohort
in_cohort = t(phenotype_names)
dimnames(in_cohort)[[2]] = in_cohort[1,]
in_cohort = in_cohort[-1,]
in_cohort = !is.na(in_cohort)

pgs_names = data.frame(cohort=c('botnia','fhs','gs','moba','finngen'),
                              ADHD=c(NA,NA,NA,'ADHD_Demontis_2023','ADHD1'),
                              age_at_first_birth=c(NA,'AFB_Repo_pgs',NA,NA,'AFB2'),
                              BMI=c('BMI_UKB','BMI_UKB_pgs','bmi','BMI_UKB','bmi'),
                              depression=c(NA,'MDD_pgs','depression','PGC_UKB_depresssion','DEP1'),
                              educational_attainment= c('EA4_2.8m','EA4_pgs','EA4_hm3','EA4','EA6'),
                              ever_smoker= c(NA,'EVSMK_UKB_pgs','ever_smoke',NA,'EVERSMOKE2'),
                              height=c('height_UKB','HGT_UKB_pgs','height','height_yengo_2022','height'),
                              number_of_children_women=c(NA,'NEB_Repo_pgs',NA,NA,'NEBwomen2'))

pgss = dimnames(pgs_names)[[2]][-1] 

## Meta-analysis results table
# Effects
effect_names =  c('population','direct','paternal_NTC','maternal_NTC','average_NTC','average_NTC_direct_ratio',
                  'maternal_minus_paternal','maternal_minus_paternal_direct_ratio',
                  'paternal','maternal','grandpaternal','grandmaternal')

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
      trait_index = match(cohort_names[1],botnia_traits)
      # Read estimates and transform
      botnia_estimates = read_gen_models(paste(botnia_dir,pgs_names_j[1],'/',trait_index,'.1.effects.txt',sep=''),
                                        paste(botnia_dir,pgs_names_j[1],'/',trait_index,'.2.effects.txt',sep=''),
                                        paste(botnia_dir,pgs_names_j[1],'/',trait_index,'.2.vcov.txt',sep=''),
                                        paste(botnia_dir,pgs_names_j[1],'/',trait_index,'.3.paternal.effects.txt',sep=''),
                                        paste(botnia_dir,pgs_names_j[1],'/',trait_index,'.3.paternal.vcov.txt',sep=''),
                                        paste(botnia_dir,pgs_names_j[1],'/',trait_index,'.3.maternal.effects.txt',sep=''),
                                        paste(botnia_dir,pgs_names_j[1],'/',trait_index,'.3.maternal.vcov.txt',sep=''))
      # Store for meta-analysis
      estimates[effect_names,1] = botnia_estimates[effect_names,1]
      estimate_ses[effect_names,1] = botnia_estimates[effect_names,2]
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
                                        paste(fhs_binary_dir,pgs_names_j[2],'/',cohort_names[2],'.3.paternal.coefficients.txt',sep=''),
                                        paste(fhs_binary_dir,pgs_names_j[2],'/',cohort_names[2],'.3.paternal.vcov.txt',sep=''),
                                        paste(fhs_binary_dir,pgs_names_j[2],'/',cohort_names[2],'.3.maternal.coefficients.txt',sep=''),
                                        paste(fhs_binary_dir,pgs_names_j[2],'/',cohort_names[2],'.3.maternal.vcov.txt',sep=''),lme4=TRUE)
      } else {
        fhs_estimates = read_gen_models(paste(fhs_dir,pgs_names_j[2],'/',cohort_names[2],'.1.effects.txt',sep=''),
                                        paste(fhs_dir,pgs_names_j[2],'/',cohort_names[2],'.2.effects.txt',sep=''),
                                        paste(fhs_dir,pgs_names_j[2],'/',cohort_names[2],'.2.vcov.txt',sep=''),
                                        paste(fhs_dir,pgs_names_j[2],'/',cohort_names[2],'.3.paternal.effects.txt',sep=''),
                                        paste(fhs_dir,pgs_names_j[2],'/',cohort_names[2],'.3.paternal.vcov.txt',sep=''),
                                        paste(fhs_dir,pgs_names_j[2],'/',cohort_names[2],'.3.maternal.effects.txt',sep=''),
                                        paste(fhs_dir,pgs_names_j[2],'/',cohort_names[2],'.3.maternal.vcov.txt',sep=''))}
      # Store for meta-analysis
      estimates[effect_names,2] = fhs_estimates[effect_names,1]
      estimate_ses[effect_names,2] = fhs_estimates[effect_names,2]
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
                                        paste0(gs_dir,pgs_names_j[3],'/',cohort_names[3],'.3.paternal.coefficients.txt'),
                                        paste0(gs_dir,pgs_names_j[3],'/',cohort_names[3],'.3.paternal.vcov.txt'),
                                        paste0(gs_dir,pgs_names_j[3],'/',cohort_names[3],'.3.maternal.coefficients.txt'),
                                        paste0(gs_dir,pgs_names_j[3],'/',cohort_names[3],'.3.maternal.vcov.txt'),
                                        lme4=TRUE)
      } else {
        gs_estimates = read_gen_models(paste0(gs_dir,pgs_names_j[3],'/',trait_index,'.1.effects.txt'),
                                        paste0(gs_dir,pgs_names_j[3],'/',trait_index,'.2.effects.txt'),
                                        paste0(gs_dir,pgs_names_j[3],'/',trait_index,'.2.vcov.txt'),
                                        paste0(gs_dir,pgs_names_j[3],'/',trait_index,'.3.paternal.effects.txt'),
                                        paste0(gs_dir,pgs_names_j[3],'/',trait_index,'.3.paternal.vcov.txt'),
                                        paste0(gs_dir,pgs_names_j[3],'/',trait_index,'.3.maternal.effects.txt'),
                                        paste0(gs_dir,pgs_names_j[3],'/',trait_index,'.3.maternal.vcov.txt'))}
      # Store for meta-analysis
      estimates[effect_names,3] = gs_estimates[effect_names,1]
      estimate_ses[effect_names,3] = gs_estimates[effect_names,2]
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
                                    paste0(moba_dir,'/',pgs_names_j[4],'_',trait_index,'.3.paternal.effects.txt',sep=''),
                                    paste0(moba_dir,'/',pgs_names_j[4],'_',trait_index,'.3.paternal.vcov.txt',sep=''),
                                    paste0(moba_dir,'/',pgs_names_j[4],'_',trait_index,'.3.maternal.effects.txt',sep=''),
                                    paste0(moba_dir,'/',pgs_names_j[4],'_',trait_index,'.3.maternal.vcov.txt',sep=''))
      # Store for meta-analysis
      estimates[effect_names,4] = moba_estimates[effect_names,1]
      estimate_ses[effect_names,4] = moba_estimates[effect_names,2]
      # Store in MoBa results
      moba_results[j,match(phenotypes[i],dimnames(moba_results)$phenotype),effect_names] = moba_estimates[,1]
      moba_results[j,match(phenotypes[i],dimnames(moba_results)$phenotype),paste(effect_names,'SE',sep='_')] = moba_estimates[,2]
    }
    # Finngen
    if (!is.na(cohort_names[5]) & !is.na(pgs_names_j[5])){
      # Find trait index
      trait_index = match(cohort_names[5],finngen_traits)
      # Read estimates and transform
      #if (pgs_names_j[5]=='EA4'){sign_flip=TRUE} else {sign_flip=FALSE}
      if (phenotypes[i]%in%binary_phens & !phenotypes[i]=='ADHD'){
        finngen_estimates = read_gen_models(paste0(finngen_dir,pgs_names_j[5],'/',cohort_names[5],'.1.coefficients.txt'),
                                        paste0(finngen_dir,pgs_names_j[5],'/',cohort_names[5],'.2.coefficients.txt'),
                                        paste0(finngen_dir,pgs_names_j[5],'/',cohort_names[5],'.2.vcov.txt'),
                                        paste0(finngen_dir,pgs_names_j[5],'/',cohort_names[5],'.3.paternal.coefficients.txt'),
                                        paste0(finngen_dir,pgs_names_j[5],'/',cohort_names[5],'.3.paternal.vcov.txt'),
                                        paste0(finngen_dir,pgs_names_j[5],'/',cohort_names[5],'.3.maternal.coefficients.txt'),
                                        paste0(finngen_dir,pgs_names_j[5],'/',cohort_names[5],'.3.maternal.vcov.txt'),
                                        lme4=TRUE)
      } else {
        finngen_estimates = read_gen_models(paste0(finngen_dir,pgs_names_j[5],'/',trait_index,'.1.effects.txt'),
                                        paste0(finngen_dir,pgs_names_j[5],'/',trait_index,'.2.effects.txt'),
                                        paste0(finngen_dir,pgs_names_j[5],'/',trait_index,'.2.vcov.txt'),
                                        paste0(finngen_dir,pgs_names_j[5],'/',trait_index,'.3.paternal.effects.txt'),
                                        paste0(finngen_dir,pgs_names_j[5],'/',trait_index,'.3.paternal.vcov.txt'),
                                        paste0(finngen_dir,pgs_names_j[5],'/',trait_index,'.3.maternal.effects.txt'),
                                        paste0(finngen_dir,pgs_names_j[5],'/',trait_index,'.3.maternal.vcov.txt'))}
      # Store for meta-analysis
      estimates[effect_names,5] = finngen_estimates[effect_names,1]
      estimate_ses[effect_names,5] = finngen_estimates[effect_names,2]
      # Store in GS results
      finngen_results[j,match(phenotypes[i],dimnames(finngen_results)$phenotype),effect_names] = finngen_estimates[,1]
      finngen_results[j,match(phenotypes[i],dimnames(finngen_results)$phenotype),paste(effect_names,'SE',sep='_')] = finngen_estimates[,2]
    }
    # Meta-analysis
    for (effect_j in effect_names){
      if (effect_j%in%c('average_NTC_direct_ratio','maternal_minus_paternal_direct_ratio')){
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
pgs_primary = data.frame(pgs=dimnames(meta_results)[[1]],
              phenotype=c('ADHD','age_at_first_birth_women',
              'BMI','depressive_symptoms','educational_attainment',
              'ever_smoker','height','number_of_children_women'),
              average_NTC_direct_ratio=NA,average_NTC_direct_ratio_SE=NA,
              average_NTC_log10P=NA,
              maternal_minus_paternal_direct_ratio=NA,
              maternal_minus_paternal_direct_ratio_SE=NA,
              maternal_minus_paternal_log10P=NA,
              paternal_Z=NA,paternal_log10P=NA,maternal_Z=NA,maternal_log10P=NA,
              grandpaternal_Z=NA,grandpaternal_log10P=NA,grandmaternal_Z=NA,grandmaternal_log10P=NA)

for (i in 1:dim(pgs_primary)[1]){
  pgs_primary[i,'average_NTC_direct_ratio'] = meta_results[pgs_primary[i,'pgs'],pgs_primary[i,'phenotype'],'average_NTC_direct_ratio']
  pgs_primary[i,'average_NTC_direct_ratio_SE'] = meta_results[pgs_primary[i,'pgs'],pgs_primary[i,'phenotype'],'average_NTC_direct_ratio_SE']
  pgs_primary[i,'average_NTC_log10P'] = meta_results[pgs_primary[i,'pgs'],pgs_primary[i,'phenotype'],'average_NTC_log10P'] 
  pgs_primary[i,'maternal_minus_paternal_direct_ratio'] = meta_results[pgs_primary[i,'pgs'],pgs_primary[i,'phenotype'],'maternal_minus_paternal_direct_ratio']
  pgs_primary[i,'maternal_minus_paternal_direct_ratio_SE'] = meta_results[pgs_primary[i,'pgs'],pgs_primary[i,'phenotype'],'maternal_minus_paternal_direct_ratio_SE']
  pgs_primary[i,'maternal_minus_paternal_log10P'] = meta_results[pgs_primary[i,'pgs'],pgs_primary[i,'phenotype'],'maternal_minus_paternal_log10P']
  pgs_primary[i,'paternal_Z'] = meta_results[pgs_primary[i,'pgs'],pgs_primary[i,'phenotype'],'paternal']
  pgs_primary[i,'paternal_log10P'] = meta_results[pgs_primary[i,'pgs'],pgs_primary[i,'phenotype'],'paternal_log10P']
  pgs_primary[i,'maternal_Z'] = meta_results[pgs_primary[i,'pgs'],pgs_primary[i,'phenotype'],'maternal']
  pgs_primary[i,'maternal_log10P'] = meta_results[pgs_primary[i,'pgs'],pgs_primary[i,'phenotype'],'maternal_log10P']
  pgs_primary[i,'grandpaternal_Z'] = meta_results[pgs_primary[i,'pgs'],pgs_primary[i,'phenotype'],'grandpaternal']
  pgs_primary[i,'grandpaternal_log10P'] = meta_results[pgs_primary[i,'pgs'],pgs_primary[i,'phenotype'],'grandpaternal_log10P']
  pgs_primary[i,'grandmaternal_Z'] = meta_results[pgs_primary[i,'pgs'],pgs_primary[i,'phenotype'],'grandmaternal']
  pgs_primary[i,'grandmaternal_log10P'] = meta_results[pgs_primary[i,'pgs'],pgs_primary[i,'phenotype'],'grandmaternal_log10P']
}

write.csv(pgs_primary,'~/Dropbox/grandparental/pgs_primary.csv',quote=F,row.names=F)

pgs_primary$ntc_lower = pgs_primary[,'average_NTC_direct_ratio']+qnorm(0.025)*pgs_primary[,'average_NTC_direct_ratio_SE']
pgs_primary$ntc_upper = pgs_primary[,'average_NTC_direct_ratio']-qnorm(0.025)*pgs_primary[,'average_NTC_direct_ratio_SE']
pgs_primary$patdiff_lower = pgs_primary[,'maternal_minus_paternal_direct_ratio']+qnorm(0.025)*pgs_primary[,'maternal_minus_paternal_direct_ratio_SE']
pgs_primary$patdiff_upper = pgs_primary[,'maternal_minus_paternal_direct_ratio']-qnorm(0.025)*pgs_primary[,'maternal_minus_paternal_direct_ratio_SE']


pgs_primary$pgs = factor(pgs_primary$pgs,levels=pgs_primary[order(abs(pgs_primary$average_NTC_direct_ratio)),'pgs'])

plot_gen = rbind(data.frame(pgs=pgs_primary$pgs,phenotype=pgs_primary$phenotype,effect='average_NTC_direct_ratio',
                            value=pgs_primary$average_NTC_direct_ratio,lower=pgs_primary$ntc_lower,upper=pgs_primary$ntc_upper),
                 data.frame(pgs=pgs_primary$pgs,phenotype=pgs_primary$phenotype,effect='maternal_minus_paternal_direct_ratio',
                 value=pgs_primary$maternal_minus_paternal_direct_ratio,lower=pgs_primary$patdiff_lower,upper=pgs_primary$patdiff_upper))

require(ggplot2)
par_plot = ggplot(plot_gen, aes(y=value, x=pgs, fill=effect)) +
    geom_errorbar(aes(ymin = lower, ymax = upper,color=effect),width=0.1,position=position_dodge(.5))+coord_flip()+
    geom_point(position=position_dodge(.5),stat="identity",aes(color=effect))+theme_minimal() +theme(axis.line = element_line(color="black"),
                          axis.ticks = element_line(color="black"),
                          panel.border = element_blank(),
                          axis.text.x = element_text(angle = 45, hjust=1))+
                          geom_hline(yintercept=0,linetype='dashed')+
                          ylab('Ratio with direct genetic effect')+
                          xlab('PGI')

ggsave('~/Dropbox/grandparental/ratio_plot.pdf',plot=par_plot,width=7,height=5)

############# QQ Plots ###############
qqplot_ci = function(pmatrix,outfile){
  # QQ-plots
  require("MASS") 
  require("reshape2") 
  require('qqplotr')
  di <- "exp"
  dp <- list(rate = log(10))
  parqq=data.frame(pmatrix)
  parqq$pgs=dimnames(pmatrix)[[1]]
  parqq = melt(parqq)
  dimnames(parqq)[[2]] = c('pgs','phenotype','value')
  qqpar <- ggplot(data = parqq, mapping = aes(sample = value)) +
  stat_qq_band(distribution = di, dparams = dp, identity=TRUE) +
  stat_qq_line(distribution = di, dparams = dp, identity=TRUE) +
  stat_qq_point(distribution = di, dparams = dp)+
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
  ggsave(filename=outfile,plot=qqpar,width=5,height=5)
  return(qqpar)
}

qq_paternal=qqplot_ci(meta_results[,,'paternal_log10P'],
'~/Dropbox/grandparental/paternal_log10P_qqplot.pdf')

qq_maternal=qqplot_ci(meta_results[,,'maternal_log10P'],
'~/Dropbox/grandparental/maternal_log10P_qqplot.pdf')

qq_maternal_minus_paternal = qqplot_ci(meta_results[,,'maternal_minus_paternal_log10P'],
'~/Dropbox/grandparental/maternal_minus_paternal_log10P_qqplot.pdf') 

qq_avg_NTC = qqplot_ci(meta_results[,,'average_NTC_log10P'],
'~/Dropbox/grandparental/average_NTC_qqplot.pdf')

qq_grandparental = qqplot_ci(meta_results[,,'grandpaternal_log10P'],
'~/Dropbox/grandparental/grandpaternal_qqplot.pdf')

qq_grandparental = qqplot_ci(meta_results[,,'grandmaternal_log10P'],
'~/Dropbox/grandparental/grandmaternal_qqplot.pdf')