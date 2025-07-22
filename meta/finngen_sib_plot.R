# Finngen
# Read gen models 1-3 function
read_sib_model = function(gen2_effects,gen2_vcov,sign_flip=FALSE,lme4=FALSE){
  ## Estimates to output 
  parameter_names = c('direct','sibling','paternal_NTC','maternal_NTC','average_NTC','maternal_minus_paternal')
  results = matrix(NA,nrow=length(parameter_names),ncol=2)
  dimnames(results)[[1]] = parameter_names
  dimnames(results)[[2]] = c('estimates','SE')
  ## Check files exist
 if (!file.exists(gen2_effects)){
    print(paste(gen2_effects,'does not exist'))
    return(results)
  } else if (!file.exists(gen2_vcov)){
    print(paste(gen2_vcov,'does not exist'))
    return(results)
  } 
  ## Get non-transmitted coefficients
  if (lme4){
   results_2gen = read.table(gen2_effects,row.names=1,skip=1,header=F) 
  } else {results_2gen = read.table(gen2_effects,row.names=1)}
  results_2gen = results_2gen[c('proband','sibling','paternal','maternal'),1:2]
  results_2gen_vcov = read.table(gen2_vcov,row.names=1)
  results_2gen_vcov = results_2gen_vcov[c('proband','sibling','paternal','maternal'),c('proband','sibling','paternal','maternal')]
  # Transform
  A_ntc = matrix(0,nrow=6,ncol=4)
  A_ntc[1:4,1:4] = diag(4)
  A_ntc[5,] = c(0,0,0.5,0.5)
  A_ntc[6,] = c(0,0,-1,1)
  results_2gen_vcov = A_ntc%*%as.matrix(results_2gen_vcov)%*%t(A_ntc)
  gen2_effects = c('direct','sibling','paternal_NTC','maternal_NTC','average_NTC','maternal_minus_paternal')
  dimnames(results_2gen_vcov)[[1]] = gen2_effects
  dimnames(results_2gen_vcov)[[2]] = gen2_effects
  results[gen2_effects,1] = A_ntc%*%as.matrix(results_2gen[,1])
  results[gen2_effects,2] = sqrt(diag(results_2gen_vcov))
  # Return
  if (sign_flip){results[,1] = -results[,1]}
  return(results)
}

setwd('finngen/pgs_out/')
finngen_traits = c('height','NC','BMI','ever_smoker',
                'depression','AAFB','AAFB_WOMEN',
                'AAFB_MEN','NC_WOMEN','NC_MEN',
                'ADHD','asthma','eczema','hypertension',
                'alcohol_use_disorder','allergic_rhinitis',
                'migraine','copd')

pgs_names = data.frame(cohort=c('botnia','fhs','gs','moba','finngen'),
                              ADHD=c(NA,NA,NA,'ADHD_Demontis_2023','ADHD1'),
                              age_at_first_birth=c(NA,'AFB_Repo_pgs',NA,NA,'AFB2'),
                              BMI=c('BMI_UKB','BMI_UKB_pgs','bmi','BMI_UKB','bmi'),
                              depression=c(NA,'MDD_pgs','depression','PGC_UKB_depresssion','DEP1'),
                              educational_attainment= c('EA4_2.8m','EA4_pgs','EA4_hm3','EA4','EA6'),
                              ever_smoker= c(NA,'EVSMK_UKB_pgs','ever_smoke',NA,'EVERSMOKE2'),
                              height=c('height_UKB','HGT_UKB_pgs','height','height_yengo_2022','height'),
                              number_of_children_women=c(NA,'NEB_Repo_pgs',NA,NA,'NEBwomen2'))

pgs_names = pgs_names[5,c(2,3,5,7)]

finngen_results = array(NA,dim=c(length(pgs_names),length(finngen_traits),6,2))
dimnames(finngen_results) = list(pgs_names,finngen_traits,
                                    c('direct','sibling','paternal_NTC','maternal_NTC','average_NTC','maternal_minus_paternal'),
                            c('estimate','SE'))

for (trait_index in 1:18){
    for (pgs_name in pgs_names){
    finngen_estimates = read_sib_model(paste0(pgs_name,'/',trait_index,'_sib.2.effects.txt'),
                                        paste0(pgs_name,'/',trait_index,'_sib.2.vcov.txt'))
      # Store in GS results
      finngen_results[match(pgs_name,pgs_names),trait_index,,] = finngen_estimates
    }
}