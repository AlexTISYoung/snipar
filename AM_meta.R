setwd('~/Google Drive/EA4_revision/AM_analysis/')
ukb = read.csv('UKB_AM_nonlinear.csv')
gs = read.csv('GS_AM_nonlinear.csv')

inv_var_meta = function(ests,ses){
  weights = ses^(-2)
  est = sum(weights*ests)/sum(weights)
  return(c(est,sqrt(1/sum(weights))))
}

ci_to_se = function(lower,upper){
  return((upper-lower)/(2*qnorm(0.975)))
}

meta = data.frame(trait=ukb$trait,
                  phenotypic_correlation=rep(NA,4),
                  phenotpic_correlation_SE=rep(NA,4),
                  r_p = rep(NA,4),
                  r_p_SE = rep(NA,4),
                  r_m = rep(NA,4),
                  r_m_SE = rep(NA,4),
                  prediction = rep(NA,4),
                  prediction_SE = rep(NA,4),
                  PGI_r = rep(NA,4),
                  PGI_r_SE = rep(NA,4),
                  PGI_trait_r = rep(NA,4),
                  PGI_trait_r_SE = rep(NA,4),
                  PGI_trait_PCs_r = rep(NA,4),
                  PGI_trait_PCs_r_SE = rep(NA,4))

for (i in 1:4){
  meta[i,2:3] = inv_var_meta(c(ukb[i,2],gs[i,2]),
                             c(ci_to_se(ukb[i,3],ukb[i,4]),
                               ci_to_se(gs[i,3],gs[i,4])))
  meta[i,4:5] = inv_var_meta(c(ukb[i,6],gs[i,6]),
                             c(ci_to_se(ukb[i,7],ukb[i,8]),
                               ci_to_se(gs[i,7],gs[i,8])))
  meta[i,6:7] = inv_var_meta(c(ukb[i,9],gs[i,9]),
                             c(ci_to_se(ukb[i,10],ukb[i,11]),
                               ci_to_se(gs[i,10],gs[i,11])))
  meta[i,8:9] = inv_var_meta(c(ukb[i,12],gs[i,12]),
                             c(ukb[i,13],gs[i,13]))
  meta[i,10:11] = inv_var_meta(c(ukb[i,14],gs[i,14]),
                               c(ci_to_se(ukb[i,15],ukb[i,16]),
                                 ci_to_se(gs[i,15],gs[i,16])))
  meta[i,12:13] = inv_var_meta(c(ukb[i,'PGI_trait_r'],gs[i,'PGI_trait_r']),
                               c(ci_to_se(ukb[i,'PGI_trait_r_lower'],ukb[i,'PGI_trait_r_upper']),
                                 ci_to_se(gs[i,'PGI_trait_r_lower'],gs[i,'PGI_trait_r_upper'])))
  meta[i,14:15] = inv_var_meta(c(ukb[i,'PGI_trait_PCs_r'],gs[i,'PGI_trait_PCs_r']),
                               c(ci_to_se(ukb[i,'PGI_trait_PCs_r_lower'],ukb[i,'PGI_trait_PCs_r_upper']),
                                 ci_to_se(gs[i,'PGI_trait_PCs_r_lower'],gs[i,'PGI_trait_PCs_r_upper'])))
}

library(xlsx)
write.xlsx(ukb,'AM_UKB_nonlinear.xlsx')
write.xlsx(gs,'AM_GS_nonlinear.xlsx')
write.xlsx(meta,'AM_meta_nonlinear.xlsx')

