setwd('~/Dropbox/intergenerational/simulations/')

### General ggplot2 scatterplot function ###
scatter_plot = function(x,x_SE,y,y_SE,pch,pch_name,color,color_name,xlab,ylab,outfile,hlines=TRUE){
  require(ggplot2)
  require(gridExtra)
  require(reshape2)
  # Plot data frame
  plot_df = data.frame(x=x,y=y,point_type=pch,colour=color)
  # CIs
  x_lower = x+qnorm(0.025)*x_SE
  x_upper = x-qnorm(0.025)*x_SE
  y_lower = y+qnorm(0.025)*y_SE
  y_upper = y-qnorm(0.025)*y_SE
  # Plot
  splot = ggplot(plot_df,aes(x = x, y = y, color=colour, pch=point_type)) +
    geom_point()+geom_errorbarh(aes(xmin = x_lower,xmax=x_upper))+
    geom_errorbar(aes(ymin = y_lower,ymax=y_upper))+
    theme_minimal() + theme(axis.line = element_line(color="black"),
                            axis.ticks = element_line(color="black"),
                            panel.border = element_blank(),
                            axis.text.x = element_text(angle = 45, hjust=1))+
    geom_abline(slope=1,intercept=0)+xlab(xlab)+ylab(ylab)+labs(color = color_name, pch=pch_name)
  if (hlines==TRUE){splot = splot+geom_hline(yintercept=c(0,1),linetype=2)}
  ggsave(outfile,plot=splot,width=7,height=5)
  return(splot)
}

########## Functions to compute correlations from pedigree ###########
compute_cousin_cor=function(pedigree,last_gen){
  ngen = as.integer(strsplit(pedigree[dim(pedigree)[1],1],'_')[[1]][1])
  parent_gen = pedigree[sapply(pedigree[,1],function(x) strsplit(x,'_')[[1]][1]==as.character(ngen-1)),]
  parent_fams = unique(parent_gen[,1])
  parent_sibs = t(sapply(parent_fams,function(x) parent_gen[parent_gen[,1]==x,'IID']))
  c1_phenotypes = t(sapply(parent_sibs[,1],function(x) last_gen[last_gen$FATHER_ID==x,'PHENO']))
  c2_phenotypes = t(sapply(parent_sibs[,2],function(x) last_gen[last_gen$MOTHER_ID==x,'PHENO']))
  return(sum(cor(c1_phenotypes,c2_phenotypes))/4)
}

compute_gpar_cor=function(pedigree,last_gen){
  gpar_ids = cbind(pedigree[match(last_gen$FATHER_ID,pedigree$IID),c('FATHER_ID','MOTHER_ID')],
                   pedigree[match(last_gen$MOTHER_ID,pedigree$IID),c('FATHER_ID','MOTHER_ID')])
  gpar_phenotypes = cbind(pedigree[match(gpar_ids[,1],pedigree$IID),'PHENO'],
                          pedigree[match(gpar_ids[,2],pedigree$IID),'PHENO'],
                          pedigree[match(gpar_ids[,3],pedigree$IID),'PHENO'],
                          pedigree[match(gpar_ids[,4],pedigree$IID),'PHENO'])
  return(sum(cor(last_gen$PHENO,gpar_phenotypes))/4)
}

compute_corrs=function(pedigree){
  ngen = strsplit(pedigree[dim(pedigree)[1],1],'_')[[1]][1]
  last_gen = pedigree[sapply(pedigree[,1],function(x) strsplit(x,'_')[[1]][1]==ngen),]
  sib_cor = cor(last_gen[seq(1,dim(last_gen)[1],2),'PHENO'],
                     last_gen[seq(2,dim(last_gen)[1],2),'PHENO'])
  #### FIX MOTHER PHENO ###
  last_gen$MOTHER_PHENO = pedigree[match(last_gen$MOTHER_ID,pedigree[,2]),'PHENO']
  last_gen$FATHER_PHENO = pedigree[match(last_gen$FATHER_ID,pedigree[,2]),'PHENO']
  po_cor = (cor(last_gen$PHENO,last_gen$FATHER_PHENO)+cor(last_gen$PHENO,last_gen$MOTHER_PHENO))/2
  cousin_cor = compute_cousin_cor(pedigree,last_gen)
  gpar_cor = compute_gpar_cor(pedigree,last_gen)
  outcors = c(sib_cor,cousin_cor,po_cor,gpar_cor)
  return(outcors)
}

##### Predicted correlation functions #####
predicted_cor=function(k,m,v_g,v_eg,c_ge,v_y,r_delta,r_eta){
  r_out=((1+r_delta)/2)^(2*k+m+1)*v_g+
    ((1+r_eta)/2)^(2*k+m)*v_eg+
    0.5*((1+r_delta)/2)^k*((1+r_eta)/2)^k*(((1+r_delta)/2)^m+((1+r_eta)/2)^m)*c_ge
  return(r_out/v_y)
}

predicted_sib_cor=function(v_g,v_eg,c_ge,v_y,r_delta){
  return(predicted_cor(0,0,v_g,v_eg,c_ge,v_y,r_delta,0))
}

predicted_cousin_cor=function(v_g,v_eg,c_ge,v_y,r_delta,r_eta){
  return(predicted_cor(1,0,v_g,v_eg,c_ge,v_y,r_delta,r_eta))
}

predicted_cousin_cor=function(v_g,v_eg,c_ge,v_y,r_delta,r_eta){
  r_cousin=((1+r_delta)/2)^3*v_g+((1+r_eta)/2)^2*v_eg+((1+r_delta)*(1+r_eta)/4)*c_ge
  return(r_cousin/v_y)
}

predicted_ancestor_cor=function(k,v_g,v_eg,c_ge,v_y,r_delta,r_eta){
  r_out=((1+r_delta)/2)^(k)*v_g+
    ((1+r_eta)/2)^(k)*v_eg+
    0.5*(((1+r_delta)/2)^k+((1+r_eta)/2)^(k-1))*c_ge
  return(r_out/v_y)
}

predicted_c_ge_eq=function(v_g,v_eg,r_delta,r_eta,r_delta_eta_c,r_delta_eta_tau){
  c_ge_eq = sqrt(2*v_g*v_eg/((1-r_delta)*(1-r_eta)))*(r_delta_eta_c+r_delta_eta_tau)
  return(c_ge_eq)
}

predicted_r_po=function(v_g,v_eg,c_ge,v_y,r_delta,r_eta,r_delta_eta_c,r_y){
  r=predicted_ancestor_cor(1,v_g,v_eg,c_ge,v_y,r_delta,r_eta)
  T_eq = (v_g+v_eg+c_ge)/v_y
  r=r+r_y*(1-T_eq)*((v_g+c_ge/2+v_eg)/(2*v_y)+r_delta_eta_c*sqrt((v_g*v_eg)/(2*v_y^2*(1+r_eta))))
  return(r)
}

###################################### Read results #########################################
statistic = c('v_g','v_g_eq','v_eg','v_eg_eq','c_ge','c_ge_eq','v_y_eq',
               'r_sib','r_cousin','r_po','r_gpar')

results = data.frame(r_y = rep(c(0, 0.25, 0.5, 0.75),4),
                     v_indir=c(rep(0,4),rep(0.25,12)),
                     r_dir_indir=c(rep(0,8),rep(0.5,4),rep(1,4)))

pgi_results_true = cbind(results,data.frame(delta=rep(NA,dim(results)[1]),delta_SE=rep(NA,dim(results)[1]),alpha=rep(NA,dim(results)[1]),alpha_SE=rep(NA,dim(results)[1]),rk=rep(NA,dim(results)[1]),rk_SE=rep(NA,dim(results)[1]),k=rep(NA,dim(results)[1]),k_SE=rep(NA,dim(results)[1]),r=rep(NA,dim(results)[1]),r_SE=rep(NA,dim(results)[1]),h2_eq=rep(NA,dim(results)[1]),h2_eq_SE=rep(NA,dim(results)[1]),rho=rep(NA,dim(results)[1]),rho_SE=rep(NA,dim(results)[1]),alpha_delta=rep(NA,dim(results)[1]),alpha_delta_SE=rep(NA,dim(results)[1]),v_eta_delta=rep(NA,dim(results)[1]),v_eta_delta_SE=rep(NA,dim(results)[1])))
pgi_results_v1 = pgi_results_true
pgi_results_v10 = pgi_results_true
pgi_results_v100 = pgi_results_true

pop_pgi_results_true = pgi_results_true
pop_pgi_results_v1 = pgi_results_true
pop_pgi_results_v10 = pgi_results_true
pop_pgi_results_v100 = pgi_results_true

results=cbind(results,data.frame(v_g=rep(NA,dim(results)[1]),
                                 v_g_eq=rep(NA,dim(results)[1]),
                                 v_eg=rep(NA,dim(results)[1]),
                                 v_eg_eq=rep(NA,dim(results)[1]),
                                 c_ge=rep(NA,dim(results)[1]),
                                 c_ge_eq=rep(NA,dim(results)[1]),
                                 v_y_eq=rep(NA,dim(results)[1]),
                                 r_delta=rep(NA,dim(results)[1]),
                                 r_eta=rep(NA,dim(results)[1]),
                                 r_delta_eta_c=rep(NA,dim(results)[1]),
                                 r_delta_eta_tau=rep(NA,dim(results)[1]),
                                 r_sib=rep(NA,dim(results)[1]),
                                 r_cousin=rep(NA,dim(results)[1]),
                                 r_po=rep(NA,dim(results)[1]),
                                 r_gpar=rep(NA,dim(results)[1])))

sim = 'r_y_0.5'

for (i in 1:dim(results)[1]){
  if (results$v_indir[i]==0){
    simname = paste('r_y',results$r_y[i],sep='_')
  } else {
    simname = paste('v_indir',results$v_indir[i],'r_dir_indir',results$r_dir_indir[i],'r_y',results$r_y[i],sep='_')
  }
  ped = read.table(paste(simname,'ped',sep='.'),sep='',header=T,stringsAsFactors = F)
  # Variance component results
  vcs = read.table(paste(simname,'_VCs.txt',sep=''),header=T,stringsAsFactors = F)
  ngen = dim(vcs)[1]
  if (results$v_indir[i]==0){
    results[i,3+1:11] = c(vcs[1,'v_g'],vcs[ngen,'v_g'],0,0,0,0,vcs[ngen,'v_y'],
                         vcs[ngen,'r_delta'],0,0,0)
  } else {
    results[i,3+1:11] = c(vcs[1,'v_g'],vcs[ngen,'v_g'],
                          vcs[1,'v_eg'],vcs[ngen,'v_eg'],
                          vcs[1,'c_ge'],vcs[ngen,'c_ge'],
                          vcs[ngen,'v_y'],
                          vcs[ngen,'r_delta'],vcs[ngen,'r_eta'],
                          vcs[ngen,'r_delta_eta_c'],vcs[ngen,'r_delta_eta_tau'])
  }
  results[i,3+12:15] = compute_corrs(ped)
  # PGI results
  pgi_results_true[i,4:dim(pgi_results_true)[2]] = c(t(read.table(paste(simname,'pgi_v0.am_adj_pars.txt',sep='_'),header=T,row.names=1)))
  pgi_results_v1[i,4:dim(pgi_results_true)[2]] = c(t(read.table(paste(simname,'pgi_v1.am_adj_pars.txt',sep='_'),header=T,row.names=1)))
  pgi_results_v10[i,4:dim(pgi_results_true)[2]] = c(t(read.table(paste(simname,'pgi_v10.am_adj_pars.txt',sep='_'),header=T,row.names=1)))
  pgi_results_v100[i,4:dim(pgi_results_true)[2]] = c(t(read.table(paste(simname,'pgi_v100.am_adj_pars.txt',sep='_'),header=T,row.names=1)))
  if (results$v_indir[i]>0){
    pop_pgi_results_true[i,4:dim(pgi_results_true)[2]] = c(t(read.table(paste(simname,'population_pgi_v0.am_adj_pars.txt',sep='_'),header=T,row.names=1)))
    pop_pgi_results_v1[i,4:dim(pgi_results_true)[2]] = c(t(read.table(paste(simname,'population_pgi_v1.am_adj_pars.txt',sep='_'),header=T,row.names=1)))
    pop_pgi_results_v10[i,4:dim(pgi_results_true)[2]] = c(t(read.table(paste(simname,'population_pgi_v10.am_adj_pars.txt',sep='_'),header=T,row.names=1)))
    pop_pgi_results_v100[i,4:dim(pgi_results_true)[2]] = c(t(read.table(paste(simname,'population_pgi_v100.am_adj_pars.txt',sep='_'),header=T,row.names=1)))
  }
  print(simname)
}
save.image('~/Dropbox/intergenerational/simulations/sim_results.RData')
write.csv(results,'~/Dropbox/intergenerational/simulations/results_table.csv',quote=F,row.names=F)
write.csv(pgi_results_true,'~/Dropbox/intergenerational/simulations/true_direct_effect_pgi_results_table.csv',quote=F,row.names=F)
write.csv(pgi_results_v1,'~/Dropbox/intergenerational/simulations/v1_direct_effect_pgi_results_table.csv',quote=F,row.names=F)
write.csv(pgi_results_v10,'~/Dropbox/intergenerational/simulations/v10_direct_effect_pgi_results_table.csv',quote=F,row.names=F)
write.csv(pgi_results_v100,'~/Dropbox/intergenerational/simulations/v100_direct_effect_pgi_results_table.csv',quote=F,row.names=F)
write.csv(pop_pgi_results_v1,'~/Dropbox/intergenerational/simulations/v1_pop_effect_pgi_results_table.csv',quote=F,row.names=F)
write.csv(pop_pgi_results_v10,'~/Dropbox/intergenerational/simulations/v10_pop_effect_pgi_results_table.csv',quote=F,row.names=F)
write.csv(pop_pgi_results_v100,'~/Dropbox/intergenerational/simulations/v100_pop_effect_pgi_results_table.csv',quote=F,row.names=F)


#################################### Plot Results ###################################
load('~/Dropbox/intergenerational/simulations/sim_results.RData')
################### Variance Components ###############
v_g_plot = scatter_plot(1/(1-results$r_delta),0,
                          results$v_g_eq/results$v_g,0,
                          pch=pgi_results_true$v_indir==0,pch_name='v_indirect==0',
                          color=pgi_results_true$r_y,color_name='Phenotype correlation',
                          xlab='Predicted inflation of genetic variance',
                          ylab='Observed inflation of genetic variance',
                          outfile='plots/v_g.pdf')

v_eg_plot = scatter_plot((1+results$r_eta[results$v_indir>0])/(1-results$r_eta[results$v_indir>0]),0,
                        results$v_eg_eq[results$v_indir>0]/results$v_eg[results$v_indir>0],0,
                        pch=as.factor(pgi_results_true$r_dir_indir[results$v_indir>0]),pch_name='r_direct_indirect',
                        color=pgi_results_true$r_y[results$v_indir>0],color_name='Phenotype correlation',
                        xlab='Predicted inflation of variance due to parental IGEs',
                        ylab='Observed inflation of variance due to parental IGEs',
                        outfile='plots/v_eg.pdf')

results$predicted_c_ge_eq=predicted_c_ge_eq(results$v_g,results$v_eg,results$r_delta,results$r_eta,results$r_delta_eta_c,
                                            results$r_delta_eta_tau)

c_ge_plot = scatter_plot(results$predicted_c_ge_eq/results$v_y_eq,0,
                         results$c_ge_eq/results$v_y_eq,0,
                         pch=pgi_results_true$r_y==0,pch_name='Random-mating',
                         color=pgi_results_true$r_dir_indir,color_name='r_direct_indirect',
                         xlab='Predicted equilibrium covariance between direct and indirect genetic effects',
                         ylab='Observed equilibrium covariance between direct and indirect genetic effects',
                         outfile='plots/c_ge.pdf')

c_ge_inflation_plot = scatter_plot(results$c_ge_eq[results$r_dir_indir>0]/results$c_ge[results$r_dir_indir>0],0,
                         (results$r_delta_eta_c[results$r_dir_indir>0]+results$r_delta_eta_tau[results$r_dir_indir>0])/(results$r_delta_eta_c[results$r_dir_indir>0]-results$r_delta_eta_tau[results$r_dir_indir>0]),0,
                         pch=pgi_results_true$v_indir[results$r_dir_indir>0]==0,pch_name='v_indirect==0',
                         color=pgi_results_true$r_y[results$r_dir_indir>0],color_name='Phenotype correlation',
                         xlab='Predicted inflation of covariance between direct and indirect genetic effects',
                         ylab='Observed inflation of covariance between direct and indirect genetic effects',
                         outfile='plots/c_ge_inflation.pdf',hlines=FALSE)

################### Relative Correlations #############
## Perform theoretical predictions of relative correlations
results$predicted_r_sib = predicted_sib_cor(results$v_g_eq,results$v_eg_eq,
                                            results$c_ge_eq,results$v_y_eq,
                                            results$r_delta)

r_sib_plot = scatter_plot(results$predicted_r_sib,0,
                         results$r_sib,0,
                         pch=pgi_results_true$v_indir==0,pch_name='v_indirect==0',
                         color=pgi_results_true$r_y,color_name='Phenotype correlation',
                         xlab='Theoretical sibling correlation',
                         ylab='Observed sibling correlation',
                         outfile='plots/r_sib.pdf')

results$predicted_r_cousin = predicted_cousin_cor(results$v_g_eq,results$v_eg_eq,
                                                 results$c_ge_eq,results$v_y_eq,
                                                 results$r_delta,results$r_eta)

r_cousin_plot = scatter_plot(results$predicted_r_cousin,0,
                          results$r_cousin,0,
                          pch=pgi_results_true$v_indir==0,pch_name='v_indirect==0',
                          color=pgi_results_true$r_y,color_name='Phenotype correlation',
                          xlab='Theoretical cousin correlation',
                          ylab='Observed cousin correlation',
                          outfile='plots/r_cousin.pdf')

results$predicted_r_po = predicted_r_po(results$v_g_eq,results$v_eg_eq,
                                                  results$c_ge_eq,results$v_y_eq,
                                                  results$r_delta,results$r_eta,
                                        results$r_delta_eta_c,results$r_y)

r_po_plot = scatter_plot(results$predicted_r_po,0,
                             results$r_po,0,
                             pch=pgi_results_true$v_indir==0,pch_name='v_indirect==0',
                             color=pgi_results_true$r_y,color_name='Phenotype correlation',
                             xlab='Theoretical parent-offspring correlation',
                             ylab='Observed parent-offspring correlation',
                             outfile='plots/r_po.pdf')


results$predicted_r_gpar = predicted_ancestor_cor(2,results$v_g_eq,results$v_eg_eq,
                                                results$c_ge_eq,results$v_y_eq,
                                                results$r_delta,results$r_eta)

r_gpar_plot = scatter_plot(results$predicted_r_gpar,0,
                         results$r_gpar,0,
                         pch=pgi_results_true$v_indir==0,pch_name='v_indirect==0',
                         color=pgi_results_true$r_y,color_name='Phenotype correlation',
                         xlab='Theoretical grandparent-grandchild correlation',
                         ylab='Observed grandparent-grandchild correlation',
                         outfile='plots/r_grandpar.pdf')


########################## PGI results ############################

results$alpha_delta = results$c_ge_eq/(2*(1+results$r_delta)*results$v_g_eq)

results$alpha_delta_theoretical = ((results$r_delta_eta_c+results$r_delta_eta_tau)/(1+results$r_delta))*sqrt((1-results$r_delta)*results$v_eg/(2*(1-results$r_eta)*results$v_g))


theoretical_alpha_delta = scatter_plot(results$alpha_delta_theoretical,0,
                                       results$alpha_delta,0,
                                       pch=pgi_results_true$v_indir==0,pch_name='v_indirect==0',
                                       color=pgi_results_true$r_dir_indir,color_name='r_direct_indirect',
                                       xlab='Theoretical alpha_delta (using random-mating variance components)',
                                       ylab='Theoretical alpha_delta (using equilibrium variance components)',
                                       outfile='plots/alpha_delta_v0.pdf')

v0_alpha_delta = scatter_plot(results$alpha_delta,0,
                                       pgi_results_true$alpha_delta,pgi_results_true$alpha_delta_SE,
                                       pch=pgi_results_true$v_indir==0,pch_name='v_indirect==0',
                                       color=pgi_results_true$r_dir_indir,color_name='r_direct_indirect',
                                       xlab='Theoretical',ylab='True PGI alpha_delta',
                                       outfile='plots/alpha_delta_v0.pdf')



v1_alpha_delta = scatter_plot(pgi_results_true$alpha_delta,pgi_results_true$alpha_delta_SE,
                       pgi_results_v1$alpha_delta,pgi_results_v1$alpha_delta_SE,
                       pch=pgi_results_true$v_indir==0,pch_name='v_indirect==0',
                       color=pgi_results_true$r_dir_indir,color_name='r_direct_indirect',
                       xlab='True alpha_delta',ylab='Estimated alpha_delta (noise=1)',
                       outfile='plots/alpha_delta_v1.pdf')

v10_alpha_delta = scatter_plot(pgi_results_true$alpha_delta,pgi_results_true$alpha_delta_SE,
                              pgi_results_v10$alpha_delta,pgi_results_v10$alpha_delta_SE,
                              pch=pgi_results_true$v_indir==0,pch_name='v_indirect==0',
                              color=pgi_results_true$r_dir_indir,color_name='r_direct_indirect',
                              xlab='True alpha_delta',ylab='Estimated alpha_delta (noise=10)',
                              outfile='plots/alpha_delta_v10.pdf')

v100_alpha_delta = scatter_plot(pgi_results_true$alpha_delta,pgi_results_true$alpha_delta_SE,
                                pgi_results_v100$alpha_delta,pgi_results_v100$alpha_delta_SE,
                                pch=pgi_results_true$v_indir==0,pch_name='v_indirect==0',
                                color=pgi_results_true$r_dir_indir,color_name='r_direct_indirect',
                                xlab='True alpha_delta',ylab='Estimated alpha_delta (noise=100)',
                                outfile='plots/alpha_delta_v100.pdf')

v1_eta_delta = scatter_plot(pgi_results_true$v_eta_delta,pgi_results_true$v_eta_delta_SE,
                              pgi_results_v1$v_eta_delta,pgi_results_v1$v_eta_delta_SE,
                              pch=pgi_results_true$v_indir==0,pch_name='v_indirect==0',
                              color=pgi_results_true$r_dir_indir,color_name='r_direct_indirect',
                              xlab='True v_eta_delta',ylab='Estimated v_eta_delta (noise=1)',
                              outfile='plots/v_eta_delta_v1.pdf')

v10_eta_delta = scatter_plot(pgi_results_true$v_eta_delta,pgi_results_true$v_eta_delta_SE,
                            pgi_results_v10$v_eta_delta,pgi_results_v10$v_eta_delta_SE,
                            pch=pgi_results_true$v_indir==0,pch_name='v_indirect==0',
                            color=pgi_results_true$r_dir_indir,color_name='r_direct_indirect',
                            xlab='True v_eta_delta',ylab='Estimated v_eta_delta (noise=10)',
                            outfile='plots/v_eta_delta_v10.pdf')

v100_eta_delta = scatter_plot(pgi_results_true$v_eta_delta,pgi_results_true$v_eta_delta_SE,
                             pgi_results_v100$v_eta_delta,pgi_results_v100$v_eta_delta_SE,
                             pch=pgi_results_true$v_indir==0,pch_name='v_indirect==0',
                             color=pgi_results_true$r_dir_indir,color_name='r_direct_indirect',
                             xlab='True v_eta_delta',ylab='Estimated v_eta_delta (noise=100)',
                             outfile='plots/v_eta_delta_v100.pdf')


v1_r = scatter_plot(pgi_results_true$r,pgi_results_true$r_SE,
                              pgi_results_v1$r,pgi_results_v1$r_SE,
                              pch=as.factor(pgi_results_true$r_dir_indir),pch_name='r_direct_indir',
                              color=pgi_results_true$r_y,color_name='phenotype correlation',
                              xlab='True r',ylab='Estimated r (noise=1)',
                              outfile='plots/r_v1.pdf')

v10_r = scatter_plot(pgi_results_true$r,pgi_results_true$r_SE,
                    pgi_results_v10$r,pgi_results_v10$r_SE,
                    pch=as.factor(pgi_results_true$r_dir_indir),pch_name='r_direct_indir',
                    color=pgi_results_true$r_y,color_name='phenotype correlation',
                    xlab='True r',ylab='Estimated r (noise=10)',
                    outfile='plots/r_v10.pdf')

v100_r = scatter_plot(pgi_results_true$r,pgi_results_true$r_SE,
                     pgi_results_v100$r,pgi_results_v100$r_SE,
                     pch=as.factor(pgi_results_true$r_dir_indir),pch_name='r_direct_indir',
                     color=pgi_results_true$r_y,color_name='phenotype correlation',
                     xlab='True r',ylab='Estimated r (noise=100)',
                     outfile='plots/r_v100.pdf')


v1_h2_eq = scatter_plot(results$v_g_eq[results$r_y>0]/results$v_y_eq[results$r_y>0],0,
                    pgi_results_v1$h2_eq[results$r_y>0],pgi_results_v1$h2_eq_SE[results$r_y>0],
                    pch=results$v_indir[results$r_y>0]==0,pch_name='v_indir==0',
                    color=results$r_y[results$r_y>0],color_name='phenotype correlation',
                    xlab='True h2_eq',ylab='Estimated h2_eq (noise=1)',
                    outfile='plots/h2_eq_v1.pdf')

v10_h2_eq = scatter_plot(results$v_g_eq[results$r_y>0]/results$v_y_eq[results$r_y>0],0,
                        pgi_results_v10$h2_eq[results$r_y>0],pgi_results_v10$h2_eq_SE[results$r_y>0],
                        pch=results$v_indir[results$r_y>0]==0,pch_name='v_indir==0',
                        color=results$r_y[results$r_y>0],color_name='phenotype correlation',
                        xlab='True h2_eq',ylab='Estimated h2_eq (noise=10)',
                        outfile='plots/h2_eq_v10.pdf')

v100_h2_eq = scatter_plot(results$v_g_eq[results$r_y>0]/results$v_y_eq[results$r_y>0],0,
                        pgi_results_v100$h2_eq[results$r_y>0],pgi_results_v100$h2_eq_SE[results$r_y>0],
                        pch=results$v_indir[results$r_y>0]==0,pch_name='v_indir==0',
                        color=results$r_y[results$r_y>0],color_name='phenotype correlation',
                        xlab='True h2_eq',ylab='Estimated h2_eq (noise=100)',
                        outfile='plots/h2_eq_v100.pdf')



#################### Population Effect PGI results ##############################

pop_direct_r = scatter_plot(pgi_results_true$rk, pgi_results_true$rk_SE, 
                            pop_pgi_results_true$rk, pop_pgi_results_true$rk_SE,
                            pch=as.factor(pgi_results_true$r_dir_indir),pch_name='r_direct_indirect',
                            color=pgi_results_true$r_y,color_name='Phenotype correlation',
                            xlab='Direct effect PGI correlation',ylab='Population effect PGI correlation',
                            outfile='plots/direct_and_population_PGI_correlations.pdf'
                            )

pop_direct_alpha = scatter_plot(pgi_results_true$alpha, pgi_results_true$alpha_SE, 
                            pop_pgi_results_true$alpha, pop_pgi_results_true$alpha_SE,
                            pch=as.factor(pgi_results_true$r_dir_indir),pch_name='r_direct_indirect',
                            color=pgi_results_true$r_y,color_name='Phenotype correlation',
                            xlab='Direct effect PGI average NTC',ylab='Population effect PGI average NTC',
                            outfile='plots/direct_and_population_PGI_avg_NTCs.pdf'
)


pop_v1_r = scatter_plot(pgi_results_true$r,pgi_results_true$r_SE,
                    pop_pgi_results_v1$r,pop_pgi_results_v1$r_SE,
                    pch=as.factor(pgi_results_true$r_dir_indir),pch_name='r_direct_indir',
                    color=pgi_results_true$r_y,color_name='phenotype correlation',
                    xlab='True r',ylab='Estimated r (population effect PGI, noise=1)',
                    outfile='plots/r_pop_v1.pdf')

pop_v1_h2_eq = scatter_plot(results$v_g_eq/results$v_y_eq,0,
                        pop_pgi_results_v1$h2_eq,pop_pgi_results_v1$h2_eq_SE,
                        pch=as.factor(pgi_results_true$r_dir_indir),pch_name='r_direct_indir',
                        color=pgi_results_true$r_y,color_name='phenotype correlation',
                        xlab='True h2_eq',ylab='Estimated h2_eq (population effect PGI, noise=1)',
                        outfile='plots/h2_eq_pop_v1.pdf')

pop_v_eta_delta = scatter_plot(pgi_results_true$v_eta_delta,pgi_results_true$v_eta_delta_SE,
                            pop_pgi_results_v1$v_eta_delta,pop_pgi_results_v1$v_eta_delta_SE,
                            pch=as.factor(pgi_results_true$r_dir_indir),pch_name='r_direct_indir',
                            color=pgi_results_true$r_y,color_name='phenotype correlation',
                            xlab='Direct effect PGI v_eta_delta',ylab='Estimated v_eta_delta (population effect PGI, noise=1)',
                            outfile='plots/v_eta_delta_pop_v1.pdf')

