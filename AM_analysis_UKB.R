########################################################################
# Script to perform assortative mating analyses using PGIs in UK Biobank
########################################################################
################################# Input files ###############################################
# RData containing raw UK Biobank traits, sample QC information, genetic principal components
# Output of process_phenotypes_UKB.R
raw_traits_and_sqc = 'raw_traits_and_sqc.RData'
# Processed phenotypes afer adjusting for covariates etc
# Output of process_phenotypes_UKB.R
processed_traits = 'processed_traits.fam'
# File containing north and east birth coordinates of UKB participants, with columns ID, north coordinate, east coordinate
north_east_file = 'north_east.fam'
# File containing UK Biobank assessment centre of UKB participants
assessment_centre_file = 'assessment_centre.txt'
# Pedigree file for UK Biobank participants
pedfile = 'ukb.ped'
# PGI files for EA, cognitive ability, height, and BMI
PGI_files = c('EA4_PGS_2.fam','cog_PGS.fam','height_PGS_2.fam','BMI_PGS_2.fam')
##########################################################################################

compute_cor = function(bpg_ped,traits,trait_name){
  pgs_bpg = matrix(NA,nrow=dim(bpg_ped),ncol=2)
  pgs_bpg[,1] = traits[match(bpg_ped[,3],traits$IID),trait_name]
  pgs_bpg[,2] = traits[match(bpg_ped[,4],traits$IID),trait_name]
  return(cor.test(pgs_bpg[,1],pgs_bpg[,2]))
}

# Read traits
load(raw_traits_and_sqc)
processed = read.table(processed_traits,header=T)

# Birth coordinates and assessment centre
north_east = read.table(north_east_file,header=T)
north_east$north.east = north_east$north*north_east$east

# Assessment centre
assess = read.table(assessment_centre_file,header=T)

# Pedigree
ped = read.table(pedfile,header=T)
bpg_ped = ped[ped[,3]%in%ped[,2] & ped[,4]%in%ped[,2],]

# Filter
in_bpg = raw_traits$IID%in%c(bpg_ped[,3],bpg_ped[,4])
traits = raw_traits[in_bpg,]
pcs = pcs[in_bpg,]
traits$EA = processed[match(traits$IID,processed$IID),'EA']

trait_names = c('EA','Cognitive.ability','height','BMI')
pgs_names = PGI_files
reg_names = c('EA','Cognitive.ability','height','BMI')

# Find unique parent-pairs
parent_pairs = apply(bpg_ped[,3:4],1,paste,collapse=',')
unique_pairs = unique(parent_pairs)
parent_pairs = t(sapply(unique_pairs,function(x) strsplit(x,',')[[1]]))

results = data.frame(trait=trait_names,phenotypic_correlation=NA,
                     phenotypic_correlation_lower = NA,
                     phenotypic_correlation_upper = NA,
                     trait_PGI_R2 = NA,
                     r_p = NA,
                     r_p_lower = NA,
                     r_p_upper = NA,
                     r_m = NA,
                     r_m_lower = NA,
                     r_m_upper = NA,
                     phenotypic_assortment_prediction=NA,
                     prediction_SE = NA,
                     PGI_r = NA,
                     PGI_r_lower = NA,
                     PGI_r_upper = NA,
                     PGI_trait_r = NA,
                     PGI_trait_r_lower = NA,
                     PGI_trait_r_upper = NA,
                     PGI_trait_spouse_r = NA,
                     PGI_trait_spouse_r_lower = NA,
                     PGI_trait_spouse_r_upper = NA,
                     PGI_trait_PCs_r = NA,
                     PGI_trait_PCs_r_lower = NA,
                     PGI_trait_PCs_r_upper = NA,
                     PGI_trait_PCs_coords_r = NA,
                     PGI_trait_PCs_coords_r_lower = NA,
                     PGI_trait_PCs_coords_r_upper = NA,
                     N=NA)

options(na.action = 'na.exclude')
for (i in 1:4){
  print(trait_names[i])
  # Phenotypic correlation
  spouses = matrix(NA,nrow=dim(parent_pairs)[1],2)
  spouses[,1] = as.numeric(traits[match(parent_pairs[,1],traits$IID),reg_names[i]])
  spouses[,2] = as.numeric(traits[match(parent_pairs[,2],traits$IID),reg_names[i]])
  spouse_corr = cor.test(spouses[,1],spouses[,2],use='pairwise.complete.obs')
  results[i,'phenotypic_correlation'] = spouse_corr$estimate
  results[i,'phenotypic_correlation_upper'] = spouse_corr$conf.int[2]
  results[i,'phenotypic_correlation_lower'] = spouse_corr$conf.int[1]
  results[i,'N']= sum(!is.na(spouses[,1]) & !is.na(spouses[,2]))
  # Read PGS
  pgs = read.table(pgs_names[i],header=T,stringsAsFactors=F)
  pgs_bpg = matrix(NA,nrow=dim(parent_pairs),ncol=2)
  pgs_bpg[,1] = pgs[match(parent_pairs[,1],pgs[,2]),3]
  pgs_bpg[,2] = pgs[match(parent_pairs[,2],pgs[,2]),3]
  # Unadjusted
  pgs_cor = cor.test(pgs_bpg[,1],pgs_bpg[,2])
  results[i,'PGI_r'] = pgs_cor$estimate
  results[i,'PGI_r_lower'] = pgs_cor$conf.int[1]
  results[i,'PGI_r_upper'] = pgs_cor$conf.int[2]
  # Trait adjusted
  father_trait = as.numeric(traits[match(parent_pairs[,1],traits$IID),reg_names[i]])
  mother_trait = as.numeric(traits[match(parent_pairs[,2],traits$IID),reg_names[i]])
  rp = lm(pgs_bpg[,1]~father_trait+I(father_trait^2)+I(father_trait^3))
  rm = lm(pgs_bpg[,2]~mother_trait+I(mother_trait^2)+I(mother_trait^3))
  results[i,'trait_PGI_R2'] = sqrt(summary(rp)$r.square*summary(rm)$r.square)
  r_p = cor.test(pgs_bpg[,1],as.numeric(traits[match(parent_pairs[,1],traits$IID),reg_names[i]]))
  r_m = cor.test(pgs_bpg[,2],as.numeric(traits[match(parent_pairs[,2],traits$IID),reg_names[i]]))
  results[i,'r_p'] = r_p$estimate
  results[i,c('r_p_lower','r_p_upper')] = r_p$conf.int
  results[i,'r_m'] = r_m$estimate
  results[i,c('r_m_lower','r_m_upper')] = r_m$conf.int
  pgs_cor = cor.test(residuals(rp),residuals(rm),use='pairwise.complete.obs')
  results[i,'PGI_trait_r'] = pgs_cor$estimate
  results[i,'PGI_trait_r_lower'] = pgs_cor$conf.int[1]
  results[i,'PGI_trait_r_upper'] = pgs_cor$conf.int[2]
  # Spouse trait adjusted 
  rp = lm(pgs_bpg[,1]~father_trait+I(father_trait^2)+I(father_trait^3)+
            as.numeric(traits[match(parent_pairs[,2],traits$IID),reg_names[i]]))
  rm = lm(pgs_bpg[,2]~mother_trait+I(mother_trait^2)+I(mother_trait^3)+
            as.numeric(traits[match(parent_pairs[,2],traits$IID),reg_names[i]]))
  pgs_cor = cor.test(residuals(rp),residuals(rm),use='pairwise.complete.obs')
  results[i,'PGI_trait_spouse_r'] = pgs_cor$estimate
  results[i,'PGI_trait_spouse_r_lower'] = pgs_cor$conf.int[1]
  results[i,'PGI_trait_spouse_r_upper'] = pgs_cor$conf.int[2]
  # Trait and PC adjusted
  rp = lm(pgs_bpg[,1]~as.numeric(traits[match(parent_pairs[,1],traits$IID),reg_names[i]])+
            as.matrix(pcs[match(parent_pairs[,1],traits$IID),]))
  rm = lm(pgs_bpg[,2]~mother_trait+I(mother_trait^2)+I(mother_trait^3)+
            as.matrix(pcs[match(parent_pairs[,2],traits$IID),]))
  pgs_cor = cor.test(residuals(rp),residuals(rm),use='pairwise.complete.obs')
  results[i,'PGI_trait_PCs_r'] = pgs_cor$estimate
  results[i,'PGI_trait_PCs_r_lower'] = pgs_cor$conf.int[1]
  results[i,'PGI_trait_PCs_r_upper'] = pgs_cor$conf.int[2]
  # Trait, coords, assessment centre
  rp = lm(pgs_bpg[,1]~father_trait+I(father_trait^2)+I(father_trait^3)+
            as.matrix(pcs[match(parent_pairs[,1],traits$IID),])+
            as.matrix(north_east[match(parent_pairs[,1],north_east[,2]),-c(1:2)])+
            factor(assess[match(parent_pairs[,1],assess[,1]),2]))
  rm = lm(pgs_bpg[,2]~mother_trait+I(mother_trait^2)+I(mother_trait^3)+
            as.matrix(pcs[match(parent_pairs[,2],traits$IID),])+
            as.matrix(north_east[match(parent_pairs[,2],north_east[,2]),-c(1:2)])+
            factor(assess[match(parent_pairs[,2],assess[,1]),2]))
  pgs_cor = cor.test(residuals(rp),residuals(rm),use='pairwise.complete.obs')
  results[i,'PGI_trait_PCs_coords_r'] = pgs_cor$estimate
  results[i,'PGI_trait_PCs_coords_r_lower'] = pgs_cor$conf.int[1]
  results[i,'PGI_trait_PCs_coords_r_upper'] = pgs_cor$conf.int[2]
  ## Prediction SE
  # Boostrap standard error
  not_NA = !is.na(spouses[,1]) & !is.na(spouses[,2])
  spouses = spouses[not_NA,]
  pgs_not_NA = pgs_bpg[not_NA,]
  n_bootstrap = 1000
  bootstrap_sample = rep(NA,n_bootstrap)
  for (j in 1:n_bootstrap){
    bsample = sample(1:dim(spouses)[1],size=dim(spouses)[1],replace=T)
    rp = lm(pgs_not_NA[bsample,1]~spouses[bsample,1])
    rm = lm(pgs_not_NA[bsample,2]~spouses[bsample,2])
    spouse_corr = cor(spouses[bsample,1],spouses[bsample,2],use='pairwise.complete.obs')
    bootstrap_sample[j] = sqrt(summary(rp)$r.square*summary(rm)$r.square)*spouse_corr
  }
  results[i,'prediction_SE'] = sqrt(var(bootstrap_sample))
}
results$phenotypic_assortment_prediction = results$phenotypic_correlation*results$trait_PGI_R2

write.csv(results,'AM_PGI_results_nonlinear.csv',quote=F,row.names=F)






