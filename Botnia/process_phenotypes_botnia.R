qt_transform = function(y,sex,remove_upper=F,remove_lower=F,remove_neg=T){
  y = as.numeric(y)
  if (remove_neg){
  y[y<0] = NA}
  if (remove_lower){
    y[y<quantile(y,0.001,na.rm=T)] = NA
  }
  if (remove_upper){
    y[y>quantile(y,0.999,na.rm=T)] = NA
  }
  # Scale within sex
  non_NA_sex = sex[!is.na(y)]
  if (sum(non_NA_sex=='Male')>0){
    y[sex=='Male'] = scale(y[sex=='Male'],center=T)}
  if (sum(non_NA_sex=='Female')>0){
    y[sex=='Female'] = scale(y[sex=='Female'],center=T)}
  return(y)
}

library(foreign)
setwd('/gpfs/gpfs0/Active_Projects/BOTNIA_AYoung_analysis/Private/phenotypes/')

# Read from SPSS produced by Mikko
d = read.spss('Botex_selected_variables_2021-07-09_FirstVisit.sav', to.data.frame=TRUE)
d$PROF = NULL
d$PROF0 = NULL
d$PROF1 = NULL

# Find IIDs
fam = read.table('../haplotypes/bedfiles/autosome.fam')
d$IID = NA
d$IID[d$PATIENT_FINAL%in%fam[,2]] = d$PATIENT_FINAL[d$PATIENT_FINAL%in%fam[,2]]
d$IID[!d$PATIENT_FINAL%in%fam[,2] & d$PATIENT2_FINAL%in%fam[,2]] = d$PATIENT2_FINAL[!d$PATIENT_FINAL%in%fam[,2] & d$PATIENT2_FINAL%in%fam[,2]]
d = d[!is.na(d$IID),] 

###### EA CODING ####
options(na.action='na.exclude')
d$SCHOOLY2[d$SCHOOLY2>30 | d$AGE<30] = NA
r_agesex = lm(SCHOOLY2~AGE*SEX+I(AGE^2)*SEX+I(AGE^3)*SEX,data=d)
d$resid_years = mean(d$SCHOOLY2,na.rm=T)+residuals(r_agesex)

school_map = sapply(levels(d$SCHOOL),function(x) median(d$resid_years[d$SCHOOL==x],na.rm=T))
school1_map = sapply(levels(d$SCHOOLY1),function(x) median(d$resid_years[d$SCHOOLY1==x],na.rm=T))
school1_map['Ingen yrkesutbildning'] = NA
school1_map['Kurser och arbetsplatsskolning'] = NA

d$EA = NA
school = sapply(d$SCHOOL,function(x) if (is.na(x)){return(NA)} else {return(school_map[x])})
school = cbind(school,sapply(d$SCHOOLY1,function(x) if (is.na(x)){return(NA)} else {return(school1_map[x])}))
d$EA = apply(school,1,max,na.rm=T)
d$EA[d$EA<0 | d$AGE<30] = NA
d$resid_years = NULL
d$FID = d$IID
d = cbind(d[,c('FID','IID')],d[,-match(c('PATIENT_FINAL','PATIENT2_FINAL','IID','FID'),dimnames(d)[[2]])])

###### PHENOTYPES ####
phenotypes = data.frame(d[,c('FID','IID','EA')])
phenotypes$glucose = qt_transform(d$PFASTGLUC,d$SEX,remove_upper=T)
phenotypes$HDL = qt_transform(d$HDLCHOL,d$SEX,remove_upper=T)
phenotypes$non_HDL = qt_transform(d$CHOL-d$HDLCHOL,d$SEX,remove_upper=T)
phenotypes$height = qt_transform(d$HEIGHT,d$SEX,remove_lower=T,remove_upper=T)
phenotypes$BMI = qt_transform(d$BMI,d$SEX,remove_lower=T,remove_upper=T)
phenotypes$DBP = qt_transform(d$DIAST,d$SEX,remove_lower=T,remove_upper=T)
phenotypes$SBP = qt_transform(d$SYST,d$SEX,remove_lower=T,remove_upper=T)

write.table(phenotypes,'processed_traits_noadj.txt',quote=F,row.names=F)

###### COVARIATES ######
pcs = read.table('/ludc/Active_Projects/BOTNIA_AYoung_analysis/Private/haplotypes/bedfiles/pca/PCs.txt',header=T,row.names=1)
covariates = d[c('FID','IID','AGE')]
covariates$sex = sapply(d$SEX,function(x) if (x=='Male'){return(0)} else if (x=='Female'){return(1)})
covariates$age_sex = scale(covariates$AGE)*scale(covariates$sex)
covariates$age2_sex = scale(covariates$AGE)^2*scale(covariates$sex)
covariates$age3_sex = scale(covariates$AGE)^3*scale(covariates$sex)
covariates$PC=pcs[match(d$IID,dimnames(pcs)[[1]]),]

write.table(covariates,'covariates.txt',quote=F,row.names=F)