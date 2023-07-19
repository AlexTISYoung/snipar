setwd('/disk/genetics/sibling_consortium/GS20k/alextisyoung/')
#### Read raw traits ####
raw_traits = read.csv('phenotypes/agesex.csv',header=T)
raw_traits = rbind(raw_traits,read.csv('../phenotype/Agesex.csv',header=T))

files1 = c('biochemistry.csv',  'body.csv',	
              'cognitive.csv',  'education.csv',  'household.csv',  
              	'smoking.csv',  'spirometry.csv','BPHR.csv')
files2 = c(     "Biochemistry.csv", "Body.csv", "Cognition.csv",        
            "Education.csv" ,"Household.csv","Smoking.csv","Spirometry.csv",'BPHR.csv')

for (i in 1:length(files1)){
  p = read.csv(paste('phenotypes',files1[i],sep='/'),header=T)
  p2 = read.csv(paste('../phenotype',files2[i],sep='/'),header=T)
  colmatch = match(dimnames(p)[[2]],dimnames(p2)[[2]])
  p = rbind(p[,!is.na(colmatch)],p2[,colmatch[!is.na(colmatch)]])
  raw_traits = data.frame(raw_traits,p[match(raw_traits[,1],p[,1]),-1])
}     

raw_traits$age2 = raw_traits$age^2
raw_traits$age3 = raw_traits$age^3
raw_traits$sex = sapply(raw_traits$sex,
                        function(x) if (x=='M'){return(1)} else if (x=='F'){return(0)} else {return(NA)})
raw_traits$agesex = raw_traits$age*raw_traits$sex
raw_traits$age2sex = raw_traits$age2*raw_traits$sex
raw_traits$age3sex = raw_traits$age3*raw_traits$sex

## Remove non-European ancestry individuals
gs_pops = read.table('HM3/pca/GS_population_assignment.txt',header=T)
gs_EUR = gs_pops[gs_pops$pop=='EUR',1]
raw_traits = raw_traits[raw_traits$id%in%gs_EUR,]

################# Process traits ################
processed_traits = data.frame(IID=raw_traits$id)

qt_transform = function(y,raw_traits,remove_upper=F,remove_lower=F,remove_neg=T){
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
  sex=raw_traits$sex
  non_NA_sex = sex[!is.na(y)]
  if (sum(non_NA_sex==0)>0){
    y[sex==0] = scale(y[sex==0],center=T)}
  if (sum(non_NA_sex==1)>0){
    y[sex==1] = scale(y[sex==1],center=T)}
  return(y)
}

## EA
years_convert = function(y){
  if (is.na(y)){return(NA)}
  if (y==0){return(1)}
  if (y==1){return(7)}
  if (y==2){return(7)}
  if (y==3){return(10)}
  if (y==4){return(13)}
  if (y==5){return(15)}
  if (y==6){return(19)}
  if (y==7){return(19)}
  if (y==8){return(22)}
  if (y==9){return(22)}
  if (y==10){return(22)}
}

years_convert_mid = function(y){
  if (is.na(y)){return(NA)}
  if (y==0){return(0)}
  if (y==1){return(2.5)}
  if (y==2){return(7)}
  if (y==3){return(10.5)}
  if (y==4){return(12.5)}
  if (y==5){return(14.5)}
  if (y==6){return(16.5)}
  if (y==7){return(18.5)}
  if (y==8){return(20.5)}
  if (y==9){return(22.5)}
  if (y==10){return(24)}
}

quals_convert = function(qual,year){
  if (is.na(qual)){return(NA)}
  if (qual==1){return(19)}
  if (qual==2){return(13)}
  if (qual==3){return(NA)}
  if (qual==4){return(13)}
  if (qual==5){return(10)}
  if (qual==6){return(10)}
  if (qual==7){return(10)}
  if (qual==8){return(NA)}
  if (qual==9){return(7)}
}

## Ever smoke
# 1 - Yes, currently smoke, 2 - Yes, but stopped within past 12 months, 3 - Yes, but stopped more than 12 months ago, 4 - No, never smoked
convert_ever_smoked = function(y){
  if (is.na(y)){return(NA)}
  if (y==4){return(0)}
  else {return(1)}
}

processed_traits$Glucose = qt_transform(raw_traits$Glucose,raw_traits,remove_upper=T,remove_lower=T)
processed_traits$Non_HDL = qt_transform(raw_traits$Total_cholesterol-raw_traits$HDL_cholesterol,raw_traits,remove_upper=T,remove_lower=T)
processed_traits$HDL = qt_transform(raw_traits$HDL_cholesterol,raw_traits,remove_upper=T,remove_lower=T)
processed_traits$height = qt_transform(raw_traits$height,raw_traits,remove_upper=T,remove_lower=T)
# FEV1
options(na.action='na.exclude')
fev_height_reg = lm(raw_traits$FEV_1~processed_traits$height*raw_traits$sex)
processed_traits$FEV1 = qt_transform(residuals(fev_height_reg),raw_traits,remove_upper=T,remove_lower=T,remove_neg=F)
#
processed_traits$BMI = qt_transform(raw_traits$bmi,raw_traits,remove_upper=T,remove_lower=T)
processed_traits$ever.smoked = sapply(raw_traits$ever_smoke,convert_ever_smoked)
raw_traits$packs_day[raw_traits$packs_day==0] = NA
processed_traits$cigarettes.per.day = qt_transform(raw_traits$packs_day,raw_traits)

# Cognitive PCA
cog = raw_traits[,c('verbal_total','logical_mem_1','logical_mem_2','digit_symbol')]
cog[,2] = cog[,2]+cog[,3]
cog = cog[,-3]
dimnames(cog)[[1]] = raw_traits[,1]
cogna = is.na(cog)
cog = cog[apply(cogna,1,sum)==0,]
cogpcs = prcomp(cog,center = TRUE, scale = TRUE)
raw_traits$cog= cogpcs$x[match(raw_traits[,1],dimnames(cog)[[1]]),1]
processed_traits$cog = qt_transform(raw_traits$cog,raw_traits,remove_neg=F)
#
processed_traits$vocab = qt_transform(raw_traits$vocabulary,raw_traits,remove_lower=T)
processed_traits$SBP = qt_transform(raw_traits$avg_sys,raw_traits,remove_upper=T,remove_lower=T)
processed_traits$DBP = qt_transform(raw_traits$avg_dia,raw_traits,remove_upper=T,remove_lower=T)
processed_traits$neuroticism = qt_transform(raw_traits$eysenck_N ,raw_traits)
# New EA based on qualifications
processed_traits$EA_years = sapply(raw_traits$years,years_convert)
processed_traits$EA_years_mid = sapply(raw_traits$years,years_convert_mid)
processed_traits$EA_quals = sapply(raw_traits$qualification,quals_convert)
processed_traits$EA_quals[is.na(processed_traits$EA_quals)] = processed_traits$EA_years[is.na(processed_traits$EA_quals)]

# Match with geno ids
fam = read.table('grandpar/haplotypes/bedfiles/autosome.fam')
fam$IID = sapply(fam[,2],function(x) strsplit(x,'_')[[1]][2])
processed_traits$IID = fam[match(processed_traits$IID,fam$IID),2]
processed_traits$FID = processed_traits$IID
processed_traits = processed_traits[,c(dim(processed_traits)[2],1:(dim(processed_traits)[2]-1))]
processed_traits[raw_traits$age<30,c('EA_years','EA_years_mid','EA_quals')] = NA
# Write
write.table(processed_traits,'grandpar/processed_traits_noadj.txt',quote=F,row.names=F)

# Agesex
agesex = data.frame(FID=processed_traits$FID,IID=processed_traits$IID,
                    age=raw_traits$age,
                    sex=sapply(raw_traits$sex,function(x) if (x==0){return('F')} else if (x==1){return('M')} else {return(NA)}))
write.table(agesex,'grandpar/agesex.txt',quote=F,row.names=F)

# Write correlation matrix
write.csv(cor(processed_traits[,c('cog','vocab','EA_years','EA_years_mid','EA_quals')],use='pairwise.complete.obs'),'grandpar/EA_cog_vocab_corrs.csv',quote=F)