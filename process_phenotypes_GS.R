###############################################################################################
# Process phenotypes in generation scotland, adjusting for covariates and normalizing within sex 
###############################################################################################

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

pcs = as.matrix(read.table('GS_PCs.eigenvec'))[,-1]
pcs = pcs[match(traits$IID,pcs[,1]),-1]

#save.image('phenotypes/raw_traits.RData')

################# Process traits ################
#load('phenotypes/raw_traits.RData')
processed_traits = data.frame(IID=raw_traits$id)

# quantitative traits
get_residuals = function(y,raw_traits,pcs,binary=F){
  options(na.action = 'na.exclude')
  if (binary){
    r = glm(y~sex+age+age2+age3+agesex+age2sex+age3sex+as.matrix(pcs),
            data=raw_traits,family='binomial')
  } else{
    r = lm(y~sex+age+agesex+age2sex+age3sex+as.matrix(pcs),
            data=raw_traits)
  }
  resid = residuals(r)
  non_NA_sex = raw_traits$sex[!is.na(y)]
  if (sum(non_NA_sex==0)>0){
    resid[raw_traits$sex==0] = scale(resid[raw_traits$sex==0])}
  if (sum(non_NA_sex==1)>0){
    resid[raw_traits$sex==1] = scale(resid[raw_traits$sex==1])}
  return(resid)
}

qt_transform = function(y,raw_traits,pcs,remove_upper=F,remove_lower=F,binary=F,remove_neg=T){
  y = as.numeric(y)
  if (remove_neg){
  y[y<0] = NA}
  if (remove_lower){
    y[y<quantile(y,0.001,na.rm=T)] = NA
  }
  if (remove_upper){
    y[y>quantile(y,0.999,na.rm=T)] = NA
  }
  return(get_residuals(y,raw_traits,pcs,binary))
}

ea_convert = function(y){
  if (is.na(y)){return(NA)}
  if (y==0 | y==1){return(1)}
  if (y==2){return(7)}
  if (y==3){return(10)}
  if (y==4){return(13)}
  if (y==5){return(15)}
  if (y==6){return(15)}
  if (y==7){return(19)}
  if (y==8){return(19)}
  if (y==9){return(22)}
  if (y==10){return(22)}
}

processed_traits$Glucose = qt_transform(raw_traits$Glucose,raw_traits,pcs,remove_upper=T,remove_lower=T)
processed_traits$Non_HDL = qt_transform(raw_traits$Total_cholesterol-raw_traits$HDL_cholesterol,raw_traits,pcs,remove_upper=T,remove_lower=T)
processed_traits$HDL = qt_transform(raw_traits$HDL_cholesterol,raw_traits,pcs,remove_upper=T,remove_lower=T)
processed_traits$FEV1 = qt_transform(raw_traits$FEV_1,raw_traits,pcs,remove_upper=T,remove_lower=T)
processed_traits$height = qt_transform(raw_traits$height,raw_traits,pcs,remove_upper=T,remove_lower=T)
processed_traits$BMI = qt_transform(raw_traits$bmi,raw_traits,pcs,remove_upper=T,remove_lower=T)
raw_traits$packs_day[raw_traits$packs_day==0] = NA
processed_traits$cigarettes.per.day = qt_transform(raw_traits$packs_day,raw_traits,pcs)
# Cognitive PCA
cog = raw_traits[,c('verbal_total','logical_mem_1','logical_mem_2','digit_symbol')]
cog[,2] = cog[,2]+cog[,3]
cog = cog[,-3]
dimnames(cog)[[1]] = raw_traits[,1]
cogna = is.na(cog)
cog = cog[apply(cogna,1,sum)==0,]
cogpcs = prcomp(cog,center = TRUE, scale = TRUE)
raw_traits$cog= cogpcs$x[match(raw_traits[,1],dimnames(cog)[[1]]),1]
processed_traits$cog = qt_transform(raw_traits$cog,raw_traits,pcs,remove_neg=F)
processed_traits$vocab = qt_transform(raw_traits$vocabulary,raw_traits,pcs)
processed_traits$EA = qt_transform(raw_traits$years,raw_traits,pcs)
processed_traits$SBP = qt_transform(raw_traits$avg_sys,raw_traits,pcs,remove_upper=T,remove_lower=T)
processed_traits$DBP = qt_transform(raw_traits$avg_dia,raw_traits,pcs,remove_upper=T,remove_lower=T)
processed_traits$neuroticism = qt_transform(raw_traits$eysenck_N ,raw_traits,pcs)
processed_traits$EA_years = qt_transform(sapply(raw_traits$years,ea_convert),raw_traits,pcs)
# Scale
ids = read.table('sibs_and_trios.txt',stringsAsFactors = F)[,1]
sfs = apply(processed_traits[processed_traits$IID%in%ids,-c(1:2)],2,sd,na.rm=T)
processed_traits[,-c(1:2)] = processed_traits[,-c(1:2)]/sfs
# FIDs
ped = read.table('pedigree.txt',header=T,stringsAsFactors = F)
processed_traits$FID = ped[match(processed_traits$IID,ped$IID),'FID']
processed_traits$FID[is.na(processed_traits$FID)] = processed_traits$IID[is.na(processed_traits$FID)]
processed_traits = processed_traits[order(processed_traits$FID),]
processed_traits = processed_traits[,c(dim(processed_traits)[2],1:(dim(processed_traits)[2]-1))]
write.table(processed_traits,'processed_traits.fam',quote=F,row.names=F)


