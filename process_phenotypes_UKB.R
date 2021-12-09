trait_names = c('Glucose',
'HDL cholesterol',
'Total cholesterol',
'Systolic blood pressure',
'Diastolic blood pressure',
'Self-rated health',
'FEV1',
'Ever-smoked',
'Cigarattes per day (current)',
'Cigarettes per day (former)',
'Drinks per week',
'height',
'BMI',
'Cognitive ability',
'Household income',
'Subjective well-being',
'Depressive symptoms',
'Neuroticism')

data_ids = c(rep('40305',3),rep('39604',7),'41490',
             rep('39604',7))

field_ids = c('30740.0.0','30760.0.0','30690.0.0',
              '4080.0.0','4079.0.0','2178.0.0','3063.0.0',
              '20116.0.0','3456.0.0','2887.0.0','20414.0.0',
              '50.0.0','21001.0.0','20016.0.0','738.0.0','4526.0.0','2050.0.0','20127.0.0')
field_ids = paste('f',field_ids,sep='.')

data_ids = rep('39604',2)
field_ids = c('129.0.0','130.0.0')
field_ids = paste('f',field_ids,sep='.')
trait_names = c('north','east')

data_ids = c('39604')
field_ids = c('54.0.0')
field_ids = paste('f',field_ids,sep='.')
trait_names = c('assessment_centre')

data_ids = c('39604')
field_ids = c('189.0.0')
field_ids = paste('f',field_ids,sep='.')
trait_names = c('townsend')

data_ids = c('39604')
field_ids = c('6138.0.0','6138.0.1','6138.0.2',
              '6138.0.3','6138.0.4','6138.0.5',
              '845.0.0')
field_ids = paste('f',field_ids,sep='.')
trait_names = c('qual_0','qual_1','qual_2','qual_3','qual_4','qual_5',
                'age_completed_edu')

traits = data.frame(names=trait_names,data_ids = data_ids,field_ids = field_ids)

#traits$qtrait = c(rep(T,5),F,T,F,T,T,F,T,T,T,F,F,F,T)
#traits$qtrait = c(T,T)
traits$qtrait=F

### Get raw traits ###

# Get sample QC
sqc = read.table('/disk/genetics2/ukb/orig/UKBv2/linking/ukb_sqc_v2_combined_header.txt',
                 header=T,stringsAsFactors = F)

# Filter sample 
sqc = sqc[sqc$het.missing.outliers==0 & sqc$putative.sex.chromosome.aneuploidy==0 &
            sqc$excess.relatives==0 & sqc$in.white.British.ancestry.subset==1,]

raw_traits = data.frame(matrix(NA,nrow=dim(sqc)[1],ncol=length(trait_names)))
dimnames(raw_traits)[[1]] = sqc$IID
dimnames(raw_traits)[[2]] = trait_names

datasets = unique(data_ids)
for (dataset in datasets){
  print(dataset)
  source(paste('ukb',dataset,'.r',sep=''))
  traits_d = traits[traits$data_ids==dataset,]
  trait_match = sapply(traits_d[,3],function(x) match(x,dimnames(bd)[[2]]))
  id_match = match(dimnames(raw_traits)[[1]],bd[,1])
  raw_traits[,as.character(traits_d[,1])] = as.matrix(bd[id_match,trait_match])
}

############# Convert EA
raw_traits$age_completed_edu = as.integer(raw_traits$age_completed_edu)-5
raw_traits$age_completed_edu[raw_traits$age_completed_edu<7] = NA
ea = matrix(NA,nrow=dim(raw_traits)[1],ncol=dim(raw_traits)[2])

for (j in 1:6){
ea[raw_traits[,j]=='NVQ or HND or HNC or equivalent' & !is.na(raw_traits[,j]),j] = raw_traits$age_completed_edu[raw_traits[,j]=='NVQ or HND or HNC or equivalent'  & !is.na(raw_traits[,j])]
ea[raw_traits[,j]=='College or University degree' & !is.na(raw_traits[,j]),j] = 20
ea[raw_traits[,j]=='A levels/AS levels or equivalent' & !is.na(raw_traits[,j]),j] = 13
ea[raw_traits[,j]=='O levels/GCSEs or equivalent' & !is.na(raw_traits[,j]),j] = 10
ea[raw_traits[,j]=='None of the above' & !is.na(raw_traits[,j]),j] = 7
ea[raw_traits[,j]=='CSEs or equivalent' & !is.na(raw_traits[,j]),j] = 10
ea[raw_traits[,j]=='Other professional qualifications eg: nursing, teaching' & !is.na(raw_traits[,j]),j] = 15
}

EA = apply(ea,1,max,na.rm=T)
EA[EA<0] = NA
EA = data.frame(IID=sqc$IID,EA=EA)
raw_traits = EA

north_east = read.table('north_east.fam',header=T)
assessment = read.table('assessment_centre.txt',header=T)
raw_traits$north = north_east[match(raw_traits$IID,north_east$IID),'north']
raw_traits$north2 = raw_traits$north^2
raw_traits$north3 = raw_traits$north^3
raw_traits$east = north_east[match(raw_traits$IID,north_east$IID),'east']
raw_traits$east2 = raw_traits$east^2
raw_traits$east3 = raw_traits$east^3
raw_traits$assessment = as.factor(assessment[match(raw_traits$IID,assessment$IID),2])
raw_traits$ea4 = ea4[match(raw_traits$IID,ea4$IID),3]

r = lm(ea$ea~sex+age+age2+age3+agesex+age2sex+age3sex+genotyping.array+as.matrix(pcs)+
         north*east+north2*east+east2*north+north3*east+east3*north+north2*east3+
         north3*east2+north3*east3+assessment*(north+north2+north3+east+east2+east3),data=raw_traits)

r = lm(ea4~as.matrix(pcs)+
         north*east+north2*east+east2*north+north3*east+east3*north+north2*east3+
         north3*east2+north3*east3,data=raw_traits)

raw_traits$ea4_resid = residuals(r)

ea$EA_agesexarraypcs_northeast_assessment = residuals(r)

write.table(ea,'EA_with_controls.fam',quote=F,row.names=F)
save.image('EA_with_controls.RData')

################

fert = read.table('fertility_traits.txt',header=T,stringsAsFactors = F)
unenthusiasm = read.table('unenthusiasm.txt',header=T)
raw_traits = data.frame(raw_traits,fert[match(dimnames(raw_traits)[[1]],fert[,1]),])
raw_traits = data.frame(raw_traits,unenthusiasm[match(dimnames(raw_traits)[[1]],dimnames(unenthusiasm)[[1]]),1])
dimnames(raw_traits)[[2]][dim(raw_traits)[2]]='unenthusiasm'

# add covariates
raw_traits = data.frame(raw_traits)
raw_traits$sex = sapply(sqc$Inferred.Gender,function(x) if (x=='M'){return(0)} 
                        else if (x=='F'){return(1)} else {return(NA)})
age = read.table('age.txt',header=T,stringsAsFactors = F)
raw_traits$age = age[match(dimnames(raw_traits)[[1]],age[,1]),2]
raw_traits$age2 = scale(raw_traits$age)^2
raw_traits$age3 = scale(raw_traits$age)^3
raw_traits$agesex= raw_traits$age*raw_traits$sex
raw_traits$age2sex= raw_traits$age2*raw_traits$sex
raw_traits$age3sex= raw_traits$age3*raw_traits$sex
raw_traits$genotyping.array = sqc$genotyping.array
raw_traits$Non_HDL = as.numeric(raw_traits[,'Total.cholesterol'])-
  as.numeric(raw_traits[,'HDL.cholesterol'])
raw_traits$cigarettes.per.day = as.numeric(raw_traits$Cigarattes.per.day..current.)
raw_traits$cigarettes.per.day[is.na(raw_traits$cigarettes.per.day)] = as.numeric(raw_traits$Cigarettes.per.day..former.[is.na(raw_traits$cigarettes.per.day)])
save(raw_traits,file='raw_traits.RData')
write.table(raw_traits,'raw_traits.txt',quote=F,sep='\t')
pcs = as.matrix(sqc[,grep('PC',dimnames(sqc)[[2]])])
save.image('raw_traits_and_sqc.RData')


################# Process traits ################
load('raw_traits_and_sqc.RData')
processed_traits = data.frame(IID=as.integer(dimnames(raw_traits)[[1]]))

# quantitative traits
get_residuals = function(y,raw_traits,pcs,binary=F){
  options(na.action = 'na.exclude')
  if (binary){
    r = glm(y~sex+age+age2+age3+agesex+age2sex+age3sex+genotyping.array+as.matrix(pcs),
            data=raw_traits,family='binomial')
  } else{
    r = lm(y~sex+age+agesex+age2sex+age3sex+genotyping.array+as.matrix(pcs),
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

qt_transform = function(y,raw_traits,pcs,remove_upper=F,remove_lower=F,binary=F){
  y = as.numeric(y)
  y[y<0] = NA
  if (remove_lower){
    y[y<quantile(y,0.001,na.rm=T)] = NA
  }
  if (remove_upper){
    y[y>quantile(y,0.999,na.rm=T)] = NA
  }
  return(get_residuals(y,raw_traits,pcs,binary))
}


convert_self_rated_health = function(x){
  if (is.na(x)){return(NA)}
  if (x=='Do not know'){return(NA)}
  if (x=="Prefer not to answer"){return(NA)}
  if (x=='Poor'){return(0)}
  if (x=='Fair'){return(1)}
  if (x=='Good'){return(2)}
  if (x=='Excellent'){return(3)}
}

convert_ever_smoked = function(x){
  if (is.na(x)){return(NA)}
  else if (x=='Prefer not to answer'){return(NA)}
  else if (x=='Never'){return(0)}
  else {return(1)}
}

convert_drinks = function(x){
  if (is.na(x)){return(NA)}
  if (x=='Prefer not to answer'){return(NA)}
  if (x=='Never' | x=='Monthly or less'){return(0)}
  if (x=='2 to 3 times a week'){return(2.5)}
  if (x=='2 to 4 times a month'){return(0.75)}
  if (x=='4 or more times a week'){return(4)}
}

convert_income = function(x){
  if (is.na(x)){return(NA)}
  if (x=='Prefer not to answer'){return(NA)}
  if (x=='Do not know'){return(NA)}
  if (x=="Less than 18,000"){return(log(18000))}
  if (x=="18,000 to 30,999"){return(mean(c(log(18000),log(30999))))}
  if (x=="31,000 to 51,999"){return(mean(c(log(31000),log(51999))))}
  if (x=="52,000 to 100,000"){return(mean(c(log(52000),log(100000))))}
  if (x=='Greater than 100,000'){return(log(100000))}
}

happiness_convert = function(x){
  if (is.na(x)){return(NA)}
  if (x=='Prefer not to answer'){return(NA)}
  if (x=='Do not know'){return(NA)}
  if (x=='Extremely unhappy'){return(1)}
  if (x=="Very unhappy"){return(2)}
  if (x=='Moderately unhappy'){return(3)}
  if (x=='Moderately happy'){return(4)}
  if (x=='Very happy'){return(5)}
  if (x=='Extremely happy'){return(6)}
}


convert_mood = function(x){
  if (is.na(x)){return(NA)}
  if (x=='Prefer not to answer'){return(NA)}
  if (x=='Do not know'){return(NA)}
  if (x=='Not at all'){return(1)}
  if (x=='Several days'){return(2)}
  if (x=='More than half the days'){return(3)}
  if (x=='Nearly every day'){return(4)}
}

#processed_traits$north = qt_transform(raw_traits$north,raw_traits,pcs,remove_upper=F,remove_lower=F)
#processed_traits$east = qt_transform(raw_traits$east,raw_traits,pcs,remove_upper=F,remove_lower=F)


processed_traits$Glucose = qt_transform(raw_traits$Glucose,raw_traits,pcs,remove_upper=T,remove_lower=T)
processed_traits$Non_HDL = qt_transform(raw_traits$Non_HDL,raw_traits,pcs,remove_upper=T,remove_lower=T)
processed_traits$HDL = qt_transform(raw_traits$HDL.cholesterol,raw_traits,pcs,remove_upper=T,remove_lower=T)
processed_traits$SBP = qt_transform(raw_traits$Systolic.blood.pressure,raw_traits,pcs,remove_upper=T,remove_lower=T)
processed_traits$DBP = qt_transform(raw_traits$Diastolic.blood.pressure,raw_traits,pcs,remove_upper=T,remove_lower=T)
processed_traits$self.rated.health = qt_transform(sapply(raw_traits$Self.rated.health,convert_self_rated_health),
                                                  raw_traits,pcs)
processed_traits$FEV1 = qt_transform(raw_traits$FEV1,raw_traits,pcs,remove_upper=T,remove_lower=T)
processed_traits$ever.smoked = qt_transform(sapply(raw_traits$Ever.smoked,convert_ever_smoked),
                                            raw_traits,pcs,binary=T)
processed_traits$cigarettes.per.day = qt_transform(raw_traits$cigarettes.per.day,raw_traits,pcs)
processed_traits$drinks.per.week = qt_transform(sapply(raw_traits$Drinks.per.week,convert_drinks),
                                                  raw_traits,pcs)
processed_traits$height = qt_transform(raw_traits$height,raw_traits,pcs,remove_upper=T,remove_lower=T)
processed_traits$BMI = qt_transform(raw_traits$BMI,raw_traits,pcs,remove_upper=T,remove_lower=T)
processed_traits$Cognitive.ability = qt_transform(raw_traits$Cognitive.ability,raw_traits,pcs,remove_upper=T,remove_lower=T)
processed_traits$Neuroticism = qt_transform(raw_traits$Neuroticism,raw_traits,pcs,remove_upper=T,remove_lower=T)
processed_traits$AAFB = qt_transform(raw_traits$AAFB,raw_traits,pcs)
processed_traits$NC_M = qt_transform(raw_traits$NC_M,raw_traits,pcs,remove_upper=T)
processed_traits$NC_F = qt_transform(raw_traits$NC_F,raw_traits,pcs,remove_upper=T)
processed_traits$household.income = qt_transform(sapply(raw_traits$Household.income,convert_income),raw_traits,pcs)
processed_traits$subjective.well.being = qt_transform(sapply(raw_traits$Subjective.well.being,happiness_convert),raw_traits,pcs)
write.table(data.frame(FID=processed_traits$IID,processed_traits),'processed_traits.fam',quote=F,row.names=F)



