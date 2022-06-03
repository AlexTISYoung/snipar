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

# Remove outliers
convert_school = function(x){
    if (is.na(x)){return(NA)}
    if (x=="Mindre  n folkskola"){return(7)}
    if (x=="Folk- eller medborgarskola"){return(7)}
    if (x=="L roverk, mindre  n mellanskola"){return(10)}
    if (x=="Mellanskola eller 9- rig grundskola"){return(10)}
    if (x=="En del av gymnasium"){return(10)}
    if (x=="Student"){return(13)}
}

convert_school_1 = function(x){
    if (is.na(x)){return(NA)}
    if (x=="Ingen yrkesutbildning"){return(7)}
    if (x=="Kurser och arbetsplatsskolning"){return(7)}
    if (x=="Skolm ssiga studier h gst 2  r"){return(10)}
    if (x=="Skolm ssiga studier  ver 2  r"){return(10)}
    if (x=="Akademisk examen"){return(19)}
}

d$EA = NA
school = sapply(d$SCHOOL,convert_school)
school = cbind(school,sapply(d$SCHOOLY1,convert_school_1))
d$EA = apply(school,1,max,na.rm=T)
d$EA[d$EA<0 | d$AGE<30] = NA
d$FID = d$IID
d = cbind(d[,c('FID','IID')],d[,-match(c('PATIENT_FINAL','PATIENT2_FINAL','IID','FID'),dimnames(d)[[2]])])
write.table(d,'processed_traits_noadj.txt',quote=F,row.names=F)
