require(lme4)
binary_outcomes = c('depression','ADHD',
                    'hypertension','alcohol_use_disorder',
                    ,'ever_smoker')


phenotypes = read.table('~/phenotypes/processed_traits.txt',header=T)
covariates = read.table('~/phenotypes/covariates_reduced.txt',header=T)
covariates = covariates[,c('FID','IID','age','sex','age_sex','PC_3','PC_4')]
covariates[,-c(1:2)] = scale(covariates[,-c(1:2)])
pgs_names = c('EA4','externalizing','ADHD1','AFB2','EVERSMOKE2'
        ,'DEP1','height','bmi','NEBwomen2')
sib = FALSE

for (pgs_name in pgs_names){
    pgs = read.table(paste0('~/pgs/',pgs_name,'.pgs.txt'),header=T)
    if (sib){pgs_sib = read.table(paste0('~/pgs/',pgs_name,'_sib.pgs.txt'),header=T)}
    print(paste('PGS:',pgs_name))
    for (phenotype in binary_outcomes){
        ## 1-3 generation analysis without sib ##
        p = pgs
        p$phenotype = phenotypes[match(p$IID,phenotypes$IID),phenotype]
        p = p[!is.na(p$phenotype),]
        covar = as.matrix(covariates[match(p$IID,covariates$IID),-c(1:2)])
        # Standardize pgs variances
        p[,c('proband','paternal','maternal','gpp','gpm','gmp','gmm')] = 
            p[,c('proband','paternal','maternal','gpp','gpm','gmp','gmm')]/sd(p$proband)
        # 1 gen logistic linear mixed model fit
        glmm = glmer(phenotype ~ covar+proband+(1|FID),data=p,family=binomial(link='logit'),nAGQ=0,control=glmerControl(optimizer='bobyqa'))
        if (!is.null(glmm@optinfo$conv$lme4)){
            print('1 gen failed to converge')
            glmm = glm(phenotype ~ covar+proband,data=p,family=binomial(link='logit'))
        }
        write.table(summary(glmm)$coefficients,paste0('~/pgs/',pgs_name,'/',phenotype,'.1.coefficients.txt'),quote=F)
        write.table(as.matrix(vcov(glmm)),paste0('~/pgs/',pgs_name,'/',phenotype,'.1.vcov.txt'),quote=F) 
        # 2 gen logistic linear mixed model fit
        g2 = try({glmm = glmer(phenotype ~ covar+proband+paternal+maternal+(1|FID),data=p,family=binomial(link='logit'),nAGQ=0,control=glmerControl(optimizer='bobyqa'))})
        if (class(g2)=='try-error'){
            print('2 gen failed')
            glmm = glm(phenotype ~ covar+proband+paternal+maternal,data=p,family=binomial(link='logit'))
        } else if (g2@optinfo$message!='Normal exit from bobyqa' | !is.null(g2@optinfo$conv$lme4)){
            print('2 gen failed to converge')
            glmm = glm(phenotype ~ covar+proband+paternal+maternal,data=p,family=binomial(link='logit'))
        }
        write.table(summary(glmm)$coefficients,paste0('~/pgs/',pgs_name,'/',phenotype,'.2.coefficients.txt'),quote=F)
        write.table(as.matrix(vcov(glmm)),paste0('~/pgs/',pgs_name,'/',phenotype,'.2.vcov.txt'),quote=F)
        # 3 gen logistic linear mixed model fit
        g3= try({glmm = glmer(phenotype ~ covar+proband+paternal+maternal+gpp+gpm+gmp+gmm+(1|FID),data=p,family=binomial(link='logit'),nAGQ=0,control=glmerControl(optimizer='bobyqa'))})
        if (class(g3)=='try-error'){
            print('3 gen failed')
            glmm = glm(phenotype ~ covar+proband+paternal+maternal+gpp+gpm+gmp+gmm,data=p,family=binomial(link='logit'))
        } else if (g3@optinfo$message!='Normal exit from bobyqa' | !is.null(g3@optinfo$conv$lme4)){
            print('3 gen failed to converge')
            glmm = glm(phenotype ~ covar+proband+paternal+maternal,data=p,family=binomial(link='logit'))
        }
        write.table(summary(glmm)$coefficients,paste0('~/pgs/',pgs_name,'/',phenotype,'.3.coefficients.txt'),quote=F)
        write.table(as.matrix(vcov(glmm)),paste0('~/pgs/',pgs_name,'/',phenotype,'.3.vcov.txt'),quote=F)
        if (sib){
            ## 2 generation analysis with sib ##
            p = pgs_sib
            p$phenotype = phenotypes[match(p$IID,phenotypes$IID),phenotype]
            p = p[!is.na(p$phenotype),]
            covar = as.matrix(covariates[match(p$IID,covariates$IID),-c(1:2)])
            # Standardize pgs variances
            p[,c('proband','sibling','paternal','maternal')] = 
                p[,c('proband','sibling','paternal','maternal')]/sd(p$proband)
            # 2 gen logistic linear mixed model fit
            glmm = glmer(phenotype ~ covar+proband+sibling+paternal+maternal+(1|FID),data=p,family=binomial(link='logit'),nAGQ=0)
            write.table(summary(glmm)$coefficients,paste0('~/pgs/',pgs_name,'/',phenotype,'.sib.coefficients.txt'),quote=F)
            write.table(as.matrix(vcov(glmm)),paste0('~/pgs/',pgs_name,'/',phenotype,'.sib.vcov.txt'),quote=F)}
        print(phenotype)
    }
}



