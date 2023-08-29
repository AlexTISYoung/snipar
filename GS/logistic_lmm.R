require(lme4)
binary_outcomes = c('ever.smoked')
gpardir='/disk/genetics/sibling_consortium/GS20k/alextisyoung/grandpar'
phenotypes = read.table(paste0(gpardir,'/processed_traits_noadj.txt'),header=T)
covariates = read.table(paste0(gpardir,'/covariates.fam'),header=T)
covariates = covariates[,c('FID','IID','age','sex','pcs.V3','pcs.V4')]
covariates[,-c(1:2)] = scale(covariates[,-c(1:2)])
pgs_names = c('EA4_hm3')

for (pgs_name in pgs_names){
    pgs = read.table(paste0(gpardir,'/pgs/',pgs_name,'.pgs.txt'),header=T)
    pgs_sib = read.table(paste0(gpardir,'/pgs/',pgs_name,'_sib.pgs.txt'),header=T)
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
        glmm = glmer(phenotype ~ covar+proband+(1|FID),data=p,family=binomial(link='logit'),nAGQ=2)
        write.table(summary(glmm)$coefficients,paste0(gpardir,'/pgs/',pgs_name,'/',phenotype,'.1.coefficients.txt'),quote=F)
        write.table(as.matrix(vcov(glmm)),paste0(gpardir,'/pgs/',pgs_name,'/',phenotype,'.1.vcov.txt'),quote=F) 
        # 2 gen logistic linear mixed model fit
        g2 = try({glmm = glmer(phenotype ~ covar+proband+paternal+maternal+(1|FID),data=p,family=binomial(link='logit'),nAGQ=2)})
        if (class(g2)=='try-error'){
            print('2 gen failed')
            glmm = glm(phenotype ~ covar+proband+paternal+maternal,data=p,family=binomial(link='logit'))
        }
        write.table(summary(glmm)$coefficients,paste0(gpardir,'/pgs/',pgs_name,'/',phenotype,'.2.coefficients.txt'),quote=F)
        write.table(as.matrix(vcov(glmm)),paste0(gpardir,'/pgs/',pgs_name,'/',phenotype,'.2.vcov.txt'),quote=F)
        # 3 gen logistic linear mixed model fit
        g3= try({glmm = glmer(phenotype ~ covar+proband+paternal+maternal+gpp+gpm+gmp+gmm+(1|FID),data=p,family=binomial(link='logit'),nAGQ=2)})
        if (class(g3)=='try-error'){
            print('3 gen failed')
            glmm = glm(phenotype ~ covar+proband+paternal+maternal+gpp+gpm+gmp+gmm,data=p,family=binomial(link='logit'))
        }
        write.table(summary(glmm)$coefficients,paste0(gpardir,'/pgs/',pgs_name,'/',phenotype,'.3.coefficients.txt'),quote=F)
        write.table(as.matrix(vcov(glmm)),paste0(gpardir,'/pgs/',pgs_name,'/',phenotype,'.3.vcov.txt'),quote=F)
        ## 2 generation analysis with sib ##
        p = pgs_sib
        p$phenotype = phenotypes[match(p$IID,phenotypes$IID),phenotype]
        p = p[!is.na(p$phenotype),]
        covar = as.matrix(covariates[match(p$IID,covariates$IID),-c(1:2)])
        # Standardize pgs variances
        p[,c('proband','sibling','paternal','maternal')] = 
            p[,c('proband','sibling','paternal','maternal')]/sd(p$proband)
        # 2 gen logistic linear mixed model fit
        glmm = glmer(phenotype ~ covar+proband+sibling+paternal+maternal+(1|FID),data=p,family=binomial(link='logit'),nAGQ=2)
        write.table(summary(glmm)$coefficients,paste0(gpardir,'/pgs/',pgs_name,'/',phenotype,'.sib.coefficients.txt'),quote=F)
        write.table(as.matrix(vcov(glmm)),paste0(gpardir,'/pgs/',pgs_name,'/',phenotype,'.sib.vcov.txt'),quote=F)
        print(phenotype)
    }
}
