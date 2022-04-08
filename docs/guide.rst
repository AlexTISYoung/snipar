=====
Guide
=====

**Introduction**

*snipar* (single nucleotide imputation of parents) is a python library for imputing missing parental genotypes from observed genotypes in a nuclear family,
and for performing family based genome-wide association and polygenic score analyses using the resulting imputed parental genotypes.

*snipar* contains command line scripts:
    #. for inferring IBD segments shared between siblings (ibd.py)
    #. imputing missing parental genotypes from observed parent/offspring genotypes and IBD segments (impute.py),
    #. performing genome-wide estimation of direct genetic effects, non-transmitted coefficients, and population effects of SNPs (gwas.py),
    #. for estimating direct effects and non-transmitted coefficients of polygenic scores (pgs.py),
    #. and for estimating genome-wide correlations between direct and population effects and direct effects and non-transmtitted coefficients (correlate.py)

