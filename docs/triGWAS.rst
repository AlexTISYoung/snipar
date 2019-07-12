Documentation for triGWAS.py script
====================================

This script uses observed genotypes of individuals and their mothers and fathers to estimate direct and indirect genetic effects/confounding effects.

The script fits models for all SNPs in common between the sibling and parental genotypes passing MAF and missingness thresholds.

The phenotype file should be a tab separate text file with columns FID, IID, Y1, Y2, ...

Siblings should have the same family id (FID), FATHER_ID, and MOTHER_ID in the pedigree file.

The script outputs a .hdf5 file with variance parameter estimates and estimates of the X^T X ('xtx') and X^T Y ('xty') matrices for each SNP as a multidimensional array,
the first index of which gives the SNP, whose corresponding ID is in the 'sid' vector. The 'N_L' vector contains the sample size for each SNP (after removal of missing genotypes).
The rows/columns of the X^T X and the columns of the X^T Y matrices correspond to: intercept, direct effect,
paternal effect, maternal effect. The .hdf5 file also contains the estimated variance
parameters: sigma2, the residual variance estimate; and tau, the estimated ratio between the within-family variance and the residual variance.

The covariate file formats is the same. The first
column is family ID, and the second column is individual ID; subsequent columns are phenotype or covariate
observations.

If covariates are also specified, it will output estimates of the covariate effects from the null model as
outprefix.null_mean_effects.txt. --no_covariate_estimates suppresses this output.

**Arguments**

Required arguments:

**gts**
    path to .bed file with genotypes of siblings observed parents

**ped**
    path to pedigree file with columns FID (family ID), IID (individual ID), FATHER_ID (ID of father), MOTHER_ID (ID of mother).

**phenofile**
    path to phenotype file

**outprefix**
    prefix of output hdf5 file

Options:

--mean_covar
   Location of mean covariate file (default no mean covariates)

--fit_covariates
   Fit covariates for each locus. Default is to fit covariates for the null model and project out the covariates'

--tau_init
   Initial value for the ratio of within family variance to residual variance. Default 1.0.

--phen_index
   If the phenotype file contains multiple phenotypes, specify the phenotype to analyse. Default is first phenotype in file.
   Index counts starting from 1, the first phenotype in the phenotye file.

--min_maf
   Ignore SNPs with minor allele frequency below min_maf (default 5%)

--missing_char
   Missing value string in phenotype file (default NA)

--max_missing
   Ignore SNPs with greater % missing calls than max_missing (default 5%)

--append
   Append results to existing output file with given outprefix (default to open new file and overwrite existing file with same name)

--no_covariate_estimates
   Suppress output of covariate effect estimates

--no_sib
    If this flag is given, no indirect genetic effects from siblings will be fit. (Fit by default.)

--fit_VC
    If this flag is given, variance components (within-family and residual) will be estimated for each SNP.
    Default is to fit for the null model and fix.