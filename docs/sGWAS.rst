Documentation for sGWAS.py script
====================================


This script uses genotypes of siblings to estimate 'within family' and 'between family' effects of SNPs.

The script fits models for all SNPs in a .bed file passing MAF and missingness thresholds.

The phenotype file should be a tab separate text file with columns FID, IID, Y1, Y2, ...

Siblings should have the same family id (FID), and non-siblings should have different FIDs.

The covariate file formats is the same. The first
column is family ID, and the second column is individual ID; subsequent columns are phenotype or covariate
observations.

Minimally, the script will output a file outprefix.models.gz, which contains a table of within family and between family
effect estimates along with their standard errors and correlation.

If covariates are also specified, it will output estimates of the covariate effects from the null model as
outprefix.null_mean_effects.txt. --no_covariate_estimates suppresses this output.

**Arguments**

Required positional arguments:

**genofile**
   Path to genotypes in BED format

**phenofile**
   Location of the phenotype (response) file with format: FID, IID, y1, y2, ...

**outprefix**
   Location to output csv file with association statistics

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


**Example Usage**

Minimal usage:

   ``python sGWAS.py genotypes.bed phenotype.fam phenotype``

This will estimate between family and within family effects for all the SNPs in genotypes.bed passing MAF and missingness thresholds, using the first phenotype in phenotype.fam. It will output
the results of fitting the models to phenotype.models.gz.