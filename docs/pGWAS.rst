Documentation for pGWAS.py script
====================================


This script uses genotypes of siblings and parental genotypes imputed from sibling genotypes to estimate direct genetic effects, indirect genetic effects from siblings,
and indirect genetic effects from parents effects/confounding effects.

The genentic data is specified in a HDF5 file. Let M be the number of sibling pairs in the data. The HDF5 file must include datasets:

--ped
    a [M x 3] array, where each row corresponds to a sibpair. The first column is the family ID (FID),
    the second column is the IID (individual ID) of the first sib in the pair, and the third column
    is the IID of the second sib in the pair.

--gts
    an [M x 3 x L] array, where each row corresponds to a sib pair. The second index gives the
    genotype of the first sib, the genotype of the second sib, and the parental genotype imputed
    from this sib pair. The third index gives the SNP.

--freqs
    an [L] vector of allele frequencies

--vnames
    an [L] vector of names of the SNPs

The scripts fits models for all SNPs in a .hdf5 file passing MAF and missingness thresholds.

The phenotype file should be a tab separate text file with columns FID, IID, Y1, Y2, ...

Siblings should have the same family id (FID), and non-siblings should have different FIDs.

The covariate file formats is the same. The first
column is family ID, and the second column is individual ID; subsequent columns are phenotype or covariate
observations.

Minimally, the script will output a file outprefix.models.gz, which contains a table of within direct effects,
 indirect effects from siblings, and parental/confounding effects, along with their correlations

If covariates are also specified, it will output estimates of the covariate effects from the null model as
outprefix.null_mean_effects.txt. --no_covariate_estimates suppresses this output.

**Arguments**

Required positional arguments:

**genofile**
   Path to genotype file in HDF5 format

**phenofile**
   Location of the phenotype (response) file with format: FID, IID, y1, y2, ...

**outprefix**
   Location to output csv file with association statistics

Options:

--mean_covar
   Location of mean covariate file (default no mean covariates)

--fit_covariates
   Fit covariates for each locus. Default is to fit covariates for the null model and project out (mean) and rescale (variance)'

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

   ``python pGWAAS.py genotypes.bed phenotype.fam phenotype``

This will estimate between effects for all the SNPs in genotypes.bed passing MAF and missingness thresholds, using the first phenotype in phenotype.fam. It will output
the results of fitting the models to phenotype.models.gz.