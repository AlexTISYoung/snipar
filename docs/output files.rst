.. _output_files
============
Output files
============

IBD segments file 
-----------------
.. _ibd_segments_file:

The :ref:`ibd.py <ibd.py>` script outputs one gzipped text file per chromosome containing the IBD segments for each sibling pair
in the input. The IBD segments file is a tab delimited text file where each row contains information on a particular IBD segment.
It has columns: 

    'ID1'
        Individual ID (IID) of the first sibling in pair
    'ID2'
        Inidividual ID (IID) of the second sibling in the pair
    'IBDType'
        The IBD type (0, 1, or 2) of the segment
    'Chr'
        The chromosome of the segment
    'start_coordinate'
        The base-pair (bp) position of the start of the segment (inclusive)
    'stop_coordinate'
        The base-pair (bp) position of the end of the segment (inclusive)
    'start_SNP'
        The ID of the variant at the start of the segment
    'stop_SNP'
        The ID of the variant at the end of the segment
    'length'
        The length of the segment in centi Morgans (cM)
    


imputed parental genotypes file 
-------------------------------
.. _imputed_file:

For each chromosome, the :ref:`impute.py <impute.py>` script outputs a HDF5 file containing the imputed parental genotypes. 
(Users typically do not need to look in this file, but if they want to, HDF5 files can be read in *R* using *rhdf5* and in Python
using *h5py*.) Consider imputing missing parental genotypes for n families on a chromosome with L genotyped variants. 
The resulting HDF5 file will contain the following datasets:

    'imputed_par_gts'
        [n x L] floating point array of imputed parental genotypes. It's the imputed missing parent if only one parent is missing and the imputed average of the both parents if both are missing.

    'pos'
        [L] vector of base pair (bp) position of variants (in the order of appearance in genotypes)

    'families'
        [n] vector of family (sibship) ids of the imputed parents (in the order of appearance in genotypes); these are used internally by *snipar* and are distinct from FIDs in input files

    'parental_status'
        [n] vector where each row shows the family status [ref] of the family of the corresponding row in families.

    'sib_ratio_backup'
        [L] vector giving the ratio of backup imputation (not using IBD information) among families with 2 or more genotyped siblings for each variant.

    'parent_ratio_backup'
        [l] vector giving the ratio of backup imputation among parent-offspring imputations for each variant.

    'mendelian_error_ratio'
        [L] vector giving ratio of mendelian errors among parent-offspring pairs for each variant

    'estimated_genotyping_error'
        [L] estimated genotyping error rate for each variant 

    'ratio_ibd0'
        [L] vector giving the fraction of sibships with an observed IBD0 pair

    'bim_columns'
        [k] vector giving the column names of the 'bim' file (table containing information similar to a PLINK .bim file)

    'bim_values'
        [L x k] matrix of variant-level information (see 'bim_columns')

    'pedigree'
        pedigree with columns family ID, individual ID, father ID, mother ID, has_father, has_mother. 
        the family ID here corresponds to the family ID in the 'families' dataset, which 
        indexes the rows of 'imputed_par_gts'. 'has_father' and 'has_mother' denotes whether a genotyped father
        or genotyped mother (respectively) was used in the imputation

    'non_duplicates'
        vector of indexes of the unique snps; imputation is restricted to these SNPs

    'standard_f'
        Whether the allele frequencies are just population average instead of MAFs estimated using PCs

    'MAF_*'
        info about the MAF estimator if MAF estimator is used.
    

text summary statistics
-----------------------
.. _sumstats_text:

The :ref:`gwas.py <gwas.py>` script outputs one gzipped text file per chromosome containing the summary statistics for variants in the input. 
Variants that have been filtered out (by having an MAF below the threshold, too much missingness, or missing IBD/genetic map information)
will appear in the output file but the summary statistics will be 'nan'. The sumstats file is a white-space delimited text file. Exactly which columns are present depends on the model used.
If using a model with proband and maternal and paternal genotypes, the sumstats file will have the following columns:

    'chromosome'
        The chromosome of the variant
    'SNP'
        The ID of the variant
    'pos'
        The base-pair (bp) position of the variant
    'A1'
        The effect allele
    'A2'
        The alternative allele
    'freq'
        The frequency of the 'A1' effect allele
    'direct_N'
        The effective sample size for estimation of the direct effect
    'direct_Beta'
        The estimated direct effect
    'direct_SE'
        The standard error of the direct effect estimate
    'direct_Z'
        The Z-score of the direct effect estimate
    'direct_log10_P'
        The negative log10 P-value for a non-zero direct effect
    'paternal_N'
        The effective sample size for estimation of the paternal non-transmitted coefficient (NTC)
    'paternal_Beta'
        The estimated paternal NTC
    'paternal_SE'
        The standard error of the paternal NTC estimate
    'paternal_Z'
        The Z-score of the paternal NTC estimate
    'paternal_log10_P'
        The negative log10 P-value for a non-zero paternal NTC
    'maternal_N'
        The effective sample size for estimation of the maternal non-transmitted coefficient (NTC)
    'maternal_Beta'
        The estimated maternal NTC
    'maternal_SE'
        The standard error of the maternal NTC estimate
    'maternal_Z'
        The Z-score of the maternal NTC estimate
    'maternal_log10_P'
        The negative log10 P-value for a non-zero maternal NTC
    'avg_NTC_N'
        The effective sample size for estimation of the average non-transmitted coefficient (NTC): average of maternal and paternal NTCs
    'avg_NTC_Beta'
        The estimated average NTC
    'avg_NTC_SE'
        The standard error of the average NTC estimate
    'avg_NTC_Z'
        The Z-score of the average NTC estimate
    'avg_NTC_log10_P'
        The negative log10 P-value for a non-zero average NTC
    'population_N'
        The effective sample size for estimation of the population effect: sum of direct effect and average NTC
    'population_Beta'
        The estimated population effect
    'population_SE'
        The standard error of the population effect estimate
    'population_Z'
        The Z-score of the population effect estimate
    'population_log10_P'
        The negative log10 P-value for a non-zero population effect
    'r_direct_avg_NTC'
        The sampling correlation between the direct effect and average NTC estimates
    'r_direct_population'
        The sampling correlation between the direct effect and population effect estimates
    'r_paternal_maternal'
        The sampling correlation between paternal and maternal NTC estimates

Note that, if using parental genotypes imputed from siblings (without any observed parents),
then separate maternal and paternal NTCs cannot be estimated, so only the average NTC will appear 
in the summary statistics output. Also, if '--fit_sib' is used to include an indirect effect from siblings,
this will be included in the output. 

HDF5 summary statistics 
-----------------------
.. _sumstats_hdf5:

pgs file
--------
.. _pgs_file: 

pgs effects
-----------
.. _pgs_effects:

pgs effects sampling covariance
-------------------------------
.. _pgs_vcov: