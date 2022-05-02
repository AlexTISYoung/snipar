.. _output_files
============
Output files
============

IBD segments file 
-----------------
.. _ibd_segments_file:

imputed parental genotypes file 
-------------------------------
.. _imputed_file:

For each chromosome i, an HDF5 file is created. The resulting HDF5 files contains the following key,values:
    'imputed_par_gts'
        imputed parental genotypes. It's the imputed missing parent if only one parent is missing and the imputed average of the both parents if both are missing.

    'pos'
        the position of SNPs(in the order of appearance in genotypes)

    'families'
        family ids of the imputed parents(in the order of appearance in genotypes)

    'parental_status'
        a numpy array where each row shows the family status of the family of the corresponding row in families.

    'sib_ratio_backup'
        An array with the size of number of snps. Show the ratio of backup imputation among offspring imputations in each snp.

    'parent_ratio_backup'
        An array with the size of number of snps. Show the ratio of backup imputation among parent-offspring imputations in each snp.

    'mendelian_error_ratio'
        Ratio of mendelian errors among parent-offspring pairs for each snp

    'estimated_genotyping_error'
        estimated for each snp using mendelian_error_ratio and maf

    'ratio_ibd0'
        ratio of families with offsprings in ibd0 to all the fams.

    'bim_columns'
        Columns of the resulting bim file

    'bim_values'
        Contents of the resulting bim file

    'pedigree'
        pedigree table Its columns are has_father, has_mother, single_parent respectively.

    'non_duplicates'
        Indexes of the unique snps. Imputation is restricted to them.

    'standard_f'
        Whether the allele frequencies are just population average instead of MAFs estimated using PCs

    'MAF_*'
        info about the MAF estimator if MAF estimator is used.
    

text summary statistics
-----------------------
.. _sumstats_text:

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