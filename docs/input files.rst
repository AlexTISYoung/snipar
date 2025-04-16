.. _input files:
===========
Input files
===========

We describe the input files here. Examples are available in the :ref:`tutorial <tutorial>` data. 

IDs
---

We have followed a typical convention that input files giving individual level information (phenotype, pedigree, etc.) 
tend to have both family ID (FID) and individual ID (IID) columns. Note, however, 
that *snipar* ignores the entries in the FID column, using only individual IDs. Therefore,
entries in the IID column need to be unique. 

Observed genotypes
------------------
.. _observed genotypes:

Observed genotypes can be provided either in PLINK .bed format or in phased .bgen format. 
(Unphased .bgen format is not currently supported).
Phased .bgen files are the recommended input since the imputation is more accurate 
when phase information can be used (see Young et al. 2022 [ref]).
If one has phased genotypes in another format (for example VCF), then these can be converted
to phased .bgen format using QCTOOL (https://www.well.ox.ac.uk/~gav/qctool_v2/documentation/examples/converting.html).

The *snipar* :ref:`workflow <workflow>` has only been tested with high-quality genotype information on bi-allelic variants. 
We recommend first filtering your observed genotypes to keep only bi-allelic variants with INFO/R-square>0.99 
(if using genotypes imputed from a reference panel) and Hardy-Weinberg Equilibrium P-value>1e-6. 
Using low-quality SNPs may result in bias in the imputed parental genotypes and direct effect estimates and is not recommended. 

Observed genotypes should be split into separate files for each chromosome in the autosome (*snipar* does not support sex-chromosomes).
For example, if you have phased genotypes in .bgen format in a directory /path/to/haplotypes/ so that the genotypes for each chromosome
are in /path/to/haplotypes/chr_1.bgen, /path/to/haplotypes/chr_2.bgen, ..., /path/to/haplotypes/chr_22.bgen, then this can be specified to
*snipar* scripts with '--bgen /path/to/haplotypes/chr_@', where @ is interpreted as a numerical wildcard character. 

Pedigree 
--------
.. _pedigree:

This is a white-space delimited text file with header "FID", "IID", "FATHER_ID", "MOTHER_ID", 
corresponding to columns for family-ID (FID), individual ID (IID), father's ID, and mother's ID. 
A '0' in 'father's ID' or 'mother's ID' implies the parent is unknown. 
If the missing value is indicated by something other than '0', be sure to specify it with --pedigree_nan option.

WARNING: monozygotic (identical) twins, when coded in the pedigree the same way as full-siblings, will cause errors
in IBD inference and imputation. We recommend filtering out one individual from each identical twin pair to
prevent downstream issues. 

Instead of providing a pedigree file, one can provide the results of KING relationship inference 
and age and sex information (see below). 

kinship file
------------
.. _kinship: 

The kinship file is as output by KING: https://www.kingrelatedness.com/manual.shtml#RELATED. WARNING: KING relationship inference behaves differently if the .fam file of the input .bed file has meaningful family IDs (FIDs) and parental IDs; i.e., if your .fam file already contains pedigree/family information, KING will produce a .kin file. We have found the output of KING in this case (the .kin file) to be unpredictable in ways that causes issues with *snipar* analyses. We therefore recommend setting the family IDs (FIDs) to the individual IDs (IIDs) in .fam file input to KING, and removing information on parents (if present), before running KING with the --related command. This will output a .kin0 file containing the sibling and parent-offspring relations needed by *snipar*. 

agesex file
-----------
.. _agesex: 

This is a white-space delimited text file with header "FID", "IID", "sex", "age".
Each row contains the family-ID, individual-ID, age, and sex of one individual. 
Male and Female sex should be represented with 'M' and 'F' respectively.
The age column is used for distinguishing between parent and child in a parent-offspring relationship inferred from the :ref:`kinship file <kinship>`.
ID1 is a parent of ID2 if there is a parent-offspring (PO) relationship between them and 'ID1' is at least 12 years older than ID2.

phenotype file
--------------
.. _phenotype:

The phenotype file is a white-space delimited text file with a header. It has columns (in order) for
family-ID, individual-ID, and phenotype values for the different phenotypes. 
To specify the k'th phenotype for analysis (relevant for :ref:`gwas.py <gwas.py>` and :ref:`pgs.py <pgs.py>`),
add '--phen_index k' to your command; by default, the first phenotype will be used.  

covariate file
--------------
.. _phenotype:

The covariate file has the same format as the phenotype file (above). It is a white-space delimited text file with a header. It has columns (in order) for
family-ID, individual-ID, and covariate values for the different covariates. (Note, covariates must be numerical). 

genetic relationship (GRM) file
-------------------------------
.. _grm:

The genetic relationship file is used for modeling correlations between relatives in the :ref:`gwas.py <gwas.py>` and :ref:`pgs.py <gwas.py>` script. 
The GRM file can be specified using the --grm option. This file should be a white space or tab delimited text file where 
each row specifies a pairwise relationship. Minimally, it should have columns ID1 and ID2, giving the IDs of the pair, and 
a column named 'PropIBD' or 'relatedness' giving the relatedness coefficient between the pair. 
This is designed to work with the output of the `KING <https://www.kingrelatedness.com/manual.shtml>`_ IBD segment inference with the --ibdseg argument.
But other ways of computing relatedness coefficient can be used, such as from a pedigree. 
Alternatively, one can be specify a `GCTA GRM <https://yanglab.westlake.edu.cn/software/gcta/#MakingaGRM>`_ in the .gz format. 
This can be very slow for large samples since it requires loading all pairs but only pairs related above the 
set threshold (default 0.05) are used to compute a sparse GRM. 

weights file
------------
.. _weights: 

This file is used to input the SNP weights to the :ref:`pgs.py <pgs.py>` script for computation of the PGS. 
The weights file is a plain-text file with columns giving (minimally) the SNP ID, the SNP weight, the 
effect allele, and the alternative allele. The script is setup to process weights files as output by LD-pred
by default. If your weights file has different column names, these can be specified through the command 
line arguments of the :ref:`pgs.py <pgs.py>` script:
    '--SNP'
        the column name for the column containing the SNP IDs
    '--beta_col'
        the column name for the column with the SNP weights
    '--A1' 
        the column name for the column with the effect allele
    '--A2'
        the column name for the column with the alternative allele