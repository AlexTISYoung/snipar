.. _input files:
===========
Input files
===========

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
when phase information can be used (see Young et al. [ref]).
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

The kinship file is as output by KING: https://www.kingrelatedness.com/manual.shtml#RELATED.

agesex file
-----------
.. _agesex: 

This is a white-space delimited text file with header "FID", "IID", "sex", "age".
Each row contains the family-ID, individual-ID, age, and sex of one individual. 
Male and Female sex should be represented with 'M' and 'F' respectively.
The age column is used for distinguishing between parent and child in a parent-offsring relationship inferred from the :ref:`kinship file <kinship>`.
ID1 is a parent of ID2 if there is a parent-offspring (PO) relationship between them and 'ID1' is at least 12 years older than ID2.

phenotype file
--------------
.. _phenotype:

The phenotype file is a white-space delimited text file with header "FID", "IID", "phenotype_1", "phenotype_2", ...,
corresponding to columns for family-ID, individual-ID, and phenotype values for the different phenotypes. 
To specify the k'th phenotype for analysis (relevant for :ref:`gwas.py <gwas.py>` and :ref:`pgs.py <pgs.py>`),
add '--phen_index k' to your command; by default, the first phenotype will be used.  

weights file
------------
.. _weights: 

