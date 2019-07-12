Documentation for poGWAS.py script
====================================

This script imputes the missing parental genotype in a families where one parent and at least one offspring are observed.

**Arguments**

Required arguments:

**chr**
    which chromosome to do imputation for

**gts**
    path to .bed file with genotypes of siblings

**ibd**
    path to file with IBD information. Has columns: IID of first sib, IID of second sib, chromosome, start position of segment,
    indicator for IBD_2 state ('True' implies IBD_2; 'False' implies IBD_1); size of segment in cM (not required); end position of segment;
    number of SNPs in segment (not required).

**ped**
    path to pedigree file with columns IID (individual ID), FID (family ID), FATHER_ID (ID of father), MOTHER_ID (ID of mother).

**out**
    prefix of output hdf5 file