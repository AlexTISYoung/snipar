Documentation for impute_po.py script
====================================

This script imputes the missing parental genotype in a families where one parent and at least one offspring are observed.

**Arguments**

Required arguments:

**gts**
    path to .bed file with genotypes of siblings observed parents

**ped**
    path to pedigree file with columns IID (individual ID), FID (family ID), FATHER_ID (ID of father), MOTHER_ID (ID of mother).

**out**
    prefix of output hdf5 file