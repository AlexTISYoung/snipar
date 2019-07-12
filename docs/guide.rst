Guide
************

**Introduction**

sibreg is a python library for performing robust GWAS using nuclear family data with a random effects
model for within family phenotypic correlations

sibreg also contains scripts for imputing missing parental genotypes from sibling pairs (provided IBD sharing information)
and from parent-offspring pairs (with one missing parent)

In the sibreg/bin subdirectory, there are two scripts for
imputing missing parental genotypes:


--'impute_from_sibs.py'(:doc:`impute_from_sibs`)
    imputes the expected sum of maternal and paternal genotypes given genotype data on the sibling
    offspring of the parents and IBD sharing between the sibling offspring

--'impute_po.py'(:doc:`impute_po`)
    imputes the expected genotype of the missing parent given one parent and a number of sibling offspring
    of the parent and the missing parent

These two scripts output expected genotypes of missing parents, and are used as input for
scripts that perform GWAS using the missing parental genotypes.

There are three scripts for performing GWAS depending on the number of missing parental genotypes:

--'pGWAS.py'(:doc:`pGWAS`)
    Performs GWAS using observed sibling genotypes and the missing parental genotypes imputed from
    the sibling genotypes (produced by 'impute_from_sibs.py'(:doc:`impute_from_sibs`))

--'poGWAS.py'(:doc:`pGWAS`)
    Performs GWAS using observed sibling genotypes, the single observed parental genotype in each family, and the imputed missing parental genotypes
    (produced by 'impute_po.py'(:doc:`impute_from_sibs`))

--'triGWAS.py'(:doc:`pGWAS`)
    Performs GWAS using observed sibling genotypes and observed maternal and paternal genotypes

All of the above scripts require provision of a pedigree file. The pedigree file is a plain text file
with header and columns: IID (individual ID), FID (family ID), FATHER_ID (ID of father), MOTHER_ID (ID of mother).

Siblings are identified through individuals that have the same FID and the same FATHER_ID and MOTHER_ID.

The core model is the sibreg model (:class:`sibreg.model`), which consists of a linear model for the mean along
with a vector of class labels that allows for correlations within-class. (The correlations within-class result
from modelling mean differences between classes as independent, normally distributed random effects.) For
the GWAS applications, the classes are the families given by the pedigree file.

The documentation for the sibreg module (:doc:`sibreg`) contains information on how to define a :class:`sibreg.model`,
how to optimise a :class:`sibreg.model`, how to predict from
a :class:`sibreg.model`, and how to simulate a :class:`sibreg.model`.

**Running tests**

To check that the code is working properly and computing likelihoods and gradients accurately, you can
run tests. In the sibreg/tests subdirectory, type

    ``python tests.py``

The output should say

    ``Ran 4 tests in...``

    ``OK``





