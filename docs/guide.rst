Guide
************

**Introduction**

sibreg is a python library for imputing missing parental genotypes from observed sibling and parental genotypes,
and for performing robust GWAS using the resulting imputed parental genotypes

In the sibreg/sibreg/bin subdirectory, there are two scripts for
imputing missing parental genotypes:


--'impute_runner.py'(:doc:`impute_runner`)
    imputes the expected sum of maternal and paternal genotypes given genotype data on the sibling
    offspring of the parents and IBD sharing between the sibling offspring

--'impute_po.py'(:doc:`sibreg.bin.impute_runner`)
    imputes the expected genotype of the missing parent given one parent and a number of sibling offspring
    of the parent and the missing parent

These two scripts output expected genotypes of missing parents, and are used as input for
scripts that perform GWAS using the missing parental genotypes.

There are three scripts for performing GWAS depending on the number of missing parental genotypes:

--'pGWAS.py'(:doc:`pGWAS`)
    Performs GWAS using observed sibling genotypes and the missing parental genotypes imputed from
    the sibling genotypes (produced by 'impute_from_sibs.py'(:doc:`impute_from_sibs`))

--'poGWAS.py'(:doc:`poGWAS`)
    Performs GWAS using observed sibling genotypes, the single observed parental genotype in each family, and the imputed missing parental genotypes
    (produced by 'impute_po.py'(:doc:`impute_po`))

--'triGWAS.py'(:doc:`triGWAS`)
    Performs GWAS using observed sibling genotypes and observed maternal and paternal genotypes

All of the above scripts require either provision of a pedigree file or the script will construct a pedigree for you if you
provide it with the KING relatedness inference (output using the --related --degree 1 options) and age & sex information.

The pedigree file is a plain text file
with header and columns: FID (family ID), IID (individual ID), FATHER_ID (ID of father), MOTHER_ID (ID of mother).

Note that individuals are assumed to have unique individual IDS (IID).

Siblings are identified through individuals that have the same FID and the same FATHER_ID and MOTHER_ID.

We recommend working through the (:doc:`tutorial`) to get an idea of the workflow required for a full analysis.

The core model is the sibreg model (:class:`sibreg.model`), which consists of a linear model for the mean along
with a vector of class labels that allows for correlations within-class. (The correlations within-class result
from modelling mean differences between classes as independent, normally distributed random effects.) For
the GWAS applications, the classes are the families given by the pedigree file.

The documentation for the sibreg module (:doc:`sibreg`) contains information on how to define a :class:`sibreg.model`,
how to optimise a :class:`sibreg.model`, how to predict from
a :class:`sibreg.model`, and how to simulate a :class:`sibreg.model`. Note that to run the imputation and GWAS scripts it
is not necessary for the user to interact directly with the sibreg module.

***Package Install Instructions**

sibreg has the following dependencies:

python 2.7

Packages:

- numpy
- scipy
- pysnptools
- pandas
- networkx
- Cython

We highly recommend using a python distribution such as Anaconda (https://store.continuum.io/cshop/anaconda/).
This will come with both numpy and scipy installed and can include an MKL-compiled distribution
for optimal speed.

To install from source, clone the git repository, and in the directory
containing the sibreg source code, at the shell type

    'python setupy.py install'

**Running tests**

To check that the code is working properly, you should
run tests. To run the tests, in the main sibreg directory enter the command:

    ``python setup.py tests``




