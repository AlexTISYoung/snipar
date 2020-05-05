Guide
************

**Introduction**

sibreg is a python library for imputing missing parental genotypes from observed sibling and parental genotypes,
and for performing robust GWAS using the resulting imputed parental genotypes

In the sibreg/sibreg/bin subdirectory, there is a script for
imputing missing parental genotypes:

--'impute_runner.py'(:doc:`impute_runner`)
    imputes the expected sum of maternal and paternal genotypes given genotype data on the sibling
    offspring of the parents and IBD sharing between the sibling offspring

The script outputs expected genotypes of missing parents, and are used as input for the fGWAS.py
script that perform GWAS using the missing parental genotypes. See tutorial for an example of use.

All of the above scripts require either provision of a pedigree file or the script will construct a pedigree for you if you
provide it with the KING relatedness inference (output using the --related --degree 1 options) and age & sex information. Providing
the script with KING output is recommended.

The pedigree file is a plain text file
with header and columns: FID (family ID), IID (individual ID), FATHER_ID (ID of father), MOTHER_ID (ID of mother).

Note that individuals are assumed to have unique individual IDS (IID).

Siblings are identified through individuals that have the same FID.

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




