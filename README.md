# sibreg
sibreg is a python library for performing robust GWAS using nuclear family data with a random effects
model for within family phenotypic correlations

sibreg also contains scripts for imputing missing parental genotypes from sibling pairs (provided IBD sharing information)
and from parent-offspring pairs (with one missing parent)

# Main features:

Impute the expected sum of maternal and paternal genotypes given genotype data on the sibling
    offspring of the parents and IBD sharing between the sibling offspring

Imputes the expected genotype of the missing parent given one parent and a number of sibling offspring
    of the parent and the missing parent

Performs GWAS using observed sibling genotypes and the missing parental genotypes imputed from
    the sibling genotypes (produced by 'impute_from_sibs.py')

Performs GWAS using observed sibling genotypes, the single observed parental genotype in each family, and the imputed missing parental genotypes
    (produced by 'impute_po.py')

Performs GWAS using observed sibling genotypes and observed maternal and paternal genotypes

# Documentation

Documentation for the modules and scripts is at: https://readthedocs.org/projects/sibreg/latest

# Package Install Instructions

sibreg has the following dependencies:

python 2.7

Packages: 

- numpy
- scipy
- pysnptools

We highly recommend using a python distribution such as Anaconda (https://store.continuum.io/cshop/anaconda/). 
This will come with both numpy and scipy installed and can include an MKL-compiled distribution
for optimal speed. 

To install from source, clone the git repository, and in the directory
containing the sibreg source code, at the shell type

    'sudo python setupy.py install'

or, on the windows command prompt, type

    'python setup.py install' 
    
# Running tests

The tests directory contains scripts for testing the computation of 
the likelihoods, gradients, and maximum likelihood solutions
To run these tests, a further dependency is required: numdifftools. 

To run the tests, first install sibreg. Change to the tests/ directory and at the shell type

    'python test.py'

Tests should run without any failures.