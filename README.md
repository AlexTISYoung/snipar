# sibreg
sibreg is a python library for performing robust GWAS using sibling pairs with a random effects
model for within family phenotypic correlations


# Main features:

sibreg class: random effects regression model allowing for intra-class correlation

pGWAS.py: script for estimating direct genetic effects, indirect genetic effects, and
confounding using sibling and imputed parental genotypes

sGWAS.py:  script for estimating 'within-family' and 'between-family' effects of SNPs
using sibling genotypes

# Documentation

Documentation for the modules and scripts is at: https://readthedocs.org/projects/sibreg

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