# sibreg
sibreg is a python library for performing robust GWAS using sibling pairs with a random effects
model for within family phenotypic correlations


# Main features:

regrnd class: given data, find the parameters that maximise
the likelihood and their sampling distribution.


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