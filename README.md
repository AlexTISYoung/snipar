# SNIPar

SNIPar (single nucleotide imputation of parents) is a python library for imputing missing parental genotypes from observed genotypes in a nuclear family, and for performing
family based genome-wide association and polygenic score analyses. 

# Main features:

Impute missing parental genotypes given the observed genotypes in a nuclear family (impute_runner.py).

Perform family based GWAS using observed and imputed parental genotypes (fGWAS.py). 

Compute polygenic scores for probands, siblings, and parents from SNP weights using observed/imputed parental genotypes, and perform family
 based analysis of polygenic scores (fPGS.py script). 

# Documentation

It is recommended to work through the tutorial: https://snipar.readthedocs.io/en/latest/tutorial.html

and read the guide: https://snipar.readthedocs.io/en/latest/guide.html

Documentation for the modules and scripts is at: https://snipar.readthedocs.io/en/latest/


# Package Install Instructions

SNIPar has the following dependencies:

python 3.7

Packages:

- h5py
- bgen-reader
- numpy
- scipy
- pysnptools
- pandas
- networkx
- Cython

We highly recommend using a python distribution such as Anaconda 3 (https://store.continuum.io/cshop/anaconda/).
This will come with both numpy and scipy installed and can include an MKL-compiled distribution
for optimal speed.

To install from source, clone the git repository, and in the directory
containing the SNIPar source code, at the shell type

    'python setup.py install'

# Running tests

To check that the code is working properly and that the C modules have compiled, you should
run tests. To run the tests, in the main SNIPar directory enter the command:

    ``python setup.py pytest``
