# SNIPar

*snipar* (single nucleotide imputation of parents) is a python library for inferring identity-by-descent (IBD) segments shared between siblings, imputing missing parental genotypes from observed genotypes and IBD segments in a nuclear family, and for performing
family based genome-wide association and polygenic score analyses using observed and/or imputed parental genotypes. 

# Main features:

Infer identity-by-descent segments shared between siblings (ibd.py). 

Impute missing parental genotypes given the observed genotypes in a nuclear family (impute.py).

Perform family based GWAS using observed and imputed parental genotypes (gwas.py). 

Compute polygenic scores for probands, siblings, and parents from SNP weights using observed/imputed parental genotypes, and perform family
 based analysis of polygenic scores (pgs.py script). 
 
 Compute genome-wide correlations between different effects estimated by gwas.py (correlate.py). 

# Documentation

It is recommended to work through the tutorial: https://github.com/AlexTISYoung/SNIPar/blob/master/docs/tutorial.rst

and read the guide: https://snipar.readthedocs.io/en/latest/guide.html

Documentation for the modules and scripts is at: https://snipar.readthedocs.io/en/latest/


# Package Install Instructions

*snipar* has the following dependencies:

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
- pooch
- numba
- statsmodels
- scikit-learn

We highly recommend using a python distribution such as Anaconda 3 (https://store.continuum.io/cshop/anaconda/).

To install from source, clone the git repository, and in the directory
containing the SNIPar source code, at the shell type

    'python setup.py install'
   
If you have problems with the installation, it may be due to package conflicts with your existing Python installation. To overcome this, you can 
try installing in a virtual environment. In a bash shell, this could be done by using the following commands in the SNIPar directory:
    
    'python -m venv env'
    'source env/bin/activate'
    'python setup.py install' 

# Running tests

To check that the code is working properly and that the C modules have compiled, you should
run tests. To run the tests, in the main SNIPar directory enter the command:

    ``python setup.py pytest``
