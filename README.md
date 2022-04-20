# snipar

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

It is recommended to work through the tutorial: https://snipar.readthedocs.io/en/latest/tutorial.html

# Installing Using Pip

*snipar* currently supports Python 3.7-3.9 on Linux, Windows, and Mac OSX. We recommend using a python distribution such as Anaconda 3 (https://store.continuum.io/cshop/anaconda/). 

The easiest way to install is using pip:

    pip install snipar

Sometimes this may not work because the pip in the system is outdated. You can upgrade your pip using:

    pip install --upgrade pip

# Virtual Environment

You may encounter problems with the installation due to Python version incompatability or package conflicts with your existing Python environment. To overcome this, you can try installing in a virtual environment. In a bash shell, this could be done by using the following commands in your directory of choice:
    
    python -m venv path-to-where-you-want-the-virtual-environment-to-be

You can activate and use the environment using

    source path-to-where-you-want-the-virtual-environment-to-be/bin/activate

# Installing From Source
To install from source, clone the git repository, and in the directory
containing the *snipar* source code, at the shell type:

    pip install .
   
# Running tests
To check that the code is working properly and that the C modules have been compiled, you should run tests. To run the tests, after the installation run this command:

    python -m unittest snipar.tests
