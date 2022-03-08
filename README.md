This is a branch of snipar was created to enable replication of the methods employed in EA4 (Okbay et al., Polygenic prediction of educational attainment within and between families from genome-wide association analyses in 3 million individuals) for computing direct and population effects of the polygenic indices (PGIs), and for analysis of assortative mating using PGIs. 

# Relevant scripts:

Preprocess phenotype data in UKB Biobank: process_phenotypes_UKB.R

Preprocess phenotype data in Generation Scotland: process_phenotypes_GS.R

Estimate direct, paternal, maternal, and population effects of PGIs using a linear mixed model: EA4_PGI_analysis.py

Meta-analyse direct, paternal, maternal PGI effect estimates from UK Biobank, Generation Scotland, and Swedish Twin Resource: meta_analysis.R

Analyse evidence for assortative mating using PGIs in UK Biobank: AM_analysis_UKB.R

Analyse evidence for assortative mating using PGIs in Generation Scotland: AM_analysis_GS.R

Meta-analyse evidence for assortative mating from UK Biobank and Generation Scotland: AM_meta.R

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
- pooch
- numba

We highly recommend using a python distribution such as Anaconda 3 (https://store.continuum.io/cshop/anaconda/).
This will come with both numpy and scipy installed and can include an MKL-compiled distribution
for optimal speed.

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
