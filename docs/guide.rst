=====
Guide
=====

**Introduction**

SNIPar (single nucleotide imputation of parents) is a python library for imputing missing parental genotypes from observed genotypes in a nuclear family,
and for performing family based genome-wide association and polygenic score analyses using the resulting imputed parental genotypes.

In the main SNIPar directory, there is a script for
imputing missing parental genotypes (impute_runner.py). 

The script outputs expected genotypes of missing parents, which are used as input for the fGWAS.py
script that perform family based GWAS using observed and imputed parental genotypes. 

The impute_runner.py script takes un-phased genotypes in .bed format, and phased haplotypes in .bgen format. The imputation becomes more efficient when using phased haplotypes, at the cost of slower runtime. The script requires IBD segments in the KING (https://people.virginia.edu/~wc9c/KING/manual.html)
format as input, where IBD segments shared between first-degree relatives are used. This can be computed by using *--ibdseg --degree 1* options in KING. 

The script will construct a pedigree for you if you
provide it with the KING relatedness inference (output using the --related --degree 1 options) and age & sex information. Alternatively, a user-input pedigree can be provided. Providing
the script with KING output is recommended since this ensures the pedigree is constructed in the correct way. 

The pedigree file is a plain text file
with header and columns: FID (family ID), IID (individual ID), FATHER_ID (ID of father), MOTHER_ID (ID of mother).

Note that individuals are assumed to have unique individual IDS (IID).

Siblings are identified through individuals that have the same FID.

We recommend working through the tutorial (https://github.com/AlexTISYoung/SNIPar/edit/master/docs/tutorial.rst) to get an idea of the workflow required for a full analysis.

Family based GWAS is performed using the fGWAS.py script, which uses a linear mixed model with a random-effect that models correlations between individuals with the same family ID to estimate direct and indirect effects for genome-wide SNPs. 

Polygenic score analyses are performed using the fPGS.py script, which computes polygenic scores for individuals and their first degree relatives based on observed/imputed genotypes. It also estimates direct and indirect effects of the polygenic score in the linear mixed model. 

***Package Install Instructions**

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

    'python setupy.py install'
    
    'python setup.py build_ext --inplace'

**Running tests**

To check that the code is working properly and that the C modules have compiled, you should
run tests. To run the tests, in the main SNIPar directory enter the command:

    ``python setup.py pytest``




