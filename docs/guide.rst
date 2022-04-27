=====
Guide
=====

Introduction
------------

*snipar* (single nucleotide imputation of parents) is a Python package for inferring identity-by-descent (IBD) segments shared between siblings, imputing missing parental genotypes, and for performing
family based genome-wide association and polygenic score analyses using observed and/or imputed parental genotypes. 

Installation
------------

*snipar* currently supports Python 3.7-3.9 on Linux, Windows, and Mac OSX. We recommend using a python distribution such as Anaconda 3 (https://store.continuum.io/cshop/anaconda/). 


Installing Using pip
~~~~~~~~~~~~~~~~~~~~

The easiest way to install is using pip:

    pip install snipar

Sometimes this may not work because the pip in the system is outdated. You can upgrade your pip using:

    pip install --upgrade pip

You may encounter problems with the installation due to Python version incompatability or package conflicts with your existing Python environment. 
To overcome this, you can try installing in a virtual environment. 
In a bash shell, this could be done by using the following commands in your directory of choice:
    
    python -m venv path-to-where-you-want-the-virtual-environment-to-be

You can activate and use the environment using

    source path-to-where-you-want-the-virtual-environment-to-be/bin/activate

*snipar* can also be installed within a conda environment using pip. 

Installing From Source
~~~~~~~~~~~~~~~~~~~~~~~

To install from source, clone the git repository (https://github.com/AlexTISYoung/snipar), and in the directory
containing the *snipar* source code, at the shell type:

    pip install .
   
Running tests
~~~~~~~~~~~~~
To check that the code is working properly and that the C modules have been compiled, you can run the tests using this command:

    python -m unittest snipar.tests

Workflow
--------

A typical *snipar* workflow for performing family-based GWAS (see flowchart below) is:

1. Inferring identity-by-descent (IBD) segments shared between siblings (ibd.py)
2. Imputing missing parental genotypes (impute.py)
3. Estimating direct effects and non-transmitted coefficients (NTCs) of genome-wide SNPs (gwas.py)

.. figure:: snipar_flowchart.png
   :scale: 30 %
   :alt: typical snipar workflow

   This illustrates the a typical workflow for performing family-based GWAS

A *snipar* workflow requires input files in certain formats. See :ref:`input files <input files>`.

The :ref:`tutorial <tutorial>` allows you to work through an example workflow before trying real data. 

Inferring identity-by-descent segments 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If your sample contains full-sibling pairs without both parents genotyped,
it is necessary to first infer the identity-by-descent (IBD) segments
shared between the siblings before imputing the missing parental genotypes. 

If your sample does not contain any full-sibling pairs, but has genotyped
parent-offspring pairs (i.e. one parent's genotype is missing), imputation
can proceed without inferring IBD. 

*snipar* contains a Hidden Markov Model (HMM) algorithm for doing this, 
which can be accessed through the command line script ibd.py. 

The ibd.py script requires the :ref:`observed genotypes <observed genotypes>` of the siblings and information on which pairs
in the sample are sibling pairs: this information can be provided through a pedigree file or through
a the output of KING relationship inference (as output using the --related --degree 1 options: see https://www.kingrelatedness.com/manual.shtml#RELATED)
along with a file 

Imputing missing parental genotypes 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Family-based genome-wide association analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Family-based polygenic score analyses
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Estimating correlations between effects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~