=====
Guide
=====

Introduction
------------

*snipar* (single nucleotide imputation of parents) is a Python package for inferring identity-by-descent (IBD) segments shared between siblings, imputing missing parental genotypes, and for performing
family based genome-wide association and polygenic score analyses using observed and/or imputed parental genotypes.

*snipar* can use any genotyped samples who have at least one genotyped full-sibling or parent

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
   :scale: 35 %
   :alt: typical snipar workflow

   Illustration of a typical workflow for performing family-based GWAS

A *snipar* workflow requires input files in certain formats. See :ref:`input files <input files>`.

The :ref:`tutorial <tutorial>` allows you to work through an example workflow before trying real data. 

Inferring identity-by-descent segments 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If your sample contains full-sibling pairs (without both parents genotyped),
it is necessary to first infer the identity-by-descent (IBD) segments
shared between the siblings before imputing the missing parental genotypes. 
If your sample does not contain any full-sibling pairs, but has genotyped
parent-offspring pairs (i.e. one parent's genotype is missing), imputation
can proceed without inferring IBD. 

*snipar* contains a Hidden Markov Model (HMM) algorithm for inferring IBD shared between siblings, 
which can be accessed through the command line script ibd.py. 

The ibd.py script requires the :ref:`observed genotypes <observed genotypes>` of the siblings and information
on the sibling and parent-offspring relations in the genotyped sample. 

To infer IBD, one can use a smaller set of genetic variants than one intends to 
use in downstream analyses (imputation, gwas, etc.). 
For example, one could use the variants on a genotyping array to
infer IBD segments, and these IBD segments could be used to impute missing parental genotypes
for all variants imputed from a reference panel. This can be useful since the accuracy of IBD
inference plateaus as the density of variants increases, so inputting millions of variants
imputed from a reference panel to ibd.py will result in a long computation time for little gain
in accuracy over using variants from a genotyping array. 

The information on the relations present in genotyped sample can be provided through a pedigree file[ref] or through
the output of KING relationship inference (as output using the --related --degree 1 options: see https://www.kingrelatedness.com/manual.shtml#RELATED)
along with a file giving the age and sex information[ref] on the genotyped sample.
(The age and sex information along with the parent-offspring and sibling relations inferred by KING are used to construct a pedigree
if a pedigree is not provided.)

The algorithm requires a genetic map to compute the probabilities of transitioning between different IBD states. 
If the genetic map positions (in cM) are provided in the .bim file (if using .bed formatted genotypes), the script will use these. 
Alternatively, the *--map* argument allows the user to specify a genetic map in the same format as used by SHAPEIT 
(https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#formats).
If no genetic map is provided, then the deCODE sex-averaged map on Hg19 coordinates (Halldorsson, Bjarni V., et al. "Characterizing mutagenic effects of recombination through a sequence-level genetic map." Science 363.6425 (2019).),
which is distributed as part of *snipar*, will be used. 

The HMM employs a genotyping error model that requires a genotyping error probability parameter. 
By default, the algorithm will estimate the per-SNP genotyping error probability from Mendelian errors
observed in parent-offspring pairs. However, if your data does not contain any genotyped parent-offspring pairs, 
then you will need to supply a genotyping error probability to ibd.py.
If you have no external information on the genotyping error rate in your data, using a value of 1e-4 has 
worked well when applied to typical genotyping array data. 

The HMM will output the IBD segments to a gzipped text file with suffix ibd.segments.gz. As part of the algorithm,
LD scores are calculated for each SNP. These can also be output in LDSC format using the --ld_out option. 


Imputing missing parental genotypes 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Family-based genome-wide association analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*snipar* includes a command-line script 

Family-based polygenic score analyses
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Estimating correlations between effects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~