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
.. _workflow:

A typical *snipar* workflow for performing family-based GWAS (see flowchart below) is:

1. Inferring identity-by-descent (IBD) segments shared between siblings (:ref:`ibd.py <ibd.py>`)
2. Imputing missing parental genotypes (:ref:`impute.py <impute.py>`)
3. Estimating direct effects and non-transmitted coefficients (NTCs) of genome-wide SNPs (:ref:`gwas.py <gwas.py>`)

.. figure:: snipar_flowchart.png
   :scale: 35 %
   :alt: typical snipar workflow

   Illustration of a typical workflow for performing family-based GWAS

A *snipar* workflow requires input files in certain formats. See :ref:`input files <input files>`.

The :ref:`tutorial <tutorial>` allows you to work through an example workflow before trying real data. 

Inputting multiple chromosomes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We recommend splitting up :ref:`observed genotype files <observed genotypes>`  by chromosome since certain
scripts in *snipar* cannot handle observed genotype files with SNPs from multiple chromosomes. 

To run scripts for all chromosomes simultaneously (recommended), the @ character can be used as a numerical wildcard.
For example, if you had observed genotype files chr_1.bed, chr_2.bed, ..., chr_22.bed, then you could specify
these as inputs to the command line scripts as "--bed chr_@". If you only want to analyse a subset of the chromosomes,
you can use the "--chr_range" argument; for example, '--bed chr_@ --chr_range 1-9' would specify analysing observed genotype
files chr_1.bed, chr_2.bed, ..., chr_9.bed. 

This will result in :ref:`output files <output_files>` that are also split by chromosome. The names of the output files
can also be specified using the numerical wildcard character, @, e.g. '--out /path/to/output/dir/chr_@'.

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

The information on the relations present in the genotyped sample can be provided through a :ref:`pedigree file <pedigree>` or through
the output of KING relationship inference (as output using the --related --degree 1 options: see https://www.kingrelatedness.com/manual.shtml#RELATED)
along with a :ref:`file giving the age and sex information <agesex>` on the genotyped sample.
(The age and sex information along with the parent-offspring and sibling relations inferred by KING are used to construct a pedigree
if a pedigree is not provided.)

The algorithm requires a genetic map to compute the probabilities of transitioning between different IBD states. 
If the genetic map positions (in cM) are provided in the .bim file (if using .bed formatted genotypes), the script will use these. 
Alternatively, the *--map* argument allows the user to specify a genetic map in the same format as used by SHAPEIT 
(https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#formats).
If no genetic map is provided, then the deCODE sex-averaged map on GRCh38 coordinates (Halldorsson, Bjarni V., et al. "Characterizing mutagenic effects of recombination through a sequence-level genetic map." Science 363.6425 (2019).),
which is distributed as part of *snipar*, will be used. 

The HMM employs a genotyping error model that requires a genotyping error probability parameter. 
By default, the algorithm will estimate the per-SNP genotyping error probability from Mendelian errors
observed in parent-offspring pairs. However, if your data does not contain any genotyped parent-offspring pairs, 
then you will need to supply a genotyping error probability.
If you have no external information on the genotyping error rate in your data, using a value of 1e-4 has 
worked well when applied to typical genotyping array data. 

The HMM will output the IBD segments to a gzipped text file with suffix ibd.segments.gz. As part of the algorithm,
LD scores are calculated for each SNP. These can also be output in LDSC format using the --ld_out option. 

Imputing missing parental genotypes 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:ref:`impute.py <impute.py>` is responsible for imputing the missing parental genotypes.
This is done for families with at least one offspring and one parent missing.
You should provide the script with information identity-by-descent (IBD) segments shared between
the siblings if there are sibling pairs present within the data. This data can be either in snipar or king format.
Use –ibd_is_king to specify which format is used.

The script needs information about family structure of the sample. You can either supply it with a :ref:`pedigree file <pedigree>` or
let it build the pedigree from :ref:`kinship <kinship>` and :ref:`agesex <agesex>` files.

If you want to do the imputation only on a subset of SNPS you can achieve  it by using -start and -end options.
This can help you with possible memory issues. -chunks option implements a similar functionaity.
When the script is run with -chunks x, the SNPs are broken into x different batches and those batches are run consecutively.

If your system has more than one processor you can take advantage of -threads and -processes for higher performance.
snipar processes chromosomes one by one each as one process. In presense of -processes x, x of the chromosomes will be
processesed at the same time. -threads specifies how many threads will be used for each of the chromosomes.
With -processes x -threads y at each point there will be x*y threads open.
This number shouldn't be much higher than available physical threads on the system.

Imputed parental genotypes and other informations about the imputation will be written to a file in HDF5 format for each chromosome.
You can see information about the outputs :ref:`here <imputed_file>`.

Family-based genome-wide association analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Family-based GWAS is performed by the gwas.py script. 
This script will estimate direct effects, non-transmitted coefficients, and population effects of input genetic variants
on the phenotype specified in the :ref:`phenotype file <phenotype>`. (If multiple phenotypes are present in the :ref:`phenotype file <phenotype>`,
the phenotype to analyse can be specified using the '--phen_index' argument, where '--phen_index 1' corresponds to the first phenotype.)


The script will use both :ref:`observed <observed genotypes>` and :ref:`imputed parental genotypes <imputed_file>` to estimate these effects. 
Note that if no imputed parental genotypes are input, gwas.py will estimate effects using individuals with both parents genotyped only,
provided that a :ref:`pedigree file <pedigree>` is also input. 
(A pedigree input is not needed when inputting :ref:`imputed parental genotypes <imputed_file>`.)

By default, for each variant, the script performs a regression of an individual's phenotype onto their genotype,
their (imputed/observed) father's genotype, and their (imputed/observed) mother's genotype. This estimates
the direct effect of the variant, and the paternal and maternal non-transmitted coefficients (NTCs). See
Young et al. [ref] 2022 for more details. 

If no parental genotypes are observed, then the imputed maternal & paternal genotypes become perfectly correlated.
In this case, to overcome collinearity, gwas.py will perform a regression of an individual's phenotype onto their genotype,
and the imputed sum of their parents' genotypes. This will estimate the direct effect of the SNP, and
the average NTC. 

If one wishes to model indirect genetic effects from siblings, one can use the '--fit_sib' option to add the genotype(s)
of the individual's sibling(s) to the regression. 

The gwas.py script first estimates a variance component model that models the phenotypic correlation between siblings, 
then does a transformation that allows the SNP effects to be estimated by simple linear regression while
accounting for correlations between siblings. 

The script outputs summary statistics in both gzipped :ref:`text format <_sumstats_text>` and
:ref:`HDF5 format <sumstats_hdf5>`.

Estimating correlations between effects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:ref:`correlate.py <correlate.py>` script estimates the genome-wide correlation between direct and population effects,
and between direct effects and average non-transmitted coefficients (NTCs). 
It takes as input the :ref:`summary statistics <_sumstats_text>` files output by :ref:`gwas.py <gwas.py>`
and LD-scores for the SNPs (as output by :ref:`ibd.py <ibd.py>` or by LDSC). 
It applies a method-of-moments based estimator that 
accouts for the known sampling variance-covariance of the effect estimates, and for the correlations
between effect estimates of nearby SNPs due to LD. (See Young et al. 2022 [ref] for more details.)

Note that this is different to genetic correlation as estimated by LDSC. LDSC attempts to use LD-scores to estimate
heritability and to separate out the contribution of population stratification. This estimator only uses
LD-scores to account for correlations between nearby SNPs, not to separate out population stratification. 
This is because we are (potentially) interested in the contribution of population stratification to population effects,
and whether direct effects are meaningfully different from population effects for particular phenotypes. 

Family-based polygenic score analyses
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Polygenic scores based on observed/imputed genotypes can be calculated and analysed using the :ref:`pgs.py <pgs.py>` script.

The :ref:`pgs.py <pgs.py>` takes similar inputs to the :ref:`gwas.py <gwas.py>` script. 
The main addition is that in order to compute a PGS, a :ref:`weights file <weights>` must be provided. 

By default, if no :ref:`phenotype file <phenotype>` is provided, the :ref:`pgs.py <pgs.py>` will compute
the PGS values of the all the genotyped individuals for whom :ref:`observed <observed genotypes>` or :ref:`imputed parental genotypes <imputed_file>`
parental genotypes are available. The script will output a :ref:`PGS file <pgs_file>`, 
including the imputed/observed PGS values for each individual's parents, facilitating family-based polygenic score analyses. 

If the '--fit_sib' argument is provided, the :ref:`PGS file <pgs_file>` will include a column corresponding to the average PGS value of the individual's sibling(s). 

To estimate the direct and population effects as well as the non-transmitted coefficients (NTCs) of the polygenic score, 
input a :ref:`phenotype file <phenotype>` to :ref:`pgs.py <pgs.py>`. 
One can first compute the PGS and write it to :ref:`file <pgs_file>`, 
and then use this as input to :ref:`pgs.py <pgs.py>` along with a :ref:`phenotype file <phenotype>`.

The direct effect and NTCs of the PGS are estimated as fixed effects in a linear mixed model that includes
a random effect that models (residual) phenotypic correlations between siblings. The population effect is estimated
from a separate linear mixed regression model that includes only the proband PGS as a fixed effect. 
The estimates and their standard errors are output to :ref:`file <pgs effects>` along with a separate
:ref:`file <pgs_vcov>` giving the sampling variance-covariance matrix of the direct effect and NTCs. 