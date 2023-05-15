.. _simulation:
===================
Simulation Exercise
===================

Exercise on simulating data using the simulate.py script and performing family-based GWAS and polygenic score analysis. 

Simulating data
--------------------

If *snipar* has been installed succesfully, the :ref:`command line scripts <scripts>` should be accessible as
executables in your terminal. *snipar* includes a script for simulating genotype-phenotype data according to 
different scenarios: direct and indirect genetic effects, with and without assortative mating. 
To simulate data, please first create a directory to store the data:

    ``mkdir sim``

Now, we are going to simulate data for 3000 families genotyped at 1000 independent SNPs. We are going to simulate 20 generations of assortative mating with parental phenotype correlation 0.5. 

    ``simulate.py 1000 0.5 sim/ --nfam 3000 --impute --n_am 20 --r_par 0.5``

Please change your working directory to sim/:

    ``cd sim``

In this directory, the file phenotype.txt is a :ref:`phenotype file <phenotype>` containing the simulated phenotype. 

The genotype data (chr_1.bed) has been simulated so that there are 3000 independent families, each with two siblings genotyped. 

Inferring IBD between siblings
------------------------------

The first step is to infer the identity-by-descent (IBD) segments shared between siblings. 
However, for the purpose of this simulation exercise (where SNPs are independent, so IBD inference doesn't work)
we have provided the true IBD states in the file chr_1.segments.gz.


Imputing missing parental genotypes
-----------------------------------

This is performed using the :ref:`impute.py <impute.py>` script. 
To impute the missing parental genotypes without using phase information, use this command:

    ``impute.py --ibd chr_@ --bed chr_@ --pedigree pedigree.txt --out chr_@ --threads 4``

The pedigree along with the IBD segments shared between siblings recorded in chr_1.segments.gz are used to impute missing parental genotypes
from the observed sibling and parental genotypes in chr_1.bed. 
The imputed parental genotypes are output to a :ref:`HDF5 file <imputed_file>`, chr_1.hdf5. 

Family based GWAS
-----------------

This is performed using the :ref:`gwas.py <gwas.py>` script. 
To compute summary statistics for direct effects, non-transmitted coefficients (NTCs), and population effects for the SNPs in the .bed file, use this command:

    ``gwas.py phenotype.txt --bed chr_@ --imp chr_@ --threads 4``

This takes the observed genotypes in chr_1.bed and the imputed parental genotypes in chr_1.hdf5 and uses
them to perform, for each SNP, a joint regression onto the proband's genotype, the father's (imputed/observed) genotype, and the mother's
(imputed/observed) genotype. This is done using a linear mixed model that models phenotypic correlations between siblings,
where sibling relations are stored in the :ref:`output of the imputation script <imputed_file>`. 
The 'family variance estimate' output is the  phenotypic variance explained by mean differences between sibships, 
and the residual variance is the remaining phenotypic variance. 

Polygenic score analyses
------------------------

In addition to family based GWAS, *snipar* provides a script (:ref:`pgs.py <pgs.py>`) for computing polygenic scores (PGS) based on observed/imputed genotypes,
and for performing family based polygenic score analyses. 
Here, we give some examples of how to use this script. The script computes a PGS
from a :ref:`weights file <weights>`. 
For the tutorial, we provide a weights file (direct_weights.txt) in `LD-pred <https://github.com/bvilhjal/ldpred>`_ format
where the weights are the true direct genetic effect of the SNP. 

To compute the PGS from the true direct genetic effects, use the following command:

    ``pgs.py direct --bed chr_@ --imp chr_@ --weights causal_effects.txt --beta_col direct``
    
It outputs the PGS to a :ref:`PGS file <pgs_file>`: direct.pgs.txt. 

To estimate direct effect and average NTC of the PGS, use the following command:

    ``pgs.py direct --pgs direct.pgs.txt --phenofile phenotype.txt``

This will output a population effect estimate (1 generation model) to direct.1.effects.txt, and 
direct effect and average NTC estimates to (2 generation model) to direct.2.effects.txt.  

To compute the PGS from the true direct genetic effects+estimation error, use the following command:

    ``pgs.py direct_v1 --bed chr_@ --imp chr_@ --weights causal_effects.txt --beta_col direct_v1``
    
It outputs the PGS to a :ref:`PGS file <pgs_file>`: direct_v1.pgs.txt. 

To estimate direct effect and average NTC of the PGS, use the following command:

    ``pgs.py direct_v1 --pgs direct_v1.pgs.txt --phenofile phenotype.txt``

This will output a population effect estimate (1 generation model) to direct_v1.1.effects.txt, and 
direct effect and average NTC estimates to (2 generation model) to direct_v2.2.effects.txt.  

Unlike when using the true direct genetic effects as weights, the direct effect of the PGS estimated
from noisy weights will be smaller than the population effect. 

