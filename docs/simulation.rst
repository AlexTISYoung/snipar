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

Polygenic score analyses
------------------------

*snipar* provides a script (:ref:`pgs.py <pgs.py>`) for computing polygenic scores (PGS) based on observed/imputed genotypes,
and for performing family based polygenic score analyses. 
The script computes a PGS from a :ref:`weights file <weights>`. 

To compute the PGS using the true direct genetic effects as weights, use the following command:

    ``pgs.py direct --bed chr_@ --imp chr_@ --weights causal_effects.txt --beta_col direct``
    
It outputs the PGS to a :ref:`PGS file <pgs_file>`: direct.pgs.txt. The pgs computation script
automatically estimates the correlation between parents PGS values (also using full-sibling offspring PGS values to do this)
and performs an adjustment for assortative mating when using the imputed parental genotypes to
compute the PGS. 

To estimate direct effect and average NTC of the PGS, use the following command:

    ``pgs.py direct --pgs direct.pgs.txt --phenofile phenotype.txt``

This will output a population effect estimate (1 generation model) to direct.1.effects.txt, and 
direct effect and average NTC estimates to (2 generation model) to direct.2.effects.txt. The 
sampling variance-covariance matrix of the estimates is output to direct.1.vcov.txt (for the 1 generation model) and
direct.2.vcov.txt (for the 2 generation model).

As we are using the true direct effects as weights, the PGS captures all of the heritability,
and the direct and population effects should both be the same (1 in expectation), and the 
average parental NTC should be zero (in expectation). To check this, read in the 
effect estimate output files in *R* or look at them using a text viewer (e.g. less -S on a unix system).

To compute the PGS from the true direct genetic effects+estimation error (such as would be obtained from a GWAS), 
use the following command:

    ``pgs.py direct_v1 --bed chr_@ --imp chr_@ --weights causal_effects.txt --beta_col direct_v1``
    
It outputs the PGS to a :ref:`PGS file <pgs_file>`: direct_v1.pgs.txt. 

To estimate direct effect and average NTC of the PGS, use the following command:

    ``pgs.py direct_v1 --pgs direct_v1.pgs.txt --phenofile phenotype.txt``

This will output a population effect estimate (1 generation model) to direct_v1.1.effects.txt, and 
direct effect and average NTC estimates to (2 generation model) to direct_v2.2.effects.txt. 

Unlike when using the true direct genetic effects as weights, the direct effect of the PGS estimated
from noisy weights (in direct_v1.1.effects.txt) will be smaller than the population effect (direct_v2.2.effects.txt).
This is because the PGS does not capture all of the heritability due to estimation error in the weights. 
The PGS has its population effect inflated (relative to its
direct effect) by assortative mating, which induces a correlation of the PGS with the component of the heritability
not captured by the PGS due to estimation error. This inflation is not captured by the direct effect of the PGS
because chromosomes segregate independently during meiosis. (In this simulation, all causal SNPs segregate independently.) 
Here, the ratio between direct and population effects of the PGS should be around 0.87. 

One should also observe a statistically significant average parental NTC (in direct_v2.2.effects.txt) of the PGS from 
the two-generation model despite the absence of parental indirect genetic effects in this simulation. Here,
the ratio between the average NTC and the direct effect should be around 0.15. This demonstrates
that statistically significant average NTC estimates cannot be interpreted as unbiased estimates of 
parental indirect genetic effects, especially for phenotypes affected by assortative mating. 
