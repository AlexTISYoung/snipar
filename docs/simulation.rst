.. _simulation:
===================
Simulation Exercise
===================

Exercise on simulating data using the simulate.py script and performing family-based polygenic score analysis. 

Simulating data
--------------------

If *snipar* has been installed succesfully, the :ref:`command line scripts <scripts>` should be accessible as
executables in your terminal. *snipar* includes a script for simulating genotype-phenotype data according to 
different scenarios: direct and indirect genetic effects, with and without assortative mating. 
To simulate data, please first create a directory to store the data:

    ``mkdir sim``

Now, we are going to simulate data for 3000 families genotyped at 1000 independent SNPs. We are going to simulate 10 generations of assortative mating with parental phenotype correlation 0.5. 

    ``simulate.py 1000 0.5 sim/ --nfam 3000 --impute --n_am 10 --r_par 0.5 --save_par_gts``

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
population and direct effect estimates are the coefficients on the proband PGS in the 1 and 2
generation models, so are indicated by the 'proband' row. The average NTC estimate is the
coefficient on the parental PGS in the two-generation model. The first column gives the name
of the covariate/PGS column, the second column gives the estimated regression coefficient,
and the third column gives the standard error of the estimate. The sampling variance-covariance matrix of the estimates is output to direct.1.vcov.txt (for the 1 generation model) and
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
from noisy weights (in direct_v1.1.effects.txt) will be smaller than the population effect (direct_v1.2.effects.txt).
This is because the PGS does not capture all of the heritability due to estimation error in the weights. 
The PGS has its population effect inflated (relative to its
direct effect) by assortative mating, which induces a correlation of the PGS with the component of the heritability
not captured by the PGS due to estimation error. This inflation is not captured by the direct effect of the PGS
because chromosomes segregate independently during meiosis. (In this simulation, all causal SNPs segregate independently.) 
Here, the ratio between direct and population effects of the PGS should be around 0.87. 

One should also observe a statistically significant average parental NTC (in direct_v1.2.effects.txt) of the PGS from 
the two-generation model despite the absence of parental indirect genetic effects in this simulation. Here,
the ratio between the average NTC and the direct effect should be around 0.15. This demonstrates
that statistically significant average NTC estimates cannot be interpreted as demonstrating
parental indirect genetic effects, especially for phenotypes affected by assortative mating. 

Adjusting for assortative mating
--------------------------------

We now show how to adjust two-generation PGI results for assortative mating. 
To do this, we will combine the offspring and parental genotype files. 
This enables us to estimate the correlation between parents' scores 
using the observed parental genotypes. (This is better than using the sibling 
genotypes because the estimate from observed parental genotypes is uncorrelated with the PGI direct effect estimate.)

    ``plink --bfile chr_1 --bmerge chr_1_par --out chr_1_combined``

We now compute the noisy PGI using the observed offspring and parental genotypes:

    ``pgs.py direct_v1_obs --bed chr_@_combined --weights causal_effects.txt --beta_col direct_v1 --pedigree pedigree.txt``

To complete the inputs to two-generation PGI analysis, we need an estimate of heritability,
as one would obtain from sib-regression, RDR, MZ-DZ twin comparisons. This estimate is 
a downard biased estimate of the equilibrium heritability by a factor of (1-r), where
r is the correlation between the parents direct genetic effect components. We can obtain
this from VCs.txt output of the simulation. Each row gives, for each generation, 
the variance of the direct genetic effect component, the phenotypic variance, and The
correlation between parents direct genetic effect components. The equilibrium heritability is
obtained by using the values in the last row: 
dividing the variance of the direct genetic effect component (first column) by the phenotypic variance
(second column). To then obtain the heritability as estimated by sib-regression, RDR, MZ-DZ twin comparisons,
we multiply the equilibrium heritability by (1-r), where r is obtained from the third column of 
the last row. The equilibrium heritability should be around 0.58, and the heritability as estimated
by sib-regression, RDR, MZ-DZ twin comparisons should be around (1-0.30)*0.58=0.41. 

We can now adjust the two-generation PGI results for assortative mating using the following command:

    ``pgs.py direct_v1_obs --pgs direct_v1_obs.pgs.txt --phenofile phenotype.txt --h2f 0.42,0.01``

This script will take the input heritability estimate (0.42) and the standard error of the estimate (0.01)
and will use this to estimate the fraction of heritability the PGI would explain in a random mating population,
k, which should be around 0.5; the correlation between parents' direct genetic effect components, r, 
which should be around 0.30; the equilibrium heritability, which should be around 0.58; 
the ratio between direct and population effects that would be expected based on assortative mating alone, rho,
which should be around 0.85; the indirect effect of true direct effect PGI, alpha_delta, which should not be
statistically significantly different from zero because there are no parental indirect genetic effects in this simulation; 
and v_eta_delta, the contribution to the phenotypic variance from indirect genetic effects correlated with direct genetic effects,
which should also be statistically indistinguishable from zero. 

