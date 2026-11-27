.. _simulation:

===================
Simulation Exercise
===================

Exercise on simulating data using the simulate.py script and performing family-based polygenic score analysis. This assumes you are working in a bash environment on a unix based machine, but similar steps should work in different environments.
Activate the environment containing *snipar* before starting. The genotype merging and
sample removal steps also require `PLINK 1.9 <https://www.cog-genomics.org/plink/1.9/>`_
available as ``plink`` on your PATH. The final binary-phenotype example requires R
and the R package ``lme4``.

Simulating data
--------------------

If *snipar* has been installed successfully, the :ref:`command line scripts <scripts>` should be accessible as
executables in your terminal. *snipar* includes a script for simulating genotype-phenotype data according to 
different scenarios: direct and indirect genetic effects, with and without assortative mating. 
To simulate data, please first create a directory to store the data:

    ``mkdir sim``

Now, we are going to simulate data for 3000 families, each with two full-siblings, genotyped at 1000 independent SNPs. 
We simulate a phenotype affected by direct genetic effects and assortative mating. 
We are going to simulate 20 generations of assortative mating with parental phenotype correlation 0.5, reaching an approximate equilibrium. 
The command for this is:

    ``simulate.py 1000 0.5 sim/ --nfam 3000 --impute --n_am 20 --r_par 0.5 --save_par_gts``

where the first argument gives the number of causal SNPs, the second argument gives the 
random mating heritability, the third gives the output directory, --nfam gives the number of families, --impute 
tells the script to impute missing parental genotypes, --n_am gives the number of generations of assortative mating,
the parental phenotype correlation is given by --r_par, and --save_par_gts tells the script to output the 
genotypes of the parents of the final generation in addition to the genotypes of the final generation.

Please change your working directory to sim/:

    ``cd sim``

The genotype data (chr_1.bed) contains 3000 families, each with two siblings genotyped.
The supplied GRM includes the sibling and cousin relationships in the final
generation.

The pedigree file (pedigree.txt) contains the pedigree for all simulated generations. The pedigree has columns:
FID, IID, FATHER_ID, MOTHER_ID, SEX, PHENO, FATHER_PHENO, MOTHER_PHENO, DIRECT, FATHER_DIRECT, MOTHER_DIRECT.
The FID and IID columns are the family and individual IDs, respectively. The FATHER_ID and MOTHER_ID columns are the IDs of the parents of the individual.
The FATHER_PHENO and MOTHER_PHENO columns are the phenotype values of the parents, and the DIRECT, FATHER_DIRECT, and MOTHER_DIRECT columns are the direct genetic effect components of the individual, father, and mother, respectively.

To enable analysis of the final generation phenotypes alone, we have placed the phenotype values for the final generation in a separate file (phenotype.txt). 

We have also included a file (last_gen_rel.seg) that records the sibling and cousin relationships in the final generation to give an example of how a sparse genetic relatedness matrix, a :ref:`GRM <GRM>`,
can be provided to the :ref:`gwas.py <gwas.py>` and :ref:`pgs.py <pgs.py>` scripts. This tab or white space delimited text file has a row for each pairwise relationship with (minimally) columns
giving ID1, ID2, and a column named 'PropIBD' or 'relatedness' giving the relatedness coefficient. 
This is designed to work with the output of the `KING <https://www.kingrelatedness.com/manual.shtml>`_ IBD segment inference with the --ibdseg argument.
But other ways of computing relatedness coefficient can be used, such as from a pedigree. 

A GRM is not required to run family-GWAS or PGS analyses. By default, a linear mixed model that accounts for correlations 
between siblings will be used. By providing a GRM, correlations between other relatives can be modelled, 
leading to increased power and more accurate standard errors when other relative types are present in the analysis sample. 

Family-based GWAS without imputed parental genotypes
----------------------------------------------------

To perform a family-based GWAS, we use the :ref:`gwas.py <gwas.py>` script. 
To perform a family-based GWAS without imputed parental genotypes, use the following command:

    ``gwas.py phenotype.txt --pedigree pedigree.txt --bed chr_@ --grm last_gen_rel.seg --out chr_@_sibdiff --cpus 1``

The first argument is the phenotype file. As we are not inputting an imputed parental genotype file,
we must specify a pedigree file using the ``--pedigree`` argument. 
The genotype data in .bed format is specified by ``--bed chr_@`` argument.
The ``@`` symbol is a placeholder for the chromosome number, so the script will read the genotype data from ``chr_1.bed``. 
The output file is specified by ``--out chr_@_sibdiff``. The script will output summary statistics to a gzipped text file: ``chr_1_sibdiff.sumstats.gz``.
The ``--cpus`` argument can be used to specify the number of processes to use to parallelize the GWAS. 

Since the genotype data of the final generation contains 3000 sibling pairs, but not parental genotypes, the script will perform sib-GWAS 
using genetic differences between siblings to estimate direct genetic effects (see `Guan et al. <https://www.nature.com/articles/s41588-025-02118-0>`_).
The summary statistics are output to a gzipped text :ref:`sumstats file <sumstats_text>`: chr_1_sibdiff.sumstats.gz.

We can combine the final two generations' genotype data into one .bed file using PLINK 1.9:

    ``plink --bfile chr_1 --bmerge chr_1_par --out chr_1_combined``

If we run the GWAS script on the combined genotype data, we can estimate the direct genetic effects using the full-sibling offspring and parental genotypes 
in a trio design:

    ``gwas.py phenotype.txt --pedigree pedigree.txt --bed chr_@_combined --grm last_gen_rel.seg --out chr_@_trio``

The summary statistics are output to a gzipped text :ref:`sumstats file <sumstats_text>`: chr_1_trio.sumstats.gz.
If you read the summary statistics file (e.g. into R or using ``zless -S``) you can see that the effective sample size for 
direct genetic effects is substantially larger from the trio design than the sib-differences design. 
Note that both designs use the same number of phenotype observations in a generalized least-squares regression, but the trio design uses more information from the parents.
In this simulation, the effective sample size from the trio design should be about 45% larger than for the sib-differences design.

By default, the :ref:`gwas.py <gwas.py>` script will perform a meta-analysis of samples with genotyped siblings but without both parents genotyped —
using the sib-difference estimator — and samples with both parents genotyped — using the trio design — when
imputed parental genotypes are not provided. This should achieve something close to optimal power for family-GWAS
without imputed parental genotypes. However, improved power can be achieved by using designs that take advantage of
imputed parental genotypes. 

Inferring IBD between siblings
------------------------------

The first step in the imputation of missing parental genotypes from siblings is to infer the identity-by-descent (IBD) segments shared between siblings.
This is done using the :ref:`ibd.py <ibd.py>` script. However, for the purpose of this simulation exercise (where SNPs are independent, so IBD inference doesn't work)
we have provided the true IBD states in the file chr_1.segments.gz.

Imputing missing parental genotypes
-----------------------------------

This is performed using the :ref:`impute.py <impute.py>` script. 
To impute the missing parental genotypes without using phase information, use this command:

    ``impute.py --ibd chr_@ --bed chr_@ --pedigree pedigree.txt --out chr_@ --threads 4``

The pedigree along with the IBD segments shared between siblings recorded in chr_1.segments.gz are used to impute missing parental genotypes
from the observed sibling and parental genotypes in chr_1.bed. 
The imputed parental genotypes are output to a :ref:`HDF5 file <imputed_file>`, chr_1.hdf5. 

Family-based GWAS with imputed parental genotypes
-------------------------------------------------

To perform a family-GWAS with the parental genotypes imputed using the above command, use the following command:

    ``gwas.py phenotype.txt --bed chr_@ --imp chr_@ --grm last_gen_rel.seg --out chr_@_imp``

The script will output summary statistics to a gzipped text file: ``chr_1_imp.sumstats.gz``.
Since the genotype data of the final generation contains 3000 sibling pairs, the script will perform family-GWAS 
using parental genotypes imputed from the siblings (see `Guan et al. <https://www.nature.com/articles/s41588-025-02118-0>`_).
The summary statistics are output to a gzipped text :ref:`sumstats file <sumstats_text>`: chr_1_imp.sumstats.gz.

However, the imputation we performed above does not use phase information. Using phased data results in improved recovery of parental genoypes,
and therefore improved power, especially for imputation from parent-offspring pairs (not simulated here). 
The simulation script outputs parental genotypes imputed from phased data, which we can use instead for family-GWAS using the following command:

    ``gwas.py phenotype.txt --bed chr_@ --imp phased_impute_chr_@ --grm last_gen_rel.seg --out chr_@_phased``

The effective sample size should be slightly higher than when using the unphased imputation. For imputation 
involving parent-offspring pairs (not simulated here), the gain in power is much larger. 

Increasing power for family-GWAS by including singletons
--------------------------------------------------------

We now demonstrate the power increase from including singletons in family-GWAS
by demonstrating the unified estimator introduced by `Guan et al. <https://www.nature.com/articles/s41588-025-02118-0>`_. 
To demonstrate this, we first remove one sibling from half of the families in the final generation using the following command:

    ``plink --bfile chr_1 --remove <(head -n $(( $(wc -l < chr_1.fam) / 2 )) chr_1.fam | awk 'NR % 2 == 0 {print $1, $2}') --make-bed --out chr_1_singletons``

This produces a new .bed file (chr_1_singletons.bed) with one sibling removed for half of the families. 
We now perform the imputation of missing parental genotypes using this reduced dataset:

    ``impute.py --ibd chr_@ --bed chr_@_singletons --pedigree pedigree.txt --out chr_@_singletons --threads 4``

This means we have imputed parental genotypes for half of the families - where we have two siblings observed. 
We have singletons for the other half of the families and no imputed parental genotypes. 
We can now perform family-GWAS including the singletons by adding the ``--impute_unrel`` argument:

    ``gwas.py phenotype.txt --bed chr_@_singletons --imp chr_@_singletons --grm last_gen_rel.seg --impute_unrel --out chr_@_unified``

The output should say 

    ``3000 individuals with imputed but no observed parental genotypes.``
    ``1500 samples without imputed or observed parental genotypes will be included through linear imputation.``

This shows we are using the 1500 samples without any imputed or observed parental genotypes (singletons) through linear imputation.

We can compare this to performing family-GWAS without including the 1500 singletons:

    ``gwas.py phenotype.txt --bed chr_@_singletons --imp chr_@_singletons --grm last_gen_rel.seg --out chr_@_no_singletons``

The median effective N for direct genetic effects should be about 18% higher from the analysis including singletons 
(chr_1_unified.sumstats.gz) than from the analysis excluding singletons (chr_1_no_singletons.sumstats.gz).
The increase in effective sample size can reach up to 50% as the singleton sample size grows relative to the 
sample with observed parental genotypes or parental genotypes imputed from siblings and/or single parents. 

Polygenic score analyses
------------------------

*snipar* provides a script (:ref:`pgs.py <pgs.py>`) for computing polygenic scores (PGS) based on observed/imputed genotypes,
and for performing family based polygenic score analyses. The script computes a PGS from a :ref:`weights file <weights>`. 

To compute the PGS using the true direct genetic effects as weights, use the following command:

    ``pgs.py direct --bed chr_@ --imp chr_@ --weights causal_effects.txt --beta_col direct``
    
It outputs the PGS to a :ref:`PGS file <pgs_file>`: direct.pgs.txt. The PGS computation script
automatically estimates the correlation between parents' PGS values (also using full-sibling offspring PGS values to do this).
It adjusts imputed parental PGS values for assortative mating when the estimated correlation
divided by its standard error is greater than 1.5. Use ``--no_am_adj`` to skip this score adjustment;
see :ref:`PGS analysis <pgs_analysis>` for details.

To estimate direct effect and average NTC of the PGS, use the following command:

    ``pgs.py direct --pgs direct.pgs.txt --phenofile phenotype.txt --grm last_gen_rel.seg``

By default, this fits the 1 and 2 generation models on the same individuals with complete
phenotype and required PGS observations. All estimates are written to
``direct.effects.txt``, a tab-separated table with columns ``statistic``,
``est``, and ``SE``. ``1.proband`` is the population effect from model 1,
``2.proband`` is the direct effect from model 2, and ``2.parental`` is the
average NTC estimated using the combined parental score. Intercepts and
covariates use the same model-prefixed names. The sampling covariance matrices
remain in ``direct.1.vcov.txt`` and ``direct.2.vcov.txt``.

The ``direct_population_ratio`` row in the same table gives the direct
coefficient divided by the population coefficient, with a standard error
accounting for their sampling covariance. Inspect the estimates with:

    ``column -t direct.effects.txt``

Check the ratio status printed in the terminal before interpreting it. In particular,
``weak_denominator`` means that the population coefficient is close to zero and the
delta-method standard error can be unreliable. The ``proband`` row and ``proband``
column of ``direct.1_2.cross_vcov.txt`` give the covariance between the population
and direct estimates. See :ref:`direct/population ratio output <pgs_ratio>` for the complete output
and status descriptions.

As we are using the true direct effects as weights, the PGS captures all of the heritability,
and the direct and population effects should both be the same (1 in expectation), and the 
average parental NTC should be zero (in expectation). The direct/population ratio should
therefore be close to 1. To check the effect estimates, read in the
``direct.effects.txt`` table in *R* or look at it using a text viewer (e.g. less -S on a unix system).

To compute the PGS from the true direct genetic effects+estimation error (such as would be obtained from a family-GWAS), 
use the following command:

    ``pgs.py direct_v1 --bed chr_@ --imp chr_@ --weights causal_effects.txt --beta_col direct_v1``
    
It outputs the PGS to a :ref:`PGS file <pgs_file>`: direct_v1.pgs.txt. (Notice also that the inferred
correlation between parents' PGSs is lower than when using the true direct genetic effects as weights due to
estimation error in the weights.)

To estimate direct effect and average NTC of the PGS, use the following command:

    ``pgs.py direct_v1 --pgs direct_v1.pgs.txt --phenofile phenotype.txt --grm last_gen_rel.seg``

This writes ``direct_v1.effects.txt`` with the population effect in
``1.proband``, the direct effect in ``2.proband``, the average NTC in
``2.parental``, and the paired-model ratio in ``direct_population_ratio``.

Unlike when using the true direct genetic effects as weights, the direct effect of the PGS estimated
from noisy weights (``2.proband``) will be smaller than the population effect (``1.proband``).
This is because the PGS does not capture all of the heritability due to estimation error in the weights. 
The PGS has its population effect inflated (relative to its direct effect) by assortative mating, 
which induces a correlation with the component of the heritability not directly captured by the PGS due to estimation error. 
This inflation is not captured by the direct effect of the PGS because the within-family variation used to estimate the direct effect is due to the random segregation of genetic material during meiosis.
Here, the ratio between direct and population effects of the PGS should be around 0.86.
Read the estimate and its standard error with:

    ``column -t direct_v1.effects.txt``

The ``2.proband`` and ``1.proband`` rows show the two coefficients used in the ratio.
You can choose a single model with ``--gen_models 1`` or ``--gen_models 2``, but the
ratio is produced only when both models are requested in the same run. A standalone
one-generation analysis can retain individuals whose parental PGS is missing; the paired
analysis keeps the two models' samples identical.

One should also observe a statistically significant average parental NTC (``2.parental`` in ``direct_v1.effects.txt``) of the PGS from
the two-generation model despite the absence of parental indirect genetic effects in this simulation. Here,
the ratio between the average NTC and the direct effect should be around 0.15. This demonstrates
that statistically significant average NTC estimates cannot be interpreted as automatically demonstrating
parental indirect genetic effects, especially for phenotypes affected by assortative mating. 

Adjusting indirect genetic effect estimates for assortative mating
------------------------------------------------------------------

We now show how to adjust two-generation PGI results for assortative mating
using the procedure outlined in `Estimation of indirect genetic effects and heritability under assortative mating <https://www.biorxiv.org/content/10.1101/2023.07.10.548458v1.abstract>`_. 
The estimation procedure is summarized in this diagram: 

.. image:: two_gen_estimation.png
    :scale: 30 %
    :align: center
    :alt: Two-generation estimation procedure accounting for assortative mating

The estimation requires as inputs: an estimate of the correlation between parents' scores, :math:`r_k`;
the regression coefficients from two-generation PGI analysis, (:math:`\delta_{\text{PGI}:k},\alpha_{\text{PGI}:k}`);
and a heritability estimate, :math:`h^2_f`,from MZ-DZ twin comparisons, `RDR <https://www.nature.com/articles/s41588-018-0178-9>`_, or sib-regression.

The estimation procedure outputs estimates of: :math:`k`, the fraction of heritability the PGI would explain in a random mating population;
:math:`r_\delta`, the correlation between parents' true direct genetic effect components; 
:math:`h^2_\text{eq}`, the equilibrium heritability, adjusting for the downward bias in heritability estimates from
MZ-DZ comparisons, RDR, and sib-regression; 
:math:`\alpha_\delta`, the indirect genetic effect of true direct genetic effect PGI;
and :math:`v_{\eta:\delta}`, the fraction of phenotypic variance contributed by the indirect genetic effect component
that is correlated with the direct genetic effect component. 

We can use *snipar* to compute the two-generation PGI estimates and the correlation between parents' scores, 
and we can input a heritability estimate into *pgs.py* script to complete the inputs, so that
*snipar* will perform the two-generation analysis adjusting for assortative mating. 

To perform the estimation, we will the combined offspring and parental genotype files. 
(This was created above using the command ``plink --bfile chr_1 --bmerge chr_1_par --out chr_1_combined``.)
This enables us to estimate the correlation between parents' scores 
using the observed parental genotypes. This is better than using the siblings 
because the correlation estimate from observed parental genotypes is uncorrelated with the PGS regression coefficients.

We now compute the noisy PGI using the observed offspring and parental genotypes:

    ``pgs.py direct_v1_obs --bed chr_@_combined --weights causal_effects.txt --beta_col direct_v1 --pedigree pedigree.txt``

To complete the inputs to two-generation PGI analysis, we need an estimate of heritability --
as one would obtain from sib-regression, RDR, or MZ-DZ twin comparisons. This estimate is 
a downard biased estimate of the equilibrium heritability, :math:`h^2_\text{eq}`, by a factor of :math:`(1-r_\delta)`, where
:math:`r_\delta` is the correlation between the parents' direct genetic effect components. 

We can obtain this from the VCs.txt output of the simulation, which can be read into R/Python/etc as table. 
Each row gives, for each generation, the variance of the direct genetic effect component, the phenotypic variance, and the
correlation between parents' direct genetic effect components. The equilibrium heritability is
obtained by using the values in the last row: 
dividing the variance of the direct genetic effect component (first column) by the phenotypic variance
(second column). To then obtain the heritability as estimated by sib-regression, RDR, and MZ-DZ twin comparisons,
we multiply the equilibrium heritability by :math:`(1-r_\delta)`, where :math:`r_\delta` is obtained from the third column of 
the last row. The equilibrium heritability should be around 0.59, and :math:`r_\delta` should be around 0.29, so the heritability as estimated
by sib-regression, RDR, MZ-DZ twin comparisons should be around :math:`h^2_f \approx (1-0.29) \times 0.59=0.42`. 

Compute this value from the final row of your own simulation output:

    ``h2f=$(awk 'END {print ($1 / $2) * (1 - $3)}' VCs.txt)``

We can now perform two-generation PGI analysis accounting for assortative mating:

    ``pgs.py direct_v1_obs --pgs direct_v1_obs.pgs.txt --phenofile phenotype.txt --h2f "$h2f,0"``

This script will take the input heritability estimate (approximately 0.42) and its standard error (here 0 since we used the simulated components)
to estimate the fraction of heritability the PGI would explain in a random mating population,
:math:`k`, which should be around 0.5; the correlation between parents' direct genetic effect components, :math:`r_\delta`, 
which should be around 0.29; the equilibrium heritability, :math:`h^2_\text{eq}`, which should be around 0.59; 
the ratio between direct and population effects that would be expected based on assortative mating alone, :math:`\rho_k`,
which should be around 0.86; the indirect genetic effect of true direct genetic effect PGI, :math:`\alpha_\delta`, which should not be
statistically significantly different from zero (with high probability) because there are no parental indirect genetic effects in this simulation; 
and :math:`v_{\eta:\delta}`, the contribution to the phenotypic variance from the indirect genetic effect component correlated with direct genetic effect component,
which should also be statistically indistinguishable from zero (with high probability).
These estimates are included in ``direct_v1_obs.effects.txt`` as rows beginning
with ``am_adj.``, for example ``am_adj.k``, ``am_adj.h2_eq``, and
``am_adj.alpha_delta``. They use the same ``est`` and ``SE`` columns as the
regression coefficients.

The ``am_adj.rho`` row is the ratio predicted from assortative mating alone.
The ``direct_population_ratio`` row in the same table is the empirical ratio
estimated from the paired regressions. These two estimates answer different questions
and can be compared in this simulation, where the indirect genetic effects are zero.
For another dataset, an externally estimated parental PGS correlation and its standard
error can be supplied with ``--rk estimate,SE`` alongside ``--h2f``. Otherwise, this
post-fit adjustment estimates the correlation using observed parental scores; it requires
individuals with both parents genotyped.

Using the scores for a binary phenotype
---------------------------------------

For a binary phenotype, the :doc:`pgs_binary` script fits one- and two-generation
logistic mixed models with a family random intercept. It uses the same complete
cases in both models and reports the direct/population coefficient ratio and
its standard error.
The standard errors cluster by family ID; cousin relationships between families
in ``last_gen_rel.seg`` are not modelled by this binary analysis.

You can reuse direct_v1.pgs.txt with a separate binary phenotype file containing
columns ``FID IID case_status``, where ``case_status`` is coded 0 or 1. Both
values must occur in the analysis sample. For this demonstration, create a
binary phenotype file by splitting the simulated continuous phenotype at its
median.

    ``Rscript -e 'p <- read.table("phenotype.txt", header=TRUE); p$case_status <- as.integer(p$phenotype > median(p$phenotype)); write.table(p[c("FID", "IID", "case_status")], "binary_phenotype.txt", quote=FALSE, row.names=FALSE)'``

With R and the R package ``lme4`` installed, run the installed ``pgs.R`` command:

    ``pgs.R direct_v1_binary --pgs direct_v1.pgs.txt --phenofile binary_phenotype.txt --phen case_status --parsum``

``--parsum`` uses the combined parental score in this example.
``direct_v1_binary.effects.txt`` collects all estimates in columns
``statistic``, ``est``, and ``SE``. The ``direct_population_ratio`` row divides
the ``2.proband`` log-odds coefficient by the ``1.proband`` log-odds coefficient.
The ``2.parental`` row gives the combined parental-score coefficient. The
terminal also shows the estimated family random-intercept variance and its
SE for each model. A larger variance means more variation across families
in baseline log odds after accounting for the fitted scores.
See :doc:`pgs_binary` for installation, input formats, diagnostics,
and guidance on interpreting the ratio.

.. _stabilizing_selection_tutorial:

Simulating stabilizing selection
--------------------------------

Breeders whose phenotypes are closer to a fixed optimum have a higher chance
of surviving to reproduce. To try this with 100 unlinked causal SNPs, run

    ``mkdir -p sim_stabilizing``

    ``simulate.py 100 0.5 sim_stabilizing/ --nfam 1500 --n_stabilizing 3 --selection_variance 0.25 --save_par_gts``

The positional arguments give the number of causal SNPs, the initial
heritability due to direct effects (0.5), and the output prefix. The trailing
slash puts the output in ``sim_stabilizing/``. Each offspring's phenotype is
its direct genetic value plus a random environmental contribution. Selection
acts on the breeders for three generations. The final sample contains two
siblings per family, measured before selection on their own phenotypes.

* ``--n_stabilizing 3`` sets the number of generations of selection. The
  default is zero. This mode models direct effects with random pairing among
  survivors. Leave ``--n_am``, ``--n_random``, ``--v_indir`` and ``--beta_vert``
  at zero.
* ``--selection_optimum`` sets the favoured phenotype relative to the initial
  mean genetic value. It defaults to zero and stays fixed across generations.
* ``--selection_variance 0.25`` sets the width of the survival curve in squared
  phenotype units. Smaller positive values give stronger selection. At 0.25, a
  candidate one phenotype unit from the optimum has about a 14% chance of
  surviving, compared with 100% at the optimum.
* ``--nfam 1500`` requests 1,500 final families. ``--selection_nfam`` sets the
  initial number of candidates of each sex (default four times ``--nfam``,
  or 6,000 here). Increase the pool if too few survive. With phased BGEN input
  (for example, ``--bgen phased_chr_@ --chr_range 1``), the haplotypes determine
  the pool size, so omit ``--selection_nfam``. SNP positions must use GRCh38
  coordinates covered by the bundled deCODE map. The simulator drops unmapped SNPs.
  For simulated genotypes, ``--maf`` fixes the allele frequency across SNPs
  and ``--min_maf`` sets the lower bound when frequencies are drawn (default
  0.05).
* ``--save_par_gts`` adds the final parental generation's PLINK files to the
  offspring genotypes saved by default.

Open ``sim_stabilizing/phenotype.txt`` for the 3,000 offspring phenotypes.
In ``selection_metadata.json``, compare ``phenotypic_variance`` before and
after selection within each generation to see how the survivors differ from
the candidate pool. ``VCs.txt`` reports genetic and phenotypic variance in
successive offspring generations.

For a neutral comparison, keep the same settings and use
``--selection_variance inf``. Every candidate survives in this run. The two
runs draw their own founders and environmental effects, so their realised
variances will vary.

    ``mkdir -p sim_neutral``

    ``simulate.py 100 0.5 sim_neutral/ --nfam 1500 --n_stabilizing 3 --selection_variance inf --save_par_gts``

Family GWAS in the final generation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Begin in the selected population. PLINK combines the observed offspring and
parental genotypes. The pedigree then links each phenotyped offspring to both
parents. The relatedness file lets the GWAS account for siblings in the sample.

    ``cd sim_stabilizing``

    ``plink --bfile chr_1 --bmerge chr_1_par --out chr_1_combined``

    ``gwas.py phenotype.txt --bed chr_@_combined --pedigree pedigree.txt --grm last_gen_rel.seg --out trio_@ --cpus 1``

The GWAS writes direct, paternal and maternal effect estimates to
``trio_1.sumstats.gz``. Use these estimates to compare direct and population
effects across SNPs. The simulated loci are unlinked, so give them equal
weight with ``--ld-weight none``. With 100 loci and 20 jackknife blocks, each
jackknife fit drops about five SNPs to estimate uncertainty.

    ``correlate.py trio_@ effects --ld-weight none --n_blocks 20``

    ``awk 'NR == 1 || $1 == "r_direct_population" || $1 == "r_direct_avg_NTC" || $1 == "sqrt_var_ratio_population_direct"' effects_corrs.txt``

Move to the neutral control and repeat the four commands from ``plink``
through ``awk``. Each directory will then contain its own ``effects_corrs.txt``.

    ``cd ../sim_neutral``

Interpreting the comparison
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Read ``sqrt_var_ratio_population_direct`` first. It is the standard
deviation of population effects divided by that of direct effects. In one
run the ratio was 0.78 under selection (SE 0.02) and 1.01 in the neutral
control (SE 0.03). A ratio below one means that population GWAS effects
vary less across SNPs than direct effects. Their directions can still agree
closely, as measured by ``r_direct_population``.

Stabilizing selection favours combinations of alleles whose effects offset
one another. In this simulation it can generate a negative correlation between
direct effects and the average non-transmitted coefficient
(``r_direct_avg_NTC``), even though the phenotype has only direct genetic
and environmental components. In the neutral control the true non-transmitted
component has zero variance, so its estimated correlation can fluctuate widely
or be undefined. See ``effects_diagnostics.json`` for boundary warnings when
interpreting jackknife SEs.

The larger simulations in the revised analysis by Tan et al. used 20
generations and five replicates. The mean population-to-direct SD ratio was
0.756 under strong stabilizing selection and 1.004 in the neutral control.
The smaller exercise illustrates the same comparison. Each simulation draws
a new sample, so your estimates will vary. See the
`Tan et al. preprint and version history <https://doi.org/10.1101/2024.10.01.24314703>`_.
