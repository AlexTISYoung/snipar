.. _tutorial:

========
Tutorial
========

Tutorial on inferring IBD between siblings, imputing missing parental genotypes, and performing family-based GWAS and polygenic score analyses. 
Before working through the tutorial, please first install the package and read the :ref:`guide <guide>`. 

**family-based GWAS can be performed without imputed parental genotypes**

**snipar will meta-analyse siblings and trios by default when imputed parental genotypes are not provided**

**this tutorial guides through the entire workflow but gives examples of analyses without imputed genotypes**

Test data
--------------------

If *snipar* has been installed succesfully, the :ref:`command line scripts <scripts>` should be accessible as
executables in your terminal. A script that should be accessible loads the tutorial example data into a specified directory.
To create a directory called 'example_data/' in the current directory and load the example data into it, use the command:

    ``snipar_example_data.py --dest example_data``

You can create the example data directory elsewhere by changing the --dest argument. Please change your working directory to example_data/:

    ``cd example_data``

In this directory, there is some example data. 
The file phenotype.txt is a :ref:`phenotype file <phenotype>` containing a simulated phenotype with direct, paternal, and maternal effects, where 80% of the phenotypic
variance is explained by the combined direct, paternal and maternal effects of the SNPs; and the
pairwise correlations between the direct, paternal, and maternal effects are 0.5. 

The genotype data has been simulated so that there are 3000 independent families, where 1000 have two siblings but no parents genotyped,
1000 have one parent genotyped and a 50% chance of having a genotyped sibling, and the final 1000 have both parents genotyped and a 50%
chance of having a genotyped sibling. The example data includes :ref:`observed genotype data <observed genotypes>` formatted in both PLINK .bed format (chr_1.bed) and phased genotype
data in .bgen format (chr_1.bgen with associated sample file chr_1.sample). The folder also contains partial data—for example, chr_1_trios_sibs.bgen, which includes individuals with 
complete parental genotypes as well as individuals without parents but with siblings. This allows users to apply various estimators to different types of data.

Inferring IBD between siblings
------------------------------

The first step is to infer the identity-by-descent (IBD) segments shared between siblings.
*snipar* contains a script, :ref:`ibd.py <ibd.py>`, that employs a Hidden Markov Model (HMM) to infer the IBD segments for the sibling pairs.
The per-SNP genotyping error probability will be inferred from parent-offspring pairs when available;
alternatively, a genotyping error probability can be provided using the :code:`--p_error` option. By default, SNPs with
genotyping error rates greater than 0.01 will be filtered out, but this threshold can be changed with the :code:`--max_error` argument.
To infer the IBD segments from the genotype data in chr_1.bed,use the following command

    ``ibd.py --bed chr_@ --king king.kin0 --agesex agesex.txt --out chr_@ --threads 4 --ld_out``

This will output the IBD segments to a :ref:`gzipped text file <ibd_segments_file>` chr_1.ibd.segments.gz. 
Genotype files split over multiple chromosomes can be specified
using '@' as a numerical wildcard character: see :ref:`here <multichrom>`. 
In this example, :code:`--bed chr_@` instructs ibd.py to search for .bed files
chr_1.bed, chr_2.bed, ..., chr_22.bed, where each bed file contains SNPs from the numbered chromosome. 
In this case, only one bed file is in example_data/, chr_1.bed. 
If bed files for multiple chromosomes are found, IBD will be inferred separately for each chromosome, with one
output file per chromosome, with the chromosome number filling in the numerical wildcard in the --out argument. 
Alternatively, you can specify the chromosomes that you want to analyse: for example, you can include :code:`--chr_range 1-10`
for chromosome 1-10, or :code:`--chr_range 1 3 5-22` for chromosome 1, 3, and 5-22.

The :code:`--king` argument requires the address of the :ref:`KING kinship file <kinship>`, 
and the :code:`--agesex` argument requires the address of the :ref:`agesex file <agesex>`.
Age and sex information is needed to determine which is the mother/father in a parent-offspring relation.

The algorithm requires a genetic map to compute the probabilities of transitioning between different IBD states. 
If the genetic map positions (in cM) are provided in the .bim file, the script will use these. 
Alternatively, the ``--map`` argument allows the user to specify a genetic map in the same format as used by SHAPEIT 
(https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#formats) an example of which is 
provided in genetic_map.txt. 

If no genetic map is provided, then the deCODE sex-averaged map on GRCh38 coordinates (Halldorsson, Bjarni V., et al. "Characterizing mutagenic effects of recombination through a sequence-level genetic map." Science 363.6425 (2019).),
which is distributed as part of *snipar*, will be used. 

The algorithm computes LD scores of SNPs in order to account for correlations between SNPs. 
The ``--ld_out`` argument writes the LD scores to file in the same format as LDSC (https://github.com/bulik/ldsc). 

The user can also input a phased .bgen file. For example, to infer IBD from chr_1.bgen using the genetic map in genetic_map.txt, use this command:

    ``ibd.py --bgen chr_@ --king king.kin0 --agesex agesex.txt --out chr_@ --threads 4 --ld_out --map genetic_map.txt``

If the user has a :ref:`pedigree file <pedigree>`, they can input that instead of the *--king* and *--agesex* arguments. 
Siblings are inferred as individuals in the pedigree that share both parents. 
Using the example pedigree in pedigree.txt, you can infer IBD using this command:

    ``ibd.py --bed chr_@ --pedigree pedigree.txt --map genetic_map.txt --out chr_@ --threads 4 --ld_out``

Imputing missing parental genotypes
-----------------------------------

This is performed using the :ref:`impute.py <impute.py>` script. 
To impute the missing parental genotypes without using phase information, use this command:

    ``impute.py --ibd chr_@.ibd --bed chr_@ --king king.kin0 --agesex agesex.txt --out chr_@ --threads 4``

The script constructs a pedigree from the output of KING's relatedness inference (king.kin0),
and age and sex information (agesex.txt). 
The pedigree along with the IBD segments shared between siblings recorded in chr_1.ibd.segments.gz are used to impute missing parental genotypes
from the observed sibling and parental genotypes in chr_1.bed. 
The imputed parental genotypes are output to a :ref:`HDF5 file <imputed_file>`, chr_1.hdf5. 

If phased haplotypes are available in .bgen format, the imputation can use these as input, which improves the accuracy of the imputation. 
To perform imputation from the phased .bgen file in example_data/, use the following command:

    ``impute.py --ibd chr_@.ibd --bgen chr_@ --king king.kin0 --agesex agesex.txt --out chr_@ --threads 4``

As with the ibd.py script, the impute_runner.py script can use a user input :ref:`pedigree file <pedigree>` (with the *--pedigree* argument) rather than the *--king* and *--agesex* arguments.

Family-based GWAS with imputed parental genotypes
-------------------------------------------------

This is performed using the :ref:`gwas.py <gwas.py>` script, which implements different family-based GWAS designs.

The default design with imputed parental genotypes is described in Young et al. (2022) (https://www.nature.com/articles/s41588-022-01085-0).
This regresses the proband's phenotype jointly onto the proband's genotype, the father's imputed or observed genotype, and the mother's imputed or observed  genotype.
The regression produces variant level summary statistics on the direct genetic effect of the proband's genotype, the non-transmitted coefficient (NTC) for the father's genotype, and the NTC for the mother's genotype.
To use this design, supply the imputed parenta genotypes (above) to the gwas.py scrict. For example:

    ``gwas.py phenotype.txt --bed chr_@ --imp chr_@ --chr_range 1 --cpus 1 --out chr_@_young``

This regression is done using a linear mixed model that models phenotypic correlations between siblings,
where sibling relations are stored in the :ref:`output of the imputation script <imputed_file>`. 
The 'sibling variance component' output is the phenotypic variance explained by mean differences between sibships. 
:code:`--cpus` allows you to distribute computation across several processes to speed up analyses.

To use the .bgen file instead, use this command:

    ``gwas.py phenotype.txt --bgen chr_@ --imp chr_@ --cpus 1 --out chr_@_young``

The gwas.py script enables the use of different regression designs described in Guan et al. (https://www.nature.com/articles/s41588-025-02118-0).
For the homogeneous ancestry samples typically used in GWAS, you can increase statistical power for estimation of direct genetic effects by inclusion
of individuals without genotyped first-degree relatives (i.e., singletons) in the analysis. This also produces estimates of population effects (as estimated by standard GWAS)
that are almost identical to performing a standard GWAS in a linear mixed model framework. This is why we called this approach the 'unified estimator' in Guan et al. 
We give an example command here: 

    ``gwas.py phenotype.txt --bgen chr_@_trios_singletons --imp chr_@ --cpus 1 --impute_unrel --out chr_@_unified``

The ``--impute_unrel`` flag instructs *snipar* to linearly impute parental genotypes of singletons and include them into the analysis.

Approaches relying on imputed parental genotypes can be biased in strongly structured (Fst>0.01) and/or admixed samples. 
In Guan et al., we develop a 'robust' estimator that maximises power in such samples without introducing bias. 
Here's an example command invoking the robust estimator: 

    ``gwas.py phenotype.txt --bgen chr_@ --imp chr_@ --cpus 1 --robust --out chr_@_robust``

By default, the script outputs summary statistics in a :ref:`gzipped text file <sumstats_text>`: chr_1.sumstats.gz;
In addition to the text summary statistics, :ref:`HDF5 format summary statistics <sumstats_hdf5>` are also output to chr_1.sumstats.hdf5.
Alternatively, you can specify the output filename using the ``--out`` command: for example, with ``--out chr_@_X``, the script will
output the results to chr_1_X.sumstats.gz and chr_1_X.sumstats.hdf5; if '@' is not in the output suffix, e.g., ``--out gwas``, the results will
be stored in gwas_chr_1.sumstats.gz and gwas_chr_1.sumstats.hdf5.

Now we have estimated SNP effects. To compare the Young et al. (or robust or unified) estimates to the true effects, run
    
    ``python estimate_sim_effects.py chr_1_young.sumstats.hdf5 phenotype.effects.txt``

This should print estimates of the bias of the effect estimates.

The bias estimates for direct, paternal NTCs, maternal NTCs, and average NTCs should not be statistically significantly different from 
zero (with high probability). Population effects (as estimated by standard GWAS) are biased estimates of direct effects for this simulated 
phenotype because they also include indirect genetic effects and other confounding factors not modelled here. 

Family-based GWAS without imputed parental genotypes
----------------------------------------------------

Family-based GWAS can also be performed without imputed parental genotypes. In this case, only probands with genotypes for both parents and/or siblings available will be used.
In order to do this, one must provide a pedigree to gwas.py, as in:

    ``gwas.py phenotype.txt --out trios_sibs --bgen chr_@_trios_sibs --pedigree pedigree.txt --cpus 1``

With the above command, the script will default to meta-analysing samples with both parents genotyped (trios) and samples with siblings but without both parents genotyped. 
Alternatively, users can supply one of the following two options (``--robust`` is not applicable since it requires information derived from the imputation procedure):

- ``--sib_diff``: individuals with sibling genotypes will be used in a 'sib-GWAS' design using genetic differences between siblings, and those without will not be considered for the analysis;
- ``--impute_unrel``: individuals with both parents' genotypes and singletons will be used; individuals with sibling genotypes but incomplete parental genotypes will be ignored.

For example:

    ``gwas.py phenotype.txt --out sibs --bgen chr_@_trios_sibs --pedigree pedigree.txt --cpus 1 --sib_diff``

    ``gwas.py phenotype.txt --out trios_sibs_singletons --bgen chr_@_trios_sibs_singletons --pedigree pedigree.txt --cpus 1 --impute_unrel``

.. _correlations_tutorial:

Correlations between effects
----------------------------

How closely do population effects track direct effects across SNPs? Are they
also similar in magnitude? We can answer both questions with the family-GWAS
results from the earlier step (``chr_1_young.sumstats.gz``) and the LD scores
from ``ibd.py`` (``chr_1.l2.ldscore.gz``). In the ``example_data`` directory, run

    ``correlate.py chr_@_young effect --ldscores chr_@ --save_delete``

The first argument, ``chr_@_young``, identifies the summary-statistic files.
The script replaces ``@`` with each chromosome number. The second argument,
``effect``, is your choice of output prefix. This run creates five files
beginning with ``effect_`` in the current directory. To save them under
another name and directory, replace ``effect``
with, for example, ``results/height``. The estimates will then be in
``results/height_corrs.txt``.

Check the input filenames printed in the terminal. In this example they are
``chr_1_young.sumstats.gz`` and ``chr_1.l2.ldscore.gz``. The filtering summary
shows how many SNPs were retained and gives the reasons for any exclusions.
The tutorial summary statistics contain direct, paternal and maternal estimates
and their pairwise sampling correlations, so the script selects ``parental``
mode automatically. It fits their covariance jointly, averages the two parental
NTCs, and adds that average to the direct effect to obtain the population
effect. Every statistic below uses the same fitted covariance matrix.

The terminal shows all 17 estimates alongside their standard errors, starting
with the three statistics reported in the revised manuscript. To retrieve
those three results from ``effect_corrs.txt`` by row name, run

    ``awk 'NR == 1 || $1 == "r_direct_population" || $1 == "r_direct_avg_NTC" || $1 == "sqrt_var_ratio_population_direct"' effect_corrs.txt``

Here is an example from a run on the 1,000 retained tutorial SNPs. The first
two rows describe alignment and the third describes relative magnitude.

.. list-table::
   :header-rows: 1

   * - Statistic
     - Estimate
     - SE
   * - ``r_direct_population``
     - 0.895
     - 0.024
   * - ``r_direct_avg_NTC``
     - 0.477
     - 0.107
   * - ``sqrt_var_ratio_population_direct``
     - 1.532
     - 0.075

Direct and population effects track one another closely across these SNPs.
Direct effects also correlate positively with average NTCs, as expected in
this simulation. The ratio of 1.532 means the population effects are about
53% larger on the root-mean-square scale. The order of the names matters:
``sqrt_var_ratio_direct_population`` is the reciprocal (0.653 here), and
``var_ratio_direct_population`` is the square of that reciprocal (0.426).

The estimates use all retained SNPs. To calculate their standard errors,
the analysis omits each of 200 genomic blocks in turn and refits the weights
and covariance matrix. Each output file answers a different question.

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - File
     - How to read it
   * - ``effect_corrs.txt``
     - The 17 effect statistics and their standard errors. Look up ``r_direct_population`` to see how closely direct and population effects track each other across SNPs.
   * - ``effect_moments.tsv``
     - The fitted variances and covariances and their standard errors. For example, find ``var_paternal`` to see how much the paternal NTC varies across SNPs in this simulation.
   * - ``effect_sampling_covariance.tsv``
     - A matrix showing how estimates vary together across jackknife refits. For the sampling covariance between the two reported correlations, read the ``r_direct_population`` row at the ``r_direct_avg_NTC`` column. Diagonal entries give squared standard errors.
   * - ``effect_diagnostics.json``
     - Input and fit details. Check the retained SNP count and input mode, then look for warnings about boundary fits or undefined standard errors.
   * - ``effect_delete.txt``
     - One row for each omitted ``block``, saved by ``--save_delete``. Read down the ``r_direct_population`` column to see how the correlation changes as each block is left out. The other statistics and the fitted input moments appear alongside it.

You can also compare the fit with raw moments that use the same weight for
every effect pair at a SNP. Add ``--common-weight`` to the command above to
include these estimates in the moments, sampling covariance and deletion tables.

If you skipped the LD-score output in the IBD step, replace ``--ldscores chr_@``
with ``--bed chr_@`` to calculate LD scores from the tutorial genotypes.
See :doc:`correlate` for alternative input formats and :ref:`correlation_outputs`
for the full output definitions.

Polygenic score analyses
------------------------

The :ref:`pgs.py <pgs.py>` script can compute PGSs for individuals and their
parents and estimate their effects on a phenotype in the same run.

The tutorial includes ``direct_weights.txt``, where the weights are the true
direct SNP effects used to simulate the phenotype. This file uses LDpred column
names, which we specify explicitly. After the imputation step above, run:

    ``pgs.py direct --bed chr_@ --imp chr_@ --weights direct_weights.txt --SNP sid --beta_col ldpred_beta --A1 nt1 --A2 nt2 --phenofile phenotype.txt``
    
This writes the family scores to ``direct.pgs.txt``. The ``--phenofile`` option
also fits both models and writes ``direct.1.effects.txt``, ``direct.2.effects.txt``
and the other regression results in the current directory. Parental scores use observed
genotypes where available and imputed genotypes otherwise. If you skipped
imputation, replace ``--imp chr_@`` with ``--pedigree pedigree.txt``; scores
will then be computed only for individuals with both parents genotyped.

If you already computed ``direct.pgs.txt`` without supplying a phenotype,
run the following command from ``example_data`` to produce the regression files
from those saved scores:

    ``pgs.py direct --pgs direct.pgs.txt --phenofile phenotype.txt``

By default, this fits two linear mixed models on the same individuals. The
one-generation model estimates the population effect from the individual's
own PGS. The two-generation model includes parental scores to estimate the
direct effect and parental NTCs. Both account for residual similarities
between siblings.

Start by inspecting these files:

* ``direct.1.effects.txt``: the ``proband`` row gives the population effect.
* ``direct.2.effects.txt``: the ``proband`` row gives the direct effect, and
  the ``paternal`` and ``maternal`` rows give the parental NTCs.
* ``direct.direct_population_ratio.txt``: the ``ratio`` column gives the
  direct effect divided by the population effect, and ``se`` gives its
  standard error, accounting for covariance between the estimates.

Each effects file has three columns without a header: coefficient name,
estimate, and standard error. The ratio file has a header; for example, view
it with ``cat direct.direct_population_ratio.txt``. Its ``n`` column reports
the common sample size. Check that ``status`` is ``ok`` before interpreting
the ratio. A ``weak_denominator`` means the population coefficient is too
close to zero for the usual ratio approximation to be reliable.

In this example, the direct coefficient should be close to 1 because the
score uses the true direct effects. The population coefficient also reflects
the association with parental indirect effects. The ratio helps compare these
two estimates; it is a ratio of PGS regression coefficients, rather than the
SNP-effect variance ratio produced by ``correlate.py``.

The script also writes covariance matrices; see :ref:`PGS output files
<pgs_effects>` and :ref:`pgs_ratio` for details. To fit only one model, add
``--gen_models 1`` or ``--gen_models 2``; a single-model run does not produce
the ratio files. If parental scores are imputed from siblings alone, use
``--parsum`` to fit their sum and estimate the average NTC. A precomputed
file containing only a ``parental`` score column selects this automatically.

For a longer exercise exploring indirect effects and assortative mating,
continue to the :ref:`simulation exercise <simulation>`.

Binary outcomes
~~~~~~~~~~~~~~~

The tutorial's ``phenotype.txt`` contains a quantitative phenotype. For your
own binary outcome, coded 0/1, use :doc:`pgs.R <pgs_binary>` with the same
precomputed scores and a separate binary phenotype file. The R script fits
logistic mixed models with a family random intercept and reports a
direct/population ratio on the log-odds coefficient scale.
