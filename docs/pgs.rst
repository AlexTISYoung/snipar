.. _pgs_analysis:

Polygenic score analysis
========================

``pgs.py`` computes scores for individuals and their relatives, then fits
linear mixed models for a quantitative phenotype. By default it fits both a
one-generation model and a two-generation model on the same complete cases,
and estimates the ratio of the direct to population PGS coefficient, including
its sampling standard error.

For a binary phenotype coded 0/1, use the :doc:`pgs_binary` workflow with
``pgs.R``. It reuses scores computed by ``pgs.py`` and fits logistic mixed
models with a family random intercept.

Compute and reuse scores
------------------------

To compute scores from observed genotypes, imputed parental genotypes, and a
:ref:`weights file <weights>`:

    ``pgs.py scores --bed chr_@ --imp chr_@ --chr_range 1-22 --weights weights.txt``

Supply genotype file prefixes without their extensions. ``@`` is replaced by
chromosome numbers; ``--chr_range`` accepts ranges and individual chromosomes.
Use ``--bgen`` instead of ``--bed`` for phased BGEN input. The weights file has
a header containing ``SNP``, ``b``, ``A1``, and ``A2`` by default; change these
names with ``--SNP``, ``--beta_col``, ``--A1``, and ``--A2``. Weights are defined
with respect to ``A1``. Whitespace is the default separator; use ``--sep`` for
another delimiter.

If parental genotypes are observed and no imputation files are supplied,
provide a :ref:`pedigree <pedigree>` instead:

    ``pgs.py scores --bed combined_chr_@ --weights weights.txt --pedigree pedigree.txt``

In this mode, scores are computed only for individuals with both parents
genotyped. ``--batch_size`` controls the number of variants read at a time
(default 10000); ``--threads`` defaults to 1.

The result is ``scores.pgs.txt``. To analyse it:

    ``pgs.py analysis --pgs scores.pgs.txt --phenofile phenotype.txt``

Alternatively, add ``--phenofile`` to the score-computation command. A
precomputed ``--pgs`` input requires a phenotype file. Use either the
weights/genotype workflow or the precomputed-score workflow in a command.

PGS input and score columns
---------------------------

The generated :ref:`PGS file <pgs_file>` is space-delimited, with a header
beginning ``FID IID FATHER_ID MOTHER_ID`` followed by score columns. The
precomputed-score reader also accepts files beginning ``FID IID`` without
parental IDs; however, parental IDs/status are needed for ``--bpg`` and for
estimating the parental correlation used by ``--h2f`` when ``--rk`` is absent.
Use a single space between header fields, as in the generated file.

``FID`` identifies a full-sibling family and ``IID`` identifies the proband.
In generated files, a parental ID is ``NA`` when that parent is not genotyped;
its score can still be available through imputation. Recognised missing score
values include ``nan``, ``NA``, and ``N/A``.

.. list-table:: Score columns
   :header-rows: 1
   :widths: 25 75

   * - Column
     - Meaning
   * - ``proband``
     - The individual's own score.
   * - ``paternal``, ``maternal``
     - Separate observed or imputed parental scores.
   * - ``parental``
     - Sum of the paternal and maternal scores, used when only their sum is available.
   * - ``sibling``
     - Mean score of the proband's siblings; computed with ``--fit_sib``.
   * - ``gpp``, ``gmp``
     - Paternal grandfather and paternal grandmother scores.
   * - ``gpm``, ``gmm``
     - Maternal grandfather and maternal grandmother scores.

``--compute_controls`` additionally writes
``.pgs.control_paternal.txt``, ``.pgs.control_maternal.txt``, and
``.pgs.control_sibling.txt`` files for control families. Compute these scores
first, then analyse the required score file with a separate ``--pgs`` command.

Choose models and samples
-------------------------

``--gen_models`` selects models by number, for example ``--gen_models 2`` or
``--gen_models 1-3``. The default is ``1-2``.

.. list-table:: Regression models
   :header-rows: 1
   :widths: 15 45 40

   * - Model
     - Score predictors
     - Interpretation
   * - 1
     - ``proband``
     - Population PGS coefficient.
   * - 2
     - ``proband`` and parental scores; optionally ``sibling``.
     - Direct PGS coefficient and parental non-transmitted coefficients (NTCs).
   * - 3
     - Model 2 predictors plus grandparental scores.
     - Proband, parental, and grandparental coefficients conditional on the included generations.

Every model includes an intercept, any supplied ``--covar`` columns, and a
random effect for residual similarity among siblings. ``--grm`` adds a
:ref:`genetic relationship matrix <GRM>` to model other relatedness;
``--sparse_thresh`` defaults to 0.05.

``--parsum`` uses the sum of paternal and maternal scores in models 2 and 3.
Its coefficient estimates the average parental NTC in model 2. When reading a
precomputed file containing ``parental`` and no ``paternal`` column, the script
selects this option automatically. When computing such scores and fitting a
phenotype in the same command, specify ``--parsum`` explicitly.

To include sibling effects, use ``--fit_sib`` both when computing scores and
when fitting models 2 or 3. Model 1 still uses only the proband score and any
covariates.

Phenotypes and covariates are matched to scores by ``IID`` alone, so each
``IID`` must be unique across families and consistent across the input files.
``FID`` defines sibling groups but is not part of this matching key. Use
``--phen`` to select a case-sensitive phenotype name, or ``--phen_index`` for
its one-based position (default 1). ``--missing_char`` defaults to ``NA``.
``--bpg`` restricts the analysis to individuals with both parents genotyped.
``--scale_phen`` scales the phenotype to unit variance; ``--scale_pgs`` divides
all score columns by the proband-score standard deviation. For regression,
scaling occurs after matching phenotype IDs, before complete-case filtering.

The covariate reader requires a row for every phenotyped individual; missing
IDs cause an error. It mean-centres covariate columns and fills missing entries
with zero, equivalent to imputing the column mean on the original scale.
After this preprocessing, rows with non-finite phenotype values, required
scores, or remaining non-finite covariates are excluded from fitting.
When models 1 and 2 are requested together, model 2 is fitted first,
and model 1 uses exactly the same individuals. A standalone model 1 can
therefore use more observations than model 1 in a paired analysis. Model 3 has
its own complete-case selection; it is not included in the direct/population
ratio calculation.

.. _pgs_ratio:

Direct/population coefficient ratio
-----------------------------------

The default paired analysis estimates

.. math::

   R = \widehat\beta_{\mathrm{direct}} / \widehat\beta_{\mathrm{population}},

where the numerator is the ``proband`` coefficient from model 2 and the
denominator is the ``proband`` coefficient from model 1. Both models use the
same phenotype, scaling, and complete cases. No additional option is needed:

    ``pgs.py analysis --pgs scores.pgs.txt --phenofile phenotype.txt --covar covariates.txt --grm relatives.seg --gen_models 1-2``

The delta-method standard error accounts for the sampling covariance between
the two coefficients. The cross-model covariance uses the shared residual
from the full, two-generation model, assuming omitted parental contributions
in the population residual are orthogonal to that residual. Each model's
own marginal coefficient variances are retained. As with the model's GLS
standard errors, uncertainty in estimated variance components is not included.

The ``direct_population_ratio`` row in ``analysis.effects.txt`` contains the
ratio in ``est`` and its standard error in ``SE``. The component estimates
are the ``1.proband`` and ``2.proband`` rows in the same table.
The tab-delimited ``analysis.1_2.cross_vcov.txt`` gives the covariance
between models, with model 1 terms as rows and model 2 terms as columns.
Its ``proband`` row, ``proband`` column entry is the covariance between the
population and direct estimates. The corresponding entries in
``analysis.1.vcov.txt`` and ``analysis.2.vcov.txt`` give their variances.

The terminal reports the ratio status and warns when it is not ``ok``.
The Python return value also retains the status and method
(``GLS_shared_full_model_residual_delta``); these are not additional columns
in the three-column estimates table. Check the reported status before
interpreting the ratio:

.. list-table:: Ratio status values
   :header-rows: 1
   :widths: 40 60

   * - Status
     - Meaning
   * - ``ok``
     - Ratio and delta-method standard error were computed.
   * - ``weak_denominator``
     - The absolute population coefficient is less than 1.96 times its standard error. A ratio and standard error are reported, but the delta approximation can be unreliable near zero.
   * - ``undefined_zero_denominator``
     - The population coefficient is zero or non-finite; ratio and standard error are ``nan``.
   * - ``invalid_direct_coefficient``
     - The direct coefficient is non-finite; ratio and standard error are ``nan``.
   * - ``invalid_covariance``
     - The joint covariance contains non-finite values or a non-positive marginal variance. The standard error is ``nan``; the ratio can still be reported.
   * - ``inconsistent_marginal_covariance``
     - The joint covariance matrix has a materially negative eigenvalue. The standard error is ``nan``; the ratio can still be reported.

Running only model 1 or model 2 writes that model's rows to the estimates
table, without the ``direct_population_ratio`` row or paired-model covariance
files.

Effect and covariance files
---------------------------

``analysis.effects.txt`` collects all requested estimates in three
tab-separated columns with a header: ``statistic``, ``est``, and ``SE``.
Each coefficient is named ``N.TERM``, where ``N`` is the generation model.
For example, ``1.proband`` is the population effect and ``2.proband`` is the
direct effect; ``2.paternal`` and ``2.maternal`` give separate NTCs, or
``2.parental`` gives the combined parental-score coefficient with ``--parsum``.
Intercepts and covariates follow the same naming convention. Requested model 3
coefficients begin with ``3.``. The ratio row is named
``direct_population_ratio``.

``analysis.N.vcov.txt`` remains the sampling variance-covariance matrix for
model ``N``, with coefficient names labelling rows and columns. Separate
per-model effect tables and a separate ratio estimate file are no longer written.

Assortative mating adjustments
------------------------------

There are two distinct adjustment steps.

**When computing scores**, the script estimates the correlation between
parents' scores using observed parental and sibling scores under equilibrium.
It adjusts imputed parental scores only when the estimated correlation
divided by its standard error is greater than 1.5. ``--no_am_adj`` skips this
step. Reading an existing ``--pgs`` file does not repeat the adjustment.

.. note::

   In this version, ``--force_am_adj`` is accepted by the command-line parser
   but is not passed to the score-adjustment routine. It currently does not
   override the 1.5 threshold.

**After fitting model 2**, ``--h2f estimate,SE`` requests equilibrium
adjustment of the two-generation results using a heritability estimate from
MZ-DZ comparisons, RDR, or sibling realised relatedness. Supply the parental
score correlation and its standard error with ``--rk estimate,SE``, or let
the script estimate it from individuals with both parents genotyped. An
external ``--rk`` is needed when no such observed parental pairs are available.
``--rk`` alone has no effect; model 2 and ``--h2f`` must also be requested.
``--no_am_adj`` does not disable this post-fit step.

    ``pgs.py adjusted --pgs scores.pgs.txt --phenofile phenotype.txt --h2f 0.42,0.02 --rk 0.10,0.01``

The adjusted estimates are added to ``adjusted.effects.txt`` with names
``am_adj.delta``, ``am_adj.alpha``, ``am_adj.rk``, ``am_adj.k``, ``am_adj.r``,
``am_adj.h2_eq``, ``am_adj.rho``, ``am_adj.alpha_delta``, and
``am_adj.v_eta_delta``. They use the same ``statistic``, ``est``, and ``SE``
columns as the regression coefficients and ratio; no separate adjusted-estimate
file is written. ``am_adj.delta`` and ``am_adj.alpha`` are the standardised
direct coefficient and average parental NTC. ``am_adj.rho`` is the
direct/population ratio predicted from assortative mating alone, distinct from
the empirical ``direct_population_ratio`` row in the same table.
See the :ref:`simulation exercise <simulation>` for interpretation of the
equilibrium quantities.

Grandparental scores and model 3
--------------------------------

Add ``--grandpar`` when computing scores to append ``gpp``, ``gmp``, ``gpm``,
and ``gmm``. These are populated for individuals with both parents genotyped.
The script reuses the parents' parental scores when those parents have rows
in the score data; otherwise it imputes grandparental scores linearly from
the observed parent. Other individuals receive missing grandparental scores.
This step requires separate ``paternal`` and ``maternal`` columns and pedigree
information. With ``--no_am_adj``, the correlation used for this linear
imputation is set to zero.

Fit model 3 from a file containing grandparental scores with, for example:

    ``pgs.py three_gen --pgs grandparent_scores.pgs.txt --phenofile phenotype.txt --gen_models 3 --gparsum``

``--gparsum`` combines ``gpp + gmp`` into ``gp`` (paternal grandparents) and
``gpm + gmm`` into ``gm`` (maternal grandparents). Without it, model 3 fits
four separate grandparental coefficients. The requested score predictors
must contain enough independent variation to fit the model: scores imputed
only as linear functions of parental scores cannot identify additional
grandparental coefficients on their own.

For all command-line options, see :ref:`pgs.py <pgs.py>` or run
``pgs.py --help``.
