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

.. code-block:: bash

   pgs.py scores --bed chr_@ --imp chr_@ --chr_range 1-22 --weights weights.txt

Supply genotype file prefixes without their extensions. ``@`` is replaced by
chromosome numbers; ``--chr_range`` accepts ranges and individual chromosomes.
Use ``--bgen`` instead of ``--bed`` for phased BGEN input. The weights file has
a header containing ``SNP``, ``b``, ``A1``, and ``A2`` by default; change these
names with ``--SNP``, ``--beta_col``, ``--A1``, and ``--A2``. Weights are defined
with respect to ``A1``. Whitespace is the default separator; use ``--sep`` for
another delimiter.

If parental genotypes are observed and no imputation files are supplied,
provide a :ref:`pedigree <pedigree>` instead:

.. code-block:: bash

   pgs.py scores --bed combined_chr_@ --weights weights.txt --pedigree pedigree.txt

In this mode, scores are computed only for individuals with both parents
genotyped. ``--batch_size`` controls the number of variants read at a time
(default 10000); ``--threads`` defaults to 1.

The result is ``scores.pgs.txt``. To analyse it:

.. code-block:: bash

   pgs.py analysis --pgs scores.pgs.txt --phenofile phenotype.txt

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

.. code-block:: bash

   pgs.py analysis --pgs scores.pgs.txt --phenofile phenotype.txt \
       --covar covariates.txt --grm relatives.seg --gen_models 1-2

The delta-method standard error accounts for the sampling covariance between
the two coefficients. The cross-model covariance uses the shared residual
from the full, two-generation model, assuming omitted parental contributions
in the population residual are orthogonal to that residual. Each model's
own marginal coefficient variances are retained. As with the model's GLS
standard errors, uncertainty in estimated variance components is not included.

Three additional tab-delimited files are written:

.. list-table:: Paired-model output files
   :header-rows: 1
   :widths: 45 55

   * - Suffix
     - Contents
   * - ``.direct_population_ratio.txt``
     - Header and one result row with the ratio, its standard error, coefficient estimates, covariance, sample size, and status.
   * - ``.1_2.cross_vcov.txt``
     - Cross-model covariance matrix; rows are model 1 terms and columns are model 2 terms, including intercept and covariates.
   * - ``.direct_population.vcov.txt``
     - The 2-by-2 covariance matrix for the population and direct proband coefficients, in that order.

The ratio file has columns ``ratio``, ``se``, ``direct``, ``population``,
``direct_se``, ``population_se``, ``covariance``, ``sampling_correlation``,
``n``, ``status``, and ``method``. ``n`` is the number of individuals used by
both models. The method is recorded as
``GLS_shared_full_model_residual_delta``.

Check ``status`` before interpreting the ratio:

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

The script emits a warning for every status other than ``ok``. Running only
model 1 or model 2 produces that model's effect files without the paired-model
ratio or cross-covariance files.

Effect and covariance files
---------------------------

For each requested generation ``N``, ``analysis.N.effects.txt`` contains
three tab-delimited columns without a header: coefficient name, estimate,
and standard error. It includes the intercept, covariates, and score terms.
``analysis.N.vcov.txt`` is the corresponding sampling variance-covariance
matrix with row and column names. For example, ``analysis.2.effects.txt``
contains the direct coefficient in its ``proband`` row and either separate
``paternal``/``maternal`` NTCs or a single ``parental`` coefficient.

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

.. code-block:: bash

   pgs.py adjusted --pgs scores.pgs.txt --phenofile phenotype.txt \
       --h2f 0.42,0.02 --rk 0.10,0.01

``adjusted.am_adj_pars.txt`` has columns ``parameter``, ``estimate``, and
``SE``. Its rows are ``delta``, ``alpha``, ``rk``, ``k``, ``r``, ``h2_eq``,
``rho``, ``alpha_delta``, and ``v_eta_delta``. Here ``delta`` and ``alpha``
are the standardised direct coefficient and average parental NTC; ``rho``
is the direct/population ratio predicted from assortative mating alone. It
is different from the empirical paired-model ``ratio`` described above.
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

.. code-block:: bash

   pgs.py three_gen --pgs grandparent_scores.pgs.txt --phenofile phenotype.txt \
       --gen_models 3 --gparsum

``--gparsum`` combines ``gpp + gmp`` into ``gp`` (paternal grandparents) and
``gpm + gmm`` into ``gm`` (maternal grandparents). Without it, model 3 fits
four separate grandparental coefficients. The requested score predictors
must contain enough independent variation to fit the model: scores imputed
only as linear functions of parental scores cannot identify additional
grandparental coefficients on their own.

For all command-line options, see :ref:`pgs.py <pgs.py>` or run
``pgs.py --help``.
