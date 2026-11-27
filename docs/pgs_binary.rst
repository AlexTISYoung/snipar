.. _binary_pgs:

Binary-outcome PGS analyses in R
================================

Use ``pgs.R`` to analyse a binary phenotype, coded 0 or 1, using precomputed
polygenic scores. The script fits one- and two-generation logistic mixed
models with a family random intercept and automatically estimates the ratio
of their proband PGS coefficients, together with its standard error.

Both models use the same complete-case sample. The one-generation proband
coefficient is called the population effect and the two-generation proband
coefficient is called the direct effect. These coefficients are on the
log-odds scale; see :ref:`Interpreting the ratio <binary_pgs_ratio>` before
comparing them with results from linear PGS models.

Run an analysis
---------------

The script requires R and the R package ``lme4``. Install the package from
within R:

    ``install.packages("lme4")``

First obtain a :ref:`PGS file <pgs_file>`, for example using ``pgs.py`` with
observed or imputed parental genotypes. To fit both models with covariates, run

    ``pgs.R results/case_status --pgs scores.pgs.txt --phenofile phenotypes.txt --phen case_status --covar covariates.txt --covar_cols age,sex,PC1,PC2``

``results/case_status`` is the output prefix; ``--out results/case_status``
is also accepted. The output directory is created if necessary. Omit the
covariate options if no covariates are needed. No additional flag is required
to calculate the direct/population ratio.

If the PGS file contains a combined ``parental`` score, add ``--parsum``:

    ``pgs.R results/case_status_parsum --pgs scores.pgs.txt --phenofile phenotypes.txt --phen case_status --parsum``

``--parsum`` can also sum separate ``paternal`` and ``maternal`` scores.
Add ``--fit_sib`` to include a ``sibling`` score in the two-generation model.
To use the standalone script from a source checkout, replace ``pgs.R`` in
the commands above with ``Rscript snipar/scripts/pgs.R``. This requires
only R, ``lme4`` and the input files.

.. _binary_pgs_inputs:

Prepare the input files
-----------------------

All inputs are whitespace-delimited tables with headers and columns named
``FID`` and ``IID``. Rows are matched by the pair of IDs. IDs are read as
strings, preserving leading zeros, and each pair must occur at most once in
each input file. ``FID`` also defines the families used for the random
intercept and sampling covariance calculation.

A PGS file with separate parental scores has this format (illustrative rows):

.. code-block:: text

   FID IID proband paternal maternal
   001 00101 0.42 0.31 0.12
   001 00102 -0.08 0.31 0.12
   002 00201 0.17 -0.10 0.38

Additional columns such as ``FATHER_ID`` and ``MOTHER_ID`` are ignored.
Required scores must be numeric. Without ``--parsum``, both ``paternal`` and
``maternal`` are required. With ``--parsum``, an existing ``parental`` column
is used, or the two parental columns are summed. ``--fit_sib`` requires an
additional numeric ``sibling`` column.

The phenotype file contains the binary outcome:

.. code-block:: text

   FID IID case_status
   001 00101 1
   001 00102 0
   002 00201 1

Select a phenotype with ``--phen NAME`` or ``--phen_index N``, using only one
of these options. The index is one-based among columns other than ``FID``
and ``IID``; the default is the first phenotype. After filtering, the
phenotype must contain both 0 and 1 and no other values.

An optional covariate file contains numeric columns after the IDs:

.. code-block:: text

   FID IID age sex PC1 PC2
   001 00101 30 0 -0.02 0.11
   001 00102 32 1 0.01 0.09
   002 00201 28 0 0.03 -0.04

By default all non-ID columns in this file are included. Use
``--covar_cols age,sex,PC1`` to select a subset. Encode categorical
covariates as numeric indicator columns beforehand. Covariate names must not
reuse score or model names such as ``proband``, ``parental`` or ``phenotype``.

The analysis starts with the PGS-file rows. Unmatched rows and rows with a
missing or non-finite required phenotype, score or covariate are excluded
from both models. Missing parental scores therefore also reduce the
one-generation sample. ``NA``, ``NaN``, ``nan`` and ``N/A`` are recognized
missing values; ``--missing_char VALUE`` adds another marker. Missing IDs,
duplicate IDs or column names, and undeclared nonnumeric values cause an
error instead of being silently removed.

Models and score scaling
-------------------------

Both models include an intercept, selected covariates and a random intercept
grouped by ``FID``. The one-generation model uses the proband score. The
two-generation model adds separate paternal and maternal scores, or their
sum with ``--parsum``, and optionally the sibling score with ``--fit_sib``.
Families are assumed independent.

Coefficients are log odds per unit of the supplied score. With
``--scale_pgs``, every included score column is divided by the proband score's
standard deviation in the complete-case sample, giving all generations a
common scale. Covariates retain their supplied units and the binary
phenotype is not standardized.

The complete-case analysis sample must contain at least three families, at
least one family with multiple participants, and a proband score with
positive finite standard deviation. Both fixed-effect designs must have full
rank: constant or collinear scores and covariates cause an error.

Fitting uses ``lme4::glmer``, the Laplace approximation (``nAGQ=1``) and the
``bobyqa`` optimizer. ``--maxfun N`` changes the optimizer evaluation limit
from its default of 200000. Convergence or likelihood-check errors stop the
script before result files are written. A random-intercept
variance estimated at zero is supported and is recorded in the diagnostics.

.. _binary_pgs_outputs:

Inspect the results
--------------------

Start with ``OUT.effects.txt`` for all coefficients and the ratio, and
``OUT.model_info.tsv`` for sample counts, ratio status, and model diagnostics.
The estimates table is tab-separated with columns ``statistic``, ``est``, and
``SE``. ``1.proband`` gives the population coefficient, ``2.proband`` gives the
direct coefficient, and ``direct_population_ratio`` gives their ratio.
Other fixed effects follow the same ``MODEL.TERM`` naming convention.

.. list-table:: Files written for output prefix OUT
   :header-rows: 1
   :widths: 42 58

   * - Suffix
     - Contents
   * - ``.effects.txt``
     - All fixed effects and the direct/population ratio in columns ``statistic``, ``est``, and ``SE``. Model coefficients are named ``1.TERM`` and ``2.TERM``; the ratio is named ``direct_population_ratio``.
   * - ``.1.vcov.txt``, ``.2.vcov.txt``
     - Sampling covariance matrices for all fixed effects.
   * - ``.1_2.cross_vcov.txt``
     - Cross-model sampling covariance, with one-generation terms as rows and two-generation terms as columns.
   * - ``.direct_population.vcov.txt``
     - Two-by-two sampling covariance matrix for the population and direct proband coefficients.
   * - ``.1.model_vcov.txt``, ``.2.model_vcov.txt``
     - Model-based ``lme4::vcov`` matrices for comparison with the primary sampling covariance estimates.
   * - ``.model_info.tsv``
     - Phenotype, sample and case/control counts, exclusions, score scale, family counts, random-intercept SD, boundary indicator, convergence and likelihood diagnostics, and ``ratio_status`` and ``ratio_method`` for both models.

In ``.model_info.tsv``, ``n`` and ``n_families`` count the shared
complete-case individuals and families. ``ratio_method`` is
``Laplace_family_cluster_sandwich_delta``. The rows and columns of
``.direct_population.vcov.txt`` are ordered ``population``, then ``direct``.

Tables are tab-delimited with headers and use ``NA`` for missing results.
Covariance matrices label both rows and columns. The primary coefficient
SEs, within-model and cross-model covariance matrices, and ratio SE all use
the same family-cluster sandwich calculation with a :math:`G/(G-1)`
correction for :math:`G` families. The ``.model_vcov.txt`` files provide a
separate model-based comparison and are not used to compute the ratio SE.

.. _binary_pgs_ratio:

Interpreting the ratio
----------------------

The ratio is the two-generation proband coefficient divided by the
one-generation proband coefficient, using the same score weights and sample
in both fits. For example, coefficients of 0.2 and 0.4 give a ratio of 0.5.
This is a ratio of **log-odds coefficients**, not a ratio of odds ratios.
Because logistic models are noncollapsible, adding parental scores changes
the interpretation of the coefficient comparison; the ratio does not have
the same decomposition as the corresponding ratio from linear models.

The standard error uses the delta method and accounts for sampling
covariance between the direct and population coefficients. Their joint
covariance is obtained from stacked family contributions to the
Laplace-approximated likelihood. For interior fits this includes uncertainty
in the estimated random-intercept SD; at the boundary (SD at most ``1e-4``)
the calculation uses the fixed-effect block.

Check ``ratio_status`` in ``OUT.model_info.tsv``:

.. list-table:: Ratio status values
   :header-rows: 1
   :widths: 35 65

   * - Status
     - Interpretation
   * - ``ok``
     - The ratio and delta-method SE were computed without a denominator or covariance flag.
   * - ``weak_denominator``
     - The population coefficient is less than 1.96 sandwich SEs from zero in absolute value. The ratio and SE are retained, but the delta approximation can be unreliable near zero.
   * - ``undefined_zero_denominator``
     - The population coefficient is zero or non-finite, so the ratio and SE are missing.
   * - ``invalid_covariance``
     - The joint sampling covariance is non-finite or a marginal variance is non-positive. The ratio may still be reported, but its SE is missing.

The script also prints a warning for a non-``ok`` status. Review
``.model_info.tsv`` alongside the ratio, particularly the retained sample,
family count, boundary indicator and likelihood diagnostics.

All command-line options
-------------------------

Run ``pgs.R --help`` for the full usage message.
In addition to the required output prefix, ``--pgs`` and ``--phenofile``, the
available options are ``--phen``, ``--phen_index``, ``--covar``,
``--covar_cols``, ``--parsum``, ``--fit_sib``, ``--scale_pgs``,
``--missing_char`` and ``--maxfun``, as described above.

Regression checks
------------------

For development checks, install ``numDeriv`` in addition to ``lme4`` and run
from the source checkout:

    ``Rscript snipar/tests/test_pgs_binary.R``

The checks compare likelihood derivatives and sampling influences with
independent numerical calculations, exercise the zero-variance boundary,
and verify the CLI's matched samples, identifier handling, covariance
outputs, parental sums and sibling effects.
