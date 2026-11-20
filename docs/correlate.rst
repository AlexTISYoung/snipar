.. _effect_correlations:

Correlations and variances of GWAS effects
==========================================

Use ``correlate.py`` to compare direct effects, average non-transmitted
coefficients (NTCs), and population effects across SNPs. The command estimates
one covariance matrix for the underlying effects, using each SNP's reported
standard errors and sampling covariances to account for estimation noise.
Every reported correlation, regression coefficient, and variance ratio follows
from this fit on the same set of SNPs.

The revised manuscript (20 November 2026, Supplementary Note §§2.1–2.3
and Supplementary Table 4) reports three quantities. The correlations
:math:`r(D,P)` and :math:`r(D,N)` describe how closely the direct effects
align across SNPs with population effects and average NTCs, respectively.
The ratio :math:`k_{P/D}=\mathrm{SD}(P)/\mathrm{SD}(D)` compares the
spread of population and direct effects. For example, :math:`k_{P/D}=1.2` means
that population effects have a standard deviation 20% larger than direct
effects. The slope of population effects regressed on direct effects is
:math:`r(D,P)\,k_{P/D}`.

Running an analysis
-------------------

With :ref:`GWAS summary statistics <sumstats_text>` from ``gwas.py`` and
LD scores from ``ibd.py`` or LDSC, run

.. code-block:: bash

   correlate.py chr_@_young effect --ldscores chr_@ --save_delete

Here ``chr_@_young`` matches files such as ``chr_1_young.sumstats.gz``, and
``chr_@`` matches files such as ``chr_1.l2.ldscore.gz``. The command searches
chromosomes 1–22 by default. Use ``--chr_range 1-5 7`` to select chromosomes.
A single plain-text or gzip file can also be supplied by its full name.

The default analysis downweights SNPs in stronger linkage disequilibrium
(LD) by the inverse of their LD score. LD tables need a SNP identifier and an
``L2`` column, with one row per SNP ID. To compute scores from PLINK genotypes
instead, supply one chromosome per BED/BIM/FAM file set. For example,

.. code-block:: bash

   correlate.py chr_@_young effect --bed chr_@ --ld_wind 1 --ld_out ld_chr@

The LD window defaults to 1 cM. This example also saves the scores to
``ld_chr1.l2.ldscore.gz``, etc. Both ``@`` and ``#`` are accepted in BED
and LD output prefixes. To give every SNP equal LD weight, use
``--ld-weight none`` in place of ``--ldscores`` or ``--bed``. The fit still
uses the sampling precision of each SNP.

.. _correlation_inputs:

Choosing the input effects
--------------------------

Summary statistics must include a SNP ID (``SNP``), base-pair position
(``pos``), and allele frequency (``freq``). Supply ``CHR`` or identify
chromosomes through the file pattern. For a single legacy file with no
available chromosome labels, put the rows in genomic order before running
the command. Use the same effect allele and phenotype units for all effect
estimates at a SNP. Common alternative column names, such as ``rsid``,
``BP``, and ``EAF``, are also accepted.

Supply a beta and a standard error for **each** input effect, and a sampling
correlation or sampling covariance for **each pair**. For example,
``r_direct_avg_NTC`` is the sampling correlation between the two estimates
at one SNP. The genome-wide correlation between the underlying effects is
estimated by the analysis.

.. list-table::
   :header-rows: 1
   :widths: 22 30 48

   * - ``--mode``
     - Input effects
     - Pairwise sampling correlations
   * - ``parental``
     - ``direct``, ``paternal``, ``maternal``
     - ``r_direct_paternal``, ``r_direct_maternal``, ``r_paternal_maternal``
   * - ``average``
     - ``direct``, ``avg_NTC``
     - ``r_direct_avg_NTC``
   * - ``direct-population``
     - ``direct``, ``population``
     - ``r_direct_population``

For each effect listed above the native columns are ``<effect>_Beta`` and
``<effect>_SE``. A ``cov_`` column can replace the corresponding ``r_``
column. With correlations, the sampling covariance is the reported correlation
multiplied by the two standard errors. Retained SNPs must have strictly
positive standard errors and positive definite sampling covariance matrices.

The default ``--mode auto`` uses the first mode in the table with all required
columns in every input file. Selection is based on headers. If parental columns
exist but are empty, select ``--mode average`` to use the reported average NTC.
In parental mode, the fitted paternal and maternal NTCs can differ, and their
mean defines the average NTC. For SNIPar models with a sibling effect, use the
reported average NTC, which includes a sibling contribution. Auto mode skips
the parental option when the sibling-effect column is present.

In all modes, the population effect is defined as direct effect plus average
NTC. A direct/NTC analysis can therefore ignore redundant population columns.
In direct-population mode the script first derives average NTC as population
minus direct and applies the same transformation to the sampling covariance.

Filtering and effect scale
--------------------------

The analysis uses complete effect vectors and their sampling covariances.
Duplicate summary-statistic SNP IDs raise an error by default.
``--duplicate-snps drop`` excludes every row for each repeated ID and records
those IDs in the diagnostics. With LD weighting, retained SNPs must have a
finite positive LD score. The default minimum minor allele frequency is 0.05
(``--min_maf``, inclusive). Sampling correlations more than six standard
deviations from their pair-specific mean are filtered by default. Use
``--corr_filter`` to change the threshold or ``--no-corr-filter`` to disable it.

For its genome-wide analysis, the revised manuscript retained autosomal SNPs
with MAF at least 0.05, positive finite LD scores, complete effect estimates
and nonsingular sampling covariances. It excluded repeated SNP identifiers
and outlying sampling correlations. Use ``--duplicate-snps drop`` to apply
that duplicate-ID rule. Check sampling covariances before running the fit,
since a singular sampling matrix on a retained SNP raises an error.

The default effect scale is per allele. With ``--scale genotype-sd``, the
script multiplies the effects at a SNP by :math:`\sqrt{2f(1-f)}`, where
:math:`f` is its allele frequency, and the sampling covariance by
:math:`2f(1-f)`. This describes effects per genotype standard deviation
and changes the across-SNP moments being estimated.

How the estimates fit together
------------------------------

Write :math:`D` for direct effect, :math:`N` for average NTC and
:math:`P=D+N` for population effect. The estimator fits the variances and
covariances jointly under a zero-mean working model for the SNP effects.
It subtracts each SNP's sampling covariance from the observed effect products
and uses feasible generalised least squares to weight those products by their
joint precision under a Gaussian model. It updates the weights as the fit
changes and constrains the signal covariance matrix to be positive
semidefinite. This ensures that every linear combination of its effects has
nonnegative fitted variance.

For example, the population variance satisfies

.. math::

   \operatorname{Var}(P) =
   \operatorname{Var}(D) + 2\operatorname{Cov}(D,N) + \operatorname{Var}(N).

Inverse LD scores account approximately for redundancy between SNPs. The
result describes marginal SNP effects and can include contributions from
population stratification. The optional ``--common-weight`` diagnostic averages
the sampling-noise-corrected products using the same scalar weight for every
entry at a SNP. Compare these raw moments with the fitted moments to see
how weighting and the covariance constraint affect the results.

.. _correlation_outputs:

Reading the results
-------------------

Start with ``effect_corrs.txt`` and ``effect_diagnostics.json``. The table
has columns ``correlation`` (statistic name), ``est`` and ``SE``. The first eight
rows summarize correlations, regression slopes and variances. Their estimates
and standard errors come from the joint fit and its genomic block deletions.

.. list-table::
   :header-rows: 1
   :widths: 48 52

   * - Statistic
     - Meaning
   * - ``r_direct_population``, ``r_direct_avg_NTC``
     - Correlation of direct effects with population effects or average NTCs.
   * - ``reg_population_direct``, ``reg_avg_NTC_direct``
     - Regression slope of the named effect on direct effects, its covariance with direct effects divided by the direct-effect variance.
   * - ``var_direct``, ``var_population``, ``var_avg_NTC``
     - Fitted signal variance for each effect.
   * - ``var_ratio_direct_population``
     - Direct-effect variance divided by population-effect variance.

Nine more rows give the three pairwise covariances and the six directed
standard deviation (SD) ratios among :math:`D,N,P`. For the manuscript's
three statistics, select ``r_direct_population``, ``r_direct_avg_NTC``,
and ``sqrt_var_ratio_population_direct``. The last of these is the
population-to-direct SD ratio :math:`k_{P/D}`. The historical
``var_ratio_direct_population`` row gives the variance ratio in the reverse
direction. For example, a direct/population variance ratio of 0.64 corresponds
to a population/direct SD ratio of :math:`1/\sqrt{0.64}=1.25`.

The other files use the same output prefix.

.. list-table::
   :header-rows: 1
   :widths: 48 52

   * - Suffix
     - Contents
   * - ``_moments.tsv``
     - Each distinct fitted variance and covariance, named in the ``moment`` column, with its estimate (``est``) and jackknife standard error (``SE``). Covers the fitted input effects and the direct, average NTC and population effects.
   * - ``_sampling_covariance.tsv``
     - One square jackknife covariance matrix for the 17 reported statistics and any additional fitted input moments. The first column (``statistic``) labels the rows. The diagonal gives squared standard errors and each off-diagonal entry gives the sampling covariance between two estimates.
   * - ``_diagnostics.json``
     - Input mode and scale, filtering counts, fit convergence and numerical backend, block sizes, boundary fits and undefined statistics.
   * - ``_delete.txt``
     - With ``--save_delete``, a tab-separated table of all leave-one-block-out estimates.

For example, the ``cov_direct_population`` row in ``_moments.tsv`` gives the
covariance of direct and population effects across SNPs. To find the sampling
variance of that estimate, look up ``cov_direct_population`` on both axes of
``_sampling_covariance.tsv``.

Add ``--common-weight`` to compare fitted and raw input moments. Look for
``common_weight_moment_var_direct`` alongside ``var_direct`` in ``_moments.tsv``,
for example. These additional estimates also appear on both axes of the
sampling covariance matrix. The diagnostics record whether the raw
covariance matrix is positive semidefinite.

Each row in ``_delete.txt`` gives the estimates obtained when the block named
in the first column (``block``) is omitted. The next 17 columns give the
statistics in the same order as the rows in ``_corrs.txt``. The remaining
columns contain the fitted input moments (``moment_<name>``). If you also pass
``--common-weight``, the table includes the raw common-weight input moments
(``common_weight_moment_<name>``) for each omitted block.

Standard errors and boundary estimates
--------------------------------------

By default, the script sorts SNPs by chromosome and position and forms
200 contiguous blocks of approximately equal SNP count. Use ``--n_blocks``
to change the count, or ``--block-column NAME`` to supply genomic block labels
in your summary statistics. Blocks should span local LD. The number of blocks
must fit the retained SNP count, and each deletion must leave at least two SNPs.

The point estimate uses all retained SNPs. The default ``--jackknife refit``
leaves out each block in turn, refits the weights and constrained covariance
matrix, and recalculates every statistic. ``--jackknife frozen`` is a faster
approximation that holds the full-sample weights fixed while solving the
constrained fit for each deletion. Optional Numba kernels accelerate the same
equations used by the NumPy implementation.

Check the diagnostics for boundary fits and undefined statistics. A correlation
requires two positive fitted variances and an SD ratio requires a positive
denominator variance. An undefined deletion makes the corresponding jackknife
standard error ``NA``. At a zero variance or a perfect correlation, jackknife
standard errors are descriptive and ordinary normal confidence intervals can
be unreliable.

See :ref:`correlate.py <correlate.py>` for all command-line options.
