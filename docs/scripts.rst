.. _scripts:

====================
Command Line Scripts
====================

.. _snipar_example_data.py:

snipar_example_data.py
----------------------
.. automodule:: snipar.example.snipar_example_data

ibd.py
------
.. _ibd.py:
.. automodule:: snipar.scripts.ibd

impute.py
---------
.. _impute.py:
.. automodule:: snipar.scripts.impute

gwas.py
-------
.. _gwas.py:
.. automodule:: snipar.scripts.gwas

pgs.py
-------
.. _pgs.py:

.. automodule:: snipar.scripts.pgs

For examples and an explanation of the results, see :doc:`pgs`.

pgs.R
-----
.. _pgs.R:

Fits one- and two-generation logistic mixed models for a binary phenotype using
precomputed PGSs and a family random intercept. Both models use the same
individuals. The script requires R and ``lme4``.

**Usage**

    ``pgs.R trait --pgs scores.pgs.txt --phenofile phenotypes.txt --phen case_status``

Input files are whitespace-delimited tables with headers. Each row is matched
by ``FID`` and ``IID``, and ``FID`` identifies the family. Both models exclude
rows with a missing required score, phenotype, or selected covariate.

Args:
    ``OUT`` or ``--out OUT`` : string, required
        Output prefix. The example above writes files beginning with ``trait.``.

    ``--pgs FILE`` : path, required
        Precomputed score file with columns ``FID``, ``IID``, ``proband``,
        ``paternal``, and ``maternal``. A combined ``parental`` score can be
        used with ``--parsum``.

    ``--phenofile FILE`` : path, required
        File containing ``FID``, ``IID``, and one or more phenotype columns.
        The selected phenotype must contain both 0 and 1 in the analysis sample.

    ``--phen NAME`` : string
        Select the phenotype by column name. Use either ``--phen`` or
        ``--phen_index``.

    ``--phen_index N`` : integer, default=1
        Select the phenotype by its one-based position among the phenotype columns.

    ``--covar FILE`` : path
        File containing ``FID``, ``IID``, and numeric covariates for both models.
        By default, both models use all covariate columns.

    ``--covar_cols LIST`` : comma-separated names
        Select columns from ``--covar``, for example ``--covar_cols age,PC1,PC2``.

    ``--parsum`` : flag
        Fit a single parental score in the two-generation model. Uses the
        ``parental`` column if present, or sums ``paternal`` and ``maternal``.

    ``--fit_sib`` : flag
        Include the ``sibling`` score in the two-generation model. Requires
        a column of that name in the PGS file.

    ``--scale_pgs`` : flag
        Divide all included scores by the standard deviation of the proband
        score in the shared analysis sample.

    ``--missing_char VALUE`` : string, default=NA
        Specify an additional missing-value marker. The reader also recognises
        ``NA``, ``NaN``, ``nan``, and ``N/A``.

    ``--maxfun N`` : positive integer, default=200000
        Maximum number of function evaluations for the optimizer.

    ``-h``, ``--help`` : flag
        Print command-line help and exit.

Results:
    ``OUT.effects.txt``
        One tab-separated table with columns ``statistic``, ``est``, and ``SE``.
        Fixed effects use names ``1.TERM`` and ``2.TERM``; ``1.proband`` is
        the population log-odds coefficient and ``2.proband`` is the direct
        log-odds coefficient. The ``direct_population_ratio`` row divides the
        latter by the former, with its delta-method standard error. Standard
        errors account for clustering within families.

    ``OUT.1.vcov.txt``, ``OUT.2.vcov.txt``
        Family-cluster sandwich covariance matrices for the fixed effects of
        each model. Their diagonals give the variances used in the effects table.

    ``OUT.1_2.cross_vcov.txt``
        Sampling covariance between models, with one-generation terms in rows
        and two-generation terms in columns.

    ``OUT.model_info.tsv``
        Fit diagnostics for each model, including case and control counts,
        family count, score scaling and convergence. The script also prints each
        model's between-family random-intercept variance (squared log-odds units)
        and its family-cluster standard error. These appear here as
        ``random_intercept_variance`` and ``random_intercept_variance_se``.
        A boundary fit has an unavailable variance SE (``NA``).
        Check ``ratio_status`` and ``ratio_method`` when interpreting the ratio.

For example input files and guidance on interpreting the ratio, see :doc:`pgs_binary`.

correlate.py
------------
.. _correlate.py:

Estimates the variances and covariances of family-GWAS effects across SNPs,
accounting for their sampling covariance. Population effects are reconstructed
as direct effects plus average non-transmitted coefficients (NTCs).

**Usage**

    ``correlate.py chr_@_young effect --ldscores chr_@ --save_delete``

Args:
    ``sumstats`` : required positional argument
        Summary-statistic file or prefix. A prefix can omit ``.sumstats.gz``.
        Use ``@`` for chromosome numbers, as in ``chr_@_young`` above.

    ``out`` : required positional argument
        Output prefix. In the example, ``effect`` produces ``effect_corrs.txt``
        and ``effect_diagnostics.json``. A prefix can include a directory,
        for example ``results/height`` produces ``results/height_corrs.txt``.

    ``--chr_range``, ``--chr-range`` : integers or ranges
        Select chromosomes, for example ``--chr_range 1 3-5``. A chromosome
        pattern searches for available files on chromosomes 1–22 by default.

    ``--mode`` : default=auto
        Choose ``parental``, ``average``, or ``direct-population``. Auto selects
        the first complete input format across files. Each needs effect
        estimates, standard errors, and a sampling correlation or covariance
        for every pair. See :ref:`correlation_inputs` for the required columns.

    ``--ldscores``, ``--bed``, ``--ld-weight``
        Supply an LDSC score file or prefix with ``--ldscores``, or compute
        scores from chromosome-specific PLINK files with ``--bed``. Either
        route uses inverse LD-score weighting by default. Alternatively pass
        ``--ld-weight none`` alone to give all SNPs equal LD weight.

    ``--ld_wind``, ``--ld-wind`` : float, default=1
        LD window in centimorgans when computing scores from genotypes.

    ``--ld_out``, ``--ld-out`` : prefix
        Save genotype-derived scores as per-chromosome ``.l2.ldscore.gz`` files.

    ``--threads`` : integer
        Set the Numba thread count. Defaults to the available thread count.

    ``--min_maf``, ``--min-maf`` : float, default=0.05
        Retain SNPs with minor allele frequency at or above this threshold.

    ``--corr_filter``, ``--corr-filter`` : float, default=6
        Filter outlying per-SNP sampling correlations in units of standard
        deviations. ``--no-corr-filter`` disables this filter.

    ``--effect-scale``, ``--scale`` : default=per-allele
        Choose ``per-allele`` or ``genotype-sd``. The latter expresses each SNP
        effect per genotype standard deviation and rescales its sampling
        covariance accordingly. ``--genotype-sd`` selects this scale directly.

    ``--duplicate-snps`` : default=error
        Duplicate SNP IDs cause an error. Select ``drop`` to exclude every row
        with a repeated ID.

    ``--n_blocks``, ``--n-blocks`` : integer, default=200
        Number of contiguous genomic jackknife blocks. To supply blocks
        explicitly, use ``--block-column NAME`` with labels in the summary statistics.

    ``--jackknife`` : default=refit
        Refit weights and the constrained covariance in each deletion.
        ``frozen`` holds the full-sample weights fixed for a faster approximation.

    ``--save_delete``, ``--save-delete`` : flag
        Save all leave-one-block-out estimates in one tab-separated file.

    ``--common-weight``, ``--common_weight`` : flag
        Compute raw input moments using the same weight for every effect
        pair at a SNP. Add these estimates to the moments and sampling
        covariance tables, and to the deletion table if requested.

    ``--max_iter`` / ``--max-iter``, ``--tol``, ``--chunk_size`` / ``--chunk-size``
        Set the iteration limit (100), relative convergence tolerance
        (1e-7), and number of SNPs processed per chunk (20000).

    ``-h``, ``--help`` : flag
        Print command-line help and exit.

Results:
    Terminal output
        Lists the input filenames, summarises SNP retention and exclusions,
        and prints the 17 effect statistics with their jackknife standard errors.

    ``OUT_corrs.txt``, ``OUT_diagnostics.json``
        Named estimates and jackknife standard errors, plus input and fit
        diagnostics. The first eight rows give correlations, regression slopes,
        variances and the direct/population variance ratio. Three pairwise
        covariances and six directed standard deviation ratios follow.

    ``OUT_moments.tsv``
        One row per distinct fitted variance or covariance, covering the
        input effects and the direct, average NTC and population effects.
        The columns are ``moment``, ``est`` and ``SE``. With ``--common-weight``,
        the table also includes rows named ``common_weight_moment_<name>``.

    ``OUT_sampling_covariance.tsv``
        A square jackknife covariance matrix covering all estimates in
        ``OUT_corrs.txt`` and ``OUT_moments.tsv``. Each statistic labels one
        row and one column. The diagonal gives its squared standard error.

    ``OUT_delete.txt`` (with ``--save_delete``)
        One row per omitted block. The ``block`` label is followed by the
        17 statistics in ``OUT_corrs.txt`` order and the fitted input moments
        (``moment_<name>``). With ``--common-weight``, the table also includes
        raw input moments (``common_weight_moment_<name>``).

For a worked example, see the :ref:`tutorial <correlations_tutorial>`.
For the full file and row definitions, see :ref:`correlation_outputs`.

simulate.py
------------
.. _simulate.py:

.. automodule:: snipar.scripts.simulate

For a worked example with a neutral control, see :ref:`stabilizing_selection_tutorial`.
