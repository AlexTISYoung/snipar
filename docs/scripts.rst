.. _scripts:

====================
Command Line Scripts
====================

Large datasets
--------------

GWAS, genotype export, scoring, IBD, and imputation size their work blocks
from available memory. The planner also uses Linux cgroup limits and
Slurm resource settings when present. Existing thread, process, and batch
options provide upper bounds. Individual workers reuse sample and
variant alignment across blocks. Imputation shares the detected CPU and
memory allowance across chromosomes and includes flanking markers when
sizing phased blocks.

The family-only GWAS and quantitative PGS mixed models solve their
sibship covariance directly. For GWAS with a general relationship matrix,
SNP regressions use sparse LU, and REML estimation uses the exact determinant.
Supplying previously fitted variance components with ``--vc_list``
skips that estimation step.

IBD uses packed genotypes and traces sibling pairs in blocks. When
needed, it places the packed genotypes and sampled LD data on scratch
storage in the system temporary directory. The free space there can
affect whether a large BGEN run fits. Matrix export and imputation write
genotypes in tiles. Disk space for those outputs still scales with the
number of samples or families times the number of retained variants.

VCF and BGEN input
------------------

``ibd.py``, ``impute.py``, ``gwas.py``, ``gts_matrix.py``, and ``pgs.py``
accept ``--vcf FILE`` in place of the observed BED or BGEN input.
``correlate.py`` accepts it for computing LD scores. Plain VCF, gzip, and
BGZF are supported. Supply one file with one or more autosomes or use
``@`` in the path for chromosome-specific files. ``--chr_range`` limits
which chromosomes are read into the prepared files.

SNIPar reads diploid, biallelic ``GT`` calls and retains both SNPs and
non-symbolic indels. It skips multiallelic, symbolic, and non-autosomal
records. Missing or partially called genotypes produce missing calls.
Records without IDs receive ``CHROM:POS:REF:ALT`` identifiers, with numeric
chromosomes. The REF allele is the first allele in prepared files.

Imputation uses phase when retained genotypes have a consistent phase frame
for each sample and chromosome. It uses the unphased algorithm for unphased
heterozygotes or when phase sets differ within a sample and chromosome.
``simulate.py`` requires phased founder haplotypes and accepts phased VCF
through ``--vcf``. The other commands use allele dosages from the VCF calls.

Phased and unphased BGEN are supported for family genotype analysis. The
dosage for an unphased BGEN genotype is the expected first-allele count
from its genotype probabilities. IBD inference and parental imputation
use hard calls with posterior probability above 0.99, treating uncertain
calls as missing. Use uniformly phased or uniformly unphased diploid,
biallelic BGEN files for those commands.

Filtering large VCFs
~~~~~~~~~~~~~~~~~~~~

SNIPar streams VCF inputs and checks chromosome, variant ID, and site
metadata before decoding genotypes. ``pgs.py`` automatically skips variants
outside the nonzero weights. Use ``--vcf-variants FILE`` for an explicit
variant inclusion list, with one ID per line, or ``--vcf-pass-only`` to
require a PASS record.

``--vcf-min-info-maf`` screens on the allele counts (``INFO/AC`` and
``INFO/AN``) or frequency (``INFO/AF``) reported in the input. Sites with no
usable annotation are still considered. ``--vcf-min-maf`` computes the
frequency from the called genotypes before writing the prepared files.
These opt-in filters use the annotation cohort or the selected VCF samples.
An analysis cohort can have different allele frequencies, and its own
``--min_maf`` rules still apply.

An analysis command writes a temporary seekable BED or phased BGEN, then
removes it when the command finishes. Available scratch space and VCF
decompression can affect performance. Prepare commonly reused input once
with ``vcf_to_genotypes.py`` and pass the reported BED or BGEN prefix to
subsequent commands. Use ``--format auto`` to preserve consistent phase or
``--format phased`` to require it.

.. code-block:: sh

   vcf_to_genotypes.py chr22.vcf.gz prepared --chr_range 22 --vcf-min-info-maf .01 --format auto
   gwas.py phenotype.txt --vcf chr22.vcf.gz --pedigree pedigree.txt --out results

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

vcf_to_genotypes.py
-------------------
.. automodule:: snipar.scripts.vcf_to_genotypes

gwas.py
-------
.. _gwas.py:
.. automodule:: snipar.scripts.gwas

``--impute_unrel`` adds individuals with neither observed nor imputed
parental genotypes to the complete-parent analysis cohort. A singleton
with one available parent needs an imputation for the other parent to
enter this model.

gts_matrix.py
-------------
Writes family genotype arrays with proband and parental dosages and, with
``--fit_sib``, sibling dosages. It accepts the same BED, BGEN, or VCF inputs and
imputation or pedigree inputs as the scoring command.

    ``gts_matrix.py --bed chr_@ --imp chr_@ --out genotypes``

It writes one ``.genoarray.hdf5`` file per chromosome. The files include sample
and variant identifiers, and ``--help`` lists the filtering options.

make_rdr_grms.py
----------------
Computes proband, parental, and symmetrized parent–offspring genomic
relatedness matrices from PLINK BED genotypes and pedigree information.
Use ``--pargts`` to include parental genotypes produced by ``impute.py``.

    ``make_rdr_grms.py chr_1.bed rdr --ped pedigree.txt --hdf5_out``

With ``--hdf5_out``, ``rdr_grms.hdf5`` contains the three matrices, their
variant count and the family and individual IDs in matrix order. The
default output is three sets of GCTA ``.grm.bin``, ``.grm.N.bin`` and
``.grm.id`` files. ``--min_maf`` and ``--max_missing`` filter variants by minor
allele frequency and percent missing calls.

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

    ``--ldscores``, ``--bed``, ``--vcf``, ``--ld-weight``
        Supply an LDSC score file or prefix with ``--ldscores``, or compute
        scores from chromosome-specific PLINK files or a VCF. These inputs
        use inverse LD-score weighting by default. Alternatively pass
        ``--ld-weight none`` alone to give all SNPs equal LD weight.

    ``--vcf-variants``, ``--vcf-pass-only``, ``--vcf-min-info-maf``, ``--vcf-min-maf``
        Filter VCF records before computing LD scores. See the VCF input
        section above for the filter definitions.

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
