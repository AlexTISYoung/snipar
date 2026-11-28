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

Imputation and ``ibd.py --phased-ibd`` use phase when retained genotypes have
a consistent phase frame for each sample and chromosome. Imputation switches
to the unphased algorithm for unphased heterozygotes or when phase sets
differ within a sample and chromosome. ``--phased-ibd`` requires consistent phase.
``simulate.py`` requires phased founder haplotypes and accepts phased VCF
through ``--vcf``. The other analysis modes use allele dosages from the VCF calls.

Phased and unphased BGEN are supported for family genotype analysis. The
dosage for an unphased BGEN genotype is the expected first-allele count
from its genotype probabilities. IBD inference and parental imputation
use hard calls with posterior probability above 0.99, treating uncertain
calls as missing. Use uniformly phased or uniformly unphased diploid,
biallelic BGEN files for those commands.

Imputation quality filters
~~~~~~~~~~~~~~~~~~~~~~~~~~

The default variant quality cutoff is 0.99, inclusive. For each VCF, SNIPar
selects the first available numeric INFO tag from ``INFO``, ``R2``, ``RSQ``,
``DR2``, and ``IMPINFO`` in that order, ignoring case when choosing the tag.
VCFs without a recognized tag and BGENs without a quality table can be read
without this filter.

Use ``--min-quality FIELD=MIN`` to choose a field and threshold. Giving just
the field, as in ``--min-quality R2``, uses 0.99. Repeat the option to require
every threshold. Explicit names must match the tag or table column exactly.
Missing, nonnumeric, nonfinite, and nonscalar values fail an active filter.
Use ``--no-quality-filter`` to disable automatic filtering.

For BGEN, supply ``--bgen-quality FILE`` with a whitespace-delimited table
containing a variant ID column and one or more numeric quality columns. A
headered table can be compressed with gzip. SNIPar also accepts the
headerless eight-column UK Biobank MFI format. Use
``--bgen-quality-id COLUMN`` to specify a custom ID column. The table must
have unique IDs. SNIPar matches the BGEN variant ID first, then a unique
rsID when the variant ID is absent from the table. A variant without a
matching row is excluded. For separate chromosome tables, put ``@`` in both
the BGEN prefix and the quality filename.

VCF filtering runs before genotype parsing. BGEN filtering uses an index
into the original file, so excluded genotype blocks are not decoded for the
analysis and the genotypes do not need a second copy on disk. For imputation,
apply BGEN quality filters with BGEN as the only variant input.

.. code-block:: sh

   gwas.py phenotype.txt --vcf chr22.vcf.gz --min-quality R2=.95 --pedigree pedigree.txt --out results
   gwas.py phenotype.txt --bgen chr_@ --bgen-quality ukb_mfi_chr@.txt --pedigree pedigree.txt --out results

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
The annotation and genotype frequencies can differ from those in the
analysis cohort. The analysis command's own ``--min_maf`` rules still apply.

When you pass ``--vcf`` directly, SNIPar converts only the samples needed
for the requested analysis. ``gwas.py`` starts with non-missing phenotypes
and applies ``--keep`` if supplied. ``pgs.py`` and ``gts_matrix.py`` start
with the pedigree and also accept ``--keep`` to restrict proband output.
These commands retain observed parents and, for sibling designs, genotyped
siblings needed by the selected probands. IBD inference retains the sibling
pairs and any parent–offspring pairs used to estimate genotype error.
Imputation retains the genotyped relatives and offspring used to estimate
allele frequencies. ``simulate.py`` with ``--nfam`` stages the requested
founders in ordinary simulations. Stabilizing-selection simulations use the
full input breeding pool. ``correlate.py`` uses the full supplied LD reference.

Direct VCF commands use the full input panel for ``--vcf-min-maf`` and
therefore stage all its samples when the threshold is positive. For a
one-time conversion with an explicit sample list, use
``vcf_to_genotypes.py --samples FILE``. Its genotype-based MAF filter uses
those selected samples.

An analysis command writes a temporary seekable BED or phased BGEN, then
removes it when the command finishes. The reader still streams the entire
compressed VCF, so decompression and scratch space can affect performance.
Prepare commonly reused input once
with ``vcf_to_genotypes.py`` and pass the reported BED or BGEN prefix to
subsequent commands. Use ``--format auto`` to preserve consistent phase or
``--format phased`` to require it.

When an analysis selects a small fraction of a BGEN panel, SNIPar can
decode just those samples from supported packed diploid, biallelic blocks.
It still reads and decompresses each retained block. Other encodings use
the installed BGEN reader.

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

``--map`` accepts a SHAPEIT-format map with ``pposition`` (bp) and
``gposition`` (cM) columns. SNIPar interpolates between the map anchors and
extrapolates linearly past either end using the slope of the nearest anchor
interval. Without ``--map``, it uses informative BIM cM values when available
or the bundled GRCh38 deCODE map, with linear extension at its endpoints.

``--posterior`` writes one HDF5 file per chromosome at
``<out>_chr_<N>.ibd.posterior.hdf5``. An ``@`` in ``--out`` is replaced by the
chromosome number. The file contains sibling ``pairs``, marker ``positions``
in bp, ``variant_ids``, ``genetic_map_cm``, and a float32 ``probabilities``
array indexed by pair, retained marker, and IBD state (0, 1, or 2). Positions
are sorted, and multiple markers can share a position. The gzipped hard
segment files are also written. With ``--batches``, each posterior file uses
the same batch number as its segment file.

``--joint`` fits an exact inheritance-vector HMM to each complete set of
three or four genotyped full siblings, using unphased genotype dosages.
The hard segments for these families come from one coherent Viterbi path
without smoothing individual pairs. The optional HDF5 output contains the
projected pair marginals. SNIPar uses the pair HMM and the usual segment
smoothing for the remaining pairs.

``--phased-ibd`` fits a phased sibling-pair HMM to statistically phased BGEN
or VCF genotypes. With a pedigree, or with KING kinship and ``--agesex``,
SNIPar first applies duoHMM to each sibling with exactly one genotyped
parent. The phased sibling HMM uses confident child-origin calls and the
observed parent's transmitted homologs as soft evidence. Families with zero
or two genotyped parents use the phased pair model without this evidence.
``--no-duohmm`` disables the pre-phasing step. ``--phase-switch-rate`` sets the
sibling model's switch probability per retained heterozygote (default 0.002).
Run ``--joint`` separately for the joint unphased sibship model.

DuoHMM uses the genetic map on the input marker grid before the IBD MAF,
missingness, and error filters. When it runs, SNIPar saves the origin and
homolog labels in ``<out>_chr_<N>.ibd.origin.hdf5``, using the same
chromosome and batch names as the segments. The sidecar records directed
parent–child IDs, marker identity, and the phased calls needed to verify
the labels when imputation loads them.

.. code-block:: sh

   ibd.py --vcf chr_@.vcf.gz --pedigree pedigree.txt --map chr22.map \
       --chr_range 22 --phased-ibd --posterior --out phased_@

Every run also writes ``<out>.ibd.relatedness.tsv.gz``, with ``@`` replaced by
``genome``. For each pair, the file records the number of analyzed chromosomes,
the span covered by retained markers, and the posterior expected lengths of
IBD0, IBD1, and IBD2, in cM and Mb. Lengths are integrated by the trapezoidal
rule between retained markers and summed across chromosomes. The sharing
proportion is ``(0.5 * IBD1 + IBD2) / span`` in each unit. A zero span gives a
missing sharing proportion.

For one phased pair, the Python function
``snipar.ibd_phased.infer_phased_pair`` takes phased alleles, a genetic map,
allele frequencies, and explicit phase-switch and allele-error probabilities.
It returns IBD state marginals and each sibling's phase-switch probabilities.
During IBD2, equal switch rates give symmetric attribution of a relative
phase switch. Unequal rates can change that attribution under the model.

The :download:`1000 Genomes calibration script <../benchmarks/calibrate_ibd_posteriors.py>`
holds out sites at positions divisible by five for Mendelian transmission
labels. It estimates allele frequencies and LD in unrelated reference samples
matched to each family's superpopulation. From the source checkout, supply
the expanded NYGC single-chromosome VCF, the pedigree and population table,
and a GRCh38 deCODE map with ``pos``, ``chr``, and ``cM`` columns.

.. code-block:: sh

   python benchmarks/calibrate_ibd_posteriors.py 1kg.chr22.vcf.gz pedigree_population.txt chr22.decode.map.gz calibration.json

impute.py
---------
.. _impute.py:
.. automodule:: snipar.scripts.impute

For phased BGEN and VCF inputs, ``impute.py`` uses a four-state duoHMM-style
model for each observed parent and child in a family with one missing parent.
The model follows transmission along the chromosome while allowing for phase
switches, crossovers, and genotype mismatches. Calls with at least 99% posterior
support inform which child allele came from the missing parent. At other sites,
imputation uses the existing genotype and local haplotype rules. The output
uses the standard imputation HDF5 layout.

Supply ``--parent_origin_map`` with ``pposition`` (bp) and ``gposition`` (cM)
columns to use a genetic map for locally inferred calls. SNIPar extends
that map linearly at its endpoints. Otherwise the model uses informative BIM cM
values when BED and phased BGEN are supplied together, or approximates 1 cM per
Mb. ``--no_parent_origin`` runs the previous phased method on the same inputs.

Use the same phased genotypes in both commands to reuse the saved origin
calls. For a SNIPar ``--ibd`` prefix such as ``phased_@.ibd``,
``impute.py`` automatically looks for ``phased_@.ibd.origin.hdf5``. Use
``--parent-origin FILE`` to give the sidecar explicitly, with ``@`` for
chromosome numbers. SNIPar matches parent–child pairs and markers, and
verifies haplotype order against the current calls. It fits duoHMM locally
for duos absent from the sidecar and uses its existing fallback at uncovered
markers for matched duos. The saved calls supply transmission and homolog
identity to both ordinary and hard-state imputation.

.. code-block:: sh

   impute.py --vcf chr_@.vcf.gz --pedigree pedigree.txt --chr_range 22 \
       --ibd phased_@.ibd --out imputed_@

Use ``--ibd-posterior PATH`` to supply sibling IBD probabilities from
``ibd.py --posterior``. An ``@`` in ``PATH`` stands for the chromosome number.
For a SNIPar ``--ibd`` prefix such as ``chr_@.ibd``, imputation automatically
uses the matching ``chr_@.ibd.posterior.hdf5`` file when present. SNIPar matches
pairs by sibling IID in either order, uses marker IDs to resolve shared
positions when possible, and linearly interpolates between covered positions
in bp.

For two-sibling families using sibling IBD, SNIPar averages the state-specific
parental dosage estimates using the posterior probabilities of compatible
states. It uses the hard-state result where this average is unavailable.
Larger families continue to use hard states. With ``--ibd``, hard states come
from the supplied segments. When supplying only posteriors, SNIPar takes the
most probable state for each covered pair and marker as the hard state.
Pass the segment files to carry coherent ``--joint`` Viterbi calls into
imputation. A combined posterior estimate for larger families requires
their joint inheritance distribution.

``imputed_par_gts`` stores the ordinary estimate. ``robust_imputed_par_gts``
and ``num_observed_parental_alleles`` store the hard-state estimate and its
measured allele count for robust analyses. The separate
``num_expected_observed_parental_alleles`` dataset stores the posterior-weighted
allele count, and ``posterior_num_sibships`` gives the number of successful
two-sibling posterior imputations per marker.

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
Use ``--keep FILE`` with one proband IID per line to limit the output rows.

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

    ``--min-quality``, ``--no-quality-filter``
        Choose a VCF quality field and minimum or disable the automatic
        0.99 cutoff. See the imputation quality filters above.

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
