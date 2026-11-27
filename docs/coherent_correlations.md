# Coherent correlations between effects

The `correlate.py` command fits one positive
semidefinite (PSD) covariance matrix from the direct and NTC estimates and their
full per-SNP sampling covariance. All moments, correlations, and standard
deviation ratios come from that matrix and share the same SNPs and genomic
jackknife deletions. The implementation is in `snipar/coherent.py`; the older
`snipar.correlate` library remains available for existing Python callers.

The revised manuscript by Tan et al., *Family-GWAS reveals effects of environment
and mating on genetic associations* (20 November 2026, Supplementary Note
§§2.1–2.3 and Supplementary Table 4), reports three rows from `_corrs.txt`.
`r_direct_population` and `r_direct_avg_NTC` measure alignment across SNPs.
`sqrt_var_ratio_population_direct` measures relative root-mean-square magnitude,
with population effects in the numerator. A value above one indicates larger
population effects on this scale. The historical `var_ratio_direct_population`
row uses the reciprocal direction and reports a variance ratio.

## Run

To analyse chromosome-specific summary statistics and LD scores, run

> `correlate.py 'path/to/trait_chr@' results/trait --chr_range 1-22 --ldscores 'path/to/chr@.l2.ldscore.gz' --n_blocks 200 --save_delete`

A single text/gz file, with SNP L2 scores in another file:

> `correlate.py trait.sumstats.gz results/trait --ldscores trait.l2.ldscore.gz`

Compute LD scores from chromosome-specific PLINK files:

> `correlate.py 'path/to/trait_chr@' results/trait --bed 'path/to/genotypes_chr@' --ld_wind 1 --ld_out results/ld_chr --threads 4`

Explicitly request no LD redundancy weighting:

> `correlate.py trait.sumstats.gz results/trait --ld-weight none`

The positional arguments follow snipar's `sumstats out` convention. The input
may be an existing file or a prefix with `.sumstats.gz` appended. Supply LD scores
computed by `ibd.py` or LDSC, or use the existing `--bed` route to compute them
from PLINK genotypes. `--bed` accepts a prefix or `.bed` filename, with `@` or
`#` expanded by chromosome. Each genotype file must contain one chromosome.
`--ld_wind` sets the LD window in cM (default 1); `--threads` controls Numba threads.
`--ld_out PREFIX` writes `PREFIX<chromosome>.l2.ldscore.gz`, or replaces `@`/`#`
in the prefix. The filenames use scalar chromosome numbers, correcting the old
writer's bracketed chromosome names. Computed SNP scores are aligned by ID and
chromosome, with unmatched rows removed by the usual LD-score filter.

The existing flags and defaults are retained: `--min_maf 0.05`,
`--corr_filter 6`, and `--n_blocks 200`. Use `--no-corr-filter` to disable
sampling-correlation outlier filtering. `--help` lists all supported options.
The summary-statistic route needs NumPy, SciPy, and pandas without importing
snipar's compiled imputation extensions; the genotype route additionally needs
SNIPar's genotype dependencies.

If Numba is installed, optional `snipar/coherent_kernels.py` accelerates the same Gaussian
GLS equations with fused compiled loops. Otherwise the NumPy implementation is
used. The compiled path uses guarded small-matrix Cholesky calculations with
NumPy fallback for difficult matrices, bounded accumulation blocks, and
compensated objective summation. It does not change fitted weights, the PSD
constraint, stopping tolerance, deletion-specific scaling or jackknife refits.
The backend is recorded in fit diagnostics. Allow the first run time to compile
and cache these kernels before launching many parallel analyses.

## Input and effect definitions

The preferred input vector is **[direct, paternal NTC, maternal NTC]**, using:

```
direct_Beta direct_SE paternal_Beta paternal_SE maternal_Beta maternal_SE
r_direct_paternal r_direct_maternal r_paternal_maternal
```

For each SNP, diagonal sampling covariances are `SE**2` and off-diagonal entries
are `r * SE_i * SE_j`. The full matrix, including paternal–maternal sampling
covariance, is used. Average NTC is `(paternal + maternal)/2`; population is
`direct + average NTC`. Paternal and maternal effects are not assumed equal.

If only the average NTC is available, **[direct, average NTC]** uses:

```
direct_Beta direct_SE avg_NTC_Beta avg_NTC_SE r_direct_avg_NTC
```

`--mode auto` prefers complete parent-specific columns, then average NTC, then
direct/population columns. `--mode parental`, `average`, and `direct-population`
select the input explicitly. Auto mode detects column availability; if an input
contains entirely empty parental columns, explicitly select `--mode average`.
Direct/population input is transformed jointly to
`N=P-D`, including its sampling covariance; population is then derived from the
fitted matrix. When D/NTC input is supplied, redundant population columns do not
enter the fit. Native snipar sibling models should use the reported average NTC;
their average includes a sibling term, so a parental-only reconstruction is invalid.

SNP identifiers, chromosome, position and frequency are used for QC and genomic
ordering. A single legacy summary file without a chromosome column is also
accepted: its existing row order is preserved for the genomic jackknife, so it
must already be genomically ordered. Selecting chromosomes with `--chr_range`
requires a chromosome column or an `@` pattern in the summary filenames.
When computing LD scores for a file without chromosome labels,
SNP IDs are matched without chromosome labels and must be unique across BED files.
The command accepts native snipar names and common aliases, excludes
incomplete effect vectors and invalid MAF/LD rows, and requires positive definite
sampling covariance matrices. It stops on materially invalid or singular sampling
covariances instead of replacing missing correlations with zero. Supply a
nonredundant basis: do not include population as a third independent component of
`[direct, average NTC]`. SNP IDs must be unique. The explicit option
`--duplicate-snps drop` excludes every row belonging to a duplicated summary-statistic
identifier and records the excluded IDs and row count. It never arbitrarily chooses
one copy. Inputs must already use consistent allele coding and phenotype units.

The default scale is per-allele, as in the existing correlate command. Optional
`--scale genotype-sd` multiplies each effect by `sqrt(2f(1-f))` and its sampling
covariance by `2f(1-f)` before fitting. This changes the target moment scale.

## Estimator and what “optimal” means

For SNP j let y_j be the observed effect vector, S_j its known sampling covariance,
and Σ the common signal covariance. Effects follow a zero-mean working model:

```
E[y_j y_j'] = Σ + S_j = Ω_j
m_j = vech(y_j y_j' - S_j)
V_j[(a,b),(c,d)] = Ω_j[a,c] Ω_j[b,d] + Ω_j[a,d] Ω_j[b,c]
```

`vech` uses NumPy upper-triangle order: DD, DN, NN in two dimensions, or
DD, DF, DM, FF, FM, MM for direct/paternal/maternal input. V is the **full**
Gaussian covariance of the products, from
[Isserlis' formula](https://doi.org/10.1093/biomet/12.1-2.134).
In particular, the variance of a squared effect is `2 * Ω[a,a]**2`;
weights proportional merely to `1/Ω[a,a]` are not optimal variance weights.

With w_j = 1/LDscore_j (or 1 when disabled), each feasible GLS step forms:

```
I = sum_j w_j V_j^-1
h = sum_j w_j V_j^-1 m_j
minimize  0.5 * theta' I theta - theta' h
subject to  Σ(theta) >= 0
```

The PSD constraint is solved in the **information metric**, not by simply
clipping eigenvalues of the final unconstrained estimate. Weights are recomputed
from the current covariance. These are constrained Fisher-scoring steps for the
Gaussian composite objective

```
0.5 * sum_j w_j [ logdet(Σ + S_j) + y_j' (Σ + S_j)^-1 y_j ].
```

Backtracking makes the objective nonincreasing. The fit must satisfy a relative
constrained-score step tolerance; failure raises an error. Numerical scaling
handles small GWAS variances. When every SNP has the same sampling covariance,
the solver whitens by that covariance and clips the whitened total-covariance
eigenvalues at one. This gives the exact Gaussian fit, including for each
refitted jackknife deletion. In the general case, computation is chunked and only small information
matrices are accumulated. This seeks a converged constrained score solution;
the heteroskedastic covariance likelihood need not be globally convex.

The product weights are optimal for linear unbiased moment equations under the
common-covariance Gaussian working model; fitted weights and a PSD constraint
do not confer finite-sample unbiasedness or universal optimality. Under independent
Gaussian SNPs and equal redundancy weights, the interior score is the usual
efficient Gaussian covariance score. A scalar inverse LD score is a redundancy
approximation, **not** a full inter-SNP covariance model or proof of global
efficiency with LD. The analysis estimates marginal SNP-effect moments, not an
LDSC regression slope, LD-deconvolved heritability, or stratification correction.

Under non-Gaussian effects, the Gaussian product covariance is a working weight;
the block jackknife allows departures in fourth moments and within-block dependence.
If true covariance systematically varies with SNP precision/annotations/LD,
the fit is a model-dependent weighted projection, not generally an unbiased
arithmetic finite-panel covariance. Sampling S_j must be calibrated: an error in
S_j also changes the subtracted noise term. Its estimation uncertainty is not
included. The target uses uncentered second moments under the zero-mean effect
convention; empirically centering allele-coded effects would change the estimator.

## Outputs and uncertainty

For the expanded covariance C of **[D, N, P]**, exactly:

```
Var(P) = Var(D) + 2 Cov(D,N) + Var(N)
Cov(D,P) = Var(D) + Cov(D,N)
Cov(N,P) = Cov(D,N) + Var(N)
```

The requested correlations are `Cov(D,P)/sqrt(Var(D)Var(P))` and
`Cov(D,N)/sqrt(Var(D)Var(N))`. A name such as
`sqrt_var_ratio_direct_population` means **sqrt[Var(D)/Var(P)]**. All three
pairwise SD ratios and their reciprocals are emitted with explicit names.

The default run writes four files under the chosen output prefix.

| Suffix | Contents |
| --- | --- |
| `_corrs.txt` | All 17 effect statistics, with columns `correlation` (statistic name), `est` and `SE` |
| `_moments.tsv` | Each distinct fitted variance and covariance for the input effects and for D, N and P, with columns `moment`, `est` and `SE` |
| `_sampling_covariance.tsv` | One square jackknife covariance matrix for the 17 effect statistics and any additional input moments. Each estimate labels a row and a column. Diagonal entries are squared SEs and off-diagonal entries are sampling covariances. |
| `_diagnostics.json` | Filtering, convergence, blocks, boundary and undefined-output diagnostics |

The first eight rows retain the following statistics from earlier releases, in order:
`r_direct_population`, `reg_population_direct`, `r_direct_avg_NTC`, `reg_avg_NTC_direct`,
`var_direct`, `var_population`, `var_avg_NTC`, and `var_ratio_direct_population`.
Three pairwise covariances and six directed SD ratios follow.
For example, `var_ratio_direct_population` gives Var(D)/Var(P), and
`sqrt_var_ratio_direct_population` gives its square root.

Compared with the earlier output layout, `v_population_uncorr_direct` and
`v_avg_NTC_uncorr_direct` are no longer emitted. Separate fitted-covariance
files are consolidated into `_moments.tsv`, and the separate sampling-covariance
files into `_sampling_covariance.tsv`. Read statistics by name when updating
code that used the earlier row positions or file suffixes.

Add `--save_delete` to write all leave-one-block-out estimates to `_delete.txt`.
Each row starts with the omitted block label (`block`).
The next 17 columns give the statistics from that refit in the same order
as the rows of `_corrs.txt`. Columns named `moment_<name>` give the fitted
input-basis moments, replacing the separate `_moment_delete.tsv` file.
Tab-separated files use `NA` for undefined values.

Add `--common-weight` (also spelled `--common_weight`) to include raw
common-weight input moments in `_moments.tsv` and on both axes of
`_sampling_covariance.tsv`. The entries use names such as
`common_weight_moment_var_direct`. With `--save_delete`, the same names
identify the extra columns in `_delete.txt`.

The common-weight diagnostic implements the distinct finite-panel target
`sum_j w_j (y_j y_j' - S_j) / sum_j w_j`
with the **same scalar weights** for every entry. Its block deletions renormalize
those weights on the retained panel. It is unbiased for that weighted finite-panel
second moment conditional on calibrated S and fixed weights/selection, but can
be indefinite because of sampling variation. It is exported without projection
and flagged if non-PSD. Compare these raw moments to the optimal model-based fit
before interpreting precision-weighting changes biologically. Algebraic closure
of the primary fit alone does not prove equality to this finite-panel target.

The default `--jackknife refit` refits weights and the PSD constraint in each
genomic deletion, then recomputes every derived output. With B blocks, the
reported joint sampling covariance is

```
(B-1)/B * sum_b (t_minus_b - mean(t_minus)) (t_minus_b - mean(t_minus))'.
```

The full-sample fit is the point estimate; delete replicates are used for SEs,
not to replace it with a potentially incoherent bias correction. By default,
blocks are contiguous, approximately equal SNP counts after genomic sorting;
choose enough blocks to cover the genome while making them wider than local LD.
An explicit block column can supply prespecified genomic blocks. Highly unequal
blocks produce a diagnostic because the ordinary equal-block formula can then
be poorly calibrated.

`--jackknife frozen` is an explicit faster approximation: product weights from
the full fit are retained and block information totals are subtracted. It still
resolves the PSD constraint per deletion but does not refit the weights.
Full refitting is recommended for final inference.

Correlations require both variances to be positive; ratios require a positive
denominator. If a required full fit or deletion is undefined, the corresponding
SE is NA, with no dropping of invalid deletions. Full-basis and deletion-fit
boundaries are recorded in the fit diagnostics.
At a zero variance or perfect-correlation boundary, ordinary jackknife/Wald
inference is nonregular: the emitted SE is descriptive, not a calibrated normal
confidence interval. If the entire covariance is zero in every deletion, even
moment SEs can be zero; this is a constrained boundary result and does not imply
precise variance estimation. See also
[Shaw and Geyer (1997)](https://doi.org/10.1093/biomet/84.1.95) on constrained
covariance estimation.

## Python API

```python
from snipar.coherent import fit_coherent

# effects: n x 2 [D,N], or n x 3 [D,paternal NTC,maternal NTC]
# sampling_covariance: n x p x p, in exactly the same basis and units
result = fit_coherent(effects, sampling_covariance, ldscores=ld,
                      n_blocks=200, jackknife="refit")
covariance = result["covariance"]
moments, moment_se = result["moments"], result["moment_se"]
estimates = dict(zip(result["metric_names"], result["estimates"]))
standard_errors = dict(zip(result["metric_names"], result["se"]))
```

Rows passed directly to the API must already be ordered genomically, unless
`blocks=` supplies prespecified block membership. It accepts complete finite
vectors; file filtering is the CLI's responsibility.

Python callers can pass `common_weight=False` to skip the common-weight
diagnostic. The Python API includes it by default.
`fit_coherent` returns the core effect moments, correlations and SD ratios.
The command additionally derives the two regression slopes and the
direct/population variance ratio, so the API's `metric_names` are not the
same as the 17 rows in the command's `_corrs.txt` file.
