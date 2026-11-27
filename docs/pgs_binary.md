# Binary-outcome PGS analyses

`snipar/scripts/pgs.R` fits one- and two-generation logistic mixed models to
precomputed proband and parental polygenic scores. Both models use a family
random intercept and the same complete-case sample. The direct to population
ratio and its sampling standard error are calculated automatically.

## Run the script

The script requires R and the R package `lme4`:

> `install.packages("lme4")`

To fit both models with covariates, run the installed command

> `pgs.R results/trait --pgs scores.pgs.txt --phenofile phenotypes.txt --phen case_status --covar covariates.txt --covar_cols age,sex,PC1,PC2`

The first argument is the output prefix; `--out PREFIX` is also accepted. No
additional ratio flag is needed.
Use `pgs.R --help` to list the options. From a source checkout, you can also
run the script with `Rscript snipar/scripts/pgs.R` and the same arguments.

## Input files

Inputs are whitespace-delimited text files with headers. All files require
`FID` and `IID`. The pair identifies each observation; IDs are read as strings,
so leading zeros are retained. Each pair must occur at most once in a file.

The PGS file uses the score column names produced by `pgs.py`:

```text
FID IID proband paternal maternal
001 00101 0.42 0.31 0.12
001 00102 -0.08 0.31 0.12
002 00201 0.17 -0.10 0.38
```

Additional columns such as `FATHER_ID` and `MOTHER_ID` are ignored. `proband`,
`paternal`, and `maternal` must be numeric. With `--parsum`, the two-generation
model uses a single `parental` score: an existing `parental` column is accepted,
or the script sums `paternal` and `maternal`. With `--fit_sib`, the PGS file must
also contain a numeric `sibling` column, which is included in the
two-generation model.

The phenotype file contains IDs and one or more phenotype columns:

```text
FID IID case_status
001 00101 1
001 00102 0
002 00201 1
```

Select the phenotype with `--phen NAME` or `--phen_index N`, where the index is
one-based among the columns other than `FID` and `IID`. The default is the
first phenotype. Observed values must be 0 or 1, and both values must remain in
the analysis sample.

The optional covariate file contains IDs followed by numeric covariates:

```text
FID IID age sex PC1 PC2
001 00101 30 0 -0.02 0.11
001 00102 32 1 0.01 0.09
002 00201 28 0 0.03 -0.04
```

By default, every column other than `FID` and `IID` is included. Use
`--covar_cols age,sex,PC1` to select columns explicitly. Encode categorical
covariates as numeric indicator columns before running the script.

The files are matched by both IDs. Rows missing any required phenotype, score,
or covariate are excluded from both models. Missing parental scores therefore
also reduce the one-generation sample used in this paired analysis. `NA` is
the default missing-value marker; `--missing_char VALUE` adds a custom marker.

## Models and scaling

The one-generation model includes the proband score and selected covariates.
The two-generation model adds paternal and maternal scores, or their sum with
`--parsum`, and optionally the sibling score. Both include an intercept and a
random intercept grouped by `FID`. Families are assumed independent.

Coefficients are reported on the log-odds scale in the input score units.
`--scale_pgs` divides all included score columns by the proband score's standard
deviation in the analysis sample. This uses a common score scale for all
generations. Covariates are used in their supplied units. The binary phenotype
is never standardized.

Models use the Laplace approximation (`nAGQ=1`) and `lme4`'s `bobyqa` optimizer.
`--maxfun N` sets the evaluation limit, which defaults to 200000. A valid fit
requires at least three families and at least one family with multiple
participants in the complete-case sample, together with identifiable fixed
effects; the script reports convergence and likelihood diagnostics.

## Output files

With output prefix `results/trait`, the script writes:

| Suffix | Contents |
| --- | --- |
| `.effects.txt` | All fixed effects and the direct/population ratio, with columns `statistic`, `est`, and `SE`. |
| `.1.vcov.txt`, `.2.vcov.txt` | Labeled sampling covariance matrices for all fixed effects. |
| `.1_2.cross_vcov.txt` | Sampling covariance between models, with one-generation rows and two-generation columns. |
| `.model_info.tsv` | Family random-intercept variance and SE, random-intercept SD, boundary indicator, sample counts, convergence and likelihood diagnostics, plus `ratio_status` and `ratio_method`. |

The tab-separated `.effects.txt` table names coefficients `1.TERM` and
`2.TERM`: `1.proband` is the population effect and `2.proband` is the direct
effect. `direct_population_ratio` gives their ratio, with its estimate in
`est` and standard error in `SE`. Intercepts and covariates follow the same
model-prefixed naming convention.

Counts in `.model_info.tsv` refer to the shared complete-case sample.
The script prints each model's estimated family random-intercept variance and
its SE. The `random_intercept_variance` column measures how much families'
baseline log odds differ after accounting for the fitted scores and covariates.
Its units are squared log odds. The family-cluster delta-method SE is stored in
`random_intercept_variance_se`. It is `NA` at the boundary (estimated SD at
most `1e-4`), where the usual SE approximation breaks down.

The `ratio_method` label is `Laplace_family_cluster_sandwich_delta`. To form
the two-by-two covariance matrix of the population and direct coefficients,
take the `proband` row and `proband` column entry from each file. The
one- and two-generation covariance files give the corresponding diagonal
entries and the cross-model file gives the shared off-diagonal entry.

Tables have headers. Covariance matrices label their rows and columns. The
primary effect SEs, marginal covariance matrices, cross-model covariance, and
ratio SE all use the same family-cluster sandwich calculation, with a
`G/(G-1)` correction for `G` families.

## Ratio and sampling uncertainty

The ratio is the proband coefficient from the two-generation model divided by
the proband coefficient from the one-generation model. It uses the same score
weights and sample in both fits. It is a ratio of log-odds coefficients, rather
than a ratio of odds ratios. Logistic noncollapsibility means it should not be
interpreted as the same decomposition obtained from linear models.

The joint sampling covariance comes from stacked family contributions to the
Laplace likelihood. It includes the random-intercept SD as an estimated
nuisance parameter for interior fits; at the zero-SD boundary it uses the
fixed-effect block. A delta-method calculation then accounts for covariance
between the direct and population coefficient estimates.

A zero population coefficient gives an undefined ratio and SE. A population
coefficient within 1.96 sandwich SEs of zero is flagged as weak, because the
delta approximation can be unreliable near zero. Consult `ratio_status`
and the other diagnostics in `.model_info.tsv` when interpreting results.

## Regression checks

The regression checks additionally require the R package `numDeriv`. From the
source checkout, run:

> `Rscript snipar/tests/test_pgs_binary.R`

The checks compare likelihood derivatives and sampling influences with
independent numerical calculations, exercise the zero random-intercept
variance limit, and run the command-line interface with temporary input files.
They check matched samples, identifier handling, covariance outputs, parental
sums, and sibling effects.
