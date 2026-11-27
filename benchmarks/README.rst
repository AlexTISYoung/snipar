SNIPar audit and performance
============================

The baseline is ``154573ad189021692418c6a4422e1bdb1fead6d3`` on ``master``.

Reproducing the command profiles
--------------------------------

Use the same Python environment for both checkouts. Compile the extensions in
both source trees with ``python setup.py build_ext --inplace``, or use unpacked
wheels as source trees. Install R, ``lme4`` and ``numDeriv`` for the binary PGS
workload. PLINK is used by the documented workflow tests and is independent of
the profiler.

.. code-block:: sh

   python benchmarks/profile_cli.py --source /path/to/master --output /tmp/snipar-master
   python benchmarks/profile_cli.py --source /path/to/audit --output /tmp/snipar-audit
   python -m pstats /tmp/snipar-audit/gwas.prof

The profiler runs every installed command and two Python example helpers. It
uses bundled tutorial inputs and generates a seeded BED trio dataset for the
RDR command. It performs a warm-up, three timed subprocess runs, and a separate
``cProfile`` or ``Rprof`` run. Each output directory contains the arguments,
timings, exit codes, logs, and profiler files. Run the two series on an idle
machine for a controlled comparison. Compiled code executes inside the caller
recorded by ``cProfile``.

Selected workload results
-------------------------

These are individual Apple Silicon measurements on Python 3.12, using the same
input, command, and threads for each row. Wall-time results are warm medians.
Uncertainty and scaling depend on the cohort, file format, genotype layout,
and model. The cold-profile entries include Python and JIT startup costs.

.. list-table::
   :header-rows: 1
   :widths: 55 15 15 15

   * - Workload
     - Baseline
     - Audit
     - Speedup
   * - IBD and LD export, 2,014 sibling pairs by 1,000 variants, 2 threads
     - 9.76 s
     - 0.62 s
     - 15.8×
   * - Unphased imputation, 2,000 families by 1,000 variants, 10 chunks
     - 1.85 s
     - 0.91 s
     - 2.0×
   * - Forward simulation, 3,000 families, 1,000 loci and five AM generations
     - 15.70 s
     - 2.04 s
     - 7.7×
   * - Sibling-difference GWAS, 4,028 people by 1,000 variants, one worker
     - 5.12 s
     - 2.35 s
     - 2.2×
   * - Genotype matrix export, 1,200 probands by 1,200 variants
     - 12.01 s
     - 1.54 s
     - 7.8×
   * - Quantitative PGS with BED, observed pedigree and the default AM fit
     - 6.97 s
     - 2.04 s
     - 3.4×
   * - Binary PGS, 1,496 complete cases
     - 1.93 s
     - 1.62 s
     - 1.2×
   * - Correlation on 15,000 variants with common sampling covariance
     - 2.03 s
     - 1.51 s
     - 1.3×
   * - Correlation including LD from real BED, 1,000 variants, cold profile
     - 47.62 s
     - 6.86 s
     - 6.9×
   * - RDR GCTA export, 1,500 observed trios by 1,000 variants
     - 1.04 s
     - 0.88 s
     - 1.2×

For the first row, run the tutorial IBD command with ``--threads 2 --p_error
.001``. The imputation workload uses its BED and IBD output with 10 chunks.
The simulation command uses ``--nfam 3000 --n_am 5 --r_par 0.5 --impute
--save_par_gts``. The GWAS row uses the bundled tutorial BED, phenotype and
pedigree files, with ``--sib_diff --vc_list 0.5 0.5 --cpus 1 --threads 1``.
For the genotype matrix row, generate 600 families with two parents and two
children at 1,200 variants, and run ``gts_matrix.py`` with the BED and the
pedigree. The PGS comparison uses the same bundled tutorial BED, pedigree,
weights and phenotype as inputs and includes the paired regression with AM
adjustment.

The branch improves the compiled heterogeneous-covariance correlation kernel.
For 900,000 variants the two-effect kernel took 49 ms on baseline and 35 ms on
the branch. The three-effect kernel took 104 ms and 97 ms. End-to-end timing on
15,000 heterogeneous variants varied under concurrent tests, with medians of
2.57 s on baseline and 2.97 s on the branch and identical estimates. The
tutorial summary-statistic workload timed at 1.68 s and 1.96 s in the full
command benchmark. The data copy command is dominated by the operating
system's native copy at about 55 ms end to end. The tutorial effect
comparison helper took about 166 ms in each run. The BED phenotype
example completed in 0.70 s after the master command failed on an obsolete
import.

Output throughput trades some disk space for CPU time. The genotype matrix
benchmark grew from 1.13 MB to 1.67 MB with gzip level 3 and shuffle. The
simulation segments in the larger benchmark grew from 9.7 MB to 12.0 MB. The
RDR HDF5 path now completes, and the GCTA path writes the variant-count files
alongside the relatedness matrices. Row-blocked Gram products bound temporary
matrix memory. A pre-audit HDF5 run failed while serializing the pedigree, so
it has no successful baseline wall time. The benchmark script also measures
the tutorial's effect comparison helper and the repaired BED phenotype example.

Audit coverage
--------------

* Genotype readers and I/O. Regressions cover sample order, ID zero padding,
  duplicate SNPs, BED chromosome metadata, imputation controls, fixed-width and
  UTF-8 strings, phased and unphased dosages, current-but-corrupt BGEN indexes,
  and native backend fancy indexing. BED reads only the needed proband and
  parent rows when the design has no sibling term.
* IBD, maps, pedigrees, LD and imputation. HMM transitions and emissions are
  reused across sibling pairs, LD pair moments avoid repeated array passes,
  and native imputation reuses IBD segment lookups and thread-local counters.
  Regressions check HMM paths and LD against independent references, shuffled
  BED/BGEN and PC alignment, missing calls, phased chunk outputs, parallel
  counters and chromosome batches.
* GWAS, sparse mixed models and genotype matrices. Tests compare GLS and the
  robust covariance to dense reference calculations, cover REML, singular
  variants, cohort and covariate selection and group alignment, and check
  batched outputs. GWAS reuses factorization and nuisance covariate projections
  where the design permits. Single-worker jobs run in-process. Affected GWAS
  coefficients can change where the previous implementation used the wrong
  projection or covariance.
* Scores and correlation. Tests cover masking, allele and SNP matching, sample
  order across batches, complete-case regression, clustered uncertainty and
  independent covariance and gradient oracles. The common-covariance
  correlation path computes its exact solution and deletion fits from
  sufficient statistics. Binary PGS retains the Laplace model and computes
  family score derivatives analytically. A generated 8,000-person,
  15-covariate run completes in about 10 seconds after the master fit failed
  during post-fit finite differences.
* Simulation and examples. Tests check packed meiosis and streamed output
  against reference calculations. The four source-only research simulation
  scripts were also profiled on small workloads. Plot-enabled research paths
  require optional Matplotlib. The largest structured simulation defaults
  require a separate large-memory timing run.
* Packaging and architecture. The existing platform dependency markers are
  retained. Distribution CI exercises help for every installed Python command.
  The audit locally exercised both BGEN backends on Apple Silicon and ran
  native extension smoke tests with x86-64 Python under Rosetta. Linux, Windows,
  and the full Python wheel matrix are covered by the CI configuration.

Dependencies considered
-----------------------

PySnpTools uses the Rust ``bed-reader`` decoder for its default BED reads. Using
its sample and variant selection avoids a second genotype reader dependency.
Both existing BGEN packages wrap compiled readers, with the installed backend
chosen according to wheel availability. The changes improve access patterns
and add coverage for both backends. See the `PySnpTools documentation
<https://fastlmm.github.io/PySnpTools/>`_, the `bed-reader Rust API
<https://docs.rs/bed-reader/latest/bed_reader/>`_, the `bgen-reader documentation
<https://bgen-reader.readthedocs.io/en/latest/>`_ and the `native BGEN reader
<https://github.com/jeremymcrae/bgen>`_.

The ``lme4`` performance notes identify the quadratic number of post-fit
finite-difference evaluations and list alternative fitters. The new code
retains ``lme4`` fitting and uses analytical family scores and an observed
information calculation with a linear number of score evaluations. The R
suite compares those calculations with independent finite-difference oracles
for interior and boundary fits. See the `lme4 performance guide
<https://cran.r-project.org/web/packages/lme4/vignettes/lmerperf.html>`_.
