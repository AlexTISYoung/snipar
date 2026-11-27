Scaling workloads
=================

The comparison baseline is ``6ac422c`` on ``development/audit-performance``.
Timings below use Python 3.12 on Apple Silicon and generated BED datasets.
Peak RSS refers to the measured process. The data sizes distinguish genotype
file rows, the analysis cohort, and the number of SNPs processed.

.. list-table::
   :header-rows: 1
   :widths: 52 16 16 16

   * - Workload
     - Baseline wall time
     - Updated wall time
     - Baseline → updated RSS
   * - GWAS, 1M probands, 2M BED rows, 48 SNPs, estimated REML and SNP effects
     - 83.56 s
     - 9.36 s, four threads
     - 1,987 → 1,374 MiB
   * - Family-design reads, 1M probands, 2M BED rows, six blocks of 8 SNPs
     - 121.17 s
     - 3.12 s
     - 1.79 → 1.13 GiB [1]_
   * - Imputation, 1M BED rows, 64 SNPs, four threads
     - 17.00 s
     - 6.58–7.39 s
     - 1,745 → 1,310–1,570 MiB
   * - Imputation, 100k BED rows, 1,024 SNPs, 16 blocks, four threads
     - 17.35 s
     - 2.29–2.88 s
     - 515 → 370–408 MiB
   * - Quantitative PGS, 60k probands, 120k BED rows, 1,500 SNPs, blocks of 64
     - 17.20 s
     - 3.88 s
     - 542 → 416 MiB
   * - Genotype export, 60k probands, 120k BED rows, 1,500 SNPs
     - 21.50 s
     - 21.90 s
     - 2,554 → 500 MiB

.. [1] RSS was measured separately on two blocks of the same width and cohort.

GWAS with one updated worker took 10.61 seconds and 1,232 MiB. Updated
one- and four-worker HDF5 outputs were bit-identical. The largest GWAS effect
estimate difference from baseline was ``2.33e-10`` on this dataset.
The prepared reader had identical genotype checksums and missing counts.
The genotype export datasets were identical. Imputation agreed across the
HDF5 datasets, with tolerance for floating-point reduction in the
genotyping-error diagnostic. A separate join of 400,000 observed and 390,000
imputed variant IDs took 1.01 → 0.63 seconds, with peak RSS of 0.545 → 0.426
GiB and identical variant order and allele flips.

The 30,000-sibling-pair BED IBD benchmark generated identical segment
boundaries, states and LD scores. Warm wall time was 2.82 → 2.36 seconds
with peak RSS of 1,226 → 417 MiB. The command decodes BED and BGEN in blocks
and keeps two-bit storage for tracing.

Reproduction
------------

Use the same environment in both checkouts and keep the data fixed. The
benchmark scripts import SNIPar from ``PYTHONPATH``. Native extensions must
be built in the selected checkout for imputation. The BED generator emits
phenotypes, pedigree and weights.

.. code-block:: sh

   python benchmarks/scale_workloads.py bed /tmp/snipar-million --families 500000 --variants 48
   PYTHONPATH=/path/to/checkout python benchmarks/profile_family_reader.py /tmp/snipar-million --width 8 --threads 4
   PYTHONPATH=/path/to/checkout python -m snipar.scripts.gwas /tmp/snipar-million/phenotype.txt --bed /tmp/snipar-million/chr_1 --pedigree /tmp/snipar-million/pedigree.txt --cpus 4 --threads 1 --out /tmp/gwas
   python benchmarks/scale_workloads.py impute-bed /tmp/impute-wide --families 50000 --variants 1024
   PYTHONPATH=/path/to/checkout python -m snipar.scripts.impute --bed /tmp/impute-wide/chr_1 --pedigree /tmp/impute-wide/pedigree.txt --out /tmp/imputed --threads 4 --chunks 16 -silent_progress

For the previous GWAS pipeline use ``--cpus 1 --batch_size 8``. For the
previous family reader add ``--unprepared``. To reproduce the medium PGS,
IBD and matrix dimensions, generate 30,000 four-person families with 1,500
variants. The PGS benchmark used ``--batch_size 64 --threads 4 --no_am_adj``.
The IBD benchmark used ``--p_error .001 --threads 4 --ld_out`` and seeded
NumPy with 72 before calling its CLI entry point. For the 1M-row imputation
fixture generate 500,000 parent-offspring families with 64 variants and use
the default chunk configuration.

Scaling and resources
---------------------

The planner reserves part of available memory for the system and native
allocations, and uses Linux cgroup and Slurm limits when available.
Schedulers account for per-worker genotype and solver buffers and bound
pending batches. Family-only GWAS inverts covariance exactly using the
family mean and within-family components. For a general relationship matrix,
memory use for exact sparse REML and LU depends on relatedness topology.
IBD scratch space and the size of genotype exports scale with the underlying
data. The largest million-sample runs above used 48 or 64 SNPs. The
million-sample by 100,000-variant imputation budget was checked separately with
synthetic dimensions. A full run at that joint size has yet to be measured.
