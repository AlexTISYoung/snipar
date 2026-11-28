Sibling IBD and parental genotype imputation
===========================================

The benchmark compares the parent of the IBD-posterior change, commit
``9885cce``, with ``339a96d``. Use separate Python environments or
executable wrappers that load each checkout and its compiled imputation
extension. The benchmark scripts themselves run from the current checkout.

Simulation
----------

The input files are the public NYGC chromosome 22 GT-only VCF, the IGSR
3,202-participant pedigree/population table, and the deCODE GRCh38 map.
See `EUR_1000G.rst <EUR_1000G.rst>`_ for the source links. First simulate
120 five-child families. Then derive 2-, 3-, and 4-child panels by retaining
the first children from those same families. Each panel contains 40
families in each parent-availability group. The one-parent group includes
20 observed fathers and 20 observed mothers.

.. code-block:: sh

   root=/path/to/evaluation
   python benchmarks/prepare_1000g_eur.py chr22.vcf.gz pedigree_population.txt chr22.b38.gmap.gz \
       "$root/data/k5" --start 17000000 --end 25000000 --families 120 --children 5 \
       --unrelated 20 --min-family-maf .01 --seed 20261127
   python benchmarks/prepare_ibd_evaluation.py "$root/data/k5" --switch-rate .001 --seed 20261129
   python benchmarks/prepare_ibd_evaluation.py "$root/data/k5" --oracle-out "$root/ibd"

Generate Beagle phase from the visible calls. Beagle 5.5 (27Feb25.75f)
received no pedigree or reference-panel input. The genotype masks and
variant set are shared across the unphased, true-phase, switch-corrupted
and Beagle conditions. The two phase-quality commands verify genotype
parity and compare orientation at successive heterozygotes.

.. code-block:: sh

   for k in 2 3 4 5; do
       mkdir -p "$root/data/k$k/beagle"
       python benchmarks/benchmark_parent_origin.py phase "$root/data/k$k" chr22.b38.gmap.gz \
           "$root/data/k$k/beagle" --java /path/to/java --beagle-jar /path/to/beagle.jar --threads 2
   done
   python benchmarks/prepare_ibd_evaluation.py "$root/data/k5" --phase-quality switch
   python benchmarks/prepare_ibd_evaluation.py "$root/data/k5" --phase-quality beagle

IBD and imputation results
--------------------------

Create executable wrappers named ``$root/legacy-python`` and
``$root/current-python`` for the named commits. The IBD runner evaluates
raw Viterbi paths with a zero minimum segment length and a fixed seed,
genotyping error, and allele-frequency filter. It checks matching markers
and LD scores within each size. The eight-state phased-pair arm uses
allele-level error of 0.0005 and the inverse LD weights by default.
For controlled switches it composes interval probabilities over skipped
markers. For Beagle it assumes a switch probability of 0.002 per retained
heterozygote. The ``--phase-ld-power`` sensitivity changes only the phase
model's weighting, using 1 for the primary analysis, 0.5 for square-root
weights, and 0 for unit weights.

.. code-block:: sh

   for k in 2 3 4 5; do
       python benchmarks/evaluate_sibling_ibd.py "$root/data/k$k" "$root/ibd/k$k" \
           --legacy-python "$root/legacy-python" --current-python "$root/current-python" \
           --threads 4 --verify-joint-fallback
   done
   python benchmarks/benchmark_ibd_imputation.py run "$root" --threads 2 --jobs 2
   python benchmarks/benchmark_ibd_imputation.py score "$root" --bootstrap 999 --label summary

The true-IBD segment files make the oracle dosage comparison possible.
They use the simulated transmissions and the same markers as every other
imputation arm. Imputation methods separate the effects of phase, the
duoHMM parent-origin step, joint or phased IBD hard calls, and posterior
averaging. The eight-state model is used end to end only for two siblings.
DuoHMM uses the observed parent during phased imputation, after the IBD
run and before calculating parental dosages. The benchmark holds the
IBD calls fixed when isolating the effect of this step.

The output has hard-call and marginal-mode accuracy, per-state
confusion and recall, multiclass Brier score and calibration error, and
joint coherence of pair hard calls. Dosage scores use the common finite
intersection separately in each size, phase and mask. The target is the
known masked-parent ALT dosage or the mean of both parents' ALT dosages
if both are absent. The tables report signed error, RMSE, pooled and
within-SNP regression calibration, variance ratios, and call coverage.
The paired percentile intervals resample families and hold SNP centering
fixed. The sampling scope is the one tested interval. The same families
are reused across sizes.

Empirical masked-parent check
-----------------------------

Apply both masks to the same five real 1000 Genomes quads. The preparation
script writes two sets of observed calls and matching parental truth.
It reserves markers whose position is divisible by five as Mendelian
anchors. At an evaluation marker the paternal and maternal sharing bits
must each be supported by agreeing informative flanks within 100 kb.
The unrelated reference samples enter Beagle as study genotypes and the
parental imputation frequency calculation as singleton candidates. A
separate sensitivity analysis uses the pooled unrelated panel for the
unphased HMM's allele frequencies and/or LD weights on the same retained
markers. Results are descriptive for five families.

.. code-block:: sh

   python benchmarks/prepare_empirical_ibd_quads.py chr22.vcf.gz pedigree_population.txt \
       chr22.b38.gmap.gz "$root/empirical"
   for mask in one_parent zero_parents; do
       mkdir -p "$root/empirical/$mask/beagle"
       python benchmarks/benchmark_parent_origin.py phase "$root/empirical/$mask" chr22.b38.gmap.gz \
           "$root/empirical/$mask/beagle" --java /path/to/java --beagle-jar /path/to/beagle.jar --threads 2
   done
   python benchmarks/evaluate_empirical_ibd_quads.py "$root/empirical" "$root/empirical-results" \
       --legacy-python "$root/legacy-python" --current-python "$root/current-python"

DuoHMM before phased IBD
------------------------

After generating the Beagle, IBD, and imputation outputs above, rerun the
phased arm with the current ``ibd.py --phased-ibd`` CLI. The driver runs both
the default duoHMM pre-phasing and ``--no-duohmm`` on the same calls. Both
arms use duoHMM for imputation, so their comparison isolates the effect of
using its calls in the sibling HMM. The driver also holds the duoHMM-informed
IBD files fixed while comparing origin-sidecar reuse to local recomputation.

.. code-block:: sh

   python benchmarks/evaluate_duohmm_prephase.py "$root" "$root/duo-prephase" \
       --python /path/to/current-python --empirical --threads 2 --jobs 2 --bootstrap 999

The output includes the exact commands and process timings, IBD accuracy
and calibration, within-SNP calibration and variance recovery for imputed
dosages, and equality checks for the zero-parent control and origin reuse.
The driver also reports changes in pairs with one or two observed parents.
The five empirical quads use the same held-out Mendelian truth as
the earlier comparison. Re-running the driver in an existing output
directory reuses completed fits.
