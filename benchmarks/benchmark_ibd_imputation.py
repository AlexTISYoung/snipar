"""Compare parental imputation on shared simulated families and markers.

The work root contains ``data/kN`` fixtures, ``ibd/kN/METHOD`` segment files
and executable ``current-python`` and ``legacy-python`` wrappers. Run and score
subcommands create ``imputation/kN/PHASE`` outputs and ``imputation/summary``
tables. Scoring uses the finite intersection separately for each family size,
phase input, and parental mask. Paired percentile intervals resample families
with the same weights for each method and hold SNP centering fixed.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor
from contextlib import ExitStack
import csv
import json
import os
from pathlib import Path

from bed_reader import open_bed
import h5py
import numpy as np

from profile_1000g_eur import timed_command
from validate_1000g_eur import _calibration_moments


PHASE_INPUTS = {"true": "chr_22.vcf.gz", "switch": "switch/chr_22.vcf.gz",
                "statistical": "beagle/beagle.vcf.gz"}
METRICS = ("mean_error", "rmse", "slope_pooled", "slope_snp_centered", "variance_ratio_pooled",
           "variance_ratio_snp_centered")


def _methods(k, phase):
    using_phase = phase != "unphased"
    yield "legacy_pair", "baseline-pair", False, True, False
    yield "pair_hard", "current-pair", False, False, False
    if using_phase:
        yield "pair_hard_origin", "current-pair", False, False, True
    if k == 2:
        yield "pair_posterior", "current-pair", True, False, using_phase
        if using_phase:
            yield "phased_hard", f"phased-{phase}", False, False, True
            yield "phased_posterior", f"phased-{phase}", True, False, True
    elif k in (3, 4):
        yield "joint_hard", "current-joint", False, False, using_phase
    elif k == 5:
        yield "joint_fallback_hard", "current-joint", False, False, using_phase
    yield "oracle_hard", "oracle", False, False, using_phase


def _genotypes(root, fixture, k, phase):
    if phase == "unphased":
        stem = fixture / "chr_22"
        return ("--bed", str(stem)) if stem.with_suffix(".bed").is_file() else None
    source = fixture / PHASE_INPUTS[phase]
    if not source.is_file():
        return None
    log = source.with_suffix("").with_suffix(".log")
    if phase == "statistical" and log.is_file() and "End time" not in log.read_text():
        return None
    destination = root / "imputation" / "genotypes" / f"k{k}" / phase
    stem = destination / "genotypes"
    if not stem.with_suffix(".bgen").is_file():
        from snipar.vcf import convert_vcf

        result = convert_vcf(source, destination, chromosomes=[22], output="phased", min_quality={})
        stem = Path(result.bgen)
    return "--bgen", str(stem)


def run(root, sizes, phases, current, legacy, threads, jobs):
    environment = os.environ | dict.fromkeys(("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                                              "VECLIB_MAXIMUM_THREADS"), "1")
    environment["NUMBA_NUM_THREADS"] = str(threads)
    tasks = []
    skipped = []
    for k in sizes:
        fixture = root / "data" / f"k{k}"
        for phase in phases:
            output = root / "imputation" / f"k{k}" / phase
            pending = []
            for method, ibd, soft, baseline, origin in _methods(k, phase):
                stem = root / "ibd" / f"k{k}" / ibd / "chr_22.ibd"
                segments = Path(f"{stem}.segments.gz")
                posterior = Path(f"{stem}.posterior.hdf5")
                if not segments.is_file() or (soft and not posterior.is_file()):
                    skipped.append((k, phase, method, str(posterior if soft else segments)))
                    continue
                if ibd.startswith("phased-") and not (stem.parent / "metrics.json").is_file():
                    skipped.append((k, phase, method, "phased IBD writer still running"))
                    continue
                target = output / f"{method}.hdf5"
                artifacts = (segments, posterior) if soft else (segments,)
                if target.is_file() and all(target.stat().st_mtime_ns >= path.stat().st_mtime_ns for path in artifacts):
                    continue
                if not soft:
                    # Imputation auto-discovers a sidecar beside the segment file.
                    # The alternate stem isolates hard calls on the same segments.
                    stem = root / "imputation" / "segments" / f"k{k}" / ibd / "chr_22.ibd"
                    stem.parent.mkdir(parents=True, exist_ok=True)
                    link = Path(f"{stem}.segments.gz")
                    if not link.exists():
                        link.symlink_to(segments)
                pending.append((method, stem, posterior if soft else None, baseline, origin))
            if not pending:
                continue
            genotype = _genotypes(root, fixture, k, phase)
            if genotype is None:
                skipped.extend((k, phase, method, "genotype missing") for method, *_ in pending)
                continue
            output.mkdir(parents=True, exist_ok=True)
            for method, stem, posterior, baseline, origin in pending:
                tool = legacy if baseline else current
                target = output / method
                command = [str(tool), "-m", "snipar.scripts.impute", *genotype, "--ibd", str(stem),
                           "--pedigree", str(fixture / "pedigree.txt"), "--out", str(target),
                           "--threads", str(threads), "-silent_progress"]
                if not baseline and phase != "unphased":
                    command += ["--parent_origin_map", str(fixture / "map.txt")] if origin else ["--no_parent_origin"]
                if posterior is not None:
                    command += ["--ibd-posterior", str(posterior)]
                tasks.append((k, phase, method, target, command))

    def one(task):
        k, phase, method, target, command = task
        try:
            measurement = timed_command(command, target.with_suffix(".log"), environment)
        except (RuntimeError, OSError) as error:
            return (k, phase, method, str(error))
        measurement["command"] = command
        target.with_suffix(".run.json").write_text(json.dumps(measurement, indent=2) + "\n")
        print(f"k{k} {phase} {method}  {measurement['seconds']:.3f} s  "
              f"{measurement['peak_rss_mib']:.1f} MiB", flush=True)
        return None

    with ThreadPoolExecutor(max_workers=jobs) as pool:
        errors = [result for result in pool.map(one, tasks) if result is not None]
    for label, records in (("missing", skipped), ("failed", errors)):
        for k, phase, method, detail in records:
            print(f"{label} k{k} {phase} {method}  {detail}", flush=True)
    return errors


def _read_columns(dataset, columns):
    order = np.argsort(columns)
    sorted_columns = columns[order]
    if np.all(np.diff(sorted_columns) == 1):
        selected = dataset[:, sorted_columns[0]:sorted_columns[-1] + 1]
    else:
        selected = dataset[:, sorted_columns]
    return selected[:, np.argsort(order)]


class ImputedDosages:
    def __init__(self, file, families, expected_status, sid, allele1, allele2):
        index = {fid: i for i, fid in enumerate(file["families"].asstr()[:])}
        self.rows = np.array([index.get(fid, -1) for fid in families])
        has_row = self.rows >= 0
        np.testing.assert_array_equal(file["parental_status"][:, :2][self.rows[has_row]], expected_status[has_row])
        columns = file["bim_columns"].asstr()[:].tolist()
        metadata = file["bim_values"].asstr()[:]
        index = {marker: i for i, marker in enumerate(metadata[:, columns.index("id")])}
        self.columns = np.array([index.get(marker, -1) for marker in sid])
        present = self.columns >= 0
        self.flip = np.zeros(len(sid), dtype=bool)
        self.dataset = file["imputed_par_gts"]
        if not present.any():
            return
        rows = metadata[self.columns[present]]
        if "allele_ids" in columns:
            first, second = np.array([alleles.split(",") for alleles in rows[:, columns.index("allele_ids")]]).T
        else:
            first, second = rows[:, columns.index("allele1")], rows[:, columns.index("allele2")]
        same = (allele1[present] == first) & (allele2[present] == second)
        flip = (allele1[present] == second) & (allele2[present] == first)
        if not np.all(same | flip):
            raise ValueError(f"Imputed alleles disagree with truth in {file.filename}")
        self.flip[present] = flip

    def read(self, block):
        columns = self.columns[block]
        has_row = self.rows >= 0
        present = columns >= 0
        result = np.full((len(self.rows), len(columns)), np.nan)
        if present.any() and has_row.any():
            result[np.ix_(has_row, present)] = _read_columns(self.dataset, columns[present])[self.rows[has_row]]
            flip = self.flip[block]
            result[:, flip] = 2 - result[:, flip]
        return result


def _family_index(fixture, bed):
    truth_ids = {label: i for i, label in enumerate(bed.iid)}
    observed = set(np.loadtxt(fixture / "chr_22.fam", dtype=str, usecols=1, ndmin=1))
    family_parents = {}
    with (fixture / "pedigree.txt").open() as source:
        columns = next(source).split()
        family_col, father_col, mother_col = [columns.index(key) for key in ("FID", "FATHER_ID", "MOTHER_ID")]
        for line in source:
            row = line.split()
            if not row:
                continue
            family, father, mother = row[family_col], row[father_col], row[mother_col]
            if father in truth_ids and mother in truth_ids and not (father in observed and mother in observed):
                pair = father, mother
                if family_parents.setdefault(family, pair) != pair:
                    raise ValueError(f"Inconsistent parents for family {family}")
    families = sorted(family_parents)
    pairs = [family_parents[family] for family in families]
    if not pairs:
        raise ValueError(f"No observed parental truth for masked families in {fixture}")
    lookup = np.array([[truth_ids[parent] for parent in pair] for pair in pairs])
    status = np.array([[parent in observed for parent in pair] for pair in pairs])
    return families, lookup, status


def _ratio(moment):
    n, sx, sy, sxx, sxy, stt = np.moveaxis(moment[..., :6], -1, 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        variance = sxx - sx**2 / n
        slope = np.where(variance > 0, 1 + (sxy - sx * sy / n) / variance, np.nan)
        truth_variance = stt - (sx + sy)**2 / n
        ratio = np.where(truth_variance > 0, variance / truth_variance, np.nan)
    return slope, ratio


def _estimates(pooled, centered):
    slope, ratio = _ratio(pooled)
    center_slope, center_ratio = _ratio(centered)
    with np.errstate(divide="ignore", invalid="ignore"):
        bias = pooled[..., 6] / pooled[..., 0]
        rmse = np.sqrt(pooled[..., 7] / pooled[..., 0])
    return np.stack((bias, rmse, slope, center_slope, ratio, center_ratio), axis=-1)


def _family_score(paths, fixture):
    with open_bed(fixture / "truth" / "chr_22.bed", count_A1=True) as bed, ExitStack() as stack:
        files = {mode: stack.enter_context(h5py.File(path)) for mode, path in paths.items()}
        families, parents, status = _family_index(fixture, bed)
        masks = {"both_parents_missing": ~status.any(axis=1), "one_parent_missing": status.any(axis=1)}
        sid, allele1, allele2 = bed.sid, bed.allele_1, bed.allele_2
        readers = {mode: ImputedDosages(file, families, status, sid, allele1, allele2) for mode, file in files.items()}
        coverage = {mode: np.zeros(len(families), dtype=np.int64) for mode in paths}
        truth_sizes = np.zeros(len(families), dtype=np.int64)
        pooled = {mode: np.zeros((len(families), 8)) for mode in paths}
        centered = {mode: np.zeros((len(families), 6)) for mode in paths}
        for start in range(0, len(sid), 2048):
            block = slice(start, min(start + 2048, len(sid)))
            father, mother = bed.read(index=(parents.ravel(), block), dtype="float64").reshape(len(families), 2, -1).transpose(1, 0, 2)
            expected = (father + mother) / 2
            expected[status[:, 0]] = mother[status[:, 0]]
            expected[status[:, 1]] = father[status[:, 1]]
            finite = np.isfinite(expected)
            truth_sizes += finite.sum(axis=1)
            common = finite.copy()
            dosages = {}
            for mode, reader in readers.items():
                dosage = reader.read(block)
                calls = np.isfinite(dosage)
                common &= calls
                coverage[mode] += (finite & calls).sum(axis=1)
                dosages[mode] = dosage
            for mode, dosage in dosages.items():
                error = np.where(common, dosage - expected, np.nan)
                pooled[mode][:, :6] += _calibration_moments(dosage, error, within_snp=False)
                pooled[mode][:, 6] += np.nansum(error, axis=1)
                pooled[mode][:, 7] += np.nansum(error**2, axis=1)
                for rows in masks.values():
                    centered[mode][rows] += _calibration_moments(dosage[rows], error[rows], within_snp=True)
    return {mask: (rows, truth_sizes, coverage, pooled, centered) for mask, rows in masks.items() if rows.any()}


def _mask_results(mask, rows, truth_sizes, coverage, pooled, centered, bootstrap_count, seed, *,
                  references=("legacy_pair", "pair_hard_origin", "pair_hard")):
    n_families = int(rows.sum())
    rng = np.random.default_rng(seed)
    weights = rng.multinomial(n_families, np.full(n_families, 1 / n_families), size=bootstrap_count)
    n_truth = int(truth_sizes[rows].sum())
    draws = {}
    results = {}
    for mode in pooled:
        p, c = pooled[mode][rows], centered[mode][rows]
        estimates = _estimates(p.sum(axis=0), c.sum(axis=0))
        draws[mode] = _estimates(weights @ p, weights @ c)
        n_finite = int(coverage[mode][rows].sum())
        results[mode] = dict(mask=mask, families=n_families, n_truth=n_truth, n_finite=n_finite,
                             n_common=int(p[:, 0].sum()), coverage=n_finite / n_truth,
                             estimates=dict(zip(METRICS, estimates.tolist())),
                             ci95=dict(zip(METRICS, np.nanquantile(draws[mode], [.025, .975], axis=0).T.tolist())))
    references = [mode for mode in references if mode in results]
    for mode, entry in results.items():
        entry["paired_differences"] = {}
        for reference in references:
            if mode == reference:
                continue
            change = draws[mode] - draws[reference]
            deltas = {key: entry["estimates"][key] - results[reference]["estimates"][key] for key in METRICS}
            entry["paired_differences"][reference] = dict(estimates=deltas, ci95=dict(
                zip(METRICS, np.nanquantile(change, [.025, .975], axis=0).T.tolist())))
    return results


def score(root, sizes, phases, bootstrap_count, seed, label):
    output = root / "imputation"
    results = []
    for k in sizes:
        fixture = root / "data" / f"k{k}"
        for phase in phases:
            directory = output / f"k{k}" / phase
            paths = {path.stem: path for path in sorted(directory.glob("*.hdf5"))}
            if not paths:
                continue
            print(f"Scoring k{k} {phase}: {', '.join(paths)}", flush=True)
            scored = _family_score(paths, fixture)
            for mask, moments in scored.items():
                for method, entry in _mask_results(mask, *moments, bootstrap_count, seed).items():
                    timing = paths[method].with_suffix(".run.json")
                    measurement = json.loads(timing.read_text()) if timing.exists() else {}
                    results.append(dict(k=k, phase=phase, method=method, **entry, **measurement))
    (output / f"{label}.json").write_text(json.dumps(results, indent=2) + "\n")
    columns = ("k", "phase", "mask", "method", "families", "n_truth", "n_finite", "n_common", "coverage", *METRICS)
    with (output / f"{label}.tsv").open("w") as file:
        writer = csv.DictWriter(file, columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for entry in results:
            writer.writerow(entry | entry["estimates"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("run", "score"))
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--sizes", type=int, nargs="+", default=(2, 3, 4, 5))
    parser.add_argument("--phases", nargs="+", choices=("unphased", *PHASE_INPUTS),
                        default=("unphased", *PHASE_INPUTS))
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--bootstrap", type=int, default=999)
    parser.add_argument("--seed", type=int, default=94821)
    parser.add_argument("--label", default="summary")
    args = parser.parse_args()
    if min(args.threads, args.jobs, args.bootstrap) < 1:
        parser.error("threads, jobs and bootstrap count must be positive")
    root = args.root.resolve()
    if args.command == "run":
        errors = run(root, args.sizes, args.phases, root / "current-python", root / "legacy-python",
                     args.threads, args.jobs)
        if errors:
            raise SystemExit(1)
    else:
        score(root, args.sizes, args.phases, args.bootstrap, args.seed, args.label)


if __name__ == "__main__":
    main()
