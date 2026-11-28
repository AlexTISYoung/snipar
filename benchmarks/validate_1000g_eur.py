"""Compare benchmark IBD, imputation, and proband scores with family truth."""

import argparse
import csv
import gzip
import json
from pathlib import Path

from bed_reader import open_bed
import h5py
import numpy as np


def check_ibd(fixture, output):
    cm, bp = np.loadtxt(fixture / "chr_22.bim", usecols=(2, 3), unpack=True)
    paths = {}
    with (fixture / "truth" / "transmissions.tsv").open() as file:
        for row in csv.DictReader(file, delimiter="\t"):
            paths[row["IID"]] = [
                (int(row[f"{parent}_start"]), np.fromstring(row[f"{parent}_crossover_cM"], sep=","))
                if row[f"{parent}_crossover_cM"] != "." else (int(row[f"{parent}_start"]), np.array([]))
                for parent in ("paternal", "maternal")
            ]
    confusion = np.zeros((3, 3), dtype=np.int64)
    pairs = set()
    with gzip.open(output, "rt") as file:
        for row in csv.DictReader(file, delimiter="\t"):
            left, right = row["ID1"], row["ID2"]
            pairs.add(tuple(sorted((left, right))))
            start = np.searchsorted(bp, int(row["start_coordinate"]), side="left")
            stop = np.searchsorted(bp, int(row["stop_coordinate"]), side="right")
            positions = cm[start:stop]
            states = np.zeros(len(positions), dtype=np.uint8)
            for a, b in zip(paths[left], paths[right]):
                left_hap = (a[0] + np.searchsorted(a[1], positions, side="right")) % 2
                right_hap = (b[0] + np.searchsorted(b[1], positions, side="right")) % 2
                states += left_hap == right_hap
            confusion[:, int(row["IBDType"])] += np.bincount(states, minlength=3)
    total = int(confusion.sum())
    if not total:
        raise ValueError("No inferred IBD segments cover fixture SNPs")
    return dict(sibling_pairs=len(pairs), covered_pair_snps=total,
                state_accuracy=float(np.trace(confusion) / total),
                counts_by_true_then_inferred_state=confusion.tolist())


def _calibration_moments(dosage, error, *, within_snp):
    """Return per-family moments for dosage calibration and variance recovery."""
    called = np.isfinite(error)
    x = np.where(called, dosage, 0.)
    y = np.where(called, -error, 0.)
    if within_snp:
        count = np.maximum(called.sum(axis=0), 1)
        x -= called * (x.sum(axis=0) / count)
        y -= called * (y.sum(axis=0) / count)
    return np.column_stack((called.sum(axis=1), x.sum(axis=1), y.sum(axis=1),
                            np.square(x).sum(axis=1), (x * y).sum(axis=1),
                            np.square(x + y).sum(axis=1)))


def _fit_calibration(moments):
    """Fit the calibration slope and cluster its uncertainty by simulated family."""
    from scipy.stats import t

    moments = moments[moments[:, 0] > 0]
    families = len(moments)
    if families < 2:
        raise ValueError("Calibration requires imputed dosages from at least two families")
    n, sx, sy, sxx, sxy, stt = moments.sum(axis=0)
    xbar = sx / n
    variance_sum = sxx - xbar * sx
    true_variance_sum = stt - (sx + sy) ** 2 / n
    if variance_sum <= 0:
        raise ValueError("Calibration requires variable imputed dosages")
    slope_error = (sxy - xbar * sy) / variance_sum
    intercept = sy / n - slope_error * xbar
    n_i, sx_i, sy_i, sxx_i, sxy_i = moments[:, :5].T
    score = (sxy_i - intercept * sx_i - slope_error * sxx_i
             - xbar * (sy_i - intercept * n_i - slope_error * sx_i))
    df = families - 1
    se = np.sqrt(families / df * np.square(score).sum()) / variance_sum
    slope = 1 + slope_error
    half_width = t.ppf(.975, df) * se
    p = 2 * t.sf(abs(slope_error) / se, df) if se else float(slope_error == 0)
    variance_ratio = float(variance_sum / true_variance_sum) if true_variance_sum > 0 else None
    return dict(families=families, slope=float(slope), intercept=float(intercept), se=float(se),
                ci95=[float(slope - half_width), float(slope + half_width)], p_slope_equals_one=float(p),
                imputed_to_true_variance_ratio=variance_ratio)


def check_imputation(fixture, output, *, paired_dosages=None, calibration=False):
    bed_path = fixture / "truth" / "chr_22.bed"
    families = np.loadtxt(bed_path.with_suffix(".fam"), usecols=0, dtype=str)[::2]
    snps, alt, ref = np.loadtxt(bed_path.with_suffix(".bim"), usecols=(1, 4, 5), dtype=str, unpack=True)
    truth_families = {fid: index for index, fid in enumerate(families)}
    truth_snps = {sid: index for index, sid in enumerate(snps)}
    totals = np.zeros((2, 4), dtype=float)
    groups = ("one_parent_missing", "both_parents_missing")
    with h5py.File(output, "r") as file:
        output_families = file["families"].asstr()[:]
        status = file["parental_status"][:, :2].astype(bool)
        rows = np.flatnonzero(np.isin(output_families, families) & ~status.all(axis=1))
        family_index = np.array([truth_families[output_families[row]] for row in rows])
        parent_rows = np.column_stack((2 * family_index, 2 * family_index + 1)).ravel()
        father, mother = status[rows].T
        both_missing = ~(father | mother)
        one_missing = father ^ mother
        statuses = (one_missing, both_missing)
        if calibration:
            moments = {group: np.zeros((2, int(status.sum()), 6)) for group, status in zip(groups, statuses)}
        columns = file["bim_columns"].asstr()[:].tolist()
        output_snps = file["bim_values"].asstr()[:]
        if "allele_ids" in columns:
            snp_column = columns.index("rsid")
            allele_1, allele_2 = np.array([alleles.split(",")
                                          for alleles in output_snps[:, columns.index("allele_ids")]]).T
        else:
            snp_column = columns.index("id")
            allele_1 = output_snps[:, columns.index("allele1")]
            allele_2 = output_snps[:, columns.index("allele2")]
        variants = np.array([truth_snps[sid] for sid in output_snps[:, snp_column]])
        same = (allele_1 == alt[variants]) & (allele_2 == ref[variants])
        flip = (allele_1 == ref[variants]) & (allele_2 == alt[variants])
        if not np.all(same | flip):
            raise ValueError("Imputation output and truth alleles do not match")
        bed = open_bed(bed_path)
        for start in range(0, len(variants), 2048):
            block = slice(start, start + 2048)
            truth = bed.read(index=(parent_rows, variants[block]), dtype="float32")
            truth[:, flip[block]] = 2 - truth[:, flip[block]]
            father_gt, mother_gt = truth[::2], truth[1::2]
            expected = (father_gt + mother_gt) / 2
            expected[father] = mother_gt[father]
            expected[mother] = father_gt[mother]
            imputed = file["imputed_par_gts"][rows, block]
            error = imputed - expected
            if paired_dosages is not None:
                error[~np.isfinite(paired_dosages[rows, block])] = np.nan
            # Report signed error with respect to the same ALT allele for every input format.
            error[:, flip[block]] *= -1
            if calibration:
                dosage = imputed.astype(np.float64)
                dosage[:, flip[block]] = 2 - dosage[:, flip[block]]
            for index, (group, status) in enumerate(zip(groups, statuses)):
                delta = error[status]
                valid = delta[np.isfinite(delta)]
                totals[index] += delta.size, valid.size, valid.sum(dtype=np.float64), np.square(valid).sum(dtype=np.float64)
                if calibration:
                    moments[group][0] += _calibration_moments(dosage[status], delta, within_snp=False)
                    moments[group][1] += _calibration_moments(dosage[status], delta, within_snp=True)
    results = {
        group: dict(n=int(n), finite=int(finite), mean_error=float(bias / finite),
                    rmse=float(np.sqrt(squared / finite)))
        for group, (n, finite, bias, squared) in zip(groups, totals) if finite
    }
    if calibration:
        for group, result in results.items():
            result["calibration"] = dict(pooled=_fit_calibration(moments[group][0]),
                                         snp_centered=_fit_calibration(moments[group][1]))
    return results


def check_pgs(fixture, output):
    with output.open() as file:
        columns = next(file).split()
        scores = np.loadtxt(file, usecols=(columns.index("IID"), columns.index("proband")),
                            dtype=str, ndmin=2)
    weights = np.loadtxt(fixture / "weights.txt", skiprows=1, dtype=str, ndmin=2)
    bed = open_bed(fixture / "chr_22.bed")
    samples = {iid: index for index, iid in enumerate(bed.iid)}
    variants = {sid: index for index, sid in enumerate(bed.sid)}
    rows = [samples[iid] for iid in scores[:, 0]]
    cols = [variants[sid] for sid in weights[:, 0]]
    genotypes = bed.read(index=(rows, cols), dtype="float64")
    expected = (genotypes - genotypes.mean(axis=0)) @ weights[:, 3].astype(float)
    observed = scores[:, 1].astype(float)
    np.testing.assert_allclose(observed, expected, rtol=2e-6, atol=2e-5)
    return dict(probands=len(rows), weighted_snps=len(cols),
                maximum_absolute_error=float(np.max(np.abs(observed - expected))))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fixture", type=Path, help="Output directory from prepare_1000g_eur.py")
    parser.add_argument("--ibd", type=Path, help="Inferred .segments.gz file")
    parser.add_argument("--imputation", type=Path, help="Imputation .hdf5 file")
    parser.add_argument("--compare-imputation", type=Path, help="Second imputation .hdf5 file for matched-call comparison")
    parser.add_argument("--calibration", action="store_true", help="Fit true on imputed dosage calibration regressions")
    parser.add_argument("--pgs", type=Path, help="Calculated .pgs.txt file")
    args = parser.parse_args()
    if args.ibd is None and args.imputation is None and args.pgs is None:
        parser.error("Supply at least one output to validate")
    if args.imputation is None and (args.compare_imputation is not None or args.calibration):
        parser.error("--compare-imputation and --calibration require --imputation")
    results = {}
    if args.ibd is not None:
        results["ibd"] = check_ibd(args.fixture, args.ibd)
    if args.imputation is not None:
        if args.compare_imputation is None:
            results["imputation"] = check_imputation(args.fixture, args.imputation, calibration=args.calibration)
        else:
            with (h5py.File(args.imputation, "r") as baseline,
                  h5py.File(args.compare_imputation, "r") as comparison):
                for key in ("families", "parental_status", "bim_columns", "bim_values"):
                    np.testing.assert_array_equal(baseline[key][:], comparison[key][:], err_msg=key)
                results["imputation"] = check_imputation(args.fixture, args.imputation,
                                                          paired_dosages=comparison["imputed_par_gts"],
                                                          calibration=args.calibration)
                results["compare_imputation"] = check_imputation(args.fixture, args.compare_imputation,
                                                                  paired_dosages=baseline["imputed_par_gts"],
                                                                  calibration=args.calibration)
    if args.pgs is not None:
        results["pgs"] = check_pgs(args.fixture, args.pgs)
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
