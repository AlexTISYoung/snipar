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


def check_imputation(fixture, output):
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
        index = np.array([truth_families[output_families[row]] for row in rows])
        parent_rows = np.column_stack((2 * index, 2 * index + 1)).ravel()
        father, mother = status[rows].T
        both_missing = ~(father | mother)
        one_missing = father ^ mother
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
            # Report signed error with respect to the same ALT allele for every input format.
            error[:, flip[block]] *= -1
            for summary, group in zip(totals, (one_missing, both_missing)):
                delta = error[group]
                valid = delta[np.isfinite(delta)]
                summary += delta.size, valid.size, valid.sum(dtype=np.float64), np.square(valid).sum(dtype=np.float64)
    return {
        group: dict(n=int(n), finite=int(finite), mean_error=float(bias / finite),
                    rmse=float(np.sqrt(squared / finite)))
        for group, (n, finite, bias, squared) in zip(groups, totals) if finite
    }


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
    parser.add_argument("--pgs", type=Path, help="Calculated .pgs.txt file")
    args = parser.parse_args()
    if args.ibd is None and args.imputation is None and args.pgs is None:
        parser.error("Supply at least one output to validate")
    results = {}
    if args.ibd is not None:
        results["ibd"] = check_ibd(args.fixture, args.ibd)
    if args.imputation is not None:
        results["imputation"] = check_imputation(args.fixture, args.imputation)
    if args.pgs is not None:
        results["pgs"] = check_pgs(args.fixture, args.pgs)
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
