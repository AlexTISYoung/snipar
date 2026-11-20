#!/usr/bin/env python
"""Estimate coherent effect moments from SNIPar text summary statistics.

The fitted basis is either (direct, paternal NTC, maternal NTC) or
(direct, average NTC). Population effects are always reconstructed as
direct + average NTC. Run with --help for input and output conventions.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


ALIASES = {
    "SNP": ["SNP", "snp", "rsid", "RSID", "ID"],
    "CHR": ["CHR", "chromosome", "CHROM", "chrom", "chr"],
    "pos": ["pos", "BP", "bp", "position", "POSITION"],
    "freq": ["freq", "f", "FRQ", "EAF", "eaf", "MAF"],
    "direct_Beta": ["direct_Beta", "theta_direct"],
    "direct_SE": ["direct_SE", "se_direct"],
    "avg_NTC_Beta": [
        "avg_NTC_Beta", "avg_parental_Beta", "theta_avgNTC", "theta_avg_NTC"
    ],
    "avg_NTC_SE": ["avg_NTC_SE", "avg_parental_SE", "se_avgNTC", "se_avg_NTC"],
    "paternal_Beta": ["paternal_Beta", "paternal_NTC_Beta", "theta_paternal"],
    "paternal_SE": ["paternal_SE", "paternal_NTC_SE", "se_paternal"],
    "maternal_Beta": ["maternal_Beta", "maternal_NTC_Beta", "theta_maternal"],
    "maternal_SE": ["maternal_SE", "maternal_NTC_SE", "se_maternal"],
    "population_Beta": ["population_Beta", "theta_population"],
    "population_SE": ["population_SE", "se_population"],
    "sib_Beta": ["sib_Beta", "sibling_Beta", "theta_sib"],
}

PAIRS = [
    ("direct", "avg_NTC"), ("direct", "population"),
    ("direct", "paternal"), ("direct", "maternal"), ("paternal", "maternal"),
]
for _left, _right in PAIRS:
    _right_aliases = (
        ["avg_NTC", "avg_parental", "avgNTC"] if _right == "avg_NTC" else [_right]
    )
    for _kind in ("r", "cov"):
        _canonical = "{}_{}_{}".format(_kind, _left, _right)
        _prefixes = ["r", "rg"] if _kind == "r" else ["cov", "sampling_cov"]
        ALIASES[_canonical] = [
            "{}_{}_{}".format(prefix, a, b)
            for prefix in _prefixes
            for right in _right_aliases
            for a, b in ((_left, right), (right, _left))
        ]

MODE_EFFECTS = {
    "parental": ["direct", "paternal", "maternal"],
    "average": ["direct", "avg_NTC"],
    "direct-population": ["direct", "population"],
}


def build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__.split("Run with --help")[0].strip(),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=(
            "Supply the complete sampling covariance through effect SEs and all "
            "pairwise r_EFFECT_EFFECT columns (or cov_EFFECT_EFFECT columns). "
            "Parental mode uses average NTC=(paternal+maternal)/2. Population "
            "columns are redundant and ignored when direct and NTCs are supplied. "
            "Direct-population mode derives NTC=population-direct, transforming "
            "the full sampling covariance. No missing parental covariance is "
            "imputed. Supply LD scores, genotypes to compute them, or --ld-weight none."
        ),
    )
    parser.add_argument("sumstats", help="Text/gzip file or prefix without .sumstats.gz; @ expands to chromosomes.")
    parser.add_argument("out", help="Output filename prefix, e.g. 'effect' writes effect_corrs.txt.")
    parser.add_argument("--chr_range", "--chr-range", nargs="*", default=None,
                        help="Chromosomes as integers/ranges, e.g. 1-22 or 1 3-5; @ defaults to available 1-22.")
    parser.add_argument("--ldscores", help="LDSC text/gzip file or prefix without .l2.ldscore.gz; @ is supported.")
    parser.add_argument("--bed", help="PLINK .bed file or prefix for computing LD scores; @ or # expands to chromosomes.")
    parser.add_argument("--ld_wind", "--ld-wind", type=float, default=1.0,
                        help="Window in cM for computing LD scores from --bed.")
    parser.add_argument("--ld_out", "--ld-out", help="Prefix for computed per-chromosome LDSC files.")
    parser.add_argument("--threads", type=int, default=None,
                        help="Numba thread count for LD computation; defaults to all available.")
    parser.add_argument("--ld-weight", choices=["inverse", "none"], default="inverse",
                        help="Use inverse LD-score weights, or explicitly disable LD weighting.")
    parser.add_argument("--mode", choices=["auto", "parental", "average", "direct-population"], default="auto",
                        help="Auto uses complete parental inputs, then average NTC, then direct/population.")
    parser.add_argument("--min_maf", "--min-maf", type=float, default=0.05,
                        help="Minimum minor allele frequency; boundary values are retained.")
    parser.add_argument("--corr_filter", "--corr-filter", type=float, default=6.0,
                        help="Maximum absolute sampling-correlation z score (legacy default: 6).")
    parser.add_argument("--no-corr-filter", dest="corr_filter", action="store_const", const=None,
                        help="Disable sampling-correlation outlier filtering.")
    parser.add_argument("--effect-scale", "--scale", choices=["per-allele", "genotype-sd"], default="per-allele",
                        help="Genotype-sd multiplies effects by sqrt(2*f*(1-f)) and covariances by 2*f*(1-f).")
    parser.add_argument("--genotype-sd", dest="effect_scale", action="store_const", const="genotype-sd",
                        help="Alias for --effect-scale genotype-sd.")
    parser.add_argument("--n_blocks", "--n-blocks", type=int, default=200,
                        help="Number of contiguous blocks after chromosome/position sorting.")
    parser.add_argument("--block-column", help="Use supplied block labels instead of constructing contiguous blocks.")
    parser.add_argument("--duplicate-snps", choices=["error", "drop"], default="error",
                        help="Reject duplicate summary-statistic SNP IDs, or explicitly exclude every row for each duplicated ID.")
    parser.add_argument("--save_delete", "--save-delete", action="store_true",
                        help="Write all leave-one-block-out estimates to OUT_delete.txt, one row per omitted block.")
    parser.add_argument("--common-weight", "--common_weight", action="store_true",
                        help="Include raw common-weight moments in the moment and sampling covariance tables "
                             "and, with --save_delete, in the per-block table.")
    parser.add_argument("--jackknife", choices=["refit", "frozen"], default="refit",
                        help="Refit optimal weights in each deletion, or freeze full-sample weights.")
    parser.add_argument("--max_iter", "--max-iter", type=int, default=100)
    parser.add_argument("--tol", type=float, default=1e-7)
    parser.add_argument("--chunk_size", "--chunk-size", type=int, default=20000)
    return parser


def parse_chromosomes(tokens):
    if tokens is None or len(tokens) == 0:
        return None
    values = []
    for token in tokens:
        for part in str(token).split(","):
            try:
                if "-" in part:
                    first, last = [int(value) for value in part.split("-")]
                    if first > last:
                        raise ValueError()
                    values.extend(range(first, last + 1))
                else:
                    values.append(int(part))
            except ValueError as exc:
                raise ValueError("Invalid chromosome range: " + part) from exc
    if not values or any(value < 1 for value in values):
        raise ValueError("Chromosomes must be positive integers.")
    return sorted(set(values))


def resolve_files(pattern, suffix, chromosomes=None):
    """Resolve existing chromosome files, retaining the legacy partial-panel behavior."""
    def existing_file(value):
        candidates = [Path(value), Path(value + "." + suffix)]
        if suffix.endswith(".gz"):
            candidates.append(Path(value + "." + suffix[:-3]))
        for candidate in candidates:
            if candidate.is_file():
                return candidate
        return None

    if "@" not in pattern:
        path = existing_file(pattern)
        if path is None:
            raise ValueError("Cannot find input file or prefix: " + pattern)
        inferred = chromosomes[0] if chromosomes and len(chromosomes) == 1 else None
        return [(path, inferred)]
    files = []
    for chromosome in chromosomes if chromosomes is not None else range(1, 23):
        path = existing_file(pattern.replace("@", str(chromosome)))
        if path is not None:
            files.append((path, chromosome))
    if not files:
        raise ValueError("No files matched chromosome pattern: " + pattern)
    return files


def resolve_columns(columns):
    lookup = {str(name).lower(): str(name) for name in columns}
    resolved = {}
    for canonical, aliases in ALIASES.items():
        for alias in aliases:
            if alias.lower() in lookup:
                resolved[canonical] = lookup[alias.lower()]
                break
    return resolved


def missing_mode_columns(columns, mode):
    effects = MODE_EFFECTS[mode]
    missing = [effect + suffix for effect in effects for suffix in ("_Beta", "_SE")
               if effect + suffix not in columns]
    for i, left in enumerate(effects):
        for right in effects[i + 1:]:
            pair = left + "_" + right
            if "r_" + pair not in columns and "cov_" + pair not in columns:
                missing.append("r_{} (or cov_{})".format(pair, pair))
    return missing


def choose_mode(resolved_headers, requested):
    has_sibling = any("sib_Beta" in header for header in resolved_headers)
    if requested == "parental" and has_sibling:
        raise ValueError("SNIPar sibling models include the sibling term in average NTC; use --mode average.")
    modes = [requested] if requested != "auto" else ["parental", "average", "direct-population"]
    if requested == "auto" and has_sibling:
        modes.remove("parental")
    reasons = []
    for mode in modes:
        missing = sorted(set(column for header in resolved_headers
                             for column in missing_mode_columns(header, mode)))
        if not missing:
            return mode
        reasons.append("{}: {}".format(mode, ", ".join(missing)))
    raise ValueError("No complete effect vector and sampling covariance available. " + "; ".join(reasons))


def _table_header(path):
    # Tabs are meaningful field boundaries in TSVs, including missing cells.
    # General whitespace parsing would collapse those cells and shift the row.
    header = pd.read_csv(path, sep="\t", nrows=0, compression="infer")
    if len(header.columns) > 1:
        return header, "\t"
    return pd.read_csv(path, sep=r"\s+", nrows=0, compression="infer"), r"\s+"


def read_sumstats(files, mode="auto", block_column=None, enforce_chromosome=False, duplicate_snps="error"):
    tables = [_table_header(path) for path, _ in files]
    headers = [header for header, _ in tables]
    resolved_headers = [resolve_columns(header.columns) for header in headers]
    selected_mode = choose_mode(resolved_headers, mode)
    frames = []
    chromosome_unavailable = False
    for (path, inferred_chromosome), (header, delimiter), resolved in zip(files, tables, resolved_headers):
        required = ["SNP", "pos", "freq"]
        missing = [name for name in required if name not in resolved]
        if missing:
            raise ValueError("{} missing required metadata: {}".format(path, ", ".join(missing)))
        if block_column and block_column not in header.columns:
            raise ValueError("{} missing block column {}".format(path, block_column))
        usecols = list(dict.fromkeys(resolved.values()))
        if block_column and block_column not in usecols:
            usecols.append(block_column)
        string_columns = {resolved[name]: str for name in ("SNP", "CHR") if name in resolved}
        if block_column:
            string_columns[block_column] = str
        frame = pd.read_csv(path, sep=delimiter, compression="infer", usecols=usecols,
                            dtype=string_columns, low_memory=False)
        blocks = frame[block_column].copy() if block_column else None
        frame = frame.rename(columns={source: canonical for canonical, source in resolved.items()})
        # Normalize covariance sources separately in each file: one chromosome may
        # supply correlations while another supplies covariances explicitly.
        effect_names = MODE_EFFECTS[selected_mode]
        for i, left in enumerate(effect_names):
            for right in effect_names[i + 1:]:
                pair = left + "_" + right
                if "cov_" + pair not in frame:
                    frame["cov_" + pair] = (_numeric(frame["r_" + pair])
                                            * _numeric(frame[left + "_SE"])
                                            * _numeric(frame[right + "_SE"]))
        if "CHR" not in frame:
            if inferred_chromosome is None:
                # Older single-file inputs omit CHR. Preserve their supplied
                # genomic row order instead of interleaving positions across chromosomes.
                chromosome_unavailable = True
                frame["CHR"] = "1"
            else:
                frame["CHR"] = str(inferred_chromosome)
        elif enforce_chromosome and inferred_chromosome is not None:
            actual_chromosomes = _chromosome_numbers(frame["CHR"])
            if np.any(np.isfinite(actual_chromosomes) & (actual_chromosomes != inferred_chromosome)):
                raise ValueError("{} contains CHR values inconsistent with chromosome {} in its @ pattern."
                                 .format(path, inferred_chromosome))
        if blocks is not None:
            frame["_block"] = blocks
        frames.append(frame)
    data = pd.concat(frames, ignore_index=True)
    rows_read = len(data)
    duplicated = data["SNP"].notna() & data["SNP"].duplicated(keep=False)
    if duplicated.any() and duplicate_snps == "error":
        examples = ", ".join(data.loc[duplicated, "SNP"].astype(str).unique()[:5])
        raise ValueError("Summary statistics must have unique SNP IDs; duplicates include " + examples)
    if duplicate_snps == "drop":
        excluded_ids = data.loc[duplicated, "SNP"].astype(str).unique().tolist()
        data = data.loc[~duplicated].copy()
        data.attrs["duplicate_id_qc"] = dict(rows_read=rows_read,
            excluded_rows=int(duplicated.sum()), excluded_ids=excluded_ids)
    elif duplicate_snps != "error":
        raise ValueError("duplicate_snps must be error or drop")
    data.attrs["chromosome_unavailable"] = chromosome_unavailable
    return data, selected_mode


def _numeric(values):
    return pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)


def build_effect_arrays(data, mode):
    effects = MODE_EFFECTS[mode]
    betas = np.column_stack([_numeric(data[effect + "_Beta"]) for effect in effects])
    ses = np.column_stack([_numeric(data[effect + "_SE"]) for effect in effects])
    covariance = np.zeros((len(data), len(effects), len(effects)))
    covariance[:, np.arange(len(effects)), np.arange(len(effects))] = ses ** 2
    for i, left in enumerate(effects):
        for j in range(i + 1, len(effects)):
            pair = left + "_" + effects[j]
            if "cov_" + pair in data:
                values = _numeric(data["cov_" + pair])
            else:
                values = _numeric(data["r_" + pair]) * ses[:, i] * ses[:, j]
            covariance[:, i, j] = covariance[:, j, i] = values
    if mode == "direct-population":
        transform = np.array([[1.0, 0.0], [-1.0, 1.0]])
        betas = betas @ transform.T
        covariance = np.einsum("ab,nbc,dc->nad", transform, covariance, transform)
        effects = ["direct", "avg_NTC"]
    return betas, covariance, ses, effects


def _chromosome_numbers(values):
    labels = values.astype("string").str.replace(r"(?i)^chr", "", regex=True).str.upper()
    labels = labels.replace({"X": "23", "Y": "24", "XY": "25", "M": "26", "MT": "26"})
    return pd.to_numeric(labels, errors="coerce").to_numpy(dtype=float, na_value=np.nan)


def read_ldscores(files):
    frames = []
    for path, _ in files:
        header, delimiter = _table_header(path)
        lookup = {str(name).lower(): str(name) for name in header.columns}
        snp = next((lookup[name] for name in ("snp", "rsid", "id") if name in lookup), None)
        score = next((lookup[name] for name in ("l2", "ldscore", "ld_score", "ldscores") if name in lookup), None)
        if snp is None or score is None:
            raise ValueError("LD-score file {} requires SNP and L2 columns.".format(path))
        frame = pd.read_csv(path, sep=delimiter, usecols=[snp, score], dtype={snp: str},
                            compression="infer", low_memory=False).rename(columns={snp: "SNP", score: "L2"})
        frames.append(frame)
    data = pd.concat(frames, ignore_index=True)
    data = data.loc[data["SNP"].notna()].copy()
    if data["SNP"].duplicated().any():
        raise ValueError("LD-score files must have unique SNP IDs.")
    data["L2"] = _numeric(data["L2"])
    return data.set_index("SNP")["L2"]


def prepare_inputs(args):
    chromosomes = parse_chromosomes(args.chr_range)
    files = resolve_files(args.sumstats, "sumstats.gz", chromosomes)
    for path, _ in files:
        print(f"Reading summary statistics from {path}", flush=True)
    data, mode = read_sumstats(files, args.mode, args.block_column,
                              enforce_chromosome="@" in args.sumstats, duplicate_snps=args.duplicate_snps)
    effects, sampling, ses, names = build_effect_arrays(data, mode)
    freq = _numeric(data["freq"])
    positions = _numeric(data["pos"])
    chrom = _chromosome_numbers(data["CHR"])
    qc = {"rows_read": len(data)}
    duplicate_qc = data.attrs.get("duplicate_id_qc")
    if duplicate_qc is not None:
        qc["rows_read"] = duplicate_qc["rows_read"]
        qc["excluded_duplicate_id_rows"] = duplicate_qc["excluded_rows"]
    keep = (np.isfinite(effects).all(axis=1) & np.isfinite(sampling).all(axis=(1, 2))
            & np.isfinite(ses).all(axis=1) & np.isfinite(freq)
            & np.isfinite(positions) & np.isfinite(chrom) & data["SNP"].notna().to_numpy())
    if args.block_column:
        keep &= data["_block"].notna().to_numpy()
    qc["excluded_incomplete"] = int((~keep).sum())
    valid_metadata = (freq >= 0) & (freq <= 1) & (positions >= 0) & (chrom >= 1)
    valid_metadata &= (positions == np.floor(positions)) & (chrom == np.floor(chrom))
    if np.any(keep & ~valid_metadata):
        raise ValueError("Complete rows contain invalid frequencies, chromosomes, or positions.")
    chromosome_pass = (np.isin(chrom, chromosomes)
                       if chromosomes is not None and not data.attrs["chromosome_unavailable"]
                       else np.ones(len(data), dtype=bool))
    qc["excluded_chr"] = int((keep & ~chromosome_pass).sum())
    keep &= chromosome_pass
    maf_pass = np.minimum(freq, 1 - freq) >= args.min_maf
    qc["excluded_maf"] = int((keep & ~maf_pass).sum())
    keep &= maf_pass
    ldscores = None
    ld_files = []
    bed_metadata = {}
    if args.ldscores:
        ld_files = resolve_files(args.ldscores, "l2.ldscore.gz", chromosomes)
        for path, _ in ld_files:
            print(f"Reading LD scores from {path}", flush=True)
        ldscores = data["SNP"].map(read_ldscores(ld_files)).to_numpy(dtype=float)
    elif args.bed:
        from snipar.coherent_ld import compute_ldscores
        ldscores = np.full(len(data), np.nan)
        ldscores[keep], bed_metadata = compute_ldscores(
            args.bed, data.loc[keep, "SNP"].to_numpy(),
            None if data.attrs["chromosome_unavailable"] else chrom[keep],
            ld_wind=args.ld_wind, ld_out=args.ld_out)
    if ldscores is not None:
        ld_pass = np.isfinite(ldscores) & (ldscores > 0)
        qc["excluded_ld"] = int((keep & ~ld_pass).sum())
        keep &= ld_pass
    else:
        qc["excluded_ld"] = 0
    if not keep.any():
        raise ValueError("No SNPs remain after complete-case, MAF, and LD-score filtering.")
    if np.any((ses[keep] <= 0).any(axis=1)):
        raise ValueError("Retained SNPs must have strictly positive sampling SEs.")
    eigenvalues = np.linalg.eigvalsh(sampling[keep])
    invalid = eigenvalues[:, 0] <= 0
    if invalid.any():
        ids = data.loc[keep, "SNP"].to_numpy()[invalid][:5]
        raise ValueError(
            "{} retained SNPs have non-positive-definite sampling covariance (e.g. {}). "
            "No repair is applied; supply valid full sampling covariance and a nonredundant effect basis."
            .format(int(invalid.sum()), ", ".join(ids))
        )
    qc["excluded_corr"] = 0
    if args.corr_filter is not None:
        indices = np.flatnonzero(keep)
        cov = sampling[indices]
        diagonal = np.diagonal(cov, axis1=1, axis2=2)
        corr_pass = np.ones(len(indices), dtype=bool)
        for i in range(cov.shape[1]):
            for j in range(i + 1, cov.shape[1]):
                correlation = cov[:, i, j] / np.sqrt(diagonal[:, i] * diagonal[:, j])
                sd = correlation.std()
                # Numerical roundoff in a constant correlation must not delete every SNP.
                if sd > 32 * np.finfo(float).eps * max(1.0, np.max(np.abs(correlation))):
                    corr_pass &= np.abs(correlation - correlation.mean()) <= args.corr_filter * sd
        qc["excluded_corr"] = int((~corr_pass).sum())
        keep[indices[~corr_pass]] = False
    indices = np.flatnonzero(keep)
    if not data.attrs["chromosome_unavailable"]:
        indices = indices[np.lexsort((data["SNP"].astype(str).to_numpy()[indices], positions[indices], chrom[indices]))]
    if len(indices) < 3:
        raise ValueError("At least three complete SNPs are required after filtering.")
    effects, sampling = effects[indices], sampling[indices]
    if args.effect_scale == "genotype-sd":
        heterozygosity = 2 * freq[indices] * (1 - freq[indices])
        if np.any(heterozygosity <= 0):
            raise ValueError("Genotype-sd scaling requires allele frequencies strictly between zero and one.")
        effects = effects * np.sqrt(heterozygosity)[:, None]
        sampling = sampling * heterozygosity[:, None, None]
    qc["rows_retained"] = len(indices)
    metadata = {
        "mode": mode, "effect_scale": args.effect_scale, "ld_weight": args.ld_weight,
        "sumstats_files": [str(path) for path, _ in files],
        "ldscore_files": [str(path) for path, _ in ld_files], "qc": qc,
        "min_maf": args.min_maf, "corr_filter": args.corr_filter,
        "block_column": args.block_column,
        "chromosome_unavailable": data.attrs["chromosome_unavailable"],
        "duplicate_snps": args.duplicate_snps,
        "excluded_duplicate_snp_ids": duplicate_qc["excluded_ids"] if duplicate_qc is not None else [],
        "population_definition": "direct + average NTC",
        "average_ntc_definition": "(paternal + maternal) / 2" if mode == "parental" else
                                  "population - direct" if mode == "direct-population" else "input average NTC",
    }
    metadata.update(bed_metadata)
    blocks = data["_block"].to_numpy()[indices] if args.block_column else None
    return effects, sampling, ldscores[indices] if ldscores is not None else None, blocks, names, metadata


def _json_safe(value):
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def with_legacy_metrics(result):
    """Retain legacy output rows, deriving each from the same joint fit/refits."""
    from snipar.coherent import _jackknife_covariance

    names = list(result["metric_names"])

    def derived(values):
        column = {name: values[..., i] for i, name in enumerate(names)}
        direct = column["var_direct"]
        scale = np.maximum.reduce([np.abs(column[name]) for name in names
                                   if name.startswith(("var_", "cov_"))])
        positive = direct > 1e-10 * np.maximum(scale, np.finfo(float).tiny)
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.stack([
                np.where(positive, column["cov_direct_population"] / direct, np.nan),
                np.where(positive, column["cov_direct_avg_NTC"] / direct, np.nan),
                column["sqrt_var_ratio_direct_population"] ** 2,
            ], axis=-1)

    extra_names = ["reg_population_direct", "reg_avg_NTC_direct",
                   "var_ratio_direct_population"]
    estimates = np.concatenate([result["estimates"], derived(result["estimates"])])
    deletions = np.concatenate([result["delete_values"], derived(result["delete_values"])], axis=1)
    all_names = names + extra_names
    legacy_order = ["r_direct_population", "reg_population_direct",
                    "r_direct_avg_NTC", "reg_avg_NTC_direct",
                    "var_direct", "var_population", "var_avg_NTC", "var_ratio_direct_population"]
    ordered_names = legacy_order + [name for name in all_names if name not in legacy_order]
    order = [all_names.index(name) for name in ordered_names]
    updated = dict(result)
    updated["metric_names"] = ordered_names
    updated["estimates"] = estimates[order]
    updated["delete_values"] = deletions[:, order]
    updated["metric_covariance"] = _jackknife_covariance(updated["delete_values"], updated["estimates"])
    updated["se"] = np.sqrt(np.maximum(np.diag(updated["metric_covariance"]), 0))
    updated["diagnostics"] = dict(result["diagnostics"])
    updated["diagnostics"]["undefined_se_metrics"] = [
        name for name, se in zip(ordered_names, updated["se"]) if not np.isfinite(se)]
    return updated


def write_outputs(result, prefix, metadata, save_delete=False):
    from snipar.coherent import _jackknife_covariance, _moment_names

    Path(prefix).parent.mkdir(parents=True, exist_ok=True)
    metric_names = result["metric_names"]
    input_names = result["moment_names"]
    common_names = []
    estimates = [result["estimates"], result["moments"]]
    deletes = [result["delete_values"], result["moment_delete"]]
    if "common_weight_moments" in result:
        common_names = ["common_weight_moment_" + name for name in input_names]
        estimates.append(result["common_weight_moments"])
        deletes.append(result["common_weight_moment_delete"])

    # Input moments can also be reported metrics. Use a single column per
    # statistic when estimating their joint covariance from the same blocks.
    names = pd.Index([*metric_names, *input_names, *common_names])
    keep = ~names.duplicated()
    names = names[keep]
    values = np.concatenate(estimates)[keep]
    block_deletes = np.column_stack(deletes)
    sampling_covariance = _jackknife_covariance(block_deletes[:, keep], values)
    summary = pd.DataFrame({"est": values, "SE": np.sqrt(np.maximum(np.diag(sampling_covariance), 0))},
                           index=names)
    summary.loc[metric_names].to_csv(prefix + "_corrs.txt", sep="\t", index_label="correlation", na_rep="NA")
    expanded_names = _moment_names(result["expanded_names"])
    moment_names = list(dict.fromkeys([*input_names, *expanded_names, *common_names]))
    summary.loc[moment_names].to_csv(prefix + "_moments.tsv", sep="\t", index_label="moment", na_rep="NA")
    pd.DataFrame(sampling_covariance, index=names, columns=names).to_csv(
        prefix + "_sampling_covariance.tsv", sep="\t", index_label="statistic", na_rep="NA")
    diagnostics = dict(result["diagnostics"])
    diagnostics["input"] = metadata
    with open(prefix + "_diagnostics.json", "w", encoding="utf-8") as handle:
        json.dump(_json_safe(diagnostics), handle, indent=2, allow_nan=False)
        handle.write("\n")
    if save_delete:
        blocks = pd.Index(diagnostics["block_labels"], name="block")
        delete_names = [*metric_names, *("moment_" + name for name in input_names), *common_names]
        pd.DataFrame(block_deletes, index=blocks, columns=delete_names).to_csv(
            prefix + "_delete.txt", sep="\t", na_rep="NA")


def print_filter_summary(metadata):
    counts = metadata["qc"]
    total, kept = counts["rows_read"], counts["rows_retained"]
    if kept == total:
        print(f"All {total:,} SNPs passed filtering.", flush=True)
        return

    print(f"SNP filtering: retained {kept:,} of {total:,} input rows.")
    reasons = {
        "excluded_duplicate_id_rows": "duplicate SNP IDs",
        "excluded_incomplete": "missing or non-finite input values",
        "excluded_chr": "outside the selected chromosomes",
        "excluded_maf": f"minor allele frequency below {metadata['min_maf']:g}",
        "excluded_ld": "missing, non-finite or non-positive LD scores",
        "excluded_corr": f"outlying sampling correlations (threshold: {metadata['corr_filter']} SD)",
    }
    print("Rows excluded by reason:")
    for key, reason in reasons.items():
        if counts.get(key, 0):
            print(f"  {counts[key]:>8,}  {reason}", flush=True)


def print_estimates(result):
    # Lead with the three manuscript statistics. Use the output keys to keep
    # estimates and jackknife SEs paired even if the file order changes.
    labels = {
        "r_direct_population": "Correlation: direct and population",
        "r_direct_avg_NTC": "Correlation: direct and average NTC",
        "sqrt_var_ratio_population_direct": "SD ratio: population / direct",
        "reg_population_direct": "Regression slope: population on direct",
        "reg_avg_NTC_direct": "Regression slope: average NTC on direct",
        "var_direct": "Variance: direct",
        "var_population": "Variance: population",
        "var_avg_NTC": "Variance: average NTC",
        "var_ratio_direct_population": "Variance ratio: direct / population",
        "cov_direct_avg_NTC": "Covariance: direct and average NTC",
        "cov_direct_population": "Covariance: direct and population",
        "cov_avg_NTC_population": "Covariance: average NTC and population",
        "sqrt_var_ratio_direct_population": "SD ratio: direct / population",
        "sqrt_var_ratio_avg_NTC_population": "SD ratio: average NTC / population",
        "sqrt_var_ratio_direct_avg_NTC": "SD ratio: direct / average NTC",
        "sqrt_var_ratio_population_avg_NTC": "SD ratio: population / average NTC",
        "sqrt_var_ratio_avg_NTC_direct": "SD ratio: average NTC / direct",
    }
    table = pd.DataFrame({"Estimate": result["estimates"], "SE": result["se"]},
                         index=result["metric_names"])
    order = [key for key in labels if key in table.index]
    order.extend(key for key in table.index if key not in labels)
    table = table.loc[order].rename(index=labels)
    print(f"\nEstimates and jackknife standard errors ({result['diagnostics']['n_blocks']:,} blocks)")
    print(table.to_string(float_format="{:.5g}".format, na_rep="NA"), flush=True)


def main(argv=None):
    parser = build_parser()
    if isinstance(argv, argparse.Namespace):
        # Existing Python callers pass the original parser's Namespace.
        args = parser.parse_args([str(argv.sumstats), str(argv.out)])
        vars(args).update(vars(argv))
    else:
        args = parser.parse_args(argv)
    if args.ld_weight == "inverse" and args.ldscores is None and args.bed is None:
        parser.error("provide --ldscores or --bed, or explicitly choose --ld-weight none")
    if args.ld_weight == "none" and (args.ldscores is not None or args.bed is not None):
        parser.error("--ldscores/--bed and --ld-weight none are mutually exclusive")
    if args.ldscores is not None and args.bed is not None:
        parser.error("provide only one of --ldscores and --bed")
    if not np.isfinite(args.ld_wind) or args.ld_wind <= 0:
        parser.error("--ld_wind must be finite and positive")
    if args.threads is not None and args.threads < 1:
        parser.error("--threads must be positive")
    if not np.isfinite(args.min_maf) or not 0 <= args.min_maf <= 0.5:
        parser.error("--min_maf must be between 0 and 0.5")
    if args.corr_filter is not None and (not np.isfinite(args.corr_filter) or args.corr_filter <= 0):
        parser.error("--corr_filter must be finite and positive")
    if args.n_blocks < 2 or args.max_iter < 1 or args.chunk_size < 1 or not np.isfinite(args.tol) or args.tol <= 0:
        parser.error("require n_blocks >= 2, max_iter >= 1, chunk_size >= 1, and finite tol > 0")
    # Keep --help usable without importing the estimator or SciPy.
    from snipar.coherent import fit_coherent
    try:
        from numba import config as numba_config, set_num_threads
    except ImportError:
        pass  # The covariance estimator also works without Numba.
    else:
        num_threads = numba_config.NUMBA_NUM_THREADS
        if args.threads is not None:
            num_threads = min(args.threads, num_threads)
        set_num_threads(num_threads)
        print("Number of threads: " + str(num_threads), flush=True)
    try:
        effects, sampling, ldscores, blocks, names, metadata = prepare_inputs(args)
        print_filter_summary(metadata)
        print("Fitting {} mode ({} effects).".format(metadata["mode"], args.effect_scale), flush=True)
        if metadata["mode"] in ("parental", "average"):
            print("Population = direct + average NTC; supplied population columns are ignored.", flush=True)
        result = fit_coherent(effects, sampling, ldscores=ldscores, n_blocks=args.n_blocks,
                              blocks=blocks, effect_names=names, max_iter=args.max_iter,
                              tol=args.tol, chunk_size=args.chunk_size, jackknife=args.jackknife,
                              common_weight=args.common_weight)
        result = with_legacy_metrics(result)
        for warning in result["diagnostics"].get("warnings", []):
            print("WARNING: " + str(warning), flush=True)
        write_outputs(result, args.out, metadata, args.save_delete)
        print_estimates(result)
        print("Saved estimates and standard errors to " + args.out + "_corrs.txt", flush=True)
        print("Variance and covariance estimates: " + args.out + "_moments.tsv", flush=True)
        print("Jackknife sampling covariance: " + args.out + "_sampling_covariance.tsv", flush=True)
        print("Filtering and fit diagnostics: " + args.out + "_diagnostics.json", flush=True)
        if args.save_delete:
            print("All leave-one-block-out estimates saved to " + args.out + "_delete.txt", flush=True)
        return result
    except (ValueError, OSError, RuntimeError, np.linalg.LinAlgError) as exc:
        parser.exit(2, "error: {}\n".format(exc))


# Preserve the public parser used by existing programmatic callers.
parser = build_parser()

if __name__ == "__main__":
    main()
