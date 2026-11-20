"""Input and command-line integration checks for coherent SNIPar estimation."""

import os
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest


from snipar.scripts import correlate as cli


_PATH = Path(cli.__file__).resolve()
_ROOT = _PATH.parents[2]


def _run_cli(arguments, cwd):
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(_ROOT)
    return subprocess.run([sys.executable, str(_PATH), *arguments],
                          cwd=cwd, env=environment, text=True,
                          capture_output=True, check=False, timeout=120)


def _frame(n=24):
    rng = np.random.default_rng(7042)
    effects = rng.multivariate_normal([0, 0], [[.4, .07], [.07, .2]], size=n)
    return pd.DataFrame({
        "chromosome": np.repeat([1, 2], [n // 2, n - n // 2]),
        "SNP": ["rs{}".format(i) for i in range(n)],
        "pos": np.arange(n) + 1, "freq": .2,
        "direct_Beta": effects[:, 0], "direct_SE": .1,
        "avg_NTC_Beta": effects[:, 1], "avg_NTC_SE": .15,
        "r_direct_avg_NTC": -.3,
        # These redundant columns deliberately contain missing cells.
        "population_Beta": np.nan, "population_SE": np.nan,
    })


def _write(frame, path, separator="\t"):
    frame.to_csv(path, sep=separator, index=False,
                 na_rep="NA" if separator == " " else "")


def _load(tmp_path, name="input", options=()):
    args = cli.build_parser().parse_args([
        str(tmp_path / name), str(tmp_path / "out"), "--ld-weight", "none", *options,
    ])
    return cli.prepare_inputs(args)


def test_native_chromosome_files_ld_qc_and_end_to_end_outputs(tmp_path):
    frame = _frame(120)
    block_labels = ["chr1_A", "chr1_B", "chr2_A", "chr2_B"]
    frame["region"] = np.repeat(block_labels, 30)
    frame.loc[0, "direct_SE"] = np.nan
    frame.loc[1, "freq"] = .01
    ld = pd.DataFrame({"SNP": frame.SNP, "L2": 1 + np.arange(len(frame)) / 120})
    ld.loc[3, "L2"] = 0
    for chromosome in (1, 2):
        selection = frame.chromosome == chromosome
        _write(frame.loc[selection].sample(frac=1, random_state=7),
               tmp_path / "chr{}.sumstats.gz".format(chromosome), " ")
        ld_selection = selection & (frame.index != 2)
        _write(ld.loc[ld_selection].sample(frac=1, random_state=8),
               tmp_path / "ld{}.l2.ldscore.gz".format(chromosome), " ")
    arguments = [str(tmp_path / "chr@"), str(tmp_path / "fit"),
                 "--chr_range", "1-2", "--ldscores", str(tmp_path / "ld@"),
                 "--block-column", "region", "--save_delete"]
    prepared = cli.prepare_inputs(cli.build_parser().parse_args(arguments))
    expected_rows = np.arange(4, len(frame))
    np.testing.assert_allclose(prepared[0], frame.loc[expected_rows, ["direct_Beta", "avg_NTC_Beta"]])
    np.testing.assert_allclose(prepared[2], ld.loc[expected_rows, "L2"])
    assert prepared[5]["qc"] == {
        "rows_read": 120, "excluded_incomplete": 1, "excluded_chr": 0,
        "excluded_maf": 1, "excluded_ld": 2, "excluded_corr": 0, "rows_retained": 116,
    }
    completed = _run_cli(arguments, tmp_path)
    assert completed.returncode == 0, completed.stdout + completed.stderr
    report = completed.stdout
    for chromosome in (1, 2):
        assert f"Reading summary statistics from {tmp_path}/chr{chromosome}.sumstats.gz" in report
        assert f"Reading LD scores from {tmp_path}/ld{chromosome}.l2.ldscore.gz" in report
    assert "retained 116 of 120 input rows" in report
    for count, reason in [(1, "missing or non-finite input values"),
                          (1, "minor allele frequency below 0.05"),
                          (2, "missing, non-finite or non-positive LD scores")]:
        assert f"{count}  {reason}" in report
    output = pd.read_csv(tmp_path / "fit_corrs.txt", sep="\t").set_index("correlation")
    assert list(output.columns) == ["est", "SE"]
    assert not output.index.str.endswith("_uncorr_direct").any()
    assert "_uncorr_direct" not in report
    assert "uncorrelated with direct" not in report.lower()
    for label, key in [("Correlation: direct and population", "r_direct_population"),
                       ("SD ratio: population / direct", "sqrt_var_ratio_population_direct"),
                       ("SD ratio: direct / population", "sqrt_var_ratio_direct_population")]:
        row = next(line for line in report.splitlines() if line.startswith(label))
        np.testing.assert_allclose(np.array(row.split()[-2:], dtype=float),
                                   output.loc[key, ["est", "SE"]], rtol=5e-5)
    for metric in ("r_direct_population", "r_direct_avg_NTC",
                   "sqrt_var_ratio_direct_population", "sqrt_var_ratio_direct_avg_NTC",
                   "sqrt_var_ratio_avg_NTC_population"):
        assert np.isfinite(output.loc[metric]).all()
    moments = pd.read_csv(tmp_path / "fit_moments.tsv", sep="\t").set_index("moment")
    assert list(moments.columns) == ["est", "SE"]
    assert moments.index.is_unique
    assert set(moments.index) == {
        "var_direct", "var_avg_NTC", "var_population",
        "cov_direct_avg_NTC", "cov_direct_population", "cov_avg_NTC_population",
    }
    np.testing.assert_allclose(moments, output.loc[moments.index])
    np.testing.assert_allclose(moments.est["var_population"], moments.est["var_direct"]
                               + moments.est["var_avg_NTC"] + 2 * moments.est["cov_direct_avg_NTC"])
    legacy_names = [
        "r_direct_population", "reg_population_direct", "r_direct_avg_NTC", "reg_avg_NTC_direct",
        "var_direct", "var_population", "var_avg_NTC", "var_ratio_direct_population",
    ]
    assert list(output.index[:len(legacy_names)]) == legacy_names
    deletes = pd.read_csv(tmp_path / "fit_delete.txt", sep="\t")
    assert list(deletes.columns) == ["block", *output.index, "moment_var_direct",
                                    "moment_cov_direct_avg_NTC", "moment_var_avg_NTC"]
    assert list(deletes["block"]) == block_labels
    metric_covariance = pd.read_csv(tmp_path / "fit_sampling_covariance.tsv", sep="\t", index_col="statistic")
    assert len(output.index) == 17
    assert list(metric_covariance.index) == list(metric_covariance.columns) == list(output.index)
    np.testing.assert_allclose(np.diag(metric_covariance), output.SE ** 2)
    # Independently derive the compatibility outputs from the common joint
    # moments, including each deletion before calculating jackknife errors.
    def legacy_values(values):
        return {
            "reg_population_direct": values["cov_direct_population"] / values["var_direct"],
            "reg_avg_NTC_direct": values["cov_direct_avg_NTC"] / values["var_direct"],
            "var_ratio_direct_population": values["var_direct"] / values["var_population"],
        }
    for name, expected in legacy_values(output.est).items():
        np.testing.assert_allclose(output.loc[name, "est"], expected)
    for name, expected in legacy_values(deletes).items():
        np.testing.assert_allclose(deletes[name], expected)
        expected_se = np.sqrt((len(expected) - 1) * np.mean((expected - expected.mean()) ** 2))
        np.testing.assert_allclose(output.loc[name, "SE"], expected_se)
    assert {p.name for p in tmp_path.glob("fit_*")} == {
        "fit_corrs.txt", "fit_moments.tsv", "fit_sampling_covariance.tsv", "fit_diagnostics.json", "fit_delete.txt",
    }
    diagnostics = json.loads((tmp_path / "fit_diagnostics.json").read_text())
    assert "common_weight_min_eigenvalue" not in diagnostics
    assert "common_weight_is_psd" not in diagnostics
    assert diagnostics["input"]["mode"] == "average"
    assert diagnostics["n_blocks"] == 4
    assert diagnostics["block_labels"] == block_labels
    assert diagnostics["converged"]


def test_direct_population_transforms_full_sampling_covariance(tmp_path):
    frame = _frame()
    _write(frame, tmp_path / "average.sumstats.gz")
    average = _load(tmp_path, "average")
    population = frame.drop(columns=["avg_NTC_Beta", "avg_NTC_SE", "r_direct_avg_NTC"]).copy()
    population["population_Beta"] = frame.direct_Beta + frame.avg_NTC_Beta
    cov_dn = frame.r_direct_avg_NTC * frame.direct_SE * frame.avg_NTC_SE
    population["population_SE"] = np.sqrt(frame.direct_SE**2 + frame.avg_NTC_SE**2 + 2*cov_dn)
    population["r_direct_population"] = (frame.direct_SE**2 + cov_dn) / (frame.direct_SE * population.population_SE)
    _write(population, tmp_path / "population.sumstats.gz")
    derived = _load(tmp_path, "population")
    assert derived[5]["mode"] == "direct-population"
    np.testing.assert_allclose(derived[0], average[0], atol=1e-14)
    np.testing.assert_allclose(derived[1], average[1], atol=1e-14)


def test_parental_and_common_weight_moments_share_joint_sampling_covariance(tmp_path):
    frame = _frame(120)
    difference = np.random.default_rng(28).normal(0, .4, len(frame))
    frame["paternal_Beta"] = frame.avg_NTC_Beta + difference
    frame["maternal_Beta"] = frame.avg_NTC_Beta - difference
    frame["paternal_SE"] = np.linspace(.1, .3, len(frame))
    frame["maternal_SE"] = .15
    frame["r_direct_paternal"] = frame["r_direct_maternal"] = -.2
    frame["r_paternal_maternal"] = .1
    _write(frame, tmp_path / "input.sumstats.gz")
    result = cli.main([str(tmp_path / "input"), str(tmp_path / "fit"), "--ld-weight", "none",
                       "--n-blocks", "4", "--save-delete", "--common-weight"])
    diagnostics = json.loads((tmp_path / "fit_diagnostics.json").read_text())
    assert diagnostics["input"]["mode"] == "parental"
    assert {"common_weight_min_eigenvalue", "common_weight_is_psd"} <= diagnostics.keys()

    input_names = result["moment_names"]
    raw_names = ["common_weight_moment_" + name for name in input_names]
    metric_names = result["metric_names"]
    additional = [name for name in input_names if name not in metric_names]
    assert len(additional) == 5
    expanded = ["var_direct", "cov_direct_avg_NTC", "cov_direct_population",
                "var_avg_NTC", "cov_avg_NTC_population", "var_population"]
    moments = pd.read_csv(tmp_path / "fit_moments.tsv", sep="\t").set_index("moment")
    assert moments.index.is_unique
    assert set(moments.index) == set(expanded + additional + raw_names)
    metrics = pd.read_csv(tmp_path / "fit_corrs.txt", sep="\t").set_index("correlation")
    np.testing.assert_allclose(moments.loc[expanded], metrics.loc[expanded])
    np.testing.assert_allclose(moments.loc[input_names, "est"], result["moments"])
    np.testing.assert_allclose(moments.loc[input_names, "SE"], result["moment_se"])
    np.testing.assert_allclose(moments.loc[raw_names, "est"], result["common_weight_moments"])
    np.testing.assert_allclose(moments.loc[raw_names, "SE"], result["common_weight_moment_se"])

    sampling = pd.read_csv(tmp_path / "fit_sampling_covariance.tsv", sep="\t", index_col="statistic")
    assert list(sampling.index) == list(sampling.columns) == metric_names + additional + raw_names
    deletes = pd.read_csv(tmp_path / "fit_delete.txt", sep="\t")
    np.testing.assert_allclose(deletes[raw_names], result["common_weight_moment_delete"])
    unique_columns = metric_names + ["moment_" + name for name in additional] + raw_names
    values = deletes[unique_columns].to_numpy()
    centered = values - values.mean(axis=0)
    np.testing.assert_allclose(sampling, (len(values) - 1) / len(values) * centered.T @ centered, atol=1e-14)
    np.testing.assert_allclose(np.diag(sampling.loc[moments.index, moments.index]), moments.SE ** 2, atol=1e-14)
    assert {p.name for p in tmp_path.glob("fit_*.tsv")} == {"fit_moments.tsv", "fit_sampling_covariance.tsv"}


def test_missing_parental_covariance_never_imputed_and_sibling_uses_average(tmp_path):
    frame = _frame()
    for name in ("paternal", "maternal"):
        frame[name + "_Beta"] = frame.avg_NTC_Beta
        frame[name + "_SE"] = .2
        frame["r_direct_" + name] = -.2
    _write(frame, tmp_path / "input.sumstats.gz")
    assert _load(tmp_path)[5]["mode"] == "average"
    with pytest.raises(ValueError, match="r_paternal_maternal"):
        _load(tmp_path, options=["--mode", "parental"])
    frame["r_paternal_maternal"] = .25
    _write(frame, tmp_path / "input.sumstats.gz")
    assert _load(tmp_path)[0].shape == (24, 3)
    frame["sib_Beta"] = .01
    _write(frame, tmp_path / "input.sumstats.gz")
    assert _load(tmp_path)[5]["mode"] == "average"
    with pytest.raises(ValueError, match="sibling"):
        _load(tmp_path, options=["--mode", "parental"])


def test_aliases_explicit_covariance_scaling_and_constant_correlation(tmp_path):
    frame = _frame()
    frame["cov_direct_avg_NTC"] = frame.r_direct_avg_NTC * frame.direct_SE * frame.avg_NTC_SE
    frame = frame.drop(columns="r_direct_avg_NTC").rename(columns={
        "SNP": "rsid", "pos": "BP", "freq": "EAF",
        "direct_Beta": "theta_direct", "direct_SE": "se_direct",
        "avg_NTC_Beta": "theta_avgNTC", "avg_NTC_SE": "se_avgNTC",
        "cov_direct_avg_NTC": "sampling_cov_avgNTC_direct",
    })
    _write(frame, tmp_path / "input.sumstats.gz")
    baseline = _load(tmp_path, options=["--corr_filter", "6"])
    assert len(baseline[0]) == 24
    scaled = _load(tmp_path, options=["--scale", "genotype-sd"])
    np.testing.assert_allclose(scaled[0], baseline[0] * np.sqrt(.32))
    np.testing.assert_allclose(scaled[1], baseline[1] * .32)
    restricted = _load(tmp_path, options=["--chr_range", "2"])
    assert restricted[5]["qc"]["excluded_chr"] == 12
    assert len(restricted[0]) == 12


def test_mixed_covariance_sources_across_chromosomes_and_pattern_validation(tmp_path):
    frame = _frame()
    for chromosome in (1, 2):
        rows = frame.loc[frame.chromosome == chromosome].copy()
        if chromosome == 2:
            rows["cov_direct_avg_NTC"] = rows.r_direct_avg_NTC * rows.direct_SE * rows.avg_NTC_SE
            rows = rows.drop(columns="r_direct_avg_NTC")
        _write(rows, tmp_path / "chr{}.sumstats.gz".format(chromosome))
    assert len(_load(tmp_path, "chr@", ["--chr_range", "1-2"])[0]) == len(frame)
    _write(frame.loc[frame.chromosome == 2], tmp_path / "chr1.sumstats.gz")
    with pytest.raises(ValueError, match="inconsistent"):
        _load(tmp_path, "chr@", ["--chr_range", "1-2"])


def test_indefinite_parental_sampling_covariance_and_duplicate_snps_rejected(tmp_path):
    frame = _frame()
    for name in ("paternal", "maternal"):
        frame[name + "_Beta"] = frame.avg_NTC_Beta
        frame[name + "_SE"] = .2
        frame["r_direct_" + name] = .9
    frame["r_paternal_maternal"] = -.9
    _write(frame, tmp_path / "input.sumstats.gz")
    with pytest.raises(ValueError, match="non-positive-definite"):
        _load(tmp_path)
    frame.loc[1, "SNP"] = frame.loc[0, "SNP"]
    _write(frame, tmp_path / "input.sumstats.gz")
    with pytest.raises(ValueError, match="unique SNP IDs"):
        _load(tmp_path)


def test_ld_weighting_must_be_explicit(tmp_path):
    completed = _run_cli(["unused", str(tmp_path / "out")], tmp_path)
    assert completed.returncode == 2
    assert "--ldscores" in completed.stderr
    assert "--ld-weight none" in completed.stderr


def test_explicit_duplicate_exclusion_drops_all_copies_and_records_ids(tmp_path):
    frame = _frame()
    frame.loc[1, "SNP"] = frame.loc[0, "SNP"]
    # Conflicting coefficients must not result in an arbitrary first/last pick.
    frame.loc[1, "direct_Beta"] = 10
    _write(frame, tmp_path / "input.sumstats.gz")
    result = _load(tmp_path, options=["--duplicate-snps", "drop"])
    np.testing.assert_allclose(result[0], frame.loc[2:, ["direct_Beta", "avg_NTC_Beta"]])
    assert result[5]["qc"]["rows_read"] == len(frame)
    assert result[5]["qc"]["excluded_duplicate_id_rows"] == 2
    assert result[5]["excluded_duplicate_snp_ids"] == ["rs0"]


def test_legacy_parser_options_and_defaults():
    defaults = cli.parser.parse_args(["input", "out"])
    assert defaults.chr_range is None
    assert defaults.min_maf == .05
    assert defaults.corr_filter == 6
    assert defaults.n_blocks == 200
    assert defaults.ld_wind == 1
    assert defaults.threads is None
    assert not defaults.save_delete
    assert not defaults.common_weight
    assert cli.parser.parse_args(["input", "out", "--common_weight"]).common_weight
    args = cli.parser.parse_args([
        "input", "out", "--chr_range", "1-2", "4", "--bed", "chr#",
        "--threads", "3", "--min_maf", ".1", "--corr_filter", "5",
        "--n_blocks", "7", "--save_delete", "--ld_wind", "2", "--ld_out", "scores",
    ])
    assert cli.parse_chromosomes(args.chr_range) == [1, 2, 4]
    assert cli.parse_chromosomes([1, 2, 4]) == [1, 2, 4]
    assert args.bed == "chr#"
    assert args.threads == 3
    assert args.min_maf == .1
    assert args.corr_filter == 5
    assert args.n_blocks == 7
    assert args.save_delete
    assert args.ld_wind == 2
    assert args.ld_out == "scores"
    assert cli.parser.parse_args(["input", "out", "--chr_range"]).chr_range == []
    assert cli.parser.parse_args(["input", "out", "--no-corr-filter"]).corr_filter is None


def test_legacy_namespace_receives_new_defaults(tmp_path):
    import argparse

    frame = _frame(80)
    _write(frame, tmp_path / "input.sumstats.gz")
    _write(pd.DataFrame({"SNP": frame.SNP, "L2": 1 + np.arange(len(frame)) / 80}),
           tmp_path / "scores.l2.ldscore.gz")
    args = argparse.Namespace(
        sumstats=str(tmp_path / "input"), out=str(tmp_path / "fit"),
        chr_range=[1, 2], ldscores=str(tmp_path / "scores"), bed=None,
        threads=1, min_maf=.05, corr_filter=6., n_blocks=4, save_delete=True,
        ld_wind=1., ld_out=None,
    )
    result = cli.main(args)
    assert result["diagnostics"]["n_blocks"] == 4
    assert result["delete_values"].shape == (4, len(result["metric_names"]))
    metadata = json.loads((tmp_path / "fit_diagnostics.json").read_text())["input"]
    assert metadata["mode"] == "average"
    assert metadata["effect_scale"] == "per-allele"
    assert metadata["ld_weight"] == "inverse"
    deletes = pd.read_csv(tmp_path / "fit_delete.txt", sep="\t", dtype={"block": str})
    assert list(deletes["block"]) == result["diagnostics"]["block_labels"]
    np.testing.assert_allclose(deletes[result["metric_names"]], result["delete_values"], atol=1e-15)
    columns = ["moment_" + name for name in result["moment_names"]]
    np.testing.assert_allclose(deletes[columns], result["moment_delete"], atol=1e-15)
    assert "common_weight_moment_delete" not in result


def test_single_file_without_chromosomes_keeps_supplied_genomic_order(tmp_path):
    frame = _frame().drop(columns="chromosome")
    frame["pos"] = np.tile([300, 200, 100], len(frame) // 3)
    _write(frame, tmp_path / "input.sumstats.gz")
    prepared = _load(tmp_path)
    np.testing.assert_allclose(prepared[0], frame[["direct_Beta", "avg_NTC_Beta"]],
                               rtol=1e-14, atol=1e-15)
    assert prepared[5]["chromosome_unavailable"]
    assert prepared[5]["qc"]["rows_retained"] == len(frame)


@pytest.mark.parametrize("requested, expected", [(None, 4), (2, 2), (100, 4)])
def test_thread_count_is_set_and_capped(tmp_path, monkeypatch, requested, expected):
    numba = pytest.importorskip("numba")
    frame = _frame(60)
    _write(frame, tmp_path / "input.sumstats.gz")
    calls = []
    monkeypatch.setattr(numba.config, "NUMBA_NUM_THREADS", 4)
    monkeypatch.setattr(numba, "set_num_threads", calls.append)
    arguments = [str(tmp_path / "input"), str(tmp_path / "fit"),
                 "--ld-weight", "none", "--n_blocks", "3"]
    if requested is not None:
        arguments += ["--threads", str(requested)]
    cli.main(arguments)
    assert calls == [expected]


def test_bed_scores_preserve_qc_alignment_and_match_precomputed_fit(tmp_path, monkeypatch):
    from snipar import coherent_ld

    frame = _frame(80).sample(frac=1, random_state=3)
    frame.loc[0, "freq"] = .001
    frame.loc[1, "direct_SE"] = np.nan
    _write(frame, tmp_path / "input.sumstats.gz")
    ld = pd.Series(1 + np.arange(len(frame)) / 80, index=frame.SNP)
    _write(ld.rename("L2").rename_axis("SNP").reset_index(),
           tmp_path / "scores.l2.ldscore.gz")
    calls = []

    def fake_bed(pattern, snp_ids, chromosomes, ld_wind, ld_out):
        calls.append((pattern, snp_ids.copy(), chromosomes.copy(), ld_wind, ld_out))
        return ld.reindex(snp_ids).to_numpy(), {"bed_files": ["chr1.bed", "chr2.bed"]}

    monkeypatch.setattr(coherent_ld, "compute_ldscores", fake_bed)
    common = [str(tmp_path / "input"), str(tmp_path / "fit"), "--n_blocks", "4"]
    from_bed = cli.main(common + ["--bed", "chr#", "--ld_wind", "2.5", "--ld_out", "computed"])
    from_scores = cli.main(common + ["--ldscores", str(tmp_path / "scores")])
    assert {p.name for p in tmp_path.glob("fit_*")} == {
        "fit_corrs.txt", "fit_moments.tsv", "fit_sampling_covariance.tsv", "fit_diagnostics.json",
    }
    assert len(calls) == 1
    pattern, snp_ids, chromosomes, window, output = calls[0]
    assert pattern == "chr#"
    expected = frame.loc[~frame.index.isin([0, 1])]
    np.testing.assert_array_equal(snp_ids, expected.SNP)
    np.testing.assert_array_equal(chromosomes, expected.chromosome)
    assert window == 2.5
    assert output == "computed"
    assert from_bed["metric_names"] == from_scores["metric_names"]
    for key in ("estimates", "se", "delete_values", "covariance"):
        np.testing.assert_allclose(from_bed[key], from_scores[key], rtol=1e-12, atol=1e-14)


def test_requested_chromosome_range_accepts_available_partial_panel(tmp_path):
    frame = _frame(24).loc[lambda data: data.chromosome == 1].copy()
    _write(frame, tmp_path / "chr1.sumstats.gz")
    _write(pd.DataFrame({"SNP": frame.SNP, "L2": np.linspace(1, 2, len(frame))}),
           tmp_path / "ld1.l2.ldscore.gz")
    args = cli.parser.parse_args([
        str(tmp_path / "chr@"), str(tmp_path / "fit"),
        "--chr_range", "1-22", "--ldscores", str(tmp_path / "ld@"),
    ])
    prepared = cli.prepare_inputs(args)
    np.testing.assert_allclose(prepared[0], frame[["direct_Beta", "avg_NTC_Beta"]])
    np.testing.assert_allclose(prepared[2], np.linspace(1, 2, len(frame)))
    assert prepared[5]["sumstats_files"] == [str(tmp_path / "chr1.sumstats.gz")]
    assert prepared[5]["ldscore_files"] == [str(tmp_path / "ld1.l2.ldscore.gz")]
    args.chr_range = ["3-22"]
    with pytest.raises(ValueError, match="No files matched"):
        cli.prepare_inputs(args)
