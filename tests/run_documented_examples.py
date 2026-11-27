"""Run the tutorial and simulation exercise against the installed SNIPar commands."""
import json
import re
import shutil
import subprocess
from pathlib import Path
from tempfile import TemporaryDirectory

import h5py
import numpy as np
import pandas as pd


DOCS = Path(__file__).resolve().parents[1] / "docs"


def commands(page):
    blocks = re.findall(r"^    ``(.+?)``$", (DOCS / page).read_text(), re.MULTILINE)
    # The simulation page also shows terminal output lines beginning with a count.
    return [line for line in blocks if not line[0].isdigit()]


def run(lines, directory):
    # Keep cd, shell variables and process substitution exactly as a bash user enters them.
    script = "set -euxo pipefail\n" + "\n".join(lines)
    subprocess.run(["bash", "-c", script], cwd=directory, check=True)


def table(path, **kwargs):
    return pd.read_csv(path, sep=r"\s+", **kwargs)


def check_gwas(directory, names):
    for name in names:
        summary = table(directory / f"{name}.sumstats.gz")
        assert len(summary) == 1000 and summary.SNP.is_unique, name
        assert np.isfinite(summary.direct_Beta).all() and summary.direct_SE.gt(0).all(), name
        with h5py.File(directory / f"{name}.sumstats.hdf5") as h5:
            assert np.isfinite(h5["estimate"][:]).all(), name
            variances = np.diagonal(h5["estimate_covariance"][:], axis1=1, axis2=2)
            ses = h5["estimate_ses"][:]
            # SEs and covariances are rounded separately when written to HDF5.
            epsilon = max(np.finfo(variances.dtype).eps, np.finfo(ses.dtype).eps)
            np.testing.assert_allclose(variances, ses ** 2, rtol=4 * epsilon)


def check_pgs(directory, prefix):
    effects = table(directory / f"{prefix}.effects.txt", index_col="statistic")
    assert list(effects.columns) == ["est", "SE"] and effects.index.is_unique, prefix
    assert np.isfinite(effects).all().all() and effects.SE.gt(0).all(), prefix
    for model in (1, 2):
        cov = table(directory / f"{prefix}.{model}.vcov.txt", index_col=0)
        np.testing.assert_allclose(np.diag(cov), effects.loc[[f"{model}.{t}" for t in cov.index], "SE"] ** 2)
    assert (directory / f"{prefix}.1_2.cross_vcov.txt").is_file()
    paired = table(directory / f"{prefix}.direct_population.vcov.txt", index_col=0)
    pop, direct = effects.loc[["1.proband", "2.proband"], "est"]
    gradient = np.array([-direct / pop**2, 1 / pop])
    np.testing.assert_allclose(effects.loc["direct_population_ratio"],
                               [direct / pop, np.sqrt(gradient @ paired.values @ gradient)])
    return effects


def check_binary(directory, prefix, expected_counts):
    info = table(directory / f"{prefix}.model_info.tsv")
    assert len(info) == 2 and info.convergence.eq("OK").all() and info.ratio_status.eq("ok").all()
    assert (info[["n", "cases", "controls", "n_families"]].values == expected_counts).all()


def tutorial(root):
    lines = commands("tutorial.rst")
    run(lines[:1], root)
    example = root / "example_data"
    # Alternative IBD/imputation paths write identical filenames.
    for branch in ("bgen", "pedigree"):
        shutil.copytree(example, root / branch)
    primary = [lines[1]]
    for line in lines[2:]:
        if line.startswith("ibd.py ") and "--pedigree " in line:
            primary.append(f"(cd ../pedigree && {line})")
        elif "--bgen " in line:
            primary.append(f"(cd ../bgen && {line})")
        else:
            primary.append(line)
        if line.startswith("pgs.py direct --bed "):
            primary.append("cp direct.effects.txt ../from_genotypes.effects.txt")
    run(primary, root)
    for branch in (example, root / "bgen", root / "pedigree"):
        assert (branch / "chr_1.ibd.segments.gz").stat().st_size > 0
        assert len(table(branch / "chr_1.l2.ldscore.gz")) == 1000
    for branch in (example, root / "bgen"):
        with h5py.File(branch / "chr_1.hdf5") as h5:
            gts = h5["imputed_par_gts"][:]
            assert gts.shape == (2000, 1000) and np.isfinite(gts).all()
    check_gwas(example, ["chr_1_young"])
    check_gwas(root / "bgen", ["chr_1_young", "chr_1_unified", "chr_1_robust", "sibs_chr_1",
                               "trios_sibs_chr_1", "trios_sibs_singletons_chr_1"])
    est = table(example / "effect_corrs.txt", index_col=0)
    diagnostics = json.loads((example / "effect_diagnostics.json").read_text())
    assert diagnostics["converged"] and diagnostics["n_snps"] == 1000
    assert table(example / "effect_delete.txt").block.nunique() == diagnostics["n_blocks"] == 200
    assert est.shape == (17, 2) and np.isfinite(est).all().all()
    reference = {"r_direct_population": .895, "r_direct_avg_NTC": .471, "sqrt_var_ratio_population_direct": 1.523}
    for label, value in reference.items():
        assert abs(est.loc[label, "est"] - value) < est.loc[label, "SE"], label
    effects = check_pgs(example, "direct")
    from_genotypes = table(root / "from_genotypes.effects.txt", index_col="statistic")
    np.testing.assert_allclose(effects, from_genotypes, rtol=1e-5, atol=1e-5)
    assert abs(effects.loc["2.proband", "est"] - 1) < .05
    binary = check_pgs(example, "direct_binary")
    check_binary(example, "direct_binary", [5014, 2494, 2520, 3000])
    np.testing.assert_allclose(binary.loc["direct_population_ratio"], [.777, .023], atol=.002)
    return len(lines)


def simulation(root):
    lines = commands("simulation.rst")
    run(lines, root)
    sim = root / "sim"
    phenotypes = table(sim / "phenotype.txt")
    assert len(phenotypes) == 6000 and phenotypes.FID.nunique() == 3000
    check_gwas(sim, [f"chr_1_{m}" for m in ("sibdiff", "trio", "imp", "phased", "unified", "no_singletons")])
    for prefix in ("direct", "direct_v1", "direct_v1_obs", "direct_v1_binary"):
        check_pgs(sim, prefix)
    check_binary(sim, "direct_v1_binary", [6000, 3000, 3000, 3000])
    for name, neutral in (("sim_stabilizing", False), ("sim_neutral", True)):
        folder = sim / name
        meta = json.loads((folder / "selection_metadata.json").read_text())
        assert meta.get("neutral_control", False) is neutral and meta["n_generations"] == 3
        assert meta["selection_variance"] == (None if neutral else 2.0)
        assert len(table(folder / "phenotype.txt")) == 400
        phases = {(row["generation"], row["phase"]): row for row in meta["history"]}
        for generation in range(3):
            before = phases[generation, "before_selection"]["n_individuals"]
            after = phases[generation, "after_selection"]["n_individuals"]
            if neutral:
                assert after == before
            else:
                assert 0 < after < before
    return len(lines)


if __name__ == "__main__":
    with TemporaryDirectory(prefix="snipar-examples-") as directory:
        root = Path(directory)
        total = tutorial(root) + simulation(root)
        print(f"Passed {total} documented commands and their output checks.")
