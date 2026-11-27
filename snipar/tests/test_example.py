import shlex
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

from snipar.tests.utils import SniparTest


class TestExample(SniparTest):
    def test_example(self):
        with TemporaryDirectory() as directory:
            example = Path(directory) / "example_data"
            stdout = None if self.log else subprocess.DEVNULL
            subprocess.run(
                [sys.executable, "-m", "snipar.example.snipar_example_data", "--dest", str(example)],
                cwd=directory, stdout=stdout, check=True,
            )
            commands = [
                "-m snipar.scripts.ibd --bed chr_@ --king king.kin0 --agesex agesex.txt --out chr_@ --threads 4 --ld_out",
                "-m snipar.scripts.impute --ibd chr_@.ibd --bed chr_@ --king king.kin0 --agesex agesex.txt --out chr_@ --threads 4 -silent_progress",
                "-m snipar.scripts.gwas phenotype.txt --bed chr_@ --imp chr_@ --threads 4 --out chr_@",
                "estimate_sim_effects.py chr_1.sumstats.hdf5 phenotype.effects.txt",
                "-m snipar.scripts.correlate chr_@ effect --ldscores chr_@",
                "-m snipar.scripts.pgs direct --bed chr_@ --imp chr_@ --weights direct_weights.txt --A1 nt1 --A2 nt2 --beta_col ldpred_beta --SNP sid",
            ]
            for command in commands:
                subprocess.run([sys.executable, *shlex.split(command)], cwd=example,
                               stdout=stdout, check=True)
            for filename in ("chr_1.sumstats.gz", "effect_corrs.txt", "direct.pgs.txt"):
                self.assertGreater((example / filename).stat().st_size, 0)
