import subprocess
from snipar.tests.utils import *
import os
class TestExample(SniparTest):
    
    def test_example(self):
        commands1 = [
            f"rm -rf {output_root}/example_data",
            f"snipar_example_data.py --dest {output_root}/example_data",
        ]
        commands2 = [            
            "ibd.py --bed chr_@ --king king.kin0 --agesex agesex.txt --out chr_@ --threads 4 --ld_out",#TODO add option of supressing logging in snipar
            "impute.py --ibd chr_@.ibd --bed chr_@ --king king.kin0 --agesex agesex.txt --out chr_@ --threads 4 -silent_progress",
            "gwas.py phenotype.txt --bed chr_@ --imp chr_@ --threads 4",
            "python estimate_sim_effects.py chr_1.sumstats.hdf5 phenotype.effects.txt",
            "correlate.py chr_@ effect --ldscores chr_@",
            "pgs.py direct --bed chr_@ --imp chr_@ --weights direct_weights.txt",
        ]
        stdout = subprocess.DEVNULL
        if self.log:
            stdout = None
        os.chdir(f"{output_root}")
        for c in commands1:
            subprocess.run(c.split(), stdout=stdout)
        os.chdir(f"{output_root}/example_data")
        for c in commands2:
            subprocess.run(c.split(), stdout=stdout)