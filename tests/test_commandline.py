import unittest
import subprocess
from tests.test_imputation import imputation_test

class TestCommanline(unittest.TestCase):
    p_value_threshold = 0.01
    def test_impute_runner_with_pedigree(self):
        command = ["python",
                   "impute_runner.py",
                   "1",
                   "3",
                   "test_data/sample.segments.gz",
                   "test_data/sample~",
                   "--pedigree", "test_data/sample.ped",
                   "--out_prefix", "outputs/tmp/test_sample_imputed",
                   ]
        subprocess.check_call(command)
    
    def test_impute_runner_with_pedigree_control(self):
        command = ["python",
                   "impute_runner.py",
                   "-c",
                   "1",
                   "3",
                   "test_data/sample.segments.gz",
                   "test_data/sample~",
                   "--pedigree", "test_data/sample.ped",
                   "--out_prefix", "outputs/tmp/test_sample_imputed",
                   ]
        subprocess.check_call(command)
        coef, z, p_value = imputation_test([1, 2],
                imputed_prefix = "outputs/tmp/test_sample_imputed",
                expected_prefix = "test_data/sample",
                )
        self.assertGreaterEqual(p_value[0], self.p_value_threshold)
        self.assertGreaterEqual(p_value[1], self.p_value_threshold)

    def test_impute_runner_with_king(self):
        command = ["python",
                   "impute_runner.py",
                   "1",
                   "3",
                   "test_data/sample.segments.gz",
                   "test_data/sample~",
                   "--king", "test_data/sample.king",
                   "--agesex", "test_data/sample.agesex",
                   "--out_prefix", "outputs/tmp/test_sample_imputed",
                   ]
        subprocess.check_call(command)

    # TODO handle nothing to impute errors
    def test_impute_runner_with_king_control(self):
        command = ["python",
                   "impute_runner.py",
                   "-c",
                   "1",
                   "3",
                   "test_data/sample.segments.gz",
                   "test_data/sample~",
                   "--king", "test_data/sample.king",
                   "--agesex", "test_data/sample.agesex",
                   "--out_prefix", "outputs/tmp/test_sample_imputed",
                   ]
        subprocess.check_call(command)
        coef, z, p_value = imputation_test([1, 2],
                imputed_prefix = "outputs/tmp/test_sample_imputed",
                expected_prefix = "test_data/sample",
                )
        self.assertGreaterEqual(p_value[0], self.p_value_threshold)
        self.assertGreaterEqual(p_value[1], self.p_value_threshold)


    def test_impute_runner_with_pedigree_control_multicore(self):
        command = ["python",
                   "impute_runner.py",
                   "-c",
                   "1",
                   "3",
                   "test_data/sample.segments.gz",
                   "test_data/sample~",
                   "--pedigree", "test_data/sample.ped",
                   "--out_prefix", "outputs/tmp/test_sample_imputed",
                   "--threads", "2",
                   ]
        subprocess.check_call(command)
        coef, z, p_value = imputation_test([1, 2],
                imputed_prefix = "outputs/tmp/test_sample_imputed",
                expected_prefix = "test_data/sample",
                )
        self.assertGreaterEqual(p_value[0], self.p_value_threshold)
        self.assertGreaterEqual(p_value[1], self.p_value_threshold)

    def test_impute_runner_with_pedigree_control_notilda(self):
        command = ["python",
                   "impute_runner.py",
                   "-c",
                   "1",
                   "2",
                   "test_data/sample.segments.gz",
                   "test_data/sample1",
                   "--pedigree", "test_data/sample.ped",
                   "--out_prefix", "outputs/tmp/test_sample_imputed",
                   "--threads", "2",
                   ]
        subprocess.check_call(command)
        coef, z, p_value = imputation_test([1],
                imputed_prefix = "outputs/tmp/test_sample_imputed",
                expected_prefix = "test_data/sample",
                )
        self.assertGreaterEqual(p_value[0], self.p_value_threshold)
        self.assertGreaterEqual(p_value[1], self.p_value_threshold)

