import unittest
import subprocess
from tests.test_imputation import imputation_test
#TODO add tests with nan
class TestCommanline(unittest.TestCase):
    p_value_threshold = 0.01
    def test_impute_runner_with_unphased_pedigree(self):
        command = ["python",
                   "impute_runner.py",
                   "test_data/sample",
                   "--bed", "test_data/sample~",
                   "--from_chr", "1",
                   "--to_chr", "3",
                   "--pedigree", "test_data/sample.ped",
                   "--threads", "10",
                   "--output_address", "outputs/tmp/test_impute_runner_with_unphased_pedigree~",
                   ]
        subprocess.check_call(command)
    
    def test_impute_runner_with_unphased_pedigree_chunks_control(self):
        command = ["python",
                   "impute_runner.py",
                   "-c",
                   "test_data/sample",
                   "--bed", "test_data/sample~",
                   "--from_chr", "1",
                   "--to_chr", "3",
                   "--pedigree", "test_data/sample.ped",
                   "--chunks", "7",
                   "--threads", "10",
                   "--output_address", "outputs/tmp/test_impute_runner_with_unphased_pedigree_control~",
                   ]
        subprocess.check_call(command)
        coef, z, p_value = imputation_test([1],
                imputed_prefix = "outputs/tmp/test_impute_runner_with_unphased_pedigree_control",
                expected_prefix = "test_data/sample",
                )
        self.assertGreaterEqual(p_value[0], self.p_value_threshold)
        self.assertGreaterEqual(p_value[1], self.p_value_threshold)

    def test_impute_runner_with_unphased_pedigree_control(self):
        command = ["python",
                   "impute_runner.py",
                   "-c",
                   "test_data/sample",
                   "--bed", "test_data/sample~",
                   "--from_chr", "1",
                   "--to_chr", "3",
                   "--pedigree", "test_data/sample.ped",
                   "--threads", "10",
                   "--output_address", "outputs/tmp/test_impute_runner_with_unphased_pedigree_control~",
                   ]
        subprocess.check_call(command)
        coef, z, p_value = imputation_test([1],
                imputed_prefix = "outputs/tmp/test_impute_runner_with_unphased_pedigree_control",
                expected_prefix = "test_data/sample",
                )
        self.assertGreaterEqual(p_value[0], self.p_value_threshold)
        self.assertGreaterEqual(p_value[1], self.p_value_threshold)

    def test_impute_runner_with_phased_pedigree_control(self):
        command = ["python",
                   "impute_runner.py",
                   "-c",
                   "test_data/sample",
                   "--bgen", "test_data/sample~",
                   "--from_chr", "1",
                   "--to_chr", "3",
                   "--pedigree", "test_data/sample.ped",
                   "--threads", "10",
                   "--output_address", "outputs/tmp/test_impute_runner_with_phased_pedigree_control~",
                   ]
        subprocess.check_call(command)
        coef, z, p_value = imputation_test([1, 2],
                imputed_prefix = "outputs/tmp/test_impute_runner_with_phased_pedigree_control",
                expected_prefix = "test_data/sample",
                )
        self.assertGreaterEqual(p_value[0], self.p_value_threshold)
        self.assertGreaterEqual(p_value[1], self.p_value_threshold)

    def test_impute_runner_with_phased_pedigree_chunks_control(self):
        command = ["python",
                   "impute_runner.py",
                   "-c",
                   "test_data/sample",
                   "--bgen", "test_data/sample~",
                   "--from_chr", "1",
                   "--to_chr", "3",
                   "--chunks", "7",
                   "--pedigree", "test_data/sample.ped",
                   "--threads", "10",
                   "--output_address", "outputs/tmp/test_impute_runner_with_phased_pedigree_control~",
                   ]
        subprocess.check_call(command)
        coef, z, p_value = imputation_test([1, 2],
                imputed_prefix = "outputs/tmp/test_impute_runner_with_phased_pedigree_control",
                expected_prefix = "test_data/sample",
                )
        self.assertGreaterEqual(p_value[0], self.p_value_threshold)
        self.assertGreaterEqual(p_value[1], self.p_value_threshold)

    def test_impute_runner_with_unphased_king(self):
        command = ["python",
                   "impute_runner.py",
                   "test_data/sample",
                   "--bed", "test_data/sample~",
                   "--from_chr", "1",
                   "--to_chr", "3",
                   "--king", "test_data/sample.king",
                   "--agesex", "test_data/sample.agesex",
                   "--threads", "10",
                   "--output_address", "outputs/tmp/test_impute_runner_with_unphased_king~",
                   ]
        subprocess.check_call(command)

    def test_impute_runner_with_unphased_king_control(self):
        command = ["python",
                   "impute_runner.py",
                   "-c",
                   "test_data/sample",
                   "--bed", "test_data/sample~",
                   "--from_chr", "1",
                   "--to_chr", "3",
                   "--king", "test_data/sample.king",
                   "--agesex", "test_data/sample.agesex",
                   "--threads", "10",
                   "--output_address", "outputs/tmp/test_impute_runner_with_unphased_king_control~",
                   ]
        subprocess.check_call(command)
        coef, z, p_value = imputation_test([1, 2],
                imputed_prefix = "outputs/tmp/test_impute_runner_with_unphased_king_control",
                expected_prefix = "test_data/sample",
                )
        self.assertGreaterEqual(p_value[0], self.p_value_threshold)
        self.assertGreaterEqual(p_value[1], self.p_value_threshold)


    def test_impute_runner_with_unphased_pedigree_control_nothread(self):
        command = ["python",
                   "impute_runner.py",
                   "-c",
                   "test_data/sample",
                   "--bed", "test_data/sample~",
                   "--from_chr", "1",
                   "--to_chr", "3",
                   "--pedigree", "test_data/sample.ped",
                   "--output_address", "outputs/tmp/test_impute_runner_with_unphased_pedigree_control_multithread~",
                   ]
        subprocess.check_call(command)
        coef, z, p_value = imputation_test([1, 2],
                imputed_prefix = "outputs/tmp/test_impute_runner_with_unphased_pedigree_control_multithread",
                expected_prefix = "test_data/sample",
                )
        self.assertGreaterEqual(p_value[0], self.p_value_threshold)
        self.assertGreaterEqual(p_value[1], self.p_value_threshold)

    def test_impute_runner_with_unphased_pedigree_control_multiprocess(self):
        command = ["python",
                   "impute_runner.py",
                   "-c",
                   "test_data/sample",
                   "--bed", "test_data/sample~",
                   "--from_chr", "1",
                   "--to_chr", "3",
                   "--pedigree", "test_data/sample.ped",
                   "--threads", "10",
                   "--output_address", "outputs/tmp/test_impute_runner_with_unphased_pedigree_control_multiprocess~",
                   "--processes", "2",
                   ]
        subprocess.check_call(command)
        coef, z, p_value = imputation_test([1, 2],
                imputed_prefix = "outputs/tmp/test_impute_runner_with_unphased_pedigree_control_multiprocess",
                expected_prefix = "test_data/sample",
                )
        self.assertGreaterEqual(p_value[0], self.p_value_threshold)
        self.assertGreaterEqual(p_value[1], self.p_value_threshold)

    def test_impute_runner_with_unphased_pedigree_control_notilda(self):
        command = ["python",
                   "impute_runner.py",
                   "-c",
                   "test_data/sample",
                   "--bed", "test_data/sample1",
                   "--pedigree", "test_data/sample.ped",
                   "--output_address", "outputs/tmp/test_impute_runner_with_unphased_pedigree_control_notilda1",
                   "--threads", "10",
                   ]
        subprocess.check_call(command)
        coef, z, p_value = imputation_test([1],
                imputed_prefix = "outputs/tmp/test_impute_runner_with_unphased_pedigree_control_notilda",
                expected_prefix = "test_data/sample",
                )
        self.assertGreaterEqual(p_value[0], self.p_value_threshold)
        self.assertGreaterEqual(p_value[1], self.p_value_threshold)

    #TODO raise error in multichromosome files