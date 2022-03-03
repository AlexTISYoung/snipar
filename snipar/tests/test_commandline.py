import unittest
import subprocess
from snipar.tests.test_imputation import imputation_test
import os
#TODO add tests with nan
tests_root = os.path.dirname(__file__)
output_root = os.path.join(tests_root, "tmp")
if not os.path.exists(output_root):
    os.mkdir(output_root)
class TestCommanline(unittest.TestCase):
    p_value_threshold = 0.01
    def test_impute_with_unphased_pedigree_control(self):
        command = ["impute.py",
                   "-c",
                   "--ibd", f"{tests_root}/test_data/sample.our",                   
                   "--bed", f"{tests_root}/test_data/sample_reduced~",
                   "--from_chr", "1",
                   "--to_chr", "3",
                   "--pedigree", f"{tests_root}/test_data/sample.ped",
                   "--threads", "10",
                   "--out", f"{output_root}/test_impute_with_unphased_pedigree_control~",
                   ]
        subprocess.check_call(command)
        coef, z, p_value = imputation_test([1],
                imputed_prefix = f"{output_root}/test_impute_with_unphased_pedigree_control",
                expected_prefix = f"{tests_root}/test_data/sample",
                )
        self.assertGreaterEqual(p_value[0], self.p_value_threshold)
        self.assertGreaterEqual(p_value[1], self.p_value_threshold)
    
    def test_impute_with_unphased_pedigree_wild_ibd_control(self):
        command = ["impute.py",
                   "-c",
                   "--ibd", f"{tests_root}/test_data/sample~.our",
                   "--bed", f"{tests_root}/test_data/sample_reduced~",
                   "--from_chr", "1",
                   "--to_chr", "3",
                   "--pedigree", f"{tests_root}/test_data/sample.ped",
                   "--threads", "10",
                   "--out", f"{output_root}/test_impute_with_unphased_pedigree_control~",
                   ]
        subprocess.check_call(command)
        coef, z, p_value = imputation_test([1],
                imputed_prefix = f"{output_root}/test_impute_with_unphased_pedigree_control",
                expected_prefix = f"{tests_root}/test_data/sample",
                )
        self.assertGreaterEqual(p_value[0], self.p_value_threshold)
        self.assertGreaterEqual(p_value[1], self.p_value_threshold)

    def test_impute_with_unphased_pedigree(self):
        command = ["impute.py",
                   "--ibd", f"{tests_root}/test_data/sample.our",
                   "--bed", f"{tests_root}/test_data/sample_reduced~",
                   "--from_chr", "1",
                   "--to_chr", "2",
                   "--pedigree", f"{tests_root}/test_data/sample.ped",
                   "--threads", "10",
                   "--out", f"{output_root}/test_impute_with_unphased_pedigree~",
                   ]
        subprocess.check_call(command)

    def test_impute_with_unphased_pedigree_chunks_control(self):
        command = ["impute.py",
                   "-c",
                   "--ibd", f"{tests_root}/test_data/sample.our",
                   "--bed", f"{tests_root}/test_data/sample_reduced~",
                   "--from_chr", "1",
                   "--to_chr", "3",
                   "--pedigree", f"{tests_root}/test_data/sample.ped",
                   "--chunks", "7",
                   "--threads", "10",
                   "--out", f"{output_root}/test_impute_with_unphased_pedigree_control~",
                   ]
        subprocess.check_call(command)
        coef, z, p_value = imputation_test([1],
                imputed_prefix = f"{output_root}/test_impute_with_unphased_pedigree_control",
                expected_prefix = f"{tests_root}/test_data/sample",
                )
        self.assertGreaterEqual(p_value[0], self.p_value_threshold)
        self.assertGreaterEqual(p_value[1], self.p_value_threshold)


    def test_impute_with_unphased_pedigree_control_legacy_ibd(self):
        command = ["impute.py",
                   "-c",
                   "--ibd", f"{tests_root}/test_data/sample.king",
                   "--ibd_is_king",
                   "--bed", f"{tests_root}/test_data/sample_reduced~",
                   "--from_chr", "1",
                   "--to_chr", "3",
                   "--pedigree", f"{tests_root}/test_data/sample.ped",
                   "--threads", "10",
                   "--out", f"{output_root}/test_impute_with_unphased_pedigree_control~",
                   ]
        subprocess.check_call(command)
        coef, z, p_value = imputation_test([1],
                imputed_prefix = f"{output_root}/test_impute_with_unphased_pedigree_control",
                expected_prefix = f"{tests_root}/test_data/sample",
                )
        self.assertGreaterEqual(p_value[0], self.p_value_threshold)
        self.assertGreaterEqual(p_value[1], self.p_value_threshold)


    def test_impute_with_phased_pedigree_control(self):
        command = ["impute.py",
                   "-c",
                   "--ibd", f"{tests_root}/test_data/sample.our",
                   "--bgen", f"{tests_root}/test_data/sample_reduced~",
                   "--from_chr", "1",
                   "--to_chr", "3",
                   "--pedigree", f"{tests_root}/test_data/sample.ped",
                   "--threads", "10",
                   "--out", f"{output_root}/test_impute_with_phased_pedigree_control~",
                   ]
        subprocess.check_call(command)
        coef, z, p_value = imputation_test([1, 2],
                imputed_prefix = f"{output_root}/test_impute_with_phased_pedigree_control",
                expected_prefix = f"{tests_root}/test_data/sample",
                )
        self.assertGreaterEqual(p_value[0], self.p_value_threshold)
        self.assertGreaterEqual(p_value[1], self.p_value_threshold)

    def test_impute_with_phased_pedigree_chunks_control(self):
        command = ["impute.py",
                   "-c",
                   "--ibd", f"{tests_root}/test_data/sample.our",
                   "--bgen", f"{tests_root}/test_data/sample_reduced~",
                   "--from_chr", "1",
                   "--to_chr", "3",
                   "--chunks", "7",
                   "--pedigree", f"{tests_root}/test_data/sample.ped",
                   "--threads", "10",
                   "--out", f"{output_root}/test_impute_with_phased_pedigree_control~",
                   ]
        subprocess.check_call(command)
        coef, z, p_value = imputation_test([1, 2],
                imputed_prefix = f"{output_root}/test_impute_with_phased_pedigree_control",
                expected_prefix = f"{tests_root}/test_data/sample",
                )
        self.assertGreaterEqual(p_value[0], self.p_value_threshold)
        self.assertGreaterEqual(p_value[1], self.p_value_threshold)

    def test_impute_with_unphased_king(self):
        command = ["impute.py",
                   "--ibd", f"{tests_root}/test_data/sample.our",
                   "--bed", f"{tests_root}/test_data/sample_reduced~",
                   "--from_chr", "1",
                   "--to_chr", "3",
                   "--king", f"{tests_root}/test_data/sample.king",
                   "--agesex", f"{tests_root}/test_data/sample.agesex",
                   "--threads", "10",
                   "--out", f"{output_root}/test_impute_with_unphased_king~",
                   ]
        subprocess.check_call(command)

    def test_impute_with_unphased_king_control(self):
        command = ["impute.py",
                   "-c",
                   "--ibd", f"{tests_root}/test_data/sample.our",
                   "--bed", f"{tests_root}/test_data/sample_reduced~",
                   "--from_chr", "1",
                   "--to_chr", "3",
                   "--king", f"{tests_root}/test_data/sample.king",
                   "--agesex", f"{tests_root}/test_data/sample.agesex",
                   "--threads", "10",
                   "--out", f"{output_root}/test_impute_with_unphased_king_control~",
                   ]
        subprocess.check_call(command)
        coef, z, p_value = imputation_test([1, 2],
                imputed_prefix = f"{output_root}/test_impute_with_unphased_king_control",
                expected_prefix = f"{tests_root}/test_data/sample",
                )
        self.assertGreaterEqual(p_value[0], self.p_value_threshold)
        self.assertGreaterEqual(p_value[1], self.p_value_threshold)

    def test_impute_with_unphased_king_control_pca(self):
        command = ["impute.py",
                   "-c",
                   "--ibd", f"{tests_root}/test_data/sample.our",
                   "--bed", f"{tests_root}/test_data/sample~",
                   "--from_chr", "1",
                   "--to_chr", "3",
                   "--king", f"{tests_root}/test_data/sample.king",
                   "--agesex", f"{tests_root}/test_data/sample.agesex",
                   "--pcs", f"{tests_root}/test_data/sample1_2_pca.eigenvec",
                   "--pc_num", "8",
                   "-find_optimal_pc",
                   "--threads", "10",
                   "--out", f"{output_root}/test_impute_with_unphased_king_control~",
                   ]
        subprocess.check_call(command)
        coef, z, p_value = imputation_test([1, 2],
                imputed_prefix = f"{output_root}/test_impute_with_unphased_king_control",
                expected_prefix = f"{tests_root}/test_data/sample",
                )
        #TODO use pvals
        self.assertAlmostEqual(coef[0], 1., delta=0.02)
        self.assertAlmostEqual(coef[1], 1., delta=0.02)

    def test_impute_with_unphased_king_control_legacy_ibd(self):
        command = ["impute.py",
                   "-c",
                   "--ibd", f"{tests_root}/test_data/sample.king",
                   "--ibd_is_king",
                   "--bed", f"{tests_root}/test_data/sample_reduced~",
                   "--from_chr", "1",
                   "--to_chr", "3",
                   "--king", f"{tests_root}/test_data/sample.king",
                   "--agesex", f"{tests_root}/test_data/sample.agesex",
                   "--threads", "10",
                   "--out", f"{output_root}/test_impute_with_unphased_king_control~",
                   ]
        subprocess.check_call(command)
        coef, z, p_value = imputation_test([1, 2],
                imputed_prefix = f"{output_root}/test_impute_with_unphased_king_control",
                expected_prefix = f"{tests_root}/test_data/sample",
                )
        self.assertGreaterEqual(p_value[0], self.p_value_threshold)
        self.assertGreaterEqual(p_value[1], self.p_value_threshold)

    def test_impute_with_unphased_king_control_legacy_tilda_ibd(self):
        command = ["impute.py",
                   "-c",
                   "--ibd", f"{tests_root}/test_data/sample~.king",
                   "--ibd_is_king",
                   "--bed", f"{tests_root}/test_data/sample_reduced~",
                   "--from_chr", "1",
                   "--to_chr", "3",
                   "--king", f"{tests_root}/test_data/sample.king",
                   "--agesex", f"{tests_root}/test_data/sample.agesex",
                   "--threads", "10",
                   "--out", f"{output_root}/test_impute_with_unphased_king_control~",
                   ]
        subprocess.check_call(command)
        coef, z, p_value = imputation_test([1, 2],
                imputed_prefix = f"{output_root}/test_impute_with_unphased_king_control",
                expected_prefix = f"{tests_root}/test_data/sample",
                )
        self.assertGreaterEqual(p_value[0], self.p_value_threshold)
        self.assertGreaterEqual(p_value[1], self.p_value_threshold)

    def test_impute_with_unphased_pedigree_control_nothread(self):
        command = ["impute.py",
                   "-c",
                   "--ibd", f"{tests_root}/test_data/sample.our",
                   "--bed", f"{tests_root}/test_data/sample_reduced~",
                   "--from_chr", "1",
                   "--to_chr", "3",
                   "--pedigree", f"{tests_root}/test_data/sample.ped",
                   "--out", f"{output_root}/test_impute_with_unphased_pedigree_control_multithread~",
                   ]
        subprocess.check_call(command)
        coef, z, p_value = imputation_test([1, 2],
                imputed_prefix = f"{output_root}/test_impute_with_unphased_pedigree_control_multithread",
                expected_prefix = f"{tests_root}/test_data/sample",
                )
        self.assertGreaterEqual(p_value[0], self.p_value_threshold)
        self.assertGreaterEqual(p_value[1], self.p_value_threshold)

    def test_impute_with_unphased_pedigree_control_multiprocess(self):
        command = ["impute.py",
                   "-c",
                   "--ibd", f"{tests_root}/test_data/sample.our",
                   "--bed", f"{tests_root}/test_data/sample_reduced~",
                   "--from_chr", "1",
                   "--to_chr", "3",
                   "--pedigree", f"{tests_root}/test_data/sample.ped",
                   "--threads", "10",
                   "--out", f"{output_root}/test_impute_with_unphased_pedigree_control_multiprocess~",
                   "--processes", "2",
                   ]
        subprocess.check_call(command)
        coef, z, p_value = imputation_test([1, 2],
                imputed_prefix = f"{output_root}/test_impute_with_unphased_pedigree_control_multiprocess",
                expected_prefix = f"{tests_root}/test_data/sample",
                )
        self.assertGreaterEqual(p_value[0], self.p_value_threshold)
        self.assertGreaterEqual(p_value[1], self.p_value_threshold)

    def test_impute_with_unphased_pedigree_control_notilda(self):
        command = ["impute.py",
                   "-c",
                   "--ibd", f"{tests_root}/test_data/sample.our",
                   "--bed", f"{tests_root}/test_data/sample_reduced1",
                   "--pedigree", f"{tests_root}/test_data/sample.ped",
                   "--out", f"{output_root}/test_impute_with_unphased_pedigree_control_notilda1",
                   "--threads", "10",
                   ]
        subprocess.check_call(command)
        coef, z, p_value = imputation_test([1],
                imputed_prefix = f"{output_root}/test_impute_with_unphased_pedigree_control_notilda",
                expected_prefix = f"{tests_root}/test_data/sample",
                )
        self.assertGreaterEqual(p_value[0], self.p_value_threshold)
        self.assertGreaterEqual(p_value[1], self.p_value_threshold)