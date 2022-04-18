from snipar.scripts import impute
from snipar.tests.test_imputation import imputation_test
from snipar.tests.utils import *
#TODO add tests with nan

class TestImpute(SniparTest):

    def test_impute_with_unphased_pedigree_control(self):
        command = [
                   "-c",
                   "--ibd", f"{tests_root}/test_data/sample.our",                   
                   "--bed", f"{tests_root}/test_data/sample_reduced@",
                   "--chr_range", "1-2",
                   "--pedigree", f"{tests_root}/test_data/sample.ped",
                   "--threads", "2",
                   "--out", f"{output_root}/test_impute_with_unphased_pedigree_control@",
                   ]
        if not self.log:
            command = ["-silent_progress"] + command
        args=impute.parser.parse_args(command)
        impute.main(args)
        coef, z, p_value = imputation_test([1],
                imputed_prefix = f"{output_root}/test_impute_with_unphased_pedigree_control",
                expected_prefix = f"{tests_root}/test_data/sample_reduced",
                )
        self.assertGreaterEqual(p_value[0], self.p_value_threshold)
        self.assertGreaterEqual(p_value[1], self.p_value_threshold)

    def test_impute_with_unphased_pedigree_control_backup(self):
        command = [
                   "-c",
                   "--start", "0",
                   "--end", f"{self.subsample_snp}",
                   "--ibd", f"{tests_root}/test_data/sample.empty",                   
                   "--bed", f"{tests_root}/test_data/sample_reduced@",
                   "--chr_range", "1-2",
                   "--pedigree", f"{tests_root}/test_data/sample.ped",
                   "--threads", "2",
                   "--out", f"{output_root}/test_impute_with_unphased_pedigree_control_backup@",
                   ]
        if not self.log:
            command = ["-silent_progress"] + command
        args=impute.parser.parse_args(command)
        impute.main(args)
        coef, z, p_value = imputation_test([1],
                imputed_prefix = f"{output_root}/test_impute_with_unphased_pedigree_control_backup",
                expected_prefix = f"{tests_root}/test_data/sample",
                start = 0,
                end = self.subsample_snp,
                )
        #TODO fix this
        # self.assertGreaterEqual(p_value[0], self.p_value_threshold)
        # self.assertGreaterEqual(p_value[1], self.p_value_threshold)
    
    def test_impute_with_unphased_pedigree_wild_ibd_control(self):
        command = [
                   "-c",
                   "--start", "0",
                   "--end", f"{self.subsample_snp}",
                   "--ibd", f"{tests_root}/test_data/sample@.our",
                   "--bed", f"{tests_root}/test_data/sample_reduced@",
                   "--chr_range", "1-2",
                   "--pedigree", f"{tests_root}/test_data/sample.ped",
                   "--threads", "2",
                   "--out", f"{output_root}/test_impute_with_unphased_pedigree_control@",
                   ]
        if not self.log:
            command = ["-silent_progress"] + command
        args=impute.parser.parse_args(command)
        impute.main(args)
        coef, z, p_value = imputation_test([1],
                imputed_prefix = f"{output_root}/test_impute_with_unphased_pedigree_control",
                expected_prefix = f"{tests_root}/test_data/sample",
                start = 0,
                end = self.subsample_snp,
                )
        self.assertGreaterEqual(p_value[0], self.p_value_threshold)
        self.assertGreaterEqual(p_value[1], self.p_value_threshold)

    def test_impute_with_unphased_pedigree(self):
        command = [
                   "--start", "0",
                   "--end", f"{self.subsample_snp}",
                   "--ibd", f"{tests_root}/test_data/sample.our",
                   "--bed", f"{tests_root}/test_data/sample_reduced@",
                   "--chr_range", "1",
                   "--pedigree", f"{tests_root}/test_data/sample.ped",
                   "--threads", "2",
                   "--out", f"{output_root}/test_impute_with_unphased_pedigree@",
                   ]
        if not self.log:
            command = ["-silent_progress"] + command
        args=impute.parser.parse_args(command)
        impute.main(args)

    def test_impute_with_unphased_pedigree_chunks_control(self):
        command = [
                   "-c",
                   "--start", "0",
                   "--end", f"{self.subsample_snp}",
                   "--ibd", f"{tests_root}/test_data/sample.our",
                   "--bed", f"{tests_root}/test_data/sample_reduced@",
                   "--chr_range", "1-2",
                   "--pedigree", f"{tests_root}/test_data/sample.ped",
                   "--chunks", "3",
                   "--threads", "2",
                   "--out", f"{output_root}/test_impute_with_unphased_pedigree_control@",
                   ]
        if not self.log:
            command = ["-silent_progress"] + command
        args=impute.parser.parse_args(command)
        impute.main(args)
        coef, z, p_value = imputation_test([1],
                imputed_prefix = f"{output_root}/test_impute_with_unphased_pedigree_control",
                expected_prefix = f"{tests_root}/test_data/sample",
                start = 0,
                end = self.subsample_snp,
                )
        self.assertGreaterEqual(p_value[0], self.p_value_threshold)
        self.assertGreaterEqual(p_value[1], self.p_value_threshold)


    def test_impute_with_unphased_pedigree_control_legacy_ibd(self):
        command = [
                   "-c",
                   "--start", "0",
                   "--end", f"{self.subsample_snp}",
                   "--ibd", f"{tests_root}/test_data/sample.king",
                   "--ibd_is_king",
                   "--bed", f"{tests_root}/test_data/sample_reduced@",
                   "--chr_range", "1-2",
                   "--pedigree", f"{tests_root}/test_data/sample.ped",
                   "--threads", "2",
                   "--out", f"{output_root}/test_impute_with_unphased_pedigree_control@",
                   ]
        if not self.log:
            command = ["-silent_progress"] + command
        args=impute.parser.parse_args(command)
        impute.main(args)
        coef, z, p_value = imputation_test([1],
                imputed_prefix = f"{output_root}/test_impute_with_unphased_pedigree_control",
                expected_prefix = f"{tests_root}/test_data/sample",
                start = 0,
                end = self.subsample_snp,
                )
        self.assertGreaterEqual(p_value[0], self.p_value_threshold)
        self.assertGreaterEqual(p_value[1], self.p_value_threshold)


    def test_impute_with_phased_pedigree_control(self):
        command = [
                   "-c",
                   "--start", "0",
                   "--end", f"{self.subsample_snp}",
                   "--ibd", f"{tests_root}/test_data/sample.our",
                   "--bgen", f"{tests_root}/test_data/sample_reduced@",
                   "--chr_range", "1-2",
                   "--pedigree", f"{tests_root}/test_data/sample.ped",
                   "--threads", "2",
                   "--out", f"{output_root}/test_impute_with_phased_pedigree_control@",
                   ]
        if not self.log:
            command = ["-silent_progress"] + command
        args=impute.parser.parse_args(command)
        impute.main(args)
        coef, z, p_value = imputation_test([1, 2],
                imputed_prefix = f"{output_root}/test_impute_with_phased_pedigree_control",
                expected_prefix = f"{tests_root}/test_data/sample",
                start = 0,
                end = self.subsample_snp,
                )
        self.assertGreaterEqual(p_value[0], self.p_value_threshold)
        self.assertGreaterEqual(p_value[1], self.p_value_threshold)        

    def test_impute_with_phased_pedigree_chunks_control(self):
        command = [
                   "-c",
                   "--start", "0",
                   "--end", f"{self.subsample_snp}",
                   "--ibd", f"{tests_root}/test_data/sample.our",
                   "--bgen", f"{tests_root}/test_data/sample_reduced@",
                   "--chr_range", "1-2",
                   "--chunks", "3",
                   "--pedigree", f"{tests_root}/test_data/sample.ped",
                   "--threads", "2",
                   "--out", f"{output_root}/test_impute_with_phased_pedigree_control@",
                   ]
        if not self.log:
            command = ["-silent_progress"] + command
        args=impute.parser.parse_args(command)
        impute.main(args)
        coef, z, p_value = imputation_test([1, 2],
                imputed_prefix = f"{output_root}/test_impute_with_phased_pedigree_control",
                expected_prefix = f"{tests_root}/test_data/sample",
                start = 0,
                end = self.subsample_snp,
                )
        self.assertGreaterEqual(p_value[0], self.p_value_threshold)
        self.assertGreaterEqual(p_value[1], self.p_value_threshold)

    def test_impute_with_unphased_king(self):
        command = ["--start", "0",
                   "--end", f"{self.subsample_snp}",
                   "--ibd", f"{tests_root}/test_data/sample.our",
                   "--bed", f"{tests_root}/test_data/sample_reduced@",
                   "--chr_range", "1-2",
                   "--king", f"{tests_root}/test_data/sample.king",
                   "--agesex", f"{tests_root}/test_data/sample.agesex",
                   "--threads", "2",
                   "--out", f"{output_root}/test_impute_with_unphased_king@",
                   ]
        if not self.log:
            command = ["-silent_progress"] + command
        args=impute.parser.parse_args(command)
        impute.main(args)

    def test_impute_with_unphased_king_control(self):
        command = ["-c",
                   "--start", "0",
                   "--end", f"{self.subsample_snp}",
                   "--ibd", f"{tests_root}/test_data/sample.our",
                   "--bed", f"{tests_root}/test_data/sample_reduced@",
                   "--chr_range", "1-2",
                   "--king", f"{tests_root}/test_data/sample.king",
                   "--agesex", f"{tests_root}/test_data/sample.agesex",
                   "--threads", "2",
                   "--out", f"{output_root}/test_impute_with_unphased_king_control@",
                   ]
        if not self.log:
            command = ["-silent_progress"] + command
        args=impute.parser.parse_args(command)
        impute.main(args)
        coef, z, p_value = imputation_test([1, 2],
                imputed_prefix = f"{output_root}/test_impute_with_unphased_king_control",
                expected_prefix = f"{tests_root}/test_data/sample",
                start = 0,
                end = self.subsample_snp,
                )
        self.assertGreaterEqual(p_value[0], self.p_value_threshold)
        self.assertGreaterEqual(p_value[1], self.p_value_threshold)

    def test_impute_with_unphased_king_control_pca(self):
        command = [
                   "-c",
                   "--start", "0",
                   "--end", f"{self.subsample_snp}",
                   "--ibd", f"{tests_root}/test_data/sample.our",
                   "--bed", f"{tests_root}/test_data/sample@",
                   "--chr_range", "1-2",
                   "--king", f"{tests_root}/test_data/sample.king",
                   "--agesex", f"{tests_root}/test_data/sample.agesex",
                   "--pcs", f"{tests_root}/test_data/sample1_2_pca.eigenvec",
                   "--pc_num", "8",
                   "-find_optimal_pc",
                   "--threads", "2",
                   "--out", f"{output_root}/test_impute_with_unphased_king_control@",
                   ]
        if not self.log:
            command = ["-silent_progress"] + command
        args=impute.parser.parse_args(command)
        impute.main(args)
        coef, z, p_value = imputation_test([1, 2],
                imputed_prefix = f"{output_root}/test_impute_with_unphased_king_control",
                expected_prefix = f"{tests_root}/test_data/sample",
                start = 0,
                end = self.subsample_snp,
                )
        #TODO use pvals
        print(coef)
        self.assertAlmostEqual(coef[0], 1., delta=0.02)
        self.assertAlmostEqual(coef[1], 1., delta=0.02)

    def test_impute_with_unphased_king_control_legacy_ibd(self):
        command = [
                   "-c",
                   "--start", "0",
                   "--end", f"{self.subsample_snp}",
                   "--ibd", f"{tests_root}/test_data/sample.king",
                   "--ibd_is_king",
                   "--bed", f"{tests_root}/test_data/sample_reduced@",
                   "--chr_range", "1-2",
                   "--king", f"{tests_root}/test_data/sample.king",
                   "--agesex", f"{tests_root}/test_data/sample.agesex",
                   "--threads", "2",
                   "--out", f"{output_root}/test_impute_with_unphased_king_control@",
                   ]
        if not self.log:
            command = ["-silent_progress"] + command
        args=impute.parser.parse_args(command)
        impute.main(args)
        coef, z, p_value = imputation_test([1, 2],
                imputed_prefix = f"{output_root}/test_impute_with_unphased_king_control",
                expected_prefix = f"{tests_root}/test_data/sample",
                start = 0,
                end = self.subsample_snp,
                )
        self.assertGreaterEqual(p_value[0], self.p_value_threshold)
        self.assertGreaterEqual(p_value[1], self.p_value_threshold)

    def test_impute_with_unphased_pedigree_control_multiprocess(self):
        command = [
                   "-c",
                   "--start", "0",
                   "--end", f"{self.subsample_snp}",
                   "--ibd", f"{tests_root}/test_data/sample.our",
                   "--bed", f"{tests_root}/test_data/sample_reduced@",
                   "--chr_range", "1-2",
                   "--pedigree", f"{tests_root}/test_data/sample.ped",
                #    "--threads", "5", TODO figure out why interaction goes wrong
                   "--out", f"{output_root}/test_impute_with_unphased_pedigree_control_multiprocess@",
                   "--processes", "2",
                   ]
        if not self.log:
            command = ["-silent_progress"] + command
        args=impute.parser.parse_args(command)
        impute.main(args)
        coef, z, p_value = imputation_test([1, 2],
                imputed_prefix = f"{output_root}/test_impute_with_unphased_pedigree_control_multiprocess",
                expected_prefix = f"{tests_root}/test_data/sample",
                start = 0,
                end = self.subsample_snp,
                )
        self.assertGreaterEqual(p_value[0], self.p_value_threshold)
        self.assertGreaterEqual(p_value[1], self.p_value_threshold)

    def test_impute_with_unphased_king_control_legacy_tilda_ibd(self):
        command = [
                   "-c",
                   "--start", "0",
                   "--end", f"{self.subsample_snp}",
                   "--ibd", f"{tests_root}/test_data/sample@.king",
                   "--ibd_is_king",
                   "--bed", f"{tests_root}/test_data/sample_reduced@",
                   "--chr_range", "1-2",
                   "--king", f"{tests_root}/test_data/sample.king",
                   "--agesex", f"{tests_root}/test_data/sample.agesex",
                   "--threads", "2",
                   "--out", f"{output_root}/test_impute_with_unphased_king_control@",
                   ]
        if not self.log:
            command = ["-silent_progress"] + command
        args=impute.parser.parse_args(command)
        impute.main(args)
        coef, z, p_value = imputation_test([1, 2],
                imputed_prefix = f"{output_root}/test_impute_with_unphased_king_control",
                expected_prefix = f"{tests_root}/test_data/sample",
                start = 0,
                end = self.subsample_snp,
                )
        self.assertGreaterEqual(p_value[0], self.p_value_threshold)
        self.assertGreaterEqual(p_value[1], self.p_value_threshold)

#TODO test without providing ibd, that should result in all backup
    def test_impute_with_unphased_pedigree_control_nothread(self):
        command = [
                   "-c",
                   "--start", "0",
                   "--end", f"{self.subsample_snp}",
                   "--ibd", f"{tests_root}/test_data/sample.our",
                   "--bed", f"{tests_root}/test_data/sample_reduced@",
                   "--chr_range", "1-2",
                   "--pedigree", f"{tests_root}/test_data/sample.ped",
                   "--out", f"{output_root}/test_impute_with_unphased_pedigree_control_multithread@",
                   ]
        if not self.log:
            command = ["-silent_progress"] + command
        args=impute.parser.parse_args(command)
        impute.main(args)
        coef, z, p_value = imputation_test([1, 2],
                imputed_prefix = f"{output_root}/test_impute_with_unphased_pedigree_control_multithread",
                expected_prefix = f"{tests_root}/test_data/sample",
                start = 0,
                end = self.subsample_snp,
                )
        self.assertGreaterEqual(p_value[0], self.p_value_threshold)
        self.assertGreaterEqual(p_value[1], self.p_value_threshold)

    def test_impute_with_unphased_pedigree_control_notilda(self):
        command = [
                   "-c",
                   "--start", "0",
                   "--end", f"{self.subsample_snp}",
                   "--ibd", f"{tests_root}/test_data/sample.our",
                   "--bed", f"{tests_root}/test_data/sample_reduced1",
                   "--pedigree", f"{tests_root}/test_data/sample.ped",
                   "--out", f"{output_root}/test_impute_with_unphased_pedigree_control_notilda1",
                   "--threads", "2",
                   ]
        if not self.log:
            command = ["-silent_progress"] + command
        args=impute.parser.parse_args(command)
        impute.main(args)
        coef, z, p_value = imputation_test([1],
                imputed_prefix = f"{output_root}/test_impute_with_unphased_pedigree_control_notilda",
                expected_prefix = f"{tests_root}/test_data/sample",
                start = 0,
                end = self.subsample_snp,
                )
        self.assertGreaterEqual(p_value[0], self.p_value_threshold)
        self.assertGreaterEqual(p_value[1], self.p_value_threshold)
