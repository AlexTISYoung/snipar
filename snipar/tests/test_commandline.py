from snipar.scripts import impute
from snipar.tests.test_imputation import imputation_test
import os
from snipar.tests.utils import *
#TODO add tests with nan

class TestCommanline(SniparTest):
    p_value_threshold = 0.01
    subsample_snp = 50
    @classmethod
    def setUpClass(cls):        
        logging.basicConfig(level=logging.ERROR, format='%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s')
        warnings.filterwarnings("ignore")

    @classmethod
    def tearDownClass(cls):
        logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s')

    def setUp(self):
        self.startTime = time.time()

    def tearDown(self):
        t = time.time() - self.startTime
        print('%s: %.3f' % (self.id(), t))

    def test_impute_with_unphased_pedigree_control(self):
        command = [
                   "-silent_progress",
                   "-c",
                   "--ibd", f"{tests_root}/test_data/sample.our",                   
                   "--bed", f"{tests_root}/test_data/sample_reduced~",
                   "--from_chr", "1",
                   "--to_chr", "3",
                   "--pedigree", f"{tests_root}/test_data/sample.ped",
                   "--threads", "2",
                   "--out", f"{output_root}/test_impute_with_unphased_pedigree_control~",
                   ]
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
                   "-silent_progress",
                   "-c",
                   "--start", "0",
                   "--end", f"{self.subsample_snp}",
                   "--ibd", f"{tests_root}/test_data/sample.empty",                   
                   "--bed", f"{tests_root}/test_data/sample_reduced~",
                   "--from_chr", "1",
                   "--to_chr", "3",
                   "--pedigree", f"{tests_root}/test_data/sample.ped",
                   "--threads", "2",
                   "--out", f"{output_root}/test_impute_with_unphased_pedigree_control_backup~",
                   ]
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
                   "-silent_progress",
                   "-c",
                   "--start", "0",
                   "--end", f"{self.subsample_snp}",
                   "--ibd", f"{tests_root}/test_data/sample~.our",
                   "--bed", f"{tests_root}/test_data/sample_reduced~",
                   "--from_chr", "1",
                   "--to_chr", "3",
                   "--pedigree", f"{tests_root}/test_data/sample.ped",
                   "--threads", "2",
                   "--out", f"{output_root}/test_impute_with_unphased_pedigree_control~",
                   ]
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
                   "-silent_progress",
                   "--start", "0",
                   "--end", f"{self.subsample_snp}",
                   "--ibd", f"{tests_root}/test_data/sample.our",
                   "--bed", f"{tests_root}/test_data/sample_reduced~",
                   "--from_chr", "1",
                   "--to_chr", "2",
                   "--pedigree", f"{tests_root}/test_data/sample.ped",
                   "--threads", "2",
                   "--out", f"{output_root}/test_impute_with_unphased_pedigree~",
                   ]
        args=impute.parser.parse_args(command)
        impute.main(args)

    def test_impute_with_unphased_pedigree_chunks_control(self):
        command = [
                   "-silent_progress",
                   "-c",
                   "--start", "0",
                   "--end", f"{self.subsample_snp}",
                   "--ibd", f"{tests_root}/test_data/sample.our",
                   "--bed", f"{tests_root}/test_data/sample_reduced~",
                   "--from_chr", "1",
                   "--to_chr", "3",
                   "--pedigree", f"{tests_root}/test_data/sample.ped",
                   "--chunks", "3",
                   "--threads", "2",
                   "--out", f"{output_root}/test_impute_with_unphased_pedigree_control~",
                   ]
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
                   "-silent_progress",
                   "-c",
                   "--start", "0",
                   "--end", f"{self.subsample_snp}",
                   "--ibd", f"{tests_root}/test_data/sample.king",
                   "--ibd_is_king",
                   "--bed", f"{tests_root}/test_data/sample_reduced~",
                   "--from_chr", "1",
                   "--to_chr", "3",
                   "--pedigree", f"{tests_root}/test_data/sample.ped",
                   "--threads", "2",
                   "--out", f"{output_root}/test_impute_with_unphased_pedigree_control~",
                   ]
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
                   "-silent_progress",
                   "-c",
                   "--start", "0",
                   "--end", f"{self.subsample_snp}",
                   "--ibd", f"{tests_root}/test_data/sample.our",
                   "--bgen", f"{tests_root}/test_data/sample_reduced~",
                   "--from_chr", "1",
                   "--to_chr", "3",
                   "--pedigree", f"{tests_root}/test_data/sample.ped",
                   "--threads", "2",
                   "--out", f"{output_root}/test_impute_with_phased_pedigree_control~",
                   ]
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
                   "-silent_progress",
                   "-c",
                   "--start", "0",
                   "--end", f"{self.subsample_snp}",
                   "--ibd", f"{tests_root}/test_data/sample.our",
                   "--bgen", f"{tests_root}/test_data/sample_reduced~",
                   "--from_chr", "1",
                   "--to_chr", "3",
                   "--chunks", "3",
                   "--pedigree", f"{tests_root}/test_data/sample.ped",
                   "--threads", "2",
                   "--out", f"{output_root}/test_impute_with_phased_pedigree_control~",
                   ]
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
        command = [
                   "-silent_progress",
                   "--start", "0",
                   "--end", f"{self.subsample_snp}",
                   "--ibd", f"{tests_root}/test_data/sample.our",
                   "--bed", f"{tests_root}/test_data/sample_reduced~",
                   "--from_chr", "1",
                   "--to_chr", "3",
                   "--king", f"{tests_root}/test_data/sample.king",
                   "--agesex", f"{tests_root}/test_data/sample.agesex",
                   "--threads", "2",
                   "--out", f"{output_root}/test_impute_with_unphased_king~",
                   ]
        args=impute.parser.parse_args(command)
        impute.main(args)

    def test_impute_with_unphased_king_control(self):
        command = [
                   "-silent_progress",
                   "-c",
                   "--start", "0",
                   "--end", f"{self.subsample_snp}",
                   "--ibd", f"{tests_root}/test_data/sample.our",
                   "--bed", f"{tests_root}/test_data/sample_reduced~",
                   "--from_chr", "1",
                   "--to_chr", "3",
                   "--king", f"{tests_root}/test_data/sample.king",
                   "--agesex", f"{tests_root}/test_data/sample.agesex",
                   "--threads", "2",
                   "--out", f"{output_root}/test_impute_with_unphased_king_control~",
                   ]
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
                   "-silent_progress",
                   "-c",
                   "--start", "0",
                   "--end", f"{self.subsample_snp}",
                   "--ibd", f"{tests_root}/test_data/sample.our",
                   "--bed", f"{tests_root}/test_data/sample~",
                   "--from_chr", "1",
                   "--to_chr", "3",
                   "--king", f"{tests_root}/test_data/sample.king",
                   "--agesex", f"{tests_root}/test_data/sample.agesex",
                   "--pcs", f"{tests_root}/test_data/sample1_2_pca.eigenvec",
                   "--pc_num", "8",
                   "-find_optimal_pc",
                   "--threads", "2",
                   "--out", f"{output_root}/test_impute_with_unphased_king_control~",
                   ]
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
                   "-silent_progress",
                   "-c",
                   "--start", "0",
                   "--end", f"{self.subsample_snp}",
                   "--ibd", f"{tests_root}/test_data/sample.king",
                   "--ibd_is_king",
                   "--bed", f"{tests_root}/test_data/sample_reduced~",
                   "--from_chr", "1",
                   "--to_chr", "3",
                   "--king", f"{tests_root}/test_data/sample.king",
                   "--agesex", f"{tests_root}/test_data/sample.agesex",
                   "--threads", "2",
                   "--out", f"{output_root}/test_impute_with_unphased_king_control~",
                   ]
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
                   "-silent_progress",   
                   "-c",
                   "--start", "0",
                   "--end", f"{self.subsample_snp}",
                   "--ibd", f"{tests_root}/test_data/sample.our",
                   "--bed", f"{tests_root}/test_data/sample_reduced~",
                   "--from_chr", "1",
                   "--to_chr", "3",
                   "--pedigree", f"{tests_root}/test_data/sample.ped",
                #    "--threads", "5", TODO figure out why interaction goes wrong
                   "--out", f"{output_root}/test_impute_with_unphased_pedigree_control_multiprocess~",
                   "--processes", "2",
                   ]
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
                   "-silent_progress",
                   "-c",
                   "--start", "0",
                   "--end", f"{self.subsample_snp}",
                   "--ibd", f"{tests_root}/test_data/sample~.king",
                   "--ibd_is_king",
                   "--bed", f"{tests_root}/test_data/sample_reduced~",
                   "--from_chr", "1",
                   "--to_chr", "3",
                   "--king", f"{tests_root}/test_data/sample.king",
                   "--agesex", f"{tests_root}/test_data/sample.agesex",
                   "--threads", "2",
                   "--out", f"{output_root}/test_impute_with_unphased_king_control~",
                   ]
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
                   "-silent_progress",
                   "-c",
                   "--start", "0",
                   "--end", f"{self.subsample_snp}",
                   "--ibd", f"{tests_root}/test_data/sample.our",
                   "--bed", f"{tests_root}/test_data/sample_reduced~",
                   "--from_chr", "1",
                   "--to_chr", "3",
                   "--pedigree", f"{tests_root}/test_data/sample.ped",
                   "--out", f"{output_root}/test_impute_with_unphased_pedigree_control_multithread~",
                   ]
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
                   "-silent_progress",
                   "-c",
                   "--start", "0",
                   "--end", f"{self.subsample_snp}",
                   "--ibd", f"{tests_root}/test_data/sample.our",
                   "--bed", f"{tests_root}/test_data/sample_reduced1",
                   "--pedigree", f"{tests_root}/test_data/sample.ped",
                   "--out", f"{output_root}/test_impute_with_unphased_pedigree_control_notilda1",
                   "--threads", "2",
                   ]
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
