from snipar.imputation.impute_from_sibs import *
from snipar.imputation.impute_from_sibs cimport *
import numpy as np
import unittest
from snipar.config import nan_integer
from snipar.tests.utils import *

class TestSibImpute(SniparTest):
    def test_get_IBD_type(self):
        #check whats the condition of start and end
        ibd1 = (10, 20)
        ibd0 = (25, 40)
        ibd2 = (50, 60)
        ibd_dict = {
            (b"jack", b"jim"):[ibd1[0], ibd1[1], 1, ibd0[0], ibd0[1], 0, ibd2[0], ibd2[1], 2],
        }

        inferred_ibd0 = get_IBD_type(b"jack", b"another thing", 1, ibd_dict)
        self.assertEqual(inferred_ibd0, nan_integer, msg="error when ids are not in the dict")
        inferred_ibd0 = get_IBD_type(b"another thing", b"jack", 1, ibd_dict)
        self.assertEqual(inferred_ibd0, nan_integer, msg="error when ids are not in the dict")
        inferred_ibd0 = get_IBD_type(b"another thing", b"another thing", 1, ibd_dict)
        self.assertEqual(inferred_ibd0, nan_integer, msg="error when ids are not in the dict")

        for i in range(60):
            inferred_ibd1 = get_IBD_type(b"jack", b"jim", i, ibd_dict)
            inferred_ibd2 = get_IBD_type(b"jim", b"jack", i, ibd_dict)
            self.assertEqual(inferred_ibd1, inferred_ibd2, msg="different id order changes the result")
            if ibd0[0] <= i <= ibd0[1]:
                self.assertEqual(inferred_ibd1, 0, msg="inferred IBD is not 0")
            elif ibd1[0] <= i <= ibd1[1]:
                self.assertEqual(inferred_ibd1, 1, msg="inferred IBD is not 1")
            elif ibd2[0] <= i <= ibd2[1]:
                self.assertEqual(inferred_ibd1, 2, msg="inferred IBD is not 2")
            else:
                self.assertEqual(inferred_ibd1, nan_integer, msg="inferred IBD is not 0")

    def test_dict_to_cmap(self):
        the_dict = {
            ("A","B"):[1,2,3,4],
            ("C","D"):[1,2,3,4],
            ("E","F"):[1,2,3,4],
            ("G","H"):[1,2,3,4],
        }
        expected_result = {
            (b"A",b"B"):[1,2,3,4],
            (b"C",b"D"):[1,2,3,4],
            (b"E",b"F"):[1,2,3,4],
            (b"G",b"H"):[1,2,3,4],
        }
        result = dict_to_cmap(the_dict)
        self.assertEqual(expected_result, result, msg="dict translation is not working")

    def test_impute_snp_from_offsprings_unphased(self):
        bed = np.array([[0],[1],[2]]).astype("b")
        sib_indexes = np.array([0, 1]).astype("i")
        snp_ibd0 = np.ones((10,2)).astype("i")
        snp_ibd1 = np.ones((10,2)).astype("i")
        snp_ibd2 = np.ones((10,2)).astype("i")
        snp = 0
        f = 0.1
        cdef cpair[double, bint] t
        parent_genotype_prob = np.array([(1-f)*(1-f), 2*f*(1-f), f*f])
        for i in range(3):
            for j in range(3):
                for count in range(10):
                    bed[0, snp] = i
                    bed[1, snp] = j
                    snp_ibd0[count] = [0, 1]
                    t = impute_snp_from_offsprings(snp, sib_indexes, 2, snp_ibd0, snp_ibd1, snp_ibd2, f, parent_genotype_prob, None, bed, None, count+1, 0, 0, False)
                    result, is_backup = t.first, t.second
                    sibsum = bed[snp_ibd0[0,0], snp] + bed[snp_ibd0[0,1], snp]
                    expected = sibsum/2
                    error_message = f"problem with type0, with sibs = {[i,j]}: result, expected = {(result, expected)}, is_backup={is_backup}"
                    self.assertAlmostEqual(result, expected, 4, msg = error_message)
                    self.assertTrue(not is_backup, error_message)

        for i in range(3):
            for j in range(3):
                for count in range(10):
                    bed[0, snp] = i
                    bed[1, snp] = j
                    snp_ibd1[count] = [0, 1]
                    t = impute_snp_from_offsprings(snp, sib_indexes, 2, snp_ibd0, snp_ibd1, snp_ibd2, f, parent_genotype_prob, None, bed, None, 0, count+1, 1, False)
                    result, is_backup = t.first, t.second
                    sibsum = bed[snp_ibd1[0,0], snp] + bed[snp_ibd1[0,1], snp]
                    expected_results = [f, 1+f, 1+2*f, 2+f, 3+f]
                    expected = expected_results[int(sibsum)]/2
                    error_message = f"problem with type1, with sibs = {[i,j]}: result, expected = {(result, expected)}, is_backup={is_backup}"
                    self.assertTrue(not is_backup, error_message)
                    if abs(i-j) == 2:
                        self.assertTrue(np.isnan(result), error_message)
                    else:
                        self.assertAlmostEqual(result, expected, 4, msg = error_message)


        for i in range(3):
            for j in range(3):
                for count in range(10):
                    bed[0, snp] = i
                    bed[1, snp] = j
                    snp_ibd2[count] = [0, 1]
                    t = impute_snp_from_offsprings(snp, sib_indexes,  2, snp_ibd0, snp_ibd1, snp_ibd2, f, parent_genotype_prob, None, bed, None, 0, 0, count+1, False)
                    result, is_backup = t.first, t.second
                    sibsum = bed[snp_ibd2[0,0], snp] + bed[snp_ibd2[0,1], snp]
                    expected = sibsum/4+f
                    error_message = f"problem with type2, with sibs = {[i,j]}: result, expected = {(result, expected)}, is_backup={is_backup}"
                    self.assertTrue(not is_backup, error_message)
                    #TODO add tests for backup imputation
                    if i!=j:
                        self.assertTrue(np.isnan(result), error_message)
                    else:
                        self.assertAlmostEqual(result, expected, 4, msg = error_message)

    def test_impute_snp_from_parent_offsprings_unphased(self):
        bed = np.array([[0],[1],[2]]).astype("b")
        sib_indexes = np.array([0, 1]).astype("i")
        snp_ibd0 = np.ones((10,2)).astype("i")
        snp_ibd1 = np.ones((10,2)).astype("i")
        snp_ibd2 = np.ones((10,2)).astype("i")
        snp = 0
        f = 0.1
        parent_genotype_prob = np.array([(1-f)*(1-f), 2*f*(1-f), f*f])
                  
        #first key is parent, second is sibs
        expected_result_IBD1 = [
            {
                (0,0) : 0.05263157894736842,
                (0,1) : 1,
                (1,0) : 1,
                (1,1) : 1.1818181818181819,
            },
            {
                (0,0):0,
                (0,1):0.18181818181818182,
                (1,0):0.18181818181818182,
                (1,1):0.024390243902439022,
                (1,2):1.0526315789473684,
                (2,1):1.0526315789473684,
                (2,2):2,

            },
            {
                (1,1):0.05263157894736841,
                (1,2):1,
                (2,1):1,
                (2,2):1.1818181818181819,
            },
        ]

        expected_result_IBD2 = [
            {
                (0,0):0.1,
                (1,1):1.1,
            },
            {
                (0,0):0.1,
                (1,1):0.2,
                (2,2):1.1,
            },
            {
                (1,1):0.1,
                (2,2):1.1,
            },
            
        ]

        for i in range(3):
            for j in range(3):
                for count in range(10):
                    for par in range(3):
                        bed[0, snp] = i
                        bed[1, snp] = j
                        bed[2, snp] = par
                        snp_ibd0[count] = [0, 1]
                        t = impute_snp_from_parent_offsprings(snp, 2, sib_indexes, 2, snp_ibd0, snp_ibd1, snp_ibd2, f, parent_genotype_prob, None, bed, None, None, count+1, 1, 1, False)
                        result, data = t.first, t.second
                        mendelian_error_count = data.first
                        is_backup = data.second
                        sibsum = bed[snp_ibd0[0,0], snp] + bed[snp_ibd0[0,1], snp]
                        expected = sibsum - par
                        error_message = f"problem with type0, with parent = {par}, and sibs = {[i,j]}: result, expected = {(result, expected)}, mendelian_erros={mendelian_error_count}, is_backup={is_backup}"
                        self.assertTrue(not is_backup, error_message)
                        if abs(par - i)>1 or abs(par - j)>1:                            
                            self.assertTrue(np.isnan(result), error_message)
                            self.assertTrue(mendelian_error_count>0, error_message)
                        elif expected > 2 or expected <0:
                            self.assertTrue(mendelian_error_count==0, error_message)
                            self.assertTrue(np.isnan(result), error_message)
                        else:                            
                            self.assertTrue(mendelian_error_count==0, error_message)
                            if np.isnan(expected):
                                self.assertTrue(np.isnan(result), msg = error_message)
                            else:
                                self.assertAlmostEqual(result, expected, 4, msg = error_message)

        for i in range(3):
            for j in range(3):
                for count in range(10):
                    for par in range(3):
                        bed[0, snp] = i
                        bed[1, snp] = j
                        bed[2, snp] = par
                        snp_ibd1[count] = [0, 1]
                        t = impute_snp_from_parent_offsprings(snp,  2, sib_indexes, 2, snp_ibd0, snp_ibd1, snp_ibd2, f, parent_genotype_prob, None, bed, None, None, 0, count+1, 1, False)
                        result, data = t.first, t.second
                        mendelian_error_count = data.first
                        is_backup = data.second
                        expected = expected_result_IBD1[par].get((i, j), np.nan)
                        error_message = f"problem with type1, with parent = {par}, and sibs = {[i,j]}: result, expected = {(result, expected)}, mendelian_erros={mendelian_error_count}, is_backup={is_backup}"
                        self.assertTrue(not is_backup, error_message)
                        if abs(par - i)>1 or abs(par - j)>1:                            
                            self.assertTrue(np.isnan(result), error_message)
                            self.assertTrue(mendelian_error_count>0, error_message)                            
                        elif (j==0 and i==2) or (j==2 and i==0):
                            self.assertTrue(mendelian_error_count==0, error_message)
                            self.assertTrue(np.isnan(result), error_message)
                        else:
                            self.assertTrue(mendelian_error_count==0, error_message)
                            if np.isnan(expected):
                                self.assertTrue(np.isnan(result), msg = error_message)
                            else:
                                self.assertAlmostEqual(result, expected, 4, msg = error_message)                            

        for i in range(3):
            for j in range(3):
                for count in range(10):
                    for par in range(3):
                        bed[0, snp] = i
                        bed[1, snp] = j
                        bed[2, snp] = par
                        snp_ibd2[count] = [0, 1]
                        t = impute_snp_from_parent_offsprings(snp, 2, sib_indexes, 2, snp_ibd0, snp_ibd1, snp_ibd2, f, parent_genotype_prob, None, bed, None, None, 0, 0, count+1, False)
                        result, data = t.first, t.second
                        mendelian_error_count = data.first
                        is_backup = data.second
                        expected = expected_result_IBD2[par].get((i, j), np.nan)
                        error_message = f"problem with type2, with parent = {par}, and sibs = {[i,j]}: result, expected = {(result, expected)}, mendelian_erros={mendelian_error_count}, is_backup={is_backup}"
                        self.assertTrue(not is_backup, error_message)
                        if abs(par - i)>1 or abs(par - j)>1:                            
                            self.assertTrue(np.isnan(result), error_message)
                            self.assertTrue(mendelian_error_count>0, error_message)
                        elif (j!=i):
                            self.assertTrue(mendelian_error_count==0, error_message)
                            self.assertTrue(np.isnan(result), error_message)
                        else:                            
                            self.assertTrue(mendelian_error_count==0, error_message)
                            if np.isnan(expected):
                                self.assertTrue(np.isnan(result), msg = error_message)
                            else:
                                self.assertAlmostEqual(result, expected, 4, msg = error_message)
    def test_get_IBD(self):
        length = 1000
        half_window = 100
        hap1 = np.array([1 for i in range(1000)]).astype("b")
        hap2 = np.array([1 for i in range(1000)]).astype("b")
        agreement_count = np.array([0 for i in range(1000)]).astype("i")
        agreement_percentage = np.array([0. for i in range(1000)])
        agreement = np.array([0 for i in range(1000)]).astype("i")
        get_IBD(hap1, hap2, length, half_window, 0.5, agreement_count, agreement_percentage, agreement)

        for i in range(length):
            self.assertEqual(agreement[i], 1)
        for i in range(length):
            self.assertAlmostEqual(agreement_percentage[i], 1.)

        hap1[[i%2==0 for i in range(length)]] = 0
        hap2[[i%2==0 for i in range(length)]] = 0
        get_IBD(hap1, hap2, length, half_window, 0.5, agreement_count, agreement_percentage, agreement)
        for i in range(length):
            self.assertEqual(agreement[i], 1)
        for i in range(length):
            self.assertAlmostEqual(agreement_percentage[i], 1.)

        hap1 = np.array([i//500 for i in range(1000)]).astype("b")
        hap2 = np.array([i//500 for i in range(1000)]).astype("b")
        get_IBD(hap1, hap2, length, half_window, 0.5, agreement_count, agreement_percentage, agreement)
        for i in range(length):
            self.assertEqual(agreement[i], 1)
        
        hap1 = np.array([i//500 for i in range(1000)]).astype("b")
        hap2 = np.array([1 for i in range(1000)]).astype("b")
        get_IBD(hap1, hap2, length, half_window, 0.9999, agreement_count, agreement_percentage, agreement)

        for i in range(500+half_window):
            self.assertEqual(agreement[i], 0)
        for i in range(500+half_window, length):
            self.assertEqual(agreement[i], 1)
