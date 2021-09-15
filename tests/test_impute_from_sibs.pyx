import sys
from sibreg.bin.impute_from_sibs import *
from sibreg.bin.impute_from_sibs cimport *
import numpy as np
import unittest
from config import nan_integer

class TestSibImpute(unittest.TestCase):
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
            if ibd0[0] <= i < ibd0[1]:
                self.assertEqual(inferred_ibd1, 0, msg="inferred IBD is not 0")
            elif ibd1[0] <= i < ibd1[1]:
                self.assertEqual(inferred_ibd1, 1, msg="inferred IBD is not 1")
            elif ibd2[0] <= i < ibd2[1]:
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
        snp_ibd0 = np.ones((10,2)).astype("i")
        snp_ibd1 = np.ones((10,2)).astype("i")
        snp_ibd2 = np.ones((10,2)).astype("i")
        snp = 0
        f = 0.1
        parent_genotype_prob = np.array([(1-f)*(1-f), 2*f*(1-f), f*f])
        for i in range(3):
            for j in range(3):
                for count in range(10):
                    snp_ibd0[count] = [i, j]
                    sib_indexes = np.array([0, 1, 2]).astype("i")
                    result = np.array(impute_snp_from_offsprings(snp, sib_indexes, 3, snp_ibd0, snp_ibd1, snp_ibd2, f, parent_genotype_prob, None, bed, None, count+1, 0, 0))
                    sibsum = bed[snp_ibd0[0,0], snp] + bed[snp_ibd0[0,1], snp]
                    self.assertAlmostEqual(result, sibsum/2, 4, msg = "problem with type0")

        for i in range(3):
            for j in range(3):
                for count in range(10):
                    snp_ibd1[count] = [i, j]
                    sib_indexes = np.array([0, 1, 2]).astype("i")
                    result = impute_snp_from_offsprings(snp, sib_indexes, 3, snp_ibd0, snp_ibd1, snp_ibd2, f, parent_genotype_prob, None, bed, None, 0, count+1, 1)
                    sibsum = bed[snp_ibd1[0,0], snp] + bed[snp_ibd1[0,1], snp]
                    expected_results = [f, 1+f, 1+2*f, 2+f, 3+f]
                    self.assertAlmostEqual(result, expected_results[int(sibsum)]/2, 4, msg = "problem with type1"+str((result, expected_results[int(sibsum)]))+str((i,j))+ str(sibsum))

        for i in range(3):
            for j in range(3):
                for count in range(10):
                    snp_ibd2[count] = [i, j]
                    sib_indexes = np.array([0, 1, 2]).astype("i")
                    result = impute_snp_from_offsprings(snp, sib_indexes,  3, snp_ibd0, snp_ibd1, snp_ibd2, f, parent_genotype_prob, None, bed, None, 0, 0, count+1)
                    sibsum = bed[snp_ibd2[0,0], snp] + bed[snp_ibd2[0,1], snp]
                    self.assertAlmostEqual(result, (sibsum/2. + 2*f)/2, 4, msg = "problem with type2")

    def test_impute_snp_from_parent_offsprings_unphased(self):
        bed = np.array([[0],[1],[2]]).astype("b")
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
                (1,1) : 1.1818181818181819,
            },
            {
                (0,0):0,
                (0,1):0.18181818181818182,
                (1,1):0.024390243902439022,
                (1,2):1.0526315789473684,
                (2,2):2,

            },
            {
                (1,1):0.05263157894736841,
                (1,2):1,
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
            for j in range(i,3):
                for count in range(10):
                    for par in range(3):
                        if par == 2 and (i == 0 or j == 0):
                            continue                        
                        if par == 0 and (i == 2 or j == 2):
                            continue
                        sib_indexes = np.array([0, 1, 2]).astype("i")
                        snp_ibd0[count] = [i, j]
                        result = impute_snp_from_parent_offsprings(snp, par, sib_indexes, 3, snp_ibd0, snp_ibd1, snp_ibd2, f, parent_genotype_prob, None, bed, None, None, count+1, 1, 1)
                        sibsum = bed[snp_ibd0[0,0], snp] + bed[snp_ibd0[0,1], snp]
                        expected = sibsum - par
                        if expected<0 or expected>2:
                            expected = np.nan
                        if not np.isnan(result) or not np.isnan(expected):
                            self.assertAlmostEqual(result, expected, 4, msg = "problem with type0, with parent = "+str(par)+", and sibs = "+str([i,j]) + " result, expected = "+str((result, expected)))

        for i in range(3):
            for j in range(i, 3):
                for count in range(10):
                    for par in range(3):
                        sib_indexes = np.array([0, 1, 2]).astype("i")
                        snp_ibd1[count] = [i, j]
                        result = impute_snp_from_parent_offsprings(snp,  par, sib_indexes, 3, snp_ibd0, snp_ibd1, snp_ibd2, f, parent_genotype_prob, None, bed, None, None, 0, count+1, 1)
                        expected = expected_result_IBD1[par].get((i, j), np.nan)
                        if not np.isnan(result) or not np.isnan(expected):
                            self.assertAlmostEqual(result, expected, 4, msg = "problem with type1, with parent = "+str(par)+", and sibs = "+str([i,j])+", expected,result = "+str((expected, result)))
                            

        for i in range(3):
            for j in range(i, 3):
                for count in range(10):
                    for par in range(3):
                        sib_indexes = np.array([0, 1, 2]).astype("i")
                        snp_ibd2[count] = [i, j]
                        result = impute_snp_from_parent_offsprings(snp, par, sib_indexes, 3, snp_ibd0, snp_ibd1, snp_ibd2, f, parent_genotype_prob, None, bed, None, None, 0, 0, count+1)
                        expected = expected_result_IBD2[par].get((i, j), np.nan)
                        if not np.isnan(result) or not np.isnan(expected):
                            self.assertAlmostEqual(result, expected, 4, msg = "problem with type2, with parent = "+str(par)+", and sibs = "+str([i,j])+", expected, result = "+str((expected, result)))
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
