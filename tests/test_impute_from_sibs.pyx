# distutils: language=c++
import sys
from sibreg.bin.impute_from_sibs import *
from sibreg.bin.impute_from_sibs cimport *
import unittest

class TestSibImpute(unittest.TestCase):
    def test_get_IBD_type(self):
        #check whats the condition of start and end
        ibd1 = (10, 20)
        ibd2 = (30, 40)
        ibd_dict = {
            ("jack", "jim"):[ibd1[0], ibd1[1], 1, ibd2[0], ibd2[1], 2],
        }

        inferred_ibd0 = get_IBD_type("jack", "another thing", 1, ibd_dict)
        self.assertEqual(inferred_ibd0, 0, msg="error when ids are not in the dict")
        inferred_ibd0 = get_IBD_type("another thing", "jack", 1, ibd_dict)
        self.assertEqual(inferred_ibd0, 0, msg="error when ids are not in the dict")
        inferred_ibd0 = get_IBD_type("another thing", "another thing", 1, ibd_dict)
        self.assertEqual(inferred_ibd0, 0, msg="error when ids are not in the dict")

        for i in range(50):
            inferred_ibd1 = get_IBD_type("jack", "jim", i, ibd_dict)
            inferred_ibd2 = get_IBD_type("jim", "jack", i, ibd_dict)
            self.assertEqual(inferred_ibd1, inferred_ibd2, msg="different id order changes the result")

            if ibd2[0] <= i <= ibd2[1]:
                self.assertEqual(inferred_ibd1, 2, msg="inferred IBD is not 2")
            elif ibd1[0] <= i <= ibd1[1]:
                self.assertEqual(inferred_ibd1, 1, msg="inferred IBD is not 1")
            else:
                self.assertEqual(inferred_ibd1, 0, msg="inferred IBD is not 0")

    def test_dict_to_cmap(self):
        the_dict = {
            ("A","B"):[1,2,3,4],
            ("C","D"):[1,2,3,4],
            ("E","F"):[1,2,3,4],
            ("G","H"):[1,2,3,4],
        }
        result = dict_to_cmap(the_dict)
        self.assertEqual(the_dict, result, msg="dict translation is not working")

    def test_impute_snp(self):
        bed = np.array([[0.],[1.],[2.]])
        snp_ibd0 = np.ones((10,2)).astype(int)
        snp_ibd1 = np.ones((10,2)).astype(int)
        snp_ibd2 = np.ones((10,2)).astype(int)
        snp = 0
        f = 0.1
        for i in range(3):
            for j in range(3):
                for count in range(10):
                    snp_ibd0[count] = [i, j]
                    result = impute_snp(snp, snp_ibd0, snp_ibd1, snp_ibd2, f, bed, count+1, 1, 1)
                    sibsum = bed[snp_ibd0[0,0], snp] + bed[snp_ibd0[0,1], snp]
                    self.assertAlmostEqual(result, sibsum, 3, msg = "problem with type0")

        for i in range(3):
            for j in range(3):
                for count in range(10):
                    snp_ibd1[count] = [i, j]
                    result = impute_snp(snp, snp_ibd0, snp_ibd1, snp_ibd2, f, bed, 0, count+1, 1)
                    sibsum = bed[snp_ibd1[0,0], snp] + bed[snp_ibd1[0,1], snp]
                    expected_results = [f, 1+f, 1+2*f, 2+f, 3+f]
                    self.assertAlmostEqual(result, expected_results[int(sibsum)], 3, msg = "problem with type1"+str((result, expected_results[int(sibsum)]))+str((i,j))+ str(sibsum))

        for i in range(3):
            for j in range(3):
                for count in range(10):
                    snp_ibd2[count] = [i, j]
                    result = impute_snp(snp, snp_ibd0, snp_ibd1, snp_ibd2, f, bed, 0, 0, count+1)
                    sibsum = bed[snp_ibd2[0,0], snp] + bed[snp_ibd2[0,1], snp]
                    self.assertAlmostEqual(result, sibsum/2. + 2*f, 3, msg = "problem with type2")
