import numpy as np
import logging
import doctest
import unittest
import os.path
import time

from pysnptools.kernelreader import *
from pysnptools.snpreader import Bed
from pysnptools.util import create_directory_if_necessary
from pysnptools.pstreader import PstReader
from pysnptools.snpreader import SnpData
import pysnptools.standardizer as stdizer
from six.moves import range

class _fortesting_JustCheckExists(object): #Implements ICopier

    def __init__(self,doPrintOutputNames=False):
        self.doPrintOutputNames = doPrintOutputNames
    
    def input(self,item):
        if isinstance(item, str):
            if not os.path.exists(item): raise Exception("Missing input file '{0}'".format(item))
        elif hasattr(item,"copyinputs"):
            item.copyinputs(self)
        # else -- do nothing

    def output(self,item):
        if isinstance(item, str):
            if not os.path.exists(item): raise Exception("Missing output file '{0}'".format(item))
            if self.doPrintOutputNames:
                print(item)
        elif hasattr(item,"copyoutputs"):
            item.copyoutputs(self)
        # else -- do nothing




class TestKernelReader(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.currentFolder = os.path.dirname(os.path.realpath(__file__))

    def test_merge_std(self):
        #unit vs beta
        for std in [stdizer.Beta(2,10),stdizer.Unit()]:
            np.random.seed(0)
            snp_count = 20
            snpreader = SnpData(iid=[["0","0"],["1","1"],["2","2"]],sid=[str(i) for i in range(snp_count)],val=np.array(np.random.randint(3,size=[3,snp_count]),dtype=np.float64,order='F'))
            kerneldata0,trained0,diag0 = SnpKernel(snpreader,std,block_size=1)._read_with_standardizing(to_kerneldata=True,return_trained=True)
            kerneldata1,trained1,diag1 = SnpKernel(snpreader,std,block_size=None)._read_with_standardizing(to_kerneldata=True,return_trained=True)
            np.testing.assert_array_almost_equal(kerneldata0.val,kerneldata1.val, decimal=10)
            np.testing.assert_array_almost_equal(trained0.stats,trained1.stats, decimal=10)
            assert abs(diag0.factor-diag1.factor) < 1e-7
    
    def test_cpp_std(self):

        #Order C vs F
        for order in ['C','F']:
            #32 vs 64
            for dtype in [np.float64,np.float32]:
                #unit vs beta
                for std in [stdizer.Unit(),stdizer.Beta(2,10)]:
                        np.random.seed(0)
                        snp_count = 20
                        snpreader0 = SnpData(iid=[["0","0"],["1","1"],["2","2"]],sid=[str(i) for i in range(snp_count)],val=np.array(np.random.randint(3,size=[3,snp_count]),dtype=dtype,order=order))
                        snpreader1 = SnpData(iid=[["3","3"],["4","4"]],sid=[str(i) for i in range(snp_count)],val=np.array(np.random.randint(3,size=[2,snp_count]),dtype=dtype,order=order))

                        #has SNC
                        for has_SNC_in_train in [False, True]:
                            if has_SNC_in_train:
                                snpreader0.val[:,1] = 0

                            #missing data
                            for has_missing_data in [False, True]:
                                if has_missing_data:
                                    snpreader0.val[0,2]=np.nan
                                    snpreader1.val[0,2]=np.nan

                                #gather stats vs not
                                cppa, stdcppa = snpreader0.read(order=order,dtype=dtype).standardize(std,return_trained=True,force_python_only=False)
                                pya, stdpya = snpreader0.read(order=order,dtype=dtype).standardize(std,return_trained=True,force_python_only=True)
                                np.testing.assert_array_almost_equal(cppa.val, pya.val, decimal=10 if dtype==np.float64 else 5)

                                np.testing.assert_array_almost_equal(stdcppa.stats,stdpya.stats, decimal=10 if dtype==np.float64 else 5)
                                assert (np.inf in stdcppa.stats[:,1]) == has_SNC_in_train
                                assert (np.inf in stdpya.stats[:,1]) == has_SNC_in_train

                                if has_SNC_in_train:
                                    assert np.array_equal(cppa.val[:,1],np.zeros([cppa.val.shape[0]]))
                                    assert np.array_equal(pya.val[:,1],np.zeros([pya.val.shape[0]]))

                                if has_missing_data:
                                    assert 0 == cppa.val[0,2]
                                    assert 0 == pya.val[0,2]
                                        
                                #uses stats
                                cppb = snpreader1.read(order=order,dtype=dtype).standardize(stdcppa,force_python_only=False)
                                pyb = snpreader1.read(order=order,dtype=dtype).standardize(stdpya,force_python_only=True)
                                np.testing.assert_array_almost_equal(cppb.val, pyb.val, decimal=10 if dtype==np.float64 else 5)
                                np.testing.assert_array_almost_equal(stdcppa.stats,stdpya.stats, decimal=10 if dtype==np.float64 else 5) #Make sure we haven't messed up the train stats

                                if has_SNC_in_train:
                                    assert np.array_equal(cppb.val[:,1],np.zeros([cppb.val.shape[0]]))
                                    assert np.array_equal(pyb.val[:,1],np.zeros([pyb.val.shape[0]]))

                                if has_missing_data:
                                    assert cppb.val[0,2]==0
                                    assert pyb.val[0,2]==0
        logging.info("done with 'test_cpp_std'")


    def test_intersection(self):

        from pysnptools.standardizer import Unit
        from pysnptools.kernelreader import SnpKernel
        from pysnptools.snpreader import Pheno
        from pysnptools.kernelreader._subset import _KernelSubset
        from pysnptools.snpreader._subset import _SnpSubset
        from pysnptools.util import intersect_apply

        snps_all = Bed(self.currentFolder + "/../examples/toydata",count_A1=False)
        k = SnpKernel(snps_all,stdizer.Identity())

        pheno = Pheno(self.currentFolder + "/../examples/toydata.phe")
        pheno = pheno[1:,:] # To test intersection we remove a iid from pheno

        k1,pheno = intersect_apply([k,pheno]) #SnpKernel is special because it standardizes AFTER intersecting.
        assert isinstance(k1.snpreader,_SnpSubset) and not isinstance(k1,_KernelSubset)

        #What happens with fancy selection?
        k2 = k[::2]
        assert isinstance(k2,SnpKernel)

        logging.info("Done with test_intersection")






    def test_respect_inputs(self):
        np.random.seed(0)
        for dtype_start,decimal_start in [(np.float32,5),(np.float64,10)]:
            for order_start in ['F','C','A']:
                for snp_count in [20,2]:
                    snpdataX = SnpData(iid=[["0","0"],["1","1"],["2","2"]],sid=[str(i) for i in range(snp_count)],val=np.array(np.random.randint(3,size=[3,snp_count]),dtype=dtype_start,order=order_start))
                    for stdx in [stdizer.Beta(1,25),stdizer.Identity(),stdizer.Unit()]:
                        for snpreader0 in [snpdataX,snpdataX[:,1:]]:
                            snpreader1 = snpreader0[1:,:]

                            refdata0, trained_standardizer = snpreader0.read().standardize(stdx,return_trained=True,force_python_only=True)
                            refval0 = refdata0.val.dot(refdata0.val.T)
                            refdata1 = snpreader1.read().standardize(trained_standardizer,force_python_only=True)
                            refval1 = refdata0.val.dot(refdata1.val.T)
                            for dtype_goal,decimal_goal in [(np.float32,5),(np.float64,10)]:
                                for order_goal in ['F','C','A']:
                                    k = snpreader0.read_kernel(standardizer=stdx,block_size=1,order=order_goal,dtype=dtype_goal)
                                    PstReader._array_properties_are_ok(k.val,order_goal,dtype_goal)
                                    np.testing.assert_array_almost_equal(refval0,k.val, decimal=min(decimal_start,decimal_goal))

    def test_fail(self):
        did_fail = True
        try:
            kd = KernelData(iid=[["0"],["1"],["2"]],val=[[1,2,3],[4,5,6],[7,8,9]]) #Wrong iid shape
            did_fail = False
        except:
            pass
        assert did_fail, "The constructor should fail because the iid is the wrong shape"

    def test_kernel2(self):
        logging.info("in kernel2")
        kd = KernelData(iid=[["0","0"],["1","1"],["2","2"]],val=[[1,2,3],[4,5,6],[7,8,9]])
        assert np.array_equal(kd.iid_to_index([["1","1"],["2","2"]]),np.array([1,2]))
        assert np.array_equal(kd.iid0_to_index([["1","1"],["2","2"]]),np.array([1,2]))
        kd = kd.standardize()
        assert np.abs(np.diag(kd.val).sum()-3)<1e-7
        assert kd.iid1_count == 3

    def test_snp_kernel2(self):
        logging.info("in test_snp_kernel2")
        snpreader = Bed(self.currentFolder + "/../examples/toydata",count_A1=False)
        snpkernel = SnpKernel(snpreader,standardizer=stdizer.Beta(1,25))
        s  = str(snpkernel)
        _fortesting_JustCheckExists().input(snpkernel)
        
    def test_npz(self):
        logging.info("in test_npz")
        snpreader = Bed(self.currentFolder + "/../examples/toydata",count_A1=False)
        kerneldata1 = snpreader.read_kernel(standardizer=stdizer.Unit())
        s = str(kerneldata1)
        output = "tempdir/kernelreader/toydata.kernel.npz"
        create_directory_if_necessary(output)
        KernelNpz.write(output,kerneldata1)
        kernelreader2 = KernelNpz(output)
        kerneldata2 = kernelreader2.read()
        np.testing.assert_array_almost_equal(kerneldata1.val, kerneldata2.val, decimal=10)
        logging.info("done with test")

    def test_subset(self):
        logging.info("in test_subset")
        snpreader = Bed(self.currentFolder + "/../examples/toydata",count_A1=False)
        snpkernel = SnpKernel(snpreader,stdizer.Unit())
        krsub = snpkernel[::2,::2]
        kerneldata1 = krsub.read()
        expected = snpreader.read_kernel(stdizer.Unit())[::2].read()
        np.testing.assert_array_almost_equal(kerneldata1.val, expected.val, decimal=10)

        krsub2 = snpkernel[::2]
        kerneldata2 = krsub2.read()
        np.testing.assert_array_almost_equal(kerneldata2.val, expected.val, decimal=10)
        logging.info("done with test")

    def test_identity(self):
        logging.info("in test_identity")
        snpreader = Bed(self.currentFolder + "/../examples/toydata",count_A1=False)
        assert snpreader.iid is snpreader.row
        kid = Identity(snpreader.row)
        assert np.array_equal(kid.row,kid.iid) and np.array_equal(kid.iid,kid.iid0) and np.array_equal(kid.iid0,kid.iid1) and np.array_equal(kid.iid1, kid.col)
        np.testing.assert_array_almost_equal(kid.read().val, np.identity(snpreader.iid_count))

        #subset
        sub1 = kid[1:5]
        np.testing.assert_array_almost_equal(sub1.read().val, np.identity(4))

        #zero size
        sub2 = kid[0:0]
        np.testing.assert_array_almost_equal(sub2.read().val, np.identity(0))

    def test_identity_sub(self):
        logging.info("in test_identity_sub")
        snpreader = Bed(self.currentFolder + "/../examples/toydata",count_A1=False)
        assert snpreader.iid is snpreader.row
        kid = Identity(snpreader.row)
        sub3 = kid[::2,1:5]
        expected = np.identity(kid.iid_count)[::2,:][:,1:5]
        np.testing.assert_array_almost_equal(sub3.read().val,expected)

        logging.info("done with test")

    def test_underscore_read1(self):
        logging.info("in test_underscore_read1")
        snpreader = Bed(self.currentFolder + "/../examples/toydata",count_A1=False)
        assert snpreader.iid is snpreader.row
        kid = Identity(snpreader.row)
        sub3 = kid[::2,::2]
        expected = np.identity(kid.iid_count)[::2,:][:,::2]
        np.testing.assert_array_almost_equal(sub3.read().val,expected)

        logging.info("done with test")

    def test_underscore_read2(self):
        logging.info("in test_underscore_read2")
        snpreader = Bed(self.currentFolder + "/../examples/toydata",count_A1=False)
        assert snpreader.iid is snpreader.row
        kid = Identity(snpreader.row)
        sub3 = kid[::2,::2]
        expected = np.identity(kid.iid_count)[::2,:][:,::2]
        np.testing.assert_array_almost_equal(sub3.read().val,expected)

        logging.info("done with test")

# We do it this way instead of using doctest.DocTestSuite because doctest.DocTestSuite requires modules to be pickled, which python doesn't allow.
# We need tests to be pickleable so that they can be run on a cluster.
class TestKrDocStrings(unittest.TestCase):
    def test_kernelreader(self):
        import pysnptools.kernelreader.kernelreader
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))
        result = doctest.testmod(pysnptools.kernelreader.kernelreader)
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__

    def test_snpkernel(self):
        import pysnptools.kernelreader.snpkernel
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))
        result = doctest.testmod(pysnptools.kernelreader.snpkernel)
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__

    def test_kerneldata(self):
        import pysnptools.kernelreader.kerneldata
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))
        result = doctest.testmod(pysnptools.kernelreader.kerneldata)
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__

    def test_identity(self):
        import pysnptools.kernelreader.identity
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))
        result = doctest.testmod(pysnptools.kernelreader.identity,optionflags=doctest.ELLIPSIS|doctest.NORMALIZE_WHITESPACE)
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__

    def test_snphdf5(self):
        import pysnptools.kernelreader.kernelhdf5
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))
        result = doctest.testmod(pysnptools.kernelreader.kernelhdf5)
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__

    def test_snpnpz(self):
        import pysnptools.kernelreader.kernelnpz
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))
        result = doctest.testmod(pysnptools.kernelreader.kernelnpz)
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__



def getTestSuite():
    """
    set up composite test suite
    """
    
    test_suite = unittest.TestSuite([])
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestKernelReader))
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestKrDocStrings))
    return test_suite

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    suites = getTestSuite()
    r = unittest.TextTestRunner(failfast=False)
    r.run(suites)
