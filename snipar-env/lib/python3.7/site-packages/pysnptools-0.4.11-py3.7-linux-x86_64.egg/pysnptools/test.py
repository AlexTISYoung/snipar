from __future__ import print_function

import numpy as np
import scipy as sp
import logging
import doctest
import shutil

from pysnptools.snpreader import Bed
from pysnptools.snpreader import SnpHdf5, SnpNpz
from pysnptools.snpreader import Dat
from pysnptools.snpreader import Dense
from pysnptools.snpreader import Pheno
from pysnptools.snpreader import Ped
from pysnptools.snpreader import DistributedBed
from pysnptools.standardizer import Unit
from pysnptools.standardizer import Beta
from pysnptools.util import create_directory_if_necessary
from pysnptools.kernelreader.test import TestKernelReader
from pysnptools.kernelreader.test import TestKrDocStrings
from pysnptools.pstreader.test import TestPstReader
from pysnptools.pstreader.test import TestPstDocStrings
from pysnptools.pstreader.pstmemmap import TestPstMemMap
from pysnptools.snpreader.snpmemmap import TestSnpMemMap
from pysnptools.snpreader.snpgen import TestSnpGen
from pysnptools.snpreader.distributedbed import TestDistributedBed
from pysnptools.util.generate import TestGenerate
from pysnptools.kernelreader.test import _fortesting_JustCheckExists
from pysnptools.util.intrangeset import TestIntRangeSet
from pysnptools.util.test import TestUtilTools
from pysnptools.util.filecache.test import TestFileCache

import unittest
import os.path
import time

from six.moves import range
   
class TestPySnpTools(unittest.TestCase):     

    def xtest_aaa_hdf5_speed(self): #!!too slow to use all the time

        #currentFolder + "/examples/toydata"
        #currentFolder + "/examples/delme.hdf5"
        bedFileName = r"c:\Source\carlk\fastlmm\tests\datasets\all_chr.maf0.001.N300" #!! local paths
        hdf5Pattern = r"c:\Source\carlk\fastlmm\tests\datasets\del.{0}.hdf5"#!!
        tt0 = time.time()
        snpreader_bed = Bed(bedFileName,count_A1=False)

        S0 = snpreader_bed.sid_count
        snp_index_list0 = range(min(S0,15000)) 

        hdf5FileName = hdf5Pattern.format(len(snp_index_list0))

        snpDataBed = snpreader_bed[:,snp_index_list0].read()
        tt1 = time.time()
        logging.info("Read bed %.2f seconds" % (tt1 - tt0))
        SnpHdf5.write(hdf5FileName, snpDataBed)
        tt2 = time.time()
        logging.info("write SnpHdf5 bed %.2f seconds" % (tt2 - tt1))

        snpreader_hdf5 = SnpHdf5(hdf5FileName)
        assert(snpreader_hdf5.iid_count == snpreader_bed.iid_count)
        S = snpreader_hdf5.sid_count
        N_original = snpreader_hdf5.iid_count
        iid_index_list = sorted(range(N_original - 1,0,-2))

        snp_index_list = sorted(range(S - 1,0,-2))#!!
        #snp_index_list = range(S//2)

        snpreader_hdf5 = snpreader_hdf5[iid_index_list,:]
        snpDataHdf5 = snpreader_hdf5[:,snp_index_list].read()
        tt3 = time.time()
        logging.info("read SnpHdf5 with reversed indexes bed %.2f seconds" % (tt3 - tt2))

        snpDataHdf5C = snpreader_hdf5[:,snp_index_list].read(order = "C")
        tt4 = time.time()
        logging.info("read SnpHdf5 C with reversed indexes bed %.2f seconds" % (tt4 - tt3))

        print("the end")

    @classmethod
    def setUpClass(self):
        self.currentFolder = os.path.dirname(os.path.realpath(__file__))
        #TODO: get data set with NANs!
        snpreader = Bed(self.currentFolder + "/examples/toydata",count_A1=False)
        self.pheno_fn = self.currentFolder + "/examples/toydata.phe"
        self.snpdata = snpreader.read(order='F',force_python_only=True)
        self.snps = self.snpdata.val

    def test_diagKtoN(self):
        """
        make sure standardization on SNPs results in sum(diag(K))=N
        """
        
        np.random.seed(42)
        m = np.random.random((100,1000))
        from pysnptools.standardizer import DiagKtoN
        s = DiagKtoN()
        s.standardize(m)
        K = m.dot(m.T)
        sum_diag = np.sum(np.diag(K))
        
        np.testing.assert_almost_equal(100, sum_diag)
        
    def test_c_reader_bed(self):
        snpreader = Bed(self.currentFolder + "/examples/toydata",count_A1=False)
        self.c_reader(snpreader)

    def test_c_reader_bed_count_A1(self):
        snpreader = Bed(self.currentFolder + "/examples/toydata",count_A1=True)
        snpdata = snpreader.read()
        snpdata.val = 2 - snpdata.val
        self.c_reader(snpdata)

    def test_p_reader_bed(self):
        snpreader = Bed(self.currentFolder + "/examples/toydata",count_A1=False).read(force_python_only=True)
        self.c_reader(snpreader)

    def test_p_reader_bed_count_A1(self):
        snpreader = Bed(self.currentFolder + "/examples/toydata",count_A1=True)
        snpdata = snpreader.read(force_python_only=True)
        snpdata.val = 2 - snpdata.val
        self.c_reader(snpdata)

    def test_scalar_index(self):
        snpreader = Bed(self.currentFolder + "/examples/toydata",count_A1=False)
        arr=np.int64(1)
        snpreader[arr,arr]

    def test_c_reader_hdf5(self):
        snpreader = SnpHdf5(self.currentFolder + "/examples/toydata.snpmajor.snp.hdf5")
        self.c_reader(snpreader)

    def test_hdf5_case3(self):
        snpreader1 = SnpHdf5(self.currentFolder + "/examples/toydata.snpmajor.snp.hdf5")[::2,:]
        snpreader2 = Bed(self.currentFolder + "/examples/toydata",count_A1=False)[::2,:]
        self.assertTrue(np.allclose(snpreader1.read().val, snpreader2.read().val, rtol=1e-05, atol=1e-05))

    def test_standardize_hdf5(self):
        snpreader = SnpHdf5(self.currentFolder + "/examples/toydata.iidmajor.snp.hdf5")
        self.standardize(snpreader)


    def test_c_reader_npz(self):
        snpreader = SnpNpz(self.currentFolder + "/examples/toydata10.snp.npz")
        snpdata = snpreader.read(order='F',force_python_only=False)
        snp_c = snpdata.val
        
        self.assertEqual(np.float64, snp_c.dtype)
        self.assertTrue(np.allclose(self.snps[:,:10], snp_c, rtol=1e-05, atol=1e-05))

        snpreader1 = SnpNpz(self.currentFolder + "/examples/toydata10.snp.npz")
        snpreader2 = Bed(self.currentFolder + "/examples/toydata",count_A1=False)[:,:10]
        self.assertTrue(np.allclose(snpreader1.read().val, snpreader2.read().val, rtol=1e-05, atol=1e-05))


        snpdata.val[1,2] = np.NaN # Inject a missing value to test writing and reading missing values
        output = "tempdir/snpreader/toydata10.snp.npz"
        create_directory_if_necessary(output)
        SnpNpz.write(output,snpdata)
        snpdata2 = SnpNpz(output).read()
        np.testing.assert_array_almost_equal(snpdata.val, snpdata2.val, decimal=10)

        snpdata3 = snpdata[:,0:0].read() #create snpdata with no sids
        output = "tempdir/snpreader/toydata0.snp.npz"
        SnpNpz.write(output,snpdata3)
        snpdata4 = SnpNpz(output).read()
        assert snpdata3 == snpdata4



    def test_standardize_npz(self):
        snpreader = SnpNpz(self.currentFolder + "/examples/toydata10.snp.npz")
        self.standardize(snpreader)


    def test_c_reader_dat(self):
        snpreader = Dat(self.currentFolder + "/examples/toydata.dat")[:,::100]
        _fortesting_JustCheckExists().input(snpreader)

        snpdata1 = snpreader.read()
        self.assertEqual(np.float64, snpdata1.val.dtype)
        self.assertTrue(np.allclose(self.snps[:,::100], snpdata1.val, rtol=1e-05, atol=1e-05))

        snpdata1.val[1,2] = np.NaN # Inject a missing value to test writing and reading missing values
        output = "tempdir/snpreader/toydata.dat"
        create_directory_if_necessary(output)
        Dat.write(output,snpdata1)
        snpdata2 = Dat(output).read()
        np.testing.assert_array_almost_equal(snpdata1.val, snpdata2.val, decimal=10)

        snpdata3 = snpdata1[:,0:0].read() #create snpdata with no sids
        output = "tempdir/snpreader/toydata3.dat"
        Dat.write(output,snpdata3)
        snpdata4 = Dat(output).read()
        assert snpdata3 == snpdata4

    @staticmethod
    def assert_match_012_210(snpdata1, snpdata2):
        for sid_index in range(snpdata1.sid_count): #Check that every row matches (except OK if 0,1,2 can be 2,1,0)
            goods1 = (snpdata1.val[:,sid_index] == snpdata1.val[:,sid_index]) # find non-missing
            goods2= (snpdata2.val[:,sid_index] == snpdata2.val[:,sid_index])  # find non-missing
            assert (goods1==goods2).all() #assert that they agree on non-missing
            is_ok = (snpdata1.val[goods1,sid_index] == snpdata2.val[goods2,sid_index]).all() or (snpdata1.val[goods1,sid_index] == snpdata2.val[goods2,sid_index]*-1+2).all()
            assert is_ok

    def test_c_reader_ped(self):
        if False: #Too slow for routine testing
            snpdata1 = Ped(self.currentFolder + "/examples/toydata.ped")[::25,::1000].read()
            self.assertEqual(np.float64, snpdata1.val.dtype)
            TestPySnpTools.assert_match_012_210(self.snpdata[::25,::1000].read(),snpdata1)
        else:
            snpdata1 = self.snpdata[::25,::1000].read()

        output = "tempdir/snpreader/toydata.ped"
        create_directory_if_necessary(output)

        snpdata1.val[1,2] = np.NaN # Inject a missing value to test writing and reading missing values
        Ped.write(output, snpdata1)
        snpreader = Ped(output)
        _fortesting_JustCheckExists().input(snpreader)
        s = str(snpreader)
        snpdata2 = snpreader.read()
        TestPySnpTools.assert_match_012_210(snpdata1,snpdata2)

    def test_c_reader_pheno(self):
        snpdata1 = Pheno(self.currentFolder + "/examples/toydata.phe").read()

        self.assertEqual(np.float64, snpdata1.val.dtype)

        snpdata1.val[1,0] = np.NaN # Inject a missing value to test writing and reading missing values
        output = "tempdir/snpreader/toydata.phe"
        create_directory_if_necessary(output)
        Pheno.write(output, snpdata1)
        snpreader = Pheno(output)
        _fortesting_JustCheckExists().input(snpreader)
        s = str(snpreader)
        snpdata2 = snpreader.read()
        np.testing.assert_array_almost_equal(snpdata1.val, snpdata2.val, decimal=10)

        snpdata1 = Pheno(self.currentFolder + "/examples/toydata.phe").read()
        import pysnptools.util.pheno as pstpheno
        dict = pstpheno.loadOnePhen(self.currentFolder + "/examples/toydata.phe",missing="")
        snpdata3 = Pheno(dict).read()
        np.testing.assert_array_almost_equal(snpdata1.val, snpdata3.val, decimal=10)


        dict = pstpheno.loadOnePhen(self.currentFolder + "/examples/toydata.phe",missing="",vectorize=True)
        assert len(dict['vals'].shape)==1, "test 1-d array of values"
        snpdata3 = Pheno(dict).read()
        np.testing.assert_array_almost_equal(snpdata1.val, snpdata3.val, decimal=10)

        snpdata4 = Pheno(None,iid_if_none=snpdata1.iid)
        assert (snpdata4.row == snpdata1.row).all() and snpdata4.col_count == 0

        snpdata5 = Pheno(self.currentFolder + "/examples/toydata.id.phe").read()
        np.testing.assert_array_almost_equal(snpdata1.val, snpdata5.val, decimal=10)
        snpdata6 = Pheno(self.currentFolder + "/examples/toydata.fid.phe").read()
        np.testing.assert_array_almost_equal(snpdata1.val, snpdata6.val, decimal=10)

    def test_c_reader_dense(self):
        snpdata1 = self.snpdata[:,::100].read()
        snpdata1.val[1,2] = np.NaN # Inject a missing value to test writing and reading missing values
        output = "tempdir/snpreader/toydata.dense.txt"
        create_directory_if_necessary(output)
        Dense.write(output, snpdata1)
        snpreader = Dense(output)
        _fortesting_JustCheckExists().input(snpreader)
        snpdata2 = snpreader.read()
        np.testing.assert_array_almost_equal(snpdata1.val, snpdata2.val, decimal=10)

    def test_c_reader_distributedbed(self):
        from pysnptools.util.filecache import LocalCache

        snpdata1 = self.snpdata[:,::100].read()
        snpdata1.val[1,2] = np.NaN # Inject a missing value to test writing and reading missing values
        output = "tempdir/snpreader/toydata.distributedbed"
        LocalCache(output).rmtree()
        DistributedBed.write(output, snpdata1, piece_per_chrom_count=5)
        snpreader = DistributedBed(output)
        _fortesting_JustCheckExists().input(snpreader)
        snpdata2 = snpreader.read()
        np.testing.assert_array_almost_equal(snpdata1.val, snpdata2.val, decimal=10)

    def test_some_std(self):
        k0 = self.snpdata.read_kernel(standardizer=Unit()).val
        from pysnptools.kernelreader import SnpKernel
        k1 = self.snpdata.read_kernel(standardizer=Unit())
        np.testing.assert_array_almost_equal(k0, k1.val, decimal=10)

        from pysnptools.snpreader import SnpData
        snpdata2 = SnpData(iid=self.snpdata.iid,sid=self.snpdata.sid,pos=self.snpdata.pos,val=np.array(self.snpdata.val))
        s = str(snpdata2)
        snpdata2.standardize()
        s = str(snpdata2)

        snpreader = Bed(self.currentFolder + "/examples/toydata",count_A1=False)
        k2 = snpreader.read_kernel(standardizer=Unit(),block_size=500).val
        np.testing.assert_array_almost_equal(k0, k2, decimal=10)

        from pysnptools.standardizer.identity import Identity
        from pysnptools.standardizer.diag_K_to_N import DiagKtoN
        for dtype in [sp.float64,sp.float32]:
            for std in [Unit(),Beta(1,25),Identity(),DiagKtoN()]:
                s = str(std)
                np.random.seed(0)
                x = np.array(np.random.randint(3,size=[60,100]),dtype=dtype)
                x2 = x[:,::2]
                x2b = np.array(x2)
                #LATER what's this about? It doesn't do non-contiguous?
                #assert not x2.flags['C_CONTIGUOUS'] and not x2.flags['F_CONTIGUOUS'] #set up to test non contiguous
                #assert x2b.flags['C_CONTIGUOUS'] or x2b.flags['F_CONTIGUOUS'] #set up to test non contiguous
                #a,b = std.standardize(x2b),std.standardize(x2)
                #np.testing.assert_array_almost_equal(a,b)
        logging.info("done")

    def c_reader(self,snpreader):
        """
        make sure c-reader yields same result
        """
        snpdata = snpreader.read(order='F',force_python_only=False)
        snp_c = snpdata.val
        
        self.assertEqual(np.float64, snp_c.dtype)
        self.assertTrue(np.allclose(self.snps, snp_c, rtol=1e-05, atol=1e-05))
        return snpdata

    def test_standardize_bed(self):
        snpreader = Bed(self.currentFolder + "/examples/toydata",count_A1=False)
        self.standardize(snpreader)

    def test_standardize_dat(self):
        snpreader = Dat(self.currentFolder + "/examples/toydata.dat")
        self.standardize(snpreader)

    def test_standardize_ped(self):
        snpreader = Ped(self.currentFolder + "/examples/toydata")
        self.standardize(snpreader)

    def standardize(self,snpreader):
        """
        make sure blocked standardize yields same result as regular standardize
        """

        for dtype in [sp.float64,sp.float32]:

            snps = snpreader.read(order='F',force_python_only=True,dtype=dtype).val
            self.assertEqual(dtype, snps.dtype)

            snp_s1 = Unit().standardize(snps.copy(), force_python_only=True)
            snp_s2 = Unit().standardize(snps.copy(), block_size=100, force_python_only=True)
            snps_F = np.array(snps, dtype=dtype, order="F")
            snp_s3 = Unit().standardize(snps_F)
            snps_C = np.array(snps, dtype=dtype, order="C")
            snp_s4 = Unit().standardize(snps_C)

            snp_beta1 = Beta(1, 25).standardize(snps.copy(), force_python_only=True)
            snps_F = np.array(snps, dtype=dtype, order="F")
            snp_beta2 = Beta(1, 25).standardize(snps_F)
            snps_C = np.array(snps, dtype=dtype, order="C")
            snp_beta3 = Beta(1, 25).standardize(snps_C)

            self.assertEqual(snp_s1.shape[0], snp_s2.shape[0])
            self.assertEqual(snp_s1.shape[1], snp_s2.shape[1])

            self.assertEqual(snp_s1.shape[0], snp_s3.shape[0])
            self.assertEqual(snp_s1.shape[1], snp_s3.shape[1])
        
            self.assertEqual(snp_s1.shape[0], snp_s4.shape[0])
            self.assertEqual(snp_s1.shape[1], snp_s4.shape[1])

            self.assertTrue(np.allclose(snp_s1, snp_s2, rtol=1e-05, atol=1e-05))
            self.assertTrue(np.allclose(snp_s1, snp_s3, rtol=1e-05, atol=1e-05))
            self.assertTrue(np.allclose(snp_s1, snp_s4, rtol=1e-05, atol=1e-05))

            self.assertEqual(snp_beta1.shape[0], snp_beta2.shape[0])
            self.assertEqual(snp_beta1.shape[1], snp_beta2.shape[1])
            self.assertEqual(snp_beta1.shape[0], snp_beta3.shape[0])
            self.assertEqual(snp_beta1.shape[1], snp_beta3.shape[1])
        
            self.assertTrue(np.allclose(snp_beta1, snp_beta2, rtol=1e-05, atol=1e-05))
            self.assertTrue(np.allclose(snp_beta1, snp_beta3, rtol=1e-05, atol=1e-05))

    def test_load_and_standardize_bed(self):
        snpreader2 = Bed(self.currentFolder + "/examples/toydata",count_A1=False)
        self.load_and_standardize(snpreader2, snpreader2)

    def too_slow_test_write_bedbig(self):
        iid_count = 100000
        sid_count = 50000
        from pysnptools.snpreader import SnpData
        iid = np.array([[str(i),str(i)] for i in range(iid_count)])
        sid = np.array(["sid_{0}".format(i) for i in range(sid_count)])
        pos = np.array([[i,i,i] for i in range(sid_count)])
        np.random.seed(0)
        snpdata = SnpData(iid,sid,np.zeros((iid_count,sid_count)),pos=pos) #random.choice((0.0,1.0,2.0,float("nan")),size=(iid_count,sid_count)))
        output = "tempdir/bedbig.{0}.{1}".format(iid_count,sid_count)
        create_directory_if_necessary(output)
        Bed.write(output, snpdata, count_A1=False)
        snpdata2 = Bed(output,count_A1=False).read()
        np.testing.assert_array_almost_equal(snpdata.val, snpdata2.val, decimal=10)

    def test_write_bed_f64cpp_0(self):
        snpreader = Bed(self.currentFolder + "/examples/toydata",count_A1=False)
        iid_index = 0
        logging.info("iid={0}".format(iid_index))
        #if snpreader.iid_count % 4 == 0: # divisible by 4 isn't a good test
        #    snpreader = snpreader[0:-1,:]
        #assert snpreader.iid_count % 4 != 0
        snpdata = snpreader[0:iid_index,:].read(order='F',dtype=np.float64)
        if snpdata.iid_count > 0:
            snpdata.val[-1,0] = float("NAN")
        output = "tempdir/toydata.F64cpp.{0}".format(iid_index)
        create_directory_if_necessary(output)
        Bed.write(output, snpdata ,count_A1=False)
        snpdata2 = Bed(output,count_A1=False).read()
        np.testing.assert_array_almost_equal(snpdata.val, snpdata2.val, decimal=10)

    def test_write_bed_f64cpp_1(self):
        snpreader = Bed(self.currentFolder + "/examples/toydata",count_A1=False)
        iid_index = 1
        logging.info("iid={0}".format(iid_index))
        #if snpreader.iid_count % 4 == 0: # divisible by 4 isn't a good test
        #    snpreader = snpreader[0:-1,:]
        #assert snpreader.iid_count % 4 != 0
        snpdata = snpreader[0:iid_index,:].read(order='F',dtype=np.float64)
        if snpdata.iid_count > 0:
            snpdata.val[-1,0] = float("NAN")
        output = "tempdir/toydata.F64cpp.{0}".format(iid_index)
        create_directory_if_necessary(output)
        Bed.write(output, snpdata ,count_A1=False)
        snpdata2 = Bed(output,count_A1=False).read()
        np.testing.assert_array_almost_equal(snpdata.val, snpdata2.val, decimal=10)

    def test_write_bed_f64cpp_5(self):
        snpreader = Bed(self.currentFolder + "/examples/toydata",count_A1=False)

        from pysnptools.kernelreader.test import _fortesting_JustCheckExists
        _fortesting_JustCheckExists().input(snpreader)

        iid_index = 5
        logging.info("iid={0}".format(iid_index))
        #if snpreader.iid_count % 4 == 0: # divisible by 4 isn't a good test
        #    snpreader = snpreader[0:-1,:]
        #assert snpreader.iid_count % 4 != 0
        snpdata = snpreader[0:iid_index,:].read(order='F',dtype=np.float64)
        if snpdata.iid_count > 0:
            snpdata.val[-1,0] = float("NAN")
        output = "tempdir/toydata.F64cpp.{0}".format(iid_index)
        create_directory_if_necessary(output)
        Bed.write(output, snpdata ,count_A1=False) #,force_python_only=True)
        snpdata2 = Bed(output,count_A1=False).read()
        np.testing.assert_array_almost_equal(snpdata.val, snpdata2.val, decimal=10)

    def test_write_bed_f64cpp_5_python(self):
        snpreader = Bed(self.currentFolder + "/examples/toydata",count_A1=False)
        iid_index = 5
        logging.info("iid={0}".format(iid_index))
        #if snpreader.iid_count % 4 == 0: # divisible by 4 isn't a good test
        #    snpreader = snpreader[0:-1,:]
        #assert snpreader.iid_count % 4 != 0
        snpdata = snpreader[0:iid_index,:].read(order='F',dtype=np.float64)
        if snpdata.iid_count > 0:
            snpdata.val[-1,0] = float("NAN")
        output = "tempdir/toydata.F64python.{0}".format(iid_index)
        create_directory_if_necessary(output)
        Bed.write(output,snpdata, force_python_only=True)
        snpdata2 = Bed(output,count_A1=False).read()
        np.testing.assert_array_almost_equal(snpdata.val, snpdata2.val, decimal=10)

    def test_write_x_x_cpp(self):
        for count_A1 in [False, True]:
            snpreader = Bed(self.currentFolder + "/examples/toydata",count_A1=count_A1)
            for order in ['C','F']:
                for dtype in [np.float32,np.float64]:
                    snpdata = snpreader.read(order=order,dtype=dtype)
                    snpdata.val[-1,0] = float("NAN")
                    output = "tempdir/toydata.{0}{1}.cpp".format(order,"32" if dtype==np.float32 else "64")
                    create_directory_if_necessary(output)
                    Bed.write(output, snpdata, count_A1=count_A1)
                    snpdata2 = Bed(output,count_A1=count_A1).read()
                    np.testing.assert_array_almost_equal(snpdata.val, snpdata2.val, decimal=10)

    def test_subset_view(self):
        snpreader2 = Bed(self.currentFolder + "/examples/toydata",count_A1=False)[:,:]
        result = snpreader2.read(view_ok=True)
        self.assertFalse(snpreader2 is result)
        result2 = result[:,:].read()
        self.assertFalse(sp.may_share_memory(result2.val,result.val))
        result3 = result[:,:].read(view_ok=True)
        self.assertTrue(sp.may_share_memory(result3.val,result.val))
        result4 = result3.read()
        self.assertFalse(sp.may_share_memory(result4.val,result3.val))
        result5 = result4.read(view_ok=True)
        self.assertTrue(sp.may_share_memory(result4.val,result5.val))

    def test_load_and_standardize_hdf5(self):
        snpreader2 = SnpHdf5(self.currentFolder + "/examples/toydata.snpmajor.snp.hdf5")
        snpreader3 = SnpHdf5(self.currentFolder + "/examples/toydata.iidmajor.snp.hdf5")
        self.load_and_standardize(snpreader2, snpreader3)
        snpreaderref = Bed(self.currentFolder + "/examples/toydata",count_A1=False)
        self.load_and_standardize(snpreader2, snpreaderref)

    def test_load_and_standardize_dat(self):
        snpreader2 = Dat(self.currentFolder + "/examples/toydata.dat")
        self.load_and_standardize(snpreader2, snpreader2)
        #snpreaderref = Bed(self.currentFolder + "/examples/toydata",count_A1=False)
        #self.load_and_standardize(snpreader2, snpreaderref)

    def test_load_and_standardize_ped(self):

        #!!Ped columns can be ambiguous
        ###Creating Ped data ...
        #currentFolder = os.path.dirname(os.path.realpath(__file__))
        #snpData = Bed(currentFolder + "/examples/toydata",count_A1=False).read()
        ##Ped.write(snpData, currentFolder + "/examples/toydata.ped")
        #fromPed = Ped(currentFolder + "/examples/toydata").read()
        #self.assertTrue(np.allclose(snpData.val, fromPed.val, rtol=1e-05, atol=1e-05))


        snpreader2 = Ped(self.currentFolder + "/examples/toydata")
        self.load_and_standardize(snpreader2, snpreader2)
        #snpreaderref = Bed(self.currentFolder + "/examples/toydata",count_A1=False)
        #self.load_and_standardize(snpreader2, snpreaderref)

    def load_and_standardize(self, snpreader2, snpreader3):
        """
        test c-version of load and standardize
        """

        S = snpreader2.sid_count
        N_original = snpreader2.iid_count

        iid_index_list = range(N_original - 1,0,-2)
        snpreader3 = snpreader3[iid_index_list,:]

        for dtype in [sp.float64,sp.float32]:

            G2 = snpreader2.read(order='F',force_python_only=True).val
            G2 = Unit().standardize(G2, block_size=10000, force_python_only=True)

            SNPs_floatF = snpreader2.read(order="F", dtype=dtype, force_python_only=False).val
            GF = Unit().standardize(SNPs_floatF)

            SNPs_floatC = snpreader2.read(order="C", dtype=dtype, force_python_only=False).val
            GC = Unit().standardize(SNPs_floatC)

            self.assertTrue(np.allclose(GF, G2, rtol=1e-05, atol=1e-05))
            self.assertTrue(np.allclose(GF, GC, rtol=1e-05, atol=1e-05))

            #testing selecting a subset of snps and iids
            snp_index_list = range(S - 1,0,-2)

            G2x = snpreader2.read(order='F',force_python_only=True).val
            G2x = G2x[iid_index_list,:][:,snp_index_list]
            G2x = Unit().standardize(G2x, block_size=10000, force_python_only=True)


            SNPs_floatFx = snpreader3[:,snp_index_list].read(order="F", dtype=dtype, force_python_only=False).val
            GFx = Unit().standardize(SNPs_floatFx)
            self.assertTrue(np.allclose(GFx, G2x, rtol=1e-05, atol=1e-05))

            SNPs_floatCx = snpreader3[:,snp_index_list].read(order="C", dtype=dtype, force_python_only=False).val
            GCx = Unit().standardize(SNPs_floatCx)
            self.assertTrue(np.allclose(GFx, G2x, rtol=1e-05, atol=1e-05))


    def test_writes(self):
        from pysnptools.snpreader import SnpData, SnpHdf5, SnpNpz, SnpMemMap

        the_class_and_suffix_list = [(DistributedBed, "distributed_bed"),(Dense,"dense"),(Bed,"bed"),(Dat,"dat"),(Ped,"ped"),(Pheno,"pheno"),
                                    (SnpHdf5,"hdf5"),(SnpNpz,"npz"),(SnpMemMap,"memmap")]
        cant_do_col_prop_none_set = {"dense","distributed_bed"}
        cant_do_col_len_0_set = {"distributed_bed"}
        cant_do_row_count_zero_set = {'dense','ped','pheno'}
        can_swap_0_2_set = {'ped'}
        can_change_col_names_set = {'pheno'}
        ignore_fam_id_set = {'dense'}
        ignore_pos_set = {'dense','pheno'}
        erase_any_write_dir = {'distributed_bed'}

        
        #===================================
        #    Starting main function
        #===================================
        logging.info("starting 'test_writes'")
        np.random.seed(0)
        output_template = "tempdir/snpreader/writes.{0}.{1}"
        create_directory_if_necessary(output_template.format(0,"npz"))
        i = 0
        for row_count in [0,5,2,1]:
            for col_count in [4,2,1,0]:
                val = np.random.random_integers(low=0,high=3,size=(row_count,col_count))*1.0
                val[val==3]=np.NaN
                row = [('0','0'),('1','1'),('2','2'),('3','3'),('4','4')][:row_count]
                col = ['s0','s1','s2','s3','s4'][:col_count]
                for is_none in [True,False]:
                    row_prop = None
                    col_prop = None if is_none else [(x,x,x) for x in range(5)][:col_count]
                    snpdata = SnpData(iid=row,sid=col,val=val,pos=col_prop,name=str(i))
                    for the_class,suffix in the_class_and_suffix_list:
                        if col_count == 0 and suffix in cant_do_col_len_0_set:
                            continue
                        if col_prop is None and suffix in cant_do_col_prop_none_set:
                            continue
                        if row_count==0 and suffix in cant_do_row_count_zero_set:
                            continue
                        filename = output_template.format(i,suffix)
                        logging.info(filename)
                        i += 1
                        if suffix in erase_any_write_dir and os.path.exists(filename):
                            shutil.rmtree(filename)
                        the_class.write(filename,snpdata)
                        for subsetter in [None, sp.s_[::2,::3]]:
                            reader = the_class(filename)
                            _fortesting_JustCheckExists().input(reader)
                            subreader = reader if subsetter is None else reader[subsetter[0],subsetter[1]]
                            readdata = subreader.read(order='C')
                            expected = snpdata if subsetter is None else snpdata[subsetter[0],subsetter[1]].read()
                            if not suffix in can_swap_0_2_set:
                                assert np.allclose(readdata.val,expected.val,equal_nan=True)
                            else:
                                for col_index in range(readdata.col_count):
                                    assert (np.allclose(readdata.val[:,col_index],expected.val[:,col_index],equal_nan=True) or
                                            np.allclose(readdata.val[:,col_index]*-1+2,expected.val[:,col_index],equal_nan=True))
                            if not suffix in ignore_fam_id_set:
                                assert np.array_equal(readdata.row,expected.row)
                            else:
                                assert np.array_equal(readdata.row[:,1],expected.row[:,1])
                            if not suffix in can_change_col_names_set:
                                assert np.array_equal(readdata.col,expected.col)
                            else:
                                assert readdata.col_count==expected.col_count
                            assert np.array_equal(readdata.row_property,expected.row_property) or (readdata.row_property.shape[1]==0 and expected.row_property.shape[1]==0)

                            if not suffix in ignore_pos_set:
                                assert np.allclose(readdata.col_property,expected.col_property,equal_nan=True) or (readdata.col_property.shape[1]==0 and expected.col_property.shape[1]==0)
                            else:
                                assert len(readdata.col_property)==len(expected.col_property)
                        try:
                            os.remove(filename)
                        except:
                            pass
        logging.info("done with 'test_writes'")

class NaNCNCTestCases(unittest.TestCase):
    def __init__(self, iid_index_list, snp_index_list, standardizer, snpreader, dtype, order, force_python_only, reference_snps, reference_dtype):
        self.iid_index_list = iid_index_list
        self.snp_index_list = snp_index_list
        self.standardizer = standardizer
        self.snpreader = snpreader
        self.dtype = dtype
        self.order = order
        self.force_python_only = force_python_only
        self.reference_snps = reference_snps
        self.reference_dtype = reference_dtype

    _testMethodName = "runTest"
    _testMethodDoc = None

    @staticmethod
    def factory_iterator():

        snp_reader_factory_bed = lambda : Bed("examples/toydata",count_A1=False)
        snp_reader_factory_snpmajor_hdf5 = lambda : SnpHdf5("examples/toydata.snpmajor.snp.hdf5")
        snp_reader_factory_iidmajor_hdf5 = lambda : SnpHdf5("examples/toydata.iidmajor.snp.hdf5")
        snp_reader_factory_dat = lambda : Dat("examples/toydata.dat")

        previous_wd = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))

        snpreader0 = snp_reader_factory_bed()
        S_original = snpreader0.sid_count
        N_original = snpreader0.iid_count

        snps_to_read_count = min(S_original, 100)

        for iid_index_list in [range(N_original), range(N_original//2), range(N_original - 1,0,-2)]:
            for snp_index_list in [range(snps_to_read_count), range(snps_to_read_count//2), range(snps_to_read_count - 1,0,-2)]:
                for standardizer in [Unit(),Beta(1,25)]:
                    reference_snps, reference_dtype = NaNCNCTestCases(iid_index_list, snp_index_list, standardizer, snp_reader_factory_bed(), sp.float64, "C", "False", None, None).read_and_standardize()
                    for snpreader_factory in [snp_reader_factory_bed, 
                                             snp_reader_factory_snpmajor_hdf5, snp_reader_factory_iidmajor_hdf5,
                                             snp_reader_factory_dat
                                              ]:
                        for dtype in [sp.float64,sp.float32]:
                            for order in ["C", "F"]:
                                for force_python_only in [False, True]:
                                    snpreader = snpreader_factory()
                                    test_case = NaNCNCTestCases(iid_index_list, snp_index_list, standardizer, snpreader, dtype, order, force_python_only, reference_snps, reference_dtype)
                                    yield test_case
        os.chdir(previous_wd)

    def __str__(self):
        iid_index_list = self.iid_index_list
        snp_index_list = self.snp_index_list
        standardizer = self.standardizer
        snpreader = self.snpreader
        dtype = self.dtype
        order = self.order
        force_python_only = self.force_python_only
        return "{0}(iid_index_list=[{1}], snp_index_list=[{2}], standardizer={3}, snpreader={4}, dtype={5}, order='{6}', force_python_only=={7})".format(
            self.__class__.__name__,
            ",".join([str(i) for i in iid_index_list]) if len(iid_index_list) < 10 else ",".join([str(i) for i in iid_index_list[0:10]])+",...",
            ",".join([str(i) for i in snp_index_list]) if len(snp_index_list) < 10 else ",".join([str(i) for i in snp_index_list[0:10]])+",...",
            standardizer, snpreader, dtype, order, force_python_only)

    def read_and_standardize(self):
        previous_wd = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))
        iid_index_list = self.iid_index_list
        snp_index_list = self.snp_index_list
        standardizer = self.standardizer
        snpreader = self.snpreader
        dtype = self.dtype
        order = self.order
        force_python_only = self.force_python_only
        
        snps = snpreader[iid_index_list,snp_index_list].read(order=order, dtype=dtype, force_python_only=force_python_only).val
        snps[0,0] = np.nan # add a NaN
        snps[:,1] = 2 # make a SNC
        snps = standardizer.standardize(snps,force_python_only=force_python_only)
        os.chdir(previous_wd)
        return snps, dtype

    def doCleanups(self):
        pass
        #return super(NaNCNCTestCases, self).doCleanups()

    def runTest(self, result = None):
        snps, dtype = self.read_and_standardize()
        self.assertTrue(snps[0,0] == 0)
        self.assertTrue(np.all(snps[:,1] == 0))
        if self.reference_snps is not None:
            self.assertTrue(np.allclose(self.reference_snps, snps, rtol=1e-04 if dtype == sp.float32 or self.reference_dtype == sp.float32 else 1e-12))

# We do it this way instead of using doctest.DocTestSuite because doctest.DocTestSuite requires modules to be pickled, which python doesn't allow.
# We need tests to be pickleable so that they can be run on a cluster.
class TestSnpDocStrings(unittest.TestCase):

    def test_bed(self):
        import pysnptools.snpreader.bed
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__))+"/snpreader")
        result = doctest.testmod(pysnptools.snpreader.bed)
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__

    def test_dat(self):
        import pysnptools.snpreader.snpdata
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__))+"/snpreader")
        result = doctest.testmod(pysnptools.snpreader.dat)
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__

    def test_dense(self):
        import pysnptools.snpreader.snpdata
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__))+"/snpreader")
        result = doctest.testmod(pysnptools.snpreader.dense)
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__

    def test_pairs(self):
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__))+"/snpreader")
        import pysnptools.snpreader.pairs
        result = doctest.testmod(pysnptools.snpreader.pairs)
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__

    def test_ped(self):
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__))+"/snpreader")
        import pysnptools.snpreader.ped
        result = doctest.testmod(pysnptools.snpreader.ped)
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__

    def test_pheno(self):
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__))+"/snpreader")
        import pysnptools.snpreader.pheno
        result = doctest.testmod(pysnptools.snpreader.pheno)
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__

    def test_snpdata(self):
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__))+"/snpreader")
        import pysnptools.snpreader.snpdata
        result = doctest.testmod(pysnptools.snpreader.snpdata)
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__

    def test_snphdf5(self):
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__))+"/snpreader")
        import pysnptools.snpreader.snphdf5
        result = doctest.testmod(pysnptools.snpreader.snphdf5)
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__

    def test_snpmemmap(self):
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__))+"/snpreader")
        import pysnptools.snpreader.snpmemmap
        result = doctest.testmod(pysnptools.snpreader.snpmemmap)
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__

    def test_snpgen(self):
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__))+"/util")
        import pysnptools.snpreader.snpgen
        result = doctest.testmod(pysnptools.snpreader.snpgen)
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__

    def test_snpnpz(self):
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__))+"/snpreader")
        import pysnptools.snpreader.snpnpz
        result = doctest.testmod(pysnptools.snpreader.snpnpz)
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__

    def test_distributedbed(self):
        import pysnptools.snpreader.distributedbed
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__))+"/snpreader")
        result = doctest.testmod(pysnptools.snpreader.distributedbed)
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__


    def test_snpreader(self):
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__))+"/snpreader")
        import pysnptools.snpreader.snpreader
        result = doctest.testmod(pysnptools.snpreader.snpreader,optionflags=doctest.ELLIPSIS|doctest.NORMALIZE_WHITESPACE)
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__

    def test_util(self):
        import pysnptools.util
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__))+"/util")
        import pysnptools.util
        result = doctest.testmod(pysnptools.util)
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__

    def test_standardize_testmod(self):
        import pysnptools.standardizer
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__))+"/standardizer")
        for mod in [
                    pysnptools.standardizer.beta,
                    pysnptools.standardizer.betatrained,
                    pysnptools.standardizer.bysidcount,
                    pysnptools.standardizer.bysqrtsidcount,
                    pysnptools.standardizer.diag_K_to_N,
                    pysnptools.standardizer.identity,
                    pysnptools.standardizer.standardizer,
                    pysnptools.standardizer.unit,
                    pysnptools.standardizer.unittrained]:
            result = doctest.testmod(mod,optionflags=doctest.ELLIPSIS|doctest.NORMALIZE_WHITESPACE)
            assert result.failed == 0, "failed doc test: " + __file__
        os.chdir(old_dir)

    def test_kernelstandardize_testmod(self):
        import pysnptools.kernelstandardizer
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__))+"/kernelstandardizer")
        for mod in [
                    pysnptools.kernelstandardizer,
                ]:
            result = doctest.testmod(mod)
            assert result.failed == 0, "failed doc test: " + __file__
        os.chdir(old_dir)

def getTestSuite():
    """
    set up composite test suite
    """

    test_suite = unittest.TestSuite([])

    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestPySnpTools))
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestDistributedBed))
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestFileCache))
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestUtilTools))
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestIntRangeSet))
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestSnpDocStrings))
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestPstDocStrings))
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestKrDocStrings))
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestSnpGen))
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestGenerate))
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestPstMemMap))
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestSnpMemMap))
    test_suite.addTests(NaNCNCTestCases.factory_iterator())
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestPstReader))
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestKernelReader))

    return test_suite

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    suites = getTestSuite()
    r = unittest.TextTestRunner(failfast=False)
    r.run(suites)
