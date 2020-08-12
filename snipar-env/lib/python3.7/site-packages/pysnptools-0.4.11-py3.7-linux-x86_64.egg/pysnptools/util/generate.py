from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import scipy as sp
import unittest
import logging
import sys
import os
import doctest
from six.moves import range
#from pysnptools.snpreader import SnpData, Bed, SnpNpz -- putting this here would cause a loop


def snp_gen(fst, dfr, iid_count, sid_count, maf_low=.05, maf_high=.5, seed=0, sibs_per_family=10, freq_pop_0=.5, chr_count=None, label_with_pop=False):
    """Generates a random :class:`.SnpData`

    :param fst: Degree of Population Structure, e.g. 0 (a special case), 0.005, 0.01, 0.05, 0.1
    :type fst: float

    :param dfr: Degree of Family Relatedness, the fraction of individuals belonging to a family, ie. fracSibs, e.g. 0.0, 0.5, 0.6, 0.7, 0.8, 0.9
    :type dfr: float

    :param iid_count: The number of individuals to generate. Because of rounding the actual number may be less.
    :type iid_count: int

    :param sid_count: The number of snps to generate.
    :type sid_count: int

    :param maf_low: (default .05) lower bound of uniformly-generated Minor allele frequency
    :type maf_low: float

    :param maf_high: (default .5) upper bound of uniformly-generated Minor allele frequency
    :type maf_high: float

    :param seed: (default 0) Random seed
    :type seed: int

    :param sibs_per_family: (default 10) number of siblings in each family
    :type sibs_per_family: int

    :param freq_pop_0: (default .5) Fraction of individuals in population 0 (the rest will be in population 1)
    :type freq_pop_0: float

    :param chr_count: (default one chromosome per SNP) Number of chromosomes to which SNPs should be assigned. The SNPs will
     be assigned as evenly as possible. Chromosome names are integers starting with 1. SNP positions within a chromosome are sequential
     integers starting with 1.
    :type chr_count: int

    :rtype: :class:`.SnpData`

    :Example:

    >>> from __future__ import print_function #Python 2 & 3 compatibility
    >>> snpdata = snp_gen(fst=.1,dfr=.5,iid_count=200,sid_count=20,maf_low=.05,seed=6)
    >>> print(int(snpdata.iid_count), int(snpdata.sid_count)) #because of rounding got 190 individuals
    190 20

    :Also See: :class:`.snpreader.SnpGen`

    """
    from pysnptools.snpreader import SnpData

    assert 0 <= freq_pop_0 and freq_pop_0 <=1.0,"assert 0 <= freq_pop_0 and freq_pop_0 <=1.0"

    if seed is not None:
        np.random.seed(int(seed % 2147483647)) #old maxint

    iid_solo_count = iid_count-iid_count*dfr
    family_count = int(iid_count*dfr/(2 * sibs_per_family))

    ancestral = np.random.uniform(maf_low, maf_high, sid_count)     #sample ancestral allele frequencies

    snp_list=[]
    for population_index, freq_pop in enumerate([freq_pop_0, 1.0-freq_pop_0]):
        logging.info("Simulating SNPs from a population %i" % population_index)
        snps_parents=_generate_snps(ancestral, fst, int(iid_solo_count*freq_pop), sid_count)
        snp_list.append(snps_parents)
        snp_list.append(_generate_kids(parent_snps=snps_parents, family_count=int(freq_pop*family_count), sibs_per_family=sibs_per_family))

    snp_list.append(_generate_kids(parent_snps=np.concatenate(snp_list), family_count=family_count, sibs_per_family=sibs_per_family))
    val = np.concatenate(snp_list)

    if not label_with_pop:
        iid = np.array([["i_{0}".format(iid_index),"f_{0}".format(iid_index)] for iid_index in range(val.shape[0])],dtype=str).reshape(-1,2)
    else:
        assert len(snp_list) == 5, "real assert"
        iid0 = [["0",str(iid_index)] for iid_index in range(len(snp_list[0])+len(snp_list[1]))] #parents and children of pop 0
        iid1 = [["1",str(iid_index)] for iid_index in range(len(snp_list[2])+len(snp_list[3]))] #parents and children of pop 1
        iid2 = [["2",str(iid_index)] for iid_index in range(len(snp_list[4]))]                  #children with parents in any pop
        iid = np.array(iid0+iid1+iid2,dtype=str).reshape(-1,2)

    sid = np.array(["snp_{0}".format(sid_index) for sid_index in range(val.shape[1])],dtype='str')

    if chr_count is None:
        chr_count = len(sid)

    assert len(sid) == 0 or chr_count > 0, "chr_count must be at least 1 (unless sid_count is 0)"
    sid_per_chrom = int(sp.ceil(float(len(sid))/max(1,chr_count)))
    pos = np.array(list([1+sid_index//sid_per_chrom, 1+sid_index%sid_per_chrom, 1+sid_index%sid_per_chrom] for sid_index in range(len(sid))))
    if len(sid) == 0: #make it work when no sids are wanted
        pos = pos.reshape(len(sid),3)

    snpdata = SnpData(iid=iid, sid=sid, val=val, pos=pos,
                      parent_string="snp_gen(fst={0}, dfr={1}, iid_count={2}, sid_count={3}, maf_low={4}, maf_high={5}, seed={6}, sibs_per_family={7}, freq_pop_0={8})"
                      .format(fst, dfr, iid_count, sid_count, maf_low, maf_high, seed, sibs_per_family, freq_pop_0)
                      )

    if snpdata.iid_count != iid_count:
        logging.warn("Because of rounding the actual number of iids is {0} rather than the requested {1}".format(snpdata.iid_count, iid_count))

    return snpdata

def _beta_mode(a,f): #!!!
    alpha,beta = ((a*(1.0-f)/f,(1.0-a)*(1.0-f)/f))
    print(alpha,beta)
    mean = alpha/(alpha+beta)
    var = alpha * beta / ((alpha+beta)**2 * (alpha+beta+1))
    return mean, var**.5

def _generate_snps(ancestral, fst, iid_count, sid_count):
    """
    Generates genotypes with a certain MAF and optionally with population structure.
    In case of no population structure, they are sampled from a binomial,
    otherwise from a Beta-Binomial (Balding and Nichols, 1995).
    """
    if fst == 0.0: 
        alpha = ancestral #special treatment if no population structure
    else:
        alpha = np.random.beta(ancestral*(1.0-fst)/fst,(1.0-ancestral)*(1.0-fst)/fst, sid_count)


    #generate from population frequencies    
    snps = np.zeros((iid_count,sid_count),dtype='int8') #.zeros not .empty because will be adding 1's to it
    for i in range(2): #"2" for diploid
        #sample each allele
        rand = np.random.random((iid_count,sid_count))
        snps[rand<alpha]+=1
    return snps


def _generate_kids(parent_snps, family_count, sibs_per_family): #!!! should it be sibs, kids, or children
    '''
    generate a single set of family members
    '''    
    parent_count, sid_count = parent_snps.shape
    assert parent_count>=2*family_count, "parent_count>=2*family_count"


    parent_permutation = np.random.permutation(parent_count)
    snps = np.zeros((family_count*sibs_per_family,sid_count),dtype=np.float64)
    for copy_index in range(2):#"2" for diploid
        sample = parent_snps[parent_permutation[copy_index*family_count:(copy_index+1)*family_count],:]         #sample each allele
        for kid_index in range(sibs_per_family):
            rand = np.random.random((family_count,sid_count))
            snps[kid_index*family_count:(kid_index+1)*family_count][rand<0.5*sample]+=1
    return snps




def _encode_snp(entry):
    if entry == 0:
        return "A A"
    elif entry == 1:
        return "A C"
    elif entry == 2:
        return "C C"


def _generate_phenotype(snp_data, causals, genetic_var, noise_var, seed=None):
    """
    generate phenotype given genotype

    'causals' can be either an array of indexes to the causal snps or the number of causal snps desired.
    """

    if seed is not None:
        np.random.seed(int(seed % 2147483647)) #!!!Better to use numpy.random.RandomState instead (look for other places in code, too) #old maxint
    
    try:
        num_causal = len(causals)
        causal_idx = causals
    except:
        num_causal = causals
        #the "if..else" is a work around because the linux version of np.random.choice doesn't like to select zero items from an empty list. We need to call random in either case so that the random seed ends up in the expected state
        causal_idx = np.random.choice(sp.arange(snp_data.sid_count if num_causal>0 else 1),size=num_causal,replace=False)

    num_phenotypes = 1
    mean = 0.0
    X = snp_data.val.copy()
    X.flags.writeable = False
    X_causal = X[:,causal_idx]
    X_causal = 1./np.sqrt(X_causal.shape[1]) * X_causal
    W = np.random.randn(num_causal, num_phenotypes) * np.sqrt(genetic_var) + mean #Weight of each causal SNP
    XW = np.dot(X_causal, W)
    noise_std = np.sqrt(noise_var)
    Z = noise_std*sp.randn(X_causal.shape[0], num_phenotypes)
    y = XW + Z
    y = y[:,0]

    return y


class TestGenerate(unittest.TestCase):     


    @classmethod
    def setUpClass(self):
        self.currentFolder = os.path.dirname(os.path.realpath(__file__))

    def gen_and_compare(self, output_file, **kwargs):
        from pysnptools.snpreader import Bed

        gen_snpdata = snp_gen(**kwargs)
        #pstutil.create_directory_if_necessary(self.currentFolder + "/tempdir/" + output_file,isfile=True)
        #Bed.write(gen_snpdata, self.currentFolder + "/tempdir/" + output_file)  #comment out
        ref_snpdata = Bed(self.currentFolder + "/../../tests/datasets/generate/" + output_file,count_A1=False).read()
        assert gen_snpdata == ref_snpdata, "Failure on "+output_file
        return gen_snpdata
        #!!! Ped doesn't seem to round trip well
        #!!! Hdf5 doesn't seem to round trip well


    def test_gen1(self):
        self.gen_and_compare("gen1", fst=0,dfr=.5,iid_count=200,sid_count=20,maf_low=.05,seed=5)

    def test_gen2(self):
        self.gen_and_compare("gen2", fst=.1,dfr=.5,iid_count=200,sid_count=20,maf_low=.05,seed=5)

    def test_gen2b(self):
        """
        Test that different seed produces different result
        """
        from pysnptools.snpreader import Bed
        gen_snpdata = self.gen_and_compare("gen2b", fst=.1,dfr=.5,iid_count=200,sid_count=20,maf_low=.05,seed=6)
        ref_snpdata = Bed(self.currentFolder + "/../../tests/datasets/generate/gen2",count_A1=False).read()
        assert gen_snpdata != ref_snpdata, "Expect different seeds to produce different results"

    def test_gen3(self):
        self.gen_and_compare("gen3", fst=.1,dfr=0,iid_count=200,sid_count=20,maf_low=.05,seed=5)

    def test_gen4(self):
        self.gen_and_compare("gen4", fst=.1,dfr=.01,iid_count=200,sid_count=20,maf_low=.1,seed=5)

    def test_gen5(self):
        from pysnptools.snpreader import Bed
        gen_snpdata = self.gen_and_compare("gen5", fst=.1,dfr=.5,iid_count=200,sid_count=20,maf_low=.05,maf_high=.4, seed=5)
        ref_snpdata = Bed(self.currentFolder + "/../../tests/datasets/generate/gen2",count_A1=False).read()
        assert gen_snpdata != ref_snpdata, "Expect different seeds to produce different results"

    def test_gen6(self):
        from pysnptools.snpreader import Bed
        gen_snpdata = self.gen_and_compare("gen6", fst=.1,dfr=.5,iid_count=200,sid_count=20,maf_low=.05,seed=5,sibs_per_family=5)
        ref_snpdata = Bed(self.currentFolder + "/../../tests/datasets/generate/gen2",count_A1=False).read()
        assert gen_snpdata != ref_snpdata, "Expect different seeds to produce different results"

    def test_gen7(self):
        from pysnptools.snpreader import Bed
        gen_snpdata = self.gen_and_compare("gen7", fst=.1,dfr=.5,iid_count=200,sid_count=20,maf_low=.05,seed=5,freq_pop_0=.75)
        ref_snpdata = Bed(self.currentFolder + "/../../tests/datasets/generate/gen2",count_A1=False).read()
        assert gen_snpdata != ref_snpdata


    def test_gen8(self):
        self.gen_and_compare("gen8a", fst=.1,dfr=.5,iid_count=200,sid_count=20,maf_low=.05,seed=5,chr_count=3)
        self.gen_and_compare("gen8b", fst=.1,dfr=.5,iid_count=200,sid_count=20,maf_low=.05,seed=5,chr_count=4)
        self.gen_and_compare("gen8c", fst=.1,dfr=.5,iid_count=200,sid_count=20,maf_low=.05,seed=5,chr_count=6)

    def test_pheno1(self):
        from pysnptools.snpreader import Bed, SnpData, SnpNpz
        some_snp_data = Bed(self.currentFolder + "/../../tests/datasets/generate/gen2",count_A1=False).read()
        gen_snpdata = SnpData(iid=some_snp_data.iid,sid=["pheno"],val=_generate_phenotype(some_snp_data, 10, genetic_var=.5, noise_var=.5, seed=5).reshape(-1,1))
        #SnpNpz.write(r'c:\deldir\pheno1.snp.npz',gen_snpdata)
        ref_snpdata = SnpNpz(self.currentFolder + "/../../tests/datasets/generate/pheno1.snp.npz").read()
        assert gen_snpdata == ref_snpdata


    def test_gensmall(self):
        #Just checking that doesn't generate errors
        for iid_count in [10, 5, 3, 2, 1, 0]:
            for sid_count in [0, 10, 5, 3, 2, 1]:
                for chr_count in [30, 10, 5, 3, 2, 1, 0]:
                    if chr_count == 0 and sid_count > 0:
                        continue # not break
                    logging.debug("{0}, {1}, {2}".format(iid_count, sid_count, chr_count))
                    snpdata = snp_gen(fst=.1,dfr=.5,iid_count=iid_count,sid_count=sid_count,maf_low=.05,seed=6,chr_count=chr_count)
                    assert snpdata.iid_count <= iid_count
                    assert snpdata.sid_count == sid_count
                    assert len(snpdata.pos) == 0 or max(snpdata.pos[:,0]) <= chr_count
                    assert len(snpdata.pos) == 0 or max(snpdata.pos[:,1]) <= int(max(1,np.ceil(float(sid_count) / chr_count)))
                    assert len(snpdata.pos) == 0 or max(snpdata.pos[:,2]) <= int(max(1,np.ceil(float(sid_count) / chr_count)))

    def test_doc_test(self):
        import sys
        import pysnptools.util.generate
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__))+"/..")
        result = doctest.testmod(pysnptools.util.generate)
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__


def getTestSuite():
    """
    set up test suite
    """
    
    test_suite = unittest.TestSuite([])
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestGenerate))

    return test_suite


if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)

    suites = getTestSuite()
    r = unittest.TextTestRunner(failfast=False)
    r.run(suites)

    import doctest
    doctest.testmod()


    print("done")

