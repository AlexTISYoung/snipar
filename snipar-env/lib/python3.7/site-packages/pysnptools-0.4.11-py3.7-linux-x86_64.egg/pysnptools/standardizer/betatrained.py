import numpy as np
import logging
import warnings
from pysnptools.standardizer import Standardizer

class BetaTrained(Standardizer):
    """A :class:`.Standardizer` to beta standardize one set of SNP data based on the mean of another set of SNP data

    See :class:`.Standardizer` for more information about standardization.

    **Constructor:**
        :Parameters: * **a** (*float*) -- The *a* parameter of the beta distribution
                     * **b** (*float*) -- The *b* parameter of the beta distribution
                     * **stats** (*ndarray of float*) -- The mean and stddev of each sid

    >>> from __future__ import print_function #Python 2 & 3 compatibility
    >>> from pysnptools.standardizer import Beta
    >>> from pysnptools.snpreader import Bed
    >>> train = Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False)[1:,:].read() # read SNP values for all but the first iid
    >>> _, betatrained = train.standardize(Beta(1,25),return_trained=True) #beta standardize and remember the mean and stddev of each sid
    >>> print(betatrained.stats[:5,:]) #Print the means and stddev of the first five sids
    [[...1.94983278  0.21828988]
     [...1.96989967  0.17086341]
     [...1.84280936  0.39057474]
     [...1.99665552  0.0577347 ]
     [...1.97658863  0.15120608]]
    >>> test = Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False)[0,:].read() # read SNP values for the first iid
    >>> test = test.standardize(betatrained) # Use the mean of the train data to beta standardize the test data.
    >>> print('{0:.6f}'.format(test.val[0,0]))
    0.681674
    """

    def __init__(self, a,b,sid,stats):
        super(BetaTrained, self).__init__()
        self.a=a
        self.b=b
        self.sid=sid
        self.stats=stats

    def __repr__(self): 
        return "{0}(a={1},b={2},stats={3},sid={4})".format(self.__class__.__name__,self.a,self.b,self.stats,self.sid)

    @property
    def is_constant(self):
        return True        

    def standardize(self, snps, block_size=None, return_trained=False, force_python_only=False):
        if block_size is not None:
            warnings.warn("block_size is deprecated (and not needed, since standardization is in-place", DeprecationWarning)

        if hasattr(snps,"val"):
            val = snps.val
            assert np.array_equal(self.sid,snps.sid), "sid in training and use must be the same and in the same order"
        else:
            warnings.warn("standardizing an nparray instead of a SnpData is deprecated", DeprecationWarning)
            val = snps

        self._standardize_unit_and_beta(val, is_beta=True, a=self.a, b=self.b, apply_in_place=True,use_stats=True,stats=self.stats,force_python_only=force_python_only)
        if return_trained:
            return snps, self
        else:
            return snps

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS|doctest.NORMALIZE_WHITESPACE)
