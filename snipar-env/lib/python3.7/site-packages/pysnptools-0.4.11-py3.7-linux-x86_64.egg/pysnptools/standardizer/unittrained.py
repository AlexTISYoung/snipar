import numpy as np
import logging
import warnings
from pysnptools.standardizer import Standardizer

class UnitTrained(Standardizer):
    """A :class:`.Standardizer` to unit standardize one set of SNP data based on the mean and stddev of another set of SNP data.
    NaN values are then filled with zero.

    See :class:`.Standardizer` for more information about standardization.

    **Constructor:**
        :Parameters: * **stats** (*ndarray of float*) -- The mean and stddev of each sid

    >>> from __future__ import print_function #Python 2 & 3 compatibility
    >>> from pysnptools.standardizer import Unit
    >>> from pysnptools.snpreader import Bed
    >>> train = Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False)[1:,:].read() # read SNP values for all but the first iid
    >>> _, unittrained = train.standardize(Unit(),return_trained=True) #Unit standardize and remember the mean and stddev of each sid
    >>> print(unittrained.stats[:5,:]) #Print the means and stddev of the first five sids
    [[...1.94983278  0.21828988]
     [...1.96989967  0.17086341]
     [...1.84280936  0.39057474]
     [...1.99665552  0.0577347 ]
     [...1.97658863  0.15120608]]
    >>> test = Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False)[0,:].read() # read SNP values for the first iid
    >>> test = test.standardize(unittrained) # Use the mean and stddev of the train data to unit standardize the test data.
    >>> print('{0:.6f}'.format(test.val[0,0]))
    0.229819
    """

    #!!might want to add code so that can check that sids are in the same order for both test and train
    def __init__(self, sid, stats):
        super(UnitTrained, self).__init__()
        self.sid = sid
        self.stats = stats
        self.sid_to_index = None

    def __repr__(self): 
        return "{0}(stats={1},sid={2})".format(self.__class__.__name__,self.stats,self.sid)

    @property
    def is_constant(self):
        return True        

    def standardize(self, snps, block_size=None, return_trained=False, force_python_only=False):
        if block_size is not None:
            warnings.warn("block_size is deprecated (and not needed, since standardization is in-place", DeprecationWarning)

        if hasattr(snps,"val"):
            val = snps.val
            if len(self.sid) == len(snps.sid) and np.array_equal(self.sid,snps.sid):
                stats = self.stats
            else:
                if self.sid_to_index is None:
                    self.sid_to_index = {sid:index for index,sid in enumerate(self.sid)}
                stats = np.array([self.stats[self.sid_to_index[sid]] for sid in snps.sid])
        else:
            warnings.warn("standardizing an nparray instead of a SnpData is deprecated", DeprecationWarning)
            val = snps
            stats = self.stats

        self._standardize_unit_and_beta(val, is_beta=False, a=np.nan, b=np.nan, apply_in_place=True,use_stats=True,stats=stats,force_python_only=force_python_only)

        if return_trained:
            return snps, self
        else:
            return snps

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS|doctest.NORMALIZE_WHITESPACE)
