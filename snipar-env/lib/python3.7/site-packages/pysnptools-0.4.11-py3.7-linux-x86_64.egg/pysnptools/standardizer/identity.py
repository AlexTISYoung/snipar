import numpy as np
import logging
from pysnptools.standardizer import Standardizer

class Identity(Standardizer):
    '''
    A :class:`.Standardizer` that does nothing to SNP data.

    See :class:`.Standardizer` for more information about standardization.

    >>> from __future__ import print_function #Python 2 & 3 compatibility
    >>> from pysnptools.standardizer import Identity
    >>> from pysnptools.snpreader import Bed
    >>> snpdata1 = Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False).read()
    >>> print(snpdata1.val[0,0])
    2.0
    >>> snpdata1 = snpdata1.standardize(Identity())
    >>> print(snpdata1.val[0,0])
    2.0
    '''

    def __init__(self):
        super(Identity, self).__init__()

    def standardize(self, snps, block_size=None, return_trained=False, force_python_only=False):
        if block_size is not None:
            warnings.warn("block_size is deprecated (and not needed, since standardization is in-place", DeprecationWarning)
        if return_trained:
            return snps, self
        else:
            return snps

    @property
    def is_constant(self):
        return True        

    def __repr__(self): 
        return "{0}()".format(self.__class__.__name__)

    def _merge_trained(self, trained_list):
        return self

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
