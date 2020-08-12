import numpy as np
import logging
import warnings
from pysnptools.standardizer import Standardizer


class Beta(Standardizer):
    '''
    A :class:`.Standardizer` to beta standardize SNP data.

    See :class:`.Standardizer` for more information about standardization.

    **Constructor:**
        :Parameters: * **a** (*float*) -- The *a* parameter of the beta distribution
                     * **b** (*float*) -- The *b* parameter of the beta distribution

    >>> from __future__ import print_function #Python 2 & 3 compatibility
    >>> from pysnptools.standardizer import Beta
    >>> from pysnptools.snpreader import Bed
    >>> snpdata1 = Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False).read().standardize(Beta(1,25))
    >>> print('{0:.6f}'.format(snpdata1.val[0,0]))
    0.680802
    '''
    def __init__(self,a,b):
        super(Beta, self).__init__()
        self.a = a
        self.b = b

    def __repr__(self): 
        return "{0}(a={1},b={2})".format(self.__class__.__name__,self.a,self.b)

    def standardize(self, snpdata, block_size=None, return_trained=False, force_python_only=False): #!!!later why is the 2nd argument called 'snpdata' here, but 'snps' in unit.py?
        if block_size is not None:
            warnings.warn("block_size is deprecated (and not needed, since standardization is in-place", DeprecationWarning)

        if hasattr(snpdata,"val"):
            val = snpdata.val
        else:
            warnings.warn("standardizing an nparray instead of a SnpData is deprecated", DeprecationWarning)
            val = snpdata

        stats = self._standardize_unit_and_beta(val, is_beta=True, a=self.a, b=self.b, apply_in_place=True, use_stats=False,stats=None,force_python_only=force_python_only)
        if return_trained:
            from pysnptools.standardizer import BetaTrained
            assert hasattr(snpdata,"val"), "return_trained=True must be used with SnpData"
            return snpdata, BetaTrained(self.a,self.b,snpdata.sid,stats)
        else:
            return snpdata

    def _merge_trained(self, trained_list):
        from pysnptools.standardizer import BetaTrained

        sid = np.concatenate([trained.sid for trained in trained_list])
        stats = np.concatenate([trained.stats for trained in trained_list])
        a_set = set([trained.a for trained in trained_list])
        b_set = set([trained.b for trained in trained_list])
        assert len(a_set) <= 1,"Expect all BetaTrained's to have the same 'a'"
        assert len(b_set) <= 1,"Expect all BetaTrained's to have the same 'b'"
        a = list(a_set)[0] if a_set else None
        b = list(b_set)[0] if b_set else None
        return BetaTrained(a, b, sid, stats)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
