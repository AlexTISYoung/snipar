import numpy as np
from pysnptools.standardizer import Standardizer
import logging
import warnings


class DiagKtoN(Standardizer):
    '''
    Both a :class:`.Standardizer` and a A :class:`.KernelStandardizer`.

    When applied to a :class:`.SnpData`, it multiplies the SNP values by the square root of a factor such that a kernel
    constructed from the SNP values will have a diagonal that sums to iid_count. This class thus
    standardizes the kernel before it is even constructed.

    When applied to a :class:`.KernelData`, it multiplies the kernel values by the a factor such that a kernel
    will have a diagonal that sums to iid_count.

    See :class:`.Standardizer` for more information about standardization.

    **Constructor:**
        :Parameters: * **deprecated_iid_count** (*int*) (optional) -- Deprecated.


    Example of DiagKtoN to :class:`.SnpData`:

    >>> from __future__ import print_function #Python 2 & 3 compatibility
    >>> from pysnptools.standardizer import DiagKtoN, Unit, Identity
    >>> from pysnptools.snpreader import Bed
    >>> snpdata1 = Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False).read().standardize(Unit()).standardize(DiagKtoN())
    >>> kernel1 = snpdata1.read_kernel(Identity())
    >>> print('{0:.6f}'.format(np.diag(kernel1.val).sum()))
    300.000000

    Example of DiagKtoN to :class:`.KernelData`:

    >>> import numpy as np
    >>> from pysnptools.standardizer import DiagKtoN, Unit, Identity
    >>> from pysnptools.snpreader import Bed
    >>> snpdata1 = Bed('../../tests/datasets/all_chr.maf0.001.N300',count_A1=False).read().standardize(Unit())
    >>> kernel1 = snpdata1.read_kernel(DiagKtoN(),block_size=None)
    >>> print('{0:.6f}'.format(np.diag(kernel1.val).sum()))
    300.000000

    '''
    """diag(K)=N standardization of the data"""
    def __init__(self, deprecated_iid_count=None):
        super(DiagKtoN, self).__init__()
        if deprecated_iid_count is not None:
            warnings.warn("'iid_count' is deprecated (and not needed, since can get iid_count from SNPs val's first dimension", DeprecationWarning)

    def _standardize_kernel(self, kerneldata, return_trained=False, force_python_only=False):

        factor = float(kerneldata.iid_count) / np.diag(kerneldata.val).sum()

        if abs(factor-1.0)>1e-15:
            kerneldata.val *= factor

        if return_trained:
            return kerneldata, DiagKtoNTrained(factor)
        else:
            return kerneldata

    def standardize(self, input, block_size=None, return_trained=False, force_python_only=False):
        from pysnptools.kernelreader import KernelData
        if block_size is not None:
            warnings.warn("block_size is deprecated (and not needed, since standardization is in-place", DeprecationWarning)
        if isinstance(input,KernelData):
            return self._standardize_kernel(input, return_trained=return_trained,force_python_only=force_python_only)
        else:
            return self._standardize_snps(input, return_trained=return_trained,force_python_only=force_python_only)

    def _standardize_snps(self, snps, return_trained=False, force_python_only=False):

        if hasattr(snps,"val"):
            val = snps.val
        else:
            warnings.warn("standardizing an nparray instead of a SnpData is deprecated", DeprecationWarning)
            val = snps

        vec = val.reshape(-1, order="A")
        # make sure no copy was made
        assert not vec.flags['OWNDATA']
        squared_sum = vec.dot(vec)
        factor = float(val.shape[0])/squared_sum

        if abs(factor-1.0)>1e-15:
            val *= np.sqrt(factor)

        if return_trained:
            return snps, DiagKtoNTrained(factor)
        else:
            return snps

    def __repr__(self): 
        return "{0}()".format(self.__class__.__name__)

#!!!move to its own file
class DiagKtoNTrained(Standardizer):
    """Both a :class:`.Standardizer` and a A :class:`.KernelStandardizer`.

    When applied to a :class:`.SnpData`, DiagKtoN standardizes one set of SNP data based on another set of SNP data.

    When applied to a :class:`.KernelData`, DiagKtoN standardizes one kernel data based on another kernel data.

    See :class:`.Standardizer` for more information about standardization.

    **Constructor:**
        :Parameters: * **factor** (*float*) -- The number what when multiplied into the kernel values will make the diagonal of the kernel values to  sum to iid_count. 

    """
    def __init__(self,factor):
        super(DiagKtoNTrained, self).__init__()
        self.factor = factor

    def standardize(self, input, block_size=None, return_trained=False, force_python_only=False):
        if block_size is not None:
            warnings.warn("block_size is deprecated (and not needed, since standardization is in-place", DeprecationWarning)

        from pysnptools.kernelreader import KernelData
        if isinstance(input,KernelData):
            return self._standardize_kernel(input, return_trained=return_trained,force_python_only=force_python_only)
        else:
            return self._standardize_snps(input, return_trained=return_trained,force_python_only=force_python_only)

    @property
    def is_constant(self):
        return abs(self.factor-1.0)<1e-15

    def _standardize_snps(self, snps, return_trained=False, force_python_only=False):
    
        if hasattr(snps,"val"):
            val = snps.val
        else:
            warnings.warn("standardizing an nparray instead of a SnpData is deprecated", DeprecationWarning)
            val = snps

        if not self.is_constant:
            val *= np.sqrt(self.factor)

        if return_trained:
            return snps, self
        else:
            return snps

    def _standardize_kernel(self, kerneldata, return_trained=False, force_python_only=False):
        if not self.is_constant:
            kerneldata.val *= self.factor

        if return_trained:
            return kerneldata, self
        else:
            return kerneldata


    def __repr__(self): 
        return "{0}({1})".format(self.__class__.__name__,self.factor)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
        