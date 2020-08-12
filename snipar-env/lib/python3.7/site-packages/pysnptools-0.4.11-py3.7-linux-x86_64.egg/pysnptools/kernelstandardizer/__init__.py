
#!!!move these to their own files

import numpy as np
import logging

class KernelStandardizer(object):
    '''
    A KernelStandardizer is a class such as :class:`.DiagKtoN` and :class:`.Identity` to be used by the :meth:`.KernelData.standardize` to standardize Kernel data.
    It always works in-place *and* returns the :class:`.KernelData` on which it works.

    :Example:

    Read and standardize Kernel data.

    >>> from __future__ import print_function #Python 2 & 3 compatibility
    >>> from six.moves import range #Python 2 & 3 compatibility
    >>> from pysnptools.kernelstandardizer import DiagKtoN
    >>> from pysnptools.kernelreader import KernelNpz
    >>> kerneldata1 = KernelNpz('../examples/toydata.kernel.npz').read()
    >>> print(np.diag(kerneldata1.val).sum())
    5000000.0
    >>> kerneldata1 = kerneldata1.standardize(DiagKtoN())
    >>> print(np.diag(kerneldata1.val).sum())
    500.0

    Can also return a constant kernel standardizer that be applied to other :class:`.KernelData`.

    >>> kernel_whole = KernelNpz('../examples/toydata.kernel.npz')
    >>> train_idx, test_idx = range(10,kernel_whole.iid_count), range(0,10)  #test on the first 10, train on the rest
    >>> kernel_train, trained_standardizer = DiagKtoN().standardize(kernel_whole[train_idx].read(),return_trained=True)
    >>> print('{0:.6f}'.format(np.diag(kernel_train.val).sum()))
    490.000000
    >>> print('{0:.6f}'.format(trained_standardizer.factor))
    0.000100
    >>> kernel_whole_test = kernel_whole[:,test_idx].read().standardize(trained_standardizer)
    >>> print('{0:.6f}'.format(kernel_whole_test.val[0,0]))
    0.992217


    Details of Methods & Properties:
    '''
    def standardize(self, kerneldata, return_trained=False, force_python_only=False):
        '''
        Applies standardization, in place, to :class:`.KernelData`. For convenience also returns the :class:`KernelData`.

        :param snps: kernel values to standardize
        :type snps: :class:`.KernelData`

        :param return_trained: If true, returns a second value containing a constant :class:`.KernelStandardizer` trained on this data.
        :type return_trained: bool

        :param force_python_only: optional -- If False (default), may use outside library code. If True, requests that the read
            be done without outside library code.
        :type force_python_only: bool

        :rtype: :class:`.KernelData`, (optional) constant :class:`.KernelStandardizer`

        '''
        raise NotImplementedError("subclass {0} needs to implement method '.standardize'".format(self.__class__.__name__))

class Identity(KernelStandardizer):
    '''
    A :class:`.KernelStandardizer` that does nothing to kernel data.

    See :class:`.KernelStandardizer` for more information about standardization.

    >>> from pysnptools.kernelstandardizer import Identity as KS_Identity
    >>> from pysnptools.kernelreader import KernelNpz
    >>> kerneldata1 = KernelNpz('../examples/toydata.kernel.npz').read()
    >>> print(np.diag(kerneldata1.val).sum())
    5000000.0
    >>> kerneldata1 = kerneldata1.standardize(KS_Identity())
    >>> print(np.diag(kerneldata1.val).sum())
    5000000.0
    '''

    def __init__(self):
        super(Identity, self).__init__()

    def standardize(self, kerneldata, return_trained=False, force_python_only=False):
        if return_trained:
            return kerneldata, self
        else:
            return kerneldata

    def __repr__(self): 
        return "{0}()".format(self.__class__.__name__)

from pysnptools.standardizer import DiagKtoN #as SN_DiagKtoN
from pysnptools.standardizer import DiagKtoNTrained #as SN_DiagKtoNTrained



if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()

