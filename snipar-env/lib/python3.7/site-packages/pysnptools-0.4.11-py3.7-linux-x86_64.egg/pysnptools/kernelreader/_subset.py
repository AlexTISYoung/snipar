from pysnptools.kernelreader import KernelReader
from pysnptools.pstreader._subset import _PstSubset

class _KernelSubset(KernelReader,_PstSubset):
    def __init__(self, *args, **kwargs):
        super(_KernelSubset, self).__init__(*args, **kwargs)

