from pysnptools.snpreader import SnpReader
from pysnptools.pstreader._subset import _PstSubset

class _SnpSubset(_PstSubset,SnpReader):
    def __init__(self, *args, **kwargs):
        super(_SnpSubset, self).__init__(*args, **kwargs)
