import numpy as np
import logging
from pysnptools.standardizer import Standardizer
import warnings

class BySqrtSidCount(Standardizer):
    '''
    '''

    def __init__(self):
        super(BySqrtSidCount, self).__init__()
        warnings.warn("BySqrtSidCount no longer supported", DeprecationWarning)


