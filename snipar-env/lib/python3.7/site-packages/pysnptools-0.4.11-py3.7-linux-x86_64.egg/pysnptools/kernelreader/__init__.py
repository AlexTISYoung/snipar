"""Tools for reading and manipulating kernels. A kernels is a matrix from iid to iid that typically tells the relatedness of pairs of people.
"""

from pysnptools.kernelreader.kernelreader import KernelReader
from pysnptools.kernelreader.kerneldata import KernelData
from pysnptools.kernelreader.snpkernel import SnpKernel
from pysnptools.kernelreader.identity import Identity
from pysnptools.kernelreader.kernelnpz import KernelNpz
from pysnptools.kernelreader.kernelhdf5 import KernelHdf5