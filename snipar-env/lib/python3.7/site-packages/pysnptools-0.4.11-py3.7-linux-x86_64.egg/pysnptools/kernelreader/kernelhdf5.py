import logging
from pysnptools.kernelreader import KernelReader
from pysnptools.pstreader import PstHdf5
import warnings

class KernelHdf5(PstHdf5,KernelReader):
    '''
    A :class:`.KernelReader` for reading \*.kernel.hdf5 files from disk.

    See :class:`.KernelReader` for general examples of using KernelReaders.

    The general HDF5 format is described in http://www.hdfgroup.org/HDF5/. The KernelHdf5 format stores
    val, iid, sid, and pos information in Hdf5 format.
   
    **Constructor:**
        :Parameters: * **filename** (*string*) -- The KernelHdf5 file to read.

        :Example:

        >>> from __future__ import print_function #Python 2 & 3 compatibility
        >>> from pysnptools.kernelreader import KernelHdf5
        >>> data_on_disk = KernelHdf5('../examples/toydata.kernel.hdf5')
        >>> print(data_on_disk.iid_count)
        500

    **Methods beyond** :class:`.KernelReader`

    '''
    def __init__(self,*args, **kwargs):
        super(KernelHdf5, self).__init__(*args, **kwargs)

    @staticmethod
    def write(filename, kerneldata, hdf5_dtype=None, sid_major=True):
        """Writes a :class:`KernelData` to KernelHdf5 format and returns the :class:`.KernelHdf5`.

        :param filename: the name of the file to create
        :type filename: string
        :param kerneldata: The in-memory data that should be written to disk.
        :type kerneldata: :class:`KernelData`
        :param hdf5_dtype: None (use the .val's dtype) or a Hdf5 dtype, e.g. 'f8','f4',etc.
        :type hdf5_dtype: string
        :param col_major: Tells if vals should be stored on disk in sid_major (default) or iid_major format.
        :type col_major: bool
        :rtype: :class:`.KernelHdf5`

        >>> from pysnptools.snpreader import Bed
        >>> from pysnptools.standardizer import Unit
        >>> import pysnptools.util as pstutil
        >>> kerneldata = Bed('../examples/toydata.bed',count_A1=False).read_kernel(Unit())     # Create a kernel from the data in the Bed file
        >>> pstutil.create_directory_if_necessary("tempdir/toydata.kernel.hdf5")
        >>> KernelHdf5.write("tempdir/toydata.kernel.hdf5",kerneldata)          # Write data in KernelHdf5 format
        KernelHdf5('tempdir/toydata.kernel.hdf5')
        """
        PstHdf5.write(filename,kerneldata,hdf5_dtype=hdf5_dtype,col_major=sid_major)
        return KernelHdf5(filename)


if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
