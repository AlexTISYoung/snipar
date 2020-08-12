import numpy as np
import subprocess, sys, os.path
from itertools import *
import pandas as pd
import logging
from pysnptools.snpreader import SnpReader
from pysnptools.snpreader import SnpData
import numpy as np
import warnings
from pysnptools.pstreader import _OneShot
from six.moves import range


class Ped(_OneShot,SnpReader):
    '''
    A :class:`.SnpReader` for reading Ped-formated (and the related Map-formated) files from disk.

    See :class:`.SnpReader` for general examples of using SnpReaders.

    This format is described in http://zzz.bwh.harvard.edu/plink/data.shtml#ped and looks like::

         FAM001  1  0 0  1  2  A A  G G  A C 
         FAM001  2  0 0  1  2  A A  A G  0 0 
         ...

    the direction of the encoding from allele pair to 0,1,2 is arbitrary. That is, if the alleles are "A" and "G", then "G G" could be 0 and "A A" could be 2 or visa versa. The pair "A G" will always be 1.

    **Constructor:**
        :Parameters: * **filename** (*string*) -- The Ped file to read.
                     * **missing** (*string*) -- The value in the file that represents missing data.


        :Example:

        >>> from __future__ import print_function #Python 2 & 3 compatibility
        >>> from pysnptools.snpreader import Ped
        >>> data_on_disk = Ped('../examples/toydata.ped')
        >>> print((data_on_disk.iid_count, data_on_disk.sid_count))
        (500, 10000)

    **Methods beyond** :class:`.SnpReader`
    '''

    def __init__(self, filename, missing = '0'):
        '''
            filename    : string of the filename of the ped file
            missing         : string indicating a missing genotype (default '0')
        '''
        super(Ped, self).__init__()
        self.filename = SnpReader._name_of_other_file(filename,remove_suffix="ped", add_suffix="ped")
        self.missing = missing

    def _read_pstdata(self):
        col, col_property = SnpReader._read_map_or_bim(self.filename,remove_suffix="ped", add_suffix="map")
        ped = np.loadtxt(self.filename, dtype='str', comments=None)
        ped = ped.reshape(-1,ped.shape[-1]) #Turns 1-d row into 2-d
        row = ped[:,0:2]
        snpsstr = ped[:,6::]
        inan=snpsstr==self.missing
        snps = np.zeros((snpsstr.shape[0],snpsstr.shape[1]//2))
        for i in range(snpsstr.shape[1]//2):
            snps[inan[:,2*i],i]=np.nan
            vals=snpsstr[~inan[:,2*i],2*i:2*(i+1)]
            if vals.shape[0] > 0:
                snps[~inan[:,2*i],i]+=(vals==vals[0,0]).sum(1)
        snpdata = SnpData(iid=row,sid=col,pos=col_property,val=snps)
        return snpdata

    def copyinputs(self, copier):
        # doesn't need to self.run_once() because only uses original inputs
        copier.input(self.filename)
        copier.input(SnpReader._name_of_other_file(self.filename,remove_suffix="ped", add_suffix="map"))


    @staticmethod
    def write(filename, snpdata):
        """Writes a :class:`SnpData` to Ped format. The values must be 0,1,2. The direction of the encoding to allele pairs is arbitrary. This means
        that if a SnpData is written in Ped format and then read back, then 0's may become 2's and 2's may become 0's. (1's will stay 1's).
        Returns the :class:`.Ped`

        :param filename: the name of the file to create
        :type filename: string
        :param snpdata: The in-memory data that should be written to disk.
        :type snpdata: :class:`SnpData`
        :rtype: :class:`.Ped`

        >>> from pysnptools.snpreader import Ped, Bed
        >>> import pysnptools.util as pstutil
        >>> snpdata = Bed('../examples/toydata.bed',count_A1=False)[:,:10].read()  # Read first 10 snps from Bed format
        >>> pstutil.create_directory_if_necessary("tempdir/toydata10.ped")
        >>> Ped.write("tempdir/toydata10.ped",snpdata)            # Write data in Ped format
        Ped('tempdir/toydata10.ped')
        """

        if isinstance(filename,SnpData) and isinstance(snpdata,str): #For backwards compatibility, reverse inputs if necessary
            warnings.warn("write statement should have filename before data to write", DeprecationWarning)
            filename, snpdata = snpdata, filename 

        SnpReader._write_map_or_bim(snpdata, filename, remove_suffix="ped", add_suffix="map")

        # The PED file is a white-space (space or tab) delimited file: the first six columns are mandatory:
        # Family ID
        # Case ID
        # Paternal ID
        # Maternal ID
        # Sex (1=male; 2=female; other=unknown)
        # Phenotype

        pedfile = SnpReader._name_of_other_file(filename, remove_suffix="ped", add_suffix="ped")
        with open(pedfile,"w") as ped_filepointer:
            for iid_index, iid_row in enumerate(snpdata.iid):
                ped_filepointer.write("{0} {1} 0 0 0 0".format(iid_row[0],iid_row[1]))
                row = snpdata.val[iid_index,:]
                for sid_index, val in enumerate(row):
                    if val == 0:
                        s = "A A"
                    elif val == 1:
                        s = "A G"
                    elif val == 2:
                        s = "G G"
                    elif np.isnan(val):
                        s = "0 0"
                    else:
                        raise Exception("Expect values for ped file to be 0,1,2, or NAN. Instead, saw '{0}'".format(val))
                    ped_filepointer.write("\t"+s)
                ped_filepointer.write("\n")
        return Ped(filename)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
