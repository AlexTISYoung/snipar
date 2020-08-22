import numpy as np
import subprocess, sys, os.path
from itertools import *
import pandas as pd
import logging
from pysnptools.snpreader import SnpReader
from pysnptools.snpreader import SnpData
import warnings
from pysnptools.pstreader import _OneShot

class Dat(_OneShot,SnpReader):
    '''
    A :class:`.SnpReader` for reading Dat/Fam/Map-formated files from disk.

    See :class:`.SnpReader` for general examples of using SnpReaders.

    This is a text format that can store any numeric values. (In contrast, Bed and Ped can only store 0,1,2, and missing). Its Dat files look like::
    
        null_0  	j	n	0.333	1	2
        null_100	j	n	2	1	1
        null_200	j	n	0	nan	1
        ...

    Its Map and Fam files are described in http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml.

    **Constructor:**
        :Parameters: * **filename** (*string*) -- The Dat file to read.

        :Example:

        >>> from __future__ import print_function #Python 2 & 3 compatibility
        >>> from pysnptools.snpreader import Dat
        >>> data_on_disk = Dat('../examples/toydata.dat')
        >>> print((data_on_disk.iid_count, data_on_disk.sid_count))
        (500, 10000)

    **Methods beyond** :class:`.SnpReader`
    '''

    def __init__(self, filename):
        '''
        filename    : string of the name of the Dat file.
        '''
        super(Dat, self).__init__()
        self.filename = SnpReader._name_of_other_file(filename,remove_suffix="dat", add_suffix="dat")

    def _read_pstdata(self):
        row = SnpReader._read_fam(self.filename,remove_suffix="dat")
        col, col_property = SnpReader._read_map_or_bim(self.filename,remove_suffix="dat", add_suffix="map")
        if len(row)==0 or len(col)==0:
            return SnpData(iid=row,sid=col,pos=col_property,val=np.empty([len(row),len(col)]))
        datfields = pd.read_csv(self.filename,delimiter = '\t',header=None,index_col=False)
        if not np.array_equal(datfields[0], col) : raise Exception("Expect snp list in map file to exactly match snp list in dat file")
        del datfields[0]
        del datfields[1]
        del datfields[2]
        assert len(row) == datfields.shape[1], "Expect # iids in fam file to match dat file"
        val = datfields.values.T
        snpdata = SnpData(iid=row,sid=col,pos=col_property,val=val)
        return snpdata

    def copyinputs(self, copier):
        # doesn't need to self.run_once() because creates name of all files itself
        copier.input(SnpReader._name_of_other_file(self.filename,remove_suffix="dat", add_suffix="dat"))
        copier.input(SnpReader._name_of_other_file(self.filename,remove_suffix="dat", add_suffix="fam"))
        copier.input(SnpReader._name_of_other_file(self.filename,remove_suffix="dat", add_suffix="map"))

    @staticmethod
    def write(filename, snpdata):
        """Writes a :class:`SnpData` to dat/fam/map format and returns the :class:`.Dat`.

        :param filename: the name of the file to create
        :type filename: string
        :param snpdata: The in-memory data that should be written to disk.
        :type snpdata: :class:`SnpData`
        :rtype: :class:`.Dat`

        >>> from pysnptools.snpreader import Dat, Bed
        >>> import pysnptools.util as pstutil
        >>> snpdata = Bed('../examples/toydata.bed',count_A1=False)[:,:10].read()  # Read first 10 snps from Bed format
        >>> pstutil.create_directory_if_necessary("tempdir/toydata10.dat")
        >>> Dat.write("tempdir/toydata10.dat",snpdata)              # Write data in dat/fam/map format
        Dat('tempdir/toydata10.dat')
        """

        if isinstance(filename,SnpData) and isinstance(snpdata,str): #For backwards compatibility, reverse inputs if necessary
            warnings.warn("write statement should have filename before data to write", DeprecationWarning)
            filename, snpdata = snpdata, filename 

        SnpReader._write_fam(snpdata, filename, remove_suffix="dat")
        SnpReader._write_map_or_bim(snpdata, filename, remove_suffix="dat", add_suffix="map")
        filename = SnpReader._name_of_other_file(filename,remove_suffix="dat", add_suffix="dat")

        snpsarray = snpdata.val
        with open(filename,"w") as dat_filepointer:
            for sid_index, sid in enumerate(snpdata.sid):
                if sid_index % 1000 == 0:
                    logging.info("Writing snp # {0} to file '{1}'".format(sid_index, filename))
                dat_filepointer.write("%s\tj\tn\t"%sid) #use "j" and "n" as the major and minor allele
                row = snpsarray[:,sid_index]
                dat_filepointer.write("\t".join((str(i) for i in row)) + "\n")
        logging.info("Done writing " + filename)
        return Dat(filename)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
