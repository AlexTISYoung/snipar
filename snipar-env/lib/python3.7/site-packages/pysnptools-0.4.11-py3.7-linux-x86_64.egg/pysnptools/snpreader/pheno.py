import numpy as np
import subprocess, sys, os.path
from itertools import *
import pandas as pd
import logging
from pysnptools.snpreader import SnpReader
from pysnptools.snpreader import SnpData
import pysnptools.util.pheno as pstpheno
import pysnptools.util as pstutil
from pysnptools.pstreader import _OneShot
from six.moves import range

class Pheno(_OneShot, SnpReader):
    '''
    A :class:`.SnpReader` for reading "alternative" phenotype files from disk.

    See :class:`.SnpReader` for general examples of using SnpReaders.

    This text format is described in http://zzz.bwh.harvard.edu/plink/data.shtml#pheno and looks like::

         FID    IID      qt1   bmi    site  
         F1     1110     2.3   22.22  2     
         F2     2202     34.12 18.23  1     
         ...

    where the heading row is optional.


    **Constructor:**
        :Parameters: * **input** (*string*) -- The phenotype file to read.
                     * **iid_if_none** (*None* or array of strings) -- An :attr:`SnpReader.iid` to use if the file is empty.
                     * **missing** (*None* or string) -- The value in the file that represents missing data.

        For backwards compatibility with older code, **input** can be a dictionary (instead of a filename) with these fields:

        * 'header' : [1] array phenotype name,
        * 'vals'   : [N*1] array of phenotype data,
        * 'iid'    : [N*2] array of family IDs and case IDs

        :Example:

        >>> from __future__ import print_function #Python 2 & 3 compatibility
        >>> from pysnptools.snpreader import Pheno, Bed
        >>> data_on_disk = Pheno('../examples/toydata.phe')
        >>> print((data_on_disk.iid_count, data_on_disk.sid_count))
        (500, 1)

    **Methods beyond** :class:`.SnpReader`

    '''

    def __init__(self, input, iid_if_none=None, missing=None):
        '''
        input    : string of the name of the file or an in-memory dictionary
        '''
        super(Pheno, self).__init__()

        self.filename = input
        self._iid_if_none = iid_if_none
        self.missing = missing

    def _read_pstdata(self):
        #LATER switch it, so the main code is here rather than in loadPhen
        if isinstance(self.filename,str):
            pheno_input = pstpheno.loadPhen(self.filename,missing=self.missing)
        elif self.filename is None:
            assert self._iid_if_none is not None, "If input is None then iid_if_none be given"
            pheno_input = {
            'header':np.empty((0),dtype='str'),
            'vals': np.empty((len(self._iid_if_none), 0)),
            'iid': self._iid_if_none
            }
        else:
            pheno_input = self.filename


        if len(pheno_input['vals'].shape) == 1:
            pheno_input = {
            'header' : pheno_input['header'],
            'vals' : np.reshape(pheno_input['vals'],(-1,1)),
            'iid' : pheno_input['iid']
            }

        if len(pheno_input['header']) > 0 and pheno_input['header'][0] is None:
            pheno_input['header'] = ["pheno{0}".format(i) for i in range(len(pheno_input['header']))] #LATER move to reader?
        elif len(pheno_input['header']) == 0:
            pheno_input['header'] = ["pheno{0}".format(i) for i in range(pheno_input['vals'].shape[1])]

        row = pheno_input['iid']
        col = np.array(pheno_input['header'],dtype='str')
        col_property = np.empty((len(col),3))
        col_property.fill(np.nan)
        val = pheno_input['vals']

        snpdata = SnpData(iid=row,sid=col,pos=col_property,val=val)
        return snpdata

    @staticmethod
    def write(filename, snpdata, missing='NaN', sep='\t'):
        """Writes a :class:`SnpData` to Pheno format and returns the :class:`.Pheno`.

        :param filename: the name of the file to create
        :type filename: string
        :param missing: value to threat as missing (default 'NaN')
        :type missing: string
        :param sep: the separator (default '\t')
        :type sep: string
        :rtype: :class:`.Pheno`


        >>> from pysnptools.snpreader import Pheno, Bed
        >>> import pysnptools.util as pstutil
        >>> snpdata = Bed('../examples/toydata.bed',count_A1=False)[:,:10].read()  # Read first 10 snps from Bed format
        >>> pstutil.create_directory_if_necessary("tempdir/toydata10.phe")
        >>> Pheno.write("tempdir/toydata10.txt",snpdata)       # Write data in Pheno format
        Pheno('tempdir/toydata10.txt')
        """
        with open(filename, 'w') as f:
            for i in range(snpdata.iid_count):
                tmpstr = snpdata.iid[i,0] + sep + snpdata.iid[i,1]
                for m in range(snpdata.sid_count):
                    v = snpdata.val[i,m]
                    if np.isnan(v):
                        vs = missing
                    else:
                        vs = str(v)
                    tmpstr += sep + vs
                tmpstr += "\n"
                f.write(tmpstr)
        return Pheno(filename,missing=missing)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()

