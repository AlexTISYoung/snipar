import numpy as np
import logging
from pysnptools.snpreader import SnpReader
from pysnptools.snpreader import SnpData
import warnings
from pysnptools.pstreader import PstData
from pysnptools.pstreader import _OneShot
from six.moves import range

def zero_family(s):
    '''Given an input id from the file, returns 0 as the family id and that input id as the case id.
    '''
    return "0",s

def zero_pos(s):
    '''Given an input id from the file, return the input id as the sid and returns (0,0,0) as the chromosome, genetic distance, basepair distance.
    '''
    return 0,s,0,0

def just_case(iid):
    '''Use just the case id (ignoring the family id) in the file.
    '''
    return iid[1]

def just_sid(sid,pos):
    '''Use just the sid (ignoring the chromosome, genetic distance, and basepair distance) in the file.
    '''
    return sid

class Dense(_OneShot,SnpReader):
    '''
    A :class:`.SnpReader` for reading \*.dense.txt (0,1,2 text files) from disk.

    See :class:`.SnpReader` for general examples of using SnpReaders.

    This format looks like::

        var	4006	9570    22899	37676	41236	41978	55159	66828...
        1-10004-rs12354060	22222222...
        1-707348-rs12184279	222222?2...
        1-724325-rs12564807	00000000...
        ...

    where rows are sids and columns are iids.

    **Constructor:**
        :Parameters: * **filename** (*string*) -- The Dense file to read.
                     * **extract_iid_function** (*string*) -- A function that breaks the row id from the file into a family id and an case id. Defaults to
                       setting the family id to 0 and using the whole row id as the iid.
                     * **extract_sid_pos_function** (*string*) -- A function that breaks the column id from the file into (chromosome, sid, genetic distance, basepair distance).
                       Defaults to setting the :attr:`.SnpReader.pos` information to 0 and using the whole column id as the sid.

        :Example:

        >>> from __future__ import print_function #Python 2 & 3 compatibility
        >>> from pysnptools.snpreader import Dense
        >>> data_on_disk = Dense('../examples/toydata100.dense.txt')
        >>> print((data_on_disk.iid_count, data_on_disk.sid_count))
        (500, 100)

    **Methods beyond** :class:`.SnpReader`

    '''
    def __init__(self, filename, extract_iid_function=zero_family, extract_sid_pos_function=zero_pos):
        '''
        filename    : string of the name of the Dat file.
        '''
        super(Dense, self).__init__()
        self.filename = filename
        self.extract_iid_function = extract_iid_function
        self.extract_sid_pos_function = extract_sid_pos_function

    def _read_pstdata(self):
        bim_list = []
        val_list_list = []
        with open(self.filename,"r") as fp:
            header = fp.readline()
            iid_string_list = header.strip().split()[1:]
            iid = np.array([self.extract_iid_function(iid_string) for iid_string in iid_string_list],dtype='str')
            val_list = []
            zerofloat = float('0'[0])
            missing_char = "?"[0]
            for line_index,line in enumerate(fp):
                if line_index % 1000 == 0:
                    logging.info("reading sid and iid info from line {0} of file '{1}'".format(line_index, self.filename))
                sid_string_rest = line.strip().split()
                sid_string = sid_string_rest[0]
                rest = [] if len(sid_string_rest)==1 else sid_string_rest[1]
                assert len(rest) == len(iid)
                bim_list.append(self.extract_sid_pos_function(sid_string))
                val_list = np.array([float(val)-zerofloat if val!=missing_char else np.NaN for val in rest])
                val_list_list.append(val_list)

        col = np.array([bim[1] for bim in bim_list],dtype='str')
        col_property = np.array([[bim[0],bim[2],bim[3]] for bim in bim_list],dtype=np.float64)

        val = np.zeros((len(iid),len(col)))
        for col_index in range(len(col)):
            val[:,col_index] = val_list_list[col_index]

        return PstData(iid,col,val,col_property=col_property,name=self.filename)

    @staticmethod
    def write(filename, snpdata, join_iid_function=just_case,join_sid_pos_function=just_sid):
        """Writes a :class:`SnpData` to Dense format. The values must be 0,1,2 (or missing). Returns a :class:`.Dense` with the default extractors.

        :param filename: the name of the file to create
        :type filename: string
        :param snpdata: The in-memory data that should be written to disk.
        :type snpdata: :class:`SnpData`
        :param join_iid_function: function to turn a family id and case id into a file id for columns.
                                  Defaults ignoring the family id and using the case id as the column id.
        :type join_iid_function: a function
        :param join_sid_pos_function: function to turn a sid and pos data into a file id for rows.
                                      Defaults ignoring the :attr:`.SnpReader.pos` information and using the sid id as the row id.
        :type join_sid_pos_function: a function
        :rtype: :class:`.Dense`

        >>> from pysnptools.snpreader import Dense, Bed
        >>> import pysnptools.util as pstutil
        >>> snpdata = Bed('../examples/toydata.bed',count_A1=False)[:,:10].read()  # Read first 10 snps from Bed format
        >>> pstutil.create_directory_if_necessary("tempdir/toydata10.dense.txt")
        >>> Dense.write("tempdir/toydata10.dense.txt",snpdata)        # Write data in Dense format
        Dense('tempdir/toydata10.dense.txt')
        """

        if isinstance(filename,SnpData) and isinstance(snpdata,'S'): #For backwards compatibility, reverse inputs if necessary
            warnings.warn("write statement should have filename before data to write", DeprecationWarning)
            filename, snpdata = snpdata, filename 

        snpsarray = snpdata.val
        with open(filename,"w") as filepointer:
            filepointer.write("var\t")
            filepointer.write("\t".join((join_iid_function(iid_pair) for iid_pair in snpdata.iid)) + "\n")

            for sid_index, sid in enumerate(snpdata.sid):
                pos = snpdata.pos[sid_index]
                if sid_index % 1000 == 0:
                    logging.info("Writing snp # {0} to file '{1}'".format(sid_index, filename))
                filepointer.write("%s\t" % join_sid_pos_function(sid,pos))
                row = snpsarray[:,sid_index]
                filepointer.write("".join((str(int(i)) if i==i else "?" for i in row)) + "\n")
        logging.info("Done writing " + filename)
        return Dense(filename)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
