"""Tools for reading and manipulating SNP data.
"""

def _snps_fixup(snp_input, iid_if_none=None,count_A1=None):
    if isinstance(snp_input, str):
        return Bed(snp_input,count_A1=count_A1)

    if isinstance(snp_input, dict):
        return SnpData(iid=snp_input['iid'],sid=snp_input['header'],val=snp_input['vals'])

    if snp_input is None:
        assert iid_if_none is not None, "snp_input cannot be None here"
        return SnpData(iid_if_none, sid=np.empty((0),dtype='str'), val=np.empty((len(iid_if_none),0)),pos=np.empty((0,3)),name="") #todo: make a static factory method on SnpData

    return snp_input


from pysnptools.snpreader.snpreader import SnpReader
from pysnptools.snpreader.snpdata import SnpData
from pysnptools.snpreader.bed import Bed
from pysnptools.snpreader.ped import Ped
from pysnptools.snpreader.dat import Dat
from pysnptools.snpreader.snphdf5 import SnpHdf5
from pysnptools.snpreader.snphdf5 import Hdf5
from pysnptools.snpreader.snpnpz import SnpNpz
from pysnptools.snpreader.dense import Dense
from pysnptools.snpreader.pheno import Pheno
from pysnptools.snpreader.snpmemmap import SnpMemMap
from pysnptools.snpreader._mergesids import _MergeSIDs
from pysnptools.snpreader._mergeiids import _MergeIIDs
from pysnptools.snpreader.snpgen import SnpGen
from pysnptools.snpreader.distributedbed import DistributedBed, _Distributed1Bed