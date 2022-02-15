import numpy as np
import pandas as pd
import logging
from snipar.utilities import make_id_dict

def get_sibpairs_from_ped(ped):
    # Remove rows with missing parents
    parent_missing = np.array([ped[i,2]=='0' or ped[i,3]=='0' for i in range(ped.shape[0])])
    #print('Removing '+str(np.sum(parent_missing))+' rows from pedigree due to missing parent(s)')
    ped = ped[np.logical_not(parent_missing),:]
    # Remove control families
    controls = np.array([x[0]=='_' for x in ped[:,0]])
    ped = ped[~controls,:]
    # Find unique parent-pairs
    parent_pairs = np.array([ped[i,2]+ped[i,3] for i in range(ped.shape[0])])
    unique_pairs, sib_counts = np.unique(parent_pairs, return_counts=True)
    # Fill in sibships
    for i in range(unique_pairs.shape[0]):
        ped[np.where(parent_pairs==unique_pairs[i])[0],0] = i 
    # Parent pairs with more than one offspring
    sib_parent_pairs = unique_pairs[sib_counts>1]
    fam_sizes = np.array([x*(x-1)/2 for x in sib_counts[sib_counts>1]])
    npairs = int(np.sum(fam_sizes))
    if npairs==0:
        return None, ped
    else:    
        # Find all sibling pairs
        sibpairs = np.zeros((npairs,2),dtype=ped.dtype)
        paircount = 0
        for ppair in sib_parent_pairs:
            sib_indices = np.where(parent_pairs==ppair)[0]
            for i in range(0,sib_indices.shape[0]-1):
                for j in range(i+1,sib_indices.shape[0]):
                    sibpairs[paircount,:] = np.array([ped[sib_indices[i],1],ped[sib_indices[j],1]])
                    paircount += 1
        return sibpairs, ped

class Person:
    """Just a simple data structure representing individuals

    Args:
        id : str
            IID of the individual.
        fid : str
            FID of the individual.
        pid : str
            IID of the father of that individual.
        mid : str
            IID of the mother of that individual.
    """
    def __init__(self, id, fid=None, pid=None, mid=None):
        self.id = id
        self.fid = fid
        self.pid = pid
        self.mid = mid


def create_pedigree(king_address, agesex_address):
    """Creates pedigree table from agesex file and kinship file in KING format.

    Args:
        king_address : str
            Address of a kinship file in KING format. kinship file is a '\t' seperated csv with columns "FID1", "ID1", "FID2", "ID2, "InfType".
            Each row represents a relationship between two individuals. InfType column states the relationship between two individuals.
            The only relationships that matter for this script are full sibling and parent-offspring which are shown by 'FS' and 'PO' respectively.
            This file is used in creating a pedigree file and can be generated using KING.
            As fids starting with '_' are reserved for control there should be no fids starting with '_'.

        agesex_address : str
            Address of the agesex file. This is a " " seperated CSV with columns "FID", "IID", "FATHER_ID", "MOTHER_ID", "sex", "age".
            Each row contains the age and sex of one individual. Male and Female sex should be represented with 'M' and 'F'.
            Age column is used for distinguishing between parent and child in a parent-offspring relationship inferred from the kinship file.
            ID1 is a parent of ID2 if there is a 'PO' relationship between them and 'ID1' is at least 12 years older than ID2.

    Returns:
        pd.DataFrame:
            A pedigree table with 'FID', 'IID', 'FATHER_ID', 'MOTHER_ID'. Each row represents an individual.
    """

    kinship = pd.read_csv(king_address, delimiter="\t").astype(str)
    logging.info("loaded kinship file")
    agesex = pd.read_csv(agesex_address, delim_whitespace=True)
    agesex["IID"] = agesex["IID"].astype(str)
    agesex["FID"] = agesex["FID"].astype(str)
    logging.info("loaded agesex file")
    agesex = agesex.set_index("IID")
    logging.info("creating age and sex dictionaries")
    kinship = pd.merge(kinship, agesex.rename(columns={"sex": "sex1", "age": "age1"}), left_on="ID1", right_index=True)
    kinship = pd.merge(kinship, agesex.rename(columns={"sex": "sex2", "age": "age2"}), left_on="ID2", right_index=True)
    logging.info("dictionaries created")
    people = {}
    fid_counter = 0
    dropouts = []
    kinship_cols = kinship.columns.tolist()
    index_id1 = kinship_cols.index("ID1")
    index_id2 = kinship_cols.index("ID2")
    index_sex1 = kinship_cols.index("sex1")
    index_sex2 = kinship_cols.index("sex2")
    index_age1 = kinship_cols.index("age1")
    index_age2 = kinship_cols.index("age2")
    index_inftype = kinship_cols.index("InfType")
    logging.info("creating pedigree objects")
    pop_size = kinship.values.shape[0]
    t = kinship.values.tolist()
    for row in range(pop_size):
        relation = t[row][index_inftype]
        id1 = t[row][index_id1]
        id2 = t[row][index_id2]
        age1 = t[row][index_age1]
        age2 = t[row][index_age2]
        sex1 = t[row][index_sex1]
        sex2 = t[row][index_sex2]
        p1 = people.get(id1)
        if p1 is None:
            p1 = Person(id1)
            people[id1] = p1

        p2 = people.get(id2)
        if p2 is None:
            p2 = Person(id2)
            people[id2] = p2

        if relation == "PO":
            if age1 > age2 + 12:
                if sex1 == "F":
                    p2.mid = p1.id
                if sex1 == "M":
                    p2.pid = p1.id

            if age2 > age1 + 12:
                if sex2 == "F":
                    p1.mid = p2.id
                if sex2 == "M":
                    p1.pid = p2.id
        if relation == "FS":
            if p1.fid is None and p2.fid is None:
                p1.fid = str(fid_counter)
                p2.fid = str(fid_counter)
                fid_counter += 1

            if p1.fid is None and p2.fid is not None:
                p1.fid = p2.fid

            if p1.fid is not None and p2.fid is None:
                p2.fid = p1.fid

    for excess in dropouts:
        people.pop(excess)

    data = []
    for p in people.values():
        if p.fid is None:
            p.fid = str(fid_counter)
            fid_counter += 1

        if p.mid is None:
            # default mother id
            p.mid = p.fid + "___M"

        if p.pid is None:
            # default father ir
            p.pid = p.fid + "___P"

        data.append((p.fid, p.id, p.pid, p.mid))

    data = pd.DataFrame(data, columns=['FID', 'IID', 'FATHER_ID', 'MOTHER_ID']).astype(str)
    return data

def find_individuals_with_sibs(ids, ped, gts_ids, return_ids_only = False):
    """
    Used in get_gts_matrix and get_fam_means to find the individuals in ids that have genotyped siblings.
    """
    # Find genotyped sibships of size > 1
    ped_dict = make_id_dict(ped, 1)
    ids_in_ped = np.array([x in ped_dict for x in gts_ids])
    gts_fams = np.zeros((gts_ids.shape[0]),dtype=gts_ids.dtype)
    gts_fams[ids_in_ped] = np.array([ped[ped_dict[x], 0] for x in gts_ids[ids_in_ped]])
    fams, counts = np.unique(gts_fams[ids_in_ped], return_counts=True)
    sibships = set(fams[counts > 1])
    # Find individuals with genotyped siblings
    ids_in_ped = np.array([x in ped_dict for x in ids])
    ids = ids[ids_in_ped]
    ids_fams = np.array([ped[ped_dict[x], 0] for x in ids])
    ids_with_sibs = np.array([x in sibships for x in ids_fams])
    ids = ids[ids_with_sibs]
    ids_fams = ids_fams[ids_with_sibs]
    if return_ids_only:
        return ids
    else:
        return ids, ids_fams, gts_fams

def get_sibpairs_from_king(kinfile):
    kin_header = np.array(open(kinfile,'r').readline().split('\t'))
    inf_type_index = np.where(np.array([x[0:7]=='InfType' for x in kin_header]))[0][0]
    id1_index = np.where(np.array(kin_header)=='ID1')[0][0]
    id2_index = np.where(np.array(kin_header)=='ID2')[0][0]
    sibpairs = np.loadtxt(kinfile,dtype=str,skiprows=1,usecols=(id1_index,id2_index,inf_type_index))
    sibpairs = sibpairs[sibpairs[:,2]=='FS',0:2]
    return sibpairs
