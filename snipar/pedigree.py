from itertools import combinations

import numpy as np
import pandas as pd
import logging
from snipar.utilities import make_id_dict

def get_sibpairs_from_ped(ped, *, return_pairs=True):
    """Normalize family IDs by their observed parent pair, optionally enumerating siblings."""
    ped = np.asarray(ped, dtype=str)
    valid = (ped[:, 2] != '0') & (ped[:, 3] != '0') & ~np.char.startswith(ped[:, 0], '_')
    ped = ped[valid]
    if not len(ped):
        return None, ped
    # Stable ties retain pedigree order within each sorted parent pair.
    order = np.lexsort((ped[:, 3], ped[:, 2]))
    pairs_of_parents = ped[order, 2:4]
    new_family = np.empty(len(ped), dtype=bool)
    new_family[0] = True
    new_family[1:] = np.any(pairs_of_parents[1:] != pairs_of_parents[:-1], axis=1)
    labels = np.cumsum(new_family) - 1
    width = len(str(labels[-1]))
    if width > ped.dtype.itemsize // 4:
        ped = ped.astype(f'U{width}')
    ped[order, 0] = labels.astype(f'U{width}')
    if not return_pairs:
        return None, ped
    starts = np.flatnonzero(new_family)
    stops = np.append(starts[1:], len(ped))
    pairs = []
    for start, stop in zip(starts, stops):
        if stop - start > 1:
            pairs.extend(combinations(ped[order[start:stop], 1], 2))
    return (np.asarray(pairs, dtype=ped.dtype) if pairs else None), ped

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


def create_pedigree(king_address, agesex_address, same_parents_in_ped=True):
    """Creates pedigree table from agesex file and kinship file in KING format.
    
    Args:
        king_address : str
            Address of a kinship file in KING format. kinship file is a '\t' separated csv with columns "FID1", "ID1", "FID2", "ID2, "InfType".
            Each row represents a relationship between two individuals. InfType column states the relationship between two individuals.
            The only relationships that matter for this script are full sibling and parent-offspring which are shown by 'FS' and 'PO' respectively.
            This file is used in creating a pedigree file and can be generated using KING.
            As fids starting with '_' are reserved for control there should be no fids starting with '_'.

        agesex_address : str
            Address of the agesex file. This is a " " separated CSV with columns "FID", "IID", "FATHER_ID", "MOTHER_ID", "sex", "age".
            Each row contains the age and sex of one individual. Male and Female sex should be represented with 'M' and 'F'.
            Age column is used for distinguishing between parent and child in a parent-offspring relationship inferred from the kinship file.
            ID1 is a parent of ID2 if there is a 'PO' relationship between them and 'ID1' is at least 12 years older than ID2.
    
    Returns:
        pd.DataFrame:
            A pedigree table with 'FID', 'IID', 'FATHER_ID', 'MOTHER_ID'. Each row represents an individual.
    """

    kinship = pd.read_csv(king_address, delimiter="\t", dtype=str, keep_default_na=False)
    logging.info("loaded kinship file")
    # Filter MZ
    mz_kin = kinship.loc[kinship['InfType']=='Dup/MZ']
    if mz_kin.shape[0]>0:
        print('Found '+str(mz_kin.shape[0])+' duplicates/MZ twin pairs. Removing one from each pair')
        mz_id1 = set(np.unique(mz_kin.loc[:,['ID1']]))
        id1_in_mz = np.array([x in mz_id1 for x in kinship.ID1])
        id2_in_mz = np.array([x in mz_id1 for x in kinship.ID2])
        in_mz = np.logical_or(id1_in_mz,id2_in_mz)
        kinship = kinship.drop(kinship[in_mz].index)
    agesex = pd.read_csv(agesex_address, sep=r"\s+", dtype={"IID": str, "FID": str},
                         keep_default_na=False)
    agesex["age"] = pd.to_numeric(agesex["age"], errors="coerce")
    logging.info("loaded agesex file")
    agesex = agesex.set_index("IID")
    logging.info("creating age and sex dictionaries")
    kinship = pd.merge(kinship, agesex.rename(columns={"sex":"sex1", "age":"age1"}), left_on="ID1", right_index=True)
    kinship = pd.merge(kinship, agesex.rename(columns={"sex":"sex2", "age":"age2"}), left_on="ID2", right_index=True)
    logging.info("dictionaries created")
    people = {}
    fid_counter = 0
    merged_fids = {}
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
            if age1 >= age2 + 12:
                if sex1 == "F":
                    p2.mid = p1.id
                if sex1 == "M":
                    p2.pid = p1.id

            if age2 >= age1 + 12:
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

            a, b = p1.fid, p2.fid
            while a in merged_fids:
                a = merged_fids[a]
            while b in merged_fids:
                b = merged_fids[b]
            if a != b:
                a, b = sorted((a, b), key=int)
                merged_fids[b] = a

    data = []
    for p in people.values():
        if p.fid is None:
            p.fid = str(fid_counter)
            fid_counter += 1

        while p.fid in merged_fids:
            p.fid = merged_fids[p.fid]

        if p.mid is None:
            #default mother id
            p.mid = p.fid + "___M"

        if p.pid is None:
            #default father ir
            p.pid = p.fid + "___P"
        
        
        data.append((p.fid, p.id, p.pid, p.mid))

    data = pd.DataFrame(data, columns = ['FID' , 'IID', 'FATHER_ID' , 'MOTHER_ID']).astype(str)
    
    # Check if parents are the same in ped for each sibship
    if same_parents_in_ped:
        parent_pairs = data[["FID", "FATHER_ID", "MOTHER_ID"]].drop_duplicates()
        inconsistent_fams = parent_pairs.loc[parent_pairs.FID.duplicated(), "FID"].unique()
        if len(inconsistent_fams):
            rows = data.FID.isin(inconsistent_fams)
            data.loc[rows, "FATHER_ID"] = data.loc[rows, "FID"] + "___P"
            data.loc[rows, "MOTHER_ID"] = data.loc[rows, "FID"] + "___M"
            print(f"Warning: found {len(inconsistent_fams)} sibships with different parents "
                  "when creating pedigree. Removing the parents from pedigree.")
    
    return data

def find_individuals_with_sibs(ids, ped, gts_ids, return_ids_only = False):
    """
    Used in get_gts_matrix and get_fam_means to find the individuals in ids that have genotyped siblings.
    """
    ped_dict = make_id_dict(ped, 1)
    gts_fams = np.full(len(gts_ids), "", dtype=ped.dtype)
    families = {}
    for i, iid in enumerate(gts_ids):
        if iid in ped_dict:
            fam, _, father, mother = ped[ped_dict[iid], :4]
            gts_fams[i] = fam
            families.setdefault(fam, []).append((father, mother))
    sibships = {fam for fam, parents in families.items()
                if len(parents) > 1 and all(p == parents[0] for p in parents[1:])}

    indices = np.array([i for i, iid in enumerate(ids)
                        if iid in ped_dict and ped[ped_dict[iid], 0] in sibships], dtype=int)
    selected = ids[indices]
    if return_ids_only:
        return selected
    ids_fams = np.asarray([ped[ped_dict[iid], 0] for iid in selected], dtype=ped.dtype)
    return selected, ids_fams, gts_fams

def get_sibpairs_from_king(kinfile):
    kinship = pd.read_csv(kinfile, sep="\t", usecols=["ID1", "ID2", "InfType"],
                          dtype=str, keep_default_na=False)
    return kinship.loc[kinship["InfType"] == "FS", ["ID1", "ID2"]].to_numpy(dtype=str)
