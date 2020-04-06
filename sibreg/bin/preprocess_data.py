"""Contains functions for preprocessing data

Classes
-------
    Person

Functions
----------
    recurcive_append
    create_pedigree
    add_control
    prepare_data
"""
import logging
import pandas as pd
import numpy as np
from pysnptools.snpreader import Bed
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

def recurcive_append(dictionary, index, element):
    """Adds an element to value of all the keys that can be reached from index with using get recursively. 

    Args:
        dictionary : dict
            A dictionary of objects to list

        index
            The start point

        element
            What should be added to values
    """
    queue = {index}
    seen_so_far = set()
    while queue:
        current_index = queue.pop()
        seen_so_far.add(current_index)
        dictionary[current_index].add(element)
        queue = queue.union(dictionary[current_index])
        queue = queue.difference(seen_so_far)

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
    agesex = pd.read_csv(agesex_address, sep = " ")
    agesex["IID"] = agesex["IID"].astype(str)
    agesex["FID"] = agesex["FID"].astype(str)
    sex = {row["IID"]:row["sex"] for index, row in agesex.iterrows()}
    age = {row["IID"]:row["age"] for index, row in agesex.iterrows()}
    #Entails all individuals, key is id and value is the person with that id
    people = {}
    #lists relatives of a person, key is id and value is relatives of the person with that id
    families = {}
    #lists sibs of a person, key is id and value is sibs of the person with that id
    sibs = {}
    #MZs that should be removed
    dropouts = []
    fid_counter = 0
    #adds and element to dict[index], dict[dict[index][0]] ...
    for index, row in kinship.iterrows():
        p1 = people.get(row["ID1"])
        if p1 is None:
            p1 = Person(row["ID1"])
            people[row["ID1"]] = p1
            families[p1] = set()
            sibs[p1] = set()
        p2 = people.get(row["ID2"])
        if p2 is None:
            p2 = Person(row["ID2"])
            people[row["ID2"]] = p2
            families[p2] = set()
            sibs[p2] = set()
        recurcive_append(families, p1, p2)
        recurcive_append(families, p2, p1)
        #As the kinship only has first degree relatives, each row represent a family relation
        #Therefore they should share the same family id
        if p1.fid is not None and p2.fid is not None:
            p1.fid = p2.fid
            for p in families[p1]:
                p.fid = p2.fid

        if p1.fid is not None and p2.fid is None:
            p2.fid = p1.fid
            for p in families[p2]:
                p.fid = p1.fid

        if p2.fid is not None and p1.fid is None:
            p1.fid = p2.fid
            for p in families[p1]:
                p.fid = p2.fid
        
        if p2.fid is None and p1.fid is None:
            p1.fid = fid_counter
            p2.fid = fid_counter            
            #not neccessary
            for p in families[p1]:
                p.fid = fid_counter
            for p in families[p2]:
                p.fid = fid_counter
            fid_counter += 1

        #sets default parent id
        if p1.mid is None:
            p1.mid = p1.id + "_M"
        if p2.mid is None:
            p2.mid = p2.id + "_M"
        if p1.pid is None:
            p1.pid = p1.id + "_P"
        if p2.pid is None:
            p2.pid = p2.id + "_P"

        if row["InfType"] == "PO":
            #After checking sex and age sets the parent
            if age[p1.id] >  age[p2.id]+12:
                if sex[p1.id] == "F":
                    p2.mid = p1.id
                if sex[p1.id] == "M":
                    p2.pid = p1.id
                for sib in sibs[p2]:
                    sib.mid = p2.mid
                    sib.pid = p2.pid
            if age[p2.id] >  age[p1.id]+12:
                if sex[p2.id] == "F":
                    p1.mid = p2.id
                if sex[p2.id] == "M":
                    p1.pid = p2.id
                for sib in sibs[p1]:
                    sib.mid = p1.mid
                    sib.pid = p1.pid
        
        if row["InfType"] == "FS":
            #since these are full siblings the should share the same parents
            recurcive_append(sibs, p1, p2)
            recurcive_append(sibs, p2, p1)
            p1.mid = p2.mid
            p1.pid = p2.pid
            for p in sibs[p1]:
                p.mid = p2.mid
                p.pid = p2.pid

        if row["InfType"] == "MZ":
            #ranodmly removes one of MZs
            excess = row[random.choice(["ID1", "ID2"])]
            dropouts.append(excess)

    #removes excess from people
    for excess in dropouts:
        people.pop(excess)
    data = [(p.fid, p.id, p.pid, p.mid) for p in people.values()]
    data = pd.DataFrame(data, columns = ['FID' , 'IID', 'FATHER_ID' , 'MOTHER_ID']).astype(str)
    return data
    
def add_control(pedigree):
    """Adds control families to the pedigree table for testing.

    For each family that has two or more siblings and both parents, creates a 3 new familes, one has no parents, one with no mother and one with no father.
    gFID of these families are x+original_fid where x is "_o_", "_p_", "_m_" for these cases: no parent, only has father, only has mother. IIDs are the same in both families.

    Args:
        pedigree : pd.DataFrame
            A pedigree table with 'FID', 'IID', 'FATHER_ID', 'MOTHER_ID'. Each row represents an individual.
            fids starting with "_" are reserved for control.
        
    Returns:
        pd.DataFrame
            A pedigree table with 'FID', 'IID', 'FATHER_ID', 'MOTHER_ID'. Each row represents an individual.
            For each family with both parents and more than one offspring, it has a control family(fids for control families start with '_')

    """

    pedigree["has_mother"] = pedigree["MOTHER_ID"].isin(pedigree["IID"])
    pedigree["has_father"] = pedigree["FATHER_ID"].isin(pedigree["IID"])
    families_with_both_parents = pedigree[pedigree["has_father"] & pedigree["has_mother"]]
    count_of_sibs_in_fam = families_with_both_parents.groupby(["FID", "FATHER_ID", "MOTHER_ID"]).count().reset_index()
    FIDs_with_multiple_sibs = count_of_sibs_in_fam[count_of_sibs_in_fam["IID"] > 1][["FID"]]
    families_with_multiple_sibs = families_with_both_parents.merge(FIDs_with_multiple_sibs, on = "FID")

    families_with_multiple_sibs["FID"] = "_o_" + families_with_multiple_sibs["FID"].astype(str)
    families_with_multiple_sibs["MOTHER_ID"] = families_with_multiple_sibs["FID"].astype(str) + "_M"
    families_with_multiple_sibs["FATHER_ID"] = families_with_multiple_sibs["FID"].astype(str) + "_P"

    keep_mother = families_with_both_parents.copy()
    keep_mother["FID"] = "_m_" + keep_mother["FID"].astype(str)
    keep_mother["FATHER_ID"] = keep_mother["FID"].astype(str) + "_P"

    keep_father = families_with_both_parents.copy()
    keep_father["FID"] = "_p_" + keep_father["FID"].astype(str)
    keep_father["MOTHER_ID"] = keep_father["FID"].astype(str) + "_M"
    pedigree = pedigree.append(families_with_multiple_sibs).append(keep_father).append(keep_mother)
    pedigree = pedigree[['FID' , 'IID', 'FATHER_ID' , 'MOTHER_ID']]
    return pedigree


def prepare_data(pedigree, genotypes_address, ibd, chr, start=None, end=None, bim_address = None):
    """Processes the required data for the imputation and returns it.
    
    Args:
        pedigree : pd.DataFrame 
            The pedigree table. It contains 'FID', 'IID', 'FATHER_ID' and, 'MOTHER_ID' columns.
        
        genotypes_address : str
            Address of the bed file (does not inlude '.bed').
        
        ibd : pd.DataFrame
            A pandas dataframe containing IBD statuses for all SNPs.
            This It has these columns: "chr", "ID1", "ID2", "IBDType", "StartSNP", "StopSNP".
            Each line states an IBD segment between a pair on individuals. This can be generated using King software.

        chr : int
            Number of the chromosome that's going to imputed.

        start : int, optional
            This function can be used for preparing a slice of a chromosome. This is the location of the start of the slice.

        end : int, optional
            This function can be used for preparing a slice of a chromosome. This is the location of the end of the slice.

        bim_address : str, optional
            Address of the bim file if it's different from the address of the bed file. Does not include '.bim'.

    Returns:
        tuple(pandas.Dataframe, dict, numpy.ndarray, pandas.Dataframe, numpy.ndarray, numpy.ndarray)
            Returns the data required for the imputation. This data is a tuple of multiple objects.
                sibships: A pandas DataFrame with columns ['FID', 'FATHER_ID', 'MOTHER_ID', 'IID', 'has_father', 'has_mother', 'single_parent'] where IID columns is a list of the IIDs of individuals in that family.
                    It only contains families that have more than one child or only one parent.
                iid_to_bed_index: A str->int dictionary mapping IIDs of people to their location in bed file.
                gts: Numpy array containing the genotype data from the bed file.
                ibd: A pandas DataFrame with columns "ID1", "ID2", 'segment'. The segments column is a list of IBD segments between ID1 and ID2.
                    Each segment consists of a start, an end, and an IBD status. The segment list is flattened meaning it's like [start0, end0, ibd_status0, start1, end1, ibd_status1, ...]
                pos: A numpy array with the position of each SNP in the order of appearance in gts.
                hdf5_output_dict: A  dictionary whose values will be written in the imputation output under its keys.
    """
    logging.info("initializing data")
    logging.info("loading and filtering pedigree file ...")
    #keeping individuals with no parents
    pedigree = pedigree.astype(str)
    pedigree["has_father"] = pedigree["FATHER_ID"].isin(pedigree["IID"])
    pedigree["has_mother"] = pedigree["MOTHER_ID"].isin(pedigree["IID"])
    no_parent_pedigree = pedigree[~(pedigree["has_mother"] & pedigree["has_father"])]
    #TODO handle sibs with one parent
    ped_ids =  set(no_parent_pedigree["IID"].tolist())
    #finding siblings in each family
    sibships = no_parent_pedigree.groupby(["FID", "FATHER_ID", "MOTHER_ID", "has_father", "has_mother"]).agg({'IID':lambda x: list(x)}).reset_index()

    sibships["sib_count"] = sibships["IID"].apply(len)
    sibships["single_parent"] = sibships["has_father"] ^ sibships["has_mother"]
    sibships = sibships[(sibships["sib_count"]>1) | sibships["single_parent"]]
    fids = set([i for i in sibships["FID"].values.tolist() if i.startswith("_")])
    logging.info("loading bim file ...")
    if bim_address is None:
        bim_file = genotypes_address+'.bim'
    else:
        bim_file = bim_address
    bim = pd.read_csv(bim_file, sep = "\t", header=None, names=["Chr", "id", "morgans", "coordinate", "allele1", "allele2"])
    #TODO what if people are related but do not have ibd on chrom
    logging.info("loading and transforming ibd file ...")
    ibd = ibd.astype(str)
    #Adding location of start and end of each 
    ibd = ibd[ibd["Chr"] == str(chr)][["ID1", "ID2", "IBDType", "StartSNP", "StopSNP"]]
    ibd["IBDType"] = ibd["IBDType"].apply(lambda x: 2 if x=="IBD2" else 1)
    temp = bim[["id", "coordinate"]].rename(columns = {"id":"StartSNP","coordinate":"StartSNPLoc"})
    ibd= ibd.merge(temp, on="StartSNP")
    temp = bim[["id", "coordinate"]].rename(columns = {"id":"StopSNP","coordinate":"StopSNPLoc"})
    ibd = ibd.merge(temp, on="StopSNP")
    ibd['segment'] = ibd[['StartSNPLoc', 'StopSNPLoc', "IBDType"]].values.tolist()
    def create_seg_list(x):
        elements = list(x)
        result = []
        for el in elements:
            result = result+el
        return result
    ibd = ibd.groupby(["ID1", "ID2"]).agg({'segment':lambda x:create_seg_list(x)}).to_dict()["segment"]
    logging.info("loading bed file ...")
    gts_f = Bed(genotypes_address+".bed")
    ids_in_ped = [(id in ped_ids) for id in gts_f.iid[:,1]]
    gts_ids = gts_f.iid[ids_in_ped]
    if end is not None:        
        gts = gts_f[ids_in_ped , start:end].read().val.astype(float)
        pos = gts_f.pos[start:end, 2]
        sid = gts_f.sid[start:end]
    else:
        gts = gts_f[ [(id in ped_ids) for id in gts_f.iid[:,1]], :].read().val.astype(float)
        pos = gts_f.pos[:, 2]
        sid = gts_f.sid
    iid_to_bed_index = {i:index for index, i in enumerate(gts_ids[:,1])}
    logging.info("initializing data done ...")
    pedigree_output = np.concatenate(([pedigree.columns.values.tolist()], pedigree.values))
    pedigree_output = pedigree_output.astype('S')
    pedigree_output = np.array(pedigree_output,dtype='S')
    hdf5_output_dict = {"sid":sid, "pedigree":pedigree_output}
    return sibships, iid_to_bed_index, gts, ibd, pos, hdf5_output_dict
