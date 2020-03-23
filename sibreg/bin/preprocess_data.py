import logging
import pandas as pd
import numpy as np
from pysnptools.snpreader import Bed
class Person:
    def __init__(self, id, fid=None, pid=None, mid=None):
        self.id = id
        self.fid = fid
        self.pid = pid
        self.mid = mid

def recurcive_append(dictionary, index, element):
    queue = {index}
    seen_so_far = set()
    while queue:
        current_index = queue.pop()
        seen_so_far.add(current_index)
        dictionary[current_index].add(element)
        queue = queue.union(dictionary[current_index])
        queue = queue.difference(seen_so_far)

def create_pedigree(king_address, agesex_address):
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
    data = pd.DataFrame(data, columns = ['FID' , 'IID', 'FATHER_ID' , 'MOTHER_ID'])
    return data
    
def add_control(pedigree):
    #For each family that has two or more siblings and both parents,
    # creates a new family with all the sibling and no parents.
    #fid of this family is "_"+original_fid
    #IIDs are the same in both families.
    pedigree["has_mother"] = pedigree["MOTHER_ID"].isin(pedigree["IID"])
    pedigree["has_father"] = pedigree["FATHER_ID"].isin(pedigree["IID"])
    has_parents = pedigree[pedigree["has_father"] & pedigree["has_mother"]]
    count_of_sibs_in_fam = has_parents.groupby(["FID", "FATHER_ID", "MOTHER_ID"]).count().reset_index()
    FIDs_with_multiple_sibs = count_of_sibs_in_fam[count_of_sibs_in_fam["IID"] > 1][["FID"]]
    families_with_multiple_sibs = pedigree.merge(FIDs_with_multiple_sibs, on = "FID")
    families_with_multiple_sibs = families_with_multiple_sibs[(families_with_multiple_sibs["has_father"]) & (families_with_multiple_sibs["has_father"])]
    families_with_multiple_sibs["FID"] = "_" + families_with_multiple_sibs["FID"].astype(str)
    families_with_multiple_sibs["MOTHER_ID"] = "_" + families_with_multiple_sibs["FID"].astype(str) + "_M"
    families_with_multiple_sibs["FATHER_ID"] = "_" + families_with_multiple_sibs["FID"].astype(str) + "_P"
    pedigree = pedigree.append(families_with_multiple_sibs)
    pedigree = pedigree[['FID' , 'IID', 'FATHER_ID' , 'MOTHER_ID']]
    return pedigree


def prepare_data(pedigree, genotypes_address, ibd, chr, start=None, end=None, bim_address = None):
    logging.info("initializing data")
    logging.info("loading and filtering pedigree file ...")
    #keeping individuals with no parents
    pedigree = pedigree.astype(str)
    pedigree["has_father"] = pedigree["FATHER_ID"].isin(pedigree["IID"])
    pedigree["has_mother"] = pedigree["MOTHER_ID"].isin(pedigree["IID"])
    no_parent_pedigree = pedigree[~(pedigree["has_mother"] | pedigree["has_father"])]
    #TODO handle sibs with one parent
    ped_ids =  set(no_parent_pedigree["IID"].tolist())
    #finding siblings in each family
    sibships = no_parent_pedigree.groupby(["FID", "FATHER_ID", "MOTHER_ID"]).agg({'IID':lambda x: list(x)}).reset_index()
    sibships["sib_count"] = sibships["IID"].apply(len)
    sibships = sibships[sibships["sib_count"]>1]
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
    hdf5_output_dict = {"sid":sid, "pedigree":pedigree_output}
    return sibships, iid_to_bed_index, gts, ibd, pos, hdf5_output_dict
