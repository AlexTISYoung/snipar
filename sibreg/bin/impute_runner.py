import logging
from impute_from_sibs import *
import time
import argparse
import h5py
import random
import pandas as pd
import os

random.seed(1567924)
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
    fid_counter = 1
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
            fid_counter += 1
            #not neccessary
            for p in families[p1]:
                p.fid = fid_counter
            for p in families[p2]:
                p.fid = fid_counter

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



#does the imputation and writes the results
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s')
parser = argparse.ArgumentParser()
parser.add_argument('-c', action='store_true')
parser.add_argument('from_chr',type=int,help='Which chromosome (<=)')
parser.add_argument('to_chr',type=int,help='Which chromosome (<)')
parser.add_argument('ibd',type=str,help='IBD file')
parser.add_argument('genotypes_prefix',type=str,help='prefix of genotypes in .bed format')
parser.add_argument('--bim',type=str,default = None, help='Bim file giving positions of SNPs in KING IBD file if different from Bim file of genotypes')
parser.add_argument('--out_prefix',type=str,default = None, help="Writes the result of imputation for chromosome i to outprefix{i}")
parser.add_argument('--start', type=int,
                    help='Start index of SNPs to perform imputation for in genotype file (starting at zero)',
                    default=None)
parser.add_argument('--end', type=int,
                    help='End index of SNPs to perform imputation for in genotype file (goes from 0 to (end-1)',
                    default=None)
parser.add_argument('--pedigree',type=str,default = None, help='Pedigree file with siblings sharing family ID')
parser.add_argument('--king',type=str,default = None, help='Address of the king file')
parser.add_argument('--agesex',type=str,default = None, help='Address of the agesex file with header "FID IID age sex"')

args=parser.parse_args()
pedigree_address = args.pedigree
#fids starting with _ are reserved for control
#Families should not have grandparents
if not args.pedigree:
    logging.info("creating pedigree ...")
    pedigree_address = "test_data/__tmp_pedigree"
    pedigree = create_pedigree(args.king, args.agesex)
    pedigree.to_csv(pedigree_address, sep = " ", index = False)
if args.c:
    logging.info("Adding control to the pedigree ...")
    pedigree = pd.read_csv(pedigree_address, sep = " ").astype(str)
    pedigree = add_control(pedigree)
    dirname, filename = os.path.split(pedigree_address)
    pedigree_address = os.path.join(dirname, "__controlled"+filename)
    pedigree.to_csv(pedigree_address, sep = " ", index = False)

consumed_time = 0
for chromosome in range(args.from_chr, args.to_chr):
    print(chromosome, " is chromosome")
    sibships, iid_to_bed_index, gts, ibd, pos, sid = prepare_data(pedigree_address, args.genotypes_prefix+str(chromosome), args.ibd, chromosome, args.start, args.end, args.bim)
    gts = gts.astype(float)
    pos = pos.astype(int)
    start_time = time.time()
    if args.out_prefix is None:
        imputed_fids, imputed_par_gts = impute(sibships, iid_to_bed_index, gts, ibd, pos, sid, "test_data/parent_imputed_chr"+str(chromosome))
    else:
        imputed_fids, imputed_par_gts = impute(sibships, iid_to_bed_index, gts, ibd, pos, sid, args.out_prefix+str(chromosome))
    end_time = time.time()
    consumed_time += (end_time-start_time)

logging.info("imputation time: "+str(consumed_time))
