"""Contains functions for preprocessing the data for the imputation

Classes
-------
    Person

Functions
----------
    recurcive_append
    create_pedigree
    add_control
    preprocess_king
    prepare_data
    compute_aics
    estimate_f
    prepare_gts
"""
import logging
import pandas as pd
import numpy as np
from pysnptools.snpreader import Bed
from bgen_reader import open_bgen, read_bgen
from snipar.config import nan_integer
from tqdm import tqdm
import statsmodels.api as sm
from scipy.stats import norm
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
    logging.info("loaded kinship file")
    agesex = pd.read_csv(agesex_address, delim_whitespace=True)
    agesex["IID"] = agesex["IID"].astype(str)
    agesex["FID"] = agesex["FID"].astype(str)
    logging.info("loaded agesex file")
    agesex = agesex.set_index("IID")
    logging.info("creating age and sex dictionaries")
    kinship = pd.merge(kinship, agesex.rename(columns={"sex":"sex1", "age":"age1"}), left_on="ID1", right_index=True)
    kinship = pd.merge(kinship, agesex.rename(columns={"sex":"sex2", "age":"age2"}), left_on="ID2", right_index=True)
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
            if age1 >  age2+12:
                if sex1 == "F":
                    p2.mid = p1.id
                if sex1 == "M":
                    p2.pid = p1.id

            if age2 > age1+12:
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
            #default mother id
            p.mid = p.fid + "___M"

        if p.pid is None:
            #default father ir
            p.pid = p.fid + "___P"
        
        data.append((p.fid, p.id, p.pid, p.mid))

    data = pd.DataFrame(data, columns = ['FID' , 'IID', 'FATHER_ID' , 'MOTHER_ID']).astype(str)
    return data
    
def add_control(pedigree):
    """Adds control families to the pedigree table for testing.

    For each family that has two or more siblings and both parents, creates a 3 new familes, one has no parents, one with no mother and one with no father.
    gFID of these families are x+original_fid where x is "_o_", "_p_", "_m_" for these cases: no parent, only has father, only has mother. IIDs are the same in both families.

    Args:
        pedigree : pd.DataFrame
            A pedigree table with 'FID', 'IID', 'FATHER_ID', 'MOTHER_ID', 'has_father', 'has_mother'. Each row represents an individual.
            fids starting with "_" are reserved for control.
        
    Returns:
        pd.DataFrame
            A pedigree table with 'FID', 'IID', 'FATHER_ID', 'MOTHER_ID'. Each row represents an individual.
            For each family with both parents and more than one offspring, it has a control family(fids for control families start with '_')

    """

    families_with_both_parents = pedigree[pedigree["has_father"] & pedigree["has_mother"]]
    count_of_sibs_in_fam = families_with_both_parents.groupby(["FID", "FATHER_ID", "MOTHER_ID"]).count().reset_index()
    FIDs_with_multiple_sibs = count_of_sibs_in_fam[count_of_sibs_in_fam["IID"] > 1][["FID"]]
    families_with_multiple_sibs = families_with_both_parents.merge(FIDs_with_multiple_sibs, on = "FID")

    families_with_multiple_sibs["FID"] = "_o_" + families_with_multiple_sibs["FID"].astype(str)
    families_with_multiple_sibs["MOTHER_ID"] = families_with_multiple_sibs["FID"].astype(str) + "_M"
    families_with_multiple_sibs["FATHER_ID"] = families_with_multiple_sibs["FID"].astype(str) + "_P"
    families_with_multiple_sibs["has_father"] = False
    families_with_multiple_sibs["has_mother"] = False

    keep_mother = families_with_both_parents.copy()
    keep_mother["FID"] = "_m_" + keep_mother["FID"].astype(str)
    keep_mother["FATHER_ID"] = keep_mother["FID"].astype(str) + "_P"
    keep_mother["has_father"] = False
    keep_mother["has_mother"] = True


    keep_father = families_with_both_parents.copy()
    keep_father["FID"] = "_p_" + keep_father["FID"].astype(str)
    keep_father["MOTHER_ID"] = keep_father["FID"].astype(str) + "_M"
    keep_father["has_father"] = True
    keep_father["has_mother"] = False

    pedigree = pedigree.append(families_with_multiple_sibs).append(keep_father).append(keep_mother)
    return pedigree

#TODO raise error if file is multi chrom
def preprocess_king(ibd, segs, bim, chromosomes, sibships):
    """Converts the ibds in king format to ibds in snipar format
        King format only saves ibd1 and ibd2s in the ibd file. The rest is ibd0 only if present in the segs file. This function finds the ibd0 sections and appends to the ibd data structure.
    
    Args:
        ibd: pd.DataFrame
            A pandas DataFrame with columns including ["ID1", "ID2", "IBDType", "Chr", "StartSNP", "StopSNP"] where IDs are individual IIDs.
        
        segs: pd.DataFrame
            A pandas DataFrame with columns including ["Segment", "Chr", "StartSNP", "StopSNP"]

        bim: pd.DataFrame
            A dataframe with these columns(dtype str) including: Chr id coordinate
        
        chromosomes: list
            list of chromosome numbers

        sibships: pandas.DataFrame
            A pandas DataFrame with columns ['FID', 'FATHER_ID', 'MOTHER_ID', 'IID', 'has_father', 'has_mother', 'single_parent'] where IID columns is a list of the IIDs of individuals in that family.
            It only contains families that have more than one child or only one parent.

    Returns:
        (str, str) -> list
            A dictionary where the keys are pairs of individual ids and the values are IBD segments. The list is flattened as has information about successive ibd segments meaning it's like [start0, end0, ibd_status0, start1, end1, ibd_status1, ...]
    """
    ibd["Chr"] = ibd["Chr"].astype(int)
    segs["Chr"] = segs["Chr"].astype(int)
    chromosomes = [int(x) for x in chromosomes]
    if len(chromosomes)>1:
        ibd = ibd[ibd["Chr"].isin(chromosomes)][["ID1", "ID2", "IBDType", "StartSNP", "StopSNP"]]
    else:
        ibd = ibd[ibd["Chr"]==chromosomes[0]][["ID1", "ID2", "IBDType", "StartSNP", "StopSNP"]]
    #TODO cancel or generalize this
    if set(ibd["IBDType"].unique().tolist()) == {"IBD1", "IBD2"}:
        ibd["IBDType"] = ibd["IBDType"].apply(lambda x: 2 if x=="IBD2" else 1)
    ibd["IBDType"] = ibd["IBDType"].astype(int)
    temp = bim[["id", "coordinate"]].rename(columns = {"id":"StartSNP","coordinate":"StartSNPLoc"})
    ibd= ibd.merge(temp, on="StartSNP")
    temp = bim[["id", "coordinate"]].rename(columns = {"id":"StopSNP","coordinate":"StopSNPLoc"})
    ibd = ibd.merge(temp, on="StopSNP")
    ibd['segment'] = ibd[['StartSNPLoc', 'StopSNPLoc', "IBDType"]].values.tolist()
    ibd = ibd.groupby(["ID1", "ID2"]).agg({'segment':sorted}).reset_index()
    if len(chromosomes)>1:
        segs = segs[segs["Chr"].isin(chromosomes)][["StartSNP", "StopSNP"]]
    else:
        segs = segs[segs["Chr"]==chromosomes[0]][["StartSNP", "StopSNP"]]
    temp = bim[["id", "coordinate"]].rename(columns = {"id":"StartSNP","coordinate":"StartSNPLoc"})
    segs= segs.merge(temp, on="StartSNP")
    temp = bim[["id", "coordinate"]].rename(columns = {"id":"StopSNP","coordinate":"StopSNPLoc"})
    segs = segs.merge(temp, on="StopSNP")
    segs = segs[['StartSNPLoc', 'StopSNPLoc']].sort_values('StartSNPLoc').values.tolist()
    flatten_seg_as_ibd0 = []
    for l in segs:
        flatten_seg_as_ibd0 = flatten_seg_as_ibd0 + l + [0]
    #TODO does this work with multichromosome in the ibd file? it won't work if snplocs are indexed from zero in each snp
    all_ibd_segs = []
    for index, row in ibd.iterrows():
        id1, id2, segments = row["ID1"], row["ID2"], row["segment"]
        seg_counter = 0
        row_ibd0 = []
        start, end = segs[seg_counter]
        prev_end = start
        for seg_start, seg_end, ibd_type in segments:  
            while seg_start>end:
                if prev_end<end:
                    row_ibd0.append([prev_end, end, 0])
                if (seg_counter+1) < len(segs):
                    seg_counter+=1
                    start, end = segs[seg_counter]
                    prev_end = start
                else:
                    raise Exception("this segments starts after all meaningfull segments")

            if seg_start<start:
                raise Exception("segment starts sooner than it should")

            if seg_start>prev_end:
                row_ibd0.append([prev_end, seg_start, 0])

            if seg_end>end:
                raise Exception("segment ends outside where it should have")
            prev_end=seg_end

        if prev_end<end:
            row_ibd0.append([prev_end, end, 0])
        row_ibd0 = row_ibd0 + [[start, end, 0] for start, end in segs[seg_counter+1:]]
        row_ibd = segments+row_ibd0
        row_ibd = sorted(row_ibd, key=lambda x:x[0])        
        flatten_row_ibd = []
        for l in row_ibd:
            for el in l:
                flatten_row_ibd.append(el)
        all_ibd_segs.append(flatten_row_ibd)
    ibd["segment"] = pd.Series(all_ibd_segs)
    ibd_dict = ibd.set_index(["ID1", "ID2"]).to_dict()["segment"]
    for index, row in sibships.iterrows():
        sibs = row["IID"]
        nsibs = len(sibs)
        if nsibs > 1:
            for i in range(1, nsibs):
                for j in range(0, i):
                    sib1 = sibs[i].decode()
                    sib2 = sibs[j].decode()
                    if not((sib1, sib2) in ibd_dict or (sib2, sib1) in ibd_dict):
                        ibd_dict[(sib1, sib2)] = flatten_seg_as_ibd0
    return ibd_dict

def prepare_data(pedigree, phased_address, unphased_address, ibd_address, ibd_is_king, bim_address = None, fam_address = None, control = False, chromosome = None, pedigree_nan = '0'):
    """Processes the non_gts required data for the imputation and returns it.

    Outputs for used for the imputation have ascii bytes instead of strings.
    
    Args:
        pedigree : pd.DataFrame 
            The pedigree table. It contains 'FID', 'IID', 'FATHER_ID' and, 'MOTHER_ID' columns.
        
        phased_address : str
            Address of the phased bgen file (does not inlude '.bgen'). Only one of unphased_address and phased_address is neccessary.

        unphased_address : str
            Address of the bed file (does not inlude '.bed'). Only one of unphased_address and phased_address is neccessary.
        
        ibd_address : str
            address of the ibd file. The king segments file should be accompanied with an allsegs file.

        ibd_is_king : boolean
            Whether the ibd segments are in king format or snipar format.


        bim_address : str, optional
            Address of the bim file if it's different from the address of the bed file. Does not include '.bim'.
        
        fam_address : str, optional
            Address of the fam file if it's different from the address of the bed file. Does not include '.fam'.
        
        control : boolean, optional
            If True, adds control families to the pedigree table for testing using snipar.imputation.preprocess_data.add_control.
        
        chromosome: str, optional
            Number of the chromosome that's going to be loaded.

        pedigree_nan: str, optional
            Value that's considered nan in the pedigree. The default is '0'

    Returns:
        tuple(pandas.Dataframe, dict, numpy.ndarray, pandas.Dataframe, numpy.ndarray, numpy.ndarray)
            Returns the data required for the imputation. This data is a tuple of multiple objects.
                sibships: pandas.DataFrame
                    A pandas DataFrame with columns ['FID', 'FATHER_ID', 'MOTHER_ID', 'IID', 'has_father', 'has_mother', 'single_parent'] where IID columns is a list of the IIDs of individuals in that family.
                    It only contains families that have more than one child or only one parent.

                ibd: pandas.DataFrame
                    A pandas DataFrame with columns "ID1", "ID2", 'segment'. The segments column is a list of IBD segments between ID1 and ID2.
                    Each segment consists of a start, an end, and an IBD status. The segment list is flattened meaning it's like [start0, end0, ibd_status0, start1, end1, ibd_status1, ...]

                bim: pandas.DataFrame
                    A dataframe with these columns(dtype str): Chr id morgans coordinate allele1 allele2

                chromosomes: str
                    A string containing all the chromosomes present in the data.

                ped_ids: set
                    Set of ids of individuals with missing parents.

                pedigree_output: np.array
                    Pedigree with added parental status.
    """    
    logging.info("For file "+str(phased_address)+";"+str(unphased_address)+": Finding which chromosomes")
    if unphased_address:
        if bim_address is None:
            bim_address = unphased_address+'.bim'
        bim = pd.read_csv(bim_address, delim_whitespace=True, header=None, names=["Chr", "id", "morgans", "coordinate", "allele1", "allele2"])
        if fam_address is None:
            fam_address = unphased_address+'.fam'
        fam = pd.read_csv(fam_address, delim_whitespace=True, header=None, names=["FID", "IID", "PID", "MID", "SEX", "Phen"])

    else:
        if bim_address is None:
            bim_address = phased_address+'.bgen'
        if fam_address is None:
            fam_address = phased_address+'.bgen'
        bgen = read_bgen(bim_address, verbose=False)
        bim = bgen["variants"].compute().rename(columns={"chrom":"Chr", "pos":"coordinate"})
        #TODO this line should be replaced
        bim["id"]=bim["rsid"]
        if chromosome is None:
            raise Exception("chromosome should be specified when using phased data") 
        bim["Chr"] = chromosome
        bgen = read_bgen(bim_address, verbose=False)
        bim = bgen["variants"].compute().rename(columns={"chrom":"Chr", "pos":"coordinate"})
        #TODO this line should be replaced
        bim["id"]=bim["rsid"]
        if chromosome is None:
            raise Exception("chromosome should be specified when using phased data") 
        bim["Chr"] = chromosome
        fam = pd.DataFrame({"IID":read_bgen(fam_address)["samples"]})

    chromosomes = bim["Chr"].unique().astype(int)
    logging.info(f"with chromosomes {chromosomes} initializing non_gts data")
    logging.info(f"with chromosomes {chromosomes} loading and filtering pedigree file ...")
    #Keep individuals in fam
    pedigree = pedigree[pedigree["IID"].isin(fam["IID"].astype(str))].copy()
    pedigree["has_father"] = pedigree["FATHER_ID"].isin(fam["IID"].astype(str))
    pedigree["has_mother"] = pedigree["MOTHER_ID"].isin(fam["IID"].astype(str))
    if control:
        logging.info("Adding control to the pedigree ...")
        pedigree = add_control(pedigree)
        logging.info("Control Added.")
    #keeping individuals with no parents
    no_parent_pedigree = pedigree[~(pedigree["has_mother"] & pedigree["has_father"])]
    #removing individual whose parents are nan
    no_parent_pedigree = no_parent_pedigree[(no_parent_pedigree["MOTHER_ID"] != pedigree_nan) & (no_parent_pedigree["FATHER_ID"] != pedigree_nan)]
    no_parent_pedigree[["FID", "IID", "FATHER_ID", "MOTHER_ID"]] = no_parent_pedigree[["FID", "IID", "FATHER_ID", "MOTHER_ID"]].astype("S")
    ped_ids =  set(no_parent_pedigree["IID"].tolist())
    #finding siblings in each family
    sibships = no_parent_pedigree.groupby(["FID", "FATHER_ID", "MOTHER_ID", "has_father", "has_mother"]).agg({'IID':lambda x: list(x)}).reset_index()

    sibships["sib_count"] = sibships["IID"].apply(len)
    sibships["single_parent"] = sibships["has_father"] ^ sibships["has_mother"]
    sibships = sibships[(sibships["sib_count"]>1) | sibships["single_parent"]]
    fids = set([i for i in sibships["FID"].values.tolist() if i.startswith(b"_")])
    logging.info(f"with chromosomes {chromosomes} loading bim file ...")      
    logging.info(f"with chromosomes {chromosomes} loading and transforming ibd file ...")
    if ibd_address is None:
        logging.info(f"with chromosomes {chromosomes} making an empty ibd dataframe")
        if (sibships["sib_count"]>1).any():
            raise Exception("Should provide ibd file in the presense of families with multiple sibs")
        ibd = {}
    else:
        ibd = pd.read_csv(f"{ibd_address}.segments.gz", delim_whitespace=True).astype(str)
        if ibd_is_king:            
            king_segs = pd.read_csv(f"{ibd_address}allsegs.txt", delim_whitespace=True).astype(str)                        
            if not {"ID1", "ID2", "IBDType", "Chr", "StartSNP", "StopSNP",}.issubset(set(ibd.columns.values.tolist())):
                raise Exception("Invalid ibd columns for king formatted ibd. Columns must include: ID1, ID2, IBDType, Chr, StartSNP, StopSNP")
            ibd = preprocess_king(ibd, king_segs, bim, chromosomes, sibships)
        else:
            if not {"ID1", "ID2", "IBDType", "Chr", "start_coordinate", "stop_coordinate",}.issubset(set(ibd.columns.values.tolist())):
                raise Exception("Invalid ibd columns for snipar formatted ibd. Columns must be include: ID1, ID2, IBDType, Chr, start_coordinate, stop_coordinate")
            ibd = ibd.astype(str)
            ibd[["IBDType", "start_coordinate", "stop_coordinate"]] = ibd[["IBDType", "start_coordinate", "stop_coordinate"]].astype(int)
            #Adding location of start and end of each
            chromosomes = chromosomes.astype(str)
            if len(chromosomes)>1:
                ibd = ibd[ibd["Chr"].isin(chromosomes)]
            else:
                ibd = ibd[ibd["Chr"]==chromosomes[0]]
            #TODO cancel or generalize this    
            ibd['segment'] = ibd[['start_coordinate', 'stop_coordinate', "IBDType"]].values.tolist()
            def create_seg_list(x):
                elements = list(x)
                result = []
                for el in elements:
                    result = result+el
                return result
            ibd = ibd.groupby(["ID1", "ID2"]).agg({'segment':lambda x:create_seg_list(x)}).to_dict()["segment"]
    if not ibd:
        logging.warning(f"with chromosomes {chromosomes} no matching ibd segments")        
    logging.info("ibd loaded.")        
    logging.info(f"with chromosomes {chromosomes} initializing non_gts data done ...")
    pedigree[["FID", "IID", "FATHER_ID", "MOTHER_ID"]] = pedigree[["FID", "IID", "FATHER_ID", "MOTHER_ID"]].astype(str)
    pedigree_output = np.concatenate(([pedigree.columns.values.tolist()], pedigree.values)).astype('S')
    return sibships, ibd, bim, chromosomes, ped_ids, pedigree_output

def compute_aics(unphased_gts, pc_scores, linear=True, sample_size = 1000):
    """Akaike information criterion of linear regressions with increasing number of PCs. Returns the number of PCs that minimizes aic.
    Args:    
        unphased_gts: np.array[signed char]
            A two-dimensional array containing genotypes for all individuals and SNPs respectively.

        pc_scores: np.array[float]
            A two-dimensional array containing pc scores for all individuals and SNPs respectively.

        linear : bool, optional
            Whether the model is linear regression or not. GLM is not implemented yet. Default is true.
        
        sample_size : int, optional
            number of snps to use for computing aics. aics from these snps are averaged. default is 1000.
    
    Returns:
        int:
            optimal number of PCs to use
    """
    if linear:
        pop_size = unphased_gts.shape[0]
        nsnps = unphased_gts.shape[1]
        if nsnps < sample_size:
            sample_size = nsnps
        individual_mask = np.random.permutation([True]*sample_size+[False]*(nsnps-sample_size))
        unphased_gts = unphased_gts[:, individual_mask]

        unphased_pc_gts = unphased_gts.copy().astype(np.float)
        unphased_pc_gts[unphased_gts==nan_integer]=0.
        number_of_pcs = pc_scores.shape[1]
        x = np.hstack((np.ones((pop_size, 1)), pc_scores))
        aic = []
        logging.info("computing aics ...")
        for i in tqdm(range(1,number_of_pcs+2)):
            coefs, _, _, _ = np.linalg.lstsq(x[:, :i], unphased_pc_gts/2)
            fs = x[:, :i]@coefs
            std = np.sqrt(np.var(unphased_pc_gts.flatten()/2)-np.var(fs.flatten()))
            fs[fs>1] = 1
            fs[fs<0] = 0            
            log_likelihoods = np.sum(np.log(norm.pdf(unphased_pc_gts/2-fs, 0, std)), axis=0)
            log_likelihood = np.mean(log_likelihoods)
            aic.append(2*(i-log_likelihood))
        logging.info("computing aics done")
        return np.argmin(aic)
    raise Exception("not implemented yet")

def estimate_f(unphased_gts, pc_scores, linear=True):
    """Estimates MAF with an ols or glm from by regressing unphased_gts on pc_scores
    
    Args:
        unphased_gts: np.array[signed char]
            A two-dimensional array containing genotypes for all individuals and SNPs respectively.

        pc_scores: np.array[float]
            A two-dimensional array containing pc scores for all individuals and SNPs respectively.

        linear, bool, optional
            Whether the model is linear regression or not. Default is true.
        
    Returns:
        np.array[float], dict
            A two-dimensional array containing estimated fs for all individuals and SNPs respectively and a dictionary containing information about the model. These include ['x', 'coefs', 'TSS', 'RSS1', 'RSS2', 'R2_1', 'R2_2', 'larger1', 'less0'].
    """

    pop_size = unphased_gts.shape[0]
    nsnps = unphased_gts.shape[1]
    unphased_pc_gts = unphased_gts.copy().astype(np.float)
    unphased_pc_gts[unphased_gts==nan_integer] = 0.
    x = np.hstack((pc_scores, np.ones((pop_size, 1))))
    data = {}
    if linear:
        y = unphased_pc_gts/2
        coefs, _, _, _ = np.linalg.lstsq(x, y)
        fs = x@coefs
        residuals1 = y-fs
        larger1 = np.sum(fs>1., axis=0)/pop_size
        less0 = np.sum(fs<0., axis=0)/pop_size
        fs[fs>1] = 1
        fs[fs<0] = 0
        residuals2 = y-fs
        TSS = np.sum(y*y, axis=0)
        RSS1 = np.sum(residuals1*residuals1, axis=0)
        RSS2 = np.sum(residuals2*residuals2, axis=0)
        data["x"] = x
        data["coefs"] = coefs
        data["TSS"] = TSS
        data["RSS1"] = RSS1
        data["RSS2"] = RSS2
        data["R2_1"] = 1-(RSS1/TSS)
        data["R2_2"] = 1-(RSS2/TSS)
        data["larger1"] = larger1
        data["less0"] = less0
    else:
        models = []
        model_results = []
        fs = unphased_pc_gts.copy()
        for i in tqdm(range(nsnps)):
            y = np.zeros((pop_size,2))
            y[:,0] = unphased_pc_gts[:,i]
            y[:,1] = 2-unphased_pc_gts[:,i]
            binom_model = sm.GLM(y, x, family=sm.families.Binomial())
            binom_result = binom_model.fit()
            models.append(binom_model)
            model_results.append(binom_result)
            fs[:,i] = binom_result.predict(x)
    return fs, data


def prepare_gts(phased_address, unphased_address, bim, pedigree_output, ped_ids, chromosomes, start=None, end=None, pcs=None, pc_ids=None, find_optimal_pc=None):
    """ Processes the gts required data for the imputation and returns it.

    Outputs for used for the imputation have ascii bytes instead of strings.

    Args:
        phased_address : str
            Address of the phased bgen file (does not inlude '.bgen'). Only one of unphased_address and phased_address is neccessary.

        unphased_address : str
            Address of the bed file (does not inlude '.bed'). Only one of unphased_address and phased_address is neccessary.

        bim: pandas.DataFrame
            A dataframe with these columns(dtype str): Chr id morgans coordinate allele1 allele2

        pedigree_output: np.array
            Pedigree with added parental status.

        ped_ids: set
            Set of ids of individuals with missing parents.
        
        chromosomes: str
                    A string containing all the chromosomes present in the data.

        start : int, optional
            This function can be used for preparing a slice of a chromosome. This is the location of the start of the slice.

        end : int, optional
            This function can be used for preparing a slice of a chromosome. This is the location of the end of the slice.

        pcs : np.array[float], optional
            A two-dimensional array containing pc scores for all individuals and SNPs respectively.

        pc_ids : set, optional
            Set of ids of individuals in the pcs.

        find_optimal_pc : bool, optional
            It will use Akaike information criterion to find the optimal number of PCs to use for MAF estimation.

    Returns:
        tuple(np.array[signed char], np.array[signed char], str->int, np.array[int], np.array[float], dict)
            phased_gts: np.array[signed char], optional
                A three-dimensional array containing genotypes for all individuals, SNPs and, haplotypes respectively.

            unphased_gts: np.array[signed char]
                A two-dimensional array containing genotypes for all individuals and SNPs respectively.

            iid_to_bed_index: str->int
                A str->int dictionary mapping IIDs of people to their location in bed file.

            pos: np.array[int]
                A numpy array with the position of each SNP in the order of appearance in gts.

            practical_f: np.array[float]
                A two-dimensional array containing estimated fs for all individuals and SNPs respectively.

            hdf5_output_dict: dict
                A  dictionary whose values will be written in the imputation output under its keys. It contains:
                    'bim_columns' : Columns of the resulting bim file
                    'bim_values' : Contents of the resulting bim file
                    'pedigree' : pedigree table Its columns are has_father, has_mother, single_parent respectively.
                    'non_duplicates' : Indexes of the unique snps. Imputation is restricted to them.
                    'standard_f' : Whether the allele frequencies are just population average instead of MAFs estimated using PCs
    """
    logging.info(f"with chromosomes {chromosomes} initializing gts data with start={start} end={end}")
    phased_gts = None
    unphased_gts = None
    unphased_pc_gts = None    
    if unphased_address:
        bim_as_csv = pd.read_csv(unphased_address+".bim", delim_whitespace=True, header=None)
        gts_f = Bed(unphased_address+".bed",count_A1 = True, sid=bim_as_csv[1].values.tolist())
        logging.info(f"with chromosomes {chromosomes} opened unphased file ...")
        if not pc_ids is None:
            ids_in_ped_pc = [(id in ped_ids) and (id in pc_ids) for id in gts_f.iid[:,1].astype("S")]
        else:
            ids_in_ped_pc = [(id in ped_ids) for id in gts_f.iid[:,1].astype("S")]
        logging.info(f"with chromosomes {chromosomes} loaded ids ...")
        gts_ids = gts_f.iid[ids_in_ped_pc]
        logging.info(f"with chromosomes {chromosomes} restrict to ids ...")
        all_sids = bim_as_csv[1].values

        if end is not None:
            unphased_gts = gts_f[ids_in_ped_pc , start:end].read().val
            logging.info(f"with chromosomes {chromosomes} loaded genotypes ...")
            pos = gts_f.pos[start:end, 2]
            logging.info(f"with chromosomes {chromosomes} loaded pos ...")
            sid = all_sids[start:end]
            logging.info(f"with chromosomes {chromosomes} loaded sid ...")
        else:
            unphased_gts = gts_f[ids_in_ped_pc, :].read().val
            logging.info(f"with chromosomes {chromosomes} loaded genotypes ...")
            pos = gts_f.pos[:, 2]
            logging.info(f"with chromosomes {chromosomes} loaded pos ...")
            sid = all_sids
            logging.info(f"with chromosomes {chromosomes} loaded sid ...")
    if phased_address:
        bgen = open_bgen(phased_address+".bgen", verbose=False)
        if not pc_ids is None:
            ids_in_ped_pc = [(id in ped_ids) and (id in pc_ids) for id in bgen.samples.astype("S")]
        else:
            ids_in_ped_pc = [(id in ped_ids) for id in bgen.samples.astype("S")]
        gts_ids = np.array([[None, id] for id in bgen.samples])
        gts_ids = gts_ids[ids_in_ped_pc]
        pos = np.array(bim["coordinate"][start: end])
        all_sids = np.array(bim["id"])
        sid = np.array(bim["id"][start: end])
        pop_size = len(gts_ids)
        probs= bgen.read((slice(0, len(bgen.samples)),slice(start, end)))        
        phased_gts = np.zeros((probs.shape[0], probs.shape[1], 2))
        phased_gts[:] = np.nan
        phased_gts[probs[:,:,0] > 0.99, 0] = 1
        phased_gts[probs[:,:,1] > 0.99, 0] = 0
        phased_gts[probs[:,:,2] > 0.99, 1] = 1
        phased_gts[probs[:,:,3] > 0.99, 1] = 0
        phased_gts = phased_gts[ids_in_ped_pc,:]
        if not unphased_gts:
            unphased_gts = phased_gts[:,:,0]+phased_gts[:,:,1]
    if not pcs is None:
        set_gts_ids = set(gts_ids[:, 1].astype("S"))
        pc_ids_gts_ids =  [id in set_gts_ids for id in pc_ids]
        pcs = pcs[pc_ids_gts_ids]
    _, indexes = np.unique(all_sids, return_index=True)
    indexes = np.sort(indexes)
    indexes = (indexes[(indexes>=start) & (indexes<end)]-start)
    sid = sid[indexes]
    pos = pos[indexes]
    unphased_gts = unphased_gts[:, indexes]
    if not phased_gts is None:
        phased_gts = phased_gts[:, indexes, :]        
    pos = pos.astype(int)
    unphased_gts_greater2 = unphased_gts>2
    num_unphased_gts_greater2 = np.sum(unphased_gts_greater2)
    if num_unphased_gts_greater2>0:
        logging.warning(f"with chromosomes {chromosomes}: unphased genotypes are greater than 2 in {num_unphased_gts_greater2} locations. Converted to NaN")  
        unphased_gts[unphased_gts_greater2] = np.nan
        
    unphased_gts_less0 = unphased_gts<0
    num_unphased_gts_less0 = np.sum(unphased_gts_less0)
    if num_unphased_gts_less0>0:
        logging.warning(f"with chromosomes {chromosomes}: unphased genotypes are less than 0 in {num_unphased_gts_less0} locations. Converted to NaN")  
        unphased_gts[unphased_gts_less0] = np.nan

    if not phased_gts is None:
        phased_gts_less0 = phased_gts<0
        num_phased_gts_less0 = np.sum(phased_gts_less0)
        if num_phased_gts_less0>0:
            logging.warning(f"with chromosomes {chromosomes}: phased genotypes are less than 0 in {num_phased_gts_less0} locations. Converted to NaN")  
            phased_gts[phased_gts_less0] = np.nan
        
        phased_gts_greater1 = phased_gts>1
        num_phased_gts_greater1 = np.sum(phased_gts_greater1)
        if num_phased_gts_greater1>0:
            logging.warning(f"with chromosomes {chromosomes}: phased genotypes are greater than 1 in {num_phased_gts_greater1} locations. Converted to NaN")  
            phased_gts[phased_gts_greater1] = np.nan
    
    standard_f = np.nanmean(unphased_gts,axis=0)/2.0
    practical_f = unphased_gts.copy()
    practical_f[:] = standard_f

    #transforming genotypes into int8
    unphased_gts[np.isnan(unphased_gts)] = nan_integer
    unphased_gts = unphased_gts.astype(np.int8)
    if not phased_gts is None:
        phased_gts[np.isnan(phased_gts)] = nan_integer
        phased_gts = phased_gts.astype(np.int8)    
    iid_to_bed_index = {i.encode("ASCII"):index for index, i in enumerate(gts_ids[:,1])}
    selected_bim = bim.iloc[indexes+start, :]
    bim_values = selected_bim.to_numpy().astype('S')
    bim_columns = selected_bim.columns.to_numpy().astype('S')
    logging.info(f"with chromosomes {chromosomes} initializing non_gts data done")
    if unphased_gts.shape[0] == 0:
        logging.warning(f"with chromosomes {chromosomes}: 0 individuals in unphased genotypes")
    if unphased_gts.shape[1] == 0:
        logging.warning(f"with chromosomes {chromosomes}: 0 snps in unphased genotypes")
    if unphased_gts.max()>2 or unphased_gts.min()<0:
        logging.warning(f"with chromosomes {chromosomes}: unphased genotypes has values out of 0-2 range")
    if not phased_gts is None:
        if phased_gts.shape[0] == 0:
            logging.warning(f"with chromosomes {chromosomes}: 0 individuals in phased genotypes")
        if phased_gts.shape[1] == 0:
            logging.warning(f"with chromosomes {chromosomes}: 0 snps in phased genotypes")
        if phased_gts.max()>1 or phased_gts.min()<0:
            logging.warning(f"with chromosomes {chromosomes}: phased genotypes has values out of 0-1 range")
    #value types should be hdf5 compatible
    hdf5_output_dict = {"bim_columns":bim_columns,
                        "bim_values":bim_values,
                        "pedigree":pedigree_output,
                        "non_duplicates":indexes,
                        "standard_f":standard_f,
                        }
    if not pcs is None:
        logging.info("estimating MAF using PCs ...")
        if find_optimal_pc:
            logging.info("finding optimal number of PCs using aic...")
            optimal_pc_num = compute_aics(unphased_gts, pcs, sample_size=1000)
            logging.info(f"optimal number of PCs is {optimal_pc_num}")
            logging.info("finding optimal number of PCs using aic done")
            practical_f, maf_estimator_data = estimate_f(unphased_gts, pcs[:, :optimal_pc_num])
        else:
            practical_f, maf_estimator_data = estimate_f(unphased_gts, pcs)
        for key, val in maf_estimator_data.items():
            hdf5_output_dict[f"maf_{key}"] = val
        logging.info("estimating MAF using PCs done")
    return phased_gts, unphased_gts, iid_to_bed_index, pos, practical_f, hdf5_output_dict