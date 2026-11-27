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
from array import array
from contextlib import nullcontext
from inspect import signature
from itertools import chain
import logging
import pandas as pd
import numpy as np
from pysnptools.snpreader import Bed
from snipar.utilities import open_bgen, read_bgen
from snipar.config import nan_integer
from tqdm import tqdm

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

def add_control(pedigree):
    """Adds control families to the pedigree table for testing.

    For each family that has two or more siblings and both parents, creates a 3 new families, one has no parents, one with no mother and one with no father.
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

    pedigree = pd.concat([pedigree, families_with_multiple_sibs, keep_father, keep_mother])
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
    ibd["IBDType"] = ibd["IBDType"].replace({"IBD1": "1", "IBD2": "2"})
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
    flatten_seg_as_ibd0 = list(chain.from_iterable([start, end, 0] for start, end in segs))
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
                if prev_end <= end:
                    row_ibd0.append([prev_end, end, 0])
                if (seg_counter+1) < len(segs):
                    seg_counter+=1
                    start, end = segs[seg_counter]
                    prev_end = start
                else:
                    raise Exception("this segments starts after all meaningful segments")

            if seg_start<start:
                raise Exception("segment starts sooner than it should")

            if seg_start>prev_end:
                row_ibd0.append([prev_end, seg_start - 1, 0])

            if seg_end>end:
                raise Exception("segment ends outside where it should have")
            prev_end = seg_end + 1

        if prev_end <= end:
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
            Address of the phased bgen file (does not include '.bgen'). Only one of unphased_address and phased_address is necessary.

        unphased_address : str
            Address of the bed file (does not include '.bed'). Only one of unphased_address and phased_address is necessary.
        
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
        bim = pd.read_csv(bim_address, sep=r"\s+", header=None, names=["Chr", "id", "morgans", "coordinate", "allele1", "allele2"],
                          dtype={"id": str, "allele1": str, "allele2": str}, keep_default_na=False)
        if fam_address is None:
            fam_address = unphased_address+'.fam'
        fam = pd.read_csv(fam_address, sep=r"\s+", header=None, names=["FID", "IID", "PID", "MID", "SEX", "Phen"],
                          usecols=["IID"], dtype={"IID": str}, keep_default_na=False)

    else:
        if bim_address is None:
            bim_address = phased_address+'.bgen'
        if fam_address is None:
            fam_address = phased_address+'.bgen'
        with open_bgen(bim_address) as reader:
            if chromosome is None:
                chroms = np.unique(reader.chromosomes)
                if len(chroms) != 1 or not str(chroms[0]).isdigit():
                    raise ValueError("Supply a single BGEN chromosome, or use a chromosome wildcard")
                chromosome = chroms[0]
            if hasattr(reader, "variant_metadata"):
                bim = reader.variant_metadata()
            elif hasattr(reader, "_vaddr"):
                # cbgen already loaded these arrays when opening the BGEN.
                # Copy them here before its mapped index is closed.
                bim = pd.DataFrame({
                    "id": reader.ids, "rsid": reader.rsids, "chrom": reader.chromosomes,
                    "pos": reader.positions, "nalleles": reader.nalleles,
                    "allele_ids": reader.allele_ids, "vaddr": reader._vaddr,
                }, copy=True)
            else:
                bim = read_bgen(bim_address)["variants"].compute()
            ids = np.array(reader.samples)
        bim = bim.rename(columns={"chrom": "Chr", "pos": "coordinate"})
        bim["id"] = bim["rsid"]
        bim["Chr"] = chromosome
        if fam_address != bim_address:
            with open_bgen(fam_address) as reader:
                ids = np.array(reader.samples)
        fam = pd.DataFrame({"IID": ids})

    chromosomes = bim["Chr"].unique().astype(int)
    if len(chromosomes) != 1:
        raise ValueError("Imputation requires one chromosome per genotype file")
    pedigree = pedigree.copy()
    id_columns = ["FID", "IID", "FATHER_ID", "MOTHER_ID"]
    pedigree[id_columns] = pedigree[id_columns].astype(str)
    genotype_ids = set(fam["IID"].astype(str))
    logging.info(f"with chromosomes {chromosomes} initializing non_gts data")
    logging.info(f"with chromosomes {chromosomes} loading and filtering pedigree file ...")
    #Keep individuals in fam
    pedigree = pedigree[pedigree["IID"].isin(genotype_ids)].copy()
    pedigree["has_father"] = pedigree["FATHER_ID"].isin(genotype_ids) & (pedigree["FATHER_ID"] != pedigree_nan)
    pedigree["has_mother"] = pedigree["MOTHER_ID"].isin(genotype_ids) & (pedigree["MOTHER_ID"] != pedigree_nan)
    del genotype_ids, fam
    if control:
        logging.info("Adding control to the pedigree ...")
        pedigree = add_control(pedigree)
        logging.info("Control Added.")
    #keeping individuals with no parents
    no_parent_pedigree = pedigree[~(pedigree["has_mother"] & pedigree["has_father"])]
    # Keep a singleton with one observed parent, even when the missing parent's ID is unknown.
    has_known_parents = ((no_parent_pedigree["MOTHER_ID"] != pedigree_nan) &
                         (no_parent_pedigree["FATHER_ID"] != pedigree_nan))
    single_observed = no_parent_pedigree["has_father"] ^ no_parent_pedigree["has_mother"]
    no_parent_pedigree = no_parent_pedigree.loc[has_known_parents | single_observed].copy()
    no_parent_pedigree[["FID", "IID", "FATHER_ID", "MOTHER_ID"]] = no_parent_pedigree[["FID", "IID", "FATHER_ID", "MOTHER_ID"]].astype("S")
    #finding siblings in each family
    group_keys = ["FID", "FATHER_ID", "MOTHER_ID", "has_father", "has_mother"]
    grouped = no_parent_pedigree.groupby(group_keys)
    sibships = grouped.size().rename("sib_count").reset_index()
    sibships["single_parent"] = sibships["has_father"] ^ sibships["has_mother"]
    eligible = (sibships["sib_count"] > 1) | sibships["single_parent"]
    group_index = grouped.ngroup().to_numpy()
    selected = eligible.to_numpy()[group_index]
    rows = np.argsort(group_index[selected], kind="stable")
    iid = no_parent_pedigree["IID"].to_numpy()[selected][rows]
    sibships = sibships.loc[eligible].copy()
    ends = np.cumsum(sibships["sib_count"].to_numpy())
    groups = np.empty(len(ends), dtype=object)
    for index, (first, end) in enumerate(zip(np.r_[0, ends], ends)):
        groups[index] = iid[first:end]
    sibships["IID"] = groups
    ped_ids = set(no_parent_pedigree["IID"])
    ped_ids.update(sibships.loc[sibships["has_father"], "FATHER_ID"])
    ped_ids.update(sibships.loc[sibships["has_mother"], "MOTHER_ID"])
    del grouped, no_parent_pedigree
    logging.info(f"with chromosomes {chromosomes} loading bim file ...")      
    logging.info(f"with chromosomes {chromosomes} loading and transforming ibd file ...")
    if ibd_address is None:
        logging.info(f"with chromosomes {chromosomes} making an empty ibd dataframe")
        if (sibships["sib_count"]>1).any():
            raise Exception("Should provide ibd file in the presence of families with multiple sibs")
        ibd = {}
    else:
        segments_path = f"{ibd_address}.segments.gz"
        if ibd_is_king:
            ibd = pd.read_csv(segments_path, sep=r"\s+", dtype=str, keep_default_na=False)
            king_segs = pd.read_csv(f"{ibd_address}allsegs.txt", sep=r"\s+", dtype=str, keep_default_na=False)
            if not {"ID1", "ID2", "IBDType", "Chr", "StartSNP", "StopSNP"}.issubset(ibd.columns):
                raise ValueError("KING IBD columns must include ID1, ID2, IBDType, Chr, StartSNP, StopSNP")
            ibd = preprocess_king(ibd, king_segs, bim, chromosomes, sibships)
        else:
            columns = ["ID1", "ID2", "start_coordinate", "stop_coordinate", "IBDType"]
            types = dict.fromkeys(columns[:2] + ["Chr"], str)
            types.update(dict.fromkeys(columns[2:], int))
            ibd = {}
            # The compact native plan keeps only IBD for eligible families. At
            # ingestion, discard segments whose samples are absent entirely.
            genotype_ids = {iid.decode() for iid in ped_ids}
            with pd.read_csv(segments_path, sep=r"\s+", dtype=types,
                             keep_default_na=False, usecols=[*columns, "Chr"], chunksize=200_000) as parts:
                for part in parts:
                    selected = part["Chr"].isin(chromosomes.astype(str))
                    selected &= part["ID1"].isin(genotype_ids) & part["ID2"].isin(genotype_ids)
                    for first, second, start, stop, ibd_type in part.loc[selected, columns].itertuples(
                            index=False, name=None):
                        key = (first, second)
                        segments = ibd.get(key)
                        if segments is None:
                            # Pack coordinates and states while reading the
                            # file so segment RAM does not scale with Python
                            # object overhead for every integer.
                            ibd[key] = segments = array('i')
                        segments.extend((start, stop, ibd_type))
    if not ibd:
        logging.warning(f"with chromosomes {chromosomes} no matching ibd segments")        
    logging.info("ibd loaded.")        
    logging.info(f"with chromosomes {chromosomes} initializing non_gts data done ...")
    pedigree_output = np.concatenate(([pedigree.columns.values.tolist()], pedigree.values)).astype('S')
    return sibships, ibd, bim, chromosomes, ped_ids, pedigree_output

def _fit_frequency_model(x, dosage, called):
    """Fit SNP regressions using only called individuals."""
    coefs = np.zeros((x.shape[1], dosage.shape[1]))
    complete = called.all(axis=0)
    if complete.any():
        coefs[:, complete] = np.linalg.lstsq(x, dosage[:, complete], rcond=None)[0]
    for marker in np.flatnonzero(~complete):
        rows = called[:, marker]
        if rows.any():
            coefs[:, marker] = np.linalg.lstsq(x[rows], dosage[rows, marker], rcond=None)[0]
        else:
            coefs[:, marker] = np.nan
    return coefs


def frequency_model_aic(genotypes, pc_scores):
    """Return additive AIC sums and the number of comparable SNPs."""
    called = genotypes != nan_integer
    dosage = np.where(called, genotypes / 2.0, 0)
    observations = called.sum(axis=0)
    x = np.column_stack((np.ones(len(pc_scores)), pc_scores))
    # Compare the same observations across all candidate models.
    compare = observations > x.shape[1]
    count = int(np.count_nonzero(compare))
    aic = np.zeros(x.shape[1])
    if count:
        for index, n in enumerate(range(1, x.shape[1] + 1)):
            design = x[:, :n]
            coefs = _fit_frequency_model(design, dosage, called)
            residuals = np.where(called, dosage - design @ coefs, 0)
            mse = np.sum(residuals[:, compare] ** 2, axis=0) / observations[compare]
            aic[index] = np.sum(observations[compare] * np.log(np.maximum(mse, np.finfo(float).eps)))
    return aic, count


def compute_aics(unphased_gts, pc_scores, linear=True, sample_size=1000):
    """Select PC count by mean Gaussian AIC across evenly spaced reference SNPs."""
    if not linear:
        raise NotImplementedError("AIC selection is implemented for linear models")
    nsnps = unphased_gts.shape[1]
    sample = np.linspace(0, nsnps - 1, min(sample_size, nsnps), dtype=int)
    aic, count = frequency_model_aic(unphased_gts[:, sample], pc_scores)
    return int(np.argmin(aic + count * 2 * np.arange(1, len(aic) + 1)))


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

    pop_size, nsnps = unphased_gts.shape
    called = unphased_gts != nan_integer
    x = np.column_stack((pc_scores, np.ones(pop_size)))
    data = {}
    if linear:
        y = np.where(called, unphased_gts / 2.0, 0)
        coefs = _fit_frequency_model(x, y, called)
        fs = x @ coefs
        residuals1 = np.where(called, y - fs, 0)
        larger1 = np.mean(fs > 1, axis=0)
        less0 = np.mean(fs < 0, axis=0)
        np.clip(fs, 0, 1, out=fs)
        residuals2 = np.where(called, y - fs, 0)
        with np.errstate(invalid="ignore", divide="ignore"):
            centered = np.where(called, y - y.sum(axis=0) / called.sum(axis=0), 0)
            tss = np.sum(centered**2, axis=0)
            rss1 = np.sum(residuals1**2, axis=0)
            rss2 = np.sum(residuals2**2, axis=0)
            r2_1 = np.where(tss > 0, 1 - rss1/tss, np.nan)
            r2_2 = np.where(tss > 0, 1 - rss2/tss, np.nan)
        data.update(x=x, coefs=coefs, TSS=tss, RSS1=rss1, RSS2=rss2,
                    R2_1=r2_1, R2_2=r2_2, larger1=larger1, less0=less0)
    else:
        import statsmodels.api as sm

        fs = np.full((pop_size, nsnps), np.nan)
        for i in tqdm(range(nsnps)):
            observed = called[:, i]
            if observed.any():
                dosage = unphased_gts[observed, i]
                y = np.column_stack((dosage, 2 - dosage))
                model = sm.GLM(y, x[observed], family=sm.families.Binomial())
                fs[:, i] = model.fit().predict(x)
    return fs, data


class ImputationGenotypes:
    """Open genotype files once and cache chromosome-level sample and allele alignment."""

    def __init__(self, phased_address, unphased_address, bim, ped_ids, pcs=None, pc_ids=None,
                 *, read_threads=1, probability_bytes=128 * 1024**2):
        self.bed = Bed(unphased_address + ".bed", count_A1=True) if unphased_address else None
        self.bgen = None
        self.n_samples = 0
        self.read_threads = read_threads
        self.probability_bytes = probability_bytes
        try:
            if phased_address:
                self.bgen = open_bgen(phased_address + ".bgen", verbose=False)
            if self.bed is not None and self.bed.sid_count != len(bim):
                raise ValueError("BIM and BED datasets must have the same number of SNPs")
            all_ids = self.bed.iid[:, 1] if self.bed is not None else self.bgen.samples
            all_ids = np.asarray(all_ids, dtype="S")
            self.sample_indices = np.flatnonzero(pd.Index(all_ids).isin(ped_ids))
            if not len(self.sample_indices):
                raise ValueError("No relevant individuals have genotypes")
            genotype_ids = all_ids[self.sample_indices]
            self.n_samples = len(genotype_ids)
            self.id_rows = {iid: i for i, iid in enumerate(genotype_ids)}
            if self.bgen is not None:
                self._prepare_bgen(bim, genotype_ids)
            self.matched_pcs = None
            if pcs is not None:
                pc_rows = {iid: i for i, iid in enumerate(pc_ids)}
                indices = np.fromiter((pc_rows.get(iid, -1) for iid in genotype_ids),
                                      dtype=np.intp, count=self.n_samples)
                self.pc_present = indices >= 0
                if not self.pc_present.any():
                    raise ValueError("No principal component scores match genotyped individuals")
                self.matched_pcs = np.ascontiguousarray(pcs[indices[self.pc_present]])
        except BaseException:
            self.close()
            raise

    def _prepare_bgen(self, bim, genotype_ids):
        self.threaded_bgen = "num_threads" in signature(self.bgen.read).parameters
        if self.bed is None:
            self.bgen_indices = self.sample_indices
            self.bgen_output_rows = np.arange(self.n_samples, dtype=np.intp)
            self.bgen_markers = np.arange(len(bim), dtype=np.intp)
            self.flipped = np.zeros(len(bim), dtype=np.uint8)
            self.alleles_match = np.ones(len(bim), dtype=bool)
        else:
            sample_index = {iid: i for i, iid in enumerate(np.asarray(self.bgen.samples, dtype="S"))}
            rows = np.fromiter((sample_index.get(iid, -1) for iid in genotype_ids),
                               dtype=np.intp, count=self.n_samples)
            self.bgen_output_rows = np.flatnonzero(rows >= 0)
            self.bgen_indices = rows[self.bgen_output_rows]
            variants = {(sid, bp): i for i, (sid, bp) in enumerate(zip(self.bgen.rsids, self.bgen.positions))}
            self.bgen_markers = np.fromiter((variants.get((sid, bp), -1)
                                            for sid, bp in zip(bim["id"], bim["coordinate"])),
                                           dtype=np.intp, count=len(bim))
            valid = self.bgen_markers >= 0
            alleles = np.asarray(self.bgen.allele_ids)[self.bgen_markers[valid]]
            first = bim["allele1"].to_numpy()[valid]
            second = bim["allele2"].to_numpy()[valid]
            self.flipped = np.zeros(len(bim), dtype=np.uint8)
            reverse = np.fromiter((pair == f"{b},{a}" for pair, a, b in zip(alleles, first, second)),
                                  dtype=bool, count=len(alleles))
            self.flipped[valid] = reverse
            self.alleles_match = np.zeros(len(bim), dtype=bool)
            self.alleles_match[valid] = reverse | np.fromiter(
                (pair == f"{a},{b}" for pair, a, b in zip(alleles, first, second)), dtype=bool, count=len(alleles),
            )

    def read_calls(self, markers, *, include_phased=True):
        phased = None
        dosages = None
        if self.bed is not None:
            dosages = self.bed[self.sample_indices, markers].read(
                order="C", dtype=np.int8, _require_float32_64=False, num_threads=self.read_threads,
            ).val
            dosages[dosages == -127] = nan_integer
        if self.bgen is not None and (include_phased or dosages is None):
            from snipar.imputation.impute_from_sibs import decode_phased_genotypes

            if np.any(self.bgen_markers[markers] < 0):
                raise ValueError("The BED and BGEN datasets must contain the same SNP IDs and positions")
            if not self.alleles_match[markers].all():
                raise ValueError("The BED and BGEN datasets disagree on SNP alleles")
            if include_phased:
                phased = np.empty((self.n_samples, len(markers), 2), dtype=np.int8)
            from_phases = dosages is None
            if from_phases:
                dosages = np.empty((self.n_samples, len(markers)), dtype=np.int8)
            # Expand probabilities in small tiles, including when the phase
            # window requires many more markers than the retained output.
            tile = max(1, self.probability_bytes // (16 * max(1, len(self.bgen_indices))))
            for start in range(0, len(markers), tile):
                stop = min(len(markers), start + tile)
                selected = markers[start:stop]
                if len(self.bgen_indices):
                    options = {"num_threads": min(self.read_threads, len(selected))} if self.threaded_bgen else {}
                    probs = self.bgen.read((self.bgen_indices, self.bgen_markers[selected]), np.float32, **options)
                    calls = decode_phased_genotypes(probs, self.n_samples, self.bgen_output_rows,
                                                    self.flipped[selected])
                    # Release tile probabilities before the next decoder read.
                    del probs
                else:
                    calls = np.full((self.n_samples, len(selected), 2), nan_integer, dtype=np.int8)
                if phased is not None:
                    phased[:, start:stop, :] = calls
                if from_phases:
                    output = dosages[:, start:stop]
                    np.add(calls[:, :, 0], calls[:, :, 1], out=output)
                    output[output < 0] = nan_integer
        return phased, dosages

    def close(self):
        if self.bgen is not None:
            close = getattr(self.bgen, "close", None)
            if close is not None:
                close()
            self.bgen = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


def prepare_gts(phased_address, unphased_address, bim, pedigree_output, ped_ids, chromosomes,
                start=None, end=None, pcs=None, pc_ids=None, find_optimal_pc=None, *,
                marker_indices=None, prepared=None, retained=slice(None)):
    """Read hard genotypes and return per-marker frequencies and HDF5 metadata."""
    start = 0 if start is None else start
    end = len(bim) if end is None else end
    if not 0 <= start < end <= len(bim):
        raise ValueError("Imputation requires 0 <= start < end <= number of SNPs")
    if marker_indices is None:
        _, unique_indices = np.unique(bim["id"], return_index=True)
        unique_indices.sort()
        marker_indices = unique_indices[(unique_indices >= start) & (unique_indices < end)]
    if not len(marker_indices):
        raise ValueError("No unique SNPs in the requested slice")
    context = (ImputationGenotypes(phased_address, unphased_address, bim, ped_ids, pcs, pc_ids)
               if prepared is None else nullcontext(prepared))
    with context as genotypes:
        phased, dosages = genotypes.read_calls(marker_indices)
        target = dosages[:, retained]
        called = target != nan_integer
        with np.errstate(invalid="ignore", divide="ignore"):
            standard_f = np.sum(target, axis=0, where=called) / (2 * called.sum(axis=0))
        del called
        # All rows share the population frequency unless a PC model is fitted.
        frequencies = np.broadcast_to(standard_f, target.shape)
        selected_bim = bim.iloc[marker_indices[retained]]
        metadata = {
            "bim_columns": selected_bim.columns.to_numpy().astype("S"),
            "bim_values": selected_bim.to_numpy().astype("S"),
            "pedigree": pedigree_output,
            "non_duplicates": marker_indices[retained] - start,
            "standard_f": standard_f,
        }
        if genotypes.matched_pcs is not None:
            logging.info("estimating MAF using PCs ...")
            matched_pcs = genotypes.matched_pcs
            pc_gts = target if genotypes.pc_present.all() else target[genotypes.pc_present]
            if find_optimal_pc:
                matched_pcs = matched_pcs[:, :compute_aics(pc_gts, matched_pcs)]
            fitted, model_data = estimate_f(pc_gts, matched_pcs)
            if genotypes.pc_present.all():
                frequencies = fitted
            else:
                frequencies = frequencies.copy()
                frequencies[genotypes.pc_present] = fitted
            metadata.update((f"maf_{key}", val) for key, val in model_data.items())
        return (phased, dosages, genotypes.id_rows, bim["coordinate"].to_numpy(dtype=int)[marker_indices],
                frequencies, metadata)
