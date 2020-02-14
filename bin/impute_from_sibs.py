#!/well/kong/users/wiw765/anaconda2/bin/python
import numpy as np
import pandas as pd
import numpy.ma as ma
import argparse#, h5py
import logging
from pysnptools.snpreader import Bed
import json
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s')
parser = argparse.ArgumentParser()
parser.add_argument('chr',type=int,help='Which chromosome (integer)')
parser.add_argument('ibd',type=str,help='IBD file')
parser.add_argument('genotypes',type=str,help='Genotypes in .bed format')
parser.add_argument('ped',type=str,help='Pedigree file with siblings sharing family ID')
parser.add_argument('out',type=str,help='Prefix of hdf5 output of imputed parental genotypes')
parser.add_argument('--king',action='store_true',default=False,help='IBD segments file in KING format (default 23andme)')
parser.add_argument('--bim',type=str,default = None, help='Bim file giving positions of SNPs in KING IBD file if different from Bim file of genotypes')
parser.add_argument('--start', type=int,
                    help='Start index of SNPs to perform imputation for in genotype file (starting at zero)',
                    default=0)
parser.add_argument('--end', type=int,
                    help='End index of SNPs to perform imputation for in genotype file (goes from 0 to (end-1)',
                    default=None)
args=parser.parse_args()
#TODO use scipy interface
def impute_snp(snp, snp_ibd0, snp_ibd1, snp_ibd2, f, bed):
    #TODO check division by two
    #check for NANs
    result = 0
    if len(snp_ibd0) > 0:        
        for i,j in snp_ibd0:
            result += (bed[i,snp]+bed[j,snp])
        result = result/len(snp_ibd0)

    elif len(snp_ibd1) > 0:
        for i,j in snp_ibd1:
            sibsum = (bed[i,snp]+bed[j,snp])
            additive = 0
            if sibsum==0:
                additive = f
            elif sibsum==1:
                additive = 1+f
            elif sibsum==2:
                additive = 1+2*f
            elif sibsum==3:
                additive = 2+f
            elif sibsum==4:
                additive = 3+f
            result += additive
        result = result/len(snp_ibd1)
    
    elif len(snp_ibd2) > 0:
        result = 0
        for i,j in snp_ibd2:
            result += (bed[i,snp]+bed[j,snp]) + 2*f
        result = result/len(snp_ibd2)
    if result is None:
        print("AAAAAAAAAAA")
        print(snp, snp_ibd0, snp_ibd1, snp_ibd2, f)
    return result
        

def get_IBD_type(id1, id2, loc, ibd_dict):
    #todo use get
    segments = []
    if (id1, id2) in ibd_dict:
        segments = ibd_dict[(id1, id2)]
    elif (id2, id1) in ibd_dict:
        segments = ibd_dict[(id2, id1)]
    for seg in segments:
        if seg[0] <= loc <= seg[1]:
            if seg[2] == "IBD2":
                return 2
            else:
                return 1
    return 0

def prepare_data(args):
    logging.info("initializing data")
    logging.info("loading and filtering pedigree file ...")
    #keeping individuals with no parents
    ped = pd.read_csv(args.ped, sep = " ")
    ped["has_father"] = ~ ped["FATHER_ID"].str.endswith("_P")
    ped["has_mother"] = ~ ped["MOTHER_ID"].str.endswith("_M")
    ped = ped[~(ped["has_mother"] | ped["has_father"])]
    #TODO handle sibs with one parent
    ped_ids =  set(ped["IID"].tolist())
    #finding siblings in each family
    sibships = ped.groupby(["FID", "FATHER_ID", "MOTHER_ID"]).agg({'IID':lambda x: list(x)}).reset_index()
    sibships = sibships[sibships["IID"].apply(len)>1]
    logging.info("loading bim file ...")
    if args.bim is None:
        bim_file = args.genotypes+'.bim'
    else:
        bim_file = args.bim
    bim = pd.read_csv(bim_file, sep = "\t", header=None, names=["Chr", "id", "morgans", "coordinate", "allele1", "allele2"])
    #tODO what if people are related but do not have ibd on chrom
    logging.info("loading and transforming ibd file ...")
    ibd = pd.read_csv(args.ibd, sep = "\t").astype(str)
    #Adding location of start and end of each 
    ibd = ibd[ibd["Chr"] == str(args.chr)][["ID1", "ID2", "IBDType", "StartSNP", "StopSNP"]]
    temp = bim[["id", "coordinate"]].rename(columns = {"id":"StartSNP","coordinate":"StartSNPLoc"})
    ibd= ibd.merge(temp, on="StartSNP")
    temp = bim[["id", "coordinate"]].rename(columns = {"id":"StopSNP","coordinate":"StopSNPLoc"})
    ibd = ibd.merge(temp, on="StopSNP")
    ibd['segment'] = ibd[['StartSNPLoc', 'StopSNPLoc', "IBDType"]].apply(tuple, axis=1)
    ibd = ibd.groupby(["ID1", "ID2"]).agg({'segment':lambda x: list(x)}).to_dict()["segment"]
    #TODO orders
    logging.info("loading bed file ...")
    gts_f = Bed(args.genotypes+".bed")
    if args.end is not None:
        ids_in_ped = [(id in ped_ids) for id in gts_f.iid[:,1]]
        gts_ids = gts_f.iid[ids_in_ped]
        gts = gts_f[ids_in_ped , args.start:args.end].read().val
        pos = gts_f.pos[args.start:args.end, 2]
        sid = gts_f.sid[args.start:args.end]
    else:
        gts = gts_f[ [(id in ped_ids) for id in gts_f.iid[:,1]], :].read().val
        pos = gts_f.pos[:, 2]
        sid = gts_f.sid
    iid_to_bed_index = {i:index for index, i in enumerate(gts_ids[:,1])}
    logging.info("initializing data done ...")
    return sibships, iid_to_bed_index, gts, ibd, pos

def impute(sibships, iid_to_bed_index, gts, ibd, pos):
    logging.info("imputing data ...")
    freqs = np.nanmean(gts,axis=0)/2.0
    imputed_par_gts = [[0 for i in range(gts.shape[1])] for j in range(sibships.shape[0]+1)]
    for index in range(sibships.shape[0]):
        row = sibships.iloc[index]
        sibs_index = [iid_to_bed_index[iid] for iid in row["IID"]]   
        for snp in range(gts.shape[1]):
            #TODO tranlate address
            loc = pos[snp]
            snp_ibd0 = []
            snp_ibd1 = []
            snp_ibd2 = []
            for i in range(1,len(sibs_index)):
                for j in range(i):
                    sib1_index = sibs_index[i]
                    sib2_index = sibs_index[j]
                    sib1_id = row["IID"][i]
                    sib2_id = row["IID"][j]
                    ibd_type = get_IBD_type(sib1_id, sib2_id, loc, ibd)
                    if ibd_type == 2:
                        snp_ibd2.append((sib1_index, sib2_index))
                    if ibd_type == 1:
                        snp_ibd1.append((sib1_index, sib2_index))
                    if ibd_type == 0:
                        snp_ibd0.append((sib1_index, sib2_index))
            imputed_par_gts[index][snp] = impute_snp(snp, snp_ibd0, snp_ibd1, snp_ibd2, freqs[snp], gts)
    return imputed_par_gts

#TODO add quality control steps
sibships, iid_to_bed_index, gts, ibd, pos = prepare_data(args)
imputed_par_gts = impute(sibships, iid_to_bed_index, gts, ibd, pos)
with open(args.out, "w") as f:
    json.dump(imputed_par_gts , f)