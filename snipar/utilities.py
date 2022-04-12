import numpy as np
from os import path
import argparse
import re
def make_id_dict(x,col=0):
    """
    Make a dictionary that maps from the values in the given column (col) to their row-index in the input array
    """
    if len(x.shape)>1:
        x = x[:,col]
    id_dict = {}
    for i in range(0,x.shape[0]):
        id_dict[x[i]] = i
    return id_dict

def convert_str_array(x):
    """
    Convert an ascii array to unicode array (UTF-8)
    """
    x = np.array(x)
    x_shape = x.shape
    x = x.flatten()
    x_out = np.array([y.decode('UTF-8') for y in x])
    return x_out.reshape(x_shape)

def encode_str_array(x):
    """
    Encode a unicode array as an ascii array
    """
    x = np.array(x)
    x_shape = x.shape
    x = x.flatten()
    x_out = np.array([y.encode('ascii') for y in x])
    return x_out.reshape(x_shape)

def parse_obsfiles(obsfiles, obsformat='bed', append = True, wildcard = '@', chromosomes=None):
    obs_files = []
    chroms = []
    if wildcard in obsfiles:
        bed_ixes = obsfiles.split(wildcard)
        if chromosomes is None:
            chromosomes = np.arange(1,23)
        for i in chromosomes:
            obsfile = bed_ixes[0]+str(i)+bed_ixes[1]+'.'+obsformat
            if path.exists(obsfile):
                if append:
                    obs_files.append(obsfile)
                else:
                    obs_files.append(obsfile[:-(len(obsformat)+1)])
                chroms.append(i)
        if len(obs_files)==0:
            raise(ValueError('Observed genotype files not found'))
        else:
            print(str(len(obs_files))+' files found')
    else:
            obsfile = obsfiles+'.'+obsformat
            if path.exists(obsfile):
                obs_files = [obsfile]
                chroms = [0]
            else:
                raise(ValueError(obsfile+' not found'))
    return np.array(obs_files), np.array(chroms,dtype=int)

def parse_filelist(obsfiles, impfiles, obsformat, chromosomes=None):
    obs_files = []
    imp_files = []
    chroms = []
    if '@' in obsfiles and impfiles:
        bed_ixes = obsfiles.split('@')
        imp_ixes = impfiles.split('@')
        if chromosomes is None:
            chromosomes = np.arange(1,23)
        for i in chromosomes:
            obsfile = bed_ixes[0]+str(i)+bed_ixes[1]+'.'+obsformat
            impfile = imp_ixes[0]+str(i)+imp_ixes[1]+'.hdf5'
            if path.exists(impfile) and path.exists(obsfile):
                obs_files.append(obsfile)
                imp_files.append(impfile)
                chroms.append(i)
        if len(imp_files)==0:
            raise(ValueError('Observed/imputed genotype files not found'))
        else:
            print(str(len(imp_files))+' matched observed and imputed genotype files found')
    else:
            obsfile = obsfiles+'.'+obsformat
            impfile = impfiles+'.hdf5'
            if path.exists(obsfile) and path.exists(impfile):
                obs_files = [obsfile]
                imp_files = [impfile]
                chroms = [0]
            else:
                if not path.exists(obsfile):
                    raise(ValueError(obsfile+' not found'))
                if not path.exists(impfile):
                    raise(ValueError(impfile+' not found'))
    return np.array(obs_files), np.array(imp_files), np.array(chroms,dtype=int)

def outfile_name(outprefix,outsuffix,chrom=None):
    if '@' in outprefix:
        if chrom is None:
            raise(ValueError('Must provide chromosome number with wildcard character'))
        outprefix = outprefix.split('@')
        outprefix = outprefix[0]+str(chrom)+outprefix[1]
        return outprefix+outsuffix
    elif chrom is not None:
        return outprefix+'chr_'+str(chrom)+outsuffix
    else:
        return outprefix+outsuffix

def parseNumRange(string):
    """reads either a int or a range"""
    match_range = re.match(r' *(\d+) *- *(\d+) *', string)
    match_list = re.match(r' *(\d+) *', string)
    if match_range:
        start = int(match_range.group(1))
        end = int(match_range.group(2))
        result = [str(n) for n in range(start, end+1)]
    elif match_list:
        result = match_list.group(0)
    else:
        raise Exception(f"{string} is neither a range of the form x-y nor a list of integers of the form x y z")
    return result

class NumRangeAction(argparse.Action):
    """flattens and sorts the resulting num range. also removes duplicates"""
    def __call__(self, parser, args, values, option_string=None):
        result = []
        for v in values:
            if isinstance(v,list):
                result += v
            if isinstance(v,str):
                result.append(v)
        result = np.unique(result).astype(int)
        result.sort()
        result = result.astype(str).tolist()
        setattr(args, self.dest, result)