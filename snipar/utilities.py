from typing import Any, Dict, Optional, Union
from pathlib import Path
import numpy as np
import pandas as pd
import bgen_reader
from os import path
import argparse
import re
from typing import Tuple

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

def coord2linear(ind1: int, ind2: int) -> int:
    row_ind, col_ind = max(ind1, ind2), min(ind1, ind2)
    return int(row_ind * (row_ind + 1) / 2 + col_ind)


def linear2coord(ind: int) -> Tuple[int, int]:
    ind += 1
    ind1 = np.ceil((2 * ind + 0.25) ** 0.5 - 0.5)
    ind2 = ind - (ind1 - 1) * ind1 / 2
    return int(ind1 - 1), int(ind2 - 1)

def get_parser_doc(parser):
    doc = ""
    for action in parser._actions:
        options = str(action.option_strings)[1:-1]
        default = action.default
        type = ""
        if action.type:
            if "'" in str(action.type):
                type = str(action.type).split("'")[1]
        
        help = action.help
        default_substring = ""
        if default:
            default_substring = f", default={default}"
        
        arg_doc = f"""    {options} : {type}{default_substring}
            {help}

"""
        doc += arg_doc
    
    return doc


class _bgen:
    """Simple wrapper class of bgen_reader.open_bgen that controls access to sample ids."""
    def __init__(self, 
                 filepath: Union[str, Path],
                 samples_filepath: Optional[Union[str, Path]] = None,
                 allow_complex: bool = False,
                 verbose: bool = True,):
        self._bgen = bgen_reader.open_bgen(filepath, allow_complex=allow_complex, verbose=verbose)
        self._samples = self._read_sample(samples_filepath)
        if self._bgen.samples.shape[0] != self._samples.shape[0]:
            raise ValueError('sample file length and bgen file length do not match.')
    
    @property
    def samples(self) -> np.ndarray:
        return self._samples
    
    def __getattr__(self, name: str) -> Any:
        """Access to attributes other than sample ids."""
        return getattr(self._bgen, name)
    
    def _read_sample(self, samples_filepath: str) -> np.ndarray:
        return pd.read_csv(samples_filepath, sep='\s+')['ID_2'][1:].to_numpy(dtype=str)


def open_bgen(filename: str, verbose: bool = False) -> bgen_reader.open_bgen:
    """Wrapper of bgen_reader.open_bgen that checks if sample ids make sense."""
    bgen = bgen_reader.open_bgen(filename, verbose=verbose)
    if bgen.samples[0] == 'sample_0':
        print('WARNING: Sample ids in bgen file are generic. Trying to read the corresponding .sample file ...')
        samples_filepath = filename[:-4] + 'sample'
        if not path.exists(samples_filepath):
            raise FileNotFoundError(f'{samples_filepath} does not exist.')
        bgen = _bgen(filename, verbose=verbose, samples_filepath=samples_filepath)
    return bgen


def read_bgen(filename: str, verbose: bool = False) -> Dict:
    """Wrapper of bgen_reader.read_bgen that checks if sample ids make sense."""
    bgen = bgen_reader.read_bgen(filename, verbose=verbose)
    if bgen['samples'][0] == 'sample_0':
        print('WARNING: Sample ids in bgen file are generic. Trying to read the corresponding .sample file ...')
        samples_filepath = filename[:-4] + 'sample'
        if not path.exists(samples_filepath):
            raise FileNotFoundError(f'{samples_filepath} does not exist.')
        bgen = bgen_reader.read_bgen(filename, verbose=verbose, samples_filepath=samples_filepath)
        bgen['samples'] = pd.read_csv(samples_filepath, sep='\s+')['ID_2'][1:].to_numpy()
    return bgen
        