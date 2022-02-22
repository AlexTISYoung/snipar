import numpy as np
from os import path

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

def parse_obsfiles(obsfiles, obsformat='bed'):
    obs_files = []
    chroms = []
    if '~' in obsfiles:
        bed_ixes = obsfiles.split('~')
        for i in range(1,23):
            obsfile = bed_ixes[0]+str(i)+bed_ixes[1]+'.'+obsformat
            if path.exists(obsfile):
                obs_files.append(obsfile)
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

def parse_filelist(obsfiles, impfiles, obsformat):
    obs_files = []
    imp_files = []
    chroms = []
    if '~' in obsfiles and impfiles:
        bed_ixes = obsfiles.split('~')
        imp_ixes = impfiles.split('~')
        for i in range(1,23):
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