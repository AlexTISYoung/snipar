'''
.. deprecated::
    Use :class:`.Pheno` instead.    
'''

import logging
import sys
import scipy as sp

def loadOnePhen(filename,  i_pheno = 0, missing = '-9', vectorize = False):
    '''
    .. deprecated::
       Use :class:`.Pheno` instead.    

    Load one column of a phenotype file. Remove any rows with missing data

    :param filename: name of the file
    :type filename: string
    :param i_pheno: column to return (default '0', the first column)
    :type i_pheno: int
    :param missing: value to threat as missing
    :type missing: string
    :param vectorize: if true, return a 1-D vector rather than a 2-D array
    :type vectorize: bool

    :rtype: An output dictionary

    The output dictionary looks like:

    * 'header' : [1] array phenotype namesv (only if header line is specified in file),
    * 'vals'   : [N*1] array of phenotype-data,
    * 'iid'    : [N*2] array of family IDs and case IDs
    '''
    missing = missing
    allColumns = loadPhen(filename, missing)
    i_present=allColumns['vals'][:,i_pheno]==allColumns['vals'][:,i_pheno]
    valsvector = allColumns['vals'][i_present,i_pheno]
    vals = sp.reshape(valsvector,(-1,1))
    iid = allColumns['iid'][i_present,:]
    #iid = iid.reshape(iid.shape[1], iid.shape[2])
    header = allColumns['header']
    if header is not None:
        header = [header[i_pheno]]

    if vectorize:
        vals = vals[:,0]

    ret = {
            'header':header,
            'vals':vals,
            'iid':iid
            }
    return ret


def loadPhen(filename, missing = '-9',famid='FID', sampid='ID'):
    '''
    .. deprecated::
       Use :class:`.Pheno` instead.    

    Load a phenotype or covariate file. Covariates have the same file format.

    :param filename: name of the file
    :type filename: string
    :param missing: value to threat as missing
    :type missing: string
    :param vectorize: if true, return a 1-D vector rather than a 2-D array
    :type vectorize: bool

    :rtype: An output dictionary

    The output dictionary looks like:

    * 'header' : [1] array phenotype names (only if header line is specified in file),
    * 'vals'   : [N*1] array of phenotype-data,
    * 'iid'    : [N*2] array of family IDs and case IDs
    '''
    if missing == '-9':
        logging.warning("loadPhen is using default missing value of '-9'.")

    data = sp.loadtxt(filename, dtype='str', comments=None)
    data = data.reshape(-1,data.shape[-1]) #Turns 1-d row into 2-d
    if data[0,0] == sampid: #One column of ids - use the single id as both the family id and the iid
        header = data[0,1::].tolist()
        iid = data[1:,[0,0]]
        valsStr = data[1:,1:]
    elif data[0,0] == famid:
        header = data[0,2::].tolist()
        iid = data[1:,0:2]
        valsStr = data[1:,2:]
    else:
        header = [None] * (data.shape[1]-2) # create a header containing a list of None's
        iid = data[:,0:2]
        valsStr = data[:,2:]

    
    if missing is not None:
        valsStr[valsStr==missing] = "NaN"
    vals = sp.array(valsStr,dtype = 'float')

    ret = {
            'header':header,
            'vals':vals,
            'iid':iid
            }
    return ret

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
