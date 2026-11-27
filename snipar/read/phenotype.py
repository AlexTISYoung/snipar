from pysnptools.snpreader import Pheno
import numpy as np
import pandas as pd
from snipar.gtarray import gtarray
from snipar.utilities import make_id_dict
def read_phenotype(file_path, column=None, column_index=None, na_values='NA'):
    """
    Read data from a text file with header structure where either:
    - First two columns are 'FID' and 'IID'
    - First column is 'IID'
    
    Parameters:
    file_path (str): Path to the text file
    column (str, optional): Name of column to extract (other than 'FID' or 'IID')
    column_index (int, optional): Index of column to extract (counting from 1 after 'IID'/'FID')
                                  Note: This is 1-based indexing
    na_values (str or list, optional): String or list of strings to recognize as NA/NaN. Default is 'NA'.
    
    Returns:
        y : :class:`snipar.gtarray`
            vector of non-missing phenotype values from specified column of phenofile along with individual IDs
    
    Note: If neither column nor column_index is provided, defaults to first column after IID/FID
    """
    with open(file_path, 'r') as file:
        header = next((line.split() for line in file if line.strip()), [])
    if header[:2] == ['FID', 'IID']:
        columns = header[2:]
    elif header[:1] == ['IID']:
        columns = header[1:]
    else:
        raise ValueError("Phenotype header must start with IID or FID IID")
    if column is None:
        column_index = 1 if column_index is None else column_index
        if not 1 <= column_index <= len(columns):
            raise ValueError(f"Column index {column_index} out of range")
        column = columns[column_index - 1]
    if column in ('FID', 'IID') or column not in columns:
        raise ValueError(f"Unknown phenotype column: {column}")
    if columns.count(column) != 1:
        raise ValueError(f"Duplicate phenotype column: {column}")
    print('Reading phenotype from column:', column)
    df = pd.read_csv(file_path, sep=r"\s+", usecols=['IID', column], dtype={'IID': str},
                     keep_default_na=False, na_values={column: na_values})
    values = pd.to_numeric(df[column], errors='coerce').to_numpy(dtype=float)
    keep = (df['IID'] != '').to_numpy() & np.isfinite(values)
    ids = df.loc[keep, 'IID'].to_numpy(dtype=str)
    if len(set(ids)) != len(ids):
        raise ValueError("Phenotype file contains duplicate individual IDs")
    return gtarray(values[keep, None], ids=ids)

def match_phenotype(G,y,pheno_ids):
    """Match a phenotype to a genotype array by individual IDs.

    Args:
        G : :class:`gtarray`
            genotype array to match phenotype to
        y : :class:`~numpy:numpy.array`
            vector of phenotype values
        pheno_ids: :class:`~numpy:numpy.array`
            vector of individual IDs corresponding to phenotype vector, y

    Returns:
       y : :class:`~numpy:numpy.array`
            vector of phenotype values matched by individual IDs to the genotype array

    """
    pheno_id_dict = make_id_dict(pheno_ids)
    return y[[pheno_id_dict[x] for x in G.ids]]

def read_covariates(covar, pheno_ids=None, missing_char = 'NA'):
    covar = Pheno(covar, missing=missing_char).read()
    X = np.array(covar.val)
    X = gtarray(X, ids=np.array(covar.iid)[:,1], sid=covar.sid)
    if pheno_ids is not None:
        in_covar = np.array([x in X.id_dict for x in pheno_ids])
        if np.sum((~in_covar))>0:
            raise(ValueError('Missing covariate values for some phenotyped individuals'))
    X.fill_NAs()
    return X
