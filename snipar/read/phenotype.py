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
    # Determine delimiter (tab or whitespace)
    with open(file_path, 'r') as file:
        first_line = file.readline()
        delimiter = '\t' if '\t' in first_line else ' '  
        header = first_line.split(delimiter)
        header[-1] = header[-1].strip()  # Remove newline character
    # Determine file format based on header
    has_fid = (len(header) > 1 and header[0] == 'FID' and header[1] == 'IID')
    # Set default column if neither is provided
    if column is None and column_index is None:
        # Default to first column after IID/FID
        column_index = 1
    # Determine the usecols parameter for pd.read_csv
    if column is not None:
        if column in ['FID', 'IID']:
            raise ValueError(f"Phenotype cannot be named FID or IID")
        # We need to read the IID column and the target column
        cols_to_use = ['IID', column]
    else:  # column_index is provided
        # Adjust column_index based on file format
        offset = 2 if has_fid else 1
        adjusted_index = column_index + offset - 1  # -1 for 0-based indexing
        if adjusted_index >= len(header):
            raise ValueError(f"Column index {column_index} out of range")
        column = header[adjusted_index]
        cols_to_use = ['IID', column]
    print('Reading phenotype from column:', column)
    # Read the data using pandas for efficiency, handling missing values
    df = pd.read_csv(file_path, 
                     sep=delimiter,
                     usecols=cols_to_use,
                     na_values=na_values)
    # Verify target column contains numeric data
    try:
        df[column] = pd.to_numeric(df[column], errors='coerce')
    except ValueError:
        raise ValueError(f"Phenotype contains non-numeric values that cannot be converted")
    # Remove rows with missing values in either IID or target column
    df = df.dropna(subset=['IID', column])
    # Return gtarray
    return gtarray(np.array(df[column].values).reshape((df.shape[0],1)), ids=np.array(df['IID'].values, dtype=str))

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
    in_G_dict = np.array([x in G.id_dict for x in pheno_ids])
    y = y[in_G_dict]
    pheno_ids = pheno_ids[in_G_dict]
    pheno_id_dict = make_id_dict(pheno_ids)
    y = y[[pheno_id_dict[x] for x in G.ids]]
    return y

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