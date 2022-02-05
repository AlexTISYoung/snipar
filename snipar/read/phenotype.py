from pysnptools.snpreader import Pheno
import numpy as np
from snipar.gtarray import gtarray
from snipar.utilities import make_id_dict

def read_phenotype(phenofile, missing_char = 'NA', phen_index = 1):
    """Read a phenotype file and remove missing values.

    Args:
        phenofile : :class:`str`
            path to plain text phenotype file with columns FID, IID, phenotype1, phenotype2, ...
        missing_char : :class:`str`
            The character that denotes a missing phenotype value; 'NA' by default.
        phen_index : :class:`int`
           The index of the phenotype (counting from 1) if multiple phenotype columns present in phenofile

    Returns:
        y : :class:`~numpy:numpy.array`
            vector of non-missing phenotype values from specified column of phenofile
        pheno_ids: :class:`~numpy:numpy.array`
            corresponding vector of individual IDs (IID)
    """
    pheno = Pheno(phenofile, missing=missing_char)[:,phen_index-1].read()
    y = np.array(pheno.val)
    y.reshape((y.shape[0],1))
    pheno_ids = np.array(pheno.iid)[:,1]
    # Remove y NAs
    y_not_nan = np.logical_not(np.isnan(y[:,0]))
    if np.sum(y_not_nan) < y.shape[0]:
        y = y[y_not_nan,:]
        pheno_ids = pheno_ids[y_not_nan]
    print('Number of non-missing phenotype observations: ' + str(y.shape[0]))
    return gtarray(y,ids=pheno_ids)

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
    X = gtarray(X, ids=np.array(covar.iid)[:,1])
    if pheno_ids is not None:
        in_covar = np.array([x in X.id_dict for x in pheno_ids])
        if np.sum((~in_covar))>0:
            raise(ValueError('Missing covariate values for some phenotyped individuals'))
    X.fill_NAs()
    return X