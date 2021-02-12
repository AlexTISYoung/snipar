import numpy as np
import sib_ldsc_z as ld
from scipy.optimize import minimize
from scipy.special import comb
from scipy.misc import derivative
import scipy.stats
from importlib import reload
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import datetime
import multiprocessing
# import numdifftools as nd
reload(ld)


import h5py
import pandas as pd

filenames = "/disk/genetics/ukb/alextisyoung/haplotypes/simulated_pops_large/from_chr1_to_chr23_start0_endNone_run0_p0-0_ab_corr0-5_vb0-25_length2/phenotype_dir_par_corr_0.5/1/chr_*.hdf5"
files = glob.glob(filenames)

file = files[0]
print("Reading in file: ", file)
hf = h5py.File(file, 'r')
metadata = hf.get("bim")[()]
chromosome = metadata[:, 0]
snp = metadata[:, 1]
bp = metadata[:, 3]
theta  = hf.get('estimate')[()]
se  = hf.get('estimate_ses')[()]
N = hf.get('N_L')[()]
S = hf.get('estimate_covariance')[()]
f = hf.get('freqs')[()]

# normalizing S
sigma2 = hf.get('sigma2')[()]
tau = hf.get('tau')[()]
phvar = sigma2+sigma2/tau

if len(files) > 1:
    for file in files[1:]:
        print("Reading in file: ", file)
        hf = h5py.File(file, 'r')
        metadata = hf.get("bim")[()]
        chromosome_file = metadata[:, 0]  
        snp_file = metadata[:, 1]
        bp_file = metadata[:, 3]
        theta_file  = hf.get('estimate')[()]
        se_file  = hf.get('estimate_ses')[()]
        S_file = hf.get('estimate_covariance')[()]
        f_file = hf.get('freqs')[()]
        N_file = hf.get('N_L')[()]

        # normalizing S
        sigma2 = hf.get('sigma2')[()]
        tau = hf.get('tau')[()]

        chromosome = np.append(chromosome, chromosome_file, axis = 0)
        snp = np.append(snp, snp_file, axis = 0)
        bp = np.append(bp, bp_file, axis = 0)
        theta = np.append(theta, theta_file, axis = 0)
        se = np.append(se, se_file, axis = 0)
        S = np.append(S, S_file, axis = 0)
        f = np.append(f, f_file, axis = 0)
        N = np.append(N, N_file, axis = 0)

# Constructing dataframe of data
zdata = pd.DataFrame({'CHR' : chromosome,
                    'SNP' : snp,
                    'BP' : bp,
                    'N' : N,
                    "f" : f,
                    'theta' : theta.tolist(),
                    'se' : se.tolist(),
                    "S" : S.tolist()})


# cleaning up a bit
zdata['CHR'] = zdata['CHR'].astype(int)
zdata['SNP'] = zdata['SNP'].astype(str).str.replace("b'", "").str[:-1]
zdata['BP'] = zdata['BP'].astype(str).str.replace("b'", "").str[:-1]
zdata['BP'] = zdata['BP'].astype('int')

# sorting by chromosome
zdata = zdata.sort_values(by = ['CHR']).reset_index(drop = True)
zdata['ordering'] = zdata.index


zdata_n_message = f"Number of Observations before merging LD-Scores, before removing low MAF SNPs: {zdata.shape[0]}"

print(zdata_n_message)

# dropping obs based on MAF
# zdata = zdata[zdata['f'] >= args.maf/100.0]

zdata_n_message = f"Number of Observations before merging LD-Scores, after removing low MAF SNPs: {zdata.shape[0]}"

print(zdata_n_message)

# == Reading in LD Scores == #
ldscore_path = "/disk/genetics/ukb/alextisyoung/haplotypes/simulated_pops_large/from_chr1_to_chr23_start0_endNone_run0_p0-0_ab_corr0-5_vb0-25_length2/ldscores/*[0-9].l2.ldscore.gz"
ldcolnames = ["CHR", "SNP", "BP", "L2"]
ldscores= ld.read_ldscores(ldscore_path, ldcolnames)
# ldscores['BP'] = ldscores['BP'].astype('int')

# Merging LD scores with main Data Frame
main_df = zdata.merge(ldscores, how = "inner", on = ["CHR", "SNP"])
main_df = main_df.sort_values(by = ['ordering'])

# dropping NAs
main_df = main_df.dropna()

S = np.array(list(main_df.S)) 
theta = np.array(list(main_df.theta))
f = np.array(list(main_df["f"]))
l = np.array(list(main_df["L2"]))
u = np.array(list(main_df["L2"]))

M = len(S)

S, theta = ld.transform_estimates("direct_plus_averageparental", S, theta)

zval = ld.theta2z(theta, S, M = M)

model = ld.sibreg(S = S, 
                z = zval, 
                l = l,
                f = f,
                u = u,
                M = M) 

output_matrix, result = model.solve(est_init = np.array([phvar/2, phvar/2, 0.0]),
                                    rbounds = True)


estimates = {'v1' : output_matrix['v1'],
            'v2' : output_matrix['v2'],
            'r' : output_matrix['r']}

std_errors = np.diag(output_matrix['std_err_mat'])

nblocks_blocksize = np.ceil(zval.shape[0]/200)
blocksize = int(nblocks_blocksize)


nblocks_blocksize = int(np.ceil(zval.shape[0]/200))
print("=================SNIPAR JKSE===============")
jkse = ld.jkse(model, output_matrix, 
        blocksize = nblocks_blocksize, num_procs=24)
print(f"{jkse}") #[0.01445773 0.01499649 0.03333096]
print("============================================")
print("\n")
print("\n")
print("\n")
print("\n")
print("============================================")

# def _check_shape(x, y):
#     '''Check that x and y have the correct shapes (for regression jackknives).'''
#     if len(x.shape) != 2 or len(y.shape) != 2:
#         raise ValueError('x and y must be 2D arrays.')
#     if x.shape[0] != y.shape[0]:
#         raise ValueError(
#             'Number of datapoints in x != number of datapoints in y.')
#     if y.shape[1] != 1:
#         raise ValueError('y must have shape (n_snp, 1)')
#     n, p = x.shape
#     if p > n:
#         raise ValueError('More dimensions than datapoints.')

#     return (n, p)

# class Jackknife():

#     '''
#     Base class for jackknife objects. Input involves x,y, so this base class is tailored
#     for statistics computed from independent and dependent variables (e.g., regressions).
#     The __delete_vals_to_pseudovalues__ and __jknife__ methods will still be useful for other
#     sorts of statistics, but the __init__ method will need to be overriden.
#     Parameters
#     ----------
#     x : np.matrix with shape (n, p)
#         Independent variable.
#     y : np.matrix with shape (n, 1)
#         Dependent variable.
#     n_blocks : int
#         Number of jackknife blocks
#     *args, **kwargs :
#         Arguments for inheriting jackknives.
#     Attributes
#     ----------
#     n_blocks : int
#         Number of jackknife blocks
#     p : int
#         Dimensionality of the independent varianble
#     N : int
#         Number of datapoints (equal to x.shape[0])
#     Methods
#     -------
#     jknife(pseudovalues):
#         Computes jackknife estimate and variance from the jackknife pseudovalues.
#     delete_vals_to_pseudovalues(delete_vals, est):
#         Converts delete values and the whole-data estimate to pseudovalues.
#     get_separators():
#         Returns (approximately) evenly-spaced jackknife block boundaries.
#     '''

#     def __init__(self, model, n_blocks=None, separators=None):
#         self.N, self.p = model.z.shape
#         if separators is not None:
#             if max(separators) != self.N:
#                 raise ValueError(
#                     'Max(separators) must be equal to number of data points.')
#             if min(separators) != 0:
#                 raise ValueError('Max(separators) must be equal to 0.')
#             self.separators = sorted(separators)
#             self.n_blocks = len(separators) - 1
#         elif n_blocks is not None:
#             self.n_blocks = n_blocks
#             self.separators = self.get_separators(self.N, self.n_blocks)
#         else:
#             raise ValueError('Must specify either n_blocks are separators.')

#         if self.n_blocks > self.N:
#             raise ValueError('More blocks than data points.')

#     @classmethod
#     def jknife(cls, pseudovalues):
#         '''
#         Converts pseudovalues to jackknife estimate and variance.
#         Parameters
#         ----------
#         pseudovalues : np.matrix pf floats with shape (n_blocks, p)
#         Returns
#         -------
#         jknife_est : np.matrix with shape (1, p)
#             Jackknifed estimate.
#         jknife_var : np.matrix with shape (1, p)
#             Variance of jackknifed estimate.
#         jknife_se : np.matrix with shape (1, p)
#             Standard error of jackknifed estimate, equal to sqrt(jknife_var).
#         jknife_cov : np.matrix with shape (p, p)
#             Covariance matrix of jackknifed estimate.
#         '''
#         n_blocks = pseudovalues.shape[0]
#         jknife_cov = np.atleast_2d(np.cov(pseudovalues.T, ddof=1) / n_blocks)
#         jknife_var = np.atleast_2d(np.diag(jknife_cov))
#         jknife_se = np.atleast_2d(np.sqrt(jknife_var))
#         jknife_est = np.atleast_2d(np.mean(pseudovalues, axis=0))
#         return (jknife_est, jknife_var, jknife_se, jknife_cov)

#     @classmethod
#     def delete_values_to_pseudovalues(cls, delete_values, est):
#         '''
#         Converts whole-data estimate and delete values to pseudovalues.
#         Parameters
#         ----------
#         delete_values : np.matrix with shape (n_blocks, p)
#             Delete values.
#         est : np.matrix with shape (1, p):
#             Whole-data estimate.
#         Returns
#         -------
#         pseudovalues : np.matrix with shape (n_blocks, p)
#             Psuedovalues.
#         Raises
#         ------
#         ValueError :
#             If est.shape != (1, delete_values.shape[1])
#         '''
#         n_blocks, p = delete_values.shape
#         if est.shape != (1, p):
#             raise ValueError(
#                 'Different number of parameters in delete_values than in est.')

#         return n_blocks * est - (n_blocks - 1) * delete_values

#     @classmethod
#     def get_separators(cls, N, n_blocks):
#         '''Define evenly-spaced block boundaries.'''
#         return np.floor(np.linspace(0, N, n_blocks + 1)).astype(int)
    
# class LstsqJackknifeSlow():

#     '''
#     Slow linear-regression block jackknife. This class computes delete values directly,
#     rather than forming delete values from block values. Useful for testing and for
#     non-negative least squares (which as far as I am aware does not admit a fast block
#     jackknife algorithm).
#     Inherits from Jackknife class.
#     Parameters
#     ----------
#     x : np.matrix with shape (n, p)
#         Independent variable.
#     y : np.matrix with shape (n, 1)
#         Dependent variable.
#     n_blocks : int
#         Number of jackknife blocks
#     nn: bool
#         Non-negative least-squares?
#     Attributes
#     ----------
#     est : np.matrix with shape (1, p)
#         FWLS estimate.
#     jknife_est : np.matrix with shape (1, p)
#         Jackknifed estimate.
#     jknife_var : np.matrix with shape (1, p)
#         Variance of jackknifed estimate.
#     jknife_se : np.matrix with shape (1, p)
#         Standard error of jackknifed estimate, equal to sqrt(jknife_var).
#     jknife_cov : np.matrix with shape (p, p)
#         Covariance matrix of jackknifed estimate.
#     delete_vals : np.matrix with shape (n_blocks, p)
#         Jackknife delete values.
#     Methods
#     -------
#     delete_values(x, y, func, s):
#         Compute delete values of func(x, y) the slow way, with blocks defined by s.
#     '''

#     def __init__(self, model, rbounds=True, n_blocks=None, separators=None):
        
#         mdl = Jackknife(model, n_blocks, separators)
        
#         self.rbounds = rbounds

#         func = model.solve
#         self.est = np.atleast_2d(np.array([model.output['v1'], model.output['v2'],
#                             model.output['r']]))
#         self.s = mdl.separators
#         self.delete_values = self.delete_values(model, func, mdl.separators)
#         self.pseudovalues = mdl.delete_values_to_pseudovalues(
#             self.delete_values, self.est)
#         (self.jknife_est, self.jknife_var, self.jknife_se, self.jknife_cov) =\
#             mdl.jknife(self.pseudovalues)
        
        
#     def delete_values(self, model, func, s):
#         '''
#         Compute delete values by deleting one block at a time.
#         Parameters
#         ----------
#         x : np.matrix with shape (n, p)
#             Independent variable.
#         y : np.matrix with shape (n, 1)
#             Dependent variable.
#         func : function (n, p) , (n, 1) --> (1, p)
#             Function of x and y to be jackknived.
#         s : list of ints
#             Block separators.
#         Returns
#         -------
#         delete_values : np.matrix with shape (n_blocks, p)
#             Delete block values (with n_blocks blocks defined by parameter s).
#         Raises
#         ------
#         ValueError :
#             If x.shape[0] does not equal y.shape[0] or x and y are not 2D.
#         '''
        
#         d = []
#         print(len(s))
#         for i in range(len(s) - 1):
#             print(f"Block Number: {i}")
#             d_in = func(z = np.vstack([model.z[0:s[i], ...], model.z[s[i + 1]:, ...]]), 
#                       S = np.vstack([model.S[0:s[i], ...], model.S[s[i + 1]:, ...]]),
#                       l = np.hstack([model.l[0:s[i], ...], model.l[s[i + 1]:, ...]]),
#                       u = np.hstack([model.u[0:s[i], ...], model.u[s[i + 1]:, ...]]),
#                       f = model.f,
#                       M = model.M,
#                     printout = False,
#                     est_init = self.est,
#                     rbounds = self.rbounds)[0]
#             dmat = np.array([d_in['v1'], d_in['v2'], d_in['r']])
            
#             d.append(dmat)

#         return np.array(d)
    
    
# ldsc_jk = LstsqJackknifeSlow(model, n_blocks=200)
# print("==============LDSC JKSE========================")
# print(ldsc_jk.jknife_se)
# # [[0.69163511 0.71471265 0.95458665]]
# print("============================================")
# actual ldsc SE for r: 0.1277