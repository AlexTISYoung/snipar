import numpy.ma as ma
from numba import njit, prange
import gzip
from snipar.ld import ldscores_from_bed
from snipar.utilities import make_id_dict
import numpy as np

class sumstats(object):
    def __init__(self,chrom,sid,pos, A1, A2, freqs, direct, direct_SE, avg_NTC, avg_NTC_SE, population, population_SE, r_direct_avg_NTC, r_direct_pop, ldscores = None, map=None):
        sizes = np.array([sid.shape[0],pos.shape[0],A1.shape[0],A2.shape[0],freqs.shape[0],direct.shape[0],
                            avg_NTC.shape[0],population.shape[0],r_direct_avg_NTC.shape[0],r_direct_pop.shape[0]])
        if np.unique(sizes).shape[0] > 1:
            raise(ValueError('All inputs to sumstats class must have same size'))
        self.chrom = np.zeros(sid.shape,dtype=int)
        self.chrom[:] = int(chrom)
        self.sid = np.array(sid,dtype=str)
        self.sid_dict = make_id_dict(self.sid)
        self.pos = np.array(pos,dtype=int)
        self.A1 = np.array(A1,dtype=str)
        self.A2 = np.array(A2,dtype=str)
        self.freqs = ma.array(freqs,dtype=float)
        self.freqs.mask = np.isnan(self.freqs)
        self.direct = ma.array(direct, dtype=float)
        self.direct.mask = np.isnan(self.direct)
        self.direct_SE = ma.array(direct_SE, dtype=float)
        self.direct_SE.mask = np.isnan(self.direct_SE)
        self.avg_NTC = ma.array(avg_NTC, dtype=float)
        self.avg_NTC.mask = np.isnan(self.avg_NTC)
        self.avg_NTC_SE = ma.array(avg_NTC_SE, dtype=float)
        self.avg_NTC_SE.mask = np.isnan(self.avg_NTC_SE)
        self.population = ma.array(population, dtype=float)
        self.population.mask = np.isnan(self.population)
        self.population_SE = ma.array(population_SE, dtype=float)
        self.population_SE.mask = np.isnan(self.population_SE)
        self.r_direct_avg_NTC = ma.array(r_direct_avg_NTC, dtype=float)
        self.r_direct_avg_NTC.mask = np.isnan(self.r_direct_avg_NTC)
        self.r_direct_pop = ma.array(r_direct_pop, dtype=float)
        self.r_direct_pop.mask = np.isnan(self.r_direct_pop)
        if ldscores is not None:
            if not ldscores.shape[0] == sid.shape[0]:
                raise(ValueError('LD scores must have same size as other sumstats'))
            self.ldscores = ma.array(ldscores,dtype=float)
            self.ldscores.mask = np.isnan(self.ldscores)
        else:
            self.ldscores = None
        if map is not None:
            if not map.shape[0] == sid.shape[0]:
                raise(ValueError('LD scores must have same size as other sumstats'))
            self.map = ma.array(map,dtype=float)
            self.map.mask = np.isnan(self.map)
        else:
            self.map = None
    def concatenate(self,s2):
        self.chrom = np.hstack((self.chrom,s2.chrom))
        self.sid = np.hstack((self.sid, s2.sid))
        self.sid_dict = make_id_dict(self.sid)
        self.pos = np.hstack((self.pos, s2.pos))
        self.A1 = np.hstack((self.A1, s2.A1))
        self.A2 = np.hstack((self.A2, s2.A2))
        self.freqs = ma.concatenate([self.freqs, s2.freqs])
        self.direct = ma.concatenate([self.direct, s2.direct])
        self.direct_SE = ma.concatenate([self.direct_SE, s2.direct_SE])
        self.avg_NTC = ma.concatenate([self.avg_NTC, s2.avg_NTC])
        self.avg_NTC_SE = ma.concatenate([self.avg_NTC_SE, s2.avg_NTC_SE])
        self.population = ma.concatenate([self.population, s2.population])
        self.population_SE = ma.concatenate([self.population_SE, s2.population_SE])
        self.r_direct_avg_NTC = ma.concatenate([self.r_direct_avg_NTC, s2.r_direct_avg_NTC])
        self.r_direct_pop = ma.concatenate([self.r_direct_pop, s2.r_direct_pop])
        if self.ldscores is not None and s2.ldscores is not None:
            self.ldscores = ma.concatenate([self.ldscores, s2.ldscores])
        if self.map is not None and s2.map is not None:
            self.map = ma.concatenate([self.map, s2.map])
    def filter(self,filter_pass):
        self.chrom = self.chrom[filter_pass]
        self.sid = self.sid[filter_pass]
        self.sid_dict = make_id_dict(self.sid)
        self.pos = self.pos[filter_pass]
        self.A1 = self.A1[filter_pass]
        self.A2 = self.A2[filter_pass]
        self.freqs = self.freqs[filter_pass]
        self.direct = self.direct[filter_pass]
        self.direct_SE = self.direct_SE[filter_pass]
        self.avg_NTC = self.avg_NTC[filter_pass]
        self.avg_NTC_SE = self.avg_NTC_SE[filter_pass]
        self.population = self.population[filter_pass]
        self.population_SE = self.population_SE[filter_pass]
        self.r_direct_avg_NTC = self.r_direct_avg_NTC[filter_pass]
        self.r_direct_pop = self.r_direct_pop[filter_pass]
        if self.ldscores is not None:
            self.ldscores = self.ldscores[filter_pass]
        if self.map is not None:
            self.map = self.map[filter_pass]
    def filter_NAs(self):
        no_NAs = (~self.freqs.mask)*(~self.direct.mask)*(~self.avg_NTC.mask)*(~self.population.mask)*(~self.r_direct_pop.mask)*(~self.r_direct_avg_NTC.mask)
        if self.ldscores is not None:
            no_NAs = no_NAs*(~self.ldscores.mask)
        if self.map is not None:
            no_NAs = no_NAs*(~self.map.mask)
        self.filter(no_NAs)
    def filter_maf(self,min_maf):
        if min_maf<0 or min_maf>1:
            raise(ValueError('MAF must be between 0 and 1'))
        good_maf = np.logical_and(min_maf < self.freqs, self.freqs < (1-min_maf))
        self.filter(good_maf)
    def filter_corrs(self,max_Z):
        r_dir_avg_NTC_Z = (self.r_direct_avg_NTC-np.mean(self.r_direct_avg_NTC))/np.std(self.r_direct_avg_NTC)
        r_dir_pop_Z = (self.r_direct_pop-np.mean(self.r_direct_pop))/np.std(self.r_direct_pop)
        good_corrs = np.logical_and(np.abs(r_dir_avg_NTC_Z) < max_Z, np.abs(r_dir_pop_Z) < max_Z)
        self.filter(good_corrs)
    def scores_from_ldsc(self,ld_files):
        self.ldscores = ma.array(np.zeros(self.sid.shape[0]), mask=np.ones(self.sid.shape[0]))
        for i in range(ld_files.shape[0]):
            print('Reading LD scores from '+ld_files[i])
            ld_chr = np.loadtxt(ld_files[i],dtype=str,usecols=(1,3))
            in_sid_dict = np.array([x in self.sid_dict for x in ld_chr[:,0]])
            if np.sum(in_sid_dict) > 0:
                ld_chr = ld_chr[in_sid_dict,:]
                ld_indices = np.array([self.sid_dict[x] for x in ld_chr[:,0]])
                self.ldscores[ld_indices] = np.array(ld_chr[:,1],dtype=float)
                self.ldscores.mask[ld_indices] = False
            else:
                raise(ValueError('No SNPs in common between LD scores and summary statistics'))
    def cor_direct_pop(self, n_blocks):
        print('Computing correlation between direct and population effects')
        r_dir_pop, r_dir_pop_SE, r_dir_pop_delete = jacknife_est(self.direct,self.population,np.power(self.direct_SE,2),
                                                                np.power(self.population_SE,2),self.r_direct_pop,self.ldscores,n_blocks)
        return r_dir_pop, r_dir_pop_SE, r_dir_pop_delete
    def cor_direct_avg_NTC(self, n_blocks):
        print('Computing correlation between direct effects and average NTCs')
        r_dir_avg_NTC, r_dir_avg_NTC_SE, r_dir_avg_NTC_delete = jacknife_est(self.direct,self.avg_NTC,np.power(self.direct_SE,2),
                                                                                np.power(self.avg_NTC_SE,2),self.r_direct_avg_NTC,self.ldscores, n_blocks)
        return r_dir_avg_NTC, r_dir_avg_NTC_SE, r_dir_avg_NTC_delete
    def compute_ld_scores(self, bedfiles, chroms, ld_wind, ld_out=None):
        self.ldscores = ma.array(np.zeros(self.sid.shape[0]), mask=np.ones(self.sid.shape[0]))
        for i in range(bedfiles.shape[0]):
            print('Computing LD scores for chromosome '+str(chroms[i]))
            ld_chr, ld_snps_chr = ldscores_from_bed(bedfiles[i], chroms[i], ld_wind, ld_out)
            in_sid_dict = np.array([x in self.sid_dict for x in ld_snps_chr])
            if np.sum(in_sid_dict) > 0:
                ld_indices = np.array([self.sid_dict[x] for x in ld_snps_chr[in_sid_dict]])
                self.ldscores[ld_indices] = np.array(ld_chr[in_sid_dict],dtype=float)
                self.ldscores.mask[ld_indices] = False
            else:
                raise(ValueError('No SNPs in common between LD scores and sumstats'))

def read_sumstats_file(sumstats_file, chrom):
    print('Reading sumstats from '+sumstats_file)
    # Read header
    sumstats_file_o = gzip.open(sumstats_file,'r')
    sumstats_header = sumstats_file_o.readline().decode('UTF-8').split(' ')
    sumstats_header[len(sumstats_header)-1] = sumstats_header[len(sumstats_header)-1].split('\n')[0]
    sumstats_header = np.array(sumstats_header)
    sumstats_file_o.close()
    ## Find columns
    # Keep backwards compatibility with old sumstats files
    if 'avg_parental_Beta' in sumstats_header:
        parname = 'avg_parental'
    else:
        parname = 'avg_NTC'
    # Find column indices
    colnames = ['SNP','pos','A1','A2','freq','direct_Beta',
                'direct_SE',parname+'_Beta',parname+'_SE',
                'population_Beta','population_SE','r_direct_'+parname,'r_direct_population']
    col_indices = tuple([np.where(sumstats_header==x)[0][0] for x in colnames])
    # Read summary statistics
    s = np.loadtxt(sumstats_file,dtype=str,skiprows=1,usecols=col_indices)
    # Return sumstats class
    return sumstats(chrom, s[:,0],s[:,1],s[:,2],s[:,3],s[:,4],
                    s[:,5],s[:,6],s[:,7],s[:,8],s[:,9],s[:,10],s[:,11],s[:,12])

def read_sumstats_files(sumstats_files, chroms):
    s = read_sumstats_file(sumstats_files[0], chroms[0])
    if chroms.shape[0]>1:
        for i in range(1,sumstats_files.shape[0]):
            s.concatenate(read_sumstats_file(sumstats_files[i], chroms[i]))
    return s

@njit
def estimate_var(z,v,w):
    return np.sum(w*(np.power(z,2)-v))

@njit
def estimate_weights(v_est,v,l):
    w = np.power((v_est+v)*l,-1)
    w = w/np.sum(w)
    return w

@njit
def reweighted_estimate_var(z, v, l, tol=1e-8, max_iter=10**3):
    # Initialize weights
    w = np.power(v*l,-1)
    w = w/np.sum(w)
    # Initialize estimate
    v_est = estimate_var(z, v, w)
    # Iterate to convergence
    v_diff = np.inf
    n_iter = 0
    while v_diff > tol and n_iter < max_iter:
        w = estimate_weights(v_est, v, l)
        v_est_new = estimate_var(z, v, w)
        v_diff = np.abs(v_est_new-v_est)
        v_est = v_est_new
        n_iter += 1
    return v_est 

@njit
def estimate_cov(z1, z2, v1, v2, r, w):
    return np.sum(w*(z1*z2-r*np.sqrt(v1*v2)))

@njit
def estimate_cov_weights(v1, v2, r, var_1, var_2, c, l):
    w = np.power(((var_1+v1)*(var_2+v2)+np.power(c+r*np.sqrt(v1*v2),2))*l,-1)
    w = w/np.sum(w)
    return w

@njit
def reweighted_estimate_cov(z1, z2, v1, v2, r, l, var_1, var_2, tol=1e-8, max_iter=10**3):
    # Initialize weights
    w = np.power(((var_1+v1)*(var_2+v2)+np.power(r,2)*v1*v2)*l,-1)
    w = w/np.sum(w)
    # Initialize estimate
    c_est = estimate_cov(z1, z2, v1, v2, r, w)
    # Iterate to convergence
    c_diff = np.inf
    n_iter = 0
    while c_diff > tol and n_iter < max_iter:
        w = estimate_cov_weights(v1, v2, r, var_1, var_2, c_est, l)
        c_est_new = estimate_cov(z1, z2, v1, v2, r, w)
        c_diff = np.abs(c_est_new-c_est)
        c_est = c_est_new
        n_iter += 1
    return c_est 

@njit
def compute_corr(z1,z2,v1,v2,r,l):
    var_1 = reweighted_estimate_var(z1,v1,l)
    var_2 = reweighted_estimate_var(z2,v2,l)
    cov_12 = reweighted_estimate_cov(z1, z2, v1, v2, r, l, var_1, var_2)
    corr = cov_12/np.sqrt(var_1*var_2)
    reg_21 = cov_12/var_1
    resid_2 = z2-reg_21*z1
    v_resid_2 = v2+(reg_21**2)*v1-2*reg_21*r*np.sqrt(v1*v2)
    v_s = reweighted_estimate_var(resid_2, v_resid_2, l)
    return np.array([var_1, var_2, cov_12, corr, reg_21, v_s/var_2],dtype=np.float_)

@njit(parallel=True)
def jacknife(z1,z2,v1,v2,r,l,n_blocks,block_size):
    # Construct jacknife blocks
    jack_delete = np.zeros((n_blocks,6),dtype=np.float_)
    mask = np.ones((n_blocks, z1.shape[0]),dtype=np.bool_)
    for i in prange(n_blocks-1):
        mask[i, (block_size*i):(block_size*(i+1))] = False
    mask[n_blocks-1, (block_size*(n_blocks-1)):z1.shape[0]] = False
    # Compute jacknife values 
    for i in prange(n_blocks):
        jack_delete[i,:] = compute_corr(z1[mask[i,:]],z2[mask[i,:]],v1[mask[i,:]],v2[mask[i,:]],r[mask[i,:]],l[mask[i,:]])
    return jack_delete

def jacknife_est(z1,z2,v1,v2,r,l,n_blocks):
    # Estimate
    est = compute_corr(z1,z2,v1,v2,r,l)
    # Calculate blocks
    block_size = int(np.floor(z1.shape[0]/n_blocks))
    # Get jacknife ests
    jack_delete = jacknife(z1,z2,v1,v2,r,l,int(n_blocks),block_size)
    n_blocks = jack_delete.shape[0]
    # Compute jacknife-variance
    jack_vars = np.zeros((jack_delete.shape[1]),dtype=np.float_)
    for i in range(jack_delete.shape[1]):
        jack_vars[i] = ((n_blocks-1)/n_blocks)*np.sum(np.power(jack_delete[:,i]-np.mean(jack_delete[:,i]),2))
    # Adjust estimate of uncorrelated variance for error in regression coefficient
    est[5] = est[5]-est[0]*jack_vars[4]
    return est[3:6], np.sqrt(jack_vars[3:6]), jack_delete[:,3:6]
