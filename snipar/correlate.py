import numpy.ma as ma
from numba import njit, prange
import gzip
from snipar.ld import ldscores_from_bed
from snipar.utilities import make_id_dict
import numpy as np

class sumstats(object):
    def __init__(self,chrom,sid,pos,A1,A2,freqs,direct,avg_NTC,population,r_direct_avg_NTC,r_direct_pop,ldscores = None, map=None):
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
        self.avg_NTC = ma.array(avg_NTC, dtype=float)
        self.avg_NTC.mask = np.isnan(self.avg_NTC)
        self.population = ma.array(population, dtype=float)
        self.population.mask = np.isnan(self.population)
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
        self.avg_NTC = ma.concatenate([self.avg_NTC, s2.avg_NTC])
        self.population = ma.concatenate([self.population, s2.population])
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
        self.avg_NTC = self.avg_NTC[filter_pass]
        self.population = self.population[filter_pass]
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
        r_dir_pop, r_dir_pop_SE, r_dir_pop_delete = jacknife_est(self.direct,self.population,self.r_direct_pop,1/self.ldscores, n_blocks)
        return r_dir_pop, r_dir_pop_SE, r_dir_pop_delete
    def cor_direct_avg_NTC(self, n_blocks):
        print('Computing correlation between direct effects and average NTCs')
        r_dir_avg_NTC, r_dir_avg_NTC_SE, r_dir_avg_NTC_delete = jacknife_est(self.direct,self.avg_NTC,self.r_direct_avg_NTC,1/self.ldscores, n_blocks)
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
    # Find columns
    sid_index = np.where(sumstats_header=='SNP')[0][0]
    pos_index = np.where(sumstats_header=='pos')[0][0]
    A1_index = np.where(sumstats_header=='A1')[0][0]
    A2_index = np.where(sumstats_header=='A2')[0][0]
    freq_index = np.where(sumstats_header=='freq')[0][0]
    direct_index = np.where(sumstats_header=='direct_Z')[0][0]
    avg_NTC_index = np.where(sumstats_header=='avg_NTC_Z')[0][0]
    population_index = np.where(sumstats_header=='population_Z')[0][0]
    r_direct_avg_NTC_index = np.where(sumstats_header=='r_direct_avg_NTC')[0][0]
    r_direct_pop_index = np.where(sumstats_header=='r_direct_population')[0][0]
    # Read summary statistics
    s = np.loadtxt(sumstats_file,dtype=str,skiprows=1)
    # Return sumstats class
    return sumstats(chrom, s[:,sid_index],s[:,pos_index],s[:,A1_index],s[:,A2_index],s[:,freq_index],
                    s[:,direct_index],s[:,avg_NTC_index],s[:,population_index],
                    s[:,r_direct_avg_NTC_index],s[:,r_direct_pop_index])

def read_sumstats_files(sumstats_files, chroms):
    s = read_sumstats_file(sumstats_files[0], chroms[0])
    if chroms.shape[0]>1:
        for i in range(1,sumstats_files.shape[0]):
            s.concatenate(read_sumstats_file(sumstats_files[i], chroms[i]))
    return s

@njit
def compute_corr(z1,z2,r,w):
    return np.sum(w*(z1*z2-r))/np.sqrt(np.sum(w*(np.power(z1,2)-1))*np.sum(w*(np.power(z2,2)-1)))

@njit(parallel=True)
def jacknife(z1,z2,r,w,n_blocks,block_size):
    jack_delete = np.zeros((n_blocks),dtype=np.float_)
    mask = np.ones((n_blocks, z1.shape[0]),dtype=np.bool_)
    for i in prange(n_blocks-1):
        mask[i, (block_size*i):(block_size*(i+1))] = False
    mask[n_blocks-1, (block_size*(n_blocks-1)):z1.shape[0]] = False 
    for i in prange(n_blocks-1):
        jack_delete[i] = compute_corr(z1[mask[i,:]],z2[mask[i,:]],r[mask[i,:]],w[mask[i,:]])
    jack_delete[n_blocks-1] = compute_corr(z1[mask[n_blocks-1,:]],z2[mask[n_blocks-1,:]],r[mask[n_blocks-1,:]],w[mask[n_blocks-1,:]])
    return jack_delete

def jacknife_est(z1,z2,r,w,n_blocks):
    est = compute_corr(z1,z2,r,w)
    block_size = int(np.floor(z1.shape[0]/n_blocks))
    jack_delete = jacknife(z1,z2,r,w,int(n_blocks),block_size)
    n_blocks = jack_delete.shape[0]
    jack_var = ((n_blocks-1)/n_blocks)*np.sum(np.power(jack_delete-np.mean(jack_delete),2))
    return est, np.sqrt(jack_var), jack_delete