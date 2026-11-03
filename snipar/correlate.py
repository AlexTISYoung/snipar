import numpy.ma as ma
from numba import njit, prange
import gzip
from snipar.ld import ldscores_from_bed
from snipar.utilities import make_id_dict
import numpy as np
import pandas as pd

SUMSTATS_COLUMN_ALIASES = {
    'SNP': ['SNP', 'snp', 'rsid', 'RSID', 'ID'],
    'pos': ['pos', 'BP', 'bp', 'position', 'POSITION'],
    'A1': ['A1', 'a1', 'allele1', 'effect_allele'],
    'A2': ['A2', 'a2', 'allele2', 'other_allele'],
    'freq': ['freq', 'f', 'FRQ', 'EAF', 'eaf'],
    'direct_Beta': ['direct_Beta', 'theta_direct'],
    'direct_SE': ['direct_SE', 'se_direct'],
    'avg_NTC_Beta': ['avg_NTC_Beta', 'avg_parental_Beta', 'theta_avgNTC', 'theta_avg_NTC'],
    'avg_NTC_SE': ['avg_NTC_SE', 'avg_parental_SE', 'se_avgNTC', 'se_avg_NTC'],
    'population_Beta': ['population_Beta', 'theta_population'],
    'population_SE': ['population_SE', 'se_population'],
    'r_direct_avg_NTC': [
        'r_direct_avg_NTC',
        'r_direct_avg_parental',
        'rg_direct_avgNTC',
        'rg_direct_avg_NTC',
    ],
    'r_direct_population': ['r_direct_population', 'rg_direct_population'],
}

SUMSTATS_CANONICAL_COLUMNS = [
    'SNP',
    'pos',
    'A1',
    'A2',
    'freq',
    'direct_Beta',
    'direct_SE',
    'avg_NTC_Beta',
    'avg_NTC_SE',
    'population_Beta',
    'population_SE',
    'r_direct_avg_NTC',
    'r_direct_population',
]
SUMSTATS_AVG_NTC_COLUMNS = [
    'avg_NTC_Beta',
    'avg_NTC_SE',
    'r_direct_avg_NTC',
]
SUMSTATS_DIRECT_POPULATION_COLUMNS = [
    col for col in SUMSTATS_CANONICAL_COLUMNS
    if col not in SUMSTATS_AVG_NTC_COLUMNS
]
SUMSTATS_STRING_COLUMNS = {'SNP', 'A1', 'A2'}

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
            try:
                ld_snps, ld_scores = read_ldscore_file_pandas(ld_files[i])
            except Exception as exc:  # pylint: disable=broad-except
                print('WARNING: Fast pandas LD-score reader failed for '+ld_files[i]+': '+str(exc))
                print('WARNING: Falling back to legacy np.loadtxt LD-score reader.')
                ld_snps, ld_scores = read_ldscore_file_loadtxt(ld_files[i])
            ld_indices = pd.Series(ld_snps).map(self.sid_dict)
            valid = ld_indices.notna().to_numpy() & np.isfinite(ld_scores)
            if np.sum(valid) > 0:
                valid_indices = ld_indices.loc[valid].to_numpy(dtype=int)
                self.ldscores[valid_indices] = ld_scores[valid]
                self.ldscores.mask[valid_indices] = False
            else:
                raise(ValueError('No SNPs in common between LD scores and summary statistics'))
    def cor_direct_pop(self, n_blocks, return_variances=False):
        print('Computing correlation between direct and population effects')
        return jacknife_est(self.direct,self.population,np.power(self.direct_SE,2),
                            np.power(self.population_SE,2),self.r_direct_pop,self.ldscores,n_blocks,
                            effect_1_name='direct effects',
                            effect_2_name='population effects',
                            return_variances=return_variances)
    def cor_direct_avg_NTC(self, n_blocks, return_variances=False):
        print('Computing correlation between direct effects and average NTCs')
        return jacknife_est(self.direct,self.avg_NTC,np.power(self.direct_SE,2),
                            np.power(self.avg_NTC_SE,2),self.r_direct_avg_NTC,self.ldscores, n_blocks,
                            effect_1_name='direct effects',
                            effect_2_name='average NTCs',
                            return_variances=return_variances)
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

def read_ldscore_file_pandas(ld_file):
    header = read_text_header(ld_file)
    if 'SNP' in header and 'L2' in header:
        data = pd.read_table(
            ld_file,
            sep=r'\s+',
            usecols=['SNP', 'L2'],
            dtype={'SNP': str},
            compression='infer',
            engine='c',
        )
        print('Interpreting LD-score columns in '+ld_file+': SNP <- SNP; L2 <- L2')
    else:
        data = pd.read_table(
            ld_file,
            sep=r'\s+',
            usecols=[1, 3],
            compression='infer',
            engine='c',
        )
        data.columns = ['SNP', 'L2']
        data['SNP'] = data['SNP'].astype(str)
        print('Interpreting LD-score columns in '+ld_file+': SNP <- column 2; L2 <- column 4')
    return data['SNP'].to_numpy(), pd.to_numeric(data['L2'], errors='coerce').to_numpy()

def read_ldscore_file_loadtxt(ld_file):
    data = np.loadtxt(ld_file,dtype=str,usecols=(1,3),skiprows=1)
    if data.ndim == 1:
        data = data.reshape((1, data.shape[0]))
    return data[:,0], np.array(data[:,1],dtype=float)

def read_text_header(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as handle:
        return handle.readline().strip().split()

def resolve_column_aliases(header, aliases, required_columns):
    exact_lookup = {name: name for name in header}
    lower_lookup = {}
    for name in header:
        lower_lookup.setdefault(name.lower(), name)

    resolved = {}
    missing = []
    for canonical in required_columns:
        selected = None
        for candidate in aliases[canonical]:
            if candidate in exact_lookup:
                selected = exact_lookup[candidate]
                break
            if candidate.lower() in lower_lookup:
                selected = lower_lookup[candidate.lower()]
                break
        if selected is None:
            missing.append(canonical)
        else:
            resolved[canonical] = selected

    if missing:
        raise(ValueError('Missing required column(s): '+', '.join(missing)))
    return resolved

def resolve_sumstats_columns(header):
    resolved = resolve_column_aliases(
        header,
        SUMSTATS_COLUMN_ALIASES,
        SUMSTATS_DIRECT_POPULATION_COLUMNS,
    )
    derived_columns = set()
    resolved_avg_ntc = {}
    missing_avg_ntc = []
    for canonical in SUMSTATS_AVG_NTC_COLUMNS:
        try:
            resolved_avg_ntc.update(
                resolve_column_aliases(header, SUMSTATS_COLUMN_ALIASES, [canonical])
            )
        except ValueError:
            missing_avg_ntc.append(canonical)

    if len(missing_avg_ntc) == 0:
        resolved.update(resolved_avg_ntc)
    else:
        derived_columns.update(SUMSTATS_AVG_NTC_COLUMNS)
        if len(missing_avg_ntc) < len(SUMSTATS_AVG_NTC_COLUMNS):
            print(
                'WARNING: Some average NTC columns are missing; deriving all average '
                'NTC columns from direct and population columns for consistency.'
            )
    return resolved, derived_columns

def print_column_interpretation(resolved_columns, file_label, derived_columns=None):
    if derived_columns is None:
        derived_columns = set()
    print('Interpreting columns in '+file_label+':')
    for canonical in SUMSTATS_CANONICAL_COLUMNS:
        if canonical in derived_columns:
            print('  '+canonical+' <- derived from direct/population columns')
        else:
            source = resolved_columns[canonical]
            if source == canonical:
                print('  '+canonical+' <- '+source)
            else:
                print('  '+canonical+' <- '+source+' (alias)')

def numeric_array(values):
    return np.asarray(pd.to_numeric(values, errors='coerce'), dtype=float)

def divide_or_nan(numerator, denominator):
    out = np.full(np.asarray(numerator).shape, np.nan, dtype=float)
    np.divide(
        numerator,
        denominator,
        out=out,
        where=np.isfinite(denominator) & (denominator != 0),
    )
    return out

def derive_avg_ntc_from_direct_population(columns):
    direct_beta = numeric_array(columns['direct_Beta'])
    direct_se = numeric_array(columns['direct_SE'])
    population_beta = numeric_array(columns['population_Beta'])
    population_se = numeric_array(columns['population_SE'])
    r_direct_population = numeric_array(columns['r_direct_population'])

    var_direct = np.power(direct_se, 2)
    var_population = np.power(population_se, 2)
    cov_direct_population = r_direct_population * direct_se * population_se
    var_avg_ntc = var_population + var_direct - 2.0 * cov_direct_population

    tol = 1e-12
    negative_var = np.isfinite(var_avg_ntc) & (var_avg_ntc < -tol)
    if np.any(negative_var):
        print(
            'WARNING: Derived average NTC sampling variance is negative for '
            +str(np.sum(negative_var))+
            ' SNP(s); setting derived average NTC SE and direct-average NTC '
            'sampling correlations to NA for those SNPs.'
        )
    var_avg_ntc = np.where(
        negative_var,
        np.nan,
        np.where(var_avg_ntc < 0, 0.0, var_avg_ntc),
    )
    avg_ntc_se = np.sqrt(var_avg_ntc)
    cov_direct_avg_ntc = cov_direct_population - var_direct

    columns['avg_NTC_Beta'] = population_beta - direct_beta
    columns['avg_NTC_SE'] = avg_ntc_se
    columns['r_direct_avg_NTC'] = divide_or_nan(
        cov_direct_avg_ntc,
        direct_se * avg_ntc_se,
    )
    return columns

def complete_sumstats_columns(columns, derived_columns):
    if any(col in derived_columns for col in SUMSTATS_AVG_NTC_COLUMNS):
        columns = derive_avg_ntc_from_direct_population(columns)
    return columns

def sumstats_from_columns(chrom, columns):
    return sumstats(chrom,
                    columns['SNP'],
                    columns['pos'],
                    columns['A1'],
                    columns['A2'],
                    columns['freq'],
                    columns['direct_Beta'],
                    columns['direct_SE'],
                    columns['avg_NTC_Beta'],
                    columns['avg_NTC_SE'],
                    columns['population_Beta'],
                    columns['population_SE'],
                    columns['r_direct_avg_NTC'],
                    columns['r_direct_population'])

def read_sumstats_file_pandas(sumstats_file, chrom, resolved_columns, derived_columns):
    read_columns = [col for col in SUMSTATS_CANONICAL_COLUMNS if col in resolved_columns]
    usecols = list(dict.fromkeys([resolved_columns[col] for col in read_columns]))
    dtype = {resolved_columns[col]: str for col in SUMSTATS_STRING_COLUMNS}
    data = pd.read_table(
        sumstats_file,
        sep=r'\s+',
        usecols=usecols,
        dtype=dtype,
        compression='infer',
        engine='c',
    )
    columns = {}
    for canonical in read_columns:
        source = resolved_columns[canonical]
        if canonical in SUMSTATS_STRING_COLUMNS:
            columns[canonical] = data[source].astype(str).to_numpy()
        else:
            columns[canonical] = pd.to_numeric(data[source], errors='coerce').to_numpy()
    columns = complete_sumstats_columns(columns, derived_columns)
    return sumstats_from_columns(chrom, columns)

def read_sumstats_file_loadtxt(sumstats_file, chrom, header, resolved_columns, derived_columns):
    read_columns = [col for col in SUMSTATS_CANONICAL_COLUMNS if col in resolved_columns]
    col_indices = tuple([header.index(resolved_columns[col]) for col in read_columns])
    data = np.loadtxt(sumstats_file,dtype=str,skiprows=1,usecols=col_indices)
    if data.ndim == 1:
        data = data.reshape((1, data.shape[0]))
    columns = {canonical: data[:,i] for i, canonical in enumerate(read_columns)}
    columns = complete_sumstats_columns(columns, derived_columns)
    return sumstats_from_columns(chrom, columns)

def read_sumstats_file(sumstats_file, chrom):
    print('Reading sumstats from '+sumstats_file)
    header = read_text_header(sumstats_file)
    resolved_columns, derived_columns = resolve_sumstats_columns(header)
    print_column_interpretation(resolved_columns, sumstats_file, derived_columns)
    try:
        return read_sumstats_file_pandas(sumstats_file, chrom, resolved_columns, derived_columns)
    except Exception as exc:  # pylint: disable=broad-except
        print('WARNING: Fast pandas sumstats reader failed for '+sumstats_file+': '+str(exc))
        print('WARNING: Falling back to legacy np.loadtxt sumstats reader.')
        return read_sumstats_file_loadtxt(sumstats_file, chrom, header, resolved_columns, derived_columns)

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

def _mask_outputs_with_negative_variance(est, se, delete, internal_est, residual_var_2,
                                         effect_1_name='first effects',
                                         effect_2_name='second effects'):
    negative_variances = []
    if np.isfinite(internal_est[0]) and internal_est[0] < 0:
        negative_variances.append((effect_1_name, internal_est[0]))
    if np.isfinite(internal_est[1]) and internal_est[1] < 0:
        negative_variances.append((effect_2_name, internal_est[1]))
    if np.isfinite(residual_var_2) and residual_var_2 < 0:
        negative_variances.append(('residual '+effect_2_name, residual_var_2))

    if len(negative_variances) == 0:
        return est, se, delete

    output_names = [
        'correlation between '+effect_1_name+' and '+effect_2_name,
        'regression coefficient of '+effect_2_name+' on '+effect_1_name,
        'proportion of variance in '+effect_2_name+' uncorrelated with '+effect_1_name,
    ]
    invalid_outputs = set()
    if np.isfinite(internal_est[0]) and internal_est[0] < 0:
        invalid_outputs.update([0, 1, 2])
    if np.isfinite(internal_est[1]) and internal_est[1] < 0:
        invalid_outputs.update([0, 2])
    if np.isfinite(residual_var_2) and residual_var_2 < 0:
        invalid_outputs.add(2)

    negative_msg = '; '.join([name+' variance = '+str(value) for name, value in negative_variances])
    invalid_msg = '; '.join([output_names[i] for i in sorted(invalid_outputs)])
    print('WARNING: Negative variance estimate(s) detected in correlate.py: '+negative_msg)
    print('WARNING: Setting the following correlate.py output estimate(s), SE(s), and delete value(s) to NA because they rely on negative variance estimate(s): '+invalid_msg)

    est = est.copy()
    se = se.copy()
    delete = delete.copy()
    for i in invalid_outputs:
        est[i] = np.nan
        se[i] = np.nan
        delete[:,i] = np.nan
    return est, se, delete

def _jacknife_variance(delete_values):
    n_blocks = delete_values.shape[0]
    return ((n_blocks-1)/n_blocks)*np.sum(
        np.power(delete_values-np.mean(delete_values, axis=0), 2),
        axis=0,
    )

def _variance_outputs(est, jack_delete, jack_vars):
    variance_est = np.array([est[0], est[1], np.nan], dtype=float)
    variance_se = np.array([
        np.sqrt(jack_vars[0]),
        np.sqrt(jack_vars[1]),
        np.nan,
    ])
    variance_delete = np.column_stack((
        jack_delete[:,0],
        jack_delete[:,1],
        np.full(jack_delete.shape[0], np.nan),
    ))

    if np.isfinite(est[0]) and est[0] >= 0 and np.isfinite(est[1]) and est[1] > 0:
        variance_est[2] = est[0]/est[1]
        valid_delete = np.logical_and(
            np.isfinite(jack_delete[:,0]) & (jack_delete[:,0] >= 0),
            np.isfinite(jack_delete[:,1]) & (jack_delete[:,1] > 0),
        )
        variance_delete[valid_delete,2] = (
            jack_delete[valid_delete,0]/jack_delete[valid_delete,1]
        )
        if np.all(valid_delete):
            variance_se[2] = np.sqrt(_jacknife_variance(variance_delete[:,2]))
        else:
            print(
                'WARNING: At least one delete-block direct/population variance '
                'ratio was undefined because a corrected variance was non-positive; '
                'setting its jackknife SE to NA.'
            )
    else:
        print(
            'WARNING: Corrected variance ratio is undefined because a full-sample '
            'corrected variance was non-positive; setting its estimate, SE, and '
            'delete values to NA.'
        )
    return variance_est, variance_se, variance_delete

def jacknife_est(z1,z2,v1,v2,r,l,n_blocks,effect_1_name='first effects',effect_2_name='second effects',return_variances=False):
    # Estimate
    est = compute_corr(z1,z2,v1,v2,r,l)
    # Calculate blocks
    block_size = int(np.floor(z1.shape[0]/n_blocks))
    # Get jacknife ests
    jack_delete = jacknife(z1,z2,v1,v2,r,l,int(n_blocks),block_size)
    # Compute jacknife-variance
    jack_vars = _jacknife_variance(jack_delete)
    # Adjust estimate of uncorrelated variance for error in regression coefficient
    residual_var_2 = est[5]*est[1]
    est[5] = est[5]-est[0]*jack_vars[4]
    correlation_outputs = _mask_outputs_with_negative_variance(
        est[3:6],
        np.sqrt(jack_vars[3:6]),
        jack_delete[:,3:6],
        est,
        residual_var_2,
        effect_1_name,
        effect_2_name,
    )
    if not return_variances:
        return correlation_outputs
    return correlation_outputs + _variance_outputs(est, jack_delete, jack_vars)
