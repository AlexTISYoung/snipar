import numpy as np
import numpy.ma as ma
from snipar.utilities import make_id_dict

class gtarray(object):
    """Define a genotype or PGS array that stores individual IDs, family IDs, and SNP information.

    Args:
        garray : :class:`~numpy:numpy.array`
            2 or 3 dimensional numpy array of genotypes/PGS values. First dimension is individuals. For a 2 dimensional array, the second dimension is SNPs or PGS values.
            For a 3 dimensional array, the second dimension indexes the individual and his/her relatives' genotypes (for example: proband, paternal, and maternal); and
            the third dimension is the SNPs.
        ids : :class:`~numpy:numpy.array`
            vector of individual IDs of same length as first dimension of garray
        sid : :class:`~numpy:numpy.array`
            vector of SNP ids, equal in length, L, to last dimension of array
        alleles : :class:`~numpy:numpy.array`
            [L x 2] matrix of ref and alt alleles for the SNPs. L must match size of sid
        pos : :class:`~numpy:numpy.array`
            vector of SNP positions; must match size of sid
        chrom : :class:`~numpy:numpy.array`
            vector of SNP chromosomes; must match size of sid
        map : :class:`~numpy:numpy.array`
            vector of SNP chromosomes; must match size of sid
        fams : :class:`~numpy:numpy.array`
            vector of family IDs; must match size of ids
        par_status : :class:`~numpy:numpy.array'
             [N x 2] numpy matrix that records whether parents have observed or imputed genotypes/PGS, where N matches size of ids.
             The first column is for the father of that individual; the second column is for the mother of that individual.
             If the parent is neither observed nor imputed, the value is -1; if observed, 0; and if imputed, 1.

    Returns:
        G : :class:`snipar.gtarray`

    """
    def __init__(self, garray, ids, sid=None, alleles=None, pos=None, chrom=None, map=None, error_probs=None, fams=None, par_status=None):
        if type(garray) == np.ndarray or type(garray) == np.ma.core.MaskedArray:
            if type(garray) == np.ndarray:
                self.gts = ma.array(garray,mask=np.isnan(garray))
            else:
                self.gts = garray
            self.shape = garray.shape
            self.ndim = garray.ndim
            self.dtype = garray.dtype
            self.freqs = None
        else:
            raise ValueError('Genotypes must be a numpy ndarray')
        if garray.shape[0] == ids.shape[0]:
            self.ids = ids
            self.id_dict = make_id_dict(ids)
        else:
            raise ValueError('Shape of genotypes and ids does not match')
        if sid is not None:
            if sid.shape[0] == garray.shape[1]:
                self.snp_index = 1
                self.sid = sid
                self.sid_dict = make_id_dict(sid)
            elif sid.shape[0] == garray.shape[2]:
                self.snp_index = 2
                self.sid = sid
                self.sid_dict = make_id_dict(sid)
            else:
                raise ValueError('Shape of SNP ids (sid) does not match shape of genotype array')
        if alleles is not None:
            if self.sid is not None:
                if alleles.shape[0] == self.sid.shape[0]:
                    self.alleles = alleles
                else:
                    raise ValueError('Size of alleles does not match size of genotypes')
            else:
                raise(ValueError('Must provide SNP ids'))
        else:
            self.alleles = None
        if pos is not None:
            if self.sid is not None:
                if pos.shape[0] == self.sid.shape[0]:
                    self.pos = pos
                else:
                    raise ValueError('Size of position vector does not match size of genotypes')
            else:
                raise(ValueError('Must provide SNP ids'))
        else:
            self.pos = None
        if chrom is not None:
            if self.sid is not None:
                if chrom.shape[0] == self.sid.shape[0]:
                    self.chrom = chrom
                else:
                    raise ValueError('Size of map does not match number of SNPs')
            else:
                raise(ValueError('Must provide SNP ids'))
        else:
            self.chrom = None
        if map is not None:
            if self.sid is not None:
                if map.shape[0] == self.sid.shape[0]:
                    self.map = map
                else:
                    raise ValueError('Size of map does not match number of SNPs')
            else:
                raise(ValueError('Must provide SNP ids'))
        else:
            self.map = None
        if error_probs is not None:
            if self.sid is not None:
                if error_probs.shape[0] == self.sid.shape[0]:
                    self.error_probs = error_probs
                else:
                    raise ValueError('Size of map does not match number of SNPs')
            else:
                raise(ValueError('Must provide SNP ids'))
        else:
            self.error_probs = None
        if fams is not None:
            if fams.shape[0] == ids.shape[0] and fams.ndim==1:
                self.fams = fams
            else:
                raise ValueError('Fams not of same length as IDs')
        else:
            self.fams = None

        if par_status is not None:
            if par_status.shape[0] == ids.shape[0] and par_status.shape[1] == 2:
                self.par_status = par_status
            else:
                raise ValueError('Incompatible par status array')
        else:
            self.par_status = None

        self.mean_normalised = False

        if np.sum(self.gts.mask)>0:
            self.has_NAs = True
        else:
            self.has_NAs = False

        self.info = None

    def compute_freqs(self):
        """
        Computes the frequencies of the SNPs. Stored in self.freqs.
        """
        if self.ndim == 2:
            self.freqs = ma.mean(self.gts,axis=0)/2.0
        elif self.ndim == 3:
            self.freqs = ma.mean(self.gts[:,0,:], axis=0) / 2.0

    def filter(self, filter_pass):
        if self.freqs is not None:
            self.freqs = self.freqs[filter_pass]
        if self.ndim == 2:
            self.gts = self.gts[:,filter_pass]
        elif self.ndim == 3:
            self.gts = self.gts[:,:,filter_pass]
        self.shape = self.gts.shape
        if self.sid is not None:
            self.sid = self.sid[filter_pass]
            self.sid_dict = make_id_dict(self.sid)
        if self.pos is not None:
            self.pos = self.pos[filter_pass]
        if self.alleles is not None:
            self.alleles = self.alleles[filter_pass]
        if self.chrom is not None:
            self.chrom = self.chrom[filter_pass]
        if self.map is not None:
            self.map = self.map[filter_pass]
        if self.error_probs is not None:
            self.error_probs = self.error_probs[filter_pass]

    def filter_maf(self, min_maf = 0.01, verbose=False):
        """
        Filter SNPs based on having minor allele frequency (MAF) greater than min_maf, and have % missing observations less than max_missing.
        """
        if self.freqs is None:
            self.compute_freqs()
        freqs_pass = np.logical_and(self.freqs > min_maf, self.freqs < (1 - min_maf))
        if verbose:
            print(str(self.freqs.shape[0] - np.sum(freqs_pass)) + ' SNPs with MAF<' + str(min_maf))
        self.filter(freqs_pass)

    def filter_missingness(self, max_missing = 5, verbose=False):
        if self.ndim == 2:
            missingness = ma.mean(self.gts.mask,axis=0)
        elif self.ndim == 3:
            missingness = ma.mean(self.gts.mask,axis = (0,1))
        missingness_pass = 100 * missingness < max_missing
        if verbose:
            print(str(self.freqs.shape[0] - np.sum(missingness_pass)) + ' SNPs with missingness >' + str(max_missing) + '%')
        self.filter(missingness_pass)

    def compute_info(self):
        if self.freqs is None:
            self.compute_freqs()
        if self.ndim == 2:
            self.variances = np.var(self.gts, axis=0)
        elif self.ndim==3:
            self.variances = np.var(self.gts[:,0,:], axis=0)
        self.info = self.variances/(2.0*self.freqs*(1-self.freqs))

    def filter_info(self, min_info = 0.99, verbose=False):
        if self.info is None:
            self.compute_info()
        info_pass = self.info > min_info
        if verbose:
            print(str(self.info.shape[0] - np.sum(info_pass)) + ' SNPs with INFO <' + str(min_info))
        self.filter(info_pass)

    def filter_ids(self,keep_ids, verbose=False):
        """
        Keep only individuals with ids given by keep_ids
        """
        in_ids = np.array([x in self.id_dict for x in keep_ids])
        n_filtered = np.sum(in_ids)
        if n_filtered==0:
            raise(ValueError('No individuals would be left after filtering'))
        else:
            if verbose:
                print('After filtering, '+str(n_filtered)+' individuals remain')
            indices = np.array([self.id_dict[x] for x in keep_ids[in_ids]])
            if self.ndim == 2:
                self.gts = self.gts[indices, :]
            elif self.ndim == 3:
                self.gts = self.gts[indices, :, :]
            self.ids = self.ids[indices]
            self.id_dict = make_id_dict(self.ids)
            self.shape = self.gts.shape
            if self.fams is not None:
                self.fams = self.fams[indices]

    def mean_normalise(self):
        """
        This normalises the SNPs/PGS columns to have mean-zero.
        """
        if not self.mean_normalised:
            if self.gts.ndim == 2:
                self.gts = self.gts - ma.mean(self.gts,axis=0)
            elif self.gts.ndim==3:
                for i in range(0, self.gts.shape[1]):
                    self.gts[:, i, :] = self.gts[:, i, :] - ma.mean(self.gts[:, i, :], axis=0)
            self.mean_normalised = True

    def scale(self):
        """
        This normalises the SNPs/PGS columns to have variance 1.
        """
        if self.gts.ndim == 2:
            self.gts = self.gts/ma.std(self.gts, axis=0)
        elif self.gts.ndim == 3:
            for i in range(0, self.gts.shape[1]):
                self.gts[:, i, :] = self.gts[:, i, :]/ma.std(self.gts[:, i, :], axis=0)

    def fill_NAs(self):
        """
        This normalises the SNP columns to have mean-zero, then fills in NA values with zero.
        """
        if not self.mean_normalised:
            self.mean_normalise()
        NAs = np.sum(self.gts.mask, axis=0)
        self.gts[self.gts.mask] = 0
        self.gts.mask = False
        self.has_NAs = False
        return NAs


    def add(self,garray):
        """
        Adds another gtarray of the same dimension to this array and returns the sum. It matches IDs before summing.
        """
        if type(garray)==gtarray:
            pass
        else:
            raise ValueError('Must add to another gtarray')

        if not self.gts.ndim == garray.gts.ndim:
            raise ValueError('Arrays must have same number of dimensions')

        if self.gts.ndim == 2:
            if not self.gts.shape[1] == garray.gts.shape[1]:
                raise ValueError('Arrays must have same dimensions (apart from first)')

        if self.gts.ndim == 3:
            if not self.gts.shape[1:3] == garray.gts.shape[1:3]:
                raise ValueError('Arrays must have same dimensions (apart from first)')

        # Match IDs
        common_ids = list(self.id_dict.keys() & garray.id_dict.keys())
        if len(common_ids) == 0:
            raise ValueError('No IDs in common')
        self_index = np.array([self.id_dict[x] for x in common_ids])
        other_index = np.array([garray.id_dict[x] for x in common_ids])

        # Out
        if self.ids.ndim == 1:
            ids_out = self.ids[self_index]
        else:
            ids_out = self.ids[self_index, :]

        if self.gts.ndim ==2:
            add_gts = self.gts[self_index, :]+garray.gts[other_index, :]
        else:
            add_gts = self.gts[self_index, :, :] + garray.gts[other_index, :, :]

        return gtarray(add_gts, ids_out, self.sid, alleles=self.alleles, fams=self.fams[self_index])

    def diagonalise(self,inv_root):
        """
        This will transform the genotype array based on the inverse square root of the phenotypic covariance matrix
        from the family based linear mixed model.
        """
        if self.fams is None:
            raise(ValueError('Family labels needed for diagonalization'))
        if not self.mean_normalised:
            self.mean_normalise()
        unique_fams, famsizes = np.unique(self.fams, return_counts = True)
        fam_indices = dict()
        # Transform
        for fam in unique_fams:
            fam_indices[fam] = np.where(self.fams == fam)[0]
            famsize = fam_indices[fam].shape[0]
            if self.ndim == 2:
                if famsize == 1:
                    self.gts[fam_indices[fam], :] = inv_root[1]*self.gts[fam_indices[fam],:]
                else:
                    self.gts[fam_indices[fam],:] = inv_root[famsize].dot(self.gts[fam_indices[fam],:])
            elif self.ndim == 3:
                if famsize == 1:
                    self.gts[fam_indices[fam], : , :] = inv_root[1]*self.gts[fam_indices[fam], : , :]
                else:
                    for j in range(self.shape[1]):
                        self.gts[fam_indices[fam],j, :] = inv_root[famsize].dot(self.gts[fam_indices[fam],j, :])
        self.fam_indices = fam_indices
