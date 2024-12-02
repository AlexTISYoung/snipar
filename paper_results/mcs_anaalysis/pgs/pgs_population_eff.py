import argparse
import numpy as np
from snipar.gtarray import gtarray
import snipar.read as read
from snipar.read.bed import get_snps
import snipar.slmm as lmm
from snipar.utilities import *
from pysnptools.snpreader import Bed
from bgen_reader import open_bgen
from snipar.read import get_gts_matrix

def write(pg,filename,scale_PGS = False):
    if scale_PGS:
        # Rescale by observed proband PGS
        pg.gts = pg.gts / np.std(pg.gts[:, 0])
    ####### Write PGS to file ########
    pg_out = np.column_stack((pg.fams,pg.ids,pg.gts))
    pg_header = np.column_stack((np.array(['FID','IID']).reshape(1,2),pg.sid.reshape(1,pg.sid.shape[0])))
    pg_out = np.row_stack((pg_header,pg_out))
    print('Writing PGS to ' + filename)
    np.savetxt(filename, pg_out, fmt='%s')
    return None


class pgs(object):
    """Define a polygenic score based on a set of SNPs with weights and ref/alt allele pairs.
    Args:
        snp_ids : :class:`~numpy:numpy.array`
            [L] vector of SNP ids
        weights : :class:`~numpy:numpy.array`
            [L] vector of weights of equal length to snp_ids
        alleles : :class:`~numpy:numpy.array`
            [L x 2] matrix of ref and alt alleles for the SNPs. L must match size of snp_ids
    Returns:
        pgs : :class:`snipar.pgs`
    """
    def __init__(self,snp_ids,weights,alleles):
        if snp_ids.shape[0] == weights.shape[0] and alleles.shape[0] == weights.shape[0] and alleles.shape[1]==2:
            self.snp_ids = snp_ids
            self.snp_dict = make_id_dict(snp_ids)
            self.weights = weights
            self.alleles = alleles
        else:
            raise ValueError('All inputs must have the same dimension')

    def compute(self, garray, cols=None):
        """Compute polygenic score values from a given genotype array. Finds the SNPs in the genotype array
        that have weights in the pgs and matching alleles, and computes the PGS based on these SNPs and the
        weights after allele-matching.
        Args:
            garray : :class:`sbreg.gtarray`
                genotype array to compute PGS values for
            cols : :class:`numpy:numpy.array`
                names to give the columns in the output gtarray
        Returns:
            pg : :class:`snipar.gtarray`
                2d gtarray with PGS values. If a 3d gtarray is input, then each column corresponds to
                the second dimension on the input gtarray (for example, individual, paternal, maternal PGS).
                If a 2d gtarray is input, then there will be only one column in the output gtarray. The
                names given in 'cols' are stored in 'sid' attribute of the output.
        """
        if type(garray) == gtarray:
            garray.fill_NAs()
        else:
            raise ValueError('Must be of gtarray class')
        if garray.alleles is None:
            raise ValueError('Alleles of genotype matrix must be provided')
        # Match SNP IDs
        in_pgs_snps = np.array([x in self.snp_dict for x in garray.sid])
        nmatch = np.sum(in_pgs_snps)
        if nmatch==0:
            print('No overlap between PGS SNPs and genotype SNPs')
            return None
        else:
            # Get weights
            matched_snps = garray.sid[in_pgs_snps]
            matched_alleles = garray.alleles[in_pgs_snps,:]
            snp_indices = np.zeros((nmatch),dtype=int)
            for i in range(0,nmatch):
                snp_indices[i] = self.snp_dict[matched_snps[i]]
            weights_compute = self.weights[snp_indices]
            alleles = self.alleles[snp_indices,:]

            # Match alleles and adjust weights
            a_match = np.logical_and(alleles[:,0] == matched_alleles[:, 0], alleles[:,1] == matched_alleles[:, 1])
            a_reverse = np.logical_and(alleles[:,0] == matched_alleles[:, 1], alleles[:,1] == matched_alleles[:, 0])
            a_nomatch = np.logical_and(np.logical_not(a_match), np.logical_not(a_reverse))
            n_nomatch = np.sum(a_nomatch)
            if n_nomatch > 0:
                print('Removing ' + str(n_nomatch) + ' SNPs due to allele mismatch between genotypes and PGS alleles')
                weights_compute[a_nomatch] = 0
            weights_compute[a_reverse] = -weights_compute[a_reverse]

            ### Compute PGS
            if garray.ndim == 2:
                pgs_val = garray.gts[:, in_pgs_snps].dot(weights_compute).reshape(-1, 1)
            elif garray.ndim == 3:
                pgs_val = np.zeros((garray.gts.shape[0], garray.gts.shape[1]), garray.dtype)
                for i in range(0, garray.gts.shape[1]):
                    pgs_val[:, i] = garray.gts[:, i, in_pgs_snps].dot(weights_compute)
            # return gtarray(pgs_val, garray.ids, sid=cols, fams=garray.fams if garray.fams is not None else garray.ids)
            return gtarray(pgs_val, garray.ids, sid=cols, fams=garray.fams if garray.fams is not None else garray.ids, par_status=garray.par_status)
        
######### Command line arguments #########
parser=argparse.ArgumentParser()
parser.add_argument('out',type=str,help='Prefix for computed PGS file and/or regression results files')
parser.add_argument('--bgen',
                    type=str,help='Address of the phased genotypes in .bgen format. If there is a @ in the address, @ is replaced by the chromosome numbers in the range of chr_range for each chromosome (chr_range is an optional parameters for this script).')
parser.add_argument('--bed',
                    type=str,help='Address of the unphased genotypes in .bed format. If there is a @ in the address, @ is replaced by the chromosome numbers in the range of chr_range for each chromosome (chr_range is an optional parameters for this script).')
parser.add_argument('--imp', type=str, help='Address of hdf5 files with imputed parental genotypes (without .hdf5 suffix). If there is a @ in the address, @ is replaced by the chromosome numbers in the range of chr_range (chr_range is an optional parameters for this script).', default = None)
parser.add_argument('--chr_range',
                    type=parseNumRange,
                    nargs='*',
                    action=NumRangeAction,
                    help='number of the chromosomes to be imputed. Should be a series of ranges with x-y format or integers.', default=None)
parser.add_argument('--pedigree',type=str,help='Address of pedigree file. Must be provided if not providing imputed parental genotypes.',default=None)
parser.add_argument('--weights',type=str,help='Location of the PGS allele weights', default = None)
parser.add_argument('--SNP',type=str,help='Name of column in weights file with SNP IDs',default='sid')
parser.add_argument('--beta_col',type=str,help='Name of column with betas/weights for each SNP',default='ldpred_beta')
parser.add_argument('--A1',type=str,help='Name of column with allele beta/weights are given with respect to',default='nt1')
parser.add_argument('--A2',type=str,help='Name of column with alternative allele',default='nt2')
parser.add_argument('--sep',type=str,help='Column separator in weights file. If not provided, an attempt to determine this will be made.',default=None)
parser.add_argument('--phenofile',type=str,help='Location of the phenotype file',default = None)
parser.add_argument('--pgs', type=str, help='Location of the pre-computed PGS file', default=None)
parser.add_argument('--fit_sib', action='store_true', default=False, help='Fit indirect effects from siblings')
parser.add_argument('--parsum',action='store_true',default = False, help='Use the sum of maternal and paternal PGS in the regression (useful when imputed from sibling data alone)')
parser.add_argument('--phen_index',type=int,help='If the phenotype file contains multiple phenotypes, which phenotype should be analysed (default 1, first)',
                    default=1)
parser.add_argument('--scale_phen',action='store_true',help='Scale the phenotype to have variance 1',default=False)
parser.add_argument('--scale_pgs',action='store_true',help='Scale the PGS to have variance 1 among the phenotyped individuals',default=False)
parser.add_argument('--compute_controls', action='store_true', default=False,
                    help='Compute PGS for control families (default False)')
parser.add_argument('--missing_char',type=str,help='Missing value string in phenotype file (default NA)',default='NA')


def main(args):
    """"Calling this function with args is equivalent to running this script from commandline with the same arguments.
    Args:
        args: list
            list of all the desired options and arguments. The possible values are all the values you can pass this script from commandline.
    """
    if args.weights is not None:
        if args.bed is None and args.bgen is None:
            raise ValueError('Weights provided but no observed genotypes provided')
        if args.bed is not None and args.bgen is not None:
            raise ValueError('Provide only one of --bedfiles and --bgenfiles')
        print('Computing PGS from weights file')
        ####### Read PGS #######
        if args.sep is None:
            weights = np.loadtxt(args.weights,dtype=str)
        else:
            weights = np.loadtxt(args.weights,dtype=args.sep)
        colnames = weights[0,:]
        weights = weights[1:weights.shape[0],:]
        print('Read weights for '+str(weights.shape[0])+' variants')
        beta = np.array(weights[:,np.where(colnames == args.beta_col)[0][0]],dtype=np.float64)
        allele_indices = np.array([np.where(colnames==args.A1)[0][0],np.where(colnames==args.A2)[0][0]])
        p = pgs(weights[:,np.where(colnames==args.SNP)[0][0]],
                beta,
                weights[:,allele_indices])
        bedfiles, chroms = parse_obsfiles(args.bed, 'bed', chromosomes=args.chr_range)
        # bedfiles, pargts_list, chroms = parse_filelist(args.bed, args.imp, 'bed', chromosomes=args.chr_range)
        pg = compute(p, bedfile=bedfiles[0])#, par_gts_f=pargts_list[0])
        for i in range(1,chroms.shape[0]):
            if bedfiles[i] is not None:
                print('Observed genotype file: '+bedfiles[i])
            pg_i = compute(p, bedfile=bedfiles[i])#, par_gts_f=pargts_list[i])
            if pg_i is not None:
                if pg is not None:
                    pg = pg.add(pg_i)
                else:
                        pg = pg_i
        write(pg, args.out + '.pgs.txt', scale_PGS=args.scale_pgs)
    elif args.pgs is not None:
        if args.phenofile is None:
            raise ValueError('Pre-computed PGS provided but no phenotype provided')
        print('Reading PGS from '+args.pgs)
        pgs_f = open(args.pgs,'r')
        pgs_header = pgs_f.readline().split(' ')
        pgs_header[len(pgs_header)-1] = pgs_header[len(pgs_header)-1].split('\n')[0]
        ncols = len(pgs_header)
        pgs_cols = tuple([x for x in range(2,ncols)])
        pg = gtarray(np.loadtxt(args.pgs,usecols = pgs_cols, skiprows=1),
                     np.loadtxt(args.pgs,usecols = 1, dtype=str, skiprows=1),
                     sid=np.array(pgs_header[2:ncols]),
                     fams=np.loadtxt(args.pgs,usecols = 0, dtype=str, skiprows=1))
    else:
        raise ValueError('Weights or PGS must be provided')

def compute(pgs, bedfile=None, bgenfile=None, par_gts_f=None, ped=None, sib=False, compute_controls=False, verbose=True):
    """Compute a polygenic score (PGS) for the individuals with observed genotypes and observed/imputed parental genotypes.
    Args:
        par_gts_f : :class:`str`
            path to HDF5 file with imputed parental genotypes
        gts_f : :class:`str`
            path to bed file with observed genotypes
        pgs : :class:`snipar.pgs`
            the PGS, defined by the weights for a set of SNPs and the alleles of those SNPs
        sib : :class:`bool`
            Compute the PGS for genotyped individuals with at least one genotyped sibling and observed/imputed parental genotypes. Default False.
        compute_controls : :class:`bool`
            Compute polygenic scores for control families (families with observed parental genotypes set to missing). Default False.
    Returns:
        pg : :class:`snipar.gtarray`
            Return the polygenic score as a genotype array with columns: individual's PGS, mean of their siblings' PGS, observed/imputed paternal PGS,
            observed/imputed maternal PGS
    """
    # Check for SNP overlap
    if bedfile is not None:
        bed = Bed(bedfile, count_A1=True)
        snp_ids = bed.sid
    if bgenfile is not None:
        bgen = open_bgen(bgenfile)
        snp_ids = bgen.ids
        if np.unique(snp_ids).shape[0] == 1:
            snp_ids = bgen.rsids
    snp_set = set(snp_ids)
    in_snp_set = np.array([x in snp_set for x in pgs.snp_ids])
    if np.sum(in_snp_set)==0:
        print('No overlap between variants in weights file and observed genotypes')
        return None
    else:
        # Get genotype matrix
        bim = bedfile.split('.bed')[0] + '.bim'
        chromosome, sid, pos, alleles, obs_sid_index = get_snps(bed, bim, snp_ids=snp_ids)
        ids = bed.iid[:, 1]
        gts = bed[:, obs_sid_index].read().val
        G = gtarray(gts, ids, sid, alleles=alleles, pos=pos, chrom=chromosome, par_status=-1*np.ones((ids.shape[0], 2), dtype=int))
        # G = gtarray(gts, ids, sid, alleles=alleles, pos=pos, chrom=chromosome)
        # G = get_gts_matrix(bedfile=bedfile, bgenfile=bgenfile, par_gts_f=par_gts_f, ped=ped, snp_ids=pgs.snp_ids, sib=sib, compute_controls=compute_controls, verbose=verbose)
        # if sib:
        #     cols = np.array(['proband', 'sibling', 'paternal', 'maternal'])
        # else:
        #     cols = np.array(['proband', 'paternal', 'maternal'])
        cols = np.array(['proband'])
        if compute_controls:
            pgs_out = [pgs.compute(x,cols) for x in G[0:3]]
            if sib:
                o_cols = np.array(['proband', 'sibling', 'parental'])
            else:
                o_cols = np.array(['proband','parental'])
            pgs_out.append(pgs.compute(G[3], o_cols))
            return pgs_out
        else:
            x = pgs.compute(G,cols)
            return x


if __name__ == "__main__":
    args=parser.parse_args()
    main(args)