#!/well/kong/users/wiw765/anaconda2/bin/python
import numpy as np
import numpy.ma as ma
from pysnptools.snpreader import Bed, Pheno
import h5py, argparse, code

def impute(g,pg,f):
    gcounts = np.array([np.sum(g==x) for x in range(0,3)])
    if pg == 2:
        if gcounts[0]>0:
            return np.nan
        elif gcounts[1] == 0:
            return (1+f*(2**g.shape[0]-1))/(1+f*(2**(g.shape[0]-1)-1))
        elif gcounts[1] == g.shape[0]:
            return f/(2**(g.shape[0]-1)-f*(2**(g.shape[0]-1)-1))
        else:
            return 1
    if pg == 1:
        if gcounts[0] == 0 and gcounts[2]>0:
            return (1 + f * (2 ** gcounts[2] - 1)) / (1 + f * (2 ** (gcounts[2] - 1) - 1))
        if gcounts[0]>0 and gcounts[2]>0:
            return 1
        if gcounts[2]==0 and gcounts[0]>0:
            e = g.shape[0]-gcounts[1]-1
            return f/(2**e-f*(2**e-1))
        if gcounts[1] == g.shape[0]:
            return 2*f
    if pg == 0:
        if gcounts[2]>0:
            return np.nan
        elif gcounts[1] == g.shape[0]:
            return (1+f*(2**g.shape[0]-1))/(1+f*(2**(g.shape[0]-1)-1))
        elif gcounts[1] == 0:
            return f/(2**(g.shape[0]-1)-f*(2**(g.shape[0]-1)-1))
        else:
            return 1

def simulate_sib(father,mother):
    return np.random.choice(father,1)+np.random.choice(mother,1)

def simulate_fam(n,f):
    father = np.random.binomial(1,f,2)
    mother = np.random.binomial(1,f,2)
    return [np.sum(father),np.sum(mother),np.array([simulate_sib(father,mother) for x in range(0,n)]).reshape((n))]

def test_impute(n,f):
    fam = simulate_fam(n,f)
    imputed = impute(fam[2],fam[0],f)
    return np.array([imputed,fam[1]])

######### Command line arguments #########
if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('gts',type=str,help='Path to bed file with parent and offspring genotypes')
    parser.add_argument('ped',type=str,help='Path to pedigree file')
    parser.add_argument('out',type=str,help='Prefix of HDF5 output file with imputed parental genotypes')
    args=parser.parse_args()

    ####################### Read in data #########################
### Read pedigree file ###
    ### Load pedigree
    ped = np.loadtxt(args.ped, dtype='S20', skiprows=1)
    #ped = np.loadtxt('relatedness/families.ped', dtype='S20', skiprows=1)

### Read genotypes ###
    #### Load genotypes
    gts_f = Bed(args.gts)
    #gts_f = Bed('genotypes/chr_22.bed')
    gts_ids = gts_f.iid
    # Build dict
    id_dict = {}
    for i in xrange(0, gts_ids.shape[0]):
        id_dict[gts_ids[i, 1]] = i

    # Find individuals with one parent genotyped
    genotyped = np.zeros((ped.shape[0],3),dtype=int)
    genotyped[:] = -1
    for i in xrange(0,ped.shape[0]):
        if ped[i,1] in id_dict:
            genotyped[i,0] = id_dict[ped[i,1]]
        if ped[i,2] in id_dict:
            genotyped[i,1] = id_dict[ped[i,2]]
        if ped[i,3] in id_dict:
            genotyped[i,2] = id_dict[ped[i,3]]

    pcount = np.array(genotyped[:,1]>=0,dtype=int)+np.array(genotyped[:,2]>=0,dtype=int)

    father_genotyped = np.logical_and(np.logical_and(genotyped[:,0]>=0,genotyped[:,1]>=0),genotyped[:,2]<0)
    mother_genotyped = np.logical_and(np.logical_and(genotyped[:,0]>=0,genotyped[:,2]>=0),genotyped[:,1]<0)
    one_parent_genotyped = np.logical_or(father_genotyped,mother_genotyped)

    index_vector = np.sort(np.unique(np.hstack((genotyped[one_parent_genotyped,0],
                                           genotyped[father_genotyped,1],
                                           genotyped[mother_genotyped,2]))))

    # Read genotypes
    print('Reading genotypes')
    gts = gts_f[index_vector, :].read().val
    pos = gts_f.pos[:, 2]
    sid = gts_f.sid
    gts = ma.array(gts,mask=np.isnan(gts),dtype=int)

    # rebuild ID dictionary
    gts_ids = gts_ids[index_vector]
    # Build dict
    id_dict = {}
    for i in xrange(0, gts_ids.shape[0]):
        id_dict[gts_ids[i, 1]] = i

    # Restricted pedigree
    ped = ped[one_parent_genotyped,:]
    ped_dict = {}
    for i in xrange(0,ped.shape[0]):
        ped_dict[ped[i,1]] = i
    # Get families
    fams = np.unique(ped[:,0])

    # Output array
    imputed_par_gts = np.zeros((ped.shape[0], gts.shape[1]), dtype=np.float32)
    imputed_par_gts[:] = np.nan

    freqs = ma.mean(gts,axis=0)/2.0

    father_genotyped = np.zeros((ped.shape[0]),dtype=bool)
    print('Imputing parental genotypes for '+str(fams.shape[0])+' families')
    for f in range(0,fams.shape[0]):
        fam = fams[f]
        pfam = ped[ped[:,0]==fam,:]
        # Identify sibships
        parent_pairs = np.array([pfam[x,2]+pfam[x,3] for x in range(0,pfam.shape[0])])
        unique_parent_pairs = np.unique(parent_pairs)
        for par in unique_parent_pairs:
            psib = pfam[parent_pairs==par,:]
            sib_indices = np.array([id_dict[x] for x in psib[:,1]])
            cgts = gts[sib_indices, :]
            sib_ped_indices = np.array([ped_dict[x] for x in psib[:,1]])
            if psib[0,2] in id_dict:
                pgts = gts[id_dict[psib[0,2]],:]
                father_genotyped[sib_ped_indices] = True
            elif psib[0,3] in id_dict:
                pgts = gts[id_dict[psib[0,3]],:]
            else:
                raise(ValueError('Missing parental genotype'))
            for j in range(0,gts.shape[1]):
                not_nan = np.logical_not(cgts.mask[:,j])
                if np.sum(not_nan)>0 and np.logical_not(pgts.mask[j]):
                    imputed_par_gts[sib_ped_indices,j] = impute(cgts[not_nan,j],pgts[j],freqs[j])

    par_gt_f = h5py.File(args.out+'.hdf5','w')
    par_gt_f.create_dataset('imputed_par_gts',imputed_par_gts.shape,dtype = 'f',chunks = True, compression = 'gzip', compression_opts=9)
    par_gt_f['imputed_par_gts'][:] = imputed_par_gts
    par_gt_f['father_genotyped'] = father_genotyped
    par_gt_f['ped'] = ped
    par_gt_f['pos'] = pos
    par_gt_f['sid'] = sid
    par_gt_f.close()
