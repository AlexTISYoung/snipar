import h5py, argparse, code
import numpy as np
import numpy.ma as ma
from pysnptools.snpreader import Bed

def convert_str_array(x):
    x_shape = x.shape
    x = x.flatten()
    x_out = np.array([y.decode('UTF-8') for y in x])
    return x_out.reshape(x_shape)

def make_id_dict(x,col=0):
    if len(x.shape)>1:
        x = x[:,col]
    id_dict = {}
    for i in range(0,x.shape[0]):
        id_dict[x[i]] = i
    return id_dict

parser=argparse.ArgumentParser()
parser.add_argument('gts',type=str,help='Path to bed file with sibling genotypes')
parser.add_argument('outprefix',type=str,help='Prefix for output files')
parser.add_argument('--ped',type=str,help='Path to pedigree file with siblings sharing a family ID and non-siblings not', default = None)
parser.add_argument('--pargts', type=str, help='Path to HDF5 file with imputed parental genotypes from families with one parent genotyped',default = None)
parser.add_argument('--min_maf',type=float,help='Ignore SNPs with minor allele frequency below min_maf (default 0.01)',default=0.01)
parser.add_argument('--max_missing',type=float,help='Ignore SNPs with greater percent missing calls than max_missing (default 5)',default=1)
parser.add_argument('--hdf5_out',action='store_true',default=False)
args=parser.parse_args()

args = parser.parse_args()

gts_f = Bed(args.gts)
#gts_f = Bed('sim_reduced.bed')
gts_ids = gts_f.iid
gts_ids = gts_ids[:,0:2]
id_dict = {}
for i in range(0,gts_ids.shape[0]):
    id_dict[gts_ids[i,1]] =i

if args.pargts is None:
    ped = convert_str_array(np.loadtxt(args.ped, dtype='S20', skiprows=1))
    # Identify individuals by number of parents genotyped
    # Find indices
    genotype_indices = []
    bpg_ped_indices = []
    for i in range(0, ped.shape[0]):
        in_genotype = np.zeros((ped.shape[1]), dtype=bool)
        in_genotype[1:4] = np.array([x in id_dict for x in ped[i, 1:4]])
        if np.sum(in_genotype[1:4]) == 3:
            genotype_indices = genotype_indices + [id_dict[x] for x in ped[i, in_genotype]]
            bpg_ped_indices.append(i)

    if len(bpg_ped_indices)>0:
        print('Found '+str(len(bpg_ped_indices))+' with both parents genotyped')
        bpg_ped = ped[bpg_ped_indices,:]
    else:
        raise(ValueError('No individuals with both parents genotyped found'))

    print('Reading genotypes')
    gts = gts_f[genotype_indices,:].read().val
    gts_ids = gts_ids[genotype_indices, 0:2]
    id_dict = make_id_dict(gts_ids,1)
    gts = ma.array(gts, mask=np.isnan(gts), dtype=np.float32)

    # Calculate missingness and frequency
    freqs = ma.mean(gts, axis=0) / 2.0
    missingness = 100.0*ma.mean(gts.mask, axis=0)
    pass_filters = np.logical_and(freqs<(1-args.min_maf),np.logical_and(freqs>args.min_maf,missingness<args.max_missing))
    gts = gts[:,pass_filters]
    freqs = freqs[pass_filters]
    snp_ses = np.sqrt(2*freqs*(1-freqs))

    ### Form proband and parental genotype matrices ###
    bpg_ped = ped[bpg_ped_indices, :]
    N = bpg_ped.shape[0]
    G = np.zeros((N, 2, gts.shape[1]), dtype=np.float32)
    # BPG
    for i in range(0, bpg_ped.shape[0]):
        G[i, 0, :] = gts[id_dict[bpg_ped[i, 1]], :]
        G[i, 1, :] = gts[id_dict[bpg_ped[i, 2]], :] + gts[id_dict[bpg_ped[i, 3]], :]
    ped_out = bpg_ped
else:
    ped = np.loadtxt(args.ped, dtype='S20', skiprows=1)
    # Identify individuals by number of parents genotyped
    # Find indices
    genotype_indices = []
    bpg_ped_indices = []
    opg_ped_indices = []
    npg_ped_indices = []
    for i in range(0, ped.shape[0]):
        in_genotype = np.zeros((ped.shape[1]), dtype=bool)
        in_genotype[1:4] = np.array([x in id_dict for x in ped[i, 1:4]])
        if np.sum(in_genotype[1:4]) == 3:
            genotype_indices = genotype_indices + [id_dict[x] for x in ped[i, in_genotype]]
            bpg_ped_indices.append(i)
        elif in_genotype[1] and np.sum(in_genotype[2:4])==1:
            genotype_indices = genotype_indices + [id_dict[x] for x in ped[i, in_genotype]]
            opg_ped_indices.append(i)
        elif in_genotype[1]:
            genotype_indices = genotype_indices + [id_dict[x] for x in ped[i, in_genotype]]
            npg_ped_indices.append(i)

    bpg_ped = ped[bpg_ped_indices,:]
    opg_ped = ped[opg_ped_indices,:]

    ### Identify siblings without genotyped parents
    npg_ped = ped[npg_ped_indices, :]
    npg_ped_fams = np.unique(npg_ped[:, 0])
    sibships = {}
    n_sibs = 0
    for f in npg_ped_fams:
        pedf = npg_ped[npg_ped[:, 0] == f, :]
        parent_pairs = np.array([pedf[x, 2] + pedf[x, 3] for x in range(0, pedf.shape[0])])
        unique_parent_pairs = np.unique(parent_pairs)
        pcount = 0
        for par in unique_parent_pairs:
            pmatch = parent_pairs == par
            if np.sum(pmatch) > 1:
                sibs = pedf[pmatch, 1]
                sibships[f] = sibs
                genotype_indices = genotype_indices + [id_dict[x] for x in sibs]
                pcount += 1
                n_sibs += sibs.shape[0]
        if pcount > 1:
            print('More than one sibship without genotyped parents in family ' + str(
                f) + '. Implies incorrect/unsupported pedigree.')

    print('Found '+str(bpg_ped.shape[0])+' individuals with both parents genotyped')
    print('Found '+str(opg_ped.shape[0])+' individuals with one parent genotyped')
    print('Found '+str(len(sibships))+' sibships without genotyped parents, comprised of '+str(n_sibs)+' total siblings')

    genotype_indices = np.unique(np.array(genotype_indices))

    # Read genotypes
    gts = gts_f[genotype_indices,:].read().val
    gts = ma.array(gts, mask=np.isnan(gts), dtype=np.float32)

    # Calculate missingness and frequency
    freqs = ma.mean(gts, axis=0) / 2.0
    missingness = 100.0*ma.mean(gts.mask, axis=0)
    pass_filters = np.logical_and(freqs<(1-args.min_maf),np.logical_and(freqs>args.min_maf,missingness<args.max_missing))
    gts = gts[:,pass_filters]
    freqs = freqs[pass_filters]
    snp_ses = np.sqrt(2*freqs*(1-freqs))
    m = np.array(np.logical_not(gts.mask),dtype=np.int8)

    # Mean impute missing genotypes
    for j in range(0,gts.shape[1]):
        gts[gts.mask[:,j],j] = 2*freqs[j]

    # Rebuild ID dict
    gts_ids = gts_ids[genotype_indices, :]
    # Build dict
    id_dict = {}
    for i in xrange(0, gts_ids.shape[0]):
        id_dict[gts_ids[i, 1]] = i

    # Read genotypes imputed from families with one parent
    imp_par = h5py.File(args.imp_par,'r')
    #imp_par = h5py.File('impute_po.hdf5','r')
    imp_par_ped = np.array(imp_par['ped'])
    imp_par_gts = np.array(imp_par['imputed_par_gts'])
    imp_par_gts = imp_par_gts[:,pass_filters]
    imp_par_gts = ma.array(imp_par_gts,mask=np.isnan(imp_par_gts))
    # Mean impute missing genotypes
    for j in range(0,gts.shape[1]):
        imp_par_gts[imp_par_gts.mask[:,j],j] = 2*freqs[j]
    imp_par_dict = {}
    for i in range(0,imp_par_ped.shape[0]):
        imp_par_dict[imp_par_ped[i,1]] = i

    # Read genotypes imputed from sibships without genotyped parents
    imp_sib = h5py.File(args.imp_sib,'r')
    #imp_sib = h5py.File('impute_from_sibs.hdf5','r')
    imp_sib_fams = np.array(imp_sib['families'])
    imp_sib_gts = np.array(imp_sib['imputed_par_gts'])
    imp_sib_gts = imp_sib_gts[:,pass_filters]
    imp_sib_dict = {}
    for i in range(0,imp_sib_fams.shape[0]):
        imp_sib_dict[imp_sib_fams[i]] = i

    ### Form proband and parental genotype matrices ###
    N = bpg_ped.shape[0]+opg_ped.shape[0]+n_sibs
    G = np.zeros((N,2,gts.shape[1]),dtype=np.float32)
    ped_out = np.zeros((N,2),dtype='S20')
    all_information = np.ones((N),dtype=bool)
    # BPG
    for i in range(0,bpg_ped.shape[0]):
        G[i,0,:] = gts[id_dict[bpg_ped[i,1]],:]
        G[i,1,:] = gts[id_dict[bpg_ped[i,2]],:]+gts[id_dict[bpg_ped[i,3]],:]
        ped_out[i,:] = bpg_ped[i,0:2]

    # OPG
    for i in range(0,opg_ped.shape[0]):
        j = bpg_ped.shape[0]+i
        G[j, 0, :] = gts[id_dict[opg_ped[i, 1]], :]
        # Find imputed parental genotypes
        if opg_ped[i,1] in imp_par_dict:
            if opg_ped[i,2] in id_dict:
                G[j,1,:] = gts[id_dict[opg_ped[i,2]],:]+imp_par_gts[imp_par_dict[opg_ped[i,1]],:]
            else:
                G[j, 1, :] = gts[id_dict[opg_ped[i, 3]], :] + imp_par_gts[imp_par_dict[opg_ped[i, 1]], :]
        else:
            all_information[j] = False
        ped_out[j,:] = opg_ped[i,0:2]

    # Sib imputation
    sibship_fams = np.array(sibships.keys())
    G_start = j+1
    for fam in sibship_fams:
        sibs = sibships[fam]
        G_end = G_start+sibs.shape[0]
        # Check for imputed parental genotypes
        if fam in imp_sib_dict:
            sib_indices = np.array([id_dict[x] for x in sibs])
            G[G_start:G_end,0,:] = gts[sib_indices,:]
            G[G_start:G_end,1,:] = imp_sib_gts[imp_sib_dict[fam],:]
        else:
            all_information[G_start:G_end] = False
        ped_out[G_start:G_end,0] = fam
        ped_out[G_start:G_end,1] = sibs
        G_start = G_end

    del gts
    del imp_par_gts
    del imp_sib_gts

    if np.sum(all_information)<N:
        print('Removing '+str(N-np.sum(all_information))+' individuals with incomplete imputed parental genotypes')
        G = G[all_information,:,:]
        ped_out = ped_out[all_information,:]
        N = np.sum(all_information)

    print('Final sample of '+str(N)+' individuals')

### Compute relatedness matrices ###
print('Computing relatedness matrices')
print('Computing proband relatedness matrix')
# standardise columns
G[:,0,:] = (G[:,0,:]-np.mean(G[:,0,:],axis=0))/np.std(G[:,0,:],axis=0)
G[:,1,:] = np.sqrt(2)*(G[:,1,:]-np.mean(G[:,1,:],axis=0))/np.std(G[:,1,:],axis=0)
R = G[:,0,:].dot(G[:,0,:].T)/G.shape[2]

if args.hdf5_out:
    h5out = h5py.File(args.outprefix+'_grms.hdf5','w')
    h5out.create_dataset('R', R.shape, dtype='f', chunks=True, compression='gzip',
                            compression_opts=9)
    h5out['R'][:] = R
    h5out['N'] = G.shape[2]
else:
    R = R[np.tril_indices(R.shape[0])]
    R = R.reshape((1,R.shape[0]))
    R.tofile(args.outprefix+'_R.grm.bin')
    #R.tofile('R.grm.bin')
#del R

print('Computing parental relatedness matrix')
R_par = G[:,1,:].dot(G[:,1,:].T)/(2.0*G.shape[2])
if args.hdf5_out:
    h5out.create_dataset('R_par', R_par.shape, dtype='f', chunks=True, compression='gzip',
                         compression_opts=9)
    h5out['R_par'][:] = R_par
else:
    R_par = R_par[np.tril_indices(R_par.shape[0])]
    R_par = R_par.reshape((1,R_par.shape[0]))
    R_par.tofile(args.outprefix+'_R_par.grm.bin')
    #R_par.tofile('R_par.grm.bin')
#del R_par
print('Computing parent-offspring relatedness matrix')
R_o_par = G[:,0,:].dot(G[:,1,:].T)
R_o_par = (R_o_par + R_o_par.T)/(G.shape[2]*2.0)
#code.interact(local=locals())
#R_o_par = R_o_par/np.mean(np.diag(R_o_par))
if args.hdf5_out:
    h5out.create_dataset('R_o_par', R_o_par.shape, dtype='f', chunks=True, compression='gzip',
                         compression_opts=9)
    h5out['R_o_par'][:] = R_o_par
else:
    R_o_par = R_o_par[np.tril_indices(R_o_par.shape[0])]
    R_o_par = R_o_par.reshape((1,R_o_par.shape[0]))
    R_o_par.tofile(args.outprefix+'_R_o_par.grm.bin')
    #R_o_par.tofile('R_o_par.grm.bin')
del R_o_par
# Save pedigree
if args.hdf5_out:
    h5out['ped'] = ped_out
else:
    grm_out = np.array([args.outprefix+'_R',args.outprefix+'_R_par',args.outprefix+'_R_o_par'])
    np.savetxt(grm_out[0]+'.grm.id',ped_out[:,0:2],fmt='%s')
    np.savetxt(grm_out[1]+'.grm.id',ped_out[:,0:2],fmt='%s')
    np.savetxt(grm_out[2]+'.grm.id',ped_out[:,0:2],fmt='%s')
    np.savetxt('RDR_GRMs.txt',grm_out,fmt='%s')