import numpy as np
import h5py, argparse, gzip
from scipy.stats import multivariate_normal

def simulate_ind(nsnp,f):
    return np.random.binomial(1,f,nsnp*2).reshape((nsnp,2))

def simulate_sibs(father,mother, blocksize=200):
    # Compute blocks without recombination
    n_blocks = np.int(np.ceil(father.shape[0]/float(blocksize)))
    blocksizes = np.zeros((n_blocks),dtype=int)
    for i in range(n_blocks-1):
        blocksizes[i] = blocksize
    blocksizes[n_blocks-1] = father.shape[0]-np.sum(blocksizes)
    meioses = [[np.zeros((father.shape[0]),dtype=np.int8),np.zeros((father.shape[0]),dtype=np.int8)],
               [np.zeros((father.shape[0]),dtype=np.int8),np.zeros((father.shape[0]),dtype=np.int8)]]
    for i in range(2):
        for j in range(2):
            block_start=0
            for k in range(n_blocks):
                block_end = block_start+blocksizes[k]
                meioses[i][j][block_start:block_end] = np.random.binomial(1,0.5)
                block_start = block_end
    ibd = np.array(meioses[0][0]==meioses[1][0],dtype=np.int8)+np.array(meioses[0][1]==meioses[1][1],dtype=np.int8)
    gts = np.zeros((2,father.shape[0],2),dtype=np.int8)
    snp_index = np.array([x for x in range(0, father.shape[0])])
    for i in range(0,2):
        gts[i, :, 0] = father[snp_index,meioses[i][0]]
        gts[i, :, 1] = mother[snp_index, meioses[i][1]]
    return [gts,ibd]

def gt_convert(g):
    if g==0:
        return 'AA'
    elif g==1:
        return 'AG'
    elif g==2:
        return 'GG'

def pedline(ped,gts):
    pline = ped[0]+' '
    for i in range(1,4):
        pline += ped[i]+' '
    pline += '1 -9 '
    for i in range(0,gts.shape[0]-1):
        pline += gt_convert(gts[i])+' '
    pline += gt_convert(gts[i+1])+'\n'
    return pline

parser = argparse.ArgumentParser()
parser.add_argument('nsnp',type=int,help='Number of causal variants')
parser.add_argument('f',type=float,help='frequency of variants')
parser.add_argument('nfam',type=int,help='Number of families')
parser.add_argument('n_sib_only',type=int,help='Number of families to observe sibs only (no parental genotypes)')
parser.add_argument('n_one_parent',type=int,help='Number of families to observe only one genotyped parent')
parser.add_argument('p_sib_missing',type=float,help='Probability that one sibling is missing')
parser.add_argument('outprefix', type=str, help='Prefix of output ped files')
parser.add_argument('--blocksize',type=int,help='Size of blocks without recombination (number of SNPs)', default=None)
parser.add_argument('--chrom',type=int,help='Prefix for all snp ids', default="")
parser.add_argument('--gens',type=int,help='number of generations', default=2)
parser.add_argument('--am_noise',type=int,help='size of am noise', default=100000)
args=parser.parse_args()

nsnp = args.nsnp
f=args.f
nfam = args.nfam
outprefix = args.outprefix
n_sib_only = args.n_sib_only
n_one_parent = args.n_one_parent
p_sib_missing = args.p_sib_missing
chrom = args.chrom
gens = args.gens
am_noise = args.am_noise
fsize = 2
if args.blocksize is None:
    blocksize = nsnp
else:
    blocksize = args.blocksize
print(f"blocksize is {blocksize}")

ibd = np.zeros((nfam,fsize*(fsize-1)//2,nsnp),dtype=np.int8)
# direct_effects = np.random.normal(size=(nsnp,1))#np.zeros((nsnp,1))+1
# indirect_effects = np.random.normal(size=(nsnp,1))#np.zeros((nsnp,1))

direct_var = 1/100#1/nsnp
indirect_var = 0#0.25/nsnp
direct_indirect_corr = 0
direct_indirect_cov = direct_indirect_corr*np.sqrt(indirect_var*direct_var)
direct_indirect_cov_matrix = [[direct_var, direct_indirect_cov],[direct_indirect_cov, indirect_var]]
direct_indirect = multivariate_normal.rvs(cov = direct_indirect_cov_matrix, size=nsnp)
direct_effects = direct_indirect[:,0].reshape((-1, 1))
indirect_effects = direct_indirect[:,1].reshape((-1, 1))
print("direct_effects", direct_effects)
print("indirect_effects", indirect_effects)
print("direct_indirect_cov_matrix", direct_indirect_cov_matrix)

v_direct = []
v_indirect = []
cov_parental_direct = []
spousal_corrs = []
sib_corrs = []

father_gts = np.random.binomial(1, f, size=(nfam,nsnp,2))
mother_gts = np.random.binomial(1, f, size=(nfam,nsnp,2))
sib_gts = np.zeros((nfam,fsize,nsnp,2))
# second fpgs can be replaced with normal joint regression 
father_sum = np.sum(father_gts, axis=2)
father_indirect = np.zeros(nfam)
father_direct = (father_sum@direct_effects).flatten()
father_noise = np.random.normal(size=father_direct.shape)
father_phen = father_direct+father_indirect+father_noise

mother_sum = np.sum(mother_gts, axis=2)
mother_indirect = np.zeros(nfam)
mother_direct = (mother_sum@direct_effects).flatten()
mother_noise = np.random.normal(size=mother_direct.shape)
mother_phen = mother_direct+mother_indirect

v_direct.append(np.var(np.hstack((father_direct, mother_direct))))
v_indirect.append(np.var(np.hstack((father_indirect, mother_indirect))))
cov_parental_direct.append(np.cov(father_direct, mother_direct)[0,1])
spousal_corrs.append(-1)
sib_corrs.append(-1)
for gen in range(gens):
    print("*"*100)
    print(f"gen is {gen}")
    
    father_indexes = np.argsort(father_phen + np.random.normal(scale=am_noise*np.std(father_phen), size=father_phen.shape))
    father_direct = father_direct[father_indexes]
    father_indirect = father_indirect[father_indexes]
    father_phen = father_phen[father_indexes]
    father_gts = father_gts[father_indexes,:,:]
    father_sum = father_sum[father_indexes,:]
    
    mother_indexes = np.argsort(mother_phen + np.random.normal(scale=am_noise*np.std(mother_phen), size=mother_phen.shape))
    mother_direct = mother_direct[mother_indexes]
    mother_indirect = mother_indirect[mother_indexes]
    mother_phen = mother_phen[mother_indexes]
    mother_gts = mother_gts[mother_indexes,:,:]
    mother_sum = mother_sum[mother_indexes,:]
    print(f"assort {np.corrcoef((mother_phen, father_phen))[0,1]}")
    for i in range(0,nfam):
        father = father_gts[i, :, :]
        mother = mother_gts[i, :, :]
        sibs = simulate_sibs(father,mother, blocksize = blocksize)
        sib_gts[i, :, :, :] = sibs[0]
        ibd[i,:,:] = sibs[1]

    tmp=((father_sum+mother_sum)@indirect_effects).flatten()

    sib_indirect = np.tile(tmp, (2,1)).T
    sib_direct = np.zeros((nfam, fsize))
    sib_phen = np.zeros((nfam, fsize))
    for i in range(fsize):
        sib_direct[:, i] = (np.sum(sib_gts[:, i, :, :], axis=2)@direct_effects).flatten()
        sib_noise = np.random.normal(size=sib_direct[:, i].shape)
        sib_phen[:, i] = sib_indirect[:, i]+sib_direct[:, i]+sib_noise
    
    v_direct.append(np.var(sib_direct))
    v_indirect.append(np.var(sib_indirect))
    cov_parental_direct.append(np.cov(father_direct, mother_direct)[0,1])
    sib_corr = np.corrcoef(sib_phen[:, 0].flatten(), sib_phen[:, 1].flatten())[0, 1]
    sib_corrs.append(sib_corr)
    print("sib corr", sib_corr)
    spousal_corr = np.corrcoef(father_phen, mother_phen)[0,1]
    spousal_corrs.append(spousal_corr)
    print("spousal corr", spousal_corr)
    print("v_direct", v_direct[-1])
    if gen < (gens-1):
        whos_father = np.array([True,False]*(nfam//2))
        np.random.shuffle(whos_father)
        father_gts = np.concatenate([sib_gts[whos_father, i, :, :] for i in range(sib_gts.shape[1])], 0)
        father_sum = np.sum(father_gts, axis=2)
        father_direct = np.concatenate([sib_direct[whos_father, i] for i in range(sib_direct.shape[1])], 0)
        father_indirect = np.concatenate([sib_indirect[whos_father, i] for i in range(sib_indirect.shape[1])], 0)
        father_phen = np.concatenate([sib_phen[whos_father, i] for i in range(sib_phen.shape[1])], 0)
        
        mother_gts = np.concatenate([sib_gts[~whos_father, i, :, :] for i in range(sib_gts.shape[1])], 0)
        mother_sum = np.sum(mother_gts, axis=2)
        mother_direct = np.concatenate([sib_direct[~whos_father, i] for i in range(sib_direct.shape[1])], 0)
        mother_indirect = np.concatenate([sib_indirect[~whos_father, i] for i in range(sib_indirect.shape[1])], 0)
        mother_phen = np.concatenate([sib_phen[~whos_father, i] for i in range(sib_phen.shape[1])], 0)
v_direct = np.array(v_direct)
v_indirect = np.array(v_indirect)
cov_parental_direct = np.array(cov_parental_direct)
spousal_corrs = np.array(spousal_corrs)
sib_corrs = np.array(sib_corrs)

def get_rel_change(x):
    return list(range(len(x)-2)), (x[2:]-x[1:-1])/x[1:-1]

from matplotlib import pyplot as plt
plt.rcParams["figure.figsize"] = (20,10)
plt.plot(*get_rel_change(spousal_corrs), label="spousal_corrs")
plt.plot(*get_rel_change(sib_corrs), label="sib_corrs")
plt.plot(*get_rel_change(v_direct), label="v_direct")
plt.axhline()
plt.grid()
plt.title(f"relative change nsnp {nsnp}, nfam {nfam}, blocksize {blocksize}")
plt.legend()
plt.savefig(f"{args.outprefix}rel_change_nsnp{nsnp}_nfam{nfam}_blocksize{blocksize}")


np.save(f"{args.outprefix}phen", {
    "v_direct":v_direct,
    "v_indirect":v_indirect,
    "cov_parental_direct":cov_parental_direct,
    "spousal_corrs":spousal_corrs,
    "sib_corrs":sib_corrs
})

# Make pedigree file
fp = fsize+2
ped = np.zeros((nfam*fp,4),dtype="U20")
for i in range(0,nfam):
    # Fam ID
    ped[i*fp:((i+1)*fp),0] = str(i)
    # Individual ID
    for j in range(0,fsize):
        ped[i*fp+j,1] = str(i)+'_'+str(j)
        ped[i*fp+j,2] = str(i) + '_' + 'P'
        ped[i * fp + j, 3] = str(i) + '_' + 'M'
    ped[i * fp + j+1, 1] = str(i) + '_' + 'P'
    ped[i * fp + j+2, 1] = str(i) + '_' + 'M'
    ped[i * fp + fsize, 2] = str(i) + '_' + 'P_P'
    ped[i * fp + fsize, 3] = str(i) + '_' + 'P_M'
    ped[i * fp + fsize + 1, 2] = str(i) + '_' + 'M_P'
    ped[i * fp + fsize + 1, 3] = str(i) + '_' + 'M_M'

ped_out = np.vstack((np.array(['FID','IID','FATHER_ID','MOTHER_ID'],dtype='S20').reshape((1,4)),ped))
## Write pedigree file
np.savetxt(outprefix+'_fams.ped',ped_out,fmt='%s')

# Write relationships and age sex king output
kin_out = open(args.outprefix+'.king.kin0','w')
kin_out.write('FID1\tID1\tFID2\tID2\tN_SNP\tHetHet\tIBS0\tHetConc\tHomIBS0\tKinship\tIBD1Seg\tIBD2Seg\tPropIBD\tInfType\n')
age_sex_out = open(args.outprefix+'.agesex','w')
age_sex_out.write('FID IID age sex\n')
for i in range(0, nfam):
    # Write king.kin0
    kin_out.write(str(i)+'\t'+str(i)+'_0'+'\t'+str(i)+'\t'+str(i)+'_1'+'\t0\t0\t0\t0\t0\t0\t0\t0\t0\tFS\n')
    kin_out.write(
        str(i) + '\t' + str(i) + '_0' + '\t' + str(i) + '\t' + str(i) + '_P' + '\t0\t0\t0\t0\t0\t0\t0\t0\t0\tPO\n')
    kin_out.write(
        str(i) + '\t' + str(i) + '_0' + '\t' + str(i) + '\t' + str(i) + '_M' + '\t0\t0\t0\t0\t0\t0\t0\t0\t0\tPO\n')
    kin_out.write(
        str(i) + '\t' + str(i) + '_1' + '\t' + str(i) + '\t' + str(i) + '_P' + '\t0\t0\t0\t0\t0\t0\t0\t0\t0\tPO\n')
    if i == nfam-1:
        kin_out.write(
            str(i) + '\t' + str(i) + '_1' + '\t' + str(i) + '\t' + str(i) + '_M' + '\t0\t0\t0\t0\t0\t0\t0\t0\t0\tPO')
    else:
        kin_out.write(
            str(i) + '\t' + str(i) + '_1' + '\t' + str(i) + '\t' + str(i) + '_M' + '\t0\t0\t0\t0\t0\t0\t0\t0\t0\tPO\n')
    # Write age_sex
    age_sex_out.write(str(i)+' '+str(i) + '_0 20 M\n')
    age_sex_out.write(str(i) + ' ' + str(i) + '_1 20 F\n')
    age_sex_out.write(str(i) + ' ' + str(i) + '_M 40 F\n')
    if i == nfam-1:
        age_sex_out.write(str(i) + ' ' + str(i) + '_P 40\tM')
    else:
        age_sex_out.write(str(i) + ' ' + str(i) + '_P 40 M\n')

age_sex_out.close()
kin_out.close()

ibd = ibd[:,0,:]


# Convert IBD to our input format
our_ibd_out = gzip.open(args.outprefix+'.our.segments.gz','wb')
our_ibd_out.write(b'ID1\tID2\tIBDType\tChr\tstart_coordinate\tstop_coordinate\n')

for i in range(0,nfam):
    start_snps = []
    end_snps = []
    ibd_type = []
    init_ibd = ibd[i,0]
    start_snps.append(0)
    ibd_type.append(ibd[i,1])
    for j in range(1,ibd.shape[1]):
        if not ibd[i,j]==init_ibd:
            end_snps.append(j-1)
            init_ibd = ibd[i,j]
            start_snps.append(j)
            ibd_type.append(ibd[i,j])
    end_snps.append(ibd.shape[1]-1)
    # Write to file
    nseg = len(start_snps)
    if nseg>0:
        for s in range(0,nseg):
            t = f"{i}_0\t{i}_1\t{ibd_type[s]}\t{chrom}\t{chrom*nsnp+start_snps[s]}\t{chrom*nsnp+end_snps[s]}"
            our_ibd_out.write(t.encode("ascii"))
            if s==(nseg-1):
                our_ibd_out.write(b'\n')
our_ibd_out.close()

# Determine IBD segments assessed
allsegs_out = open(args.outprefix+'.kingallsegs.txt','w')
allsegs_out.write('Segment\tChr\tStartMB\tStopMB\tLength\tN_SNP\tStartSNP\tStopSNP\n')
# Ignore first and last SNP
insegs = np.ones((nsnp),dtype=bool)
insegs[[0,insegs.shape[0]-1]] = False
allsegs_out.write(f'1\t{chrom}\t0\t1\t1\t'+str(nsnp-2)+f'\t{chrom}_rs1\t{chrom}_rs'+str(nsnp-2)+'\n')
allsegs_out.close()

# Convert IBD to KING format
king_out = gzip.open(args.outprefix+'.king.segments.gz','wb')
king_out.write(b'FID1\tID1\tFID2\tID2\tIBDType\tChr\tStartMB\tStopMB\tStartSNP\tStopSNP\tN_SNP\tLength\n')

for i in range(0,nfam):
    start_snps = []
    end_snps = []
    ibd_type = []
    init_ibd = ibd[i,1]
    if init_ibd>0:
        start_snps.append(f'{chrom}_rs'+str(1))
        ibd_type.append(ibd[i,1])
    for j in range(2,ibd.shape[1]-1):
        if not ibd[i,j]==init_ibd:
            if init_ibd>0:
                end_snps.append(f'{chrom}_rs'+str(j-1))
            init_ibd = ibd[i,j]
            if init_ibd>0:
                start_snps.append(f'{chrom}_rs'+str(j))
                ibd_type.append(ibd[i,j])
    if init_ibd>0:
        end_snps.append(f'{chrom}_rs'+str(j))
    # Write to file
    nseg = len(start_snps)
    if nseg>0:
        for s in range(0,nseg):
            t = str(i)+'\t'+str(i)+'_0\t'+str(i)+'\t'+str(i)+'_1\tIBD'+str(ibd_type[s])+f'\t{chrom}\t0.0\t0.0\t'+start_snps[s]+'\t'+end_snps[s]+'\t0\t0\n'
            king_out.write(t.encode("ascii"))
king_out.close()


## Write full haps file
haps_out = open(outprefix+'.haps','w')
haps_matrix = np.zeros((nsnp,2*(2+fsize)*nfam),dtype=np.int8)
nhaps = 0
for i in range(0,nfam):
    # Sib haplotypes
    for j in range(0,fsize):
        haps_matrix[:,nhaps] = sib_gts[i,j,:,0]
        nhaps += 1
        haps_matrix[:,nhaps] = sib_gts[i,j,:,1]
        nhaps += 1
    # Paternal haplotypes
    haps_matrix[:,nhaps] = father_gts[i,:,0]
    nhaps += 1
    haps_matrix[:,nhaps] = father_gts[i, :, 1]
    nhaps += 1
    # maternal haplotypes
    haps_matrix[:,nhaps] = mother_gts[i,:,0]
    nhaps += 1
    haps_matrix[:,nhaps] = mother_gts[i, :, 1]
    nhaps += 1

haps_snps = np.column_stack(([f"{chrom}"]*nsnp,[f'{chrom}_rs'+str(x) for x in range(nsnp)],
                             [f"{chrom*nsnp+x}" for x in range(nsnp)],['A' for x in range(nsnp)],['G' for x in range(nsnp)]))

haps_out = np.column_stack((haps_snps,haps_matrix))
np.savetxt(args.outprefix+'.haps',haps_out,fmt='%s')

# Write sample file
sample_out = open(outprefix+'.sample','w')
sample_out.write('ID_1 ID_2 missing\n')
sample_out.write('0 0 0\n')
for i in range(0,nfam):
    for j in range(0,fsize):
        ID = str(i)+'_'+str(j)
        sample_out.write(ID+' '+ID+' 0\n')
    father_ID = str(i)+'_P'
    sample_out.write(father_ID+' '+father_ID+' 0\n')
    mother_ID = str(i)+'_M'
    if i < (nfam-1):
        sample_out.write(mother_ID + ' ' + mother_ID + ' 0\n')
    else:
        sample_out.write(mother_ID + ' ' + mother_ID + ' 0')

sample_out.close()

# Write list of individuals to remove
remove = open(outprefix+'_remove.txt','w')
for i in range(0,n_sib_only):
    # Remove parents
    remove.write(ped[fp * i + fsize, 0]+' '+ped[fp * i + fsize, 1]+'\n')
    remove.write(ped[fp * i + fsize+1, 0]+' '+ped[fp * i + fsize+1, 1]+'\n')
    if fsize>2:
        for j in range(2,fsize):
            remove_sib = np.random.binomial(1, p_sib_missing, 1)
            if remove_sib==1:
                remove.write(ped[fp * i + j, 0] + ' ' + ped[fp * i + j, 1]+'\n')

for i in range(n_sib_only, n_sib_only + n_one_parent):
    remove_father = np.random.binomial(1, 0.5, 1)
    if remove_father==1:
        remove.write(ped[fp * i + fsize, 0] + ' ' + ped[fp * i + fsize, 1]+'\n')
    else:
        remove.write(ped[fp * i + fsize+1, 0] + ' ' + ped[fp * i + fsize+1, 1]+'\n')
    # Remove sib
    for j in range(1, fsize):
        remove_sib = np.random.binomial(1, p_sib_missing, 1)
        if remove_sib==1:
            remove.write(ped[fp * i + j, 0] + ' ' + ped[fp * i + j, 1] + '\n')

for i in range(n_sib_only+n_one_parent,nfam):
    for j in range(1, fsize):
        remove_sib = np.random.binomial(1, p_sib_missing, 1)
        if remove_sib==1:
            remove.write(ped[fp * i + j, 0] + ' ' + ped[fp * i + j, 1] + '\n')

remove.close()

np.savetxt(args.outprefix+'.direct_effects',direct_effects)
np.savetxt(args.outprefix+'.indirect_effects',indirect_effects)
np.savetxt(args.outprefix+'.father_phen',father_phen)
np.savetxt(args.outprefix+'.mother_phen',mother_phen)
np.savetxt(args.outprefix+'.sib_phen',sib_phen.flatten())