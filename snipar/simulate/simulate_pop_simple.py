import numpy as np
import h5py, argparse, gzip
from snipar.utilities import *
from snipar.ibd import write_segs_from_matrix

def simulate_ind(nsnp,f):
    return np.random.binomial(1,f,nsnp*2).reshape((nsnp,2))

def simulate_sibs(father,mother, blocksize=200):
    # Compute blocks without recombination
    n_blocks = int(np.ceil(father.shape[0]/float(blocksize)))
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
args=parser.parse_args()

nsnp = args.nsnp
f=args.f
nfam = args.nfam
outprefix = args.outprefix
n_sib_only = args.n_sib_only
n_one_parent = args.n_one_parent
p_sib_missing = args.p_sib_missing
chrom = args.chrom
fsize = 2
if args.blocksize is None:
    blocksize = nsnp
else:
    blocksize = args.blocksize

father_gts = np.zeros((nfam,nsnp,2),dtype=np.int8)
mother_gts = np.zeros((nfam,nsnp,2),dtype=np.int8)
ibd = np.zeros((nfam,fsize*(fsize-1)//2,nsnp),dtype=np.int8)
sib_gts = np.zeros((nfam,fsize,nsnp,2),dtype=np.int8)

for i in range(0,nfam):
    father = simulate_ind(nsnp,f)
    father_gts[i,:,:] = father
    mother = simulate_ind(nsnp,f)
    mother_gts[i, :, :] = mother
    sibs = simulate_sibs(father,mother, blocksize = blocksize)
    sib_gts[i, :, :] = sibs[0]
    ibd[i,:,:] = sibs[1]

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
np.savetxt(outprefix+'.ped',ped_out,fmt='%s')

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

# Save segs in snipar format
sibpairs = np.column_stack(([str(x)+'_0' for x in range(nfam)],
                            [str(x)+'_1' for x in range(nfam)]))
segs = write_segs_from_matrix(ibd, sibpairs, [f'{chrom}_rs'+str(x) for x in range(nsnp)], 
                                [f"{chrom*nsnp+x}" for x in range(nsnp)], 
                                np.arange(0,1,1/float(nsnp)), chrom,
                                args.outprefix+'.chr_'+str(chrom)+'.ibd.segments.gz')

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
    remove.write(ped[fp * i + fsize, 1]+' '+ped[fp * i + fsize, 1]+'\n')
    remove.write(ped[fp * i + fsize+1, 1]+' '+ped[fp * i + fsize+1, 1]+'\n')
    if fsize>2:
        for j in range(2,fsize):
            remove_sib = np.random.binomial(1, p_sib_missing, 1)
            if remove_sib==1:
                remove.write(ped[fp * i + j, 1] + ' ' + ped[fp * i + j, 1]+'\n')

for i in range(n_sib_only, n_sib_only + n_one_parent):
    #remove_father = np.random.binomial(1, 0.5, 1)
    remove_father = 1
    if remove_father==1:
        remove.write(ped[fp * i + fsize, 1] + ' ' + ped[fp * i + fsize, 1]+'\n')
    else:
        remove.write(ped[fp * i + fsize+1, 1] + ' ' + ped[fp * i + fsize+1, 1]+'\n')
    # Remove sib
    for j in range(1, fsize):
        remove_sib = np.random.binomial(1, p_sib_missing, 1)
        if remove_sib==1:
            remove.write(ped[fp * i + j, 1] + ' ' + ped[fp * i + j, 1] + '\n')

for i in range(n_sib_only+n_one_parent,nfam):
    for j in range(1, fsize):
        remove_sib = np.random.binomial(1, p_sib_missing, 1)
        if remove_sib==1:
            remove.write(ped[fp * i + j, 1] + ' ' + ped[fp * i + j, 1] + '\n')

remove.close()

# Remove from KING Kinship file
remove = np.loadtxt(outprefix+'_remove.txt',dtype=str)
remove_dict = make_id_dict(remove,1)
kin = np.loadtxt(args.outprefix+'.king.kin0',dtype=str,delimiter='\t')
in_remove = np.logical_or(np.array([x in remove_dict for x in kin[:,1]]),
                          np.array([x in remove_dict for x in kin[:,3]]))
kin = kin[~in_remove,:]
np.savetxt(args.outprefix+'.king.kin0',kin,fmt='%s',delimiter='\t')

# Write map file
ppos = [f"{chrom*nsnp+x}" for x in range(nsnp)]
gpos = np.arange(0,1,1/float(nsnp))
map = np.vstack((np.array(['pposition','gposition']),
                 np.vstack((ppos,gpos)).T))
np.savetxt(args.outprefix+'.genetic_map.txt',map,fmt='%s')