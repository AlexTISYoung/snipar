import numpy as np
import h5py, argparse

def simulate_ind(nsnp,f):
    return np.random.binomial(1,f,nsnp*2).reshape((nsnp,2))


def simulate_sibs(father,mother,n = 2):
    meioses = [[np.array(np.random.binomial(1, 0.5, father.shape[0])),np.array(np.random.binomial(1, 0.5, father.shape[0]))]
                for x in range(0,n)]
    ibd = np.zeros((n*(n-1)/2,father.shape[0]),dtype=np.int8)
    pcount = 0
    for i in range(1,n):
        for j in range(0,i):
            ibd[pcount,:] = np.array(meioses[i][0]==meioses[j][0],dtype=np.int8)+np.array(meioses[i][1]==meioses[j][1],dtype=np.int8)
            pcount+=1
    snp_index = np.array([x for x in range(0,father.shape[0])])
    gts = np.zeros((n,father.shape[0]),dtype=np.int8)
    for i in range(0,n):
        gts[i,:] = father[snp_index,meioses[i][0]]+mother[snp_index,meioses[i][1]]
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
parser.add_argument('--fsize',type=int,help='Number of children in each family (default 2)',default=2)
args=parser.parse_args()

nsnp = args.nsnp
f=args.f
nfam = args.nfam
outprefix = args.outprefix
n_sib_only = args.n_sib_only
n_one_parent = args.n_one_parent
p_sib_missing = args.p_sib_missing
fsize = args.fsize

father_gts = np.zeros((nfam,nsnp),dtype=np.int8)
mother_gts = np.zeros((nfam,nsnp),dtype=np.int8)
ibd = np.zeros((nfam,fsize*(fsize-1)/2,nsnp),dtype=np.int8)
sib_gts = np.zeros((nfam,fsize,nsnp),dtype=np.int8)

for i in range(0,nfam):
    father = simulate_ind(nsnp,f)
    father_gts[i,:] = np.sum(father,axis=1)
    mother = simulate_ind(nsnp,f)
    mother_gts[i, :] = np.sum(mother, axis=1)
    sibs = simulate_sibs(father,mother,fsize)
    sib_gts[i, :, :] = sibs[0]
    ibd[i,:,:] = sibs[1]

# Make pedigree file
fp = fsize+2
ped = np.zeros((nfam*fp,4),dtype='S20')
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
    ped[i * fp + fsize, 2] = str(i) + '_' + 'PP'
    ped[i * fp + fsize, 3] = str(i) + '_' + 'PM'
    ped[i * fp + fsize + 1, 2] = str(i) + '_' + 'MP'
    ped[i * fp + fsize + 1, 3] = str(i) + '_' + 'MM'

## Write pedigree file
ped = np.hstack((np.array(['FID','IID','FATHER_ID','MOTHER_ID'],dtype='S20'),ped))
np.savetxt(outprefix+'_fams.ped',ped,fmt='%s')

# Save in HDF5 file
hf = h5py.File(outprefix+'.hdf5','w')
hf['sib_gts'] = sib_gts
hf['ibd'] = ibd
hf['father_gts'] = father_gts
hf['mother_gts'] = mother_gts
hf['ibd_fams'] = np.array([x for x in range(0,nfam)])
hf['ped'] = ped
hf.close()

## Write full PED file
pout = open(outprefix+'.ped','w')
for i in range(0,nfam):
    for j in range(0,fsize):
        pout.write(pedline(ped[fp*i+j,:],sib_gts[i,j,:]))
    pout.write(pedline(ped[fp * i + fsize, :], father_gts[i, :]))
    pout.write(pedline(ped[fp * i + fsize+1, :], mother_gts[i, :]))

pout.close()

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

## Make MAP file
mout = file(outprefix+'.map','w')
for i in range(0,sib_gts.shape[2]):
    mout.write('1 rs'+str(i)+' '+str(i/float(100))+' '+str(i)+'\n')

mout.close()