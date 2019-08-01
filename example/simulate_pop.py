import numpy as np
import h5py, argparse

def simulate_ind(nsnp,f):
    return np.random.binomial(1,f,nsnp*2).reshape((nsnp,2))

def simulate_sib_gt(father_gt,mother_gt):
    return np.random.choice(father_gt,1)+np.random.choice(mother_gt,1)

def simulate_sibs(father,mother):
    meiosis_1_p = np.random.binomial(1, 0.5, father.shape[0])
    meiosis_1_m = np.random.binomial(1, 0.5, father.shape[0])
    meiosis_2_p = np.random.binomial(1, 0.5, father.shape[0])
    meiosis_2_m = np.random.binomial(1, 0.5, father.shape[0])
    ibd = np.array(meiosis_1_m==meiosis_2_m,dtype=np.int8)+np.array(meiosis_1_p==meiosis_2_p,dtype=np.int8)
    snp_index = np.array([x for x in range(0,father.shape[0])])
    gt_1 = father[snp_index,meiosis_1_p]+mother[snp_index,meiosis_1_m]
    gt_2 = father[snp_index, meiosis_2_p] + mother[snp_index, meiosis_2_m]
    return [gt_1,gt_2,ibd]

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
args=parser.parse_args()

nsnp = args.nsnp
f=args.f
nfam = args.nfam
outprefix = args.outprefx
n_sib_only = args.n_sib_only
n_one_parent = args.n_one_parent
f_sib_missing = args.p_sib_missing

father_gts = np.zeros((nfam,nsnp),dtype=np.int8)
mother_gts = np.zeros((nfam,nsnp),dtype=np.int8)
ibd = np.zeros((nfam,nsnp),dtype=np.int8)
sib1_gts = np.zeros((nfam,nsnp),dtype=np.int8)
sib2_gts = np.zeros((nfam,nsnp),dtype=np.int8)

for i in range(0,nfam):
    father = simulate_ind(nsnp,f)
    father_gts[i,:] = np.sum(father,axis=1)
    mother = simulate_ind(nsnp,f)
    mother_gts[i, :] = np.sum(mother, axis=1)
    sibs = simulate_sibs(father,mother)
    sib1_gts[i,:] = sibs[0]
    sib2_gts[i,:] = sibs[1]
    ibd[i,:] = sibs[2]

# Make pedigree file
ped = np.zeros((nfam*4,4),dtype='S20')
for i in range(0,nfam):
    # Fam ID
    ped[i*4:((i+1)*4),0] = str(i)
    # Individual ID
    ped[i * 4, 1] = str(i) + '_' + str(1)
    ped[i * 4+1, 1] = str(i) + '_' + str(2)
    ped[i * 4 + 2, 1] = str(i) + '_' + 'P'
    ped[i * 4 + 3, 1] = str(i) + '_' + 'M'
    # Father ID
    ped[i * 4, 2] = str(i) + '_' + 'P'
    ped[i * 4+1, 2] = str(i) + '_' + 'P'
    ped[i * 4 + 2, 2] = str(i) + '_' + 'PP'
    ped[i * 4 + 3, 2] = str(i) + '_' + 'MP'
    # Mother ID
    ped[i * 4, 3] = str(i) + '_' + 'M'
    ped[i * 4+1, 3] = str(i) + '_' + 'M'
    ped[i * 4 + 2, 3] = str(i) + '_' + 'PM'
    ped[i * 4 + 3, 3] = str(i) + '_' + 'MM'

## Write pedigree file
np.savetxt(outprefix+'_fams.ped',ped,fmt='%s')

# Save in HDF5 file
hf = h5py.File(outprefix+'.hdf5','w')
hf['sib1_gts'] = sib1_gts
hf['sib2_gts'] = sib2_gts
hf['ibd'] = ibd
hf['father_gts'] = father_gts
hf['mother_gts'] = mother_gts
p1 = np.array([4*x for x in range(0,nfam)])
p2 = np.array([4*x+1 for x in range(0,nfam)])
hf['ibd_sibs'] = np.hstack((ped[p1,1].reshape((nfam,1)),ped[p2,1].reshape((nfam,1))))
hf.close()

# # Simulate phenotype
# b = np.random.multivariate_normal(np.zeros((3)),np.array([[1,0.7,0.7],[0.7,1,0.7],[0.7,0.7,1]]),(sib1_gts.shape[1]))
# # additive genetic component
# a1 = sib1_gts.dot(b[:,0])
# a2 = sib2_gts.dot(b[:,0])
# a_par = (father_gts+mother_gts).dot(b[:,1])
# a1_sib = sib2_gts.dot(b[:,2])
# a2_sib = sib1_gts.dot(b[:,2])
#
# A1 = a1+a1_sib+a_par
# A2 = a2+a2_sib+a_par
# A = np.hstack((A1,A2))
# a_factor = np.sqrt(0.6)*np.power(np.std(A,axis=0),-1)
#
# b = b*a_factor
#
# e = np.random.randn((A.shape[0]))
# y = A*a_factor+np.sqrt(0.4)*e
#
#
# yout = np.hstack((np.vstack((ped[p1,0:2],ped[p2,0:2])),np.array(y,dtype=str).reshape((y.shape[0],1))))
# np.savetxt(outprefix+'_phenotype.ped',yout,fmt='%s')
#
# b_out = np.hstack((np.array(b[:,0],dtype=str).reshape(b.shape[0],1),
#                    np.array(b[:,1],dtype=str).reshape(b.shape[0],1),
#                    np.array(b[:,2],dtype=str).reshape(b.shape[0],1)))
# np.savetxt(outprefix+'.effects.txt',b_out,fmt='%s')

## Write full PED file
pout = open(outprefix+'.ped','w')
for i in range(0,nfam):
    pout.write(pedline(ped[4*i,:],sib1_gts[i,:]))
    pout.write(pedline(ped[4 * i+1, :], sib2_gts[i, :]))
    pout.write(pedline(ped[4 * i + 2, :], father_gts[i, :]))
    pout.write(pedline(ped[4 * i + 3, :], mother_gts[i, :]))

pout.close()

# Write list of individuals to remove
remove = open(outprefix+'_remove.txt','w')
for i in range(0,n_sib_only):
    # Remove parents
    remove.write(ped[4 * i + 2, 0]+' '+ped[4 * i + 2, 1]+'\n')
    remove.write(ped[4 * i + 3, 0]+' '+ped[4 * i + 3, 1]+'\n')

for i in range(n_sib_only, n_sib_only + n_one_parent):
    remove_father = np.random.binomial(1, 0.5, 1)
    if remove_father==1:
        remove.write(ped[4 * i + 2, 0] + ' ' + ped[4 * i + 2, 1]+'\n')
    else:
        remove.write(ped[4 * i + 3, 0] + ' ' + ped[4 * i + 3, 1]+'\n')
    # Remove sib
    remove_sib = np.random.binomial(1, f_sib_missing, 1)
    if remove_sib==1:
        remove.write(ped[4 * i + 1, 0] + ' ' + ped[4 * i + 1, 1]+'\n')

for i in range(n_sib_only+n_one_parent,nfam):
    remove_sib = np.random.binomial(1, f_sib_missing, 1)
    if remove_sib==1:
        remove.write(ped[4 * i + 1, 0] + ' ' + ped[4 * i + 1, 1]+'\n')

remove.close()

## Make MAP file
mout = file(outprefix+'.map','w')
for i in range(0,sib1_gts.shape[1]):
    mout.write('1 rs'+str(i)+' '+str(i/float(100))+' '+str(i)+'\n')

mout.close()