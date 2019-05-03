from pysnptools.snpreader import Bed
import numpy as np

cols = tuple([x for x in xrange(3,103)])

weights = np.loadtxt('white_british_pc_weights.txt',usecols=cols)

snp_ids = np.loadtxt('white_british_pc_weights.txt',usecols=(0),dtype='S20')

chr = np.loadtxt('white_british_pc_weights.txt',usecols=(1),dtype=int)

pcs = np.zeros((51774,100))

for c in xrange(1,23):
    snps = snp_ids[chr==c]
    gts = Bed('../genotypes/chr_'+str(c)+'.bed')
    sid = gts.sid
    snp_dict = {}
    for s in xrange(0,sid.shape[0]):
        snp_dict[sid[s]] = s
    snp_indices = []
    for s in snps:
        snp_indices.append(snp_dict[s])
    snp_indices = np.array(snp_indices)
    g = gts[:,snp_indices].read()
    g = g.val
    g_notnan = np.logical_not(np.isnan(g))
    for j in xrange(0,g.shape[1]):
        g[:,j] = g[:,j]-np.mean(g[g_notnan[:,j],j])
        g[:,j] = g[:,j]/np.std(g[g_notnan[:,j],j])
        g[np.logical_not(g_notnan[:,j]),j] = 0
    pcs += g.dot(weights[chr==c,:])
    print(c)

np.savetxt('pcs.txt',pcs)