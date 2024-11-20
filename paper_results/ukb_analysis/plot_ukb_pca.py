# PCA
import matplotlib.pyplot as plt
import pandas as pd
f = '/disk/genetics/ukb/jguan/ukb_analysis/output'
PCs = pd.read_csv(f + '/ukb_sqc_v2_combined_header.txt', sep='\s+')
rob = pd.read_csv(f + '/robust_ids.txt', sep='\s+')
sib = pd.read_csv(f + '/sibdiff_ids.txt', sep='\s+')


plt.figure(figsize=(10,5))
# plt.figure()
plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.spines.top'] = False
plt.subplot(121)
plt.scatter(rob[rob['in.white.British.ancestry.subset'] == 0]['PC1'], rob[rob['in.white.British.ancestry.subset'] == 0]['PC2'], color='grey', alpha=0.5, label='Non-\'white British\'') # (11.83%)')
plt.scatter(rob[rob['in.white.British.ancestry.subset'] == 1]['PC1'], rob[rob['in.white.British.ancestry.subset'] == 1]['PC2'], color='red', alpha=0.5, label='\'White British\'') # (88.17%)')
plt.legend(loc='upper left')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.title('robust')

plt.subplot(122)
plt.scatter(sib[sib['in.white.British.ancestry.subset'] == 0]['PC1'], sib[sib['in.white.British.ancestry.subset'] == 0]['PC2'], color='grey', alpha=0.5, label='Non-\'white British\'') # (11.84%)')
plt.scatter(sib[sib['in.white.British.ancestry.subset'] == 1]['PC1'], sib[sib['in.white.British.ancestry.subset'] == 1]['PC2'], color='red', alpha=0.5, label='\'White British\'') # (88.16%)')
plt.xlabel('PC1')

plt.title('sib-difference')
plt.tight_layout()
plt.savefig('PC_full_ukb.pdf')
plt.savefig('PC_full_ukb.png')