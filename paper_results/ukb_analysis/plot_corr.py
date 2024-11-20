import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.lines import Line2D

sib_corr = []
gain_imp = []
gain_noimp = []
gain_rob = []
names = []
for f in glob.glob('../output/gwas/v2/*/vc.npz'):
    b = np.load(f)
    f = f[:-6]
    name = f.split('/')[-2]
    # if name == 'Glucose': continue
        
    sibdiff = pd.read_csv(f'../output/gwas/sibdiff/with_grm/{name}/22.sumstats.gz', sep='\s+')
    robust = pd.read_csv(f'../output/gwas/robust/with_grm/{name}/22.sumstats.gz', sep='\s+')
    xx = robust.merge(sibdiff, on='SNP', how='inner', suffixes=['_rob', '_sib'])
    xx = xx.dropna()
    gain_rob.append(np.median(xx['direct_N_rob'] / xx['direct_N_sib']))

    noimp = pd.read_csv(f'../output/gwas/v2/paper_results/young/{name}/22.sumstats.gz', sep='\s+')
    imp = pd.read_csv(f'../output/gwas/v2/paper_results/unified/{name}/22.sumstats.gz', sep='\s+')
    sibdiff = pd.read_csv(f'../output/gwas/sibdiff/with_grm/europeans/{name}/22.sumstats.gz', sep='\s+')
    xx = imp.merge(sibdiff, on='SNP', how='inner', suffixes=['_imp', '_sibdiff'])
    xx = xx.dropna()
    gain_imp.append(np.median(xx['direct_N_imp'] / xx['direct_N_sibdiff']))
    yy = noimp.merge(sibdiff, on='SNP', how='inner', suffixes=['_noimp', '_sibdiff'])
    yy = yy.dropna()
    gain_noimp.append(np.median(yy['direct_N_noimp'] / yy['direct_N_sibdiff']))
    v = b['varcomps']
    sib_corr.append((v[0]/2+v[1])/sum(v))
    f = f[:-1]
    names.append(f[18:])

names = [x for _,x in sorted(zip(sib_corr, names))]
gain_imp = [x for _,x in sorted(zip(sib_corr, gain_imp))]
gain_noimp = [x for _,x in sorted(zip(sib_corr, gain_noimp))]
gain_rob = [x for _,x in sorted(zip(sib_corr, gain_rob))]
sib_corr = sorted(sib_corr)


#https://stackoverflow.com/questions/14542232/annotate-several-points-with-one-text-in-matplotlib
#textcoords='axes fraction' or textcoords='offset point'
types = ['o', 'v', 'h', '^', '<', 'x', '>', 's', '1', 'p', '*', 'd', 'P', 'X', 'D', 'H', '2', '+', '8']
matplotlib.rcParams.update({'font.size': 15})
plt.figure(figsize=(8,5))
plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.spines.top'] = False
for n, (i,j, k) in enumerate(zip(sib_corr, gain_imp, gain_noimp)):
    if types[n] in ['x', '1', '2', '+', 2]:
        plt.scatter(i, j, s=200, marker=types[n], facecolors='darkgoldenrod')
        plt.scatter(i, k, s=200, marker=types[n], facecolors='lightseagreen')
    else: 
        plt.scatter(i, j, s=200, marker=types[n], facecolors='none', edgecolors='darkgoldenrod')
        plt.scatter(i, k, s=200, marker=types[n], facecolors='none', edgecolors='lightseagreen')
# plt.scatter(sib_corr, gain_imp, label='unified vs sib-difference', c='darkgoldenrod', alpha=0.5, s=150, marker=types) #c=gain)c='darkgoldenrod', 
# plt.scatter(sib_corr, gain_noimp, label='Young et al. 2022 vs sib-difference', c='lightseagreen', alpha=0.5, s=150, marker=types) #c=gain)c='lightseagreen', 
plt.xlabel('phenotypic correlation of siblings')
plt.ylabel('relative effective sample size')
legend_elements = [Line2D([0], [0], color='darkgoldenrod', label='unified vs sib-difference', markersize=15),
                   Line2D([0], [0], color='lightseagreen', label='Young et al. 2022 vs sib-difference',
                          markerfacecolor='lightseagreen', markersize=15)]

# Create the figure
plt.legend(handles=legend_elements, loc='lower left')
plt.ylim(bottom=1.,)
plt.tight_layout()
# plt.savefig('empirical_gain_unified_young_vs_sibdiff_2.pdf')
plt.savefig('empirical_gain_unified_young_vs_sibdiff.pdf')





import pandas as pd
import numpy as np
import os
sib_corr = []
gain = []
name = []

for p in ['subjective.well.being', 'NC_M', 'BMI', 'self.rated.health', #'Glucose',
    'myopia', 'cigarettes.per.day','Non_HDL','AAFB','NC_F','EA',
    'DBP',
    'HDL','household.income','ever.smoked','Cognitive.ability','drinks.per.week','Neuroticism',
    'SBP','height']:
    sibdiff = pd.read_csv(f'../output/gwas/sibdiff/with_grm/{p}/22.sumstats.gz', sep='\s+')
    robust = pd.read_csv(f'../output/gwas/new_robust/with_grm/{p}/22.sumstats.gz', sep='\s+')
    xx = robust.merge(sibdiff, on='SNP', how='inner', suffixes=['_rob', '_sib'])
    xx = xx.dropna()
    gain.append(np.median(xx['direct_N_rob'] / xx['direct_N_sib']))
    name.append(p)
    f = f'../output/gwas/v2/{p}/vc.npz'
    try:
        b = np.load(f)
        v = b['varcomps']
        print(v)
        sib_corr.append((v[0]/2+v[1])/sum(v))
        print(p, (v[0]/2+v[1])/sum(v))
    except:
        f = f'../output/gwas/v2/{p}/vc_nogrm.npz'
        b= np.load(f)
        v = b['varcomps']
        sib_corr.append((v[0])/sum(v))
        print(p, (v[0])/sum(v))


name = [x for _,x in sorted(zip(sib_corr, name))]
gain = [x for _,x in sorted(zip(sib_corr, gain))]
sib_corr = sorted(sib_corr)


plt.figure(figsize=(8,5))
plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.spines.top'] = False
for n, (i,j) in enumerate(zip(sib_corr, gain)):
    if types[n] in ['x', '1', '2', '+', 2]:
        plt.scatter(i, j, s=200, marker=types[n], facecolors='mediumslateblue')
    else: 
        plt.scatter(i, j, s=200, marker=types[n], facecolors='none', edgecolors='mediumslateblue')
legend_elements = [Line2D([0], [0], color='mediumslateblue', label='robust vs sib-difference', markersize=15),]

# Create the figure
plt.legend(handles=legend_elements, loc='lower left')


plt.xlabel('phenotypic correlation of siblings')
plt.ylabel('relative effective sample size')
plt.ylim(bottom=1.)
plt.yticks([1, 1.1, 1.2, 1.3])
plt.tight_layout()
# plt.legend(loc='lower left')
plt.savefig('empirical_gain_linear_imputation_sibdiff_new_robust.pdf')




name_dict = {'subjective.well.being': 'subjective well-being',
             'NC_M': 'number of children (male)',
             'Neuroticism': 'neuroticism',
             'self.rated.health': 'self-rated health',
             'NC_F': 'number of children (female)',
             'DBP': 'diastolic blood pressure',
             'SBP': 'systolic blood pressure',
             'drinks.per.week': 'drinks-per-week',
             'Non_HDL': 'non-HDL cholesterol',
             'cigarettes.per.day': 'cigarettes-per-day',
             'ever.smoked': 'ever-smoked',
             'myopia': 'myopia',
             'household.income': 'household Income',
             'BMI': 'BMI',
             'HDL': 'hDL cholesterol',
             'AAFB':'age-at-first-birth (women)',
             'Cognitive.ability':'cognitive ability',
             'EA': 'educational attainment (years)',
             'height': 'height'
             }

names_ = [name_dict[i] for i in names]
legend_elements = [
                   Line2D([0], [0], marker=types[i], markeredgecolor='black', label=names_[i],linestyle = 'None',
                          markerfacecolor='none', markersize=15) if types[i] not in ['x', '1', '2', '+', 2] else Line2D([0], [0], marker=types[i], color='black', label=names_[i],linestyle = 'None', markersize=15)
                          for i in range(len(types))]
# figlegend = plt.figure()
figlegend = plt.figure(figsize=(4.4, 6))
figlegend.legend(handles=legend_elements, loc='center')
figlegend.tight_layout()
figlegend.savefig('legend.pdf')