import numpy as np
from sklearn.model_selection import LeaveOneOut
import pandas as pd


def plot(alpha, alpha_cov, alpha_max, alpha_cov_max, non_sampling=False):
    s = np.zeros(alpha.shape[0])
    loo = LeaveOneOut()
    loo.get_n_splits(alpha)
    for i, (index, _) in enumerate(loo.split(alpha)):
        if alpha_cov.ndim == 1:
            if not non_sampling:
                s[i] = np.var(alpha[index, 0] / np.sqrt(alpha_cov)[index]) / (np.var(alpha_max[index, 0] / np.sqrt(alpha_cov_max)[index,0,0]))
            else:
                s[i] = np.mean(alpha[index,0] ** 2 - alpha_cov[index]) / (np.mean(alpha_max[index,0] ** 2) - np.mean(alpha_cov_max[index,0,0]))
                # s[i] = (np.mean(alpha[index,0] ** 2) - np.mean(alpha_cov[index])) / (np.mean(alpha_max[index,0] ** 2) - np.mean(alpha_cov_max[index,0,0]))
        else:
            if not non_sampling:
                s[i] = np.var(alpha[index, 0] / np.sqrt(np.diagonal(alpha_cov, axis1=1, axis2=2))[index,0]) / (np.var(alpha_max[index, 0] / np.sqrt(alpha_cov_max)[index,0,0]))
            else:
                s[i] = (np.mean(alpha[index,0] ** 2) - np.mean(alpha_cov[index,0,0])) /(np.mean(alpha_max[index,0] ** 2) - np.mean(alpha_cov_max[index,0,0]))
    mean = np.mean(s)
    std = np.sqrt(np.sum(np.power(mean - s, 2)) * (s.shape[0] - 1) / s.shape[0])
    if alpha_cov.ndim == 1:
        return mean, std, alpha_cov.mean()
    else:
        return mean, std, alpha_cov[:, 0, 0].mean()


Zvar_std = []
Zvar_sib = []
Zvar_lin = []
Zvar_robust = []
Fst = [0, 0.001, 0.01, 0.1]
df = []
print('Non sampling var allsibs')
npz_ref = np.load(f'sp_pairs/allsibs/Fst_0.001_h2pop_0.5_h2quad_0.0_standard.npz')
for F in [0, 0.001, 0.01, 0.1]:
    print('\n', F)
    if F == 0:
        tmp = 'tmp/'
    else:
        tmp = ''
    print('robust')
    npz = np.load(f'sp_pairs/allsibs/{tmp}Fst_{F}_h2pop_0.5_h2quad_0.0_robust.npz')
    m, std, var = plot(npz['alpha'], npz['alpha_cov'], npz_ref['alpha'], npz_ref['alpha_cov'],non_sampling=True)
    df.append(pd.DataFrame({'Fst': [F], 
                            'method': ['robust'],
                            'non_sampling_var': m, 'non_sampling_var_std': std, 'var': var}))

    print('young')
    npz = np.load(f'sp_pairs/allsibs/{tmp}Fst_{F}_h2pop_0.5_h2quad_0.0.npz')
    m, std, var = plot(npz['alpha'], npz['alpha_cov'], npz_ref['alpha'], npz_ref['alpha_cov'],non_sampling=True)
    df.append(pd.DataFrame({'Fst': [F], 
                            'method': ['Young et al. 2022/unified'],
                            'non_sampling_var': m, 'non_sampling_var_std': std, 'var': var}))

    print('standard')
    npz = np.load(f'sp_pairs/allsibs/{tmp}Fst_{F}_h2pop_0.5_h2quad_0.0_standard.npz')
    m, std, var = plot(npz['alpha'], npz['alpha_cov'], npz_ref['alpha'], npz_ref['alpha_cov'],non_sampling=True)
    df.append(pd.DataFrame({'Fst': [F], 
                            'method': ['standard'],
                            'non_sampling_var': m, 'non_sampling_var_std': std, 'var': var}))

    print('sibdiff')
    npz = np.load(f'sp_pairs/allsibs/{tmp}Fst_{F}_h2pop_0.5_h2quad_0.0_sibdiff.npz')
    m, std, var = plot(npz['alpha'], npz['alpha_cov'].squeeze(), npz_ref['alpha'], npz_ref['alpha_cov'],non_sampling=True)
    df.append(pd.DataFrame({'Fst': [F], 
                            'method': ['sib diff'],
                            'non_sampling_var': m, 'non_sampling_var_std': std, 'var': var}))
    
    
    
print(Zvar_std)
print(Zvar_sib)
print(Zvar_lin)
print(Zvar_robust)
pd.concat(df).to_csv('paper_materials/sp_non_sampling_var_ref_to_standardGWAS_at_Fst0.001_allsibs.csv', index=False, sep='\t')



print('Z var allsibs')
for F in [0, 0.001, 0.01, 0.1]:
    print('\n', F)
    if F == 0:
        tmp = 'tmp/'
    else:
        tmp = ''
    print('robust')
    npz = np.load(f'sp_pairs/allsibs/{tmp}Fst_{F}_h2pop_0.5_h2quad_0.0_robust.npz')
    m, std, var = plot(npz['alpha'], npz['alpha_cov'], npz_ref['alpha'], npz_ref['alpha_cov'])
    df.append(pd.DataFrame({'Fst': [F], 
                            'method': ['robust'],
                            'non_sampling_var': m, 'non_sampling_var_std': std, 'var': var}))

    print('young')
    npz = np.load(f'sp_pairs/allsibs/{tmp}Fst_{F}_h2pop_0.5_h2quad_0.0.npz')
    m, std, var = plot(npz['alpha'], npz['alpha_cov'], npz_ref['alpha'], npz_ref['alpha_cov'])
    df.append(pd.DataFrame({'Fst': [F], 
                            'method': ['Young et al. 2022/unified'],
                            'non_sampling_var': m, 'non_sampling_var_std': std, 'var': var}))

    print('standard')
    npz = np.load(f'sp_pairs/allsibs/{tmp}Fst_{F}_h2pop_0.5_h2quad_0.0_standard.npz')
    # Zvar_std.append(plot(standard['alpha'], standard['alpha_cov']))
    m, std, var = plot(npz['alpha'], npz['alpha_cov'], npz_ref['alpha'], npz_ref['alpha_cov'])
    df.append(pd.DataFrame({'Fst': [F], 
                            'method': ['standard'],
                            'non_sampling_var': m, 'non_sampling_var_std': std, 'var': var}))

    print('sibdiff')
    npz = np.load(f'sp_pairs/allsibs/{tmp}Fst_{F}_h2pop_0.5_h2quad_0.0_sibdiff.npz')
    m, std, var = plot(npz['alpha'], npz['alpha_cov'].squeeze(), npz_ref['alpha'], npz_ref['alpha_cov'])
    df.append(pd.DataFrame({'Fst': [F], 
                            'method': ['sib diff'],
                            'non_sampling_var': m, 'non_sampling_var_std': std, 'var': var}))
    
pd.concat(df).to_csv('paper_materials/sp_Z_var_ref_to_standardGWAS_at_Fst0.001_allsibs.csv', index=False, sep='\t')


df = []
print('Non sampling var with unrel')
npz_ref = np.load(f'sp_pairs/Fst_0.001_h2pop_0.5_h2quad_0.0_standard_with_unrel.npz')
for F in [0, 0.001, 0.01, 0.1]:
    print('\n', F)
    print('robust')
    npz = np.load(f'sp_pairs/Fst_{F}_h2pop_0.5_h2quad_0.0_robust.npz')
    m, std, var = plot(npz['alpha'], npz['alpha_cov'], npz_ref['alpha'], npz_ref['alpha_cov'],non_sampling=True)
    df.append(pd.DataFrame({'Fst': [F], 
                            'method': ['robust'],
                            'non_sampling_var': m, 'non_sampling_var_std': std, 'var': var}))

    print('young')
    npz = np.load(f'sp_pairs/Fst_{F}_h2pop_0.5_h2quad_0.0_no_unrel.npz')
    m, std, var = plot(npz['alpha'], npz['alpha_cov'],npz_ref['alpha'], npz_ref['alpha_cov'],non_sampling=True)
    df.append(pd.DataFrame({'Fst': [F], 
                            'method': ['Young et al. 2022'],
                            'non_sampling_var': m, 'non_sampling_var_std': std, 'var': var}))

    print('unified')                        
    npz= np.load(f'sp_pairs/Fst_{F}_h2pop_0.5_h2quad_0.0_with_unrel.npz')
    m, std, var = plot(npz['alpha'], npz['alpha_cov'],npz_ref['alpha'], npz_ref['alpha_cov'],non_sampling=True)
    df.append(pd.DataFrame({'Fst': [F], 
                            'method': ['unified'],
                            'non_sampling_var': m, 'non_sampling_var_std': std, 'var': var}))

    print('standard')
    npz = np.load(f'sp_pairs/Fst_{F}_h2pop_0.5_h2quad_0.0_standard_with_unrel.npz')
    m, std, var = plot(npz['alpha'], npz['alpha_cov'],npz_ref['alpha'], npz_ref['alpha_cov'],non_sampling=True)
    df.append(pd.DataFrame({'Fst': [F], 
                            'method': ['standard'],
                            'non_sampling_var': m, 'non_sampling_var_std': std, 'var': var}))

    print('sibdiff') 
    npz = np.load(f'sp_pairs/Fst_{F}_h2pop_0.5_h2quad_0.0_sibdiff.npz')
    m, std, var = plot(npz['alpha'], npz['alpha_cov'].squeeze(), npz_ref['alpha'], npz_ref['alpha_cov'],non_sampling=True)
    df.append(pd.DataFrame({'Fst': [F], 
                            'method': ['sib diff'],
                            'non_sampling_var': m, 'non_sampling_var_std': std, 'var': var}))
pd.concat(df).to_csv('paper_materials/sp_non_sampling_var_ref_to_standardGWAS_at_Fst0.001.csv', index=False, sep='\t')






df = []
print('Z var with unrel')
npz_ref = np.load(f'sp_pairs/Fst_0.001_h2pop_0.5_h2quad_0.0_standard_with_unrel.npz')
for F in [0, 0.001, 0.01, 0.1]:
    print('\n', F)
    print('robust')
    npz = np.load(f'sp_pairs/Fst_{F}_h2pop_0.5_h2quad_0.0_robust.npz')
    m, std, var = plot(npz['alpha'], npz['alpha_cov'], npz_ref['alpha'], npz_ref['alpha_cov'])
    df.append(pd.DataFrame({'Fst': [F], 
                            'method': ['robust'],
                            'non_sampling_var': m, 'non_sampling_var_std': std, 'var': var}))

    print('young')
    npz = np.load(f'sp_pairs/Fst_{F}_h2pop_0.5_h2quad_0.0_no_unrel.npz')
    m, std, var = plot(npz['alpha'], npz['alpha_cov'],npz_ref['alpha'], npz_ref['alpha_cov'])
    df.append(pd.DataFrame({'Fst': [F], 
                            'method': ['Young et al. 2022'],
                            'non_sampling_var': m, 'non_sampling_var_std': std, 'var': var}))

    print('unified')                        
    npz= np.load(f'sp_pairs/Fst_{F}_h2pop_0.5_h2quad_0.0_with_unrel.npz')
    m, std, var = plot(npz['alpha'], npz['alpha_cov'],npz_ref['alpha'], npz_ref['alpha_cov'])
    df.append(pd.DataFrame({'Fst': [F], 
                            'method': ['unified'],
                            'non_sampling_var': m, 'non_sampling_var_std': std, 'var': var}))

    print('standard')
    npz = np.load(f'sp_pairs/Fst_{F}_h2pop_0.5_h2quad_0.0_standard_with_unrel.npz')
    m, std, var = plot(npz['alpha'], npz['alpha_cov'],npz_ref['alpha'], npz_ref['alpha_cov'])
    df.append(pd.DataFrame({'Fst': [F], 
                            'method': ['standard'],
                            'non_sampling_var': m, 'non_sampling_var_std': std, 'var': var}))
    
    print('sibdiff')
    npz = np.load(f'sp_pairs/Fst_{F}_h2pop_0.5_h2quad_0.0_sibdiff.npz')
    m, std, var = plot(npz['alpha'], npz['alpha_cov'].squeeze(), npz_ref['alpha'], npz_ref['alpha_cov'])
    df.append(pd.DataFrame({'Fst': [F], 
                            'method': ['sib diff'],
                            'non_sampling_var': m, 'non_sampling_var_std': std, 'var': var}))
pd.concat(df).to_csv('paper_materials/sp_Z_var_ref_to_standardGWAS_at_Fst0.001.csv', index=False, sep='\t')


########### parent offspring

print('Parent offspinrg')
df = []
for F in [0, 0.001, 0.01, 0.1]:
    print('\n', F)
    npz = np.load(f'po_pairs/Fst_{F}_h2pop_0.5_h2quad_0.0_no_unrel.npz')
    m, std = plot(npz['alpha'], npz['alpha_cov'])
    df.append(pd.DataFrame({'Fst': [F], 
                            'method': ['Young et al. 2022/robust'],
                            'non_sampling_var': m, 'non_sampling_var_std': std}))

    npz= np.load(f'po_pairs/Fst_{F}_h2pop_0.5_h2quad_0.0_with_unrel.npz')
    m, std = plot(npz['alpha'], npz['alpha_cov'])
    df.append(pd.DataFrame({'Fst': [F], 
                            'method': ['unified'],
                            'non_sampling_var': m, 'non_sampling_var_std': std}))

    npz = np.load(f'po_pairs/Fst_{F}_h2pop_0.5_h2quad_0.0_standard_with_unrel.npz')
    m, std = plot(npz['alpha'], npz['alpha_cov'])
    df.append(pd.DataFrame({'Fst': [F], 
                            'method': ['standard'],
                            'non_sampling_var': m, 'non_sampling_var_std': std}))
    
    
    
    
pd.concat(df).to_csv('paper_materials/po_Z_var.csv', index=False, sep='\t')


df = []
for F in [0, 0.001, 0.01, 0.1]:
    print('\n', F)

    npz = np.load(f'po_pairs/Fst_{F}_h2pop_0.5_h2quad_0.0_no_unrel.npz')
    m, std = plot(npz['alpha'], npz['alpha_cov'],non_sampling=True)
    df.append(pd.DataFrame({'Fst': [F], 
                            'method': ['Young et al. 2022/robust'],
                            'non_sampling_var': m, 'non_sampling_var_std': std}))
                            
    npz= np.load(f'po_pairs/Fst_{F}_h2pop_0.5_h2quad_0.0_with_unrel.npz')
    m, std = plot(npz['alpha'], npz['alpha_cov'],non_sampling=True)
    df.append(pd.DataFrame({'Fst': [F], 
                            'method': ['unified'],
                            'non_sampling_var': m, 'non_sampling_var_std': std}))

    npz = np.load(f'po_pairs/Fst_{F}_h2pop_0.5_h2quad_0.0_standard_with_unrel.npz')
    m, std = plot(npz['alpha'], npz['alpha_cov'],non_sampling=True)
    df.append(pd.DataFrame({'Fst': [F], 
                            'method': ['standard'],
                            'non_sampling_var': m, 'non_sampling_var_std': std}))
    

pd.concat(df).to_csv('paper_materials/po_non_sampling_var.csv', index=False, sep='\t')