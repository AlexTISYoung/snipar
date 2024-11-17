import sys
import pandas as pd
import statsmodels.formula.api as sm
import seaborn as sns


f = '/disk/genetics/data/mcs/private/latest/raw/downloaded/NCDS_SFTP_1TB_1/imputed/'
eur = pd.read_csv(f+'filter_extract/eur_samples_mz_removed.txt', sep='\s+', header=None)
sas = pd.read_csv(f+'filter_extract/sas_samples_mz_removed.txt', sep='\s+', header=None)
def string(sample, phen_str): 
    return \
f"""
phen= pd.read_csv(f+'phen/{phen_str}/{sample}/pheno.pheno', sep='\s+', header=None)
# phen={sample}.merge(phen.dropna(), on=[0,1])
king = pd.read_csv(f+'king/mcs_{sample}.kin0', sep='\s+')
print('{phen_str}')
# print(king.loc[(king['ID1'].isin(phen[1])) & (king['ID2'].isin(phen[1])), 'InfType'].value_counts())
# print( king.loc[(king['ID1'].isin(phen[1])) & (king['ID2'].isin(phen[1])), ['ID1', 'ID2', 'IBD1Seg', 'InfType']])
tmp = king.loc[(king['ID1'].isin(phen[1])) & (king['ID2'].isin(phen[1])), ['ID1', 'ID2', 'IBD1Seg', 'InfType']]
# tmp = king.loc[(king['ID1'].isin(phen[1])) & (king['ID2'].isin(phen[1])) & (king['InfType'] == 'PO'), ['ID1', 'ID2', 'IBD1Seg', 'InfType']]
print(tmp)
lis_{sample}_{phen_str} = []
for row in tmp.itertuples():
    if row.InfType == "PO":
        lis_{sample}_{phen_str}.append(row.ID1)
        lis_{sample}_{phen_str}.append(row.ID2)
    else: lis_{sample}_{phen_str}.append(row.ID1)
print()
"""
exec(string('eur', 'ea'))
exec(string('eur', 'bmi'))
exec(string('eur', 'height'))
exec(string('sas', 'ea'))
exec(string('sas', 'bmi'))
exec(string('sas', 'height'))



# Family PGS
def main(method: str, phen_: str, ancestry: str, population=False):
    f = f'/disk/genetics/data/mcs/private/latest/raw/downloaded/NCDS_SFTP_1TB_1/imputed/'
    f_scores = f'/disk/genetics/ukb/jguan/ukb_analysis/output/prs_analysis/mcs/nofilter/{phen_}/{method}/{ancestry}/' #direct.pgs.txt'

    import pandas as pd

    if population and method == 'unified':
        pgs = pd.read_csv(f_scores + 'population.pgs.txt', sep='\s+')
    # else: pgs = pd.read_csv(f_scores + 'direct/scoresout.sscore', sep='\s+')
    else: pgs = pd.read_csv(f_scores + 'direct.pgs.txt', sep='\s+')
    phen = pd.read_csv(f'/var/genetics/data/mcs/private/latest/raw/downloaded/NCDS_SFTP_1TB_1/imputed/phen/phenotypes_{ancestry}.txt', sep='\s+')[['FID', 'IID', phen_]].dropna()
    # phen = pd.read_csv(f'/var/genetics/data/mcs/private/latest/raw/genotyped/NCDS_SFTP_1TB_1/imputed/phen/bmi_height_sweep7.txt', sep='\s+')[['FID', 'IID', phen_+'7']].dropna()
    query = f'IID not in @lis_{ancestry}_{phen_}'
    # phen = pd.read_csv(f + f'/phen/{phen}/pheno.pheno', sep='\s+', header=None)
    phen.columns = ['FID', 'IID', 'pheno']
    phen = phen.query(query)
    cov_ = '_' + ancestry if ancestry == 'sas' else ''
    # covar = pd.read_csv(f + f'/phen/covar{cov_}.txt', sep='\s+')
    covar = pd.read_csv(f + f'/phen/PCs{cov_}.txt', sep='\s+')
    if 'FATHER_ID' in pgs.columns:
        pgs= pgs.drop('FATHER_ID', axis=1)
    if 'MOTHER_ID' in pgs.columns:
        pgs= pgs.drop('MOTHER_ID', axis=1)
    assert pgs.shape == pgs.dropna().shape
    assert phen.shape == phen.dropna().shape
    assert covar.shape == covar.dropna().shape

    X = phen.merge(covar, on=['FID', 'IID'], how='inner')
    X = X.drop('FID', axis=1)
    pgs = pgs.drop('FID', axis=1)
    X = X.merge(pgs, on=['IID'], how='inner')
    # if ancestry == 'eur':
    #     X = X.sample(n=658, replace=False, random_state=1)
    X['proband'] /= X['proband'].std()
    X['maternal'] /= X['maternal'].std()
    X['paternal'] /= X['paternal'].std()
    X['pheno'] /= X['pheno'].std()
    string = '+'.join([col for col in X.columns if col != 'FID' and col != 'IID' and col != 'pheno'])
    result = sm.ols(formula="pheno ~" + string, data=X).fit()
    print('population' if population and method == 'unified' else method, 
          phen_.upper() if phen_ != 'height' else 'Height', ancestry.upper(), result.params['proband'], result.bse['proband'], sep='\t') #[result.params['proband'] - 1.96 * result.bse['proband'], result.params['proband'] + 1.96 * result.bse['proband']])


    import sys
with open('mcs_pgi_fam_reg_results_PCadjusted_new_robust.txt', 'w'): pass # empty the file
original_stdout = sys.stdout
with open('mcs_pgi_fam_reg_results_PCadjusted_new_robust.txt', 'w') as f:
    sys.stdout = f # redirect the output to file
    # filename = 'mcs_pgi_simple_reg_results.txt'
    print('method', 'Phen', 'Ancestry', 'Beta', 'SE')
    main('unified', 'bmi', 'eur')
    main('young', 'bmi', 'eur')
    main('new_robust', 'bmi', 'eur')
    main('sibdiff', 'bmi', 'eur')
    main('unified', 'bmi', 'eur', True)

    main('unified', 'height', 'eur')
    main('young', 'height', 'eur')
    main('new_robust', 'height', 'eur')
    main('sibdiff', 'height', 'eur')
    main('unified', 'height', 'eur', True)


    main('unified', 'ea', 'eur')
    main('young', 'ea', 'eur')
    main('new_robust', 'ea', 'eur')
    main('sibdiff', 'ea', 'eur')
    main('unified', 'ea', 'eur', True)


    sys.stdout = original_stdout




    x = pd.read_csv('tammy_eur_mcs_results.txt', sep='\s+')


# Population PGS

def main(method: str, phen_: str, ancestry: str, population=False):
    f = f'/disk/genetics/data/mcs/private/latest/raw/downloaded/NCDS_SFTP_1TB_1/imputed/'
    f_scores = f'/disk/genetics/ukb/jguan/ukb_analysis/output/prs_analysis/mcs/nofilter/{phen_}/{method}/{ancestry}/pop_pgs/' #direct.pgs.txt'

    import pandas as pd

    if population and method == 'unified':
        pgs = pd.read_csv(f_scores + 'population.pgs.txt', sep='\s+')
    else: pgs = pd.read_csv(f_scores + 'direct.pgs.txt', sep='\s+')
    phen = pd.read_csv(f'/var/genetics/data/mcs/private/latest/raw/downloaded/NCDS_SFTP_1TB_1/imputed/phen/phenotypes_{ancestry}.txt', sep='\s+')[['FID', 'IID', phen_]].dropna()

    query = f'IID not in @lis_{ancestry}_{phen_}'
    # phen = pd.read_csv(f + f'/phen/{phen}/pheno.pheno', sep='\s+', header=None)
    phen.columns = ['FID', 'IID', 'pheno']
    phen = phen.query(query)
    cov_ = '_' + ancestry if ancestry == 'sas' else ''
    # covar = pd.read_csv(f + f'/phen/covar{cov_}.txt', sep='\s+')
    covar = pd.read_csv(f + f'/phen/PCs{cov_}.txt', sep='\s+')
    assert pgs.shape == pgs.dropna().shape
    assert phen.shape == phen.dropna().shape
    assert covar.shape == covar.dropna().shape

    X = phen.merge(covar, on=['FID', 'IID'], how='inner')
    X = X.merge(pgs, on=['FID', 'IID'], how='inner')

    X['proband'] /= X['proband'].std()
    X['pheno'] /= X['pheno'].std()

    string = '+'.join([col for col in X.columns if col != 'FID' and col != 'IID' and col != 'pheno'])
    result = sm.ols(formula="pheno ~" + string, data=X).fit()

    print('population' if population and method == 'unified' else method, 
          phen_.upper() if phen_ != 'height' else 'Height', ancestry.upper(), result.params['proband'], result.bse['proband'], sep='\t') #[result.params['proband'] - 1.96 * result.bse['proband'], result.params['proband'] + 1.96 * result.bse['proband']])



import sys
with open('mcs_pgi_simple_reg_results_PCadjusted_new_robust.txt', 'w'): pass # empty the file
original_stdout = sys.stdout
with open('mcs_pgi_simple_reg_results_PCadjusted_new_robust.txt', 'w') as f:
    sys.stdout = f # redirect the output to file
    print('method', 'Phen', 'Ancestry', 'Beta', 'SE')
    main('unified', 'bmi', 'eur')
    main('young', 'bmi', 'eur')
    main('new_robust', 'bmi', 'eur')
    main('sibdiff', 'bmi', 'eur')
    main('unified', 'bmi', 'eur', True)

    main('unified', 'bmi', 'sas')
    main('young', 'bmi', 'sas')
    main('new_robust', 'bmi', 'sas')
    main('sibdiff', 'bmi', 'sas')
    main('unified', 'bmi', 'sas', True)

    main('unified', 'height', 'eur')
    main('young', 'height', 'eur')
    main('new_robust', 'height', 'eur')
    main('sibdiff', 'height', 'eur')
    main('unified', 'height', 'eur', True)

    main('unified', 'height', 'sas')
    main('young', 'height', 'sas')
    main('new_robust', 'height', 'sas')
    main('sibdiff', 'height', 'sas')
    main('unified', 'height', 'sas', True)

    main('unified', 'ea', 'eur')
    main('young', 'ea', 'eur')
    main('new_robust', 'ea', 'eur')
    main('sibdiff', 'ea', 'eur')
    main('unified', 'ea', 'eur', True)

    main('unified', 'ea', 'sas')
    main('young', 'ea', 'sas')
    main('new_robust', 'ea', 'sas')
    main('sibdiff', 'ea', 'sas')
    main('unified', 'ea', 'sas', True)

    sys.stdout = original_stdout


import matplotlib.pyplot as plt
f, axes = plt.subplots(1, 3, figsize=(15,5))




z = pd.read_csv('mcs_pgi_fam_reg_results_PCadjusted_new_robust.txt', sep='\s+')
import seaborn as sns
z.replace({'method': {'sibdiff': 'sib-difference', 'young': 'Young et al. 2022', 'population': 'standard GWAS', 'new_robust': 'robust'}}, inplace=True)
z.replace({'Phen': {'EA': 'average GCSE grade\n(EA PGIs)', 'bmi': 'BMI', 'Height': 'height'}}, inplace=True)
ax = sns.barplot(x=z['Phen'], y=z['Beta'], hue=z['method'], hue_order=['sib-difference', 'robust', 'Young et al. 2022', 'unified', 'standard GWAS'], ax=axes[0], order=['BMI', 'height', 'average GCSE grade\n(EA PGIs)'])
sns.move_legend(ax, "upper left",)
x_coords = [p.get_x() + 0.5 * p.get_width() for p in ax.patches]
y_coords = [p.get_height() for p in ax.patches]
ax.errorbar(x=x_coords, y=y_coords, yerr=1.96*z['SE'], fmt="none", c="k")
ax.set(xlabel=None, ylabel='partial correlation')
ax.set_ylim(-0.1, 0.55)
ax.set_yticks([0., 0.1, 0.2, 0.3, 0.4, 0.5])
ax.set_title('EUR (direct effects)')




x = pd.read_csv('mcs_pgi_simple_reg_results_PCadjusted_new_robust.txt', sep='\s+')
# z = x[(x['Population'] ==False) ] #& (x['Phenotype'] != 'EA')]
z = x[(x['Ancestry'] == 'EUR')]
z.replace({'Phen': {'EA': 'average GCSE grade\n(EA PGIs)', 'bmi': 'BMI', 'Height': 'height'}}, inplace=True)
z.replace({'method': {'sibdiff': 'sib-difference', 'young': 'Young et al. 2022', 'population': 'standard GWAS', 'new_robust': 'robust'}}, inplace=True)
ax = sns.barplot(x=z['Phen'], y=z['Beta'], hue=z['method'], hue_order=['sib-difference', 'robust', 'Young et al. 2022', 'unified', 'standard GWAS'], ax=axes[1])
ax.legend_.remove()
# sns.move_legend(ax, "upper left",)
x_coords = [p.get_x() + 0.5 * p.get_width() for p in ax.patches]
y_coords = [p.get_height() for p in ax.patches]
ax.errorbar(x=x_coords, y=y_coords, yerr=1.96*z['SE'], fmt="none", c="k")
ax.set(xlabel=None, ylabel='correlation')
ax.set_ylim(-0.1, 0.55)
ax.set_yticks([0., 0.1, 0.2, 0.3, 0.4, 0.5])
ax.set_title('EUR (population effects)')



z = x[(x['Ancestry'] == 'SAS')]
z.replace({'Phen': {'EA': 'average GCSE grade\n(EA PGIs)',  'bmi': 'BMI', 'Height': 'height'}}, inplace=True)
z.replace({'method': {'sibdiff': 'sib-difference', 'young': 'Young et al. 2022', 'population': 'standard GWAS', 'new_robust': 'robust'}}, inplace=True)
ax = sns.barplot(x=z['Phen'], y=z['Beta'], hue=z['method'], hue_order=['sib-difference', 'robust', 'Young et al. 2022', 'unified', 'standard GWAS'], ax=axes[2])
ax.legend_.remove()
x_coords = [p.get_x() + 0.5 * p.get_width() for p in ax.patches]
y_coords = [p.get_height() for p in ax.patches]
ax.errorbar(x=x_coords, y=y_coords, yerr=1.96*z['SE'], fmt="none", c="k")
ax.set(xlabel=None, ylabel='correlation')
ax.set_ylim(-0.1, 0.55)
ax.set_yticks([0., 0.1, 0.2, 0.3, 0.4, 0.5])
ax.set_title('SAS (population effects)')
plt.tight_layout()
plt.savefig('mcs_pgs_with_population_PCadjusted_nofilter_new_robust.pdf')