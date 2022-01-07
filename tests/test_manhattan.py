import h5py, code
import numpy as np
import pandas as pd
from pysnptools.snpreader import Bed
from scipy.stats import norm
from pandas import DataFrame
import matplotlib.pyplot as plt
from tqdm import tqdm
import statsmodels.api as sm
from scipy.stats import linregress
def imputation_test(chromosomes,
                   imputed_prefix = '',
                   expected_prefix = "../UKBioRDE_revision/data/tmp/filtered_ukb_chr",
                   start = None,
                   end = None,
                   backup_threshold=0.01,
                   ):
    #Data files for chromosome i should be named in this fashion: "prefix{i}"
    chromosomes_expected_genes_o = []
    chromosomes_expected_genes_pm = []
    chromosomes_imputed_genes_o = []
    chromosomes_imputed_genes_pm = []
    chrom_sizes = []
    for chromosome in tqdm(chromosomes):
        with h5py.File(imputed_prefix+str(chromosome)+".hdf5",'r') as f:
            gts = np.array(f["imputed_par_gts"]).astype(float)
            fids = np.array(f["families"]).astype(str)
            ped_array = np.array(f["pedigree"]).astype(str)
            ped = pd.DataFrame(ped_array[1:], columns = ped_array[0])
            non_duplicates = np.array(f["non_duplicates"])
            parent_ratio_backup = np.array(f["parent_ratio_backup"])
            sib_ratio_backup = np.array(f["sib_ratio_backup"])
            #Fix backup ratio start end
            if start is not None:
                non_duplicates = non_duplicates+start
        expected = Bed(expected_prefix+str(chromosome)+".bed", count_A1 = True)
        expected_gts = expected[:, non_duplicates].read().val.astype(float)
        chrom_sizes.append(gts.shape[1])
        expected_ids = expected.iid
        iid_to_bed_index = {i:index for index, i in enumerate(expected_ids[:,1])}
        #fids of control families start with _
        #this has the predix _*_
        index_of_families_in_imputation = {fid:index for index,fid in enumerate(fids)}
        # no parent control starts with _o_
        # only has father control starts with _p_
        # only has father control starts with _m_
        control_o_families = list({row["FID"][3:] for index, row in ped.iterrows() if row["FID"].startswith("_o_")})
        #for each family select id of the parents
        parent_ids = ped.groupby("FID").agg({
                                    'FATHER_ID':lambda x: ([a for a in list(x) if a in ped["IID"].tolist()]+[None])[0],
                                    'MOTHER_ID':lambda x: ([a for a in list(x) if a in ped["IID"].tolist()]+[None])[0],
                                    })
        parents_of_control_o_families = parent_ids.loc[control_o_families]
        mother_indexes_control_o = [iid_to_bed_index[parents_of_control_o_families.loc[i, "MOTHER_ID"]] for i in control_o_families]
        father_indexes_control_o = [iid_to_bed_index[parents_of_control_o_families.loc[i, "FATHER_ID"]] for i in control_o_families]
        expected_parent_gts_control_o = (expected_gts[mother_indexes_control_o,:]+expected_gts[father_indexes_control_o,:])/2
        expected_genes_o = expected_parent_gts_control_o
        index_of_control_families_in_imputation_o = [index_of_families_in_imputation["_o_"+i] for i in control_o_families]
        imputed_genes_o = gts[index_of_control_families_in_imputation_o,:]
        control_p = list({row["FID"][3:] for index, row in ped.iterrows() if row["FID"].startswith("_p_")})
        control_m = list({row["FID"][3:] for index, row in ped.iterrows() if row["FID"].startswith("_m_")})
        control_pm_families = control_p + control_m
        parent_of_control_m = parent_ids.loc[control_m]
        parent_of_control_p = parent_ids.loc[control_p]
        father_indexes_control_m = [iid_to_bed_index[parent_of_control_m.loc[i, "FATHER_ID"]] for i in control_m]
        mother_indexes_control_p = [iid_to_bed_index[parent_of_control_p.loc[i, "MOTHER_ID"]] for i in control_p]
        expected_parent_gts_control_pm = expected_gts[mother_indexes_control_p + father_indexes_control_m, :]
        expected_genes_pm = expected_parent_gts_control_pm
        index_of_control_families_in_imputation_pm = [index_of_families_in_imputation["_p_" + i] for i in control_p] + [index_of_families_in_imputation["_m_" + i] for i in control_m]
        imputed_genes_pm = gts[index_of_control_families_in_imputation_pm,:]
        chromosomes_expected_genes_o.append(expected_genes_o)
        chromosomes_expected_genes_pm.append(expected_genes_pm)
        chromosomes_imputed_genes_o.append(imputed_genes_o)
        chromosomes_imputed_genes_pm.append(imputed_genes_pm)
    
    whole_expected_genes_o = np.hstack(chromosomes_expected_genes_o)
    whole_imputed_genes_o = np.hstack(chromosomes_imputed_genes_o)
    whole_expected_genes_pm = np.hstack(chromosomes_expected_genes_pm)
    whole_imputed_genes_pm = np.hstack(chromosomes_imputed_genes_pm)
    pm_coefs = []
    pm_pvals = []
    pop_size, nsnps = whole_imputed_genes_pm.shape
    for i in tqdm(range(nsnps)):
        x = whole_imputed_genes_pm[:,i]
        y = whole_expected_genes_pm[:,i]
        mask = ~(np.isnan(x) | np.isnan(y))
        result = sm.OLS(y[mask], x[mask]).fit()
        pm_coefs.append(result.params[0])
        pm_pvals.append(result.t_test(([1],1)).pvalue)
    pm_pvals_log10 = -np.log10(pm_pvals)
    numerical = ~(np.isnan(pm_pvals_log10) | np.isinf(pm_pvals_log10) | (parent_ratio_backup > backup_threshold))
    counter = 0
    pm_chrom_sizes = [i for i in chrom_sizes]
    for i in range(len(chrom_sizes)):
        tmp = chrom_sizes[i]
        pm_chrom_sizes[i] = np.sum(numerical[counter:counter+chrom_sizes[i]])
        counter += tmp
    pm_pvals_log10 = pm_pvals_log10[numerical]
    pm_coefs = np.array(pm_coefs)[numerical]
    pm_pvals = np.array(pm_pvals)[numerical]
    print("pm_coefs.mean()", pm_coefs.mean())
    tmp = np.random.uniform(size=pm_pvals_log10.shape)
    qs = -np.log10(tmp)
    print("coefs", min(pm_coefs), max(pm_coefs))
    print([np.min(tmp), np.max(tmp)])
    print([np.min(pm_pvals), np.max(pm_pvals)])
    print([np.min(qs), np.max(qs)])
    print([np.min(pm_pvals_log10), np.max(pm_pvals_log10)])
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    ax.plot(np.sort(qs), np.sort(pm_pvals_log10), '-b')
    ax.plot(np.sort(qs), np.sort(pm_pvals_log10), '.r')
    ax.plot([np.min(qs), np.max(qs)], [np.min(qs), np.max(qs)], 'y')
    plt.xlabel("expected(-log10(p-value))")
    plt.ylabel("observed(-log10(p-value))")
    plt.title(f"{imputed_prefix}\nqqplot")
    print(f"{imputed_prefix}_pm_qq.png")

    plt.savefig(f"{imputed_prefix}_pm_qq")
    plt.clf()



    # sample data
    nsnps = len(pm_coefs)
    all_chrom_names = []
    for i in range(len(chromosomes)):
        all_chrom_names += [f'ch-{chromosomes[i]}']*pm_chrom_sizes[i]
    df = DataFrame({'gene' : [f'gene-{i}' % i for i in np.arange(nsnps)],
    'pm_coef' : pm_coefs,
    'pm_pvals' : pm_pvals,
    'pm_pvals_log10': pm_pvals_log10,
    'chromosome' : all_chrom_names})

    for col in ['pm_coef', 'pm_pvals', 'pm_pvals_log10']:
        df.chromosome = df.chromosome.astype('category')
        df.chromosome = df.chromosome.cat.set_categories([f'ch-{i}' for i in chromosomes], ordered=True)
        df = df.sort_values('chromosome')

        # How to plot gene vs. -log10(pvalue) and colour it by chromosome?
        df['ind'] = range(len(df))
        df_grouped = df.groupby(('chromosome'))

        # manhattan plot
        fig = plt.figure(figsize=(40, 8)) # Set the figure size
        ax = fig.add_subplot(111)
        colors = ['darkred','darkgreen','darkblue', 'darkslategray']
        x_labels = []
        x_labels_pos = []
        for num, (name, group) in enumerate(df_grouped):
            # group.plot(kind='scatter', x='ind', y='coef',color=colors[num % len(colors)], ax=ax, size=1)
            ax.scatter(group["ind"], group[col], color=colors[num % len(colors)], s=1)
            x_labels.append(name) 
            x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))
        ax.set_xticks(x_labels_pos)
        ax.set_xticklabels(x_labels)

        # set axis limits
        ax.set_xlim([0, len(df)])
        ax.set_ylim([0, np.ma.masked_invalid(df[col]).max()+0.2])
        ax.grid(True)

        # x axis label
        ax.set_xlabel('Chromosome')
        ax.set_title(f"{imputed_prefix}\n{col}")
        plt.savefig(f"{imputed_prefix}_pm_manhattan_{col}")
        plt.clf()
        print(f"{imputed_prefix}_pm_manhattan_{col}.png")


        #------------------------------------------------------------------
    o_coefs = []
    o_pvals = []
    pop_size, nsnps = whole_imputed_genes_o.shape
    for i in tqdm(range(nsnps)):
        x = whole_imputed_genes_o[:,i]
        y = whole_expected_genes_o[:,i]
        mask = ~(np.isnan(x) | np.isnan(y))
        result = sm.OLS(y[mask], x[mask]).fit()
        o_coefs.append(result.params[0])
        o_pvals.append(result.t_test(([1],1)).pvalue)
        # result = linregress(x[mask],y[mask])
        # o_coefs.append(result.slope)
        # t=(result.slope-1)/result.stderr
        # o_pvals.append(2*(1-norm.cdf(abs(t))))
    o_pvals_log10 = -np.log10(o_pvals)
    numerical = ~(np.isnan(o_pvals_log10) | np.isinf(o_pvals_log10) | (sib_ratio_backup > backup_threshold))
    counter = 0
    o_chrom_sizes = [i for i in chrom_sizes]
    for i in range(len(chrom_sizes)):
        tmp = chrom_sizes[i]
        o_chrom_sizes[i] = np.sum(numerical[counter:counter+chrom_sizes[i]])
        counter += tmp
    o_pvals_log10 = o_pvals_log10[numerical]
    o_coefs = np.array(o_coefs)[numerical]
    o_pvals = np.array(o_pvals)[numerical]
    print("o_coefs.mean()", o_coefs.mean())

    tmp = np.random.uniform(size=o_pvals_log10.shape)
    qs = -np.log10(tmp)
    print([np.min(tmp), np.max(tmp)])
    print([np.min(o_pvals), np.max(o_pvals)])
    print([np.min(qs), np.max(qs)])
    print([np.min(o_pvals_log10), np.max(o_pvals_log10)])
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)    
    ax.plot(np.sort(qs), np.sort(o_pvals_log10), '-b')
    ax.plot(np.sort(qs), np.sort(o_pvals_log10), '.r')
    ax.plot([np.min(qs), np.max(qs)], [np.min(qs), np.max(qs)], 'y')
    plt.xlabel("expected(-log10(p-value))")
    plt.ylabel("observed(-log10(p-value))")
    plt.title(f"{imputed_prefix}\nqqplot")
    print(f"{imputed_prefix}_o_qq.png")
    plt.savefig(f"{imputed_prefix}_o_qq")
    plt.clf()



    # sample data
    nsnps = len(o_coefs)
    all_chrom_names = []
    for i in range(len(chromosomes)):
        all_chrom_names += [f'ch-{chromosomes[i]}']*o_chrom_sizes[i]
    print(nsnps)
    print(o_coefs.shape, o_pvals.shape, o_pvals_log10.shape)
    print(len(all_chrom_names))
    df = DataFrame({'gene' : [f'gene-{i}' % i for i in np.arange(nsnps)],
    'o_coef' : o_coefs,
    'o_pvals' : o_pvals,
    'o_pvals_log10': o_pvals_log10,
    'chromosome' : all_chrom_names})

    for col in ['o_coef', 'o_pvals', 'o_pvals_log10']:
        df.chromosome = df.chromosome.astype('category')
        df.chromosome = df.chromosome.cat.set_categories([f'ch-{i}' for i in chromosomes], ordered=True)
        df = df.sort_values('chromosome')

        # How to plot gene vs. -log10(pvalue) and colour it by chromosome?
        df['ind'] = range(len(df))
        df_grouped = df.groupby(('chromosome'))

        # manhattan plot
        fig = plt.figure(figsize=(40, 8)) # Set the figure size
        ax = fig.add_subplot(111)
        colors = ['darkred','darkgreen','darkblue', 'darkslategray']
        x_labels = []
        x_labels_pos = []
        for num, (name, group) in enumerate(df_grouped):
            # group.plot(kind='scatter', x='ind', y='coef',color=colors[num % len(colors)], ax=ax, size=1)
            ax.scatter(group["ind"], group[col], color=colors[num % len(colors)], s=1)
            x_labels.append(name) 
            x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))
        ax.set_xticks(x_labels_pos)
        ax.set_xticklabels(x_labels)

        # set axis limits
        ax.set_xlim([0, len(df)])
        ax.set_ylim([0, np.ma.masked_invalid(df[col]).max()+0.2])
        ax.grid(True)

        # x axis label
        ax.set_xlabel('Chromosome')
        ax.set_title(f"{imputed_prefix}\n{col}")
        plt.savefig(f"{imputed_prefix}_o_manhattan_{col}")
        plt.clf()
        print(f"{imputed_prefix}_o_manhattan_{col}.png")
