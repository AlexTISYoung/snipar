import matplotlib.pyplot as plt
import pandas as pd
from functools import partial

def plot_confidence_interval_(x, mean, stdev, z=1.96, color='#2187bb', horizontal_line_width=0.02, label=''):


    left = x - horizontal_line_width / 2
    top = mean - 2 * stdev
    right = x + horizontal_line_width / 2
    bottom = mean + 2 * stdev
    plt.plot([x, x], [top, bottom], color=color, alpha=0.8)
    plt.plot([left, right], [top, top], color=color,alpha=0.8)
    plt.plot([left, right], [bottom, bottom], color=color,alpha=0.8)
    plt.plot(x, mean, color=color, label=label,alpha=0.8)

f = pd.read_csv('sp_non_sampling_var_ref_to_standardGWAS_at_Fst0.001.csv', sep='\t')


plot_confidence_interval = partial(plot_confidence_interval_, horizontal_line_width=1)
max_non_sampling_var = 1#f[(f['Fst'] == 0.1) & (f['method'] == 'standard')]['non_sampling_var'].values[0]
for F in [0, 0.001, 0.01, 0.1]:
    plt.figure(figsize=(4,3))
    # plt.figure()
    pop = (f[(f['Fst'] == F) & (f['method'] == 'sib diff')]['var']).values[0]
    # print(pop)
    # plt.ylim(bottom=-0.1, top=1.1)
    plt.ylim(bottom=-5, top=100)
    plt.xlim(right=14, left=0.1)
    # color = iter(cm.rainbow(np.linspace(0, 1, 5)))
    color = iter(['b', 'g', 'r', 'orange', 'm'])
    y = f[(f['Fst'] == F) & (f['method'] == 'sib diff')]['non_sampling_var']/max_non_sampling_var
    x = pop/f[(f['Fst'] == F) & (f['method'] == 'sib diff')]['var']
    std = f[(f['Fst'] == F) & (f['method'] == 'sib diff')]['non_sampling_var_std']
    c = next(color)
    plt.scatter(y=y,
        x=x, label='robust / sib-difference', c=[c], s=10)
    plot_confidence_interval(x, y, std, color=c)
    
    y = f[(f['Fst'] == F) & (f['method'] == 'robust')]['non_sampling_var']/max_non_sampling_var
    x = pop/f[(f['Fst'] == F) & (f['method'] == 'robust')]['var']
    std = f[(f['Fst'] == F) & (f['method'] == 'robust')]['non_sampling_var_std']
    c = next(color)
    plt.scatter(y=y,
        x=x, label='non-transmitted', c=[c], s=10)
    plot_confidence_interval(x, y, std, color=c)
    
    y = f[(f['Fst'] == F) & (f['method'] == 'Young et al. 2022')]['non_sampling_var']/max_non_sampling_var
    x = pop/f[(f['Fst'] == F) & (f['method'] == 'Young et al. 2022')]['var']
    std = f[(f['Fst'] == F) & (f['method'] == 'Young et al. 2022')]['non_sampling_var_std']
    c = next(color)
    plt.scatter(y=y,
        x=x, label='Young et al. 2022', c=[c], s=10)
    plot_confidence_interval(x, y, std, color=c)

    y = f[(f['Fst'] == F) & (f['method'] == 'unified')]['non_sampling_var']/max_non_sampling_var
    x = pop/f[(f['Fst'] == F) & (f['method'] == 'unified')]['var']
    std = f[(f['Fst'] == F) & (f['method'] == 'unified')]['non_sampling_var_std']
    c = next(color)
    plt.scatter(y=y,
        x=x, label='unified', c=[c], s=10)
    plot_confidence_interval(x, y, std, color=c)

    y = f[(f['Fst'] == F) & (f['method'] == 'standard')]['non_sampling_var']/max_non_sampling_var
    x = pop/f[(f['Fst'] == F) & (f['method'] == 'standard')]['var']
    std = f[(f['Fst'] == F) & (f['method'] == 'standard')]['non_sampling_var_std']
    c = next(color)
    plt.scatter(y=y,
        x=x, label='standard GWAS', c=[c], s=10)
    plot_confidence_interval(x, y, std, color=c)
    
    plt.axhline(y=0, linestyle='--', alpha=0.5,linewidth=1, c='grey')
    plt.axvline(x=1, linestyle='--', alpha=0.5,linewidth=1, c='grey')
    if F ==0:
        plt.legend(loc='upper right')
    # plt.legend()
    plt.ylabel('bias')
    plt.xlabel('relative effective sample size')
    # plt.show()
    plt.xticks([1, 5,9 ,13])
    # plt.yticks([0, 0.5, 1])
    plt.yticks([0, 45, 90])
    plt.text(1.3, 70, '$F_{st}$' + f'={F}')
    plt.tight_layout()
    plt.savefig(f'paper_materials/tradeoff/{F}.pdf')



plot_confidence_interval = partial(plot_confidence_interval_, horizontal_line_width=0.0333333333)
for F in [0, 0.001, 0.01, 0.1]:
    plt.figure(figsize=(4,3))

    pop = (f[(f['Fst'] == F) & (f['method'] == 'sib diff')]['var']).values[0]

    plt.ylim(bottom=-0.05, top=1.2)
    plt.xlim(right=1.4, left=0.97)

    color = iter(['b', 'g', 'r', 'orange', 'm'])
    y = f[(f['Fst'] == F) & (f['method'] == 'sib diff')]['non_sampling_var']/max_non_sampling_var
    x = pop/f[(f['Fst'] == F) & (f['method'] == 'sib diff')]['var']
    std = f[(f['Fst'] == F) & (f['method'] == 'sib diff')]['non_sampling_var_std']
    c = next(color)
    plt.scatter(y=y,
        x=x, label='robust / sib-difference', c=[c], s=10)
    plot_confidence_interval(x, y, std, color=c)
    
    y = f[(f['Fst'] == F) & (f['method'] == 'robust')]['non_sampling_var']/max_non_sampling_var
    x = pop/f[(f['Fst'] == F) & (f['method'] == 'robust')]['var']
    std = f[(f['Fst'] == F) & (f['method'] == 'robust')]['non_sampling_var_std']
    c = next(color)
    plt.scatter(y=y,
        x=x, label='non-transmitted', c=[c], s=10)
    plot_confidence_interval(x, y, std, color=c)
    
    y = f[(f['Fst'] == F) & (f['method'] == 'Young et al. 2022')]['non_sampling_var']/max_non_sampling_var
    x = pop/f[(f['Fst'] == F) & (f['method'] == 'Young et al. 2022')]['var']
    std = f[(f['Fst'] == F) & (f['method'] == 'Young et al. 2022')]['non_sampling_var_std']
    c = next(color)
    plt.scatter(y=y,
        x=x, label='Young et al. 2022', c=[c], s=10)
    plot_confidence_interval(x, y, std, color=c)

    y = f[(f['Fst'] == F) & (f['method'] == 'unified')]['non_sampling_var']/max_non_sampling_var
    x = pop/f[(f['Fst'] == F) & (f['method'] == 'unified')]['var']
    std = f[(f['Fst'] == F) & (f['method'] == 'unified')]['non_sampling_var_std']
    c = next(color)
    plt.scatter(y=y,
        x=x, label='unified', c=[c], s=10)
    plot_confidence_interval(x, y, std, color=c)

    
    plt.axhline(y=0, linestyle='--', alpha=0.5,linewidth=1, c='grey')
    plt.axvline(x=1, linestyle='--', alpha=0.5,linewidth=1, c='grey')
    if F == 0: plt.legend(loc='upper right')
    plt.ylabel('bias')
    plt.xlabel('relative effective sample size')
    plt.yticks([0, 0.5, 1])
    plt.text(1.01, 1, '$F_{st}$' + f'={F}')
    plt.tight_layout()
    plt.savefig(f'paper_materials/tradeoff/{F}_zoomed.pdf')





f = pd.read_csv('sp_non_sampling_var_ref_to_standardGWAS_at_Fst0.001_allsibs.csv', sep='\t')


plot_confidence_interval = partial(plot_confidence_interval_, horizontal_line_width=0.10852713178)
max_non_sampling_var = 1#f[(f['Fst'] == 0.1) & (f['method'] == 'standard')]['non_sampling_var'].values[0]
for F in [0, 0.001, 0.01, 0.1]:
    plt.figure(figsize=(4,3))
    pop = (f[(f['Fst'] == F) & (f['method'] == 'sib diff')]['var']).values[0]
    plt.ylim(bottom=-4, top=110)
    plt.xlim(right=2.3, left=0.9)
    color = iter(['b', 'g', 'r', 'm'])
    y = f[(f['Fst'] == F) & (f['method'] == 'sib diff')]['non_sampling_var']/max_non_sampling_var
    x = pop/f[(f['Fst'] == F) & (f['method'] == 'sib diff')]['var']
    std = f[(f['Fst'] == F) & (f['method'] == 'sib diff')]['non_sampling_var_std']
    c = next(color)
    plt.scatter(y=y,
        x=x, label='robust / sib-difference', c=[c], s=10)
    plot_confidence_interval(x, y, std, color=c)
    
    y = f[(f['Fst'] == F) & (f['method'] == 'robust')]['non_sampling_var']/max_non_sampling_var
    x = pop/f[(f['Fst'] == F) & (f['method'] == 'robust')]['var']
    std = f[(f['Fst'] == F) & (f['method'] == 'robust')]['non_sampling_var_std']
    c = next(color)
    plt.scatter(y=y,
        x=x, label='non-transmitted', c=[c], s=10)
    plot_confidence_interval(x, y, std, color=c)
    
    y = f[(f['Fst'] == F) & (f['method'] == 'Young et al. 2022/unified')]['non_sampling_var']/max_non_sampling_var
    x = pop/f[(f['Fst'] == F) & (f['method'] == 'Young et al. 2022/unified')]['var']
    std = f[(f['Fst'] == F) & (f['method'] == 'Young et al. 2022/unified')]['non_sampling_var_std']
    c = next(color)
    plt.scatter(y=y,
        x=x, label='Young et al. 2022', c=[c], s=10)
    plot_confidence_interval(x, y, std, color=c)


    y = f[(f['Fst'] == F) & (f['method'] == 'standard')]['non_sampling_var']/max_non_sampling_var
    x = pop/f[(f['Fst'] == F) & (f['method'] == 'standard')]['var']
    std = f[(f['Fst'] == F) & (f['method'] == 'standard')]['non_sampling_var_std']
    c = next(color)
    plt.scatter(y=y,
        x=x, label='standard GWAS', c=[c], s=10)
    plot_confidence_interval(x, y, std, color=c)
    
    plt.axhline(y=0, linestyle='--', alpha=0.5,linewidth=1, c='grey')
    plt.axvline(x=1, linestyle='--', alpha=0.5,linewidth=1, c='grey')
    if F ==0:
        plt.legend(loc='upper right')
    # plt.legend()
    plt.ylabel('bias')
    plt.xlabel('relative effective sample size')
    # plt.show()
    plt.xticks([1, 1.4, 1.8, 2.2])
    # plt.yticks([0, 0.5, 1])
    plt.annotate('$F_{st}$' + f'={F}', xy=(0.1, 0.9), xycoords='axes fraction',)
                # horizontalalignment='right', verticalalignment='bottom')
    # plt.text(1.01, 100, '$F_{st}$' + f'={F}')
    plt.tight_layout()
    plt.savefig(f'paper_materials/tradeoff/{F}_allsibs.pdf')


plot_confidence_interval = partial(plot_confidence_interval_, horizontal_line_width=0.01240310077)
for F in [0, 0.001, 0.01, 0.1]:
    plt.figure(figsize=(4,3))
    # plt.figure()
    pop = (f[(f['Fst'] == F) & (f['method'] == 'sib diff')]['var']).values[0]

    plt.ylim(bottom=-0.01, top=0.25)
    plt.xlim(right=1.15, left=0.99)
    color = iter(['b', 'g', 'r', 'm'])
    y = f[(f['Fst'] == F) & (f['method'] == 'sib diff')]['non_sampling_var']/max_non_sampling_var
    x = pop/f[(f['Fst'] == F) & (f['method'] == 'sib diff')]['var']
    std = f[(f['Fst'] == F) & (f['method'] == 'sib diff')]['non_sampling_var_std']
    c = next(color)
    plt.scatter(y=y,
        x=x, label='robust / sib-difference', c=[c], s=10)
    plot_confidence_interval(x, y, std, color=c)
    
    y = f[(f['Fst'] == F) & (f['method'] == 'robust')]['non_sampling_var']/max_non_sampling_var
    x = pop/f[(f['Fst'] == F) & (f['method'] == 'robust')]['var']
    std = f[(f['Fst'] == F) & (f['method'] == 'robust')]['non_sampling_var_std']
    c = next(color)
    plt.scatter(y=y,
        x=x, label='non-transmitted', c=[c], s=10)
    plot_confidence_interval(x, y, std, color=c)
    
    y = f[(f['Fst'] == F) & (f['method'] == 'Young et al. 2022/unified')]['non_sampling_var']/max_non_sampling_var
    x = pop/f[(f['Fst'] == F) & (f['method'] == 'Young et al. 2022/unified')]['var']
    std = f[(f['Fst'] == F) & (f['method'] == 'Young et al. 2022/unified')]['non_sampling_var_std']
    c = next(color)
    plt.scatter(y=y,
        x=x, label='Young et al. 2022', c=[c], s=10)
    plot_confidence_interval(x, y, std, color=c)
    
    plt.axhline(y=0, linestyle='--', alpha=0.5,linewidth=1, c='grey')
    plt.axvline(x=1, linestyle='--', alpha=0.5,linewidth=1, c='grey')
    if F == 0: plt.legend(loc='upper right')
    plt.ylabel('bias')
    plt.xlabel('relative effective sample size')
    plt.xticks([1, 1.05, 1.1, 1.15])
    plt.yticks([0, 0.1, 0.2])
    plt.text(1.01, 0.22, '$F_{st}$' + f'={F}')
    plt.tight_layout()
    plt.savefig(f'paper_materials/tradeoff/{F}_allsibs_zoomed.pdf')