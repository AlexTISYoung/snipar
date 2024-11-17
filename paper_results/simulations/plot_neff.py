import matplotlib.pyplot as plt
import numpy as np

def compute_neff_ratio(sib_corr: float, f: float, n0: int, n1: int, phased: bool = False) -> float:
    """Compute theoretical effective sample size ratio of combined sample analysis versus related sample analysis.

    Args:
        sib_corr (float): sibling correlation.
        f (float): allele frequency.
        n0 (int): number of sib pairs in related sample.
        n1 (int): unrelated sample size.

    Returns:
        float: effective sample size ratio.
    """
    if phased:
        nu = 3 / 4
    else:
        nu = 3 / 4 - f * (1 - f) / 8.
    var_direct = np.linalg.inv(
        n0 / (1 - sib_corr ** 2) * np.array(
            [[2 - sib_corr, 2 * (1 - sib_corr)],
            [2 * (1 - sib_corr), 4 * nu * (1 - sib_corr)]]
        ) + n1
    )[0, 0]
    var_direct0 = 4 * nu * (1 - sib_corr ** 2) / (4 * n0 * (nu * (2 - sib_corr) + sib_corr - 1))
    return var_direct0 / var_direct

plt.rcParams.update({'font.size': 12})
plt.figure(figsize=(5,5))
n0 = 21341
n0 = 20_000
sib_corr = np.linspace(0, 0.99, 1000)
for n1 in np.arange(n0, 420_000, 70_000):
    neff = np.vectorize(compute_neff_ratio)(sib_corr, None, n0, n1, True)
    plt.plot(sib_corr, neff, label='$n_1$=' + f'{n1:,}')
plt.ylabel('Relative effective sample size (direct)')
plt.xlabel('Correlation between siblings\' residuals')
plt.legend()
plt.tight_layout()
plt.ylim([0.98,1.5])
plt.savefig('neff_vs_sibcorr.pdf')
plt.savefig('neff_vs_sibcorr.png')
plt.figure(figsize=(5,5))
f = np.linspace(0, 0.5, 1000)
for n1 in np.arange(n0, 420_000, 70_000):
    neff = np.vectorize(compute_neff_ratio)(0.52, f, n0, n1, False)
    plt.plot(f, neff, label='$n_1$=' + f'{n1:,}')
plt.ylabel('Relative effective sample size (direct)')
plt.xlabel('Minor Allele frequency')
plt.ylim([0.98,1.5])
plt.legend()
plt.tight_layout()
plt.savefig('neff_vs_maf.pdf')