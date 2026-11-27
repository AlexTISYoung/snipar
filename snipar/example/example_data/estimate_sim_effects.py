"""Compare tutorial GWAS estimates to the simulated SNP effects in variant order."""

import argparse

import h5py
import numpy as np


def slope_and_se(truth, estimates):
    valid = np.isfinite(truth) & np.isfinite(estimates)
    truth, estimates = truth[valid], estimates[valid]
    if len(truth) < 3:
        raise ValueError('At least three variants are required to estimate bias and its standard error')
    x = truth - truth.mean()
    y = estimates - estimates.mean()
    sum_squares = x @ x
    if sum_squares == 0:
        raise ValueError('The simulated effects must vary across variants')
    slope = (x @ y) / sum_squares
    residual = y - slope * x
    return slope, np.sqrt((residual @ residual) / ((len(x) - 2) * sum_squares))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('estimated_effects', help='HDF5 file with estimated effects')
    parser.add_argument('true_effects', help='Text file with simulated effects in the same variant order')
    args = parser.parse_args(argv)
    with h5py.File(args.estimated_effects, 'r') as source:
        estimates = source['estimate'][:].astype(float)
    truth = np.loadtxt(args.true_effects, ndmin=2)
    if estimates.shape[0] != truth.shape[0]:
        parser.error('Estimated and true effects must cover the same variants in the same order')
    comparisons = [('direct genetic effects', truth[:, 0], estimates[:, 0])]
    if estimates.shape[1] == 3:
        if truth.shape[1] < 3:
            parser.error('The true effects file must contain direct, paternal and maternal effects')
        avg = (estimates[:, 1] + estimates[:, 2]) * .5
        comparisons.extend([
            ('paternal NTCs', truth[:, 1], estimates[:, 1]),
            ('maternal NTCs', truth[:, 2], estimates[:, 2]),
            ('average NTCs', (truth[:, 1] + truth[:, 2]) * .5, avg),
            ('population effects as estimates of direct effects', truth[:, 0], estimates[:, 0] + avg),
        ])
    for label, expected, fitted in comparisons:
        slope, se = slope_and_se(expected, fitted)
        print(f'Bias for {label}: {round(slope - 1, 4)}; S.E. {round(se, 4)}')


if __name__ == '__main__':
    main()
