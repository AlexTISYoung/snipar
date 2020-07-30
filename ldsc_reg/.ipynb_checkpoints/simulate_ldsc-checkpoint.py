import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("whitegrid")
from scipy.special import factorial
from scipy.optimize import minimize, fmin_l_bfgs_b
from numba import jit, njit, prange
import pprint


def vparams2Vvec(Vparams):
    # function for transforming
    # parameters into appropriate matrix
    # V = np.array([[Vparams[0], Vparams[1]],
    #               [Vparams[1], Vparams[2]]])

    V = Vparams.reshape((2, 2), order = 'F')
    return V


def gvec2gparam(gvec):
    # converts gradient vector
    # to conform with dimensions of V

    # G = np.array([gvec[0, 0], gvec[0, 1], gvec[1, 1]])
    G = gvec.flatten(order = "F")

    return G


def simdata(V, S, N):
    # Simulated data (theta hats) as per section 7.1
    # V = varcov matrix of true effects
    # S = array of variance covariance matrices (each one
    # for a given snp)
    # N = Number of obs/SNPs to generate
    # Make sure S has as man
    θhat_vec = []
    # make sure they are np arrays
    for i in range(N):
        Si = S[i]
        V = np.array(V)
        Si = np.array(Si)
        # get shape of V
        d = V.shape[0]
        zeromat = np.zeros(d)
        # generate true effect vector
        θ = np.random.multivariate_normal(zeromat, V)
        sim = np.random.multivariate_normal(θ, Si)
        # Append to vector of effects
        θhat_vec.append(sim)
    θhat_vec = np.array(θhat_vec)
    return θhat_vec

N = 1000
S = np.array([np.array([[5, 0], [0, 5]])]*N)# 50 = N/2
V = np.identity(2) * 10.0


hat_vec = simdata(V, S, N)

def logll(V, θ, S):
    # calculate negative log likelihood
    # Unflatten V into a matrix
    d = S[0].shape[0]
    V = vparams2Vvec(V)
    N = len(S)
    log_ll = 0
    for i in range(N):
        Si = S[i]
        θi = θ[i]
        d, ddash = Si.shape
        assert d == ddash  # Each S has to be a square matrix
        log_ll += -(d / 2) * np.log(2 * np.pi)
        log_ll += -(1 / 2) * np.log(np.linalg.det(Si + V))
        log_ll += -(1 / 2) * np.trace(np.outer(θi, θi) @ np.linalg.inv(Si + V))
    return -log_ll

Vvec = np.array([10., 0., 0., 10.])
logll(Vvec, hat_vec, S)


def grad_logll(V, θ, S):
    # the gradient of the log
    # likelihood function
    # Unflatten V into a matrix
    d = S[0].shape[0]
    V = vparams2Vvec(V)
    N = len(S)
    Gvec = np.zeros((d, d))
    for i in range(N):
        Si = S[i]
        θi = θ[i, :]
        SV_inv = np.linalg.inv(Si + V)
        G = -(1 / 2) * SV_inv
        G += (1 / 2) * np.dot(SV_inv,np.dot(np.outer(θi, θi),SV_inv))
        Gvec += G
    #Gvec = Gvec.flatten(order = 'F')
    Gvec = gvec2gparam(Gvec)
    return -Gvec


def neg_logll_grad(V, θ, S):
    # returns negative log likelihood and negative
    # of the gradient
    # Unflatten V into a matrix
    d = S[0].shape[0]
    V = vparams2Vvec(V)
    Gvec = np.zeros((d, d))
    N = len(S)
    log_ll = 0
    for i in range(N):
        Si = S[i]
        θi = θ[i]
        d, ddash = Si.shape
        assert d == ddash  # Each S has to be a square matrix
        SV_inv = np.linalg.inv(Si + V)
        # calculate log likelihood
        log_ll += -(d / 2) * np.log(2 * np.pi)
        log_ll += -(1 / 2) * np.log(np.linalg.det(Si + V))
        log_ll += -(1 / 2) * np.trace(np.outer(θi, θi) @ SV_inv)
        # calculate gradient
        G = -(1 / 2) * np.linalg.inv(Si + V)
        G += (1 / 2) * np.dot(SV_inv,np.dot(np.outer(θi, θi),SV_inv))
        Gvec += G
    #     Gvec = Gvec.flatten(order = 'F')
    Gvec = gvec2gparam(Gvec)
    return -log_ll, -Gvec

Vinit = np.array([0., 0., 0., 0.])

result = fmin_l_bfgs_b(neg_logll_grad, Vinit,
            fprime = None,
            args = (hat_vec, S),
            bounds = [(1e-5, None),
                    (None, None),
                    (None, None),
                    (1e-5, None)])

# Print Results
print("=============================")
print("Estimated Results:\n")
pprint.pprint(vparams2Vvec(result[0]))
print("=============================")
print("True Vector:\n")
pprint.pprint(vparams2Vvec(Vvec))
