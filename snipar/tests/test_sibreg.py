import unittest
import numpy as np
from numpy import testing
from snipar import lmm
from snipar.tests.utils import *

def random_design(labels):
    unique_labels = np.unique(labels)
    n_labels = unique_labels.shape[0]
    label_dict = {}
    for i in range(0,n_labels):
        label_dict[unique_labels[i]] = i
    Z = np.zeros((labels.shape[0],n_labels))
    for i in range(0,labels.shape[0]):
        Z[i,label_dict[labels[i]]]=1
    return Z

def Sigma_make(labels,sigma2,tau):
    Z = random_design(labels)
    sigmau = sigma2/tau
    return sigmau * Z.dot(Z.T) + sigma2 * np.identity(Z.shape[0])

def safe_likelihood(y,X,Sigma):
    alpha = safe_alpha_mle(y,X,Sigma)
    alpha = alpha.reshape((alpha.shape[0],))
    slogdet = np.linalg.slogdet(Sigma)
    slogdet = slogdet[0]*slogdet[1]
    resid = y-X.dot(alpha)
    Sigma_inv = np.linalg.inv(Sigma)
    L = slogdet+np.dot(resid.T.dot(Sigma_inv),resid)
    return L

def safe_alpha_mle(y,X,Sigma):
    Sigma_inv = np.linalg.inv(Sigma)
    X_T_Sigma_inv = np.dot(X.T,Sigma_inv)
    X_T_X = np.dot(X_T_Sigma_inv,X)
    X_T_y = np.dot(X_T_Sigma_inv,y)
    return np.linalg.solve(X_T_X,X_T_y)

class test_regrnd_functions(SniparTest):

    def test_alpha_mle(self):
        sigma2 = float(1)
        sigmau = float(5.5)
        n = 10 ** 2
        for i in range(0, 100):
            alpha = np.random.randn((2))
            tau = sigma2 / sigmau
            m = lmm.simulate(n, alpha, sigma2, tau)
            Z = random_design(m.labels)
            # code.interact(local=locals())
            Sigma = sigmau * Z.dot(Z.T) + sigma2 * np.identity(n)
            safe_alpha = safe_alpha_mle(m.y, m.X, Sigma)
            alpha_mle = m.alpha_mle(tau,sigma2)
            testing.assert_almost_equal(alpha_mle, safe_alpha, decimal=5)

    def test_likelihood(self):
        sigma2 = float(1)
        sigmau = float(5.5)
        n = 10 ** 2
        for i in range(0,100):
            alpha=np.random.randn((2))
            tau = sigma2/sigmau
            m = lmm.simulate(n,alpha,sigma2,tau)
            Z = random_design(m.labels)
            #code.interact(local=locals())
            Sigma = sigmau*Z.dot(Z.T)+sigma2*np.identity(n)
            safe_lik=safe_likelihood(m.y,m.X,Sigma)
            lik, grad = m.likelihood_and_gradient(sigma2,tau)
            testing.assert_almost_equal(lik,safe_lik/float(n),decimal=5)

    def test_grad_sigma2(self):
        n = 10 ** 3
        c = 2
        sigma2 = float(1)
        sigmau = float(10)
        for i in range(0, 100):
            alpha = np.random.randn((c))
            tau = sigma2 / sigmau
            m = lmm.simulate(n, alpha, sigma2, tau)
            lik, grad = m.likelihood_and_gradient(sigma2, tau)

            # Numerical gradient
            # Compute gradient numerically
            def likelihood(sigma2):
                lik, grad = m.likelihood_and_gradient(sigma2, tau)
                return lik
            num_grad = (likelihood(sigma2 + 10**(-6)) - likelihood(sigma2 - 10**(-6))) / (2 * 10 ** (-6))
            testing.assert_almost_equal(grad[0], num_grad, decimal=5)

    def test_grad_tau(self):
        n = 10 ** 3
        c = 2
        sigma2 = float(10)
        sigmau = float(100)
        for i in range(0, 100):
            alpha = np.random.randn((c))
            tau = sigma2 / sigmau
            m = lmm.simulate(n, alpha, sigma2, tau)
            lik, grad = m.likelihood_and_gradient(sigma2, tau)

            # Numerical gradient
            # Compute gradient numerically
            def likelihood(tau):
                lik, grad = m.likelihood_and_gradient(sigma2, tau)
                return lik
            num_grad = (likelihood(tau + 10**(-6)) - likelihood(tau - 10**(-6))) / (2 * 10 ** (-6))
            testing.assert_almost_equal(grad[1], num_grad, decimal=5)


if  __name__=='__main__':
    unittest.main()