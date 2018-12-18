import unittest
import numpy as np
from numpy import testing
import sibreg


def random_design(labels):
    unique_labels = np.unique(labels)
    n_labels = unique_labels.shape[0]
    label_dict = {}
    for i in xrange(0,n_labels):
        label_dict[unique_labels[i]] = i
    Z = np.zeros((labels.shape[0],n_labels))
    for i in xrange(0,labels.shape[0]):
        Z[i,label_dict[labels[i]]]=1
    return Z

def safe_likelihood(y,X,Sigma,alpha,l):
    slogdet = np.linalg.slogdet(Sigma)
    slogdet = slogdet[0]*slogdet[1]
    resid = y-X.dot(alpha)
    Sigma_inv = np.linalg.inv(Sigma)
    L = slogdet+np.dot(resid.T.dot(Sigma_inv),resid)+l*np.sum(np.square(alpha))
    return L

class test_regrnd_functions(unittest.TestCase):

    def test_likelihood(self):
        sigma2 = float(1)
        sigmau = float(5.5)
        n = 10 ** 2
        l = 10
        for i in xrange(0,100):
            alpha=np.random.randn((2))
            tau = sigma2/sigmau
            m = regrnd.simulate(n,alpha,sigma2,tau)
            Z = random_design(m.labels)
            #code.interact(local=locals())
            Sigma = sigmau*Z.dot(Z.T)+sigma2*np.identity(n)
            safe_lik=safe_likelihood(m.y,m.X,Sigma,alpha,l)
            lik, grad = m.likelihood_and_gradient(alpha,sigma2,tau,l)
            testing.assert_almost_equal(lik,safe_lik/float(n),decimal=5)

    def test_grad_alpha(self):
        n = 10 ** 2
        l = 1
        c = 2
        sigma2 = float(10)
        sigmau = float(1)
        for i in xrange(0,100):
            alpha = np.random.randn((c))
            tau = sigma2 / sigmau
            m = regrnd.simulate(n, alpha, sigma2, tau)
            lik, grad = m.likelihood_and_gradient(alpha, sigma2, tau, l)
            # Numerical gradient
            # Compute gradient numerically
            def likelihood(alpha):
                lik, grad = m.likelihood_and_gradient(alpha, sigma2,tau,l)
                return lik
            num_grad = np.zeros((c))
            diffs = np.identity(c) * 10 ** (-6)
            for i in xrange(0, c):
                num_grad[i] = (likelihood(alpha + diffs[i, :]) - likelihood(alpha - diffs[i, :])) / (2 * 10 ** (-6))
            testing.assert_almost_equal(grad[0:c], num_grad , decimal=5)

    def test_grad_sigma2(self):
        n = 10 ** 3
        l = 1
        c = 2
        sigma2 = float(1)
        sigmau = float(10)
        for i in xrange(0, 100):
            alpha = np.random.randn((c))
            tau = sigma2 / sigmau
            m = regrnd.simulate(n, alpha, sigma2, tau)
            lik, grad = m.likelihood_and_gradient(alpha, sigma2, tau, l)

            # Numerical gradient
            # Compute gradient numerically
            def likelihood(sigma2):
                lik, grad = m.likelihood_and_gradient(alpha, sigma2, tau, l)
                return lik
            num_grad = (likelihood(sigma2 + 10**(-6)) - likelihood(sigma2 - 10**(-6))) / (2 * 10 ** (-6))
            testing.assert_almost_equal(grad[c], num_grad, decimal=5)

    def test_grad_tau(self):
        n = 10 ** 3
        l = 1
        c = 2
        sigma2 = float(10)
        sigmau = float(100)
        for i in xrange(0, 100):
            alpha = np.random.randn((c))
            tau = sigma2 / sigmau
            m = regrnd.simulate(n, alpha, sigma2, tau)
            lik, grad = m.likelihood_and_gradient(alpha, sigma2, tau, l)

            # Numerical gradient
            # Compute gradient numerically
            def likelihood(tau):
                lik, grad = m.likelihood_and_gradient(alpha, sigma2, tau, l)
                return lik
            num_grad = (likelihood(tau + 10**(-6)) - likelihood(tau - 10**(-6))) / (2 * 10 ** (-6))
            testing.assert_almost_equal(grad[c+1], num_grad, decimal=5)


if  __name__=='__main__':
    unittest.main()