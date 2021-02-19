import unittest
import numpy as np
from numpy import testing
import ldsc_reg.inferz.sib_ldsc_z as ld
import scipy.stats

def get_logll_scipy(V, z, r, S):

    '''
    Gets log likelihood function from scipy
    Used to test if our log likelihood function 
    is defined correctly
    '''
    
    S_inv_root = ld.calc_inv_root(S)
    dist = scipy.stats.multivariate_normal(mean = None,
                                          cov = np.eye(V.shape[0]) + r * S_inv_root @ V @ S_inv_root)

    nlogll = dist.logpdf(z)
    return nlogll
    
    
class test_regrnd_functions(unittest.TestCase):
    
    def test_logll(self):
        '''
        Testing if log likelihood function is 
        proper
        '''
    
        # Simulating data
        np.random.seed(123)

        N = int(100)
        S_size = int(N/2)
        S = np.array([np.array([[.5, 0], [0, .8]]),
            np.array([[0.5, 0], [0, 0.8]])] * S_size )
        V = np.identity(2) * 0.5
        f = np.random.uniform(0, 1, N)
        
        model = ld.sibreg(S = S/N)
        model.simdata(V/N, N, simr = True)
        
        scipy_logll = np.empty(N)
        a_logll = np.empty(N)
        
        for i in range(N):
            scipy_logll[i] = get_logll_scipy(V, model.z[i, :], model.r[i], model.S[i])
            a_logll[i] = model._log_ll(V, model.z[i, :], model.S[i], model.r[i])
            
        
        np.allclose(scipy_logll, a_logll)
    
    def test_grad_logll(self):
    
        '''
        Testing if gradient function is proper
        '''
        
        # Simulating data
        np.random.seed(123)

        N = int(100)
        S_size = int(N/2)
        S = np.array([np.array([[.5, 0], [0, .8]]),
            np.array([[0.5, 0], [0, 0.8]])] * S_size )
        V = np.identity(2) * 0.5
        f = np.random.uniform(0, 1, N)
        
        model = ld.sibreg(S = S/N)
        model.simdata(V/N, N, simr = True)
        
        aderiv = np.empty((N, 2, 2))
        nderiv = np.empty((N, 2, 2))
        
        for i in range(N):
            aderiv[i, :, :] = model._grad_ll_v(V, model.z[i, :], model.S[i],  model.r[i])
            nderiv[i, :, :] = model._num_grad_V(V, model.z[i, :], model.S[i],  model.r[i])

        np.allclose(aderiv, nderiv)
        
        
if  __name__=='__main__':
    unittest.main()
    