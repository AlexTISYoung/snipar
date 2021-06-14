import unittest
import numpy as np
from numpy import testing
import ldsc_reg.sib_ldsc_z as ld
import scipy.stats

def get_logll_scipy(V, z, l, S, M):

    '''
    Gets log likelihood function from scipy
    Used to test if our log likelihood function 
    is defined correctly
    '''
    
    V = ld.V2Vmat(V, M)

    Vnew, Snew = ld.standardize_mat(V, S, M)

    dist = scipy.stats.multivariate_normal(mean = None,
                                          cov = Snew + l * Vnew)

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
        model.simdata(V/N, N, simld = True)
        
        scipy_logll = np.empty(N)
        a_logll = np.empty(N)
        
        for i in range(N):
            scipy_logll[i] = get_logll_scipy(ld.Vmat2V(V, model.M), model.z[i, :], model.l[i], model.S[i], model.M)
            a_logll[i] = ld._log_ll(ld.Vmat2V(V, model.M), model.z[i, :], model.S[i], model.l[i], model.M)
            
        
        self.assertTrue(np.allclose(scipy_logll, a_logll))
    
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
        model.simdata(V/N, N, simld = True)
        
        aderiv = np.empty((N, 3))
        nderiv = np.empty((N, 3))
        
        for i in range(N):
            aderiv[i, :] = ld._grad_ll_v(ld.Vmat2V(V, model.M), model.z[i, :], model.S[i],  model.l[i], model.M)
            nderiv[i, :] = ld._num_grad_V(ld.Vmat2V(V, model.M), model.z[i, :], model.S[i], model.l[i], model.M)

        self.assertTrue(np.allclose(aderiv, nderiv))

    def test_modified_chol(self):

        A = np.random.rand(3,3)
        G = A.T @ A

        Gmod = ld.modified_cholesky(G)

        self.assertTrue(np.allclose(G, Gmod))

        
if  __name__=='__main__':
    np.random.seed(123)
    unittest.main()
    