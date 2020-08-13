import numpy as np
from helperfuncs import *
from scipy.optimize import fmin_l_bfgs_b
from scipy.special import comb

class sibreg_72():
    
    def __init__(self, S, theta = None, u = None, r = None):

        for s in S:
            n, m = s.shape
            assert n == m

        if theta is None:
            print("Warning there is no value for theta. Maybe consider simulating it")
        if u is None:
            print("No value for U given. Generating a vector of ones (all SNPs weighted equally)")
            u = np.ones(S.shape[0])
        if r is None:
            print("No value for r given. Generating a vector of ones for r")
            r = np.ones(S.shape[0])

        self.theta = theta
        self.S = S
        self.u = u
        self.r = r
        

    def simdata(self, V, r, N):
        
        # Simulated data (theta hats) as per section 7.1
        # V = varcov matrix of true effects
        # N = Number of obs/SNPs to generate
        
        S = self.S
        theta = self.theta
        r = self.r

        thetahat_vec = []
        
        # make sure they are np arrays
        for i in range(N):
            
            Si = S[i]
            ri = r[i]
            
            V = np.array(V)
            Si = np.array(Si)

            # get shape of V
            d = V.shape[0]
            zeromat = np.zeros(d)

            # generate true effect vector
            theta = np.random.multivariate_normal(zeromat, ri * V)

            sim = np.random.multivariate_normal(theta, Si)
            
            # Append to vector of effects
            thetahat_vec.append(sim)
        
        thetahat_vec = np.array(thetahat_vec)

        print("Effect Vectors Simulated!")
        
        self.theta = thetahat_vec

    def neg_logll_grad(self, V, theta = None, S = None, u = None, r = None):
        
        # ============================================ #
        # returns negative log likelihood and negative
        # of the gradient
        # ============================================ #
        
        theta = self.theta if theta is None else theta
        S = self.S if S is None else S
        u = self.u if u is None else u
        r = self.r if r is None else r

        # Unflatten V into a matrix
        d = S[0].shape[0]
        V = return_to_symmetric(V, d)
        Gvec = np.zeros((d, d))
        
        N = len(S)
        log_ll = 0
        
        for i in range(N):
            
        
            Si = S[i]
            thetai = theta[i, :]
            ui = u[i]
            ri = r[i]

            d, ddash = Si.shape
            assert d == ddash # Each S has to be a square matrix
      
            # calculate log likelihood
            log_ll += -(d/2) * np.log(2 * np.pi)
            dit_sv = np.linalg.det(Si + ri * V)
            log_ll += -(1/2) * np.log(dit_sv)
            log_ll += -(1/2) * np.trace(np.outer(thetai, thetai) @ np.linalg.inv(Si + ri * V))
            log_ll *= 1/ui
            
            
            # calculate gradient
            SV_inv = np.linalg.inv(Si + ri * V)
            G = -(1 / 2) * SV_inv
            G += (1 / 2) * np.dot(SV_inv,np.dot(np.outer(thetai, thetai),SV_inv))
            G *= 1/ui
            
            Gvec += G

        Gvec = extract_upper_triangle(Gvec)
        print("Log Likelihood: ", log_ll)
        print("Gradient: ", Gvec)
        return -log_ll, -Gvec


    def solve(self,
              theta = None, 
              S = None,
              u = None,
              r = None,
              neg_logll_grad = None,
              est_init = None,
              printout = True):
        
        # inherit parameters from the class if they aren't defined
        theta = self.theta if (theta is None) else theta
        S = self.S if S is None else S
        u = self.u if u is None else u
        r = self.r if r is None else r
        neg_logll_grad = self.neg_logll_grad if neg_logll_grad is None else neg_logll_grad

        # == Solves our MLE problem == #
        n, m = theta.shape
        
        if est_init is not None:
            # Shape of initial varcov guess
            rowstrue = est_init.shape[0] == m
            colstrue = est_init.shape[1] == m

            if rowstrue & colstrue:
                pass
            else:
                if printout == True:
                    print("Warning: Initial Estimate given is not of the proper dimension")
                    print("Making a matrix of 0s as the initial estimate")
                    print("=================================================")
                    
                est_init = np.zeros((m, m))
        else:
            if printout == True:
                print("No initial guess provided.")
                print("Making a matrix of 0s as the initial estimate")
                print("=================================================")
            
            est_init = np.zeros((m, m))
            
        
        # extract array from est init
        est_init_array = extract_upper_triangle(est_init) 
        
        bounds = extract_bounds(m)

        result = fmin_l_bfgs_b(
            neg_logll_grad, 
            est_init_array,
            fprime = None,
            args = (theta, S, u, r),
            bounds = bounds
        )
        
        output_matrix = return_to_symmetric(result[0], m)
        
        if printout == True:
            print("Final Estimate:\n", output_matrix)
            print("Convergence Flag: ", result[2]['task'])
            print("Number of Iterations: ", result[2]['nit'])
            print("Final Gradient: ", result[2]['grad'])
        
        return output_matrix, result 

    def jackknife_se(self,
                  theta  = None, S = None,
                  r = None, u = None,
                  blocksize = 1):

        # Simple jackknife estimator for SE
        # Source: https://www.stat.berkeley.edu/~hhuang/STAT152/Jackknife-Bootstrap.pdf
        # Default value of blocksize = 1 is the normal jackknife

        theta = self.theta if (theta is None) else theta
        S = self.S if (S is None) else S
        r = self.r if (r is None) else r
        u = self.u if (u is None) else u

        assert theta.shape[0] == S.shape[0]

        nobs = theta.shape[0]
        
        estimates_jk = []
        
        start_idx = 0
        while True:
            
            end_idx = start_idx + blocksize
            end_idx_cond = end_idx <= theta.shape[0]
            
            # remove blocks of observations

            vars_jk = []

            for var in [theta, S, r, u]:

                var_jk = delete_obs_jk(var, start_idx, end_idx,
                                       end_idx_cond)
                vars_jk.append(var_jk)
            
            if start_idx < theta.shape[0]:
                # Get our estimate
                output_matrix, _ = self.solve(theta = vars_jk[0],
                                              S = vars_jk[1],
                                              r = vars_jk[2],
                                              u = vars_jk[3],
                                              printout = False)

                estimates_jk.append(output_matrix)

                start_idx += 1
            else:
                break
            
        estimates_jk = np.array(estimates_jk)
        
        # calculate the SE
        estimate_jk_mean = estimates_jk.mean(axis = 0)
        estimate_jk_mean = np.array([estimate_jk_mean] * estimates_jk.shape[0])
        estimates_jk_dev = estimates_jk - estimate_jk_mean
        estimates_jk_devsq = estimates_jk_dev ** 2
        estimates_jk_devsq_sum = estimates_jk_devsq.sum(axis = 0)
        
        se_correction = (nobs - blocksize)/(blocksize *comb(nobs, blocksize))

        se = se_correction * estimates_jk_devsq_sum
        se = np.sqrt(se)
        
        return se  




