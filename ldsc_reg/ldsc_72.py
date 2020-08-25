import numpy as np
from helperfuncs import *
from scipy.optimize import fmin_l_bfgs_b, minimize
from scipy.special import comb

class sibreg():
    
    def __init__(self, S, theta = None, u = None, r = None, f = None):
        
        if S.ndim > 1:
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
        if f is None:
            print("Warning: No value given for allele frequencies. Some parameters won't be noramlized.")
        
        self.theta = None if theta is None else theta[~np.any(np.isnan(theta), axis = 1)]
        self.S = S[~np.any(np.isnan(S), axis = (1, 2))]
        self.u = u[~np.isnan(u)]
        self.r = r[~np.isnan(r)]
        self.f = None if f is None else f[~np.isnan(f)]
    

    def simdata(self, V,  N, simr = False):
        
        # Simulated data (theta hats) as per section 7.1
        # V = varcov matrix of true effects
        # N = Number of obs/SNPs to generate
        
        S = self.S
        
        if simr:
            self.r = np.random.uniform(low=1, high=5, size=N)
            print("Simulated LD scores!")
        
        r = self.r

        thetahat_vec = np.empty((N, V.shape[1]))
        
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
            if d > 1:
                sim = np.random.multivariate_normal(zeromat, Si + ri * V/N)
                if i > (2/3)*N:
                    sim = np.random.multivariate_normal(zeromat, Si + ri * V)
            else:
                sim = np.random.normal(zeromat, Si + ri * V/N)
                if i > (2/3)*N:
                    sim = np.random.normal(zeromat, Si + ri * V)
            
            # Append to vector of effects
            thetahat_vec[i, :] = sim
        

        print("Effect Vectors Simulated!")
        
        self.snp = np.arange(1, N+1, 1)
        self.pos = np.arange(1, N+1, 1)
        self.theta = thetahat_vec

    def neg_logll_grad(self, V, theta = None, S = None, u = None, r = None, f = None):
        
        # ============================================ #
        # returns negative log likelihood and negative
        # of the gradient
        # ============================================ #
        
        theta = self.theta if theta is None else theta
        S = self.S if S is None else S
        u = self.u if u is None else u
        r = self.r if r is None else r
        f = self.f if f is None else f

        # Unflatten V into a matrix
        d = S[0].shape[0]
        V = return_to_symmetric(V, d)
        Gvec = np.zeros((d, d))
        
        N = len(S)
        log_ll = 0
        
        # Normalizing variables
        # V = V * N
        V_norm = V/N
        for i in range(N):
            
            Si = S[i]
            thetai = theta[i, :]
            ui = u[i]
            ri = r[i]
            
            
            fi = f[i]  if f is not None else None

            d, ddash = Si.shape
            assert d == ddash # Each S has to be a square matrix
            
            # normalizing variables using allele frequency
            normalizer = 2 * fi  * (1 - fi) if fi is not None else 1.0
            thetai = np.sqrt(normalizer) * thetai
            Si = normalizer * Si
      
            # calculate log likelihood
            log_ll_add = -(d/2) * np.log(2 * np.pi)
            dit_sv = np.linalg.det(Si + ri * V_norm)
            dit_sv = 1e-6 if dit_sv < 0 else dit_sv
            log_ll_add += -(1/2) * np.log(dit_sv)
            log_ll_add += -(1/2) * np.trace(np.outer(thetai, thetai) @ np.linalg.inv(Si + ri * V_norm))
            log_ll_add *= 1/ui
            
            if np.isnan(log_ll_add) == False:
                log_ll += log_ll_add
            
            # calculate gradient
            SV_inv = np.linalg.inv(Si + ri * V_norm)
            G = -(1 / 2) * SV_inv
            G += (1 / 2) * np.dot(SV_inv,np.dot(np.outer(thetai, thetai),SV_inv))
            G *= 1/ui
            
            if np.any(np.isnan(G)) == False:
                Gvec += G

        Gvec = extract_upper_triangle(Gvec)
        return -log_ll, -Gvec


    def solve(self,
              theta = None, 
              S = None,
              u = None,
              r = None,
              f = None,
              neg_logll_grad = None,
              est_init = None,
              printout = True):
        
        # inherit parameters from the class if they aren't defined
        theta = self.theta if (theta is None) else theta
        S = self.S if S is None else S
        u = self.u if u is None else u
        r = self.r if r is None else r
        f = self.f if f is None else f
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
                    print("Making 'optimal' matrix")
                    print("=================================================")
                
                theta_full = theta
                S_full = S/n
                
                theta_var = np.cov(theta_full.T)
                S_hat = np.mean(S_full, axis = 0)
                est_init = n * (theta_var - S_hat)
        else:
            if printout == True:
                print("No initial guess provided.")
                print("Making 'optimal' matrix")
                print("=================================================")
            
            theta_full = theta
            S_full = S/n
            
            theta_var = np.cov(theta_full.T)
            S_hat = np.mean(S_full, axis = 0)
            est_init = n * (theta_var - S_hat)
            
        
        # exporting for potential later reference
        self.est_init = est_init

        # extract array from est init
        est_init_array = extract_upper_triangle(est_init) 
        
        bounds = extract_bounds(m)

        result = minimize(
            neg_logll_grad, 
            est_init_array,
            jac = True,
            args = (theta, S, u, r, f),
            bounds = bounds,
            method = 'L-BFGS-B'
        )
        
        output_matrix = return_to_symmetric(result.x, m)
        
        # re-normnalizing output matrix
        output_matrix = output_matrix / n
        
        self.output_matrix = output_matrix
        
        return output_matrix, result 

    def jackknife_se(self,
                  theta  = None, S = None,
                  r = None, u = None,
                  blocksize = 1):

        # Simple jackknife estimator for SE
        # Ref: https://github.com/bulik/ldsc/blob/aa33296abac9569a6422ee6ba7eb4b902422cc74/ldscore/jackknife.py#L231
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
                                              printout = False,
                                              est_init = self.est_init)

                estimates_jk.append(output_matrix)

                start_idx += blocksize
            else:
                break
            
        estimates_jk = np.array(estimates_jk)
        full_est = self.output_matrix
        
        # calculate pseudo-values
        n_blocks = nobs/blocksize
        pseudovalues = n_blocks * full_est - (n_blocks - 1) * estimates_jk
        
        # calculate jackknife se
        pseudovalues = pseudovalues.reshape((nobs, theta.shape[1] * theta.shape[1]))
        jknife_cov = np.cov(pseudovalues.T, ddof=1) / n_blocks
        jknife_var = np.diag(jknife_cov)
        jknife_se = np.sqrt(jknife_var)
    
        jknife_se  = jknife_se.reshape((theta.shape[1], theta.shape[1]))
        
        return jknife_se  




