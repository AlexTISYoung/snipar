import numpy as np
import numdifftools as nd

def delta_method(f, x, v_x):
    """Delta method for estimating standard error of a function of a vector of parameters

    Args:
        f (_type_): function of vector of parameters
        x (_type_): vector of parameters
        v_x (_type_): variance-covariance matrix of parameters

    Returns:
        _type_: _description_
    """
    if np.isnan(x).any():
        return np.nan
    else:
        grad_f = nd.Gradient(f)(x)
        v_f = grad_f@v_x@grad_f.T
        return np.sqrt(v_f) 

def estimate_k(delta, se_delta, h2f, se_h2f, rk, bound=True):
    """Estimate k, the fraction of random-mating heritability explained by PGI

    Args:
        delta (_type_): direct effect estimate 
        se_delta (_type_): standard error of direct effect estimate
        h2f (_type_): family-based heritability estimate
        se_h2f (_type_): standard error of family-based heritability estimate
        rk (_type_): estimated correlation between parents for PGI

    Returns:
        _type_: _description_
    """
    if h2f<=0:
        return np.nan
    zf_inv = se_h2f/h2f
    k =(1-rk)*(1-zf_inv**(2))*(delta**2-se_delta**2)/h2f
    if k<0 and bound:
        return 0
    elif k>1 and bound:
        return 1
    else:
        return k

def k_se(delta, se_delta, h2f, se_h2f, rk, se_rk):
    """Standard error of estimate of k, the fraction of random-mating heritability explained by PGI

    Args:
        delta (_type_): direct effect estimate 
        se_delta (_type_): standard error of direct effect estimate
        h2f (_type_): family-based heritability estimate
        se_h2f (_type_): standard error of family-based heritability estimate
        rk (_type_): estimated correlation between parents for PGI

    Returns:
        _type_: _description_
    """
    def est_k(x):
        return estimate_k(x[0], se_delta, x[1], se_h2f, x[2], bound=False)
    return delta_method(est_k, np.array([delta, h2f, rk]), np.diag([se_delta**2, se_h2f**2, se_rk**2]))

def estimate_r(delta, se_delta, h2f, se_h2f, rk, return_all=True, bound=True, allow_neg_r=False):
    k = estimate_k(delta, se_delta, h2f, se_h2f, rk, bound=bound)
    r = rk/(k+(1-k)*rk)
    if r>1 and bound:
        r=1
    elif r<(-1) and bound:
        r=-1
    return {'r':r, 'k':k} if return_all else r

def r_se(delta, se_delta, h2f, se_h2f, rk, se_rk):
    def est_r(x):
        return estimate_r(x[0], se_delta, x[1], se_h2f, x[2],return_all=False, bound=False)
    return delta_method(est_r, np.array([delta, h2f, rk]), np.diag([se_delta**2, se_h2f**2, se_rk**2]))

def estimate_h2eq(delta, se_delta, h2f, se_h2f, rk, return_all=True):
    r = estimate_r(delta, se_delta, h2f, se_h2f, rk)
    if r['r'] >= 1:
        h2eq = np.nan
    else: 
        h2eq = h2f/(1-r['r'])
    return {'h2eq':h2eq, 'r':r['r'], 'k':r['k']} if return_all else h2eq

def h2eq_se(delta, se_delta, h2f, se_h2f, rk, se_rk):
    def est_h2eq(x):
        return estimate_h2eq(x[0], se_delta, x[1], se_h2f, x[2],return_all=False)
    return delta_method(est_h2eq, np.array([delta, h2f, rk]), np.diag([se_delta**2, se_h2f**2, se_rk**2]))

def estimate_rho(delta, se_delta, h2f, se_h2f, rk, return_all=True, bound=True):
    r = estimate_r(delta, se_delta, h2f, se_h2f, rk)
    rho = 1-(1-r['k'])*r['r']
    return {'rho':rho, 'r':r['r'], 'k':r['k']} if return_all else rho

def rho_se(delta, se_delta, h2f, se_h2f, rk, se_rk):
    def est_rho(x):
        return estimate_rho(x[0], se_delta, x[1], se_h2f, x[2],return_all=False)
    return delta_method(est_rho, np.array([delta, h2f, rk]), np.diag([se_delta**2, se_h2f**2, se_rk**2]))

def alpha_from_alpha(delta, se_delta, alpha, h2f, se_h2f, rk, return_all=True):
    rho = estimate_rho(delta, se_delta, h2f, se_h2f, rk)
    alpha = (rho['rho']+rho['k']*rho['r'])*(alpha/delta)-(1-rho['rho'])
    alpha = alpha/(1+rho['r'])
    if return_all:
        return {'alpha_delta':alpha, 'rho':rho['rho'], 'r':rho['r'], 'k':rho['k']}
    else:
        return alpha

def se_alpha_from_alpha(delta, se_delta, alpha, se_alpha, r_delta_alpha, h2f, se_h2f, rk, se_rk):
    def est_alpha(x):
        return alpha_from_alpha(x[0], se_delta, x[1], x[2], se_h2f, x[3], return_all=False)
    v_x = np.diag([se_delta**2, se_alpha**2, se_h2f**2, se_rk**2])
    v_x[0,1] = r_delta_alpha*se_delta*se_alpha
    v_x[1,0] = v_x[0,1]
    return delta_method(est_alpha, np.array([delta, alpha, h2f, rk]), v_x)

def alpha_from_beta(delta, se_delta, beta, h2f, se_h2f, rk, return_all=True):
    rho = estimate_rho(delta, se_delta, h2f, se_h2f, rk)
    alpha = (rho['rho']-delta/beta)/((delta/beta)*(1+rho['r']))
    if return_all:
        return {'alpha_delta':alpha, 'rho':rho['rho'], 'r':rho['r'], 'k':rho['k']}
    else:
        return alpha

def se_alpha_from_beta(delta, se_delta, beta, se_beta, r_delta_beta, h2f, se_h2f, rk, se_rk):
    def est_alpha(x):
        return alpha_from_beta(x[0], se_delta, x[1], x[2], se_h2f, x[3], return_all=False)
    v_x = np.diag([se_delta**2, se_beta**2, se_h2f**2, se_rk**2])
    v_x[0,1] = r_delta_beta*se_delta*se_beta
    v_x[1,0] = v_x[0,1]
    return delta_method(est_alpha, np.array([delta, beta, h2f, rk]), v_x)

def v_eta_delta(delta, se_delta, h2f, se_h2f, rk, alpha=None, beta=None, return_all=True):
    if alpha is not None:
        alpha_all = alpha_from_alpha(delta, se_delta, alpha, h2f, se_h2f, rk, return_all=True)
    elif beta is not None:
        alpha_all = alpha_from_beta(delta, se_delta, beta, h2f, se_h2f, rk, return_all=True)
    else:
        raise(ValueError('Must provide average NTC or population effect'))
    if alpha_all['r'] >= 1:
        h2_eq, v_ed = np.nan, np.nan
    else:
        h2_eq = h2f/(1-alpha_all['r'])
        v_ed = 2*(1+alpha_all['r'])*alpha_all['alpha_delta']*(1+alpha_all['alpha_delta'])*h2_eq
    if return_all == True:
        return {'v_eta_delta':v_ed, 
                'alpha_delta':alpha_all['alpha_delta'], 
                'rho':alpha_all['rho'], 
                'r':alpha_all['r'], 
                'k':alpha_all['k'],
                'h2_eq':h2_eq}
    else:
        return v_ed

def se_v_eta_delta(delta, se_delta, ab, se_ab, r_delta_ab, h2f, se_h2f, rk, se_rk, is_beta=False):
    if is_beta:
        def est_v_ed(x):
            return v_eta_delta(x[0], se_delta, x[1], se_h2f, x[2], beta=x[3], return_all=False)
    else:
        def est_v_ed(x):
            return v_eta_delta(x[0], se_delta, x[1], se_h2f, x[2], alpha=x[3], return_all=False)
    v_x = np.diag([se_delta**2, se_h2f**2, se_rk**2, se_ab**2])
    v_x[0,3] = r_delta_ab*se_delta*se_ab
    v_x[3,0] = v_x[0,3]
    return delta_method(est_v_ed, np.array([delta, h2f, rk, ab]), v_x)

# Simulate sample estimates
def simulate_ests(n, delta, alpha, beta, v_alpha_delta_beta, h2f, h2f_se, rk, rk_se):
    # Simulate delta, alpha, beta
    delta_alpha = np.random.multivariate_normal([delta, alpha, beta], v_alpha_delta_beta, n)
    delta_se = np.sqrt(v_alpha_delta_beta[0,0])
    alpha_se = np.sqrt(v_alpha_delta_beta[1,1])
    beta_se = np.sqrt(v_alpha_delta_beta[2,2])
    r_delta_alpha = v_alpha_delta_beta[0,1]/(delta_se*alpha_se)
    r_delta_beta = v_alpha_delta_beta[0,2]/(delta_se*beta_se)
    # Simulate h2f
    h2f = h2f+h2f_se*np.random.randn(n)
    # Simulate rk
    rk = rk+rk_se*np.random.randn(n)
    # Estimates
    k_ests = np.zeros(n)
    k_ests[:] = np.nan
    r_ests = np.zeros(n)
    r_ests[:] = np.nan
    h2eq_ests = np.zeros(n)
    h2eq_ests[:] = np.nan
    rho_ests = np.zeros(n)
    rho_ests[:] = np.nan
    alpha_from_beta_ests = np.zeros(n)
    alpha_from_beta_ests[:] = np.nan
    alpha_from_alpha_ests = np.zeros(n)
    alpha_from_alpha_ests[:] = np.nan
    v_eta_delta_ests = np.zeros(n)
    v_eta_delta_ests[:] = np.nan
    # Standard error estimates
    se_k_ests = np.zeros(n)
    se_k_ests[:] = np.nan
    se_r_ests = np.zeros(n)
    se_r_ests[:] = np.nan
    se_h2eq_ests = np.zeros(n)
    se_h2eq_ests[:] = np.nan
    se_rho_ests = np.zeros(n)
    se_rho_ests[:] = np.nan
    se_alpha_from_beta_ests = np.zeros(n)
    se_alpha_from_alpha_ests = np.zeros(n)
    se_alpha_from_alpha_ests = np.zeros(n)
    se_alpha_from_alpha_ests[:] = np.nan
    se_v_eta_delta_ests = np.zeros(n)
    se_v_eta_delta_ests[:] = np.nan
    for i in range(n):
        if h2f[i]>0:
            # Sample estimates
            k_ests[i] = estimate_k(delta_alpha[i,0], delta_se, h2f[i], h2f_se, rk[i])
            r_ests[i] = estimate_r(delta_alpha[i,0], delta_se, h2f[i], h2f_se, rk[i], return_all=False)
            h2eq_ests[i] = estimate_h2eq(delta_alpha[i,0], delta_se, h2f[i], h2f_se, rk[i], return_all=False)
            rho_ests[i] = estimate_rho(delta_alpha[i,0], delta_se, h2f[i], h2f_se, rk[i], return_all=False)
            alpha_from_beta_ests[i] = alpha_from_beta(delta_alpha[i,0], delta_se, delta_alpha[i,2], h2f[i], h2f_se, rk[i], return_all=False)
            alpha_from_alpha_ests[i] = alpha_from_alpha(delta_alpha[i,0], delta_se, delta_alpha[i,1], h2f[i], h2f_se, rk[i], return_all=False)
            v_eta_delta_ests[i] = v_eta_delta(delta_alpha[i,0], delta_se, h2f[i], h2f_se, rk[i], alpha=delta_alpha[i,1], return_all=False)
            # Standard error estimates      
            se_k_ests[i] = k_se(delta_alpha[i,0], delta_se, h2f[i], h2f_se, rk[i], rk_se)
            se_r_ests[i] = r_se(delta_alpha[i,0], delta_se, h2f[i], h2f_se, rk[i], rk_se)
            se_h2eq_ests[i] = h2eq_se(delta_alpha[i,0], delta_se, h2f[i], h2f_se, rk[i], rk_se)
            se_rho_ests[i] = rho_se(delta_alpha[i,0], delta_se, h2f[i], h2f_se, rk[i], rk_se)
            se_alpha_from_beta_ests[i] = se_alpha_from_beta(delta_alpha[i,0], delta_se, delta_alpha[i,2], beta_se, r_delta_beta, h2f[i], h2f_se, rk[i], rk_se)
            se_alpha_from_alpha_ests[i] = se_alpha_from_alpha(delta_alpha[i,0], delta_se, delta_alpha[i,1], alpha_se, r_delta_alpha, h2f[i], h2f_se, rk[i], rk_se)
            se_v_eta_delta_ests[i] = se_v_eta_delta(delta_alpha[i,0], delta_se, delta_alpha[i,1], alpha_se, r_delta_alpha, h2f[i], h2f_se, rk[i], rk_se, is_beta=False)
    return k_ests, se_k_ests, r_ests, se_r_ests, h2eq_ests, se_h2eq_ests, rho_ests, se_rho_ests, alpha_from_beta_ests, se_alpha_from_beta_ests, alpha_from_alpha_ests, se_alpha_from_alpha_ests, v_eta_delta_ests, se_v_eta_delta_ests

def se_check(ests,se_ests):
    return np.nanmedian(se_ests)/np.nanstd(ests)

def check_se_calc(n, delta, delta_se, alpha, alpha_se, h2f, h2f_se, rk, rk_se):
    # Simulate sample estimates
    k_ests, se_k_ests, r_ests, se_r_ests, h2eq_ests, se_h2eq_ests, rho_ests, se_rho_ests, alpha_from_beta_ests, se_alpha_from_beta_ests, alpha_from_alpha_ests, se_alpha_from_alpha_ests, v_eta_delta_ests, se_v_eta_delta_ests = simulate_ests(n, delta, delta_se, alpha, alpha_se, h2f, h2f_se, rk, rk_se)
    # Check standard error estimates
    return np.nanmean(k_ests), se_check(k_ests,se_k_ests), np.nanmean(r_ests), se_check(r_ests,se_r_ests), np.nanmean(h2eq_ests), se_check(h2eq_ests,se_h2eq_ests), np.nanmean(rho_ests), se_check(rho_ests,se_rho_ests), np.nanmean(alpha_from_beta_ests), se_check(alpha_from_beta_ests,se_alpha_from_beta_ests), np.nanmean(alpha_from_alpha_ests), se_check(alpha_from_alpha_ests,se_alpha_from_alpha_ests), np.nanmean(v_eta_delta_ests), se_check(v_eta_delta_ests,se_v_eta_delta_ests)
    
def simulation_check(n, h2_eq, r_y, k, N_fam, h2f_se):
    # Set parameters
    r_delta = r_y*h2_eq
    print('r_delta: '+str(round(r_delta,3)))
    h2f = (1-r_delta)*h2_eq
    print('h2f: '+str(round(h2f,3)))
    rho = 1-(1-k)*r_delta
    print('rho: '+str(round(rho,3)))
    rk = k*r_delta/rho
    rk_se = np.sqrt((1-rk**2)**2/(N_fam-2))
    print('rk: '+str(round(rk,3)))
    delta = np.sqrt(k*h2f/(1-rk))
    print('delta: '+str(round(delta,3)))
    beta = delta/rho
    print('beta: '+str(round(beta,3)))
    alpha = (beta-delta)/(1+rk)
    print('alpha: '+str(round(alpha,3)))
    # Calculate standard errors
    r2 = delta**2+2*(1+rk)*alpha*(alpha+delta)
    v_alpha_delta = ((1-r2)/(1-rk**2))*np.array([[2*(1+rk),-(1+rk)],[-(1+rk),1]])/N_fam
    A = np.array([[1,0],[0,1],[1,(1+rk)]])
    v_alpha_delta_beta = A @ v_alpha_delta @ A.T
    # Run simulation
    k_mean, k_se_ratio, r_mean, r_se_ratio, h2_eq_mean, h2_se_ratio, rho_mean, rho_se_ratio, alpha_1_mean, alpha_1_se_ratio, alpha_2_mean, alpha_2_se_ratio, mean_v_eta_delta, v_eta_delta_se_ratio = check_se_calc(n, delta, alpha, beta, v_alpha_delta_beta, h2f, h2f_se, rk, rk_se)
    return r_y, N_fam, delta, np.sqrt(v_alpha_delta[0,0]), alpha, np.sqrt(v_alpha_delta[1,1]), rk, rk_se, h2f, h2f_se, k, k_mean, k_se_ratio, r_delta, r_mean, r_se_ratio, h2_eq, h2_eq_mean, h2_se_ratio, rho, rho_mean, rho_se_ratio, alpha_1_mean, alpha_1_se_ratio, alpha_2_mean, alpha_2_se_ratio, mean_v_eta_delta, v_eta_delta_se_ratio