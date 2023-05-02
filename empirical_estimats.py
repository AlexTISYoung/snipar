import snipar.pgs as pgs

def r_from_ratio_se(x1,x2,se1,se2,se_ratio):
    return ((x1*x2)/(2*se1*se2))*(se1**2/x1**2+se2**2/x2**2-se_ratio**2*(x2/x1)**2)

### Height ###
rk = 0.106
rk_se = 0.020
delta = 0.531
delta_se = 0.0056
beta = 0.5835
beta_se = 0.0043 
ratio_se = 0.0094
r_delta_beta = r_from_ratio_se(delta, beta, delta_se, beta_se, ratio_se)
## RDR ##
h2f = 0.554
h2f_se = 0.044
height_rdr_ests, height_rdr_ses = pgs.am_adj_2gen_calc(delta, delta_se, beta, beta_se, r_delta_beta, h2f, h2f_se, rk, rk_se, is_beta=True)
## Twin ##
h2f = 0.729
h2f_se = 0.018
height_twin_ests, height_twin_ses = pgs.am_adj_2gen_calc(delta, delta_se, beta, beta_se, r_delta_beta, h2f, h2f_se, rk, rk_se, is_beta=True)

### EA ###
rk = 0.175
rk_se = 0.020
delta = 0.1833
delta_se = 0.0073
beta = 0.330
beta_se = 0.0044 
r_delta_beta = 0.391
## Twin ##
h2f = 0.400
h2f_se = 0.024
ea_twin_ests, ea_twin_ses = pgs.am_adj_2gen_calc(delta, delta_se, beta, beta_se, r_delta_beta, h2f, h2f_se, rk, rk_se, is_beta=True)