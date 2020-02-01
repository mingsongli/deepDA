"""
"""
import numpy as np

try:
    import bayspline
except ImportError as e1:
    print('Warning:', e1)
try:
    import bayspar
except ImportError as e2:
    print('Warning:', e2)
try:
    import bayfox
except ImportError as e3:
    print('Warning:', e3)
try:
    import baymag
except ImportError as e4:
    print('Warning:', e4)

def d18o_localsw(x):
    '''
    local sea water d18o
    x: latitude
    '''
    d18osw = 0.576 + 0.041 * x - 0.0017 * x ** 2 + 1.35e-5 * x ** 3
    return d18osw

def d18o_linear(d18oc,d18olocalsw,poly):
    
    # T = 16.1–4.64 * (δ18Oc–δ18Osw) + 0.09 * (δ18Oc−δ18Osw)2
    # by Bemis et al., 1998
    
    alpha = 16.1
    beta  = -4.64
    gamma = 0.09 
    T = np.nan
    if poly == 2:
        T = alpha + beta * (d18oc - d18olocalsw) + gamma * (d18oc - d18olocalsw) ** 2
    elif poly == 1:
        T = alpha + beta * (d18oc - d18olocalsw)
    else:
        print('poly must be either 1 or 2')
    return T

def d18oc_linear(d18oc,poly):
    # corrected d18oc
    # T = 16.1–4.64 * (δ18Oc–δ18Osw) + 0.09 * (δ18Oc−δ18Osw)2
    # by Bemis et al., 1998
    
    alpha = 16.1
    beta  = -4.64
    gamma = 0.09 
    T = np.nan
    if poly == 2:
        T = alpha + beta * (d18oc) + gamma * (d18oc) ** 2
    elif poly == 1:
        T = alpha + beta * (d18oc)
    else:
        print('poly must be either 1 or 2')
    return T

def tex86_linear(tex86):
    result = tex86 * 56.2 - 10.8
    return result


def mgca_anand03(mgca, a, b, h, mgca_sw, mgca_swt):
    
    #    Based on Dunkley Jones et al., 2013 ESR
    #    T = 1/a * ln( MGCA )
    #    MGCA = mgca / b * mgca_sw ^ h / mgca_swt ^ h
    #, where
    #    mgca = Mg/Ca ratio in foraminiferal calcite
    #    a = 0.09, standard Mg/Catemperature calibration constants
    #    b = 0.38, standard Mg/Catemperature calibration constants
    #    h = 0.42, a power law dependence of test Mg/Ca with changing Mg/Casw with power component, h 
    #            (Evans and Müller, 2012)
    #    mgca_sw = 5.15 mol/mol, Mg/Ca ratios of modern seawater
    #    mgca_swt, Mg/Ca ratios of ancient seawater, mgca_swt = 2 mol/mol for the latest Paleocene (Dunkley Jones et al., 2013)

    ratio = mgca_sw ** h / mgca_swt ** h
    MGCA = mgca / b * ratio
    T = 1 / a * np.log(MGCA)
    
    return T

def obs_estimate_r_uk37(obs):
    # return R for the observation using PSMs
    y = bayspline.predict_uk(sst=obs)
    return np.max(np.var(y.ensemble, axis=1))

def obs_estimate_r_tex86(obs, temptype, search_tol):
    # return R for the observation using PSMs
    y = bayspar.predict_tex_analog(seatemp=obs, temptype = temptype, search_tol = search_tol)
    return np.max(np.var(y.ensemble, axis=1))
    
def obs_estimate_r_d18o(obs, d18o_localsw):
    # return R for the observation using PSMs
    y = bayfox.predict_d18oc(obs,d18o_localsw) # pool model for bayfox
    return np.max(np.var(y.ensemble, axis=1))
    
def obs_estimate_r_mgca_pooled(obs, cleaning, salinity, ph, omega, spp, age):
    # return R for the observation using PSMs
    prediction_mgca = baymag.predict_mgca(obs, cleaning, salinity, ph, omega, spp) # pool model for baymag reductive
    y = baymag.sw_correction(prediction_mgca, np.array([age]))
    return np.max(np.var(y.ensemble, axis=1))
    
def obs_qc(Ye, obs, obs_err, proxy_qc):
    # quality control
    # Ye: ye, an x by 1 ensemble of PSM-generated data
    # obs: observation value, single data
    # obs_err:  observation variance, single data
    # proxy_qc: settings
    innov = obs - np.mean(Ye)
    if not proxy_qc:
        return True
    else:
        obs_qc = np.sqrt( obs_err ) * proxy_qc
        ye_qc = np.std(Ye) * proxy_qc
        
        print('obs_qc {}, ye_qc {}'.format(obs_qc, ye_qc))
        
        if obs_qc > ye_qc:
            max_qc = obs_qc
        else:
            max_qc = ye_qc
            
        if np.abs(innov) > max_qc:
            return False
        else:
            return True