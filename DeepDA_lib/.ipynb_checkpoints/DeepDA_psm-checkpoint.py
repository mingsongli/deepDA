"""
"""
import numpy as np
from DeepDA_lib import modules_nc

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

def obs_estimate_r_fixed_uk37(obs):
    # return R for the observation using PSMs
    y = bayspline.predict_uk(sst=np.array([31.0]))
    return np.max(np.var(y.ensemble, axis=1))

def obs_estimate_r_tex86(obs, temptype, search_tol):
    # return R for the observation using PSMs
    y = bayspar.predict_tex_analog(seatemp=obs, temptype = temptype, search_tol = search_tol)
    return np.max(np.var(y.ensemble, axis=1))

def obs_estimate_r_fixed_tex86(obs):
    # return R for the observation using PSMs
    y = bayspar.predict_tex_analog(seatemp=np.array([1.0,15.0,31.0]), temptype = 'sst', search_tol = 15.0)
    return np.max(np.var(y.ensemble, axis=1))
    
def obs_estimate_r_d18o(obs, d18o_localsw):
    # return R for the observation using PSMs
    y = bayfox.predict_d18oc(obs,d18o_localsw) # pool model for bayfox
    return np.max(np.var(y.ensemble, axis=1))

def obs_estimate_r_fixed_d18o(obs):
    # return R for the observation using PSMs
    y = bayfox.predict_d18oc(np.array([15.0]),np.array([0.0])) # pool model for bayfox
    return np.max(np.var(y.ensemble, axis=1))
    
def obs_estimate_r_mgca_pooled(obs, cleaning, salinity, ph, omega, spp, age):
    # return R for the observation using PSMs
    prediction_mgca = baymag.predict_mgca(obs, cleaning, salinity, ph, omega, spp) # pool model for baymag reductive
    y = baymag.sw_correction(prediction_mgca, np.array([age]))
    return np.max(np.var(y.ensemble, axis=1))

def obs_estimate_r_fixed_mgca_pooled(obs, cleaning, salinity, ph, omega, spp, age):
    # return R for the observation using PSMs
    prediction_mgca = baymag.predict_mgca(np.array([1.0,15.0,16.0,31.0]), cleaning, salinity, ph, omega, spp) # pool model for baymag reductive
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
        
        #print('    obs_qc {}, ye_qc {}'.format(obs_qc, ye_qc))
        
        if obs_qc > ye_qc:
            max_qc = obs_qc
        else:
            max_qc = ye_qc
            
        if np.abs(innov) > max_qc:
            return False
        else:
            return True
        
        
# calculate Ye for cGENIE prior
def cal_ye_cgenie(yml_dict,proxies,j,Xb,proxy_assim2,proxy_psm_type,dum_lon_offset,dum_imax,dum_jmax):
    '''
    calculate ye for d18o and TEX86
    INPUT:
    
    proxies: proxy dataframe
    j: index of selected proxy
    Xb: prior, background state
    
    OUTPUT:
        precalculated ye
    '''
    # read lon lat for each line of proxy
    dum_lat = proxies['Lat'][j]  # (paleo)latitude of this site
    dum_lon = proxies['Lon'][j]  # (paleo)longitude of this site
    Filei = proxies['File'][j]    
    #yo_all[proi,:] = np.array([dum_lon, dum_lat])  # save location of this site
    
    lonlat = modules_nc.cal_find_ij(dum_lon,dum_lat,dum_lon_offset,dum_imax,dum_jmax) 
    # output [lon, lat], 
    # lon ranges from 0 (-180) to 35 (180), lat ranges from 0 (-90) to 35 (90)

    ######################## TO DO: adjusted to include d13C or other proxies ##############
    # find 1d grid location
    lonlati = lonlat[1] * dum_jmax + lonlat[0]
    # read prior
    prior_1grid = np.copy(Xb[lonlati,:])   # prior
    #print(prior_1grid.shape)
    #print(prior_1grid)
    ######################## TO DO: add  dum_ijmax * j etc. ##############
    
    # Read proxy type from the database
    data_psm_type = proxies['Proxy'][j]
    
    # Read allowed proxy from the DTDA-config.yml
    data_psm_type_find = 0
    for key, value in proxy_assim2.items():
        #print(key,value)
        # find this proxy type exist or not, how many times it occurrs
        if data_psm_type in proxy_assim2[key]:
            data_psm_type_find = data_psm_type_find + 1
    if data_psm_type_find == 1:
        for key, value in proxy_psm_type.items():
            if data_psm_type in proxy_assim2[key]:
                data_psm_key = key
        proxy_psm_type_i = proxy_psm_type[data_psm_key]
        #print('')
        #print('>>  {}. {}, grid [lon lat] {}, grid id {}'.format(j,Filei,lonlat,lonlati))
        #print('>>  PSM for {} is {}, prior mean is {}, variance is {}'.format(data_psm_type,proxy_psm_type_i, np.mean(prior_1grid), np.var(prior_1grid)))
    elif data_psm_type_find == 0:
        print('Warning, this proxy type in database is not find in DTDA-config.yml dictionary')
    else:
        print('Warning, this proxy type in database appears more than 1 time in DTDA-config.yml dictionary')
    
    # Now PSM type has been found. Let's cal Ye
    
    if proxy_psm_type_i in ['bayesreg_d18o_pooled']:
        #try:
            # bayfox
        #d18o_localsw = DeepDA_psm.d18o_localsw(abs(dum_lat))
        x = abs(dum_lat)
        d18o_localsw = 0.576 + 0.041 * x - 0.0017 * x ** 2 + 1.35e-5 * x ** 3
        psm_d18osw_adjust = yml_dict['psm']['bayesreg_d18o_pooled']['psm_d18osw_adjust']
        # total d18osw = d18o_localsw + d18o_adj + psm_d18osw_adjust
        # d18o_adj has been included in the bayfox model
        #print('>>  Prior is {}'.format(prior_1grid))
        prediction_d18O = bayfox.predict_d18oc(prior_1grid,d18o_localsw + psm_d18osw_adjust) # pool model for bayfox
        #print('>>  prediction_d18O.ensemble shape {}'.format(prediction_d18O.ensemble.shape))
        Ye = np.mean(prediction_d18O.ensemble, axis = 1)
        #yo_all[proi,:] = np.array([dum_lon, dum_lat])
        #print('>>  Ye is {}'.format(Ye))
        #print('>>   {}. Mean of Ye is {:.6f}, variance is {:.6f} '.format(proxy_psm_type_i, np.mean(Ye), np.var(Ye,ddof=1)))
        
    elif proxy_psm_type_i in ['bayesreg_tex86']:
        # bayspar
        search_tol_i = yml_dict['psm']['bayesreg_tex86']['search_tol']
        nens_i = yml_dict['psm']['bayesreg_tex86']['nens']
        try:
            prediction = bayspar.predict_tex_analog(prior_1grid, temptype = 'sst', search_tol = search_tol_i, nens=nens_i)
        except:
            print('  Warning. search_tol may be too small. try a larger number + 5')
            prediction = bayspar.predict_tex_analog(prior_1grid, temptype = 'sst', search_tol = search_tol_i + 5, nens=nens_i)
        Ye = np.mean(prediction.ensemble, axis = 1)
        #print('>>   {}. Mean of Ye is {:.6f}, variance is {:.6f} '.format(proxy_psm_type_i, np.mean(Ye), np.var(Ye,ddof=1)))
    return Ye


# calculate Ye for cGENIE prior
def cal_ye_cgenie_mgca(yml_dict,proxies,j,Xb,proxy_psm_type,dum_lon_offset,dum_imax,dum_jmax,Xb_sal,Xb_ph,Xb_omega,geologic_age):
    '''
    INPUT:
    
    OUTPUT:
        calculated ye
    '''
    # read lon lat for each line of proxy
    dum_lat = proxies['Lat'][j]  # (paleo)latitude of this site
    dum_lon = proxies['Lon'][j]  # (paleo)longitude of this site
    Filei = proxies['File'][j]    
    spp = 'all'
    prior_len = Xb.shape[1]
    # ``1`` for reductive, ``0`` for BCP (Barker).
    cleaningr = np.tile(np.array([1]),prior_len)
    cleaningb = np.tile(np.array([0]),prior_len)
    lonlat = modules_nc.cal_find_ij(dum_lon,dum_lat,dum_lon_offset,dum_imax,dum_jmax) 
    ######################## TO DO: adjusted to include d13C or other proxies ##############
    # find 1d grid location
    lonlati = lonlat[1] * dum_jmax + lonlat[0]
    # read prior
    prior_1grid = np.copy(Xb[lonlati,:])   # prior

    if proxy_psm_type in ['bayesreg_mgca_pooled_red', 'bayesreg_mgca_pooled_bcp']:
        #print('... bayesreg_mgca_pooled_red: To be done ...')
        if proxy_psm_type in ['bayesreg_mgca_pooled_red']:
            clearning_one = cleaningr
            proxy_explain = 'reductive'
        elif proxy_psm_type in ['bayesreg_mgca_pooled_bcp']:
            clearning_one = cleaningb
            proxy_explain = 'barker'

        salinity =  np.copy(Xb_sal[lonlati,:])
        ph       =  np.copy(Xb_ph[lonlati,:])
        omega    =  np.copy(Xb_omega[lonlati,:])

        prediction_mgca = baymag.predict_mgca(prior_1grid, clearning_one, salinity, ph, omega, spp) # pool model for baymag reductive
        pred_mgca_adj = baymag.sw_correction(prediction_mgca, np.array([geologic_age]))
        Ye = np.mean(pred_mgca_adj.ensemble, axis = 1)
    return Ye