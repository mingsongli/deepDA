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
    testn = 20 # run 20 times, output mean of these 20 runs
    outn = np.full(20,np.nan)
    for i in range(testn):
        prediction_mgca = baymag.predict_mgca(np.array([1.0,15.0,16.0,31.0]), cleaning, salinity, ph, omega, spp) # pool model for baymag reductive
        y = baymag.sw_correction(prediction_mgca, np.array([age]))
        outn[i] = np.max(np.var(y.ensemble, axis=1))
    return np.nanmean(outn)
    
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
        
def proxy_frac_4da_eval(proxy_select,proxy_frac):
    from random import sample
    # INPUT
    #    proxy_select: dataframe including all sites
    #    proxy_frac: fraction of proxy data to be assimilated and evalution
    site_len = len(proxy_select)
    site_len_assim = int(site_len*proxy_frac)
    index_assim = sample(list(range(0,site_len)), site_len_assim)
    index_eval  = list(set(range(0,site_len)) - set(index_assim)) # list indices of sites not chosen
    print('>>  Selected index: {}'.format(index_assim))
    print('>>  Unselected index: {}'.format(index_eval))
    sites_eval = []
    assim_i = 0
    eval_i = 0
    for j in range(len(index_assim)):
        if assim_i == 0:
            sites_assim = proxy_select.iloc[[index_assim[j]]]
            sites_assim = sites_assim.reset_index(drop=True) # reset_index, avoid index error
            assim_i = 1
        else:
            sites_assim = sites_assim.append(proxy_select.iloc[[index_assim[j]]], ignore_index=True)
    for k in range(len(index_eval)):
        if eval_i == 0:
            sites_eval = proxy_select.iloc[[index_eval[k]]]
            sites_eval = sites_eval.reset_index(drop=True) # reset_index, avoid index error
            eval_i = 1
        else:
            sites_eval = sites_eval.append(proxy_select.iloc[[index_eval[k]]], ignore_index=True)
    #sites_assim = [proxy_select.iloc[[j]] for j in index_assim]
    #sites_eval  = [proxy_select.iloc[[j]] for j in index_eval]
    #sites_eval  = [proxy_select[p] for p in index_eval]
    return sites_assim, sites_eval

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
    # Read proxy type from the database
    data_psm_type = proxies['Proxy'][j]
    # Read allowed proxy from the DTDA-config.yml
    data_psm_type_find = 0
    for key, value in proxy_assim2.items():
        if data_psm_type in proxy_assim2[key]:
            data_psm_type_find = data_psm_type_find + 1
    if data_psm_type_find == 1:
        for key, value in proxy_psm_type.items():
            if data_psm_type in proxy_assim2[key]:
                data_psm_key = key
        proxy_psm_type_i = proxy_psm_type[data_psm_key]
    
    ###### Read Prior dic ####
    # save prior variable list
    prior_variable_dict = []  # variable list
    prior_nc_file_list = []  # nc file list
    prior_variable_dict_3d = []  # variable list
    prior_nc_file_list_3d = []  # nc file list
    prior_source = yml_dict['prior']['prior_source'] #
    prior_state_variable = yml_dict['prior'][prior_source]['state_variable']  # note: ['2d': xxx; '3d': xxx]
    for key, value in prior_state_variable.items():
        nc_keyvalue = prior_state_variable[key]['ncname']  # note: 2d dict
        #print('      nc_keyvalue {}...'.format(nc_keyvalue))
        for key1, value1 in nc_keyvalue.items():
            print('      {}: {}'.format(key1,value1))
            for i in range(len(prior_state_variable[key][value1])):
                if key in ['2d']:
                    prior_variable_dict.append(prior_state_variable[key][value1][i])
                    #prior_nc_file_list.append(key1+'/'+value1+'.nc')
                elif key in ['3d']:
                    prior_variable_dict_3d.append(prior_state_variable[key][value1][i])
                    #prior_nc_file_list_3d.append(key1+'/'+value1+'.nc')

                
    psm_required_variable_key = list(yml_dict['psm'][proxy_psm_type_i]['psm_required_variables'].keys())[0]
    if psm_required_variable_key in prior_variable_dict:
        psm_required_variable_key_index = prior_variable_dict.index(psm_required_variable_key)
    else:
        psm_required_variable_key_index = 0
    
    lonlat = modules_nc.cal_find_ij(dum_lon,dum_lat,dum_lon_offset,dum_imax,dum_jmax) 
    ######################## TO DO: adjusted to include d13C or other proxies ##############
    # find 1d grid location
    lonlati = lonlat[1] * dum_jmax + lonlat[0] + psm_required_variable_key_index * dum_imax * dum_jmax
    # read prior
    prior_1grid = np.copy(Xb[lonlati,:])   # prior
    ######################## TO DO: add  dum_ijmax * j etc. ##############
    
    # Now PSM type has been found. Let's cal Ye
    if proxy_psm_type_i in ['bayesreg_d18o_pooled']:
        x = abs(dum_lat)
        d18o_localsw = 0.576 + 0.041 * x - 0.0017 * x ** 2 + 1.35e-5 * x ** 3
        psm_d18osw_adjust = yml_dict['psm']['bayesreg_d18o_pooled']['psm_d18osw_adjust']
        prediction_d18O = bayfox.predict_d18oc(prior_1grid,d18o_localsw + psm_d18osw_adjust) # pool model for bayfox
        Ye = np.mean(prediction_d18O.ensemble, axis = 1)
        
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
        
    elif proxy_psm_type_i in ['cgenie_caco3']:
        Ye = np.copy(prior_1grid)
        
    elif proxy_psm_type_i in ['cgenie_caco3_13c']:
        Ye = np.copy(prior_1grid)
        
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

def CE_NS70(data, model):
    '''
    Inputs:
        data: observation
        model: model, the same size as data
    Outputs:
    CE: CE statistic calculated following Nash & Sutcliffe (1970)
    Borrowed from LMR_utils.py by Greg Hakim & Robert Tardif, 2015
    '''
    difference = model - data
    numer = np.nansum( np.power(difference,2), axis = 0 )
    denom = np.nansum( np.power(data - np.nanmean(data, axis=0),2), axis = 0 )
    CE = 1. - np.divide(numer, denom)
    
    return CE