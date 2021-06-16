"""
Updated Oct 31, 2020: add forward model for d18O, TEX86
Updated May 6, add pH correction for d18O glassy forams

By Mingsong Li
   Penn State

Functions
    
    d18o_localsw
    
    d18o_linear
        follow Bemis et al., 1998
        
    d18oc_linear
        follow Bemis et al., 1998
        
    d18oc_linear_forward
        forward model, follow Bemis et al., 1998
    
    tex86_linear
    
    tex86_linear_forward
        forward model
        
    
    tex86h_forward
    
    mgca_anand03
    
    mgca_evans18
        Evans et al., 2018
        
    mgca_evans18_forward
        forward model
        
    mgca_sal_corr
        salinity correction
        
    mgca_sal_corr_forward
        forward model
        
    obs_estimate_r_uk37
    
    obs_estimate_r_fixed_uk37    
    
    obs_estimate_r_tex86
    
    obs_estimate_r_fixed_tex86
    
    obs_estimate_r_d18o
    
    obs_estimate_r_fixed_d18o
    
    obs_estimate_r_mgca_pooled
    
    obs_estimate_r_fixed_mgca_pooled    
    
    obs_qc
    
    proxy_frac_4da_eval
    
    cal_ye_cgenie
    
    cal_ye_cgenie_d18O
    
    cal_ye_cgenie_mgca
    
    CE_NS70
    
    rmse

"""
import numpy as np
from DeepDA_lib import modules_nc

#try:
#    import bayspline
#except ImportError as e1:
#    print('Warning:', e1)
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
    
    # T = 16.1–4.64 * (δ18Oc–δ18Osw) + 0.09 * (δ18Oc−δ18Osw)^2
    # by Bemis et al., 1998
    # Hollis et al., 2019 GMD DeepMIP protocol
    
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
    # Hollis et al., 2019 GMD DeepMIP protocol
    
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

def d18oc_linear_forward(T,d18olocalsw):
    # corrected d18oc
    # T = 16.1–4.64 * (δ18Oc–δ18Osw) + 0.09 * (δ18Oc−δ18Osw)2
    # by Bemis et al., 1998
    # Hollis et al., 2019 GMD DeepMIP protocol
    
    a = 0.09
    b = -4.64
    c = 16.1 - T
    
    sq = np.sqrt( b * b - 4 * a * c)
    d18o = (-1 * b - sq)/(2*a)
    d18oc = d18o + d18olocalsw
    return d18oc

def tex86_linear(tex86):
    '''
    TEX86 = (GDGT-2 + GDGT-3 + cren')/(GDGT-1 + GDGT-2 + GDGT-3 + cren')
    '''
    T = tex86 * 56.2 - 10.8
    return T

def tex86_linear_forward(T):
    tex86 = (T + 10.8) / 56.2
    return tex86

def tex86h(tex86):
    '''
    TEX86H = log10((GDGT-2 + GDGT-3 + cren')/(GDGT-1 + GDGT-2 + GDGT-3 + cren'))
    
    return sea surface temperature
    Ref: Kim et al., 2010
    By : Mingsong Li (Penn State, Nov 2020)
    '''
    tex86h = np.log10(tex86)
    return 68.4 * tex86h + 38.6

def tex86h_forward(sst):
    '''
    TEX86H = log10((GDGT-2 + GDGT-3 + cren')/(GDGT-1 + GDGT-2 + GDGT-3 + cren'))
    TEX86H = log10(TEX86)
    TEX86 = 10 ** TEX86H
    return TEX86
    Ref: Kim et al., 2010
    By : Mingsong Li (Penn State, Nov 2020)
    '''
    tex86h = ( sst - 38.6 )/68.4
    tex86 = 10**tex86h
    return tex86

def mgca_anand03(mgca, a, b, h, mgca_sw, mgca_swt):
    '''
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
    '''
    ratio = mgca_sw ** h / mgca_swt ** h
    MGCA = mgca / b * ratio
    T = 1 / a * np.log(MGCA)
    return T

def mgca_evans18(mgca,ph,mgcasw):
    '''
    calculate temprature for Mg/Ca proxy
    
    INPUT
        mgca  : size Nx1; mg/ca measured raw data
        ph    : size Nx1; surface seawater pH
        mgcasw: size Nx1; Mg/Casw seawater Mg/Ca ratio
    OUTPUT
        T     : size Nx1; Temperature
    
    Example #1
        mgca = np.array([3,4,5,6])
        ph = np.array([7.7,7.6,7.7,7.7])
        mgcasw = np.array([4.3,4.5,2.5,3.0])
        T = DeepDA_psm.mgca_evans18(mgca,ph,mgcasw)
    Example #2
        mgca = 3.85
        ph = 7.654
        mgcasw = 1.9216  # mean of seawater Mg/Ca at 56.0 Ma
        T = DeepDA_psm.mgca_evans18(mgca,ph,mgcasw)
        
    By Mingsong Li, Penn State, May 2, 2020 Matlab
        updated Nov 1, 2020 for python
        
    Reference:
        Evans, D., Brierley, C., Raymo, M.E., Erez, J., Müller, W., 
        2016. Planktic foraminifera shell chemistry response to seawater 
        chemistry: Pliocene–Pleistocene seawater Mg/Ca, temperature and 
        sea level change. Earth and Planetary Science Letters 438, 139-148.
      Evans, D., Sagoo, N., Renema, W., Cotton, L.J., Müller, W., Todd, J.A., 
       Saraswati, P.K., Stassen, P., Ziegler, M., Pearson, P.N., 2018. 
       Eocene greenhouse climate revealed by coupled clumped isotope-Mg/Ca 
       thermometry. Proceedings of the National Academy of Sciences, 201714744.
    
    '''
    #eq. 5 in Evans et al., 2018
    denominator = 0.76 + 0.66/(1 + np.exp(6.9 * (ph-8.0)))
    mgcanorm = mgca/denominator
    #eq. 6 in Evans et al., 2018
    B = 0.019 * mgcasw ** 2 - 0.16 * mgcasw + 0.804
    A = -0.0029 * mgcasw ** 2 + 0.032 * mgcasw
    
    t1 = 1/A
    t2 = np.log(mgcanorm/B)
    T = t1 * t2
    return T

def mgca_evans18_forward(T,ph,mgcasw):
    '''
    Mg/Ca forward model
    
    T: sst
    ph: ph
    mgcasw: Mg/Casw sea water Mg/Ca
    Ref: Evans et al., 2018
    By : Mingsong Li (Penn State, Nov 1, 2020)
    '''
    #eq. 6 in Evans et al., 2018
    B = 0.019 * mgcasw ** 2 - 0.16 * mgcasw + 0.804
    A = -0.0029 * mgcasw ** 2 + 0.032 * mgcasw
    mgcanorm = B * np.exp(A * T)
    
    denominator = 0.76 + 0.66/(1 + np.exp(6.9 * (ph-8.0)))
    
    mgca = mgcanorm * denominator
    
    return mgca
    
def mgca_sal_corr(mgca,salinity):
    '''
    Salinity correction
    
    INPUT:
      mgca: Mg/Ca_MEASURED
      salinity: salinity
      
    OUTPUT:
      Mg/Ca_CORRECTED
    Ref: Hollis et al., 2019 GMD
    By : Mingsong Li (Penn State, Nov 1, 2020)
    '''
    return (1- (salinity - 35) * 0.042) * mgca

def mgca_sal_corr_forward(mgcacorr,salinity):
    '''
    Input
        mgcacorr: corrected mg/ca
    OUTPUT
        mgca: measured Mg/Ca before salinity correction
    Ref: Hollis et al., 2019 GMD
    By : Mingsong Li (Penn State, Nov 1, 2020)
    '''
    param = 1 - (salinity - 35) * 0.042
    return mgcacorr/param
    
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
        
def proxy_frac_4da_eval(proxy_select,proxy_frac,log_level):
    from random import sample
    # INPUT
    #    proxy_select: dataframe including all sites
    #    proxy_frac: fraction of proxy data to be assimilated and evalution
    site_len = len(proxy_select)
    site_len_assim = int(site_len*proxy_frac)
    index_assim = sample(list(range(0,site_len)), site_len_assim)
    index_eval  = list(set(range(0,site_len)) - set(index_assim)) # list indices of sites not chosen
    if log_level > 1:
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
    calculate ye for caco3 and TEX86 (and d18o w/o pH correction, zeebe 2001)
    INPUT:
    
    proxies: proxy dataframe
    j: index of selected proxy
    Xb: prior, background state
    
    OUTPUT:
        precalculated ye
    '''
    # read lon lat for each line of proxy
    lon_label = yml_dict['proxies'][yml_dict['proxies']['use_from'][0]]['lon_label']
    lat_label = yml_dict['proxies'][yml_dict['proxies']['use_from'][0]]['lat_label']
    dum_lat = proxies[lat_label][j]  # (paleo)latitude of this site
    dum_lon = proxies[lon_label][j]  # (paleo)longitude of this site
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
            #print('      {}: {}'.format(key1,value1))
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
    
    ######################## TO DO: adjusted to include CaCO3, d13C or other proxies ##############
    # find 1d grid location
    if Xb.shape[0] == dum_imax * dum_jmax:
        # if size of Xb is dum_imax x dum_jmax x nens, lonlati should be <= dum_imax * dum_jmax
        lonlati = lonlat[1] * dum_jmax + lonlat[0]
    else:
        # if size of Xb is varn x dum_imax x dum_jmax x nens, lonlati should be dum_imax * dum_jmax + varn*dum_imax x dum_jmax
        lonlati = lonlat[1] * dum_jmax + lonlat[0] + psm_required_variable_key_index * dum_imax * dum_jmax
    
    #lonlati = lonlat[1] * dum_jmax + lonlat[0]
    # read prior
    prior_1grid = np.copy(Xb[lonlati,:])   # prior
    ######################## TO DO: add fit to 3d var list ##############
    
    # Now PSM type has been found. Let's cal Ye
    if proxy_psm_type_i in ['bayesreg_d18o_pooled']:
        psm_d18osw_adjust = yml_dict['psm']['bayesreg_d18o_pooled']['psm_d18osw_adjust']
        d18osw_local_choice = yml_dict['psm']['bayesreg_d18o_pooled']['d18osw_local_choice']
        d18osw_icesm_pco2 = yml_dict['psm']['bayesreg_d18o_pooled']['d18osw_icesm_pco2']
        
        if d18osw_local_choice in ['zachos94']:
            # d18o_localsw using method by Zachos et al., 1994 PALEOCEANOGRAPHY
            #d18o_localsw = DeepDA_psm.d18o_localsw(abs(dum_lat))
            x = abs(dum_lat)
            #d18o_localsw = 0.576 + 0.041 * x - 0.0017 * x ** 2 + 1.35e-5 * x ** 3
            d18o_localsw = d18o_localsw(x)
            prediction_d18O = bayfox.predict_d18oc(prior_1grid,d18o_localsw + psm_d18osw_adjust) # pool model for bayfox
        else:
            if d18osw_icesm_pco2 == 1.0:
                proxy_col_d18osw = 'd18osw_1x'
            elif d18osw_icesm_pco2 == 6.0:
                proxy_col_d18osw = 'd18osw_6x'
            elif d18osw_icesm_pco2 == 9.0:
                proxy_col_d18osw = 'd18osw_9x'
            else:
                proxy_col_d18osw = 'd18osw_3x'
            d18o_localsw = proxies[proxy_col_d18osw][j]
            prediction_d18O = bayfox.predict_d18oc(prior_1grid,d18o_localsw) # pool model for bayfox
        
        Ye = np.mean(prediction_d18O.ensemble, axis = 1)
        
    elif proxy_psm_type_i in ['deepmip_d18o']:
        psm_d18osw_adjust = yml_dict['psm']['deepmip_d18o']['psm_d18osw_adjust']
        d18osw_local_choice = yml_dict['psm']['deepmip_d18o']['d18osw_local_choice']
        d18osw_icesm_pco2 = yml_dict['psm']['deepmip_d18o']['d18osw_icesm_pco2']
        if d18osw_local_choice in ['zachos94']:
            x = abs(dum_lat)
            d18o_localsw = d18o_localsw(x)
            Ye = d18oc_linear_forward(prior_1grid,d18o_localsw + psm_d18osw_adjust)
        else:
            if d18osw_icesm_pco2 == 1.0:
                proxy_col_d18osw = 'd18osw_1x'
            elif d18osw_icesm_pco2 == 6.0:
                proxy_col_d18osw = 'd18osw_6x'
            elif d18osw_icesm_pco2 == 9.0:
                proxy_col_d18osw = 'd18osw_9x'
            else:
                proxy_col_d18osw = 'd18osw_3x'
            d18o_localsw = proxies[proxy_col_d18osw][j]
            Ye = d18oc_linear_forward(prior_1grid,d18o_localsw)
    
    elif proxy_psm_type_i in ['bayesreg_tex86']:
        # bayspar
        search_tol_i = yml_dict['psm']['bayesreg_tex86']['search_tol']
        nens_i = yml_dict['psm']['bayesreg_tex86']['nens']
        try:
            prediction = bayspar.predict_tex_analog(prior_1grid, temptype = 'sst', search_tol = search_tol_i, nens=nens_i)
        except:
            print('  bayspar Warning. search_tol may be too small. try a larger number + 10')
            prediction = bayspar.predict_tex_analog(prior_1grid, temptype = 'sst', search_tol = search_tol_i + 10, nens=nens_i)
        Ye = np.mean(prediction.ensemble, axis = 1)
        
    elif proxy_psm_type_i in ['tex86h_forward']:
        Ye = tex86h_forward(prior_1grid)
        
    elif proxy_psm_type_i in ['cgenie_caco3']:
        Ye = np.copy(prior_1grid)
        
    elif proxy_psm_type_i in ['cgenie_caco3_13c']:
        Ye = np.copy(prior_1grid)
        
    return Ye

# calculate Ye for cGENIE prior using pH correction
def cal_ye_cgenie_d18O(yml_dict,proxies,j,Xb,Xb_ph,proxy_assim2,proxy_psm_type,dum_lon_offset,dum_imax,dum_jmax):
    '''
    calculate ye for d18o w pH correction following Zeebe (2001)
    INPUT:
    
    proxies: proxy dataframe
    j: index of selected proxy
    Xb: prior, background state
    
    OUTPUT:
        precalculated ye
    '''
    # read lon lat for each line of proxy
    lon_label = yml_dict['proxies'][yml_dict['proxies']['use_from'][0]]['lon_label']
    lat_label = yml_dict['proxies'][yml_dict['proxies']['use_from'][0]]['lat_label']
    dum_lat = proxies[lat_label][j]  # (paleo)latitude of this site
    dum_lon = proxies[lon_label][j]  # (paleo)longitude of this site
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
            #print('      {}: {}'.format(key1,value1))
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
    
    ######################## TO DO: adjusted to include CaCO3, d13C or other proxies ##############
    # find 1d grid location
    if Xb.shape[0] == dum_imax * dum_jmax:
        # if size of Xb is dum_imax x dum_jmax x nens, lonlati should be <= dum_imax * dum_jmax
        lonlati = lonlat[1] * dum_jmax + lonlat[0]
    else:
        # if size of Xb is varn x dum_imax x dum_jmax x nens, lonlati should be dum_imax * dum_jmax + varn*dum_imax x dum_jmax
        lonlati = lonlat[1] * dum_jmax + lonlat[0] + psm_required_variable_key_index * dum_imax * dum_jmax
    
    #lonlati = lonlat[1] * dum_jmax + lonlat[0]
    # read prior
    prior_1grid = np.copy(Xb[lonlati,:])   # prior
    
    # read pH for d18o correction following Zeebe (2001)
    ph          =  np.copy(Xb_ph[lonlati,:])
    d18o_cor = -1.42 * (ph - 8.1)
    #print('d18o correction: ph shape is {}, mean is {}, cor factor = {}'.format(ph.shape, np.mean(ph), d18o_cor))
    ######################## TO DO: add fit to 3d var list ##############
    
    # Now PSM type has been found. Let's cal Ye
    if proxy_psm_type_i in ['bayesreg_d18o_pooled']:
        psm_d18osw_adjust = yml_dict['psm']['bayesreg_d18o_pooled']['psm_d18osw_adjust']
        d18osw_local_choice = yml_dict['psm']['bayesreg_d18o_pooled']['d18osw_local_choice']
        d18osw_icesm_pco2 = yml_dict['psm']['bayesreg_d18o_pooled']['d18osw_icesm_pco2']
        
        if d18osw_local_choice in ['zachos94']:
            # d18o_localsw using method by Zachos et al., 1994 PALEOCEANOGRAPHY
            #d18o_localsw = DeepDA_psm.d18o_localsw(abs(dum_lat))
            x = abs(dum_lat)
            #d18o_localsw = 0.576 + 0.041 * x - 0.0017 * x ** 2 + 1.35e-5 * x ** 3
            d18o_localsw = d18o_localsw(x)
            prediction_d18O = bayfox.predict_d18oc(prior_1grid,d18o_localsw + psm_d18osw_adjust) # pool model for bayfox
        else:
            if d18osw_icesm_pco2 == 1.0:
                proxy_col_d18osw = 'd18osw_1x'
            elif d18osw_icesm_pco2 == 6.0:
                proxy_col_d18osw = 'd18osw_6x'
            elif d18osw_icesm_pco2 == 9.0:
                proxy_col_d18osw = 'd18osw_9x'
            else:
                proxy_col_d18osw = 'd18osw_3x'
            d18o_localsw = proxies[proxy_col_d18osw][j]
            prediction_d18O = bayfox.predict_d18oc(prior_1grid,d18o_localsw) # pool model for bayfox
        
        Ye = np.mean(prediction_d18O.ensemble, axis = 1)
        
        #print('d18o correction: Ye shape is {}, mean is {}'.format(Ye.shape, np.mean(Ye)))
        
    elif proxy_psm_type_i in ['deepmip_d18o']:
        psm_d18osw_adjust = yml_dict['psm']['deepmip_d18o']['psm_d18osw_adjust']
        d18osw_local_choice = yml_dict['psm']['deepmip_d18o']['d18osw_local_choice']
        d18osw_icesm_pco2 = yml_dict['psm']['deepmip_d18o']['d18osw_icesm_pco2']
        if d18osw_local_choice in ['zachos94']:
            x = abs(dum_lat)
            d18o_localsw = d18o_localsw(x)
            Ye = d18oc_linear_forward(prior_1grid,d18o_localsw + psm_d18osw_adjust)
        else:
            if d18osw_icesm_pco2 == 1.0:
                proxy_col_d18osw = 'd18osw_1x'
            elif d18osw_icesm_pco2 == 6.0:
                proxy_col_d18osw = 'd18osw_6x'
            elif d18osw_icesm_pco2 == 9.0:
                proxy_col_d18osw = 'd18osw_9x'
            else:
                proxy_col_d18osw = 'd18osw_3x'
            d18o_localsw = proxies[proxy_col_d18osw][j]
            Ye = d18oc_linear_forward(prior_1grid,d18o_localsw)
        
    return Ye + d18o_cor

# calculate Ye for cGENIE prior
def cal_ye_cgenie_mgca(yml_dict,proxies,j,Xb,proxy_psm_type_i,dum_lon_offset,dum_imax,dum_jmax,Xb_sal,Xb_ph,Xb_omega,geologic_age):
    '''
    INPUT:
    OUTPUT:
        calculated ye
    '''
    proxy_assim2      = yml_dict['proxies'][yml_dict['proxies']['use_from'][0]]['proxy_assim2']
    # read lon lat for each line of proxy
    lon_label = yml_dict['proxies'][yml_dict['proxies']['use_from'][0]]['lon_label']
    lat_label = yml_dict['proxies'][yml_dict['proxies']['use_from'][0]]['lat_label']
    dum_lat = proxies[lat_label][j]  # (paleo)latitude of this site
    dum_lon = proxies[lon_label][j]  # (paleo)longitude of this site
    Filei = proxies['File'][j]
    
    # Read proxy type from the database
    data_psm_type = proxies['Proxy'][j]
    # Read allowed proxy from the DTDA-config.yml
    data_psm_type_find = 0
    for key, value in proxy_assim2.items():
        if data_psm_type in proxy_assim2[key]:
            data_psm_type_find = data_psm_type_find + 1
    
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
            #print('      {}: {}'.format(key1,value1))
            for i in range(len(prior_state_variable[key][value1])):
                if key in ['2d']:
                    prior_variable_dict.append(prior_state_variable[key][value1][i])
                    #prior_nc_file_list.append(key1+'/'+value1+'.nc')
                elif key in ['3d']:
                    prior_variable_dict_3d.append(prior_state_variable[key][value1][i])
                    #prior_nc_file_list_3d.append(key1+'/'+value1+'.nc')                    
                      
    spp = 'all'
    prior_len = Xb.shape[1]
    # ``1`` for reductive, ``0`` for BCP (Barker).
    cleaningr = np.tile(np.array([1]),prior_len)
    cleaningb = np.tile(np.array([0]),prior_len)
    lonlat = modules_nc.cal_find_ij(dum_lon,dum_lat,dum_lon_offset,dum_imax,dum_jmax) 
    ######################## TO DO: adjusted to include d13C or other proxies ##############
    
    psm_required_variable_key = list(yml_dict['psm'][proxy_psm_type_i]['psm_required_variables'].keys())[0]
    if psm_required_variable_key in prior_variable_dict:
        psm_required_variable_key_index = prior_variable_dict.index(psm_required_variable_key)
    else:
        psm_required_variable_key_index = 0
    # find 1d grid location
    if Xb.shape[0] == dum_imax * dum_jmax:
        lonlati = lonlat[1] * dum_jmax + lonlat[0]
    else:
        lonlati = lonlat[1] * dum_jmax + lonlat[0] + psm_required_variable_key_index * dum_imax * dum_jmax
    # read prior
    prior_1grid = np.copy(Xb[lonlati,:])   # prior

    salinity =  np.copy(Xb_sal[lonlati,:])
    ph       =  np.copy(Xb_ph[lonlati,:])
    if proxy_psm_type_i in ['bayesreg_mgca_pooled_red', 'bayesreg_mgca_pooled_bcp']:
        
        psm_baymag_ln =  yml_dict['psm']['bayesreg_mgca_pooled_red']['psm_baymag_ln']
        
        if proxy_psm_type_i in ['bayesreg_mgca_pooled_red']:
            clearning_one = cleaningr
            proxy_explain = 'reductive'
        elif proxy_psm_type_i in ['bayesreg_mgca_pooled_bcp']:
            clearning_one = cleaningb
            proxy_explain = 'barker'
        omega    =  np.copy(Xb_omega[lonlati,:])
        
        if psm_baymag_ln in ['no']:
            prediction_mgca = baymag.predict_mgca(prior_1grid, clearning_one, salinity, ph, omega, spp) # pool model for baymag reductive
            pred_mgca_adj = baymag.sw_correction(prediction_mgca, np.array([geologic_age]))
        if psm_baymag_ln in ['yes']:
            pred_mgca_adj = baymag.predict_mgca_ln_dt(geologic_age,prior_1grid, clearning_one, salinity, ph, omega, spp) # 
            
        Ye = np.mean(pred_mgca_adj.ensemble, axis = 1)
        
    else:
        mgcasw = yml_dict['psm'][proxy_psm_type_i]['mgcasw']
        mgcacorr = mgca_evans18_forward(prior_1grid,ph,mgcasw)
        Ye = mgca_sal_corr_forward(mgcacorr,salinity)
    return Ye

def CE_NS70(data, model, axis):
    '''
    Inputs:
        data: observation; m x n; m is value; n is time
        model: model, the same size as data
        axis: = 0 : n x m; m is value; n is time
    Outputs:
        CE: The Nash-Sutcliffe model efficiency coefficient statistic calculated following Nash & Sutcliffe (1970)
    Borrowed from LMR_utils.py by Greg Hakim & Robert Tardif, 2015
    '''
    if axis == 0:
        model = np.swapaxes(model,0,1)
        data = np.swapaxes(data,0,1)
    
    difference = model - data
    numer = np.nansum( np.power(difference,2), axis = 0 )
    denom = np.nansum( np.power(data - np.nanmean(data, axis=0), 2), axis = 0 )
    CE = 1. - np.divide(numer, denom)
        
    
    return CE

def rmse(predictions, targets):
    return np.sqrt(np.nanmean((predictions - targets) ** 2))