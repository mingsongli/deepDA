'''
Data assimilation for deep time
Stage 1:    Prior: cGENIE only
            Proxy: petmproxy3slices format database
            PSM: bayesian proxy system model
            DA: Mingsong Li, with LMR DA Core
            
            Mingsong Li
            1/15/2020
Stage 2:    Proxy confirmed: TEX86, 
            Updated: Feb. 10, 2020
            
'''
# Package
import h5py
from DeepDA_lib import LMR_DA
from DeepDA_lib import modules_nc
from DeepDA_lib import DeepDA_psm

from netCDF4 import Dataset
import os
import numpy as np
import numpy.ma as ma
import numpy.matlib as mat
import scipy.stats as stats
import pandas
from sys import platform as sys_pf
import yaml
import matplotlib.pyplot as plt
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("TkAgg")
    import matplotlib.pyplot as plt
%matplotlib inline
from mpl_toolkits.basemap import Basemap, shiftgrid, cm

print('>>  OKAY.')

print('>>  Load config and prepare prior')
config_name = "DeepDA_config.yml"
#config_name = "petmproxy3slices_v0.0.10gt1.csvexp_petm78_og1_qc_obs_20200203_test2.yml"

# read DTDA-config.yml
f = open(config_name, 'r')
yml_dict = yaml.load(f, Loader=yaml.FullLoader)
f.close()

########## Prior #########
prior_source = yml_dict['prior']['prior_source'] #
prior_state_variable = yml_dict['prior'][prior_source]['state_variable']  # note: ['2d': xxx; '3d': xxx]
locRad = yml_dict['core']['local_rad'] #
dum_lon_offset = yml_dict['prior'][prior_source]['dum_lon_offset'] # longitude offset

# save prior variable list
prior_variable_dict = []  # variable list
prior_nc_file_list = []  # nc file list
prior_variable_dict_3d = []  # variable list
prior_nc_file_list_3d = []  # nc file list

for key, value in prior_state_variable.items():
    nc_keyvalue = prior_state_variable[key]['ncname']  # note: 2d dict
    
    print('>>  nc_keyvalue {}...'.format(nc_keyvalue))
    for key1, value1 in nc_keyvalue.items():
        print('>>  {}: {}'.format(key1,value1))
        
        for i in range(len(prior_state_variable[key][value1])):
            if key in ['2d']:
                prior_variable_dict.append(prior_state_variable[key][value1][i])
                prior_nc_file_list.append(key1+'/'+value1+'.nc')
            elif key in ['3d']:
                prior_variable_dict_3d.append(prior_state_variable[key][value1][i])
                prior_nc_file_list_3d.append(key1+'/'+value1+'.nc')
                
# variable list
prior_variable_len = len(prior_variable_dict)
prior_variable3d_len = len(prior_variable_dict_3d)
print('>>  Number of prior variables is: {}. List:'.format(prior_variable_len))
print('      {}'.format(prior_variable_dict))

dir_prior = yml_dict['core']['prior_dir']
dir_prior_full = os.listdir(dir_prior)
try:
    #x0 = Dataset(dir_prior+'/'+dir_prior_full[0]+'/'+ nc_file_2d).variables[prior_variable_dict[0]][0,:,:]
    x1 = Dataset(dir_prior+'/'+dir_prior_full[0]+'/'+ prior_nc_file_list_3d[0]).variables[prior_variable_dict_3d[0]][0,:,:,:]
    zt = Dataset(dir_prior+'/'+dir_prior_full[0]+'/'+ prior_nc_file_list_3d[0]).variables['zt'][:]
    print('    Shape of prior 3d nc file {}'.format(x1.shape))
    #print(zt)
    dum_dmax = x1.shape[0] # depth
    dum_imax = x1.shape[1]  # lon
    dum_jmax = x1.shape[2]  # lat
except:
    try:
        x0 = Dataset(dir_prior+'/'+dir_prior_full[0]+'/'+ prior_nc_file_list[0]).variables[prior_variable_dict[0]][0,:,:]
        dum_imax = x0.shape[0]  # lon
        dum_jmax = x0.shape[1]  # lat
        dum_dmax = 16
        print('    Shape of prior 2d nc file {}'.format(x0.shape))
    except:
        dum_dmax = 16
        dum_imax = 36
        dum_jmax = 36
# prepare 2d Xb for lon-lat state 
dum_ijmax = dum_imax*dum_jmax  # lonn * latn
print('>>  Shape of dum_dmax {}, dum_imax {}, dum_jmax {}, dum_ijmax {}'.format(dum_dmax,dum_imax,dum_jmax,dum_ijmax))
######## 

nexp = yml_dict['core']['nexp']
nens = yml_dict['core']['nens']
dir_data_save = yml_dict['core']['wrkdir']
recon_period = yml_dict['core']['recon_period']
recon_timescale = yml_dict['core']['recon_timescale_interval']
recon_period_full = np.arange(recon_period[0],recon_period[1]+1,recon_timescale)
recon_period_len = recon_period_full.shape[0]
save_ens_full = yml_dict['core']['save_ens_full']
proxy_assim2 = yml_dict['proxies'][yml_dict['proxies']['use_from'][0]]['proxy_assim2']
proxy_psm_type    = yml_dict['proxies'][yml_dict['proxies']['use_from'][0]]['proxy_psm_type']
proxy_frac      = yml_dict['proxies']['proxy_frac']

# NetCDF file name
nc_filename = dir_data_save + '/' + yml_dict['proxies'][yml_dict['proxies']['use_from'][0]]['dbversion'] + '.' + nexp + '.nc'
# read preprior HDF5 file
dir_proxy_data = dir_data_save +'/'+ yml_dict['proxies'][yml_dict['proxies']['use_from'][0]]['dbversion']
proxy_err_eval = yml_dict['proxies'][yml_dict['proxies']['use_from'][0]]['proxy_err_eval']
hdf5name = dir_proxy_data + '.' + nexp + '_precal_ye.hdf5'
kcov_saving = 0

# for saving DA product Xa
if prior_variable_len > 0:
    Xa_output   = np.full((dum_ijmax * prior_variable_len, nens, recon_period_len),np.nan)
    Xa_output_all = Xa_output
    if prior_variable3d_len > 0:
        Xa3d_output   = np.full((dum_ijmax * dum_dmax * prior_variable_len, nens, recon_period_len),np.nan)
        Xa_output_all = np.concatenate((Xa_output, Xa3d_output), axis=0)
    else:
        print('>>  No 3d variable listed in {}'.format(config_name))
elif prior_variable_len == 0:
    if prior_variable3d_len > 0:
        Xa3d_output   = np.full((dum_ijmax * dum_dmax * prior_variable_len, nens, recon_period_len),np.nan)
        Xa_output_all = Xa3d_output
    print('>>  No 2d variable listed in {}'.format(config_name))
else:
    print('>>  Error! No 3d or 2d variables are listed in {}'.format(config_name))


# ========= dataset for plot =========
cGENIEGrid = yml_dict['core']['proj_dir'] + '/data_misc/cGENIEGrid.csv'
cGENIEGrid = pandas.read_csv(cGENIEGrid)

#print(cGENIEGrid)
cGENIEGridB_lat36 = cGENIEGrid['lat']
cGENIEGridB_lon36 = cGENIEGrid['lon']
#print(cGENIEGridB_lat36.shape)
cGENIEGrid = cGENIEGrid.to_numpy()
#for i in range(36):
#        if cGENIEGrid[i,1] < 0:
#            cGENIEGrid[i,1] = 360 + cGENIEGrid[i,1]
#cGENIEGrid = cGENIEGrid.reshape((dum_imax,2))
#print(cGENIEGrid)
print('>>  OKAY.')


# DA core script

proxies=pandas.read_hdf(hdf5name, 'proxies')
proxy_psm_type_dict_df = pandas.read_hdf(hdf5name, 'proxy_psm_type_dict_df')
proxy_psm_type_dict_list = proxy_psm_type_dict_df[0].values.tolist()

with h5py.File(hdf5name, 'r') as f:
    Xb = f.get('Xb')  # read Xb, background 2d field data
    Xb3d = f.get('Xb3d')  # read Xb, background 3d field data order: lon-lat-depth
    if Xb and Xb3d:
        Xball = np.concatenate((Xb, Xb3d), axis=0)
    elif Xb and Xb3d is None:
        Xball = Xb
    elif Xb is None and Xb3d:
        Xball = Xb3d
    else:
        print('>>  Error! No 3d or 2d variables are listed in {}'.format(config_name))
    
    Xb0 = np.copy(Xball)  # default Xb
    obvalue_full = f.get('obvalue')
    Ye_full = f.get('Ye')
    ob_err_full = f.get('ob_err')
    ob_err0_full = f.get('ob_err0')
    ob_err_comb = f.get('ob_err_comb')
    yo_all = f.get('yo_all')  # read location data
    
    if 'bayesreg_mgca_pooled_bcp' in proxy_psm_type_dict_list or 'bayesreg_mgca_pooled_red' in proxy_psm_type_dict_list:
        Xb_sal = f.get('Xb_sal')
        Xb_omega = f.get('Xb_omega')
        Xb_ph = f.get('Xb_ph')
        geologic_age = yml_dict['core']['geologic_age']
        print('>>  Mg/Ca proxy found. Load salinity, pH and omega')
    
    Xa_output_all = np.full((Xball.shape[0], Xball.shape[1], recon_period_len),np.nan)
    ob_len = obvalue_full.shape[0]
    
    print('>>  recon intervals: {}, obser number {}'.format(recon_period_len,ob_len))
    for reconi in range(recon_period_len):
        Xball = Xb0.copy()  # initialize Xball
        for obi in range(ob_len):
            print('>>  Recon ID: {}, obser ID {}'.format(reconi,obi))
            yo_loc = yo_all[obi,:]  # read location
            obvalue  = obvalue_full[obi, reconi]  # read observation value
            if proxy_err_eval in ['proxy_err_psm_mp']:
                ob_err = ob_err_comb[obi, reconi]  # read observation error, use PSM model + interval data uncertainty
            else:
                ob_err = ob_err0_full[obi, reconi] # read observation error, use PSM model only
                
            # proxy type
            proxy_psm_type_i = proxy_psm_type_dict_df[0][obi]
            if proxy_psm_type_i in ['bayesreg_tex86', 'bayesreg_d18o_pooled', 'cgenie_caco3', 'cgenie_caco3_13c']:
                Ye = DeepDA_psm.cal_ye_cgenie(yml_dict,proxies,obi,Xball,proxy_assim2,proxy_psm_type,dum_lon_offset,dum_imax,dum_jmax)
            elif proxy_psm_type_i in ['bayesreg_mgca_pooled_bcp', 'bayesreg_mgca_pooled_red']:
                Ye = DeepDA_psm.cal_ye_cgenie_mgca(yml_dict,proxies,obi,Xball,proxy_psm_type_i,dum_lon_offset,dum_imax,dum_jmax,Xb_sal,Xb_ph,Xb_omega,geologic_age)
            
            if ~np.isnan(obvalue) and ~np.isnan(ob_err_comb[obi, reconi]):
                print('>>                           Loc: {}. Mean of Ye {:.6f}, var {:.6f}, obs {:.6f}, obs_err {:.6f}'.format(yo_loc,np.mean(Ye),np.var(Ye,ddof=1), obvalue, ob_err))
                if locRad:
                    covloc = modules_nc.covloc_eval(locRad, yo_loc, dum_jmax, dum_imax, cGENIEGrid)
                    covlocext = int(Xball.shape[0] / covloc.shape[0])
                    covloc = np.matlib.repmat(covloc, covlocext, 1).reshape((Xball.shape[0],))
                else:
                    covloc = np.full((Xball.shape[0],),1)
                #print('>>  Shape of Xball {}, ye {}, ob_err {}, covloc {}'.format(Xball.shape, Ye.shape, ob_err.shape, covloc.shape))
                Xa = LMR_DA.enkf_update_array(Xball, obvalue, Ye, ob_err, loc = covloc)
                #XaMean = np.ma.MaskedArray(Xa, np.matlib.repmat(np.copy(xbm) >= 9.9692e+36, 150,1))
                #print('>>    mean of Xa is {}'.format(np.nanmean(Xa)))
                
                if reconi == 0 and obi == 0:
                    kcov_saving = 1
                    ye = np.subtract(Ye, np.mean(Ye))
                    xbm = np.mean(Xball,axis=1)
                    Xbp = np.subtract(Xball,xbm[:,None])  # "None" means replicate in this dimension
                    kcov = np.dot(Xbp,np.transpose(ye)) / (nens-1)
                # update Xb using Xa, to assimilate next observation
                Xball = np.copy(Xa)
            else:
                print('>>                           No valid observation, skip ...')
        #print('>>  ... global mean is {}'.format(np.nanmean(Xa)))
        Xa_output_all[:,:,reconi] = np.copy(Xa) # for each reconi, all observations were assimilated. Save final result for this reconi
        
    if Xb is not None:
        lenn1 = f.get('Xb').shape[0]
        Xa_output_2d = Xa_output_all[0:lenn1,:,:]
        if Xb3d:
            lenn2 = f.get('Xb3d').shape[0]
            Xa_output_3d = Xa_output_all[lenn1:lenn2+lenn1,:,:]
    elif Xb is None:
        if Xb3d:
            lenn2 = f.get('Xb3d').shape[0]
            Xa_output_3d = Xa_output_all[0:lenn2,:,:]
    else:
        print('>>  Error! No 3d or 2d variables are listed in {}'.format(config_name))
print('>>  All Done')



# DA save output
with h5py.File(hdf5name, 'r') as f:

    print('')
    print('>>  Start writing netCDF ...')
    
    # save netCDF file
    nf = Dataset(nc_filename, 'w', format='NETCDF4')
    nf.description = 'DeepDA' + nc_filename
    #Specifying dimensions
    nf.createDimension('lon', len(cGENIEGridB_lon36))
    nf.createDimension('lat', len(cGENIEGridB_lat36))
    z = np.arange(0,1,1) # level 2d
    nf.createDimension('z', len(z))  # level
    nf.createDimension('nens', nens)  # number of ens
    nf.createDimension('time', recon_period_len)
    # Building variables
    longitude = nf.createVariable('Longitude', 'f4', 'lon')
    # Passing data into variables
    longitude[:] = cGENIEGridB_lon36.values

    latitude = nf.createVariable('Latitude', 'f4', 'lat')
    latitude[:] = cGENIEGridB_lat36.values

    levels = nf.createVariable('Levels', 'i4', 'z')
    levels[:] = z  # 2d level
    if Xb3d is not None:
        nf.createDimension('zt', len(zt))
        levels = nf.createVariable('zt', 'f4', 'zt')
        levels[:] = zt
        
    if locRad:
        #nf.createDimension('prior_var', prior_variable_len)  # level
        covloc_nc = nf.createVariable('covloc', 'f4', ('lat', 'lon'))
        covloc_nc[:,:] = np.copy(covloc[0:dum_ijmax].reshape(dum_jmax,dum_imax))
        
    if Xb is not None:
        for nc_var_i in range(prior_variable_len):
            nc_var_name = prior_variable_dict[nc_var_i]

            j0 = dum_ijmax * nc_var_i
            j1 = dum_ijmax * (nc_var_i+1)
            print('>>    id from {} to {}: {}'.format(j0, j1,nc_var_name))

            Xb0_i = np.copy(f.get('Xb')[j0:j1,:])
            
            Xa_output_i = np.copy(Xa_output_2d[j0:j1,:,:])
            Xa_outputi = Xa_output_i.reshape(dum_imax,dum_jmax,nens,recon_period_len)

            XbNC_mean = nf.createVariable(nc_var_name+'_Xb_mean', 'f4', ('lat', 'lon','z'))
            xbm = np.mean(Xb0_i,axis=1)
            XbNC_mean[:,:,:] = np.copy(xbm.reshape(dum_jmax,dum_imax,1))

            XbNC_variance = nf.createVariable(nc_var_name+'_Xb_variance', 'f4', ('lat', 'lon','z'))
            Xb_temp = np.copy(np.var(Xb0_i,axis=1).reshape(dum_jmax,dum_imax,1))
            Xb_temp = np.ma.MaskedArray(Xb_temp, np.copy(xbm.reshape(dum_jmax,dum_imax,1)) >= 9.9692e+36)
            XbNC_variance[:,:,:] = Xb_temp
            print('>>    Xb mean is {}, std is {}, var is {}'.format(np.nanmean(XbNC_mean), np.sqrt(np.nanmean(Xb_temp)), np.nanmean(Xb_temp)))

            XaNC_mean = nf.createVariable(nc_var_name+'_Xa_mean', 'f4', ('lat', 'lon','z','time'))
            Xam_temp = np.copy(np.nanmean(Xa_outputi,axis=2).reshape(dum_jmax,dum_imax,1,recon_period_len))
            XaNC_mean[:,:,:,:] = Xam_temp

            XaNC_variance = nf.createVariable(nc_var_name+'_Xa_variance', 'f4', ('lat', 'lon','z','time'))
            #print(Xa_outputi[0,0:36,0,0])
            Xa_temp = np.copy(np.ma.var(Xa_outputi,axis=2).reshape(dum_jmax,dum_imax,1,recon_period_len))
            Xa_temp = np.ma.MaskedArray(Xa_temp, Xam_temp >= 9.9692e+36)
            XaNC_variance[:,:,:,:] = Xa_temp
            
            for reconii in range(recon_period_len):
                XaNC_mean_i = XaNC_mean[:,:,:,reconii]
                XaNC_var_i = XaNC_variance[:,:,:,reconii]
                print('>>    recon {}. Xa mean is {}, std is {}, var is {}'.format(reconii, np.nanmean(XaNC_mean_i),np.sqrt(np.nanmean(XaNC_var_i)), np.nanmean(XaNC_var_i)))
            
            if save_ens_full:
                XaNC_full = nf.createVariable(nc_var_name+'_Xa_full', 'f4', ('lat', 'lon', 'nens', 'z','time'))
                XaNC_full[:,:,:,:,:] = np.copy(Xa_outputi.reshape(dum_jmax,dum_imax,nens,1,recon_period_len))
                
                XbNC_full = nf.createVariable(nc_var_name+'_Xb_full', 'f4', ('lat', 'lon', 'nens', 'z'))
                XbNC_full[:,:,:,:] = np.copy(Xb0_i.reshape(dum_jmax,dum_imax,nens,1))

            if kcov_saving > 0:
                kcov_i = np.copy(kcov[j0:j1]).reshape(dum_imax,dum_jmax,1)
                kcov_i = np.ma.MaskedArray(kcov_i, np.copy(xbm.reshape(dum_jmax,dum_imax,1)) >= 9.9692e+36)
                cov_ob0 = nf.createVariable(nc_var_name+'_obs0'+'_cov', 'f4', ('lat', 'lon','z'))
                cov_ob0[:,:,:] = kcov_i

            #Add local attributes to variable instances
            longitude.units = '°'
            latitude.units = '°'
            levels.units = 'm'
            XbNC_mean.units = '°C'
            XbNC_variance.units = '°C^2'
            if save_ens_full:
                XaNC_full.units = '°C'
                XbNC_full.units = '°C'

            #variance.warning = 'test ...'
    if Xb3d is not None:
        for nc_var_i in range(prior_variable3d_len):
            nc_var_name = prior_variable_dict_3d[nc_var_i]

            j0 = dum_ijmax * dum_dmax * nc_var_i
            j1 = dum_ijmax * dum_dmax * (nc_var_i+1)
            print('>>  Writing 3d field. ID from {} to {}: {}'.format(j0, j1,nc_var_name))

            Xb0_i = np.copy(f.get('Xb3d')[j0:j1,:])
            Xa_output_i = np.copy(Xa_output_3d[j0:j1,:,:])
            Xa_outputi = Xa_output_i.reshape(dum_imax, dum_jmax,dum_dmax, nens,recon_period_len)
            
            XbNC_mean = nf.createVariable(nc_var_name+'_Xb_3d_mean', 'f4', ( 'zt', 'lat','lon'))
            xbm = np.mean(Xb0_i,axis=1)
            XbNC_mean[:,:,:] = np.copy(xbm.reshape(dum_dmax,dum_jmax,dum_imax))

            XbNC_variance = nf.createVariable(nc_var_name+'_Xb_3d_variance', 'f4', ( 'zt', 'lat','lon'))
            Xb_temp = np.copy(np.var(Xb0_i,axis= 1).reshape(dum_dmax,dum_jmax,dum_imax))
            Xb_temp = np.ma.MaskedArray(Xb_temp, np.copy(xbm.reshape(dum_dmax,dum_jmax,dum_imax)) >= 9.9692e+36)
            XbNC_variance[:,:,:] = Xb_temp
            print('>>    Xb mean is {}, std is {}, var is {}'.format(np.nanmean(XbNC_mean),np.sqrt(np.nanmean(Xb_temp)), np.nanmean(Xb_temp)))

            XaNC_mean = nf.createVariable(nc_var_name+'_Xa_3d_mean', 'f4', ('zt','lat', 'lon','time'))
            Xam_temp = np.copy(np.nanmean(Xa_outputi,axis=3).reshape(dum_dmax,dum_jmax,dum_imax,recon_period_len))
            XaNC_mean[:,:,:,:] = Xam_temp

            XaNC_variance = nf.createVariable(nc_var_name+'_Xa_3d_variance', 'f4', ('zt','lat', 'lon','time'))
            Xa_temp = np.copy(np.ma.var(Xa_outputi,axis=3).reshape(dum_dmax,dum_jmax,dum_imax,recon_period_len))
            Xa_temp = np.ma.MaskedArray(Xa_temp, Xam_temp >= 9.9692e+36)
            XaNC_variance[:,:,:,:] = Xa_temp
            
            for reconii in range(recon_period_len):
                XaNC_mean_i = XaNC_mean[:,:,:,reconii]
                XaNC_var_i = XaNC_variance[:,:,:,reconii]
                print('>>    recon {}. Xa mean is {}, std is {}, var is {}'.format(reconii, np.nanmean(XaNC_mean_i), np.sqrt(np.nanmean(XaNC_var_i)), np.nanmean(XaNC_var_i)))

            if save_ens_full:
                XaNC_full = nf.createVariable(nc_var_name+'_Xa_3d_full', 'f4', ('zt','lat', 'lon', 'nens', 'time'))
                XaNC_full[:,:,:,:,:] = np.copy(Xa_outputi.reshape(dum_dmax,dum_jmax,dum_imax,nens,recon_period_len))
                
                XbNC_full = nf.createVariable(nc_var_name+'_Xb_3d_full', 'f4', ('zt','lat', 'lon', 'nens'))
                XbNC_full[:,:,:,:] = np.copy(Xb0_i.reshape(dum_dmax,dum_jmax,dum_imax,nens))

            if kcov_saving > 0:
                kcov_i = np.copy(kcov[lenn1:lenn1+dum_ijmax*dum_dmax]).reshape(dum_dmax,dum_jmax,dum_imax)
                kcov_i = np.ma.MaskedArray(kcov_i, np.copy(xbm.reshape(dum_dmax,dum_jmax,dum_imax)) >= 9.9692e+36)
                cov_ob0 = nf.createVariable(nc_var_name+'_3d_obs0'+'_cov', 'f4', ( 'zt', 'lat','lon'))
                cov_ob0[:,:,:] = kcov_i

            #Add local attributes to variable instances
            longitude.units = '°'
            latitude.units = '°'
            levels.units = 'm'
            XbNC_mean.units = '°C'
            XbNC_variance.units = '°C^2'
            if save_ens_full:
                XaNC_full.units = '°C'
                XbNC_full.units = '°C'
    # Closing the dataset
    nf.close()  # close the new file
    print('>>  End writing netCDF')
    config_save_name = dir_proxy_data + nexp + '.yml'
    configos = 'cp ' + config_name + ' ' +  config_save_name
    os.system(configos)
print('')    
print(config_save_name)
print(nc_filename)
print('')  
print('************  All saved  ************')