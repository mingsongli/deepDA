'''
    Modules for prior
    
    cal_find_ij
        : find loc_i, and loc_j, return dum_ij
        : need dum_lon, dum_lat, dum_lon_offset, dum_imax, dum_jmax
        : corrected by Mingsong Li, Penn State, Dec. 5, 2018

    get_gdep
        : return
        : need kmax, ez0, maxdep

    calc_find_k
        : return k
        : need dum_depth, dum_kmax
        : calls for get_gdep

    cal_find_t
        : return t
        : need dum_time

    is_number
        : return True or False
        
    haversine
    
    cov_local_i
    
    covloc_eval

EXAMPLE:

    find_ij = cal_find_ij(dum_lon,dum_lat,dum_lon_offset,dum_imax,dum_jmax)
    find_k  = calc_find_k(dum_depth,dum_kmax)

    By Mingsong Li, Penn State, 2017

'''
import numpy as np
from math import sin, pi, pow

# *********************************************************************** %
#  A function 
def cal_find_ij(dum_lon,dum_lat,dum_lon_offset,dum_imax,dum_jmax):
# *** INITIALIZE PARAMETERS & VARIABLES ********************************* %
# dum_lon: log in degree
# dum_lat: lat in degree
# dum_lon_offset: -180 for cGENIE
# dum_imax: grid length, 36 for cGENIE
# dum_jmax: grid length, 36 for cGENIE
# output [lon, lat], 
#      lon ranges from 0 (-180) to 35 (180), lat ranges from 0 (-90) to 35 (90)
# *********************************************************************** %
# set passed parameters
    lon_offset = dum_lon_offset
    imax = dum_imax
    jmax = dum_jmax
# LOCAL VARIABLES
    loc_dlon = 360.0/imax

# *** CALCULATE (i,j) *************************************************** %
# CALCULATE 'i'
# precondition lon
    if dum_lon >= (360.0 + lon_offset):
        dum_lon = dum_lon - 360.0

    if dum_lon < lon_offset:
        dum_lon = dum_lon + 360.0

# calculate 'i'
    #loc_i = int((dum_lon - lon_offset)/loc_dlon + 0.5)
    loc_i = int((dum_lon - lon_offset)/loc_dlon) # ML corrected for cGENIE
    if loc_i > imax:
         loc_i = imax

# CALCULATE 'j'
# check lat
    if (dum_lat > 90) or (dum_lat < -90):
        disp(['warning. lat out-of-range ...'])

# calculate 'j'
    loc_sinlat = sin(pi*dum_lat/180.0)
    #loc_j = int(jmax*0.5*(1.0 + loc_sinlat) + 0.5)
    loc_j = int(jmax*0.5*(1.0 + loc_sinlat)) # ML corrected for cGENIE
    if (loc_j > jmax):
        loc_j = jmax
    dum_ij = [loc_i, loc_j]
#              lon,   lat
    return dum_ij

# *********************************************************************** %
def get_gdep(kmax, ez0, maxdep):
    z1 = ez0*(pow((1.0 + 1/ez0),(1.0/kmax)) - 1.0)
    tv4 = ez0*(pow((z1/ez0+1),0.5)-1)
    tv2 = 0
    tv1 = 0
    zro = np.zeros((kmax,1))
    zw = np.zeros((kmax+1,1))
    dz = np.zeros((kmax,1))
    dza = np.zeros((kmax,1))
    zro[kmax-1] = -tv4
    zw[kmax-1] = tv2

    for k in np.arange(1,kmax):
        if ez0 > 0:
            tv3 = ez0*(pow((z1/ez0+1),k)-1)
            dz[kmax-k] = tv3 - tv2
            tv2 = tv3
            tv5 = ez0*(pow((z1/ez0+1),(k+0.5))-1)
            if k < kmax:
                dza[kmax-k-1] = tv5 - tv4
            tv4 = tv5
            tv1 = tv1 + dz[kmax-k]
        else:
            dz[k-1] = 1.0/kmax
            dza[k-1] = 1.0/kmax

    for k in np.arange(kmax,1,-1):
        if k > 1:
            zro[k-2] = zro[k-1] - dza[k-2]
        zw[k-1] = zw[k] - dz[k-1]
    a = zw * maxdep
    a[0] = -1*maxdep # ensure the deepest one is 5km
    return a

# *********************************************************************** %
def calc_find_k(dum_depth, dum_kmax):
    par_D_max = 5000.0          # max depth (m)
    par_ez0 = 0.1               # GOLDSTEIN ez0 parameter
# set depth array
    D = get_gdep(dum_kmax,par_ez0,par_D_max) * -1
    #print(D)
    # set max k as default value - surface ocean
    loc_n_k = dum_kmax
    for n_k in np.arange(dum_kmax,2,-1):
        if dum_depth > D[n_k-1]:
            loc_n_k = n_k - 1
    if dum_depth > D[0]:
        loc_n_k = 1
    dum_k = [loc_n_k]
    return dum_k

# *********************************************************************** %
from netCDF4 import Dataset
#from operator import itemgetter
def cal_find_t(name_nc,dum_t):
        # open nc file as dataset
    dataset = Dataset(name_nc)
    time = dataset.variables['time']
#    val = 0
    mylist = [i-dum_t for i in time]
    val, idx = min((val, idx) for (idx, val) in enumerate(my_list))
    if val == 0:
        return idx
    else:
        return idx
        print('No such value, the nearest value is {}, index is {}'.format(val,idx))

# *********************************************************************** %
def is_number(n):
    try:
        float(n)
    except ValueError:
        return False
    return True
# *********************************************************************** %


def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """

    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = list(map(np.radians, [lon1, lat1, lon2, lat2]))
    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = np.sin(dlat/2.)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2.)**2
    c = 2 * np.arcsin(np.sqrt(a))
    km = 6367.0 * c
    return km


def cov_local_i(locRad, lon1, lat1, lon2, lat2):
    '''
    locRad: covariance localization distance, in km
    lon1:
    lat1:
    lon2:
    lat2:
    '''
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = list(map(np.radians, [lon1, lat1, lon2, lat2]))
    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = np.sin(dlat/2.)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2.)**2
    c = 2 * np.arcsin(np.sqrt(a))
    km = 6367.0 * c
    #
    hlr = 0.5*locRad; # work with half the localization radius
    r = km/hlr;
    
    # Gaspari-Cohn function
    # for pts within 1/2 of localization radius
    if km <= hlr:
        covLoc= (((-0.25*r+0.5)*r+0.625)* \
                 r-(5.0/3.0))*(r**2)+1.0
    # for pts between 1/2 and one localization radius
    if km >  hlr:
        covLoc = ((((r/12. - 0.5) * r + 0.625) * \
                          r + 5.0/3.0) * r - 5.0) * \
                          r + 4.0 - 2.0/(3.0*r)
    # Impose zero for pts outside of localization radius
    if km >  2.*hlr:
        covLoc = 0.0
    return covLoc

def covloc_eval(locRad, yo_loc,dum_jmax,dum_imax,cGENIEGrid):
    dum_ijmax = dum_jmax * dum_imax
    lon1 = yo_loc[0] #  longitude
    lat1 = yo_loc[1] # latitude
    LAT, LON = np.meshgrid(cGENIEGrid[:,0], cGENIEGrid[:,1])
    #LON, LAT = np.meshgrid(cGENIEGrid['lon'], cGENIEGrid['lat'])
    
    LON = LON.T.reshape(dum_ijmax,1)
    LAT = LAT.T.reshape(dum_ijmax,1)
    covloc = np.full((dum_ijmax,1),np.nan)
    
    for ii in range(dum_ijmax):
        #if ii % 500 == 0:
        #    print('  covariance localization estimation: {} of {}'.format(ii, dum_ijmax))
        lon2 = LON[ii]
        lat2 = LAT[ii]
        covloc[ii] = cov_local_i(locRad, lon1, lat1, lon2, lat2) # lon-lat
        #if ii < 50:
        #    print('      lon2 {}, lat2 {}, covloc {}'.format(lon2,lat2,covloc[ii]))
    return covloc