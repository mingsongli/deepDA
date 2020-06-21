#==========================================================================================
# Data assimilation for deep time tools. 
#
#==========================================================================================

import numpy as np

def deepda_hard_limit(Xa, yml_dict, prior_variable_dict, dum_ijmax,verbose):
    """
    Function to ensure the Xa value within the hard limit

    Originator: Mingsong Li, Penn State

    Revisions:

    June 17 2020: 
                  
    
    -----------------------------------------------------------------
     Inputs:
          Xa: ensemble estimates of state (Nx x Nens) 
          yml_dict: loaded yaml file
          prior_variable_dict: dict of variable for prior
    """
    prior_source = yml_dict['prior']['prior_source'] #
    limit_hard_keys = list(yml_dict['prior'][prior_source]['limit_hard'].keys())
    #print(limit_hard_keys)
    
    # 2d variable
    prior_variable_len = len(prior_variable_dict)
    
    for Xa2d_vari in range(prior_variable_len):
        if prior_variable_dict[Xa2d_vari] in limit_hard_keys:
            
            # Read hard limit. Some variables have hard limitation: e.g., CaCO3 = [0, 100]
            lim_min = yml_dict['prior'][prior_source]['limit_hard'][prior_variable_dict[Xa2d_vari]]['lim_min']
            lim_max = yml_dict['prior'][prior_source]['limit_hard'][prior_variable_dict[Xa2d_vari]]['lim_max']
            #print(' {} min {} max {}'.format(prior_variable_dict[Xa2d_vari], lim_min, lim_max))
            # read Xa index
            # for 2d variables only!! Need to include 3D variables
            j0 = dum_ijmax * Xa2d_vari
            j1 = dum_ijmax * (Xa2d_vari+1)
            
            if lim_min is not None or lim_max is not None:
                # read selected variable
                Xa_full_vari = np.copy(Xa[j0:j1,:])
                
                if lim_min is not None:
                    #Xa_full_vari[np.logical_and(Xa_full_vari<lim_min,Xa_full_vari<9.9692e+36)] = lim_min
                    if np.any(Xa_full_vari < lim_min):
                        Xa_full_vari[Xa_full_vari<lim_min] = lim_min
                        if verbose:
                            print(' set variable {} : id {} - {} min limit to {}'.format(prior_variable_dict[Xa2d_vari],j0,j1,lim_min))
                if lim_max is not None:
                    if np.any(np.logical_and(Xa_full_vari>lim_max,Xa_full_vari<9.9692e+36)):
                        Xa_full_vari[np.logical_and(Xa_full_vari>lim_max,Xa_full_vari<9.9692e+36)] = lim_max
                        #if np.any(Xa_full_vari>lim_max):
                        #    Xa_full_vari[Xa_full_vari>lim_max] = lim_max
                        if verbose:
                            print(' set variable {} : id {} - {} max limit to {}'.format(prior_variable_dict[Xa2d_vari],j0,j1,lim_max))
                # save back
                Xa[j0:j1,:] = Xa_full_vari
                #Xa = np.ma.MaskedArray(Xa, Xa >= 9.9692e+36)
    # Return the full state
    return Xa

