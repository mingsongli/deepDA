# DeepDA

---
`DeepDA is a data assimilation framework for deep time paleoclimate project.`


Contact: Mingsong Li
         Peking University

Email: msli{at}pku.edu.cn

## Highlights
> 
## Structure

- DeepDA_config.yml    # Configuration parameters for running a data assimilation reconstruction using cGENIE priors.
                       # It defines almost all parameters.
- DeepDA_allMC.ipynb   # Jupyter notebook. main function. Usually users have to adjust the code for their own project
                       # 
- DeepDA_lib
    - DeepDA_psm.py    # functions useful for the proxy system models
    - DeepDA_tools.py  # miscellaneous functions
    - LMR_DA.py        # DA core code; borrowed from LMR project by G. J. Hakim with code borrowed from L. Madaus
    - modules_nc.py    # functions useful for the prior
    
- misc
    - deepda_pyenv.yml    # python environment dependencies
    
- mlwrk
    - data_misc
        - cGENIEGrid.csv                 # grid definition for cGENIE
        - cGENIEGridBound.csv            # grid definition for cGENIE
        - oceanmask_p0055c_Atlantic.csv  # grid definition for cGENIE
        - ...
    - proxy
        - petmproxy3slices_v0.1.csv         # proxy database
        - PETMTEX_baysparSettings_v0.2.csv  # tolerance settings for BAYSPAR
        - ...
    - wrk
        -       # directory for saving outputs
- src
    - bayfox      # bayfox python package,  unused?
    - baymagpy    # baymag python package,  unused?
    
- utils
    - CESMreadNC4R18O.ipynb    # read d18Osw from netCDF file of CESM simulations
    - correct_cgenie_biogem_2d.ipynb   # Generate surface field for cGENIE biogem_2d.nc file
    - correct_cgenie_carb_ohm_cal.ipynb # Pick bottom water calcite saturation field 'carb_ohm_cal_ben' using user-defined bathymetry
                                        # Update cGENIE output file of sedgem fields_sedgem_2d.nc
    - correct_cgenie_sed_caco3_13c.ipynb # Correct cGENIE field of sed_CaCO3_13C
    - ReadcGENIEresViaName.ipynb  # Read cGENIE .res files via names
    - ReadDASummaryXlsxPlot-prePETM_SST.ipynb  # Read wrk folder for writing *.summary.xlsx; extract warming and cooling data for all variables
    - ReadDASummaryXlsxPlot-SST.ipynb          # Read wrk folder for writing *.summary.xlsx; extract warming and cooling data for all variables 
    - ReadDASummaryXlsxPlot.ipynb              # Read wrk folder for writing *.summary.xlsx; extract warming and cooling data for all variables 
    
- verification
    - DeepDA_verify_proxyunit-jobs.ipynb  # DeepDA_verify is able to verify DA outputs
                                          # It reads proxy, prior, and posterior from DA output files and configuration files.
                                          # Then, it calculates statistics (corrcoef and CE) of the DA results and saves outputs
    - DeepDA_verify_proxyunit-testR.ipynb # DeepDA_verify is able to verify DA outputs
                                          # It reads proxy, prior, and posterior from DA output files and configuration files.
                                          # Then, it calculates statistics (corrcoef and CE) of the DA results and saves outputs
                                          
                                          
                                          

## How to run

###  0. Ensure that all python packages have been installed.
        See /misc/deepda_pyenv.yml  & /misc/lmr_py3EnvCondaList.docx
###  1. Revise settings in the DeepDA_config.yml accordingly
        May need to adjust the code for their own project
###  2. If the prior lacking required fileds, rerun the simulation to update the prior. 
        The code entilted correc_xxx.ipynb within the 'utils' folder may be helpful for some minor corrections.
###  3. Run DeepDA_allMC.ipynb
        The DA output will be saved at the user-defined directory
        
## Data



#### Useful sources:

    https://atmos.washington.edu/~hakim/lmr/docs/index.html
    
    