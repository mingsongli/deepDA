# DeepDA

---
`DeepDA is a data assimilation framework for deep time paleoclimate project.`


## Warning

**This project is under heavy development.**

---

## Contact:

>  **Mingsong Li** 
>  
>  Peking University
>  
>  Email: msli {at} pku.edu.cn
> 
>  April 7, 2022
 
---

## DA Frame Structure

### DeepDA_config.yml
> Configuration parameters for running a data assimilation reconstruction using cGENIE priors.
> 
> It defines **almost** all parameters.

### DeepDA_allMC.ipynb   
> Jupyter notebook. 
> 
> Main function. 
> 
> Users have to adjust the code for their own projects.

### DeepDA_lib/
* **DeepDA_psm.py**
	
> Functions useful for the proxy system models

* **DeepDA_tools.py**
	
> Miscellaneous functions
	
* **LMR_DA.py**

> DA core code; borrowed from LMR project by G. J. Hakim with code borrowed from L. Madaus
	
* **modules_nc.py**   

> Functions useful for the prior
    
### misc/
* **deepda_pyenv.yml**

	*python environment dependencies*
    
### mlwrk/

#### data_misc/
* **cGENIEGrid.csv**

	*grid definition for cGENIE*
        
* **cGENIEGridBound.csv**
	
	*grid definition for cGENIE*
	
* **oceanmask_p0055c_Atlantic.csv**

	*grid definition for cGENIE*
	
#### proxy/
* **petmproxy3slices_v0.1.csv**

	*proxy database*
	
* **PETMTEX_baysparSettings_v0.2.csv**

	*tolerance settings for BAYSPAR*

#### wrk/

*directory for saving outputs*
    
### utils/

*utilities for the analysis of the outputs*

* CESMreadNC4R18O.ipynb    

> 	Read d18Osw from netCDF file of CESM simulations

* correct_cgenie_biogem_2d.ipynb

> 	Generate surface field for cGENIE `biogem_2d.nc` file

* correct_cgenie_carb_ohm_cal.ipynb

> 	Pick bottom water calcite saturation field '`carb_ohm_cal_ben`' using user-defined bathymetry.
> 
>   Update cGENIE output file of sedgem `fields_sedgem_2d.nc`

* correct_cgenie_sed_caco3_13c.ipynb

> 	Correct cGENIE field of `sed_CaCO3_13C`
	
* ReadcGENIEresViaName.ipynb

> 	Read cGENIE `.res` files via names
	
* ReadDASummaryXlsxPlot-prePETM_SST.ipynb

> 	Read `wrk`/ folder for writing `*.summary.xlsx`; extract warming and cooling data for all variables

* ReadDASummaryXlsxPlot-SST.ipynb

> 	Read `wrk`/ folder for writing `*.summary.xlsx`; extract warming and cooling data for all variables 


* ReadDASummaryXlsxPlot.ipynb

> 	Read `wrk`/ folder for writing `*.summary.xlsx`; extract warming and cooling data for all variables 


    
### verification/

* DeepDA_verify_proxyunit-jobs.ipynb

> 	DeepDA_verify is able to verify DA outputs.
> 	
> 	It reads proxy, prior, and posterior from DA output files and configuration files.
> 	
> 	Then, it calculates statistics (corrcoef and CE) of the DA results and saves outputs.
	
* DeepDA_verify_proxyunit-testR.ipynb

> 	DeepDA_verify is able to verify DA outputs.
> 	
> 	It reads proxy, prior, and posterior from DA output files and configuration files.
> 	
> 	Then, it calculates statistics (corrcoef and CE) of the DA results and saves outputs.
                                          
                                          
                                          
---
        
## Prior

### Directory

	../../ML.petm/
				ML.petm029.ID.1/    # folder, exp 1
				ML.petm029.ID.2/    # folder, exp 2
				ML.petm029.ID.3/    # folder, exp 3
				...
				ML.petm029.ID.n/    # folder, exp n
                 
Each folder of experiment:

	ML.petm029.ID.1/
				archive/
				atchem/
				biogem/
				ecogem/
				...
				rokgem/
				sedgem/
					
See [cGENIE.muffin Earth system model](https://www.seao2.info/mycgenie.html) for the structure of cGENIE model.

	
---
---

# How to run deepDA?

### 1. Create an environment:

	conda create --name deepda
        
 Activate this environment, use

	conda activate deepda

 To deactivate an active environment, use

	conda deactivate

### 2. Clone the deepDA frame code
#### git clone

    git clone https://github.com/mingsongli/deepDA.git

#### Change directory

    cd misc

### 3. Setup environement

#### Install packages using `deepda_pyenv.yml` environment file:

    conda env update --file deepda_pyenv.yml

#### Install PSMs: BAYSPAR, BAYFOX & BAYMAG

##### Clone

    cd ../..
    git clone https://github.com/mingsongli/bayfox.git
    git clone https://github.com/mingsongli/baysparpy.git
    git clone https://github.com/mingsongli/BAYMAG.git

> I did some corrections/adjustments for the DA project. The original code may not work as expected.

##### Install

Go to each package, install the package from the local directory

###### BAYSPAR

    cd baysparpy
    python setup.py install

###### BAYFOX

    cd ../bayfox
    python setup.py install
    
###### BAYMAG
    
    cd ../baymagpy
    python setup.py install
    
 Read more [python setup.py](https://stackoverflow.com/questions/19048732/python-setup-py-develop-vs-install)

### 4. DA configuration

Go to the `deepDA` folder

    cd ../deepDA

Run `jupyterlab`

    jupyter lab

Revise settings in the **`DeepDA_config.yml`** accordingly

Double click **`DeepDA_config.yml`** to revise the settings.

Double click **`DeepDA_allMC.ipynb`** to open the notebook.

*May need to adjust the code for their own project*

### 5. Reconstruction

Click the triangle to run the notebook.

The DA output will be saved in the user-defined directory.

Good luck!

---

## Reminder

### If the above steps do not work, you may want to:

#### 1. ensure that all required python packages have been installed.

 See
 
 	/misc/deepda_pyenv.yml
 	
 	/misc/lmr_py3EnvCondaList.docx

####  2. Revise settings in the `DeepDA_config.yml` accordingly

May need to adjust the code for your own project.

####  3. If the prior lacks required fields, rerun the simulation to update the prior. 

The code entitled

	correct_xxx.ipynb
	
within the `utils/` folder may be helpful for some minor corrections.


## Useful resource

[LMR project](https://atmos.washington.edu/~hakim/lmr/docs/index.html)
