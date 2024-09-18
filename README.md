# DeepDA

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1234567.svg)](https://doi.org/10.5281/zenodo.1234567)


---
`DeepDA is a data assimilation framework for paleoclimate projects focusing on deep time.`


## Warning

**This project is under active development.**

---

## Contact:

>  **Mingsong Li** 
>  
>  Peking University
>  
>  Email: msli {at} pku.edu.cn
> 
>  Sept 18, 2024
 
---

## Data Assimilation Framework Structure

### DeepDA_config.yml
> Configuration parameters for running a data assimilation reconstruction using cGENIE priors.
> 
> Defines **almost** all parameters.

### DeepDA_allMC.ipynb   
> Jupyter notebook. 
> 
> Main function. 
> 
> Users must adjust the code for their own projects.

### DeepDA_lib/
* **DeepDA_psm.py**
	
> Functions useful for proxy system models

* **DeepDA_tools.py**
	
> Miscellaneous functions
	
* **LMR_DA.py**

> DA core code borrowed from the LMR project by G. J. Hakim, with additional code by L. Madaus
	
* **modules_nc.py**   

> Functions useful for working with the prior data
    
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

	*Grid mask definition for the Atlantic in cGENIE*
	
#### proxy/
* **petmproxy3slices_v0.1.csv**

	*proxy database*
	
* **PETMTEX_baysparSettings_v0.2.csv**

	*Tolerance settings for BAYSPAR*

#### wrk/

* A directory for saving outputs*
    
### utils/

*utilities for analyzing the outputs*

* CESMreadNC4R18O.ipynb    

> 	Reads d18Osw from CESM simulation netCDF files

* correct_cgenie_biogem_2d.ipynb

> 	Generates surface fields for cGENIE biogem_2d.nc file

* correct_cgenie_carb_ohm_cal.ipynb

> 	Selects bottom water calcite saturation fields ('carb_ohm_cal_ben') using user-defined bathymetry
> 
>   Updates cGENIE output files for sedgem (fields_sedgem_2d.nc)

* correct_cgenie_sed_caco3_13c.ipynb

> 	Corrects the cGENIE sed_CaCO3_13C field
	
* ReadcGENIEresViaName.ipynb

> 	Reads cGENIE .res files by name
	
* ReadDASummaryXlsxPlot-prePETM_SST.ipynb

> 	Reads the wrk folder and writes *.summary.xlsx; extracts warming and cooling data for all variables

* ReadDASummaryXlsxPlot-SST.ipynb

> 	Same as above, but for different SST-related data


* ReadDASummaryXlsxPlot.ipynb

> 	Reads the wrk folder and writes *.summary.xlsx; extracts warming and cooling data for all variables


    
### verification/

* DeepDA_verify_proxyunit-jobs.ipynb

> 	Verifies DA outputs by reading proxy, prior, and posterior data from DA output files and config files.
> 	
> 	Calculates statistics (corrcoef and CE) and saves outputs.
	
* DeepDA_verify_proxyunit-testR.ipynb

> 	Similar functionality as above for verifying DA outputs.  
                                          
---
        
## Prior

### Directory Structure

	../../ML.petm/
				ML.petm029.ID.1/    # folder, exp 1
				ML.petm029.ID.2/    # folder, exp 2
				ML.petm029.ID.3/    # folder, exp 3
				...
				ML.petm029.ID.n/    # folder, exp n
                 
Each experiment folder contains:

	ML.petm029.ID.1/
				archive/
				atchem/
				biogem/
				ecogem/
				...
				rokgem/
				sedgem/
					
Refer to the [cGENIE.muffin Earth system model](https://www.seao2.info/mycgenie.html) for more details on the model structure.

	
---
---

# How to Run DeepDA?

### 1. Create an environment:

	conda create --name deepda
        
 Activate the environment:

	conda activate deepda

 To deactivate, use:

	conda deactivate

### 2. Clone the DeepDA repository
#### git clone

    git clone https://github.com/mingsongli/deepDA.git

#### Navigate to the misc directory:

    cd misc

### 3. Setup environment

#### Install packages using the `deepda_pyenv.yml` file:

    conda env update --file deepda_pyenv.yml

#### 4. Install PSMs: BAYSPAR, BAYFOX & BAYMAG

##### Clone the repositories:

    cd ../..
    git clone https://github.com/mingsongli/bayfox.git
    git clone https://github.com/mingsongli/baysparpy.git
    git clone https://github.com/mingsongli/BAYMAG.git

> Note: Some adjustments were made for the DA project, and the original code may not work as expected.

##### To install:

###### For BAYSPAR:

    cd baysparpy
    python setup.py install

###### For BAYFOX:

    cd ../bayfox
    python setup.py install
    
###### For BAYMAG
    
    cd ../baymagpy
    python setup.py install
    
 Read more about [python setup.py](https://stackoverflow.com/questions/19048732/python-setup-py-develop-vs-install)

### 5. Data Assimilation Configuration

Navigate to the `deepDA` folder:

    cd ../deepDA

Run `jupyterlab`:

    jupyter lab

Revise settings in **`DeepDA_config.yml`** accordingly.

Open **`DeepDA_allMC.ipynb`** to start the notebook.

*Adjust the code as needed for your project.*

### 6. Run the Reconstruction

Click the triangle button to execute the notebook.

The output will be saved in the user-defined directory.

Good luck!

---

## Reminder

### If the above steps do not work:

#### 1. Ensure that all required Python packages are installed.

 See
 
 	/misc/deepda_pyenv.yml
 	
 	/misc/lmr_py3EnvCondaList.docx

####  2. Revise settings in `DeepDA_config.yml` to suit your project.

####  3. If the prior lacks necessary fields, rerun the simulation to update the prior.

Use the

	correct_xxx.ipynb
	
scripts in the `utils/` folder for minor corrections.


## Useful Resource

[LMR project Documentation](https://atmos.washington.edu/~hakim/lmr/docs/index.html)
