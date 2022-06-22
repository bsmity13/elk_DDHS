# Code and data for, "Density-dependent habitat selection alters drivers of population distribution in northern Yellowstone elk."

_Authors_:  

  - Brian J. Smith <a itemprop="sameAs" content="https://orcid.org/0000-0002-0531-0492" href="https://orcid.org/0000-0002-0531-0492" target="orcid.widget" rel="me noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" alt="ORCID iD icon" style="width:1em;margin-right:.5em;"/></a>

  - Dan MacNulty
  - Daniel R. Stahler
  - Douglas W. Smith
  - Tal Avgar

## Manuscript status
To be published to bioRxiv week of June 22, 2022.  

To be submitted for review to *Ecology Letters*.

## About Repository

Repository version 0.1 -- Prior to peer review

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6687905.svg)](https://doi.org/10.5281/zenodo.6687905)

This repository contains the processed data and code needed to fit the model described in the manuscript. The models are fit using MCMC in R with the package `nimble`. The posterior samples are stored in directory `models/`, but they are almost 3 GB in total. 

To save space in this repository, the posterior samples are not included; however, a very small sample of the posterior is included in the directory `demo_models` to help demo the code. Scripts, such as `02_figures.R`, contain a variable `demo_data` that can be set to `TRUE` to load the sample posterior for testing the code.

Users can reproduce the full posterior by running the scripts in this repo, but be aware that it takes approximately 10 days to fit the full model. In this case, `demo_data` should be set to `FALSE` to use the real posterior samples, *e.g.*, to create figures.

### Data

All data are contained in the directory `data/`. Data have been pre-processed from raw counts to a format used to fit the models.

- `all_data_test.csv` -- comma-delimited file; testing dataset used for model validation
- `all_data.csv` -- comma-delimited file; primary dataset used to fit full model
- `good_cells.rds` -- R binary data file; a vector of the raster cell indices that fall within the northern range and should be used for analysis
- `NR_elk_SSM_1988-2020.csv` -- comma-delmited file; results of a state-space model used to interpolate elk counts for all years from 1988 -- 2020; updated version of model presented in appendix of Tallian et al. 2017 *Func Ecol* (see main text for full citation)
- `r_template.rds` -- R binary data file; an object of class `RasterLayer` from the `raster` package with the resolution and extent used to process all other raster layers
- `scale_df.csv` -- comma-delimited file; contains means and standard deviations for all variables in `all_data.csv`; used to scale and center all variables before analyses or for model prediction

### Scripts

- `01_fit_model.R` -- uses package `nimble` and custom functions to fit the full model
- `01_fit_null_model.R` -- uses package `nimble` and custom functions to fit the null model
- `02_figures.R` -- loads posterior samples from full model and makes the majority of diagnostic plots and manuscript figures
- `03_spatial_corr.R` -- checks model residuals for spatial autocorrelation using spline correlograms
- `04_pseudoR2.R` -- calculates pseudo-R^2^ by comparing posteriors of full model to null model
- `05_oos_valid.R` -- calculates out-of-sample (OOS) validation metrics
- `06_open_sens_fit.R` -- fits full model to reduced datasets to check potential sensitivity to imperfect detection due to forest cover
- `07_open_sens_fig.R` -- evaluates fitted models from script `06_open_sens_fit.R` to check potential sensitivity to imperfect detection due to forest cover
- `99_demo_chains.R` -- subsamples full posterior samples for inclusion of demo data in this repo
- `99_fun.R` -- helper functions
- `99_nimble_NB_model_function.R` -- contains the model fitting function that can be run in parallel in scripts `01` and `06` for model fitting
