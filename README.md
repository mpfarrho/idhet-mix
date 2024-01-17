# Description
These files are to be used with R and did run in version 4.3.0. 

## Data files (sources are described in the paper)
The folder 'data' collects all required datafiles alongside the collected output of the estimation procedure across horizons:

- 'dailydata.rda' is the main daily data
- 'exoshocks.rda' are the additionally observed high-frequency shocks that act as controls
- 'mpu_indexes.rda' is the monetary policy uncertainty series
- 'VARdata.rda' is the monthly data used to estimate the auxiliary VAR model
- 'FACTOR_mix_pub_exo_J30Q2sv_lev_facsv.rda' are posteriors required for the VAR

## Code files
The individual R-source files are the following:

- 'main_estimation.R' runs the model for all required specifications (set in 'spec_grid.R')
- 'main_estVAR.R' runs the auxiliary VAR model
- 'spec_grid.R' sets up a grid for all required estimations and defines algorithm settings
- 'utils.R' contains several helper functions alongside the main estimation procedure
- 'utils_conjVAR.R' contains the VAR codes
- 'main_output.R' produces all charts shown in the paper from the collected MCMC draws

The required MCMC estimation output 'mix_pub_exo_J30Q2sv_lev_facsv.rda' is available [here (Dropbox)](https://www.dropbox.com/scl/fi/bvmris93grldxn8y3tihx/mix_pub_exo_J30Q2sv_lev_facsv.rda?rlkey=bnwaqyq0vh9dlj1va05yo2qix&dl=0).

This code comes without technical support of any kind. Please report any typos or errors to: mpfarrho@gmail.com. The code is free to use, provided that the paper is cited properly.
