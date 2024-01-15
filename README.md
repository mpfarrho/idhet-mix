These files are to be used with R and did run in version 4.3.0.

The folder 'data' collects all required datafiles alongside the collected output of the estimation procedure across horizons:

- 'dailydata.rda' is the main daily data
- 'exoshocks.rda' are the additionally observed high-frequency shocks that act as controls
- 'mpu_indexes.rda' is the monetary policy uncertainty series
- 'VARdata.rda' is the monthly data used to estimate the auxiliary VAR model

- 'mix_pub_exo_J30Q2sv_lev_facsv.rda' is the collected MCMC estimation output
- 'FACTOR_mix_pub_exo_J30Q2sv_lev_facsv.rda' are posteriors required for the VAR

The individual R-source files are the following:

- 'main_estimation.R' runs the model for all required specifications (set in 'spec_grid.R')
- 'main_estVAR.R' runs the auxiliary VAR model
- 'spec_grid.R' sets up a grid for all required estimations and defines algorithm settings
- 'utils.R' contains several helper functions alongside the main estimation procedure
- 'ugils_conjVAR.R' contains the VAR codes

- 'main_output.R' produces all charts shown in the paper from the collected MCMC draws
