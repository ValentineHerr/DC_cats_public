# B_Random_Thinning_Analysis

This folder contains the data and script to reproduce the SCR Random thinning model with individual covariates analysis of the paper.

Because this analysis is so computer intensive the study area was divided in 4 blocks. The data for each block is given in files `data_block[1-4].Rdata`. In those files, the location of the traps have been altered because it is sensitive information. 

The files containing the MCMC results and trace of these blocks (too large) are available upon request to Valentine Herrmann HerrmannV@si.edu


# RandomThinning_Scripts

The files in RandomThinning_Scripts are adapted from Ben Augustine's code found [here](https://github.com/benaug/RandomThinIDCov).

The functions to save and restart the MCMC states came from [here](http://danielturek.github.io/public/restartingMCMC/restartingMCMC.html)
