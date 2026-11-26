
# Pooling after multiple imputation under posterior skewness

Implements the analyses in the paper: \`\`Pooling after multiple imputation
under posterior skewness”, Phillips,
Christodoulou and Steinsaltz (in preparation).

Code to apply approximate Bayesian pooling, and the diagnostic workflow for Rubin's rules, are available in the R package `abpool`.
The package is available at https://github.com/danphillipsstats/abpool.

Files and folders:
 - HPV - This contains the HPV and cervical cancer correlation code, outputs and instructions for how to download the data
 - SUPPORT - This contains the SUPPORT code and outputs
 - COVID - This contains the COVID-19 vaccine efficacy code and outputs
 - Simulations - This contains the code required to run the simulations, and the outputs from the simulations, including outputted R objects containing more details of the results.
 - Error_bounds.R - This generates Figure S1: Approximate error in uncertainty interval bounds due to posterior skewness or kurtosis.

 ## Data and notes
Data are not included in the repository, although instructions of where to access the data are.
### HPV example
The data is available in the supplementary material of the paper:
Plummer, M. Cuts in Bayesian graphical models. Stat Comput 25, 37–43 (2015). 
https://doi.org/10.1007/s11222-014-9503-z
The data can be extracted from the file using the following links: 
https://static-content.springer.com/esm/art%3A10.1007%2Fs11222-014-9503-z/MediaObjects/11222_2014_9503_MOESM2_ESM.rnw
which defines the data objects nhpv, Npart, ncases and Npop.
### SUPPORT
The data is freely available as part of the Hmisc package, so the code will run without needing to first download the data.
Since we used a large number of imputations, the imputation model takes a long time (hours on my laptop). 
Outputs from this model and ABpool and Rubin's rules are saved, so that the code can be reproduced without rerunning the imputation model.
### COVID-19 vaccine efficacy example
The COV002 trial data is available via request on the Vivli platform https://doi.org/10.25934/PR00010299.
We are yet to verify whether this dataset matches exactly the data we used, so the code may still not run even with this dataset.
Since intermediate R objects are too large to save to GitHub, and for privacy reasons, we do not provide these.
