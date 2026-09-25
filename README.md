
# Pooling after multiple imputation under posterior skewness

Implements the analyses in the paper: \`\`Pooling after multiple imputation
under posterior skewness”, Phillips,
Christodoulou and Steinsaltz (in preparation).

Code to apply approximate Bayesian pooling, and the diagnostic workflow for Rubin's rules, are available in the R package `abpool`.
The package is available at https://github.com/danphillipsstats/abpool.

Files and folders:
 - HPV - This contains the HPV and cervical cancer correlation code, and instructions for how to download the data
 - SUPPORT - This contains the SUPPORT code
 - Simulations - This contains the code required to run the simulations. Currently, we don't include simulation outputs due to size constraints, we hope to add these in future.
 - Error_bounds.R - This generates Figure S1: Approximate error in uncertainty interval bounds due to posterior skewness or kurtosis.
 
Not currently included in this repository are:
The COVID-19 example. We hope to add this in future, but need to first check whether any sensitive data is contained in the code files.
The R outputs from the simulations. We plan to add these in future, but the outputs are large, so we may generate a smaller summary and add this.