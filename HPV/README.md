# HPV and cervical cancer correlation
This recreates an analysis from Plummer (2014) Cuts in Bayesian graphical models. Statistics and Computing.
The analysis estimates the association between human papillomavirus prevalence and cervical cancer incidence, 
itself recreating part of an analysis from 
Maucort-Boulch, D., Franceschi, S., Plummer, M.: International correlation between human papillomavirus prevalence and cervical cancer incidence. Cancer Epidemiol. Biomark. Prev. 17, 717–720 (2008).
Plummer compares Rubin's rules to various methods to fit the Bayesian cut model (which we refer to as Bayesian inference after multiple imputation).
In our analysis, we compare Rubin's rules and ABpool in their ability to recreate the results from Bayesian inference after multiple imputation.

## Computation time
Fitting Bayesian inference after multiple imputation takes a few minutes to run, so an intermediate object (Plummer_BayesMI_thetas.RDS) is provided to speed up re-analysis.
