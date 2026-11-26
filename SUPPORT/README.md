# SUPPORT
The data is freely available as part of the Hmisc package, so the code will run without needing to first download the data.
Since we used a large number of imputations, the imputation model takes a long time (hours on my laptop) - this is from line 53 in the code.
Outputs from this model and ABpool and Rubin's rules are saved, so that the code can be reproduced without rerunning the imputation model.
Saved intermediate objects:
imp_data_support_age_numco_scoma_alb_bili_pafi_crea_meanbp_hrt_resp_temp_adlsc_dzgroup_NACH_100000.RDS
^contains the imputations, the output of mice::mice()
cox_mice_imp_data_support_age_numco_scoma_alb_bili_pafi_crea_meanbp_hrt_resp_temp_adlsc_dzgroup_NACH_100000.RDS
^contains the fitted cox models, from running with() on the above object
SUPPORT_output_objects.RDS 
^contains the results from ABpool and Rubin's rules, as well as the estimates and variances for each imputation.
