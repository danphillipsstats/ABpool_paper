#####
# These create a list of data to be inputted into the stan models

# Creates a list of data to input into Longitudinal_hierarchical_model_normalised_ri_rs_exp_slope_t.stan
stan_long_model <- function(long_data, log_y_var, obs_time_var, individual_data, 
                            longitudinal_covariates, ind_var, t_center = "mean"){
  n <- nrow(individual_data)
  
  X_l_df <- as.data.frame(individual_data[,longitudinal_covariates])
  X_l <- as.matrix(model.matrix(~.,data=X_l_df)[,-1],ncol=length(longitudinal_covariates))
  ll <- match(long_data[,c(ind_var)],
              individual_data[,c(ind_var)],
              incomparables = NA) # ll[i] is the individual corresponding to the i th antibody observation
  n_l <- nrow(long_data)
  p_l <- ncol(X_l)
  if (t_center=="mean") {t_center <- mean(long_data[,obs_time_var])}
  
  data <- list("n" =n, # number of individuals
               "n_l" = n_l, # number of antibody observations
               "p_l" = p_l, # number of covariates in X_l
               "ll" = ll,
               "X_l" = X_l, # covariate matrix
               "log_y" = long_data[,log_y_var], # the observed log antibodies
               "t_l" = long_data[,obs_time_var], # the timings of the observations
               "t_center" = t_center) # the transformed time 0 in the running of the stan model to improve mixing
  return(data)
}

# Gauss-Kronrod nodes and weights
q_nodes.unadjusted <- c(
  -0.991455371120813,
  -0.949107912342758,
  -0.864864423359769,
  -0.741531185599394,
  -0.586087235467691,
  -0.405845151377397,
  -0.207784955007898,
  0,
  0.207784955007898,
  0.405845151377397,
  0.586087235467691,
  0.741531185599394,
  0.864864423359769,
  0.949107912342758,
  0.991455371120813
)
q_weights.unadjusted <- c(
  0.022935322010529,
  0.063092092629979,
  0.10479001032225,
  0.140653259715526,
  0.169004726639268,
  0.190350578064785,
  0.204432940075299,
  0.209482141084728,
  0.204432940075299,
  0.190350578064785,
  0.169004726639268,
  0.140653259715526,
  0.10479001032225,
  0.063092092629979,
  0.022935322010529
)
q <- length(q_nodes.unadjusted)

myremove_characters <- function(x){
  x <- str_replace_all(x,"[^[:alnum:]]", "")
  return(x)
}

# Takes a vector with a date in (almost) any format and returns it in the ISO
# format
myfix_dates <- function(x, years = c("2020","2021","2022","2023")){
  Datevec <- x
  Datevec <- myremove_characters(Datevec)
  Datevec[which(Datevec=="")] <- NA
  Datevec_start <- substring(Datevec[which(!is.na(Datevec))[1]], first = 1, last = 4)
  if(is.element(Datevec_start, years)){
    Years <- substring(Datevec, first = 1, last = 4) 
    Months <- substring(Datevec, first = 5, last = 6)
    Days <- substring(Datevec, first = 7, last = 8)
    Datevec <- paste(Years,Months,Days,
                     sep = "-")
    Datevec[which(Datevec=="NA-NA-NA")] <- NA
  } else{
    Years <- substring(Datevec, first = 5, last = 8) 
    Months <- substring(Datevec, first = 3, last = 4)
    Days <- substring(Datevec, first = 1, last = 2)
    Datevec <- paste(Years,Months,Days,
                     sep = "-")
    Datevec[which(Datevec=="NA-NA-NA")] <- NA
  }
  return(Datevec)
}
