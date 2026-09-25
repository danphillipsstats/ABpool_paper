// Exponential model in Stan with missing values in covariate X
data {
  int<lower=1> n;
  int<lower=0,upper=n> n_obs;
  int<lower=n-n_obs,upper=n-n_obs>    n_mis;
  
  int<lower=1,upper=n> obs[n_obs];
  int<lower=1,upper=n> mis[n_mis];
  
  vector[n_obs] X_obs;
  vector[n] t; // event times
  int<lower=0, upper=1> event[n]; // event indicator
  vector<lower=0, upper=1>[n] Z;
}
transformed data {
  // Center Z and divide by sd to reduce correlation with the intercept
  vector[n] Zc = (Z - 0.5)*4;
}
parameters {
  real log_h_0;
  real beta_X;
  real beta_Z;
  vector[n_mis] X_mis;
}
model {
  vector[n] X;
  X[mis] = X_mis;
  X[obs] = X_obs;
  
  vector[n] log_haz; // log sigma_i
  log_haz = log_h_0 + beta_X * X + beta_Z * Zc; // constant baseline hazard
  vector[n] hazard = exp(log_haz);
  
  // Priors
  log_h_0 ~ normal(0, 1);
  beta_X ~ normal(0, 1);
  beta_Z ~ normal(0, 1);
  
  // Likelihood
  X_obs ~ normal(0,1);
  X_mis ~ normal(0,1);
  
  // Likelihood with right-censoring
  for (i in 1:n) {
    if (event[i] == 1) {
      target += log_haz[i] - hazard[i] * t[i];
    } else {
      target += - hazard[i] * t[i];
    }
  }
}
