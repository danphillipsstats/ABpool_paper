// Exponential model in Stan with missing values in covariate X
data {
  int<lower=1> n;
  int<lower=0,upper=n> n_obs;
  int<lower=n-n_obs,upper=n-n_obs>    n_mis;
  
  int<lower=1,upper=n> obs[n_obs];
  int<lower=1,upper=n> mis[n_mis];
  
  vector[n_obs] log_X_obs_unstd;
  vector[n] t; // event times
  int<lower=0, upper=1> event[n]; // event indicator
  vector<lower=0, upper=1>[n] Z;
}
transformed data {
  // Center Z and divide by sd to reduce correlation with the intercept
  vector[n] Zc = (Z - 0.5)*4;
}
parameters {
  // X lognormal sd
  real log_ln_sigma; // standard deviation of lognormal variable
  // parameters
  real log_h_0;
  real beta_X;
  real beta_Z;
  vector[n_mis] std_log_X_mis_unstd;  // predictions for (log( X_mis*sd(ln(mu,sigma)) + E(ln(mu,sigma)) ) - mu)/sigma = (log( X_mis*sd(ln(0,sigma)) + E(ln(0,sigma)) ))/sigma
}
transformed parameters{
  real ln_sigma = exp(log_ln_sigma);
  real log_mean_ln = 0.5*square(ln_sigma);  // log of the mean of a lognormal(0,sigma^2) variable (since we standardised X ~ ln(0,sigma^2) by subtracting mean and dividing by sd)
  real<lower=0> denominator = sqrt(expm1(square(ln_sigma))); // sd of a lognormal(0,sigma^2) variable
  vector[n_mis]  X_mis = (expm1(std_log_X_mis_unstd*ln_sigma - log_mean_ln))/denominator;
}
model {
  vector[n] X;
  vector[n] log_haz; // log sigma_i
  vector[n] hazard;
  X[obs] = expm1(log_X_obs_unstd-log_mean_ln)/denominator;
  X[mis] = X_mis;
  
  log_haz = log_h_0 + beta_X * X + beta_Z * Zc; // constant baseline hazard
  hazard = exp(log_haz);
  
  // Priors
  log_h_0 ~ normal(0, 1);
  beta_X ~ normal(0, 1);
  beta_Z ~ normal(0, 1);
  // Priors for X lognormal sd
  log_ln_sigma ~ normal(log(0.5),1); // 
  
  // Likelihood
  target += normal_lpdf(log_X_obs_unstd | 0, ln_sigma);
  std_log_X_mis_unstd ~ normal(0,1); // Prior for X_mis (the generating distribution before Y is generated)
  
  // Likelihood with right-censoring
  for (i in 1:n) {
    if (event[i] == 1) {
      target += log_haz[i] - hazard[i] * t[i];
    } else {
      target += - hazard[i] * t[i];
    }
  }
}
