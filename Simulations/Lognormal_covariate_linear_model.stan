// A linear regression with a partly missing-completely-at-random continuous covariate, and a binary covariate.
data {
	int<lower=1>            n;                      // number of individuals
	int<lower=0,upper=n>    n_obs;                  // number of individuals with observed X
	int<lower=n-n_obs,upper=n-n_obs>    n_mis;      // number of individuals with missing X
	
  vector[n_obs] log_X_obs_unstd;
	                        
	int<lower=1,upper=n>    obs[n_obs];             // Indices of individuals with observed X
	int<lower=1,upper=n>    mis[n_mis];             // Indices of individuals with missing X
	
	vector<lower=0, upper=1>[n]   Z;                   // Indicator
	vector[n]               Y;                      // Outcome
}
transformed data {
  // Center Z and divide by sd to reduce correlation with the intercept
  vector[n] Zc = (Z - 0.5)*4;
}
parameters {
  // X lognormal sd
  real log_ln_sigma; // standard deviation of lognormal variable
	// Parameters
	real                          beta_0;      // Intercept in regression of Y ~ 1 + X + Z
	real                          beta_X;      // X coefficient in regression of Y ~ 1 + X + Z
	real                          beta_Z;      // Z coefficient in regression of Y ~ 1 + X + Z
	real<lower=0>                 sigma;     // Variance parameter in Y ~ 1 + X + Z
	// Missing values
  vector[n_mis] std_log_X_mis_unstd;  // predictions for (log( X_mis*sd(ln(mu,sigma)) + E(ln(mu,sigma)) ) - mu)/sigma
}
transformed parameters{
  real ln_sigma = exp(log_ln_sigma);
  real log_mean_ln = 0.5*square(ln_sigma);  // log of the mean of a lognormal(0,sigma^2) variable (since we standardised X ~ ln(0,sigma^2) by subtracting mean and dividing by sd)
  real<lower=0> denominator = sqrt(expm1(square(ln_sigma))); // sd of a lognormal(0,sigma^2) variable
  vector[n_mis]  X_mis = (expm1(std_log_X_mis_unstd*ln_sigma - log_mean_ln))/denominator;
}
model {
  vector[n] X;
  X[obs] = expm1(log_X_obs_unstd-log_mean_ln)/denominator;
  X[mis] = X_mis;
  
  // Priors
  beta_0 ~ normal(0,1);
  beta_X ~ normal(0,2);
  beta_Z ~ normal(0,1);
  sigma ~ normal(0,1);
  // Priors for X lognormal sd
  log_ln_sigma ~ normal(log(0.5),1); // 
    
  // Predict X_mis
  target += normal_lpdf(log_X_obs_unstd | 0, ln_sigma);
  std_log_X_mis_unstd ~ normal(0,1); // Prior for X_mis (the generating distribution before Y is generated)
  
  Y ~ normal(beta_0 + X*beta_X + Zc*beta_Z, sigma);
}
