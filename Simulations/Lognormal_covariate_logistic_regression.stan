data {
  int<lower=1> n;
  int<lower=0,upper=n> n_obs;
	int<lower=n-n_obs,upper=n-n_obs>    n_mis;

  int<lower=1,upper=n> obs[n_obs];
  int<lower=1,upper=n> mis[n_mis];

  vector[n_obs] log_X_obs_unstd; // log of the unstandardised observed X
  
  array[n] int<lower=0, upper=1> Y; // Standardised
  vector<lower=0, upper=1>[n] Z;
}
transformed data {
  // Center Z and divide by sd to reduce correlation with the intercept
  vector[n] Zc = (Z - 0.5)*4;
  real mean_ln01 = exp(0.5);  // mean of a lognormal(0,1) variable (since we standardised X ~ ln(0,1) by subtracting mean and dividing by sd)
  real sd_ln01 = sqrt(exp(1)*(exp(1)-1)); // sd of a lognormal(0,1) variable
}
parameters {
  // X lognormal sd
  real log_ln_sigma; // standard deviation of lognormal variable
  // Parameters
  real beta_0;
  real beta_X;
  real beta_Z;
	// Missing values
  vector[n_mis] std_log_X_mis_unstd;  // predictions for (log( X_mis*sd(ln(mu,sigma)) + E(ln(mu,sigma)) ) - mu)/sigma
}
transformed parameters {
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
  beta_0 ~ student_t(5,0,1);
  beta_X ~ student_t(5,0,1);
  beta_Z ~ student_t(5,0,1);
  // Priors for X lognormal sd
  log_ln_sigma ~ normal(log(0.5),1); // 
  
  // Predict X_mis
  target += normal_lpdf(log_X_obs_unstd | 0, ln_sigma);
  std_log_X_mis_unstd ~ normal(0,1); // Prior for X_mis (the generating distribution before Y is generated)
  
  Y ~ bernoulli_logit(beta_0 + beta_X * X + beta_Z * Zc);
}
