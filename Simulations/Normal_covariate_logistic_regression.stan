data {
  int<lower=1> n;
  int<lower=0,upper=n> n_obs;
	int<lower=n-n_obs,upper=n-n_obs>    n_mis;

  int<lower=1,upper=n> obs[n_obs];
  int<lower=1,upper=n> mis[n_mis];

  vector[n_obs] X_obs;
  array[n] int<lower=0, upper=1> Y; // Standardised
  vector<lower=0, upper=1>[n] Z;
}
transformed data {
  // Center Z and divide by sd to reduce correlation with the intercept
  vector[n] Zc = (Z - 0.5)*4;
}
parameters {
  real beta_0;
  real beta_X;
  real beta_Z;
  vector[n_mis] X_mis;
}
model {
  vector[n] X;
  X[mis] = X_mis;
  X[obs] = X_obs;
  // Priors
  beta_0 ~ student_t(5,0,1);
  beta_X ~ student_t(5,0,1);
  beta_Z ~ student_t(5,0,1);
  
  // Likelihood
  X_obs ~ normal(0,1);
  X_mis ~ normal(0,1);
  Y ~ bernoulli_logit(beta_0 + beta_X * X + beta_Z * Zc);
}
