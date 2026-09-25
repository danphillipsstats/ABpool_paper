// A linear regression model, where Y = X beta_X_std + Z beta_Z
// and X is sampled from N(0,1) (assumed known), and missing completely at random.
data {
  int<lower=1> n;
  int<lower=0,upper=n> n_obs;
  int<lower=n-n_obs,upper=n-n_obs> n_mis;

  int<lower=1,upper=n> obs[n_obs];
  int<lower=1,upper=n> mis[n_mis];

  vector[n_obs] X_obs;
  vector[n] Y; // Standardised
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
  real<lower=0> sigma;
}
model {
  // Priors
  beta_0 ~ normal(0, 1);
  beta_X ~ normal(0, 2); // Wider prior to avoid over-shrinking key parameter
  beta_Z ~ normal(0, 1);
  sigma  ~ normal(0, 1);
  
  // Likelihood for observed-X rows
  {
    X_obs ~ normal(0,1);
    vector[n_obs] mu_obs = beta_0 + beta_X * X_obs + beta_Z * Zc[obs];
    Y[obs] ~ normal(mu_obs, sigma);
  }
  // Likelihood for missing-X rows (marginalised over missing X)
  {
    vector[n_mis] mu_mis = beta_0 + beta_Z * Zc[mis];
    real sigma_mis = sqrt(square(sigma) + square(beta_X));
    Y[mis] ~ normal(mu_mis, sigma_mis);
  }
}
generated quantities {
  // Posterior draws of X for missing rows (imputations)
  vector[n_mis] X_mis;
  for (k in 1:n_mis) {
    int i = mis[k];
    real s2 = 1 / (1 + square(beta_X) / square(sigma));
    real m  = s2 * (beta_X / square(sigma)) *
               (Y[i] - beta_0 - beta_Z * Z[i]); 
    X_mis[k] = normal_rng(m, sqrt(s2));
  }
}
