// A joint model for observations of a Binomial count variable
// and a Poisson count variable, whose rate depends on the Binomial probability
data {
	int<lower=0>            n;                      // number of individuals
	vector<lower=0>[n]      t;                   // time at-risk for each individual
	
  int<lower=0>            Y[n];                   // Poisson count
  // Imputed random effects
	vector<lower=0,upper=1>[n] 		phi;			// probability in Binomial distribution for each individual
}
transformed data{
  vector[n] 		phi_std;			// standardised
  vector<lower=0>[n]      t_std;                   // standardised
	
  phi_std = (phi-mean(phi))/sd(phi);
  t_std = t/mean(t);
}

parameters {
	// Regression coefficients
	real    		theta_1_std;		               	// intercept in regression for Poisson rate
	real    		theta_2_std;		               	// phi effect in regression for Poisson rate
	
}

transformed parameters{
  vector<lower=0>[n] 		lambda_t;			// Poisson mean for each individual
	lambda_t = exp(theta_1_std + theta_2_std*phi_std + log(t_std));
}

model {
        // Priors
        theta_1_std ~ normal(0,sqrt(1000));
        theta_2_std ~ normal(0,sqrt(1000));
  
        // Observations
        Y ~ poisson(lambda_t);    // Poisson model for the number of cancer cases, Y
}

generated quantities{
  real theta_1;
  real theta_2;
  
  theta_1 = theta_1_std - log(mean(t)) - theta_2_std*mean(phi)/sd(phi);
  theta_2 = theta_2_std/sd(phi);
}
