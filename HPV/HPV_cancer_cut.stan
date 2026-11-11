// A joint model for observations of a Binomial count variable
// and a Poisson count variable, whose rate depends on the Binomial probability
data {
	int<lower=0>            n;                      // number of individuals
	int<lower=0>            N[n];                   // size of binomial for each individual
	
	int<lower=0, upper=N>   Z[n];                   // Binomial count
}

parameters {
	// Random effects
	vector<lower=0,upper=1>[n] 		phi;			// probability in Binomial distribution for each individual
}

model {
        // Individual probabilities
        phi ~ uniform(0,1);
  
        // Observations
        Z ~ binomial(N,phi);      // Binomial model for the number of HPV infections, Z
}
