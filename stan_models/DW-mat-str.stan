// DW-Mat model fitted using RStan

data { // data block
  int<lower=1> N;     // rows of data (number of observations)
  int<lower=0, upper=1> mat[N];  // Binary outcome variable (0 or 1)
  int<lower=0> DW[N];  // response (independent variable)
  
   int <lower = 0> K; //  number of probabilistic predictions required from the model
   vector <lower = 0> [K] DW_pred; // specifying DW values to predict for

   real Lmat; // known length at 50% maturity
}

parameters { // parameters block
  real beta; // slope (coefficient for DW)
  real alpha; // intercept 
}

model { // model block
  for (i in 1:N) {
    mat[i] ~ bernoulli_logit(alpha + beta * DW[i]);
  }
  
  // priors
  beta ~ normal(0,10); // prior for b
  // Informative prior on alpha based on known length at maturity
  alpha ~ normal(logit(0.5) + beta * Lmat, 10); // Assuming a logistic link function for the prior and adjusting alpha based on the known length
}

generated quantities { // need to create posterior predictive samples for maturity and transform values from the logit scale to the probability scale
  vector[K]  mat_prob_pred; // predicted probabilities for specified number of predictions
  
  for (k in 1:K) {
    mat_prob_pred[k] = inv_logit(alpha + beta * DW_pred[k]); // probability maturity values that can be used to create plot
  }
}


