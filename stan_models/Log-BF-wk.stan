// Logistic growth model fitted using RStan

data { // data block
  int<lower=1> N;     // rows of data
  int<lower=0> age[N];  // predictor x
  int<lower=0> DW[N];  // response y
}

transformed data { // optional transformed data block
  real logDW[N];
  for (i in 1:N) {
    logDW[i] = log(DW[i]);
  }
}

parameters { // parameters block
  real<lower=0> sigma; // error term
  real<lower=0, upper=1> k;
  real<lower=0> Linf;
  real<lower=0> L0;
  real<lower=0.7, upper=1.3> kappa;  // distribution of multiplier of Linf
}

transformed parameters { // optional transformed parameters block
  real<lower=0> Linfprop;             // max size multiplied by proportion
  Linfprop = 1970 * kappa; // updated to max size of 3500m instead of 3100mm
}

model { // model block
  real mu[N];
  for (i in 1:N) {
    mu[i] = log((Linf*L0*exp(k * age[i])) / (Linf + L0 * (exp(k * age[i]) - 1))); // check correct equation
    logDW[i] ~ normal(mu[i], sigma);
  }
  
  // priors
  sigma ~ cauchy(0, 30000); // halfCauchy distribution 30000
  k ~ uniform(0, 2);
  kappa ~ gamma(200, 196);       // hyperprior 
  Linf ~ normal(Linfprop, 400);
  L0 ~ normal(700, 300);
}

generated quantities { // add in generated quantities block to compute log likelihood info
  real log_lik[N]; 
  for (i in 1:N) {
    log_lik[i] = normal_lpdf(logDW[i] | log((Linf*L0*exp(k * age[i])) / (Linf + L0 * (exp(k * age[i]) - 1))), sigma);
  }
}


