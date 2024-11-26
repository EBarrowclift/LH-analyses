// Gompertz growth model

data { // data block
  int<lower=1> N;     // rows of data
  int<lower=0> age[N];  // predictor
  int<lower=0> DW[N];  // response
}

transformed data { // optional transformed data block
  real logDW[N];
  for (i in 1:N) {
    logDW[i] = log(DW[i]);
  }
}

parameters { // parameters block
  real<lower=0> sigma; // error term
  real<lower=0, upper=1> k; // constrained between 0 and 1
  real<lower=0> Linf;
  real<lower=0> L0;
  real<lower=0.7, upper=1.3> kappa;  // distribution of multiplier of Linf
}

transformed parameters { // optional transformed parameters block
  real<lower=0> Linfprop;             // max size multiplied by proportion
  Linfprop = 3500 * kappa; // updated to max size of 3500m instead of 3100mm
}

model { // model block
  real mu[N];
  for (i in 1:N) {
    mu[i] = log(L0*exp(log(Linf/L0)*(1-exp(-k * age[i])))); // check equation correct
    logDW[i] ~ normal(mu[i], sigma);
  }
  
  // priors
  sigma ~ cauchy(0, 30000); // halfCauchy distribution 30000
  k ~ uniform(0, 2);
  kappa ~ gamma(200, 198);       // hyperprior 
  Linf ~ normal(Linfprop, 400);
  L0 ~ normal(900, 300);
}

generated quantities { // add in generated quantities block to compute log likelihood info
  real log_lik[N]; 
  for (i in 1:N) {
    log_lik[i] = normal_lpdf(logDW[i] | log(L0*exp(log(Linf/L0)*(1-exp(-k * age[i])))), sigma);
  }
}


