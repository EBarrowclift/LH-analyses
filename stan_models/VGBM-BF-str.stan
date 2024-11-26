// von Bertalanffy growth model fitted using RStan

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
  real<lower=0, upper=1> k;
  real<lower=0> Linf;
  real<lower=0> L0;
  real<lower=0.7, upper=1.3> kappa;  // distribution of multiplier of Linf
}

transformed parameters { // optional transformed parameters block
  real<lower=0> Linfprop;             // max size multiplied by proportion
  Linfprop = 1970 * kappa; // updated to max size of 1970mm
}

model { // model block
  real mu[N];
  for (i in 1:N) {
    mu[i] = log(Linf * (1 - ((1 - (L0/Linf)) * exp(-k * age[i]))));
    logDW[i] ~ normal(mu[i], sigma);
  }
  
  // priors
  sigma ~ cauchy(0, 30000); // halfCauchy distribution 30000
  k ~ uniform(0, 2);
  kappa ~ gamma(1000, 980);       // updated hyperprior based around mean of 1.02
  Linf ~ normal(Linfprop, 100);
  L0 ~ normal(700, 200); // update mean to min offspring size of 700 mm
}

generated quantities { // add in generated quantities block to compute log likelihood info
  real log_lik[N]; 
  for (i in 1:N) {
    log_lik[i] = normal_lpdf(logDW[i] | log(Linf * (1 - ((1 - (L0/Linf)) * exp(-k * age[i])))), sigma);
  }
}


