// DW-Wt model fitted using RStan

data { // data block
  int<lower=1> N;     // rows of data
  int<lower=0> Wgt[N];  // predictor
  int<lower=0> DW[N];  // response
}

transformed data { // optional transformed data block for log transformation
  real logWgt[N];
  for (i in 1:N) {
    logWgt[i] = log(Wgt[i]);
  }
}

parameters { // parameters block
  real<lower=0> sigma; // error term
  real b; // slope (parameter b)
  real loga; // intercept (parameter log(a))
}

model { // model block
  real mu[N];
  for (i in 1:N) {
    mu[i] = loga + (b * log(DW[i]));
    logWgt[i] ~ normal(mu[i], sigma);
  }
  
  // priors
  loga ~ normal(-5,3); // prior for a
  b ~ normal(3,1); // prior for b
  sigma ~ cauchy(0, 30000); // halfCauchy distribution 30000
}


