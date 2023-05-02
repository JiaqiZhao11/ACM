
// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  vector[N] Source1;
  vector[N] Source2;
  vector[N] second;
}


transformed data{
  vector[N] l_Source1;
  vector[N] l_Source2;
  l_Source1 = logit(Source1) * 0.5;
  l_Source2 = logit(Source2) * 0.5;
  
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real bias;
  real<lower=0> p_1;
  real<lower=0> p_2;
}

// The model to be estimated.
model {
  target += beta_proportion_lpdf(bias | 0.5, 2);
  target += lognormal_lpdf(p_1 | 0.5, 1);
  target += lognormal_lpdf(p_2 | 0.5, 1);
  target += beta_proportion_lpdf(Source1 | bias, p_1);
  target += beta_proportion_lpdf(second | inv_logit(l_Source1 + l_Source2), p_2);
  
}

