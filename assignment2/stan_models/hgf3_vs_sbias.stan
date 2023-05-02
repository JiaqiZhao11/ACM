data {
  int<lower=1> n;
  
  array[n] int yhat;
  
  array[n] int u;
  
  int <lower = 0, upper = 1> prior;
  
}



// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  
  real <lower = 0> theta;
  real omega;
  real <lower = 0> kappa;

}


transformed parameters{
  array[n] real mu2;
  array[n] real sa2;
  array[n] real mu1hat;
  array[n] real da;
  array[n] real sa1hat;
  array[n] real sa2hat;
  array[n] real mu3;
  array[n] real sa3;
  array[n] real pi3hat;
  array[n] real pi3;
  array[n] real da2;
  array[n] real r2;
  array[n] real w2;
  
  mu2[1] = 0;
  sa2[1] = 2;
  mu1hat[1] = 0.5;
  sa1hat[1] = 0;
  sa2hat[1] = 0;
  da[1] = 0;
  
  mu3[1] = 0;
  sa3[1] = 2;
  pi3hat[1] = 0;
  pi3[1] = 1/sa3[1];
  
  da2[1] = 0;
  r2[1] = 0;
  w2[1] = 0;
  
  for  (t in 2:n){
    
    sa2hat[t] = sa2[t-1]+exp(kappa*mu3[t-1]+omega);
    mu1hat[t] = 1/(1+exp(-mu2[t-1]));
    
    
    sa1hat[t] = mu1hat[t]*(1-mu1hat[t]);
    
    da[t] = u[t]-mu1hat[t];
    
    
    sa2[t] = 1/((1/sa2hat[t])+sa1hat[t]);
    
    mu2[t] = mu2[t-1]+da[t]*sa2[t];
    
    
    da2[t] = ((sa2[t]+(mu2[t]-mu2[t-1])^2)/(sa2[t-1]+exp(kappa*mu3[t-1]+omega)))-1;
    
    r2[t] = (exp(kappa*mu3[t-1]+omega)-sa2[t-1])/(sa2[t-1]+exp(kappa*mu3[t-1]+omega));
    
    w2[t] = exp(kappa*mu3[t-1]+omega)/(sa2[t-1]+exp(kappa*mu3[t-1]+omega));
    
    pi3hat[t] = 1/(sa3[t-1]+theta);
    
    pi3[t] = pi3hat[t]+(kappa^2/2)*w2[t]*(w2[t]+r2[t]*da2[t]);
    
    
    sa3[t] = 1/pi3[t];
    
    mu3[t] = mu3[t-1]+sa3[t]*(kappa/2)*w2[t]*da2[t];
    
  }
  
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  
  target += lognormal_lpdf(theta | -3,0.5);
  target += normal_lpdf(omega | -3,3);
  target +=lognormal_lpdf(kappa | 0,1);
  
  
  if(prior == 0){
  
    for (t in 1:n){
    
      target += bernoulli_lpmf(yhat[t] | mu1hat[t]);
    }
  }
  
}

generated quantities{

  array[n] int sim_yhatt;
  
  real <lower = 0 > theta_prior = lognormal_rng(0,1);
  
  real omega_prior = normal_rng(0,1);
  
  real <lower = 0> kappa_prior = lognormal_rng(0,1);
  
 
}


