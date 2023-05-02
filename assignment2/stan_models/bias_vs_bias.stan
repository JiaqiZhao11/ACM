data {
  int<lower=1> n;
  int<lower=0, upper = 1> prior;
  array[n] int rw1;
  array[n] int rw2;
  
}

// The parameters accepted by the model. Our model
// accepts two parameters 'muwww' and 'sigma'.
parameters {
  real  bias_1;
  real  bias_2;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  target +=normal_lpdf(bias_1 | 0,1);
  target +=normal_lpdf(bias_2 | 0,1);
  

  if(prior == 0){  
    for (i in 1:n){
  
        target +=bernoulli_lpmf(rw1[i] |inv_logit(bias_1));
        target +=bernoulli_lpmf(rw2[i] |1-inv_logit(bias_2));
      }
    
  }
  
}


generated quantities{

  int sim_rw1;
  int sim_rw2;
 
  real <lower = 0, upper = 1> theta1_prior = inv_logit(bias_1);
  real <lower = 0, upper = 1> theta2_prior = inv_logit((bias_2));
  

  
  sim_rw1 = binomial_rng(n,theta1_prior);
  sim_rw2 = binomial_rng(n,(1-theta2_prior));


}
