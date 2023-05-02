data {
  int<lower=1> n;
  array[n] int rw1;
  array[n] int rw2;
  array[n] int fb_rw1;
  array[n] int fb_rw2;
  int <lower = 0, upper = 1> prior;
  
  
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  
  real  bias_1;
  real  bias_2;
  
  real  alpha_1;
  real  alpha_2;

}


transformed parameters{
  array[n] real <lower = 0, upper = 1> belief_1;
  array[n] real <lower = 0, upper = 1> belief_2;
  

  belief_1[1] = inv_logit(bias_1);
  belief_2[1] = inv_logit(bias_2);
  
  for (i in 2:n){
    belief_1[i] = belief_1[i-1]+inv_logit(alpha_1)*(rw2[i-1]-belief_1[i-1]);
    belief_2[i] = belief_2[i-1]+inv_logit(alpha_2)*(rw1[i-1]-belief_2[i-1]);
  }
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  
  target += normal_lpdf(bias_1 | 0,1);
  target += normal_lpdf(bias_2 | 0,1);
  
  target +=normal_lpdf(alpha_1 | 0,1);
  target +=normal_lpdf(alpha_2 | 0,1);
  
  
  if(prior == 0){
    
    for (i in 1:n){
    
      target +=bernoulli_lpmf(rw1[i] | belief_1[i]);
      target +=bernoulli_lpmf(rw2[i] | (1-belief_2[i]));
      
    }
  }
}


generated quantities{

  array[n] int sim_rw1t;
  array[n] int sim_rw2t;
  
  real <lower = 0, upper = 1> theta1_prior = inv_logit(bias_1);
  real <lower = 0, upper = 1> theta2_prior = inv_logit(bias_2);
  
  real <lower = 0, upper = 1> alpha1_prior = inv_logit(alpha_1);
  real <lower = 0, upper = 1> alpha2_prior = inv_logit(alpha_2);
  
  
  for (i in 1:n){
    sim_rw1t[i] = bernoulli_rng(belief_1[i]);
    sim_rw2t[i] = bernoulli_rng(1-belief_2[i]);
  }
  
  int sim_rw1 = sum(sim_rw1t);
  int sim_rw2 = sum(sim_rw2t);
  
}

