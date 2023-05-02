data {
  int<lower=1> n;
  array[n] int rw1;
  array[n] int rw2;
  array[n] int fb_rw1;
  array[n] int fb_rw2;
  int <lower = 0, upper = 1> prior;
  real bias1_mean; 
  real bias2_mean; 
  real bias1_sd; 
  real bias2_sd; 

  real alpha_1w_mean; 
  real alpha_2w_mean; 
  real alpha_1w_sd; 
  real alpha_2w_sd; 
  
  
  real alpha_1l_mean; 
  real alpha_2l_mean; 
  real alpha_1l_sd; 
  real alpha_2l_sd; 
  

  
}


// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  
  
  real  bias_1;
  real  bias_2;
  
  real alpha_1w;
  real alpha_1l;
  
  real alpha_2w;
  real alpha_2l;
}


transformed parameters{
  array[n] real <lower = 0, upper = 1> belief_1;
  array[n] real <lower = 0, upper = 1> belief_2;
  

  belief_1[1] = inv_logit(bias_1);
  belief_2[1] = inv_logit(bias_2);
  
  for (i in 2:n){
    if(fb_rw1[i-1])
      belief_1[i] = belief_1[i-1]+inv_logit(alpha_1w)*(rw2[i-1]-belief_1[i-1]);
    else
      belief_1[i] = belief_1[i-1]+inv_logit(alpha_1l)*(rw2[i-1]-belief_1[i-1]);
      
    if(fb_rw2[i-1])
      belief_2[i] = belief_2[i-1]+inv_logit(alpha_2w)*(rw1[i-1]-belief_2[i-1]);
    else
      belief_2[i] = belief_2[i-1]+inv_logit(alpha_2l)*(rw1[i-1]-belief_2[i-1]);
  }
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  
  
  target +=normal_lpdf(bias_1 | bias1_mean,bias1_sd);
  target +=normal_lpdf(bias_2 | bias2_mean,bias2_sd);
  
  target +=normal_lpdf(alpha_1w | alpha_1w_mean,alpha_1w_sd);
  target +=normal_lpdf(alpha_1l | alpha_1l_mean,alpha_1l_sd);
  
  target +=normal_lpdf(alpha_2w | alpha_2w_mean,alpha_2w_sd);
  target +=normal_lpdf(alpha_2l | alpha_2l_mean,alpha_2l_sd);
  
  
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
  
  real <lower = 0, upper = 1> alpha1l_prior = inv_logit(alpha_1l);
  real <lower = 0, upper = 1> alpha1w_prior = inv_logit(alpha_1w);
  
  real <lower = 0, upper = 1> alpha2l_prior = inv_logit(alpha_2l);
  real <lower = 0, upper = 1> alpha2w_prior = inv_logit(alpha_2w);
  
  
  for (i in 1:n){
    sim_rw1t[i] = bernoulli_rng(belief_1[i]);
    sim_rw2t[i] = bernoulli_rng(1-belief_2[i]);
  }
  
  int sim_rw1 = sum(sim_rw1t);
  int sim_rw2 = sum(sim_rw2t);
  
}
