data {
  int<lower=1> trials; //number of trials
  int<lower=1> subjects; //number of subjects
  int<lower=0, upper = 1> prior;
  array[trials,subjects] int rw1;
  array[trials,subjects] int rw2;
  array[trials,subjects] int fb_rw1;
  array[trials,subjects] int fb_rw2;
  
  
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real bias_1_mu;
  real bias_1_sd;
  real bias_2_mu;
  real bias_2_sd;
  array[subjects] real bias_1;
  array[subjects] real bias_2;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  target +=normal_lpdf(bias_1_mu | 0,1);
  target +=normal_lpdf(bias_1_sd | 0,1);
  
  target +=normal_lpdf(bias_2_mu | 0,1);
  target +=normal_lpdf(bias_2_sd | 0,1);
    
  
  for (s in 1:subjects){
    target +=normal_lpdf(bias_1[s] | bias_1_mu, bias_1_sd);
    target +=normal_lpdf(bias_2[s] | bias_2_mu, bias_2_sd);
    
    if (prior == 0){
      for (t in 1:trials){
        target +=bernoulli_lpmf(rw1[t,s] | inv_logit(bias_1[s]));
        target +=bernoulli_lpmf(rw2[t,s] | 1-inv_logit(bias_2[s]));
        
      }
    }
  }
}


generated quantities{
  
  
  array[subjects] int sim_rw1;
  array[subjects] int sim_rw2;
  
  real <lower = 0, upper = 1> theta1_prior_p = inv_logit(bias_1_mu);
  real <lower = 0, upper = 1> theta2_prior_p = inv_logit(bias_2_mu); 
  
  array[subjects] real <lower = 0, upper = 1> theta1_prior;
  array[subjects] real <lower = 0, upper = 1> theta2_prior; 
 
  for (s in 1:subjects){
    theta1_prior[s] = inv_logit(bias_1[s]);
    theta2_prior[s] = inv_logit(bias_2[s]);
      
    sim_rw1[s] = binomial_rng(trials,theta1_prior[s]);
    sim_rw2[s] = binomial_rng(trials,1-theta2_prior[s]);
        
  
    
  }




}

