data {
  int<lower=1> trials; //number of trials
  int<lower=1> subjects; //number of subjects
  array[trials,subjects] int rw1;
  array[trials,subjects] int rw2;
  array[trials,subjects] int fb_rw1;
  array[trials,subjects] int fb_rw2;
  
  
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  
  array[subjects] real <lower = 0, upper = 1> bias_1;
  array[subjects] real <lower = 0, upper = 1> bias_2;
  
  array[subjects] real <lower = 0, upper = 1> alpha_1;
  array[subjects] real <lower = 0, upper = 1> alpha_2;
  
  real <lower = 0, upper = 1> bias_1_mu;
  real <lower = 0, upper = 1> bias_1_sd;
  
  real <lower = 0, upper = 1> bias_2_mu;
  real <lower = 0, upper = 1> bias_2_sd;
  
  real <lower = 0, upper = 1> alpha_1_mu;
  real <lower = 0, upper = 1> alpha_1_sd;
  
  real <lower = 0, upper = 1> alpha_2_mu;
  real <lower = 0, upper = 1> alpha_2_sd;
  
}


transformed parameters{
  array[trials,subjects] real <lower = 0, upper = 1> belief_1;
  array[trials,subjects] real <lower = 0, upper = 1> belief_2;
  
  for (s in 1:subjects){
    belief_1[1,s] = bias_1[s];
    belief_2[1,s] = bias_2[s];
      
    for (t in 2:trials){
      belief_1[t,s] = belief_1[t-1,s]+alpha_1[s]*(rw2[t-1,s]-belief_1[t-1,s]);
      belief_2[t,s] = belief_2[t-1,s]+alpha_2[s]*(rw1[t-1,s]-belief_2[t-1,s]);
  }
}
}
// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  target += beta_lpdf(alpha_1_mu | 1,1);
  target += beta_lpdf(alpha_1_sd | 1,1);
  target += beta_lpdf(alpha_2_mu | 1,1);
  target += beta_lpdf(alpha_2_sd | 1,1);
  
  target +=beta_lpdf(bias_1_mu | 1,1);
  target +=beta_lpdf(bias_1_sd | 1,1);
  
  target +=beta_lpdf(bias_2_mu | 1,1);
  target +=beta_lpdf(bias_2_sd | 1,1);
    
  
  for (s in 1:subjects){
    target +=beta_lpdf(bias_1[s] | bias_1_mu, bias_1_sd);
    target +=beta_lpdf(bias_2[s] | bias_2_mu, bias_2_sd);
    
    target +=beta_lpdf(alpha_1[s] | alpha_1_mu,alpha_1_sd);
    target +=beta_lpdf(alpha_2[s] | alpha_2_mu,alpha_2_sd);

    for (t in 1:trials){
    
      target +=bernoulli_lpmf(rw1[t,s] | belief_1[t,s]);
      target +=bernoulli_lpmf(rw2[t,s] | (1-belief_2[t,s]));
      
    }
    
}

}



