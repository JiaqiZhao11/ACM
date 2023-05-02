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
  real <lower = 0, upper = 1> bias_1_mu;
  real <lower = 0, upper = 1> bias_1_sd;
  real <lower = 0, upper = 1> bias_2_mu;
  real <lower = 0, upper = 1> bias_2_sd;
  array[subjects] real <lower = 0, upper = 1> bias_1;
  array[subjects] real <lower = 0, upper = 1> bias_2;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  target +=beta_lpdf(bias_1_mu | 1,1);
  target +=beta_lpdf(bias_1_sd | 1,1);
  
  target +=beta_lpdf(bias_2_mu | 1,1);
  target +=beta_lpdf(bias_2_sd | 1,1);
    
  
  for (s in 1:subjects){
    target +=beta_lpdf(bias_1[s] | bias_1_mu,bias_1_sd);
    target +=beta_lpdf(bias_2[s] | bias_2_mu,bias_2_sd);
    
    for (t in 1:trials){
    
      target +=bernoulli_lpmf(rw1[t,s] | bias_1[s]);
      target +=bernoulli_lpmf(rw2[t,s] | (1-bias_2[s]));
      
    }
  }
}
