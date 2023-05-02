SimpleBayes_f <- function(Source1, Source2){
  
  outcome <- inv_logit_scaled(0.5 * logit_scaled(Source1) + 0.5 * logit_scaled(Source2))
  
  return(outcome)
  
}



# set up the social conformatity game 
socialconform_f = function(ntrials, bias, p_1, p_2) {

  ntrials = ntrials
  rating1 = array(NA, ntrials)
  Source1 = array(NA, ntrials)
  first_rating = array(NA, ntrials)
  difference = array(NA, ntrials)
  Source2 = array(NA, ntrials)
  group_rating = array(NA, ntrials)
  
  outcome = array(NA, ntrials)
  second = array(NA, ntrials)
  second_rating = array(NA, ntrials)
  
  # bias and precision for the first rating 
  bias = bias
  p_1 = p_1
  # precision for the second rating
  p_2 = p_2

  
  for (i in seq(ntrials)){
    
    rating1[i] = rprop(1, mean = bias, size = p_1)
    first_rating[i] = round(rating1[i] * 7 + 1)
    Source1[i] = first_rating[i] / 9
  
    
    # Generate a random difference from the set {-3, -2, 0, 2, 3}
    difference[i] = sample(c(-3, -2, 0, 2, 3), 1)
    
    # Calculate the random number by adding the difference to the reference number
    group_rating[i] = first_rating[i] + difference[i]
    
    # Ensure the random number falls within the range [1, 8]
    range_max = 8
    range_min = 1
    while (group_rating[i] < range_min || group_rating[i] > range_max) {
      difference[i] = sample(c(-3, -2, 0, 2, 3), 1)
      group_rating[i] = first_rating[i] + difference[i]
    }
    
    Source2[i] = (group_rating[i] ) / 9
    
    outcome[i] = SimpleBayes_f(Source1 = Source1[i], Source2 = Source2[i])
    second[i] = rprop(1, mean = outcome[i], size = p_2)
    second_rating[i] = round(second[i]) / 9
    
  }
  
  df = list(first_rating = first_rating, Source1 = Source1, difference = difference, group_rating = group_rating, Source2 = Source2, outcome = outcome, second = second, second_rating = second_rating)
  
  return(df)
  
  
}



