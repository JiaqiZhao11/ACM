# Distance 
distance <- function(vect1, vect2, w) {
  return(sum(w * abs(vect1 - vect2)))
}

# Similarity
similarity <- function(distance, c) {
  return(exp(-c * distance))
}


stimulus_f = function(){
  arms = c(0,1)
  legs = c(0,1)
  eyes = c(0,1)
  spots = c(0,1)
  color = c(0,1)
  
  stimuli = expand.grid(arms_up = arms, legs_thick = legs, eyes_stalks = eyes, spots_on = spots, color_green = color) %>% 
    mutate(danger = ifelse(spots_on == 1 & eyes_stalks == 1, 1, 0)) 
  
  perm = sample(nrow(stimuli))
  stimuli = stimuli[perm,]
  
  
}

### generative model ###
gcm <- function(w, c, obs, cat_one, quiet = TRUE) {
  # create an empty list to save probability of saying "1" for each trial
  r <- c()
  
  ntrials <- nrow(obs)
  
  for (i in 1:ntrials) {
    # If quiet is FALSE, print every ten trials
    if (!quiet && i %% 10 == 0) {
      print(paste("i =", i))
    }
    # if this is the first trial, or there any category with no exemplars seen yet, set the choice to random
    if (i == 1 || sum(cat_one[1:(i - 1)]) == 0 || sum(cat_one[1:(i - 1)]) == (i - 1)) {
      r <- c(r, .5)
    } else {
      similarities <- c()
      # for each previously seen stimulus assess distance and similarity
      for (e in 1:(i - 1)) {
        sim <- similarity(distance(obs[i, ], obs[e, ], w), c)
        similarities <- c(similarities, sim)
      }
      # Calculate prob of saying "1" by dividing similarity to 1 by the sum of similarity to 1 and to 2
      numerator <- 0.5 * sum(similarities[cat_one[1:(i - 1)] == 1])
      denominator <- 0.5 * sum(similarities[cat_one[1:(i - 1)] == 1]) + 0.5 * sum(similarities[cat_one[1:(i - 1)] == 0])
      r <- c(r, numerator / denominator)
    }
  }
  
  return(rbinom(ntrials, 1, r))
}