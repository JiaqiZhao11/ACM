
# Agent: shifting bias along trials
ShiftBiasAgent_f <- function(ntrials, biases, nbiases, i) {
  for (n in 1:nbiases){
    if ((i / (ntrials/nbiases)) > (n-1) & (i / (ntrials/nbiases)) <= n){bias = biases[n]}
  }
  
  return(bias)
  
}

# Agent: HGF_3
sigmoid = function(x) {
  1 / (1 + exp(-x))
}

hgf_agent_3level = function(ntrials,u,omega,theta,kappa){
  
  
  mu2 = array(NA,ntrials)
  sa2 = array(NA,ntrials)
  mu1hat = array(NA,ntrials)
  da = array(NA,ntrials)
  sa1hat = array(NA,ntrials)
  sa2hat = array(NA,ntrials)
  yhat= array(NA,ntrials)
  exp_p= array(NA,ntrials)
  
  mu3 = array(NA,ntrials)
  sa3 = array(NA,ntrials)
  pi3hat = array(NA,ntrials)
  pi3 = array(NA,ntrials)
  da2 = array(NA,ntrials)
  r2 = array(NA,ntrials)
  w2 = array(NA,ntrials)
  
  
  
  mu2[1] = 0
  sa2[1] = 2
  mu1hat[1] = 0.5
  sa1hat[1] = 0
  sa2hat[1] = 0
  da[1] = 0
  yhat[1] = 0
  
  mu3[1] <- 0
  sa3[1] <- 2
  pi3hat[1] <- 0
  pi3[1] <- 1/sa3[1]
  
  
  
  da2[1] <- 0
  r2[1] <- 0
  w2[1] <- 0
  
  
  for  (t in 2:ntrials){
    
    sa2hat[t] <- sa2[t-1]+exp(kappa*mu3[t-1]+omega)
    mu1hat[t] = sigmoid(mu2[t-1])
    
    
    sa1hat[t] = mu1hat[t]*(1-mu1hat[t])
    
    da[t] = u[t]-mu1hat[t]
    
    
    sa2[t] = 1/((1/sa2hat[t])+sa1hat[t])
    
    mu2[t] = mu2[t-1]+da[t]*sa2[t]
    
    
    da2[t] <- ((sa2[t]+(mu2[t]-mu2[t-1])^2)/(sa2[t-1]+exp(kappa*mu3[t-1]+omega)))-1
    
    r2[t] <-(exp(kappa*mu3[t-1]+omega)-sa2[t-1])/(sa2[t-1]+exp(kappa*mu3[t-1]+omega))
    
    w2[t] <-exp(kappa*mu3[t-1]+omega)/(sa2[t-1]+exp(kappa*mu3[t-1]+omega))
    
    pi3hat[t] <- 1/(sa3[t-1]+theta)
    
    pi3[t] <- pi3hat[t]+(kappa^2/2)*w2[t]*(w2[t]+r2[t]*da2[t])
    
    
    sa3[t] = 1/pi3[t]
    
    mu3[t] <- mu3[t-1]+sa3[t]*(kappa/2)*w2[t]*da2[t]
    
    
    
    yhat[t] = rbinom(n = 1,size = 1, prob = mu1hat[t])
    
  }
  
  data = data.frame(u,mu1hat,da,sa1hat,mu2,sa2,sa2hat,mu3,sa3,pi3hat,w2,r2,da2,yhat)
  
}



# game simulation: HGF_3 as matcher vs ShiftBiasAgent as non_matcher
hgf3_vs_shiftBias = function(ntrials, biases, nbiases, omega, theta, kappa) {
  ntrials = ntrials
  choice1 = array(NA,ntrials)
  choice2 = array(NA,ntrials)
  u = array(NA,ntrials)
  bias = array(NA, ntrials)
  
  bias[1] = biases[1]
  
  choice1[1] = rbinom(1,1,bias[1])
  
  for (i in 2:ntrials) {
    bias[i] = ShiftBiasAgent_f(ntrials = ntrials,
                            biases = biases,
                            nbiases = nbiases,
                            i = i)
    choice1[i] = rbinom(1,1,bias[i])
  }
  
    datahgf = hgf_agent_3level(ntrials = ntrials,
                                  u = choice1,
                                  omega = omega,
                                  theta = theta,
                                  kappa = kappa)
  
  
  #return(list(choice1 = choice1, choice2 = choice2))
  return(list(datahgf = datahgf, choice1 = choice1 ))
}




rw_vs_rw = function(ntrials,alpha1_l,alpha1_w,alpha2_l,alpha2_w,bias1, bias2, incentive1, incentive2){
  
  
  ntrials = ntrials
  rw1 = array(NA,ntrials)
  rw2 = array(NA,ntrials)
  expectation1 = array(NA,ntrials)
  expectation2 = array(NA,ntrials)
  
  feedback_rw1 = array(NA,ntrials)
  feedback_rw2 = array(NA,ntrials)
  
  expectation1[1] = bias1
  expectation2[1] = bias2
  
  
  rw1[1] = rbinom(1,1,expectation1[1])
  rw2[1] = rbinom(1,1,expectation2[1])
  
  if(rw1[1] == rw2[1]){
    feedback_rw1[1] = 1
    feedback_rw2[1] = 1-feedback_rw1[1]
  } else{
    feedback_rw1[1] = 0
    feedback_rw2[1] = 1-feedback_rw1[1]
  }
  
  for (i in 2:ntrials){
    
    
    expectation1[i] = rw_agent2(previous_expected = expectation1[i-1],
                                previous_them = rw2[i-1],
                                alpha_w = alpha1_w,
                                alpha_l = alpha1_l,
                                feedback = feedback_rw1[i-1])
    if(incentive1 == 0){
    rw1[i] = rbinom(1,1,expectation1[i])
    }else{
      rw1[i] = rbinom(1,1,(1-expectation1[i]))
    }
    expectation2[i] = rw_agent2(previous_expected = expectation2[i-1],
                                previous_them = rw1[i-1],
                                alpha_w = alpha2_w,
                                alpha_l = alpha2_l,
                                feedback = feedback_rw2[i-1])
    
    
    if(incentive2 == 0){
    rw2[i] = rbinom(1,1,(1-expectation2[i]))
    #this overwrites the second agents decision as he would otherwise play to match which he shouldn't he tries to not match.
    }else{
      rw2[i] = rbinom(1,1,(1-(1-expectation2[i])))
      #this overwrites the second agents decision as he would otherwise play to match which he shouldn't he tries to not match.
    }
    
    
    
    if(rw1[i] == rw2[i]){
      feedback_rw1[i] = 1
      feedback_rw2[i] = 1-feedback_rw1[i]
    } else{
      feedback_rw1[i] = 0
      feedback_rw2[i] = 1-feedback_rw1[i]
    }
    
    
  }
  
  return(list(rw1 = rw1, rw2 = rw2, feedback_rw1 = feedback_rw1,feedback_rw2 = feedback_rw2, expectation1 = expectation1, expectation2 = expectation2))
  
}


rw_agent2 = function(previous_expected, previous_them, alpha_w, alpha_l, feedback){
  if(feedback == 1){
    expected = previous_expected + alpha_w * (previous_them-previous_expected)
  }
  if(feedback == 0){
    expected = previous_expected + alpha_l * (previous_them-previous_expected)
  }
  
  return(expected)
}



#plot choices:


plot_choices = function(df){
  
  return(df %>% ggplot()+theme_classic()+geom_line(color = "red",aes(1:ntrials, df$rw1))+geom_line(color = "blue",aes(1:ntrials, df$rw2))+xlab("Trial")+ylab("Choice")+
    ggtitle("Plot of choices of matcher (red) and non-matcher (blue)"))
  
}



plot_cumwin = function(df){
  
  df = df %>% mutate(trials = 1:nrow(df)) %>%  mutate(cumrw1 = cumsum(feedback_rw1)/seq_along(feedback_rw1),
                                                   cumrw2 = cumsum(feedback_rw2)/seq_along(feedback_rw2))
  
  plot = df %>% ggplot()+
    theme_classic()+
    geom_line(color = "red",aes(trials, cumrw1))+
    geom_line(color = "blue",aes(trials, cumrw2))+
    xlab("Trial")+
    ylab("Procent of wins")+
    ggtitle("Plot of Procent of wins of matcher (red) and non-matcher (blue)")
  
  return(list(plot = plot, data = df))
  
  
}




rw_vs_rw_times = function(times, ntrials,alpha1_l,alpha1_w,alpha2_l,alpha2_w,bias1, bias2, incentive1, incentive2){
  
  agg = data.frame()
  for (i in 1:times){
    
    df = rw_vs_rw(ntrials = ntrials,
             alpha1_l = alpha1_l,
             alpha1_w = alpha1_w,
             alpha2_l = alpha2_l,
             alpha2_w = alpha2_w,
             bias1 = bias1,
             bias2 = bias2,
             incentive1 = incentive1,
             incentive2 = incentive2)
    
    df = data.frame(df)
    
    cumm = plot_cumwin(df)    
    data = cumm$data
    data$trial = 1:ntrials
    agg = rbind(agg,data)
    
    
  }
  
  agg = agg %>% 
    group_by(trials) %>% 
    summarize(n = n(),meanrw1 = mean(cumrw1),meanrw2 = mean(cumrw2), serw1 = sd(cumrw1)/sqrt(n),serw2 = sd(cumrw2)/sqrt(n))
  
  return(agg)
  
}


plot_agg = function(aggregate_data){
  plot_agg = aggregate_data %>% ggplot()+theme_classic()+
    geom_line(color = "red",aes(trials, meanrw1))+
    geom_line(color = "blue",aes(trials, meanrw2))+
    xlab("Trial")+
    ylab("Procent of wins")+
    ggtitle("Plot of Procent of wins of matcher (red) and non-matcher (blue)")+
    geom_ribbon(aes(x = trials, ymin = meanrw1-2*serw1, ymax = meanrw1+2*serw1), fill = "red", alpha = 0.2)+
    geom_ribbon(aes(x = trials, ymin = meanrw2-2*serw2, ymax = meanrw2+2*serw2), fill = "blue", alpha = 0.2)+theme(text = element_text(size = 20))
  
  return(plot_agg)
}





# MixRLAgent uses a strategy that combines two strategies together 
# on the first half trials and second half trials respectively

MixRLAgent_f <- function(previous_expectation, previous_other, feedback, alpha1, alpha2,i) {
  #probmodel = rcat(1,c(p1,p2,p3,p4)) not sure if this can be applied
  
  if (i < 60){
    expectation =  previous_expectation + alpha1 * (previous_other - previous_expectation)
  }else{
    expectation =  previous_expectation + alpha2 * (previous_other - previous_expectation)
  }
  
  return(expectation)
  
}

RL_vs_MIX = function(ntrials, alpha_w, alpha_l, alpha1, alpha2, bias_self, bias_other){
  # self refers to the matcher, other refers to the opponent
  self = rep(NA, ntrials) # matcher
  other = rep(NA, ntrials) # opponent
  expectation_self = rep(NA, ntrials)
  expectation_other = rep(NA, ntrials)
  feedback_self = rep(NA, ntrials)
  feedback_other = rep(NA, ntrials)
  
  expectation_self[1] = bias_self
  expectation_other[1] = bias_other
  self[1] = rbinom(1,1,expectation_self[1])
  other[1] = rbinom(1,1,expectation_other[1])
  
  if (self[1] == other[1]){
    feedback_self[1] = 1
    feedback_other[1] = 0} else{
      feedback_self[1] = 0
      feedback_other[1] = 1
    }
  
  for (i in 2:ntrials){
    
    expectation_self[i] = MixRLAgent_f(alpha1 = alpha1, 
                                       alpha2 = alpha2, 
                                       feedback = feedback_self[i-1],
                                       previous_expectation = expectation_self[i-1],
                                       previous_other = other[i-1],
                                       i = i)
    
    self[i] = rbinom(1,1,expectation_self[i])
    
    
    expectation_other[i] = RLAgent2_f(alpha_w = alpha_w,
                                      alpha_l = alpha_l,
                                      feedback = feedback_other[i-1],
                                      previous_expectation = expectation_other[i-1],
                                      previous_other = self[i-1])
    
    other[i] = rbinom(1,1,expectation_other[i])
    other[i] = 1 - other[i]
    
    if (self[i] == other[i]){
      feedback_self[i] = 1
      feedback_other[i] = 0} else{
        feedback_self[i] = 0
        feedback_other[i] = 1
      }
  }
  
  return(list(self = self, other = other, feedback_self = feedback_self, feedback_other = feedback_other))
}



RLAgent2_f <- function(previous_expectation, previous_other, feedback, alpha_w, alpha_l){
  
  if (feedback == 1){
    expectation =  previous_expectation + alpha_w * (previous_other - previous_expectation)
  } else {
    expectation =  previous_expectation + alpha_l * (previous_other - previous_expectation)
  }
  
  return(expectation)
}







zib = function(mean, sigma){
  kappa = 1/sigma
  
  alpha = mean*(kappa-1)
  beta = (kappa-1)*(1-mean)
  
  return(list(alpha = alpha, beta = beta))
}






rw_vs_rw_hier = function(subjects, ntrials,alpha1_l_mu,alpha1_l_sd,alpha1_w_mu,alpha1_w_sd,alpha2_l_mu,alpha2_l_sd,alpha2_w_mu,alpha2_w_sd,bias1_mu,bias1_sd, bias2_mu,bias2_sd, incentive1, incentive2){
  #making participants from a hiercial model:
  
  alpha1_l = rbeta(subjects, zib(alpha1_l_mu, alpha1_l_sd)$alpha,zib(alpha1_l_mu, alpha1_l_sd)$beta)
  alpha1_w = rbeta(subjects, zib(alpha1_w_mu, alpha1_w_sd)$alpha,zib(alpha1_w_mu, alpha1_w_sd)$beta)
  
  alpha2_l = rbeta(subjects, zib(alpha2_l_mu, alpha2_l_sd)$alpha,zib(alpha2_l_mu, alpha2_l_sd)$beta)
  alpha2_w = rbeta(subjects, zib(alpha2_w_mu, alpha2_w_sd)$alpha,zib(alpha2_w_mu, alpha2_w_sd)$beta)
  
  bias1 = rbeta(subjects, zib(bias1_mu, bias1_sd)$alpha,zib(bias1_mu, bias1_sd)$beta)
  bias2 = rbeta(subjects, zib(bias2_mu, bias2_sd)$alpha,zib(bias2_mu, bias2_sd)$beta)
  
  alpha2_l = rbeta(subjects, zib(alpha2_l_mu, alpha2_l_sd)$alpha,zib(alpha2_l_mu, alpha2_l_sd)$beta)
  alpha2_w = rbeta(subjects, zib(alpha2_w_mu, alpha2_w_sd)$alpha,zib(alpha2_w_mu, alpha2_w_sd)$beta)
  
  
  
  true = data.frame(pair = 1:subjects, alpha1_l = alpha1_l, alpha1_w = alpha1_w, alpha2_l = alpha2_l, alpha2_w = alpha2_w, bias1 = bias1, bias2 = bias2, incentive1 = incentive1, incentive2 = incentive2, bias1 = bias1, bias2 = bias2)
  
  agg = data.frame()
  rw1 = data.frame(1:ntrials)
  rw2 = data.frame(1:ntrials)
  rw1_fb = data.frame(1:ntrials)
  rw2_fb = data.frame(1:ntrials)
  
  for (i in 1:subjects){
    
    df = rw_vs_rw(ntrials = ntrials,
                  alpha1_l = alpha1_l[i],
                  alpha1_w = alpha1_w[i],
                  alpha2_l = alpha2_l[i],
                  alpha2_w = alpha2_w[i],
                  bias1 = bias1[i],
                  bias2 = bias2[i],
                  incentive1 = incentive1,
                  incentive2 = incentive2)
    
    data = data.frame(df)
    
    
    data$pair = i
    agg = rbind(agg,data)
    rw1 = cbind(rw1, data$rw1)
    rw2 = cbind(rw2, data$rw1)
    rw1_fb = cbind(rw1_fb, data$feedback_rw1)
    rw2_fb = cbind(rw2_fb, data$feedback_rw2)
    
    
  }
  
  rw1[,1] = NULL
  colnames(rw1) = 1:subjects
  rw2[,1] = NULL
  colnames(rw2) = 1:subjects
  rw1_fb[,1] = NULL
  colnames(rw1_fb) = 1:subjects
  rw2_fb[,1] = NULL
  colnames(rw2_fb) = 1:subjects
  
  
  return(list(triallevel = agg, pairlevel = true, rw1 = rw1, rw2 = rw2, rw1_fb = rw1_fb, rw2_fb = rw2_fb))
  
}








