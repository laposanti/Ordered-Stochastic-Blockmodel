
#-------------------------------------------------------------------------------
############################ MCMC functions ####################################
#-------------------------------------------------------------------------------

tuning_proposal<- function(iteration, acceptance_count, sigma, acceptanceTarget, min_sigma){
  #compute acceptance rate
  acceptanceRate <- acceptance_count / iteration
  #setting the change in the variance
  
  #passing top the log scale, to have a finer scale
  lsi= log(sigma) + (iteration**(-0.8))*(acceptanceRate - acceptanceTarget)
  
  # Update the proposal standard deviations
  sigma_updated <- min(max(min_sigma,exp(lsi)),0.3)
  return(sigma_updated)
}


label_probability <- function(distance, sigma_z) {
  prob <- dnorm(distance, mean = 0, sd = sigma_z)
  return(prob)
}

upper.tri.extractor = function(x){ x[upper.tri(x,diag = T)]}

inverse_logit_f = function(x){
  y= exp(x)/(1+exp(x))
  return(y)
}








#-------------------------------------------------------------------------------
#-------------------------- P within the model ---------------------------------
#-------------------------------------------------------------------------------

#---------------------- Computing the probabilities ----------------------------


# Likelihood

llik_over_blocks_f_binomial = function(lamdabar, ybar, mbar, theta, K, t=1){
  P = inverse_logit_f(theta)
  A_current = t*sum(lamdabar + ybar * log(P)+ (mbar)*log(1 -  P))
  
  return(A_current)
}

# Density for P

theta_prior_probability = function(theta,K,mu_vec, sigma_squared, model){
  p_theta = matrix(0,K,K)
  if(model == "WST"){
    for(diag_i in 0:(K-1)){
      p_theta[col(p_theta)-row(p_theta)==diag_i] <- dunif(x = theta[col(theta)-row(theta)==diag_i], 
                                                          min= mu_vec[diag_i+1]-sigma_squared,
                                                          max = mu_vec[diag_i+2]+sigma_squared)
    }
  }else if (model == "Simple"){
    P = inverse_logit_f(theta)
    P_upper.tri =upper.tri.extractor(P)
    
    p_theta[upper.tri(p_theta,diag=T)] = dbeta(P_upper.tri,1,1)*(exp(P_upper.tri)/((1+exp(P_upper.tri))**2))
  }else if(model =='SST'){
    
    for(diag_i in 0:(K-1)){
      p_theta[col(p_theta)-row(p_theta)==diag_i] <- dunif(x = theta[col(theta)-row(theta)==diag_i], min= 
                                                            mu_vec[diag_i+1],
                                                          max = mu_vec[diag_i+2])
    }
  }
  log_p_theta = sum(log(p_theta[upper.tri(p_theta,diag = T)]))
  
  return(log_p_theta)
}

#mu prior distribution
d_sA_mu = function(K,mu_vec){
  # mu_0_lp = dnorm(mu_vec[1],0,1,log=T)
  fac <- lfactorial(K+1)
  joint_density<- sum(log(dtruncnorm(mu_vec[1:(K+1)],a = 0,mean = 0,sd = 3)) - 
                        log(1- pnorm(0,mean = 0,sd = 3)))
  return(joint_density+fac)
}

order_stat_truncnorm = function(K, mu, mean, sd, lb, ub){
  fac <- lfactorial(K+1)
  joint_density<- sum(log(dtruncnorm(mu,a = lb,b = ub, mean = mean,sd = sd)) - 
                        log(pnorm(ub,mean = mean,sd = sd) - pnorm(lb,mean = mean,sd = sd)))
  return(fac+joint_density)
}




# Proportional posterior


lprop_posterior <- function(lamdabar, ybar,mbar,theta, 
                            alpha_vec, n_k,sigma_squared, mu_vec,K, model,t){
  
  if(model=='Simple'){
    
    #log likelihood
    log_lik <- llik_over_blocks_f_binomial(lamdabar = lamdabar,ybar =  ybar,  mbar = mbar,theta = theta, K=K, t=t)
    #log prior on z
    prior_z <- ddirichlet_multinomial(N = sum(n_k),K = K,n_k = n_k,my_alpha =  alpha_vec)
    #log prior on P
    prior_theta<- theta_prior_probability(theta = theta,K=K,mu_vec = mu_vec, sigma_squared = sigma_squared, model=model)
    #computing the whole log proportional posterior
    results<- log_lik+ prior_theta + prior_z 
    
  }else if (model=='WST'){
    
    #log likelihood
    log_lik <- llik_over_blocks_f_binomial(lamdabar = lamdabar,ybar =  ybar,  mbar = mbar,theta = theta,K=K,t=t)
    #log prior on z
    prior_z <- ddirichlet_multinomial(N = sum(n_k), K = K,n_k = n_k, my_alpha = alpha_vec)
    #log prior on P
    prior_theta<- theta_prior_probability(theta = theta, K=K,mu_vec = mu_vec,  sigma_squared = sigma_squared,model = model)
    #log prior on mu
    hyperprior_mu <- d_sA_mu(K = K, mu_vec = mu_vec)
    #log prior on sigma^2
    hyperprior_sigmasquared <- LaplacesDemon::dinvgamma(x = sigma_squared,shape = 0.001,scale = 0.001,log = T)
    #computing the whole log proportional posterior
    results<- log_lik+ prior_theta + prior_z + hyperprior_mu + hyperprior_sigmasquared 
    
  }else if(model =='SST'){
    
    #log likelihood
    log_lik <- llik_over_blocks_f_binomial(lamdabar = lamdabar,ybar =  ybar,  mbar = mbar,theta = theta,K=K,t=t)
    #log prior on P
    prior_theta<- theta_prior_probability(theta = theta, K = K,mu_vec = mu_vec, sigma_squared = sigma_squared,model = model)
    #log prior on z
    prior_z <- ddirichlet_multinomial(N = sum(n_k), K = K,n_k = n_k, my_alpha = alpha_vec)
    #log prior on mu
    hyperprior_mu <- d_sA_mu(K = K,mu_vec =  mu_vec)
    #computing the whole log proportional posterior
    results <- log_lik + prior_z + hyperprior_mu  + prior_theta
    
  }
  return(results)
}

#-------------------------- MCMC steps -----------------------------------------

theta_update_f = function(lamdabar,ybar,mbar,theta, alpha_vec, n_k,
                          sigma_squared, mu_vec,K, tau_theta,
                          acc.count_theta, model,t){
  
  theta_current<- theta
  #Updating each entry of P, one at the time
  ut <- upper.tri(theta,diag = T) # get the logical matrix for upper triangular elements
  theta_combn = which(ut, arr.ind = TRUE) # get the indices of the upper triangular elements
  uo<- data.frame(theta_combn[sample(nrow(theta_combn)), ])# permuting the order of the rows
  
  for(i_th in 1:nrow(uo)){
    theta_prime<-theta_current
    i_star<- uo$row[i_th]
    j_star<- uo$col[i_th]
    if(model=="SST"){
      lower.bound = mu_vec[j_star - i_star + 1]
      upper.bound = mu_vec[j_star - i_star + 2]
    }else if(model == 'WST'){
      
      lower.bound = mu_vec[j_star - i_star + 1] - sigma_squared
      upper.bound = mu_vec[j_star - i_star + 2] + sigma_squared
    }else if(model=='Simple'){
      lower.bound = -5
      upper.bound = +10
    }
    
    theta_prime[i_star,j_star]<- rtruncnorm(1, a = lower.bound , b = upper.bound,
                                            mean = theta_current[i_star,j_star], 
                                            sd =  .25)
    
    theta_prime[lower.tri(theta_prime)] = - t(theta_prime)[lower.tri(theta_prime)]
    
    #computing the proportional posterior in P'
    prop_posterior_prime <-  lprop_posterior(lamdabar = lamdabar, ybar = ybar, mbar = mbar, 
                                             theta = theta_prime,
                                             K = K,mu_vec = mu_vec,
                                             sigma_squared = sigma_squared,
                                             alpha_vec = alpha_vec,n_k = n_k,model = model,t = t)
    
    #evaluating the proportional posterior in P^(t)
    prop_posterior_current<- lprop_posterior(lamdabar = lamdabar, ybar = ybar, 
                                             mbar = mbar, theta = theta_current,
                                             K = K,mu_vec = mu_vec,
                                             sigma_squared = sigma_squared,
                                             alpha_vec = alpha_vec,n_k = n_k,model = model,t = t)
    
    #evaluating the proposal density g(P'| P^(t)) 
    log_proposal_prime <- log(dtruncnorm(theta_prime[i_star,j_star],
                                         mean = theta_current[i_star,j_star],
                                         sd = tau_theta[i_star,j_star], 
                                         a =   lower.bound, 
                                         b  = upper.bound))
    
    #evaluating the proposal density g(P^(t)| P') 
    log_proposal_current <- log(dtruncnorm(theta_current[i_star,j_star],
                                           mean = theta_prime[i_star,j_star],
                                           sd = tau_theta[i_star,j_star], 
                                           a =  lower.bound,
                                           b  = upper.bound))
    
    #acceptance ratio
    log_r=  prop_posterior_prime + log_proposal_current - 
      prop_posterior_current - log_proposal_prime
    
    
    #create statements that check conditiond to accept move
    MH_condition_P_update= min(log_r,0)>=log(runif(1))
    if(MH_condition_P_update){
      acc.count_theta[i_star,j_star] =acc.count_theta[i_star,j_star] +1
      theta_current = theta_prime
    }
    
  }
  
  theta_current= theta_prime
  
  return(list(acc.moves = acc.count_theta,
              theta= theta_current))
  
} 





mu_update_f = function(lamdabar, ybar,mbar,theta, alpha_vec, n_k,
                       sigma_squared, mu_vec,tau_mu_vec, K,
                       acc.count_mu_vec,model,t){
  
  #computing the proportional posterior in mu' ~ g(mu^(t), tau_mu_vec)
  #mu_1_K_prime <- truncnorm::rtruncnorm(K+1,a = 0,b = 10, mean = mu_vec[1:(K+1)],sd = .05)
  # mu_0_prime = truncnorm::rtruncnorm(1,a = -Inf,b = min(mu_1_K_prime), mean = mu_vec[1],sd = .2)
  # mu_vec_prime <- c(mu_0_prime,sort(mu_1_K_prime))
  
  
  split = diag_split_matrix(theta)
  mins = vector()
  maxs = vector()
  for(splittino in 1:length(split)){
    mins <- append(mins,max(split[[splittino]]))
    maxs <- append(maxs,min(split[[splittino]]))
  }
  lbs = c(0, mins)
  ubs = c(maxs, 10)
  mu_vec[1] =0 
  for(mu in 2:length(mu_vec)){
    mu_1_K_prime <- rtruncnorm(1,a = lbs[mu],b = ubs[mu],mean = mu_vec[mu],sd = 1)

    mu_vec_prime= mu_vec
    mu_vec_prime[mu] <- mu_1_K_prime
    
    #computing the proportional posterior in mu'
    prop_posterior_prime <- lprop_posterior(lamdabar = lamdabar,ybar = ybar,mbar = mbar,theta = theta, 
                                            alpha_vec = alpha_vec, 
                                            n_k = n_k,sigma_squared = sigma_squared
                                            ,mu_vec = mu_vec_prime,K = K,model=model,t = t)
    
    #evaluating the proportional posterior in mu^(t)
    prop_posterior_current <- lprop_posterior(lamdabar = lamdabar,ybar = ybar,mbar = mbar,
                                              theta = theta,alpha_vec = alpha_vec,
                                              n_k = n_k,sigma_squared = sigma_squared,
                                              mu_vec = mu_vec,K = K,model = model,t = t)
    
    #evaluating the proposal density g(mu'| mu^(t)) 
    p_proposal_prime = dtruncnorm(mu_1_K_prime,a = lbs[mu],b = ubs[mu], mean = mu_vec[mu],sd = 1)
    p_proposal_current = dtruncnorm(mu_vec[mu],a = lbs[mu],b = ubs[mu], mean = mu_1_K_prime,sd = 1)
    #evaluating the proposal density g(mu^(t)| mu') 
    # log_proposal_mu_1_K_current <-order_stat_truncnorm(K, mu = mu_vec[1:(K+1)], mean = mu_1_K_prime,
    #                                                    sd =tau_mu_vec,lb = 0,ub = 10)
    # log_proposal_mu0_current <- log(truncnorm::dtruncnorm(mu_vec[1] ,a= -Inf,b = min(mu_1_K_prime), 
    #                                                       mean = mu_0_prime,sd = tau_mu_vec))
    # 
    # log_proposal_mu_current = log_proposal_mu_1_K_current
    
    
    log_r =  prop_posterior_prime  - prop_posterior_current + log(p_proposal_current) - log(p_proposal_prime)
    
    #create statements that check conditionion to accept move
    MH_condition_mu_vec= min(log_r,0)>=log(runif(1))
    if(MH_condition_mu_vec){
      acc.count_mu_vec[mu] = acc.count_mu_vec[mu] +1
      mu_vec <- mu_vec_prime
    }
  }
  return(list(acc.moves = acc.count_mu_vec,
              mu_vec=mu_vec))
  
} 





sigma_squared_update_f= function(lamdabar,ybar,mbar,theta, alpha_vec, n_k,
                                 sigma_squared, mu_vec,K, tau_sigma_squared,
                                 acc.count_sigma_squared, model,t){
  
  
  #simulating (sigma^2)' from a g ~ truncated normal
  sigma_squared_prime <- rtruncnorm(1,a = 0,b = 1 ,
                                    mean = sigma_squared,
                                    sd = tau_sigma_squared)
  
  #computing the proportional posterior in (sigma^2)'
  prop_posterior_prime <- lprop_posterior(lamdabar, ybar,mbar,theta, 
                                          alpha_vec, n_k,sigma_squared_prime, mu_vec,K,model = model,t = t)
  #evaluating the proportional posterior in (sigma^2)^(t)
  prop_posterior_current<- lprop_posterior(lamdabar, ybar,mbar,theta,
                                           alpha_vec, n_k,sigma_squared, mu_vec,K,model=model,t = t)
  
  #evaluating the proposal density g(sigma^2)') 
  log_proposal_prime <- log(dtruncnorm(sigma_squared_prime,a = 0,b = 1,
                                       mean = sigma_squared,
                                       sd = tau_sigma_squared))
  
  #evaluating the proposal density g(sigma^2)^(t)) 
  log_proposal_current <- log(dtruncnorm(sigma_squared,a = 0,b = 1,
                                         mean = sigma_squared_prime,
                                         sd = tau_sigma_squared))
  #acceptance ratio
  log_r=  prop_posterior_prime + log_proposal_current - 
    prop_posterior_current - log_proposal_prime
  
  #create statements that check conditionion to accept move
  MH_condition_sigma_squared= min(log_r,0)>=log(runif(1))
  if(MH_condition_sigma_squared){
    acc.count_sigma_squared = acc.count_sigma_squared+1
    sigma_squared <- sigma_squared_prime
  }
  
  
  return(list(acc.moves = acc.count_sigma_squared,
              sigma_squared=sigma_squared))
  
} 


z_update_f = function(N_ij, Y_ij, z,lamdabar,ybar,mbar, theta, alpha_vec, n_k,
                      K, 
                      acc.count_z,labels_available,model,t){
  
  P = inverse_logit_f(theta)
  #redefining prior quantities
  n_prime = n_k
  
  B_prime<- ddirichlet_multinomial(N = n,K = K,n_k = n_k, my_alpha =  alpha_vec)
  
  
  upper.tri.Y_ij = Y_ij*upper.tri(Y_ij)
  upper.tri.N_ij = N_ij*upper.tri(N_ij)
  coef1 = lchoose(N_ij, Y_ij)*upper.tri(N_ij)
  n_minus_y1 <- (N_ij-Y_ij)*upper.tri(N_ij)
  
  z_prime = z
  z_mat_prime = vec2mat_0_P(z_prime, theta)
  
  A_prime = llik_over_blocks_f_binomial(lamdabar = lamdabar,
                                        ybar = ybar,
                                        mbar = mbar,
                                        theta = theta,
                                        K = K,
                                        t = t)
  
  scanning_order = sample(1:n,n, replace=F)
  # full sweep
  for(i_th_turn in scanning_order){
    #save current label of z_ii
    k_prime <- z_prime[i_th_turn]
    
    # Sample a new label using the adjusted probabilities
    if(model != 'Simple'){
      labels_to_sample = c(min(k_prime+1, K), max(k_prime-1, 1))
      
      #checking which adjacent labels are available, and if k_prime is among those, exclude it
      avail_options_k_prime = setdiff(labels_to_sample, k_prime)
      
      #sample one from those available
      k_scanning <- sample(x = avail_options_k_prime, size = 1, replace = F)
      #the probability of extracting k_scanning is 1/number of labels available (either 1 or two)
      probability_k_prime_given_k_current = 1/length(avail_options_k_prime)
      #given k_scanning, what is the porbability of getting k_current?
      #1) get the adjacent labels to k_scanning
      labels_to_sample_current = c(min(k_scanning+1, K), max(k_scanning-1, 1))
      #2) compute how many they are
      avail_options_k_current = setdiff(labels_to_sample_current, k_scanning)
      #3) let's compute the conditional probability to go from k_scanning to k_current
      probability_k_current_given_k_prime = 1/length(avail_options_k_current)
    }else{
      #checking which adjacent labels are available, and if k_prime is among those, exclude it
      avail_options_k_prime = setdiff(labels_available, k_prime)
      
      #sample one from those available
      k_scanning <- sample(x = avail_options_k_prime, size = 1, replace = F)
      
      probability_k_prime_given_k_current = 1/(K-1)
      probability_k_current_given_k_prime = 1/(K-1)
    }
    
    
    #assign the new label
    z_mat_scanning = z_mat_prime
    z_mat_scanning[i_th_turn,] <- rep(0,K) + (labels_available == k_scanning)
    
    
    
    
    #----------
    # Subtracting the contribution of the current i-th item from ybar, mbar and lambdabar
    
    #-------ybar
    i_th_victoriesvs_all_clusters = t(z_mat_prime)%*%(upper.tri.Y_ij)[i_th_turn,]
    all_clusters_victories_vs_ith = (upper.tri.Y_ij)[,i_th_turn]%*%z_mat_prime
    
    subtract_matrix = matrix(0,K,K)
    subtract_matrix[k_prime,]= i_th_victoriesvs_all_clusters
    subtract_matrix[,k_prime]= subtract_matrix[,k_prime] + all_clusters_victories_vs_ith
    
    #-------mbar
    i_th_games_vs_all_clusters = t(z_mat_prime)%*%(n_minus_y1)[i_th_turn,]
    all_clusters_games_vs_ith = (n_minus_y1)[,i_th_turn]%*%z_mat_prime
    
    
    subtract_matrix_N_ij = matrix(0,K,K)
    subtract_matrix_N_ij[k_prime,]= i_th_games_vs_all_clusters
    subtract_matrix_N_ij[,k_prime] = subtract_matrix_N_ij[,k_prime] + all_clusters_games_vs_ith
    
    #-------lamdabar
    i_th_possibilities_vs_all_clusters = t(z_mat_prime)%*%(coef1)[i_th_turn,]
    all_clusters_possibilities_vs_ith = (coef1)[,i_th_turn]%*%z_mat_prime
    
    subtract_matrix_lamda_ij = matrix(0,K,K)
    subtract_matrix_lamda_ij[k_prime,]= i_th_possibilities_vs_all_clusters
    subtract_matrix_lamda_ij[,k_prime] = subtract_matrix_lamda_ij[,k_prime] + all_clusters_possibilities_vs_ith
    
    #subtracting contribution to the likelihood
    A_minus = sum(subtract_matrix_lamda_ij+ subtract_matrix * log(P)+(subtract_matrix_N_ij)*log(1 - P))
    
    #----------
    # Adding the contribution of the new i-th item from ybar, mbar and lambdabar
    
    #-------ybar
    i_th_victoriesvs_all_clusters_scanning = t(z_mat_scanning)%*%(upper.tri.Y_ij)[i_th_turn,]
    all_clusters_victories_vs_ith_scanning = (upper.tri.Y_ij)[,i_th_turn]%*%z_mat_scanning
    
    
    to_add_matrix = matrix(0,K,K)
    to_add_matrix[k_scanning,] = i_th_victoriesvs_all_clusters_scanning
    to_add_matrix[,k_scanning] = to_add_matrix[,k_scanning] + all_clusters_victories_vs_ith_scanning
    
    
    #-------mbar
    i_th_gamesvs_all_clusters_scanning = t(z_mat_scanning)%*%(n_minus_y1*upper.tri(n_minus_y1))[i_th_turn,]
    all_clusters_games_vs_ith_scanning = (n_minus_y1*upper.tri(n_minus_y1))[,i_th_turn]%*%z_mat_scanning
    
    to_add_matrix_N_ij = matrix(0,K,K)
    to_add_matrix_N_ij[k_scanning,] = i_th_gamesvs_all_clusters_scanning
    to_add_matrix_N_ij[,k_scanning] = to_add_matrix_N_ij[,k_scanning] + all_clusters_games_vs_ith_scanning
    
    #-------lambdabar
    i_th_possibilities_all_clusters_scanning = t(z_mat_scanning)%*%coef1[i_th_turn,]
    all_clusters_possibilities_vs_ith_scanning = coef1[,i_th_turn]%*%z_mat_scanning
    
    to_add_matrix_lamda_ij = matrix(0,K,K)
    to_add_matrix_lamda_ij[k_scanning,] = i_th_possibilities_all_clusters_scanning
    to_add_matrix_lamda_ij[,k_scanning] = to_add_matrix_lamda_ij[,k_scanning] + all_clusters_possibilities_vs_ith_scanning
    
    #adding contribution to the likelihood
    A_plus = sum(to_add_matrix_lamda_ij+ to_add_matrix * log(P)+(to_add_matrix_N_ij)*log(1 - P))
    
    A_scanning = A_prime - A_minus +A_plus
    
    
    n_scanning<- n_prime
    n_scanning[c(k_prime, k_scanning)] <- n_prime[c(k_prime, k_scanning)] + c(-1, 1)
    
    B_scanning<- ddirichlet_multinomial(N = n,K = K,n_k = n_scanning,my_alpha = alpha_vec)
    
    log_r= t*A_scanning - t*A_prime + B_scanning - B_prime + 
      log(probability_k_current_given_k_prime) - log(probability_k_prime_given_k_current)
    
    #create statements that check conditiond to accept move
    GS_condition= min(log_r,0)>=log(runif(1))
    if(GS_condition){
      acc.count_z[i_th_turn]= acc.count_z[i_th_turn]+1
      z_mat_prime <- z_mat_scanning
      A_prime<- A_scanning
      B_prime<- B_scanning
      n_prime <- n_scanning
    }
  }
  
  z = z_mat_prime %*% matrix(labels_available,K,1)
  
  # number of victories between block p and block q
  ybar = t(z_mat_prime)%*%(Y_ij*upper.tri(Y_ij))%*%z_mat_prime
  # number of missed victories between block p and block q
  mbar<- t(z_mat_prime)%*%n_minus_y1%*%z_mat_prime
  
  lamdabar <- t(z_mat_prime)%*%(coef1)%*%z_mat_prime
  
  return(list(acc.moves = acc.count_z, z_current= z, ybar =ybar, mbar=mbar, lamdabar= lamdabar))
} 





