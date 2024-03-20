
#-------------------------------------------------------------------------------
############################ MCMC functions ####################################
#-------------------------------------------------------------------------------

tuning_proposal<- function(iteration, acceptance_count, sigma, acceptanceTarget, min_sigma){
  #compute acceptance rate
  acceptanceRate <- acceptance_count / iteration
  #setting the change in the variance

  #passing top the log scale, to have a finer scale
  lsi= log(sigma) + (iteration**(-1/2))*(acceptanceRate - acceptanceTarget)

  # Update the proposal standard deviations
  sigma_updated <- max(min_sigma,exp(lsi))
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
#-------------------------- P integrated out -----------------------------------
#-------------------------------------------------------------------------------


#---------------------- Computing the probabilities ----------------------------

log_lik_f_binom = function(N,Y,z,P, directed=T,t=1){
  z_P<- vec2mat_0_P(z,P)
  P_nbyn<- calculate_victory_probabilities(z_P, P)
  if(directed==T){
    #computing the pairwise log-probabilitiees
    bigM = lchoose(N,Y)+(Y* log(P_nbyn)+(N-Y)*log(1 - P_nbyn))
    #remember to subtract the diagonal
    log_lik= t* sum(bigM) - sum(diag(bigM))
  }else if(directed==F){
    bigM = lchoose(N,Y)+(Y* log(P_nbyn)+(N-Y)*log(1 - P_nbyn))
    #remember to subtract the diagonal
    log_lik= t*sum(bigM*upper.tri(bigM))
  }
  return(log_lik)
}


llike_integrated_out = function(lamdabar, ybar,mbar,K,U_vec, sigma_squared){
  
  beta_params = beta_mean_var(U_vec,rep(sigma_squared,K) )
  alpha_k = beta_params$alpha
  beta_k = beta_params$beta
  
  llik= matrix(0,K,K)
  for(p in 1:K){
    for(q in p:K){
      if(q-p != 0){
        llik[p,q]<-  log(lamdabar[p,q]) + log(lamdabar[q,p]) +
          lbeta(alpha_k[q-p] + ybar[p,q]+ mbar[q,p],beta_k[q-p] + 
                  ybar[q,p]+ mbar[q,p]) - lbeta(alpha_k[q-p], beta_k[q-p])  
      }else if (p==q){
        llik[p,p]<-  log(lamdabar[p,p]) +lbeta(1 + ybar[p,p],1 
                                               + ybar[p,p]) - lbeta(1, 1)  
      }
    }
  }
  
  return(sum(llik))
}











#-------------------------------------------------------------------------------
#-------------------------- P within the model ---------------------------------
#-------------------------------------------------------------------------------

#---------------------- Computing the probabilities ----------------------------


# Likelihood

llik_over_blocks_f_binomial = function(lamdabar, ybar, mbar, P, K, t=1){
  #llik<- sum(bin_coef+ ybar*log(P) + (mbar)*log(1-P))
  llik = matrix(0,K,K)
  for(diag_iii in 0:(K-1)){
    llik[col(P)-row(P)==diag_iii] <- lamdabar[col(P)-row(P)==diag_iii]+ ybar[col(P)-row(P)==diag_iii]*P[col(P)-row(P)==diag_iii] -
      (mbar[col(P)-row(P)==diag_iii]+ybar[col(P)-row(P)==diag_iii])*log(1+exp(P[col(P)-row(P)==diag_iii]))
  }
  
  return(t*sum(llik))
}

# Density for P

P_prior_probability = function(P,K,mu_vec, sigma_squared, model){
  p_P = matrix(0,K,K)
  if(model == "WST"){
    sigma= sqrt(sigma_squared)
    for(diag_i in 0:(K-1)){
      p_P[col(p_P)-row(p_P)==diag_i] <- dunif(x = P[col(P)-row(P)==diag_i], 
                                              min= mu_vec[diag_i+1]-sigma,
                                              max = mu_vec[diag_i+2]+sigma)
    }
  }else if (model == "Simple"){
    theta = inverse_logit_f(P)
    theta=upper.tri.extractor(theta)
    
    p_P[upper.tri(p_P,diag=T)] = dbeta(theta,1,1)*(exp(theta)/((1+exp(theta))**2))
  }else if(model =='SST'){
    
    for(diag_i in 0:(K-1)){
      p_P[col(p_P)-row(p_P)==diag_i] <- dunif(x = P[col(P)-row(P)==diag_i], min= 
                                                mu_vec[diag_i+1],
                                              max = mu_vec[diag_i+2])
    }
  }
  log_p_P = sum(log(p_P[upper.tri(p_P,diag = T)]))
  
  return(log_p_P)
}

#mu prior distribution
d_sA_mu = function(K,mu_vec){
  mu_0_lp = dnorm(mu_vec[1],0,1,log=T)
  fac <- lfactorial(K)
  joint_density<- sum(log(dtruncnorm(mu_vec[2:(K+1)],a = 0,mean = 0,sd = 1)) - log(1- pnorm(0,mean = 0,sd = 1)))
  return(joint_density+fac+mu_0_lp)
}

order_stat_truncnorm = function(K, mu, mean, sd, lb, ub){
  fac <- lfactorial(K)
  joint_density<- sum(log(dtruncnorm(mu,a = lb,b = ub, mean = mean,sd = sd)) - 
                        log(pnorm(ub,mean = mean,sd = sd)- pnorm(lb,mean = mean,sd = sd)))
  return(fac+joint_density)
}




# Proportional posterior


lprop_posterior_withP <- function(lamdabar, ybar,mbar,P, 
                                  alpha_vec, n_k,sigma_squared, mu_vec,K, model,t){
  
  if(model=='Simple'){
    
    #log likelihood
    log_lik <- llik_over_blocks_f_binomial(lamdabar = lamdabar,ybar =  ybar,  mbar = mbar,P =  P, K=K, t=t)
    #log prior on z
    prior_z <- ddirichlet_multinomial(N = sum(n_k),K = K,n_k = n_k,my_alpha =  alpha_vec)
    #log prior on P
    prior_P<- P_prior_probability(P = P,K=K,mu_vec = mu_vec, sigma_squared = sigma, model=model)
    #computing the whole log proportional posterior
    results<- log_lik+ prior_P + prior_z 
    
  }else if (model=='WST'){
    
    #log likelihood
    log_lik <- llik_over_blocks_f_binomial(lamdabar = lamdabar,ybar =  ybar,  mbar = mbar,P =  P,K=K,t=t)
    #log prior on z
    prior_z <- ddirichlet_multinomial(N = sum(n_k), K = K,n_k = n_k, my_alpha = alpha_vec)
    #log prior on P
    prior_P<- P_prior_probability(P = P, K=K,mu_vec = mu_vec, sigma_squared = sigma_squared,model = model)
    #log prior on mu
    hyperprior_mu <- d_sA_mu(K = K, mu_vec = mu_vec)
    #log prior on sigma^2
    hyperprior_sigmasquared <- LaplacesDemon::dinvgamma(x = sigma_squared,shape = 0.001,scale = 0.001,log = T)
    #computing the whole log proportional posterior
    results<- log_lik+ prior_P + prior_z + hyperprior_mu + hyperprior_sigmasquared 
    
  }else if(model =='SST'){
    
    #log likelihood
    log_lik <- llik_over_blocks_f_binomial(lamdabar = lamdabar,ybar =  ybar,  mbar = mbar,P =  P,K=K,t=t)
    #log prior on P
    prior_P<- P_prior_probability(P = P, K = K,mu_vec = mu_vec, sigma_squared = sigma,model = model)
    #log prior on z
    prior_z <- ddirichlet_multinomial(N = sum(n_k), K = K,n_k = n_k, my_alpha = alpha_vec)
    #log prior on mu
    hyperprior_mu <- d_sA_mu(K = K,mu_vec =  mu_vec)
    #computing the whole log proportional posterior
    results <- log_lik + prior_z + hyperprior_mu  + prior_P
    
  }
  return(results)
}

#-------------------------- MCMC steps -----------------------------------------

P_update_f = function(lamdabar,ybar,mbar,P, alpha_vec, n_k,
                      sigma_squared, mu_vec,K, tau_P,
                      acc.count_P, model,t){
  
  P_current<- P
  #Updating each entry of P, one at the time
  ut <- upper.tri(P,diag = T) # get the logical matrix for upper triangular elements
  Pcombn = which(ut, arr.ind = TRUE) # get the indices of the upper triangular elements
  uo<- data.frame(Pcombn[sample(nrow(Pcombn)), ])# permuting the order of the rows
  n_P_entries<- nrow(uo)
  
  for(i_th in 1:n_P_entries){
    P_prime<-P_current
    i_star<- uo$row[i_th]
    j_star<- uo$col[i_th]
    if(model=="SST"){
      lower.bound = mu_vec[j_star - i_star + 1]
      upper.bound = mu_vec[j_star - i_star + 2]
    }else if(model == 'WST'){
      sigma=sqrt(sigma_squared)
      lower.bound = mu_vec[j_star - i_star + 1] - sigma
      upper.bound = mu_vec[j_star - i_star + 2] + sigma
    }else if(model=='Simple'){
      lower.bound = -5
      upper.bound = +10
    }
    
    P_prime[i_star,j_star]<- rtruncnorm(1, a = lower.bound , b = upper.bound,
                                        mean = P_current[i_star,j_star], 
                                        sd =  tau_P[i_star,j_star])
    
    P_prime[lower.tri(P_prime)] = - t(P_prime)[lower.tri(P_prime)]
    
    #computing the proportional posterior in P'
    prop_posterior_prime <-  lprop_posterior_withP(lamdabar = lamdabar, ybar = ybar, mbar = mbar, 
                                                   P = P_prime,
                                                   K = K,mu_vec = mu_vec,
                                                   sigma_squared = sigma_squared,
                                                   alpha_vec = alpha_vec,n_k = n_k,model,t)
    
    #evaluating the proportional posterior in P^(t)
    prop_posterior_current<- lprop_posterior_withP(lamdabar = lamdabar, ybar = ybar, 
                                                   mbar = mbar, P = P_current,
                                                   K = K,mu_vec = mu_vec,
                                                   sigma_squared = sigma_squared,
                                                   alpha_vec = alpha_vec,n_k = n_k,model,t)
    
    #evaluating the proposal density g(P'| P^(t)) 
    log_proposal_prime <- log(dtruncnorm(P_prime[i_star,j_star],
                                         mean = P_current[i_star,j_star],
                                         sd = tau_P[i_star,j_star], 
                                         a =   lower.bound, 
                                         b  = upper.bound))
    
    #evaluating the proposal density g(P^(t)| P') 
    log_proposal_current <- log(dtruncnorm(P_current[i_star,j_star],
                                           mean = P_prime[i_star,j_star],
                                           sd = tau_P[i_star,j_star], 
                                           a =  lower.bound,
                                           b  = upper.bound))
    
    #acceptance ratio
    log_r=  prop_posterior_prime + log_proposal_current - 
      prop_posterior_current - log_proposal_prime
    
    
    #create statements that check conditiond to accept move
    MH_condition_P_update= min(log_r,0)>=log(runif(1))
    if(MH_condition_P_update){
      acc.count_P[i_star,j_star] =acc.count_P[i_star,j_star] +1
      P_current = P_prime
    }
    
  }
  
  P_current= P_prime
  
  return(list(acc.moves = acc.count_P,
              P= P_current))
  
} 





mu_update_f_withP = function(lamdabar, ybar,mbar,P, alpha_vec, n_k,
                             sigma_squared, mu_vec,tau_mu_vec, K,
                             acc.count_mu_vec,model,t){
  
  #computing the proportional posterior in mu' ~ g(mu^(t), tau_mu_vec)
  mu_1_K_prime <- truncnorm::rtruncnorm(K,a = 0,b = 10, mean = mu_vec[2:(K+1)],sd = tau_mu_vec)
  mu_0_prime = truncnorm::rtruncnorm(1,a = -Inf,b = min(mu_1_K_prime), mean = mu_vec[1],sd = tau_mu_vec)
  mu_vec_prime <- c(mu_0_prime,sort(mu_1_K_prime))
  
  
  #computing the proportional posterior in mu'
  prop_posterior_prime <- lprop_posterior_withP(lamdabar = lamdabar,ybar = ybar,mbar = mbar,P = P, 
                                                alpha_vec = alpha_vec, 
                                                n_k = n_k,sigma_squared = sigma_squared
                                                ,mu_vec = mu_vec_prime,K = K,model=model,t = t)
  
  #evaluating the proportional posterior in mu^(t)
  prop_posterior_current <- lprop_posterior_withP(lamdabar = lamdabar,ybar = ybar,mbar = mbar,
                                                  P = P,alpha_vec = alpha_vec,
                                                  n_k = n_k,sigma_squared = sigma_squared,
                                                  mu_vec = mu_vec,K = K,model = model,t = t)
  
  #evaluating the proposal density g(mu'| mu^(t)) 
  log_proposal_mu_1_K_prime <- order_stat_truncnorm(K, mu = mu_1_K_prime, 
                                                    mean = mu_vec[2:(K+1)],sd =tau_mu_vec,lb = 0,ub = 10)
  log_proposal_mu0_prime <- log(truncnorm::dtruncnorm(mu_0_prime ,a= -Inf,b = min(mu_1_K_prime), 
                                                      mean = mu_vec[1],sd = tau_mu_vec))
  log_proposal_mu_prime = log_proposal_mu_1_K_prime+ log_proposal_mu0_prime
  
  #evaluating the proposal density g(mu^(t)| mu') 
  log_proposal_mu_1_K_current <-order_stat_truncnorm(K, mu = mu_vec[2:(K+1)], mean = mu_1_K_prime,sd =tau_mu_vec,lb = 0,ub = 10)
  log_proposal_mu0_current <- log(truncnorm::dtruncnorm(mu_vec[1] ,a= -Inf,b = min(mu_1_K_prime), 
                                                        mean = mu_0_prime,sd = tau_mu_vec))
  log_proposal_mu_current = log_proposal_mu_1_K_current+ log_proposal_mu0_current
  
  
  log_r =  prop_posterior_prime  - prop_posterior_current + 
    log_proposal_mu_current - log_proposal_mu_prime
  
  #create statements that check conditionion to accept move
  MH_condition_mu_vec= min(log_r,0)>=log(runif(1))
  if(MH_condition_mu_vec){
    acc.count_mu_vec = acc.count_mu_vec+1
    mu_vec <- mu_vec_prime
  }
  
  return(list(acc.moves = acc.count_mu_vec,
              mu_vec=mu_vec))
  
} 





sigma_squared_update_f_withP = function(lamdabar,ybar,mbar,P, alpha_vec, n_k,
                                        sigma_squared, mu_vec,K, tau_sigma_squared,
                                        acc.count_sigma_squared, model,t){
  
  
  #simulating (sigma^2)' from a g ~ truncated normal
  sigma_squared_prime <- rtruncnorm(1,a = 0,b = 1 ,
                                    mean = sigma_squared,
                                    sd = tau_sigma_squared)
  
  #computing the proportional posterior in (sigma^2)'
  prop_posterior_prime <- lprop_posterior_withP(lamdabar, ybar,mbar,P, 
                                                alpha_vec, n_k,sigma_squared_prime, mu_vec,K,model = model,t = t)
  #evaluating the proportional posterior in (sigma^2)^(t)
  prop_posterior_current<- lprop_posterior_withP(lamdabar, ybar,mbar,P,
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


z_update_f_withP = function(N_ij, Y_ij, z,lamdabar,ybar,mbar, P, alpha_vec, n_k,
                            K, 
                            acc.count_z,labels_available,model,t){
  P<- inverse_logit_f(P)
  n<- nrow(N_ij)
  A_prime<- log_lik_f_binom(N = N_ij,Y = Y_ij,z =z,P = P,directed = T)
  B_prime<- ddirichlet_multinomial(N = n,K = K,n_k = n_k, my_alpha =  alpha_vec)
  z_prime= z
  P_NbyN_prime<- calculate_victory_probabilities(vec2mat_0_P(z_prime,P),P)
  n_prime = matrix(0,nrow(P),1)
  for(h in 1:K){
    n_prime[h] = sum(length(which(z_prime==h)))
  }
  
  scanning_order = sample(1:n,n, replace=F)
  
  # full sweep
  for(i_th_turn in scanning_order){
    
    z_scanning = z_prime
    #save current label of z_ii
    k_prime <- z_prime[i_th_turn]
    
    # Sample a new label using the adjusted probabilities
    if(model != 'Simple'){
      labels_to_sample = c(min(k_prime+1, K), max(k_prime-1, 1))
    }else{
      labels_to_sample = labels_available
    }
    k_scanning <- sample(x = setdiff(labels_to_sample, k_prime), size = 1, replace = F)
    
    z_scanning[i_th_turn] <- k_scanning
    
    
    #compute the likelihood of the data with the current assignment just for i_th_turn
    A_minus = sum(dbinom(Y_ij[i_th_turn,], N_ij[i_th_turn,], P_NbyN_prime[i_th_turn,], log=T)) + 
      sum(dbinom(Y_ij[,i_th_turn], N_ij[,i_th_turn], P_NbyN_prime[,i_th_turn], log=T)) 
    
    #update P_NbyN
    P_NbyN_scanning = P_NbyN_prime
    for(nodes in 1:n){
      P_NbyN_scanning[i_th_turn,nodes]<- P[k_scanning,z_scanning[nodes]]
      P_NbyN_scanning[nodes,i_th_turn]<- P[z_scanning[nodes],k_scanning]
    }
    #compute the likelihood of the same points with the new assignment
    A_plus = sum(dbinom(Y_ij[i_th_turn,], N_ij[i_th_turn,], P_NbyN_scanning[i_th_turn,], log=T)) + sum(dbinom(Y_ij[,i_th_turn], N_ij[,i_th_turn], P_NbyN_scanning[,i_th_turn], log=T)) 
    
    #Updating the likelihood
    A_scanning = A_prime - A_minus + A_plus
    
    
    n_scanning<- n_prime
    n_scanning[c(k_prime, k_scanning)] <- n_prime[c(k_prime, k_scanning)] + c(-1, 1)
    
    B_scanning<- ddirichlet_multinomial(N = n,K = K,n_k = n_scanning,my_alpha = alpha_vec)
    
    log_r= t*A_scanning - t*A_prime + B_scanning - B_prime
    #create statements that check conditiond to accept move
    GS_condition= min(log_r,0)>=log(runif(1))
    if(GS_condition){
      acc.count_z[i_th_turn]=acc.count_z[i_th_turn]+1
      z_prime<-z_scanning
      A_prime<- A_scanning
      B_prime<- B_scanning
      P_NbyN_prime <- P_NbyN_scanning
      n_prime <- n_scanning
    }
    #labels_available are the same
    #else, z_prime[ii] stays equal to z_current[ii]
  }
  z<- z_prime
  z_P<- vec2mat_0_P(z,P)
  # number of victories between block p and block q
  ybar = t(z_P)%*%(Y_ij*upper.tri(Y_ij))%*%z_P
  # number of missed victories between block p and block q
  n_minus_y1 <- (N_ij-Y_ij)*upper.tri(N_ij)
  # number of missed victories between block p and block q
  mbar<- t(z_P)%*%n_minus_y1%*%z_P
  
  coef1 = lchoose(N_ij, Y_ij)*upper.tri(N_ij)
  lamdabar <- t(z_P)%*%(coef1)%*%z_P
  
  return(list(acc.moves = acc.count_z, z_current= z, ybar =ybar, mbar=mbar, lamdabar= lamdabar))
} 








#-------------------------------------------------------------------------------
#-------------------------- UNORDERED ---------------------------------
#-------------------------------------------------------------------------------

#---------------------- Computing the probabilities ----------------------------


# Prior on P

P_prior_probability_UNORDERED = function(P,K){
  
  p_P = matrix(0,K,K)
  for(p_th in 1:K){
    for(q_th in p_th:K){
      p_P[p_th,q_th]<-  dbeta(P[p_th,q_th],1,1,log = TRUE) 
    }
  }
  
  return(sum(p_P))
}

# Proportional posterior


lprop_posterior_withP_UNORDERED <- function(lamdabar, ybar,mbar,P, alpha_vec, n_k,K){
  
  sum1 <- llik_over_blocks_f_binomial(lamdabar, ybar, mbar, P)
  #prior on P
  sum0<- P_prior_probability_UNORDERED(P,K)
  #prior on z
  sum2 <- ddirichlet_multinomial(sum(n_k),K,n_k, alpha_vec)
  #prior on K
  sum5 <- -lfactorial(K) - log(exp(1) -1)
  result <- sum0+ sum1 + sum2 + sum5 
  return(result)
}

#---------------------------- MCMC steps ---------------------------------------



P_update_f_UNORDERED = function(lamdabar,ybar,mbar,P, a, alpha_vec, n_k,
                                sigma_squared, U_vec,K, tau_P,
                                acc.count_P){
  
  P_current = P
  
  for(p_th in 1:K){
    for(q_th in p_th:K){
      
      P_prime<-P_current
      P_prime[p_th,q_th]<- rtruncnorm(1, mean = P_current[p_th,q_th],sd = tau_P[p_th,q_th], a =   0, b  = 1)
      if(p_th!=q_th){
        P_prime[q_th,p_th] <- 1 - P_prime[p_th,q_th]
      }
      #computing the proportional posterior in a'
      prop_posterior_prime <-  lprop_posterior_withP_UNORDERED(lamdabar = lamdabar,
                                                               ybar = ybar, mbar = mbar, 
                                                               P = P_prime, 
                                                               alpha_vec = alpha_vec,
                                                               n_k = n_k,K = K)
      
      #evaluating the proportional posterior in a^(t)
      prop_posterior_current <-  lprop_posterior_withP_UNORDERED(lamdabar = lamdabar,
                                                                 ybar = ybar, mbar = mbar, 
                                                                 P = P_current, 
                                                                 alpha_vec = alpha_vec,
                                                                 n_k = n_k,K = K) 
      #evaluating the proposal density g(a') 
      log_proposal_prime <- log(dtruncnorm(P_prime[p_th,q_th],mean = P_current[p_th,q_th],sd = tau_P[p_th,q_th], a =   0, b  = 1))
      
      #evaluating the proposal density g(sigma^2)^(t)) 
      log_proposal_current <- log(dtruncnorm(P_current[p_th,q_th],mean = P_prime[p_th,q_th],sd = tau_P[p_th,q_th], a =   0, b  = 1))
      
      #acceptance ratio
      log_r=  prop_posterior_prime + log_proposal_current - 
        prop_posterior_current - log_proposal_prime
      
      
      #create statements that check conditiond to accept move
      MH_condition_P_update= min(log_r,0)>=log(runif(1))
      if(MH_condition_P_update){
        acc.count_P[p_th,q_th] =acc.count_P[p_th,q_th] +1
        P_current = P_prime
      }
    }
  }
  P_current= P_prime
  
  return(list(acc.moves = acc.count_P,
              P= P_current))
  
} 


