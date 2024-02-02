
#-------------------------------------------------------------------------------
############################ MCMC functions ####################################
#-------------------------------------------------------------------------------

tuning_proposal<- function(iteration, acceptance_count, sigma, acceptanceTarget, min_sigma){
  #compute acceptance rate
  acceptanceRate <- acceptance_count / iteration
  #setting the change in the variance
  delta= min(0.01, iteration**(-1/2))
  #passing top the log scale, to have a finer scale
  lsi= log(sigma)
  #if we are accepting too much ==> increase variance, otherwise, reduce it
  if (acceptanceRate > acceptanceTarget) {
    lsi <- lsi + delta
  } else {
    lsi <- lsi - delta}
  
  # Update the proposal standard deviations
  sigma_updated <- max(min_sigma,exp(lsi))
  return(sigma_updated)
}


label_probability <- function(distance, sigma_z) {
  prob <- dnorm(distance, mean = 0, sd = sigma_z)
  return(prob)
}


inverse_logit_f = function(x){
  y= exp(x)/(1+exp(x))
  return(y)
}
#-------------------------------------------------------------------------------
#-------------------------- P integrated out -----------------------------------
#-------------------------------------------------------------------------------


#---------------------- Computing the probabilities ----------------------------

log_lik_f_binom = function(N,Y,z,P, directed=T){
  z_P<- vec2mat_0_P(z,P)
  P_nbyn<- calculate_victory_probabilities(z_P, P)
  if(directed==T){
    #computing the pairwise log-probabilitiees
    bigM = lchoose(N,Y)+(Y* log(P_nbyn)+(N-Y)*log(1 - P_nbyn))
    #remember to subtract the diagonal
    log_lik= sum(bigM) - sum(diag(bigM))
  }else if(directed==F){
    bigM = lchoose(N,Y)+(Y* log(P_nbyn)+(N-Y)*log(1 - P_nbyn))
    #remember to subtract the diagonal
    log_lik= sum(bigM*upper.tri(bigM))
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

llik_over_blocks_f_binomial = function(lamdabar, ybar, mbar, P){
  #llik<- sum(bin_coef+ ybar*log(P) + (mbar)*log(1-P))
  llik = matrix(0,K,K)
  for(diag_iii in 0:(K-1)){
    llik[col(P)-row(P)==diag_iii] <- lamdabar[col(P)-row(P)==diag_iii]+ ybar[col(P)-row(P)==diag_iii]*P[col(P)-row(P)==diag_iii] -
      (mbar[col(P)-row(P)==diag_iii]+ybar[col(P)-row(P)==diag_iii])*log(1+exp(P[col(P)-row(P)==diag_iii]))
  }
  
  return(sum(llik))
}

# Density for P

P_prior_probability = function(P,K,mu_vec, sigma_squared, model){
  p_P = matrix(0,K,K)
  if(model == "WST"){
    p_P[col(P)-row(P)==0] <- dnorm(x = P[col(P)-row(P)==0],mean = mu_vec[1],sd = sqrt(sigma_squared))
    for(diag_i in 1:(K-1)){
      p_P[col(P)-row(P)==diag_i] <- truncnorm::dtruncnorm(x = P[col(P)-row(P)==diag_i],a = 0,b=Inf,mean = mu_vec[diag_i+1],sd = sqrt(sigma_squared))
    }
  }else if (model == "Simple"){
    theta = inverse_logit_f(P)
    p_P = dbeta(theta,1,1)*exp(theta)/((1+exp(theta))**2)
  }
  log_p_P = sum(log(p_P[upper.tri(p_P,diag = T)]))
  
  return(log_p_P)
}

#mu prior distribution
d_sA_mu = function(K,mu_vec){
  fac <- lfactorial(K)
  joint_density<- sum(log(dtruncnorm(mu_vec,a = 0,mean = -1,sd = 1)) - log(1- pnorm(0,mean = -1,sd = 1)))
  return(joint_density+fac)
}


# Proportional posterior


lprop_posterior_withP <- function(lamdabar, ybar,mbar,P, 
                                  alpha_vec, n_k,sigma_squared, mu_vec,K, model){
  
  if(model=='Simple'){
    #log likelihood
    log_lik <- llik_over_blocks_f_binomial(lamdabar = lamdabar,ybar =  ybar,  mbar = mbar,P =  P)
    #log prior on z
    prior_z <- ddirichlet_multinomial(sum(n_k),K,n_k, alpha_vec)
    #log prior on P
    prior_P<- P_prior_probability(P,K,mu_vec, sigma_squared,model)
    #computing the whole log proportional posterior
    results<- log_lik+ prior_P + prior_z 
  }else if (model=='WST'){
    #log likelihood
    log_lik <- llik_over_blocks_f_binomial(lamdabar = lamdabar,ybar =  ybar,  mbar = mbar,P =  P)
    #log prior on z
    prior_z <- ddirichlet_multinomial(sum(n_k),K,n_k, alpha_vec)
    #log prior on P
    prior_P<- P_prior_probability(P,K,mu_vec, sigma_squared,model)
    #log prior on mu
    hyperprior_mu <- d_sA_mu(K, mu_vec)
    #log prior on sigma^2
    hyperprior_sigmasquared <- LaplacesDemon::dinvgamma(x = sigma_squared,shape = 0.001,scale = 0.001,log = T)
    #computing the whole log proportional posterior
    results<- log_lik+ prior_P + prior_z + hyperprior_mu + hyperprior_sigmasquared 
  }else if(model =='SST'){
    #log likelihood
    log_lik <- llik_over_blocks_f_binomial(lamdabar = lamdabar,ybar =  ybar,  mbar = mbar,P =  P)
    #log prior on z
    prior_z <- ddirichlet_multinomial(sum(n_k),K,n_k, alpha_vec)
    #log prior on mu
    hyperprior_mu <- d_sA_mu(K, mu_vec)
    #computing the whole log proportional posterior
    results <- log_lik + prior_z + hyperprior_mu  
  }
  return(results)
}

#-------------------------- MCMC steps -----------------------------------------

P_update_f = function(lamdabar,ybar,mbar,P, alpha_vec, n_k,
                      sigma_squared, mu_vec,K, tau_P,
                      acc.count_P, model){
  tau_P = matrix(.2,K,K)
  P_current<- P
  #Updating each entry of P, one at the time
  
  ut <- upper.tri(P,diag = T) # get the logical matrix for upper triangular elements
  Pcombn = which(ut, arr.ind = TRUE) # get the indices of the upper triangular elements
  uo<- data.frame(Pcombn[sample(nrow(Pcombn)), ])# permuting the order of the rows
  
  for(i_th in 1:nrow(uo)){
    
    P_prime<-P_current
    #P' ~ g(P^(t), tau_P) 
    P_prime[uo$row[i_th],uo$col[i_th]]<- rtruncnorm(1, mean = P_current[uo$row[i_th],uo$col[i_th]],
                                                    sd = tau_P[uo$row[i_th],uo$col[i_th]], 
                                                    a = ifelse(uo$row[i_th]==uo$col[i_th],-10,0), b  = 10)
    
    if(uo$row[i_th] != uo$col[i_th]){
      p_ij<- inverse_logit_f(P_prime[uo$row[i_th],uo$col[i_th]])
      P_prime[uo$col[i_th],uo$row[i_th]] <- log((1-p_ij)/p_ij)                  
    }
    
    
    #computing the proportional posterior in P'
    prop_posterior_prime <-  lprop_posterior_withP(lamdabar = lamdabar, ybar = ybar, mbar = mbar, P = P_prime,
                                                   K = K,mu_vec = mu_vec,sigma_squared = sigma_squared,
                                                   alpha_vec = alpha_vec,n_k = n_k,model)
    
    #evaluating the proportional posterior in P^(t)
    prop_posterior_current<- lprop_posterior_withP(lamdabar = lamdabar, ybar = ybar, mbar = mbar, P = P_current,
                                                   K = K,mu_vec = mu_vec,sigma_squared = sigma_squared,
                                                   alpha_vec = alpha_vec,n_k = n_k,model)
    
    #evaluating the proposal density g(P'| P^(t)) 
    log_proposal_prime <- log(dtruncnorm(P_prime[uo$row[i_th],uo$col[i_th]],
                                         mean = P_current[uo$row[i_th],uo$col[i_th]],
                                         sd = tau_P[uo$row[i_th],uo$col[i_th]], 
                                         a =   ifelse(uo$row[i_th]==uo$col[i_th],-10,0), 
                                         b  = 10))
    
    #evaluating the proposal density g(P^(t)| P') 
    log_proposal_current <- log(dtruncnorm(P_current[uo$row[i_th],uo$col[i_th]],
                                           mean = P_prime[uo$row[i_th],uo$col[i_th]],
                                           sd = tau_P[uo$row[i_th],uo$col[i_th]], 
                                           a =  ifelse(uo$row[i_th]==uo$col[i_th],-10,0),
                                           b  = 10))
    
    #acceptance ratio
    log_r=  prop_posterior_prime + log_proposal_current - 
      prop_posterior_current - log_proposal_prime
    
    
    #create statements that check conditiond to accept move
    MH_condition_P_update= min(log_r,0)>=log(runif(1))
    if(MH_condition_P_update){
      acc.count_P[uo$row[i_th],uo$col[i_th]] =acc.count_P[uo$row[i_th],uo$col[i_th]] +1
      P_current = P_prime
    }
    
  }
  P_current= P_prime
  
  return(list(acc.moves = acc.count_P,
              P= P_current))
  
} 





mu_update_f_withP = function(lamdabar, ybar,mbar,P, alpha_vec, n_k,
                             sigma_squared, mu_vec,tau_mu_vec, K,
                             acc.count_mu_vec,model){
  
  #computing the proportional posterior in mu' ~ g(mu^(t), tau_mu_vec)
  mu_prime <- truncnorm::rtruncnorm(K,a = -1,b = 10, mean = mu_vec,sd = tau_mu_vec)
  mu_vec_prime <- sort(mu_prime)
  P_prime = P
  
  if(model =='SST'){
    P_prime<-matrix(NA, K,K)
    for(i in 0:(K-1)){
      P_prime[col(P_prime)-row(P_prime)==(i)] <- rep(mu_vec_prime[i+1],(K-i))
      p_ij<- inverse_logit_f(rep(mu_vec_prime[i+1],(K-i)))
      P_prime[col(P_prime)-row(P_prime)==(-i)] <- log((1-p_ij)/p_ij)   
    }
    
    
    
  }
  
  #computing the proportional posterior in mu'
  prop_posterior_prime <- lprop_posterior_withP(lamdabar = lamdabar,ybar = ybar,mbar = mbar,P = P_prime,alpha_vec = alpha_vec,
                                                n_k = n_k,sigma_squared = sigma_squared,mu_vec = mu_vec_prime,K = K,model=model)
  
  #evaluating the proportional posterior in mu^(t)
  prop_posterior_current <- lprop_posterior_withP(lamdabar = lamdabar,ybar = ybar,mbar = mbar,P = P,alpha_vec = alpha_vec,
                                                  n_k = n_k,sigma_squared = sigma_squared,mu_vec = mu_vec,K = K,model = model)
  
  #evaluating the proposal density g(mu'| mu^(t)) 
  log_proposal_prime <- log(dtruncnorm(mu_prime,mean = mu_vec,sd =tau_mu_vec, a =   -1, b  = Inf))
  
  #evaluating the proposal density g(mu^(t)| mu') 
  log_proposal_current <- log(dtruncnorm(mu_vec,mean = mu_prime,sd = tau_mu_vec, a =   -1, b  = Inf))
  
  log_r =  prop_posterior_prime  - prop_posterior_current + log_proposal_current - log_proposal_prime
  
  #create statements that check conditionion to accept move
  MH_condition_mu_vec= min(log_r,0)>=log(runif(1))
  if(MH_condition_mu_vec){
    acc.count_mu_vec = acc.count_mu_vec+1
    # U_vec[t] <- U_prime
    mu_vec <- mu_vec_prime
    if(model=='SST'){
      P<-P_prime
    }
  }
  
  return(list(acc.moves = acc.count_mu_vec,
              mu_vec=mu_vec,
              P= P))
  
} 





sigma_squared_update_f_withP = function(lamdabar,ybar,mbar,P, alpha_vec, n_k,
                                        sigma_squared, mu_vec,K, tau_sigma_squared,
                                        acc.count_sigma_squared, model){
  
  
  #simulating (sigma^2)' from a g ~ truncated normal
  sigma_squared_prime <- rtruncnorm(1,a = 0,b = 1 ,
                                    mean = sigma_squared,
                                    sd = tau_sigma_squared)
  
  #computing the proportional posterior in (sigma^2)'
  prop_posterior_prime <- lprop_posterior_withP(lamdabar, ybar,mbar,P, 
                                                alpha_vec, n_k,sigma_squared_prime, mu_vec,K,model = model)
  #evaluating the proportional posterior in (sigma^2)^(t)
  prop_posterior_current<- lprop_posterior_withP(lamdabar, ybar,mbar,P,
                                                 alpha_vec, n_k,sigma_squared, mu_vec,K,model=model)
  
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
                            acc.count_z,labels_available,model){
  P<- inverse_logit_f(P)
  n<- nrow(N_ij)
  A_prime<- log_lik_f_binom(N = N_ij,Y = Y_ij,z =z,P = P,directed = T)
  B_prime<- ddirichlet_multinomial(N = n,K = K,n_k = n_k, my_alpha =  alpha_vec)
  z_prime= z
  P_NbyN_prime<- calculate_victory_probabilities(vec2mat_0_P(z_prime,P),P)
  n_prime = matrix(0,nrow(P),1)
  for(h in 1:nrow(P)){
    n_prime[h] = sum(length(which(z_prime==h)))
  }
  
  scanning_order = sample(1:n,n, replace=F)
  
  # full sweep
  for(ii in scanning_order){
    
    z_scanning = z_prime
    #save current label of z_ii
    k_prime <- z_prime[ii]
    
    # Sample a new label using the adjusted probabilities
    labels_to_sample = c(min(k_prime+1, K), max(k_prime-1, 1))
    k_scanning <- sample(x = setdiff(labels_to_sample, k_prime), size = 1, replace = F)
    
    z_scanning[ii] <- k_scanning
    
    
    #compute the likelihood of the data with the current assignment just for ii
    A_minus = sum(dbinom(Y_ij[ii,], N_ij[ii,], P_NbyN_prime[ii,], log=T)) + sum(dbinom(Y_ij[,ii], N_ij[,ii], P_NbyN_prime[,ii], log=T)) 
    
    #update P_NbyN
    P_NbyN_scanning = P_NbyN_prime
    for(nodes in 1:n){
      P_NbyN_scanning[ii,nodes]<- P[k_scanning,z_scanning[nodes]]
      P_NbyN_scanning[nodes,ii]<- P[z_scanning[nodes],k_scanning]
    }
    #compute the likelihood of the same points with the new assignment
    A_plus = sum(dbinom(Y_ij[ii,], N_ij[ii,], P_NbyN_scanning[ii,], log=T)) + sum(dbinom(Y_ij[,ii], N_ij[,ii], P_NbyN_scanning[,ii], log=T)) 
    
    #Updating the likelihood
    A_scanning = A_prime - A_minus + A_plus
    
    
    n_scanning<- n_prime
    n_scanning[c(k_prime, k_scanning)] <- n_prime[c(k_prime, k_scanning)] + c(-1, 1)
    
    B_scanning<- ddirichlet_multinomial(N = n,K = K,n_k = n_scanning,my_alpha = alpha_vec)
    
    log_r= A_scanning - A_prime + B_scanning - B_prime
    #create statements that check conditiond to accept move
    GS_condition= min(log_r,0)>=log(runif(1))
    if(GS_condition){
      acc.count_z[ii]=acc.count_z[ii]+1
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


