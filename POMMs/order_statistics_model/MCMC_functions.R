
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


lprop_posterior_P_out <- function(lamdabar, ybar,mbar, a, 
                                  alpha_vec, n_k,sigma_squared, U_vec,K){
  
  sum1 <- llike_integrated_out(lamdabar, ybar,mbar,K,U_vec, sigma_squared)
  #prior on z
  sum2 <- ddirichlet_multinomial(N,K,n_k, alpha_vec)
  #prior on U
  sum3 <- lfactorial(K)+K*(-log(a-0.5))+log(ifelse(all(U_vec<a)&all(U_vec>0.5),1,0))
  #prior on sigma^2
  sum4 <- -sigma_squared - log(1-exp(-max(U_vec)*(1-max(U_vec))))
  #prior on K
  sum5 <- -lfactorial(K) - log(exp(1) -1)
  #pior on a
  sum6 <- 1/2
  result <- sum1 + sum2 + sum3 + sum4 + sum5 + sum6
  return(result)
}

#-------------------------- MCMC steps -----------------------------------------
a_update_f_P_out = function(lamdabar, ybar,mbar, a, alpha_vec, n_k,
                            sigma_squared, U_vec,K, tau_a,
                            acc.count_a){
  
  
  
  #simulating (sigma^2)' from a g ~ truncated normal
  a_prime <- rtruncnorm(1,a = 0.5,b = 1,
                        mean = a,
                        sd = sqrt(tau_a))
  
  #computing the proportional posterior in a'
  prop_posterior_prime <- lprop_posterior(lamdabar, ybar,mbar, a_prime, 
                                          alpha_vec, n_k,sigma_squared, 
                                          U_vec,K)
  
  #evaluating the proportional posterior in a^(t)
  prop_posterior_current<- lprop_posterior(lamdabar, ybar,
                                           mbar, a, alpha_vec, n_k,
                                           sigma_squared, U_vec,K)
  
  #evaluating the proposal density g(a') 
  log_proposal_prime <- log(dtruncnorm(a_prime,a = 0.5,b = 1,
                                       mean = a,
                                       sd = tau_a**2))
  
  #evaluating the proposal density g(sigma^2)^(t)) 
  log_proposal_current <- log(dtruncnorm(a,a = 0.5,b = 1,
                                         mean = a_prime,
                                         sd = tau_a**2))
  #acceptance ratio
  log_r=  prop_posterior_prime + log_proposal_current - 
    prop_posterior_current - log_proposal_prime
  
  #create statements that check conditionion to accept move
  MH_condition_a= min(log_r,0)>=log(runif(1))
  if(MH_condition_a){
    acc.count_a = acc.count_a + 1
    a <- a_prime
  }
  
  
  return(list(acc.moves = acc.count_a,
              a=a))
  
} 




U_update_f_P_out = function(lamdabar, ybar,mbar, a, alpha_vec, n_k,
                            sigma_squared, U_vec,K, tau_U_vec,
                            acc.count_a){
  # 
  # for(k in 2:(K+1)){
  
  lb<- max(0.5*(1-sqrt(1-4*sigma_squared)),0.5)
  ub<- min(0.5*(1+sqrt(1-4*sigma_squared)),a)
  # U_ext <- c(lb,U_vec,ub)
  # U_vec_prime <- U_vec
  
  
  #simulating (sigma^2)' from a g ~ truncated normal
  # U_prime <- rtruncnorm(1,a = U_ext[k-1],b = U_ext[k+1],
  #                       mean = U_ext[k],
  #                       sd = tau_U_vec[k-1]**2)
  
  U_prime <- runif(K,lb,ub)
  U_vec_prime <- sort(U_prime)
  
  #computing the proportional posterior in a'
  prop_posterior_prime <- lprop_posterior(lamdabar, ybar,mbar, a, 
                                          alpha_vec, n_k,sigma_squared, 
                                          U_vec_prime,K)
  
  #evaluating the proportional posterior in a^(t)
  prop_posterior_current<- lprop_posterior(lamdabar, ybar,
                                           mbar, a, alpha_vec, n_k,
                                           sigma_squared, U_vec,K)
  
  #evaluating the proposal density g(a') 
  log_proposal_prime <- log(dtruncnorm(U_prime,a = U_ext[k-1],b = U_ext[k+1],
                                       mean = U_ext[k],
                                       sd = tau_U_vec[k-1]**2))
  
  #evaluating the proposal density g(sigma^2)^(t)) 
  log_proposal_current <- log(dtruncnorm(U_ext[k],a = U_ext[k-1],b = U_ext[k+1],
                                         mean = U_ext[k],
                                         sd = tau_U_vec[k-1]**2))
  #acceptance ratio
  # log_r=  prop_posterior_prime + log_proposal_current - 
  #   prop_posterior_current - log_proposal_prime
  # 
  log_r =  prop_posterior_prime  - prop_posterior_current 
  
  #create statements that check conditionion to accept move
  MH_condition_U_vec= min(log_r,0)>=log(runif(1))
  if(MH_condition_U_vec){
    acc.count_a = acc.count_a+1
    # U_vec[k-1] <- U_prime
    U_vec <- U_prime
  }
  # }
  
  return(list(acc.moves = acc.count_a,
              U_vec=U_vec))
  
} 



sigma_squared_update_f_P_out = function(lamdabar,ybar,mbar, a, alpha_vec, n_k,
                                        sigma_squared, U_vec,K, tau_sigma_squared,
                                        acc.count_sigma_squared){
  
  
  ub<- min(U_vec*(1-U_vec))
  #simulating (sigma^2)' from a g ~ truncated normal
  sigma_squared_prime <- rtruncnorm(1,a = 0,b = ub ,
                                    mean = sigma_squared,
                                    sd = sqrt(tau_sigma_squared))
  
  #computing the proportional posterior in (sigma^2)'
  prop_posterior_prime <- lprop_posterior(lamdabar, ybar,mbar, a, alpha_vec, 
                                          n_k,sigma_squared_prime, U_vec,K)
  
  #evaluating the proportional posterior in (sigma^2)^(t)
  prop_posterior_current<- lprop_posterior(lamdabar, ybar,mbar, a, alpha_vec, 
                                           n_k,sigma_squared, U_vec,K)
  
  #evaluating the proposal density g(sigma^2)') 
  log_proposal_prime <- log(dtruncnorm(sigma_squared_prime,a = 0,b = ub,
                                       mean = sigma_squared,
                                       sd = tau_sigma_squared**2))
  
  #evaluating the proposal density g(sigma^2)^(t)) 
  log_proposal_current <- log(dtruncnorm(sigma_squared,a = 0,b = ub,
                                         mean = sigma_squared,
                                         sd = tau_sigma_squared**2))
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













#-------------------------------------------------------------------------------
#-------------------------- P within the model ---------------------------------
#-------------------------------------------------------------------------------

#---------------------- Computing the probabilities ----------------------------


# Likelihood

llik_over_blocks_f_binomial = function(lamdabar, ybar, mbar, P){
  #llik<- sum(bin_coef+ ybar*log(P) + (mbar)*log(1-P))
  llik = matrix(0,K,K)
  for(p in 1:K){
    for(q in p:K){
      if(q-p==0){
        llik[p,q] = lamdabar[p,q]+ ybar[p,q]*log(P[p,q]) + mbar[p,q]*log(1-P[p,q])
      }else{
        llik[p,q] = lamdabar[q,p]+ lamdabar[p,q]+(ybar[p,q]+mbar[q,p])*log(P[p,q]) + (mbar[p,q]+ybar[q,p])*log(1-P[p,q]) 
      }
    }
  }
  return(sum(llik))
}

# Prior on P

P_prior_probability = function(P,K,U_vec, sigma_squared){
  
  beta_params = beta_mean_var(mu = U_vec,var = rep(sigma_squared,K-1) )
  alpha_k = beta_params$alpha
  beta_k = beta_params$beta
  
  p_P = matrix(0,K,K)
  for(p_th in 1:K){
    for(q_th in p_th:K){
      if(q_th-p_th != 0){
        p_P[p_th,q_th]<-  dbeta(P[p_th,q_th],alpha_k[q_th-p_th],beta_k[q_th-p_th],log = TRUE) 
      }else if (p_th==q_th){
        p_P[p_th,p_th]<-  dbeta(P[p_th,q_th],1,1,log = TRUE) 
      }
    }
  }
  
  return(sum(p_P))
}

# Proportional posterior


lprop_posterior_withP <- function(lamdabar, ybar,mbar,P, a, 
                                  alpha_vec, n_k,sigma_squared, U_vec,K){
  
  sum1 <- llik_over_blocks_f_binomial(lamdabar, ybar, mbar, P)
  #prior on P
  sum0<- P_prior_probability(P,K,U_vec, sigma_squared)
  #prior on z
  sum2 <- ddirichlet_multinomial(sum(n_k),K,n_k, alpha_vec)
  #prior on U
  sum3 <- lfactorial(K-1)+(K-1)*(-log(a-0.5))+log(ifelse(all(U_vec<a)&all(U_vec>0.5),1,0))
  #prior on sigma^2
  sum4 <- -sigma_squared - log(1-exp(-min(U_vec*(1-U_vec))))
  #prior on K
  sum5 <- -lfactorial(K) - log(exp(1) -1)
  #pior on a
  sum6 <- dtruncnorm(x = a,mean = 0.75,sd =0.25,a=0.5,b = 1)
  result <- sum0+ sum1 + sum2 + sum3 + sum4 + sum5 + sum6
  return(result)
}

#-------------------------- MCMC steps -----------------------------------------

P_update_f = function(lamdabar,ybar,mbar,P, a, alpha_vec, n_k,
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
      prop_posterior_prime <-  lprop_posterior_withP(lamdabar = lamdabar, ybar = ybar, mbar = mbar, P = P_prime,
                                                     K = K,U_vec = U_vec,sigma_squared = sigma_squared,a = a,alpha_vec = alpha_vec,n_k = n_k)
      
      #evaluating the proportional posterior in a^(t)
      prop_posterior_current<- lprop_posterior_withP(lamdabar = lamdabar, ybar = ybar, mbar = mbar, P = P_current,
                                                     K = K,U_vec = U_vec,sigma_squared = sigma_squared,a = a,alpha_vec = alpha_vec,n_k = n_k)
      
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


a_update_f_withP = function(lamdabar, ybar,mbar, P,a, alpha_vec, n_k,
                            sigma_squared, U_vec,K, tau_a,
                            acc.count_a){
  
  
  #simulating (sigma^2)' from a g ~ truncated normal
  a_prime <- rtruncnorm(1,a = 0.5,b = 1,
                        mean = a,
                        sd = tau_a)
  
  #computing the proportional posterior in a'
  prop_posterior_prime <-  lprop_posterior_withP(lamdabar, ybar,mbar,P, a_prime, 
                                                 alpha_vec, n_k,sigma_squared, U_vec,K)
  
  
  #evaluating the proportional posterior in a^(t)
  prop_posterior_current<- lprop_posterior_withP(lamdabar, ybar,mbar,P, a, 
                                                 alpha_vec, n_k,sigma_squared, U_vec,K)
  
  #evaluating the proposal density g(a') 
  log_proposal_prime <- log(dtruncnorm(a_prime,a = 0.5,b = 1,
                                       mean = a,
                                       sd = tau_a))
  
  #evaluating the proposal density g(sigma^2)^(t)) 
  log_proposal_current <- log(dtruncnorm(a,a = 0.5,b = 1,
                                         mean = a_prime,
                                         sd = tau_a))
  #acceptance ratio
  log_r=  prop_posterior_prime + log_proposal_current - 
    prop_posterior_current - log_proposal_prime
  
  #create statements that check conditionion to accept move
  MH_condition_a= min(log_r,0)>=log(runif(1))
  if(MH_condition_a){
    acc.count_a = acc.count_a + 1
    a <- a_prime
  }
  
  
  return(list(acc.moves = acc.count_a,
              a=a))
  
} 




U_update_f_withP = function(lamdabar, ybar,mbar,P, a, alpha_vec, n_k,
                            sigma_squared, U_vec,K, tau_U_vec,
                            acc.count_U){
  # 
  # for(k in 2:(K+1)){
  
  lb<- max(0.5*(1-sqrt(1-4*sigma_squared)),0.5)
  ub<- min(0.5*(1+sqrt(1-4*sigma_squared)),a)
  # U_ext <- c(lb,U_vec,ub)
  # U_vec_prime <- U_vec
  
  
  #simulating (sigma^2)' from a g ~ truncated normal
  # U_prime <- rtruncnorm(1,a = U_ext[k-1],b = U_ext[k+1],
  #                       mean = U_ext[k],
  #                       sd = tau_U_vec[k-1]**2)
  
  U_unif_prime <- runif((K-1),lb,ub)
  U_vec_prime <- sort(U_unif_prime)
  
  #computing the proportional posterior in a'
  prop_posterior_prime <- lprop_posterior_withP(lamdabar,ybar,mbar,P,
                                                a,alpha_vec,n_k,sigma_squared,
                                                U_vec = U_vec_prime,K)
  
  #evaluating the proportional posterior in a^(t)
  prop_posterior_current <- lprop_posterior_withP(lamdabar,ybar,mbar,
                                                  P,a,alpha_vec,
                                                  n_k,sigma_squared,U_vec,K)
  
  #evaluating the proposal density g(a') 
  # log_proposal_prime <- log(dtruncnorm(U_prime,a = U_ext[k-1],b = U_ext[k+1],
  #                                      mean = U_ext[k],
  #                                      sd = tau_U_vec[k-1]**2))
  # 
  # #evaluating the proposal density g(sigma^2)^(t)) 
  # log_proposal_current <- log(dtruncnorm(U_ext[k],a = U_ext[k-1],b = U_ext[k+1],
  #                                        mean = U_ext[k],
  #                                        sd = tau_U_vec[k-1]**2))
  #acceptance ratio
  # log_r=  prop_posterior_prime + log_proposal_current - 
  #   prop_posterior_current - log_proposal_prime
  # 
  log_r =  prop_posterior_prime  - prop_posterior_current 
  
  #create statements that check conditionion to accept move
  MH_condition_U_vec= min(log_r,0)>=log(runif(1))
  if(MH_condition_U_vec){
    acc.count_U = acc.count_U+1
    # U_vec[k-1] <- U_prime
    U_vec <- U_vec_prime
  }
  # }
  
  return(list(acc.moves = acc.count_U,
              U_vec=U_vec))
  
} 





sigma_squared_update_f_withP = function(lamdabar,ybar,mbar,P, a, alpha_vec, n_k,
                                        sigma_squared, U_vec,K, tau_sigma_squared,
                                        acc.count_sigma_squared){
  
  
  ub<- min(U_vec*(1-U_vec))
  #simulating (sigma^2)' from a g ~ truncated normal
  sigma_squared_prime <- rtruncnorm(1,a = 0,b = ub ,
                                    mean = sigma_squared,
                                    sd = tau_sigma_squared)
  
  #computing the proportional posterior in (sigma^2)'
  prop_posterior_prime <- lprop_posterior_withP(lamdabar, ybar,mbar,P, a, 
                                                alpha_vec, n_k,sigma_squared_prime, U_vec,K)
  #evaluating the proportional posterior in (sigma^2)^(t)
  prop_posterior_current<- lprop_posterior_withP(lamdabar, ybar,mbar,P, a, 
                                                 alpha_vec, n_k,sigma_squared, U_vec,K)
  
  #evaluating the proposal density g(sigma^2)') 
  log_proposal_prime <- log(dtruncnorm(sigma_squared_prime,a = 0,b = ub,
                                       mean = sigma_squared,
                                       sd = tau_sigma_squared))
  
  #evaluating the proposal density g(sigma^2)^(t)) 
  log_proposal_current <- log(dtruncnorm(sigma_squared,a = 0,b = ub,
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
                            K, tau_z,
                            acc.count_z,labels_available){
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





z_update_f_withP1 = function(N_ij, Y_ij, z,lamdabar,ybar,mbar, P, a, alpha_vec, n_k,
                             sigma_squared, U_vec,K, tau_z,
                             acc.count_z,labels_available){
  
  n<- nrow(N_ij)
  scanning_order = sample(1:n,n, replace=F)
  
  # full sweep
  for(ii in scanning_order){
    
    z_prime<- z
    z_ii <- z_prime[ii]
    
    z_prime[ii]<- sample(x = setdiff(labels_available,z_ii),size = 1)
    
    z_P<- vec2mat_0_P(z_prime,P)
    # number of victories between block p and block q
    ybar_prime = t(z_P)%*%(Y_ij*upper.tri(Y_ij))%*%z_P
    # number of missed victories between block p and block q
    n_minus_y1 <- (N_ij-Y_ij)*upper.tri(N_ij)
    # number of missed victories between block p and block q
    mbar_prime<- t(z_P)%*%n_minus_y1%*%z_P
    
    coef1 = lchoose(N_ij, Y_ij)*upper.tri(N_ij)
    lamdabar_prime <- t(z_P)%*%(coef1)%*%z_P
    
    #computing the proportional posterior in (sigma^2)'
    prop_posterior_prime <- lprop_posterior_withP(lamdabar_prime, ybar_prime,mbar_prime,P, a, 
                                                  alpha_vec, n_k,sigma_squared, U_vec,K)
    #evaluating the proportional posterior in (sigma^2)^(t)
    prop_posterior_current<- lprop_posterior_withP(lamdabar, ybar,mbar,P, a, 
                                                   alpha_vec, n_k,sigma_squared, U_vec,K)
    #acceptance ratio
    log_r=  prop_posterior_prime  - prop_posterior_current 
    
    #create statements that check conditionion to accept move
    MH_condition_z= min(log_r,0)>=log(runif(1))
    if(MH_condition_z){
      lamdabar = lamdabar_prime
      mbar=mbar_prime
      ybar = ybar_prime
    }
    
  }
  
  
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


