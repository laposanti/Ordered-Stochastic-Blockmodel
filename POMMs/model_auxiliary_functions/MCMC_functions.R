
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
# 
# theta_prior_probability = function(theta,K, mu_vec, model,diag0.5=F){
#   
#   if(model == "WST"){
#     
#     
#     if(diag0.5 == F){
#       P = inverse_logit_f(theta)
#       P_upper.tri = P[upper.tri(P,diag = F)]
#       
#       p_prod = dunif(P_upper.tri,0.5,0.9999)*
#         (exp(P_upper.tri)/((1+exp(P_upper.tri))**2))
#       
#       f_diag = dunif(diag(P),0.0001,0.9999)*
#         (exp(diag(P))/((1+exp(diag(P)))**2))
#       
#       log_p_sum = sum(log(p_prod) + log(f_diag))
#     }else if(diag0.5 == T){
#       P = inverse_logit_f(theta)
#       P_upper.tri = P[upper.tri(P,diag = F)]
#       
#       p_prod = dunif(P_upper.tri,0.5,0.9999)*
#         (exp(P_upper.tri)/((1+exp(P_upper.tri))**2))
#       
#       log_p_sum = sum(log(p_prod))
#       
#     }
#   }else if (model == "Simple"){
#     P = inverse_logit_f(theta)
#     P_upper.tri = P[upper.tri(P,diag = F)]
#     
#     p_prod = dunif(P_upper.tri,0.0001,0.9999)*
#       (exp(P_upper.tri)/((1+exp(P_upper.tri))**2))
#     
#     log_p_sum = sum(log(p_prod))
#     
#   }else if(model =='SST'){
#     if(diag0.5==T){
#       fac <- lfactorial(K-1)
#       joint_density<- sum(log(dtruncnorm(mu_vec[2:(K)],a = 0,mean = 0,sd = 1)) - 
#                             log(1- pnorm(0,mean = 0,sd = 1)))  
#       log_p_sum = joint_density + fac 
#     }else if(diag0.5==F){
#       fac <- lfactorial(K-1)
#       joint_density<- sum(log(dtruncnorm(mu_vec[2:(K)],a = 0,mean = 0,sd = 1)) - 
#                             log(1- pnorm(0,mean = 0,sd = 1)))  
#       P = inverse_logit_f(theta)
#       f_diag = dunif(diag(P),0.0001,0.9999)*
#         (exp(diag(P))/((1+exp(diag(P)))**2))
#       log_p_sum = joint_density + fac + sum(log(f_diag))
#     }
#     
#     
#   }
#   
#   
#   return(log_p_sum)
# }








lprop_posterior <- function(Y_ij, N_ij, z,theta,
                            alpha_vec,mu_vec, 
                            n_k,K,t,common_indices){
  
  # Compute log likelihood only if llik is not provided
  log_lik <- ll_naive(z = z, theta = theta, Y_ij = Y_ij, N_ij = N_ij,common_indices) * t
  
  #log prior on P
  prior_theta<- theta_prior_probability(theta = theta, K = K, mu_vec = mu_vec)
  
  #log prior on z
  prior_z <- ddirichlet_multinomial(N = sum(n_k), K = K,n_k = n_k, my_alpha = alpha_vec)
  
  #computing the whole log proportional posterior
  results <- log_lik + prior_z  + prior_theta
  
  return(results)
}


#-------------------------- MCMC wrappers --------------------------------------
########################### MCMC-wrappers ######################################
#-------------------------- MCMC wrappers --------------------------------------

theta_update_f = function(Y_ij, N_ij,z, theta, alpha_vec, n_k, mu_vec,K, tau_theta,
                          acc.count_theta,common_indices, model,t,diag0.5){
  
  #number of nodes
  n<- nrow(N_ij)
  
  #storing current values
  theta_prime <- theta
  P_prime<- inverse_logit_f(theta_prime)
  
  P_NbyN_prime<- calculate_victory_probabilities(vec2mat_0_P(z,P_prime),P_prime)
  
  A_prime<-   sum(dbinom(Y_ij[common_indices],
                         size = N_ij[common_indices], 
                         P_NbyN_prime[common_indices],log = T))
  
  C_prime<- theta_prior_probability(theta = theta_prime, K=K,
                                    mu_vec = mu_vec)
  
  #Updating each entry of P, one at the time
  
  ut <- upper.tri(theta, diag= !diag0.5)
  theta_combn = which(ut, arr.ind = TRUE) # get the indices of the upper triangular elements
  n_entries = nrow(theta_combn)
  uo<- data.frame(row = theta_combn[,1], col = theta_combn[,2] )# permuting the order of the rows
  
  for(i_th in 1:nrow(uo)){
    
    theta_scanning <- theta_prime
    P_NbyN_scanning = P_NbyN_prime
    
    i_star<- uo$row[i_th]
    j_star<- uo$col[i_th]
    # print(paste0('theta_prime = ', theta_prime[i_star, j_star]))
    
    if(model == 'WST'& i_star != j_star){
      lower.bound = 0
      upper.bound = 9.21023
    }else{
      lower.bound = -9.21023
      upper.bound = 9.21023
    }
    #saving for convenience, to avoid multiple computations
    theta_ij_prime<- theta_prime[i_star,j_star]
    
    #proposing a new value for theta p_q
    theta_ij_scanning <- rtruncnorm(1, a = lower.bound , b = upper.bound,
                                    mean = theta_ij_prime, 
                                    sd =  tau_theta[i_star,j_star])
    
    #the items in cluster i_star
    Z_i_star = which(z == i_star)
    #the items in cluster j_star
    Z_j_star = which(z == j_star)
    
    
    logical_matrix1 = matrix(FALSE, n,n)
    logical_matrix1[Z_i_star,]<-TRUE
    logical_matrix1[,Z_j_star]<-TRUE
    
    logical_matrix2 = matrix(FALSE, n,n)
    logical_matrix2[Z_j_star,]<-TRUE
    logical_matrix2[,Z_i_star]<-TRUE
    # Get the upper triangular indices for the relevant clusters
    
    
    # Create a matrix that is TRUE only at the positions that are both in the upper triangle and in the selected clusters
    filtering_matrix = (logical_matrix1|logical_matrix2)*common_indices == T
    # Filter the relevant entries
    Y_ij_upper <- Y_ij[filtering_matrix]
    
    N_ij_upper <- N_ij[filtering_matrix]
    P_NbyN_prime_upper <- P_NbyN_prime[filtering_matrix]
    
    #---------------------------------------------------------------------------
    #compute the likelihood of the the interactions between items in clusters i_star,j_star. All the rest are unaffected
    
    #First, subtract the contribution to the likelihood of those items
    A_minus = sum(dbinom(Y_ij_upper, N_ij_upper,P_NbyN_prime_upper,log = T))
    
    #Second, subtract also the contribution of the theta[i_star,j_star] entry to the prior
    C_minus = dunif(theta_ij_prime,min = lower.bound, max = upper.bound, log = T)
    
    #recompute the interaction success probabilities, just those that were affected
    
    p_ij_scanning = inverse_logit_f( theta_ij_scanning )
    
    P_NbyN_scanning[Z_i_star,Z_j_star]<- p_ij_scanning
    if(i_star != j_star){
      P_NbyN_scanning[Z_j_star,Z_i_star]<- 1 - p_ij_scanning
    }
    P_NbyN_scanning_upper = P_NbyN_scanning[filtering_matrix]
    
    #Third,add the contribution to the likelihood of those items with new interaction probabilities
    A_plus= sum(dbinom(x = Y_ij_upper,size = N_ij_upper,prob = P_NbyN_scanning_upper,log = T))
    #Fourth,add the contribution to the prior
    C_plus = dunif(theta_ij_scanning,min = lower.bound, max = upper.bound, log = T)
    
    #Updating the likelihood
    A_scanning = A_prime - A_minus + A_plus
    
    
    #Update the prior
    C_scanning = C_prime - C_minus + C_plus
    
    #evaluating the proposal density g(P'| P^(t)) 
    log_proposal_scanning <- log(dtruncnorm(theta_ij_scanning,
                                            mean = theta_ij_prime,
                                            sd = tau_theta[i_star,j_star], 
                                            a =   lower.bound, 
                                            b  = upper.bound))
    
    #evaluating the proposal density g(P^(t)| P') 
    log_proposal_prime <- log(dtruncnorm(theta_ij_prime,
                                         mean =theta_ij_scanning,
                                         sd = tau_theta[i_star,j_star], 
                                         a =  lower.bound,
                                         b  = upper.bound))
    
    
    log_r= A_scanning*t - A_prime*t + C_scanning - C_prime + log_proposal_scanning - log_proposal_prime
    
    
    
    #create statements that check conditiond to accept move
    MH_condition_P_update= min(log_r,0)>=log(runif(1))
    if(MH_condition_P_update){
      acc.count_theta[i_star,j_star] =acc.count_theta[i_star,j_star] +1
      A_prime = A_scanning
      C_prime = C_scanning
      theta_prime[i_star,j_star] <- theta_ij_scanning
      if(i_star != j_star){
        theta_prime[j_star,i_star] <- -theta_ij_scanning
      }
      P_NbyN_prime = P_NbyN_scanning
    }
    
  }
  theta_current = theta_prime
  A_current = A_prime
  
  return(list(acc.moves = acc.count_theta,
              theta = theta_current,
              llik  = A_current))
  
} 





mu_update_f = function(Y_ij, N_ij,z,theta, alpha_vec, n_k, mu_vec,tau_mu_vec, K,
                       acc.count_mu_vec,model,t, common_indices, llik = NULL, diag0.5){
  
  #computing the proportional posterior in mu' ~ g(mu^(t), tau_mu_vec)
  #mu_1_K_prime <- truncnorm::rtruncnorm(K+1,a = 0,b = 10, mean = mu_vec[1:(K+1)],sd = .05)
  # mu_0_prime = truncnorm::rtruncnorm(1,a = -Inf,b = min(mu_1_K_prime), mean = mu_vec[1],sd = .2)
  # mu_vec_prime <- c(mu_0_prime,sort(mu_1_K_prime))
  
  
  theta_prime = theta
  mu_vec_prime = mu_vec
  
  #if diag0.5 ==T, then we do not need to update the main diagonal mu, set equal to 0.5 indeed
  for(mu in (1+diag0.5):length(mu_vec) ){
    
    #if it's the last diagonal, we have no boundary above it. We fix arbitrarily 9.21
    ub = ifelse(test = mu != length(mu_vec), 
                yes = mu_vec_prime[mu+1],
                no = 9.21)
    if(mu == 1){
      ub = 9.21
    }
    lb= ifelse(test = mu == 1, 
               yes = -9.21,
               no = mu_vec_prime[mu-1])
    
    mu_k_scanning <- rtruncnorm(1,a = lb,b = ub,
                                mean = mu_vec_prime[mu],sd = tau_mu_vec[mu])
    
    theta_scanning = theta_prime
    theta_scanning[which(col(theta_prime) - row(theta_prime) == (mu-1))] <- mu_k_scanning
    theta_scanning[lower.tri(theta_scanning)] = -t(theta_scanning)[lower.tri(theta_scanning)]
    
    #computing the proportional posterior in mu'
    prop_posterior_scanning <- lprop_posterior(Y_ij = Y_ij, N_ij = N_ij, z= z, 
                                               theta = theta_scanning,
                                               alpha_vec = alpha_vec, 
                                               mu_vec = mu_vec,
                                               n_k = n_k,K = K,
                                               common_indices=common_indices,
                                               t = t)
    
    #evaluating the proportional posterior in mu^(t)
    prop_posterior_current <- lprop_posterior(Y_ij = Y_ij, N_ij = N_ij, z= z,
                                              theta = theta_prime,
                                              alpha_vec = alpha_vec,
                                              mu_vec = mu_vec,
                                              n_k = n_k,K = K,
                                              common_indices=common_indices,
                                              t = t)
    
    #evaluating the proposal density g(mu'| mu^(t)) 
    p_proposal_scanning = dtruncnorm(mu_k_scanning, 
                                     a =lb,
                                     b =ub,
                                     mean = mu_vec_prime[mu],
                                     sd = tau_mu_vec[mu])
    
    p_proposal_current = dtruncnorm(mu_vec_prime[mu],
                                    a = lb,
                                    b = ub, 
                                    mean = mu_k_scanning,
                                    sd = tau_mu_vec[mu])
    
    log_r =  prop_posterior_scanning  - prop_posterior_current + 
      log(p_proposal_current) - log(p_proposal_scanning)
    
    #create statements that check conditionion to accept move
    MH_condition_mu_vec= min(log_r,0)>=log(runif(1))
    if(MH_condition_mu_vec){
      acc.count_mu_vec[mu] = acc.count_mu_vec[mu] +1
      mu_vec_prime[mu] <- mu_k_scanning
      theta_prime <- theta_scanning
    }
  }
  mu_vec = mu_vec_prime
  theta =theta_prime
  return(list(acc.moves = acc.count_mu_vec,
              mu_vec = mu_vec,
              theta = theta))
  
} 






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

ll_naive = function(z,theta,Y_ij, N_ij,common_indices){
  P = inverse_logit_f(theta)
  z_mat = vec2mat_0_P(z,  P)
  P_ij = calculate_victory_probabilities(z_mat,  P)
  A_cur = dbinom(  Y_ij[common_indices], N_ij[common_indices], P_ij[common_indices],log = T)
  return(sum(A_cur))
}



z_update_f = function(N_ij, Y_ij, z,lamdabar,ybar,mbar, theta,
                      common_indices, alpha_vec, n_k, K, 
                      acc.count_z,labels_available,model,t){
  
  P<- inverse_logit_f(theta)
  n<- nrow(N_ij)
  z_prime= z
  P_NbyN_prime<- calculate_victory_probabilities(vec2mat_0_P(z_prime,P),P)
  
  
  A_prime<-  sum(dbinom(Y_ij[common_indices],
                        size = N_ij[common_indices], 
                        P_NbyN_prime[common_indices],log = T))
  
  B_prime<- ddirichlet_multinomial(N = n,K = K,n_k = n_k, my_alpha =  alpha_vec)
  
  
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
    k_scanning <- sample(x = setdiff(labels_available, k_prime), size = 1, replace = F)
    
    z_scanning[i_th_turn] <- k_scanning
    
    #the items in cluster i_star
    logical_matrix1 = matrix(FALSE, n,n)
    logical_matrix1[i_th_turn,]<-TRUE
    logical_matrix1[,i_th_turn]<-TRUE
    
    
    # Get the upper triangular indices for the relevant clusters
    
    
    # Create a matrix that is TRUE only at the positions that are both in the upper triangle and in the selected clusters
    filtering_matrix = (logical_matrix1)*common_indices == T
    # Filter the relevant entries
    Y_ij_upper <- Y_ij[filtering_matrix]
    N_ij_upper <- N_ij[filtering_matrix]
    P_NbyN_prime_upper <- P_NbyN_prime[filtering_matrix]
    
    #compute the likelihood of the data with the current assignment just for i_th_turn
    A_minus = sum(dbinom(Y_ij_upper, N_ij_upper, P_NbyN_prime_upper, log=T)) 
    
    #update P_NbyN
    P_NbyN_scanning = P_NbyN_prime
    for(nodes in 1:n){
      P_NbyN_scanning[i_th_turn,nodes]<- P[k_scanning,z_scanning[nodes]]
      P_NbyN_scanning[nodes,i_th_turn]<- P[z_scanning[nodes],k_scanning]
    }
    
    P_NbyN_scanning_upper <- P_NbyN_scanning[filtering_matrix]
    #compute the likelihood of the same points with the new assignment
    A_plus = sum(dbinom(Y_ij_upper, N_ij_upper, P_NbyN_scanning_upper, log=T)) 
    
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
  # ybar = t(z_P)%*%(Y_ij*upper.tri(Y_ij))%*%z_P
  # # number of missed victories between block p and block q
  # n_minus_y1 <- (N_ij-Y_ij)*upper.tri(N_ij)
  # # number of missed victories between block p and block q
  # mbar<- t(z_P)%*%n_minus_y1%*%z_P
  # 
  # coef1 = lchoose(N_ij, Y_ij)*upper.tri(N_ij)
  # lamdabar <- t(z_P)%*%(coef1)%*%z_P
  
  return(list(acc.moves = acc.count_z, z_current= z, A_prime = A_prime))
} 

