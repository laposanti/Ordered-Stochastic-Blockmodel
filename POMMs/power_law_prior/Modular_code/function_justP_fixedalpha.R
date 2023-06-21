

#fix alpha and get P: function


P_POMM_update2_fixed_alpha = function(z_current, p_current, A_current,C_current,y_ij,n_ij,labels_available,upper.tri.non.zero,K,alpha_current,beta_max){
  z_current_mat<- vec2mat(z_current)
  acc.value_p<- 0
  p_current <-p_current
  

  truncations_current<-  improper_prior5(K,beta_max,alpha = alpha_current,diag0.5 = T)
  
  #proposing a new P
  p_prime <-  simulating_POMM_powerlaw_norm(K, alpha_current, truncations_current,beta_max = beta_max, diag0.5 = T)
  
  #computing full probabilities
  p_ij_prime_nbyn <- calculate_victory_probabilities(z_current_mat, P = p_prime)
  p_ij_prime <- p_ij_prime_nbyn[upper.tri.non.zero]
  
  C_prime <- l_like_p_ij_normal(K = K, P_matrix = p_prime,truncations = truncations_current,diag0.5 = T) 
  A_prime <- sum(dbinom(y_ij, n_ij, p_ij_prime, log = T))
  
  log_r=  C_prime  - C_current
  
  #create statements that check conditiond to accept move
  MH_condition= min(log_r,0)>=log(runif(1))
  if(MH_condition){
    acc.value_p=acc.value_p+1
    C_current<- C_prime
    p_current<- p_prime
    A_current <- A_prime
  }
  
  
  return(list(acc.moves = acc.value_p,
              p_current=p_current,
              C_current = C_current, 
              A_current = A_current))
  
} 


