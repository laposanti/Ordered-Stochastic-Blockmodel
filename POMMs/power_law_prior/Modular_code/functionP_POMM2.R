

P_POMM_update2 = function(z_current, p_current, A_current,C_current,y_ij,n_ij,labels_available,upper.tri.non.zero,K,alpha_current,beta_max){
  z_current_mat<- vec2mat(z_current)
  acc.value_p<- 0
  p_current <-p_current
  
  #proposing a new alpha
  alpha_prime <- sample_norm_trunc(1, alpha_current,s =sigma_prime,a = 0.01,b = 3)
  truncations_prime <-  improper_prior5(K,beta_max,alpha = alpha_prime,diag0.5 = T)
  
  #proposing a new P
  p_prime <-  simulating_POMM_powerlaw_norm(K, alpha_prime, truncations_prime,beta_max = beta_max, diag0.5 = T)
  
  #computing full probabilities
  p_ij_prime_nbyn <- calculate_victory_probabilities(z_current_mat, P = p_prime)
  p_ij_prime <- p_ij_prime_nbyn[upper.tri.non.zero]
  
  C_prime <- l_like_p_ij_normal(K = K, P_matrix = p_prime,truncations = truncations_prime,diag0.5 = T) + dlnorm_param(alpha_prime)
  A_prime <- sum(dbinom(y_ij, n_ij, p_ij_prime, log = T))
  
  log_r= A_prime + C_prime - A_current- C_current
  
  #create statements that check conditiond to accept move
  MH_condition= min(log_r,0)>=log(runif(1))
  if(MH_condition){
    acc.value_p=acc.value_p+1
    alpha_current <- alpha_prime
    C_current<- C_prime
    p_current<- p_prime
    A_current <- A_prime
  }
  
  
  return(list(acc.moves = acc.value_p,
              p_current=p_current,
              C_current = C_current, 
              A_current = A_current,
              alpha_current = alpha_current))
  
} 


