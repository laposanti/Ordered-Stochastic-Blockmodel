
alpha_POMM_update1 = function(z_current, A_current,C_current,y_ij,n_ij,p_current,labels_available,upper.tri.non.zero,K,alpha_current, truncations_current,beta_max){
  z_current_mat=vec2mat(z_current)
  acc.value_p=0
  
  alpha_prime = sample_norm_trunc(1, alpha_current,s =sigma_prime,a = 0.01,b = 3)
  truncations_prime =  improper_prior5(K,beta_max,alpha = alpha_prime,diag0.5 = T)
  
  
  C_prime<- l_like_p_ij_normal(K = K, P_matrix = p_current,truncations = truncations_prime,diag0.5 = T) + dlnorm_param(alpha_prime)
  
  log_r= C_prime  - C_current
  
  #create statements that check conditiond to accept move
  GS_condition= min(log_r,0)>=log(runif(1))
  if(GS_condition){
    acc.value_p=acc.value_p+1
    truncations_current = truncations_prime
    alpha_current <- alpha_prime
    C_current<- C_prime
  }


return(list(acc.moves = acc.value_p,
            C_current = C_current, 
            alpha_current = alpha_current))

} 


