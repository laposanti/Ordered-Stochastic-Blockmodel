

P_POMM_update_given_overlap1 = function(z_current, p_current, 
                                       A_current,C_current,y_ij,n_ij,labels_available,
                                       upper.tri.non.zero,K,alpha_current,truncations_current,beta_max, overlap_current, diag0.5,acc.count_p,sigma_p){
  
  z_current_mat<- vec2mat(z_current)
  
  

  j_start = ifelse(diag0.5, yes = 1, no = 0)
  K_stop = ifelse(diag0.5, yes = K-1, no = K)
  
  
  for( ii in 1:K_stop){
    for(jj in (ii+j_start):K){
      

      p_prime = p_current
      p_prime[ii,jj] <- rnormTrunc(1, p_current[ii,jj],sd = sigma_p[ii,jj], min = 0.5, max = beta_max)
      p_prime[jj,ii] <- 1 - p_prime[ii,jj]
      
      #computing full probabilities
      p_ij_prime_nbyn <- calculate_victory_probabilities(z_current_mat, P  = p_prime)
      p_ij_prime <- p_ij_prime_nbyn[upper.tri.non.zero]
      
      C_prime <- l_like_p_ij_normal_overlap(K = K, P_matrix = p_prime,overlap = 
                                              overlap_current, 
                                            truncations = truncations_current,
                                            diag0.5 = T) + dlnorm_param(alpha_current) + dlnorm_param(overlap_current,z.mu = 0.3,z.sigma  = 0.5)         
      
      A_prime <- sum(dbinom(y_ij, n_ij, p_ij_prime, log = T))
      
      log_r= A_prime + C_prime - A_current- C_current
      
      #create statements that check conditiond to accept move
      MH_condition= min(log_r,0)>=log(runif(1))
      if(MH_condition){
        acc.count_p[ii,jj] =acc.count_p[ii,jj] +1
        acc.count_p[jj,ii] =acc.count_p[jj,ii] +1
        C_current<- C_prime
        p_current<- p_prime
        A_current <- A_prime
      }
    }}
  
  
  return(list(acc.moves = acc.count_p,
              sigma_p=sigma_p,
              p_current=p_current,
              C_current = C_current, 
              A_current = A_current))
  
} 


