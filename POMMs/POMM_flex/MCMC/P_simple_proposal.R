

P_simple_update_adaptive = function(z_current, p_current, 
                                        A_current,C_current,y_ij,n_ij,labels_available,
                                        upper.tri.non.zero,K, diag0.5,acc.count_p,sigma_p,beta_max){
  
  z_current_mat<- vec2mat(z_current)
  
  

  j_start = ifelse(diag0.5, yes = 1, no = 0)
  K_stop = ifelse(diag0.5, yes = K-1, no = K)
  
  
  for( ii in 1:K_stop){
    for(jj in (ii+j_start):K){
      
      
      p_prime = p_current
      p_prime[ii,jj] <- rnormTrunc(1, p_current[ii,jj],sd = sigma_p[ii,jj], min = 0.01, max = .99)
      p_prime[jj,ii] <- 1 - p_prime[ii,jj]
      
      #computing full probabilities
      p_ij_prime_nbyn <- calculate_victory_probabilities(z_current_mat, P  = p_prime)
      p_ij_prime <- p_ij_prime_nbyn[upper.tri.non.zero]
      
      C_prime <- sum(dunif(p_prime[upper.tri(p_prime)],min = 1-beta_max, max = beta_max, log = T))
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


