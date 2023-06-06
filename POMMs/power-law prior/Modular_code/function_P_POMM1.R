

P_POMM_update1 = function(z_current, A_current,C_current,y_ij,n_ij,p_current,labels_available,upper.tri.non.zero,K,alpha_current, truncations_current,beta_max){
  z_current_mat=vec2mat(z_current)
  acc.value_p=0
 
 
  # full sweep
  for(i in 1:(K-1)){
    for(j in (i+1):K){
  
      alpha_prime = sample_norm_trunc(1, alpha_current,s =sigma_prime,a = 0.01,b = 3)
      truncations_prime =  improper_prior5(K,beta_max,alpha = alpha_prime,diag0.5 = T)
      #p_prime = simulating_POMM_powerlaw2(K = K,alpha = alpha_prime,truncations = truncations_current,beta_max = beta_max,diag0.5 = T)
      
      
      truncations_scanning = truncations_current
      #p_scanning = p_current
      level_set <- abs(j-i)
      alpha_scanning = alpha_prime
      #updating truncations
      
      truncations_scanning[c(level_set,level_set + 1)] = truncations_prime[c(level_set,level_set + 1)]
      #updating p_scanning
      
      #p_scanning[i,j] <-  p_prime[i,j]
      #p_scanning[j,i] <- p_prime[j,i]
      
      

      #computing p_ij_probabilities
      #p_ij_scanning= calculate_victory_probabilities(z_current_mat, p_current)
      #p_ij_scanning = p_ij_scanning[upper.tri.non.zero]
      
      #A_prime<- sum(dbinom(y_ij, n_ij, p_ij_scanning, log=T))
      C_prime<- l_like_p_ij_normal(K = K, P_matrix = p_current,truncations = truncations_scanning,diag0.5 = T) + dlnorm_param(alpha_scanning)
      
      log_r= C_prime  - C_current

      #create statements that check conditiond to accept move
      GS_condition= min(log_r,0)>=log(runif(1))
      if(GS_condition){
        acc.value_p=acc.value_p+1
        #p_prime <-p_scanning
        truncations_prime = truncations_scanning
        alpha_current <- alpha_prime
        #A_current<- A_prime
        C_current<- C_prime
      }
    }
  }
  #labels_available are the same
  #else, Z_prime[ii] stays equal to Z_current[ii]

#p_current<-p_prime
alpha_current <- alpha_prime
truncations_current<-truncations_prime
return(list(acc.moves = acc.value_p,
            #A_current = A_current,
            C_current = C_current, 
            p_current = p_current,
            alpha_current = alpha_prime,
            truncations_current = truncations_prime))

} 

