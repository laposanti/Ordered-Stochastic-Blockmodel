
P_POMM_update = function(z_current, P_matrix, K, n_ij,y_ij,A_current, C_current, upper.tri.non.zero, alpha_current, truncations_current, beta_max){
  z_mat_current <- vec2mat(z_current)
  acc.count_p=0
  
  
  sampled_order = sample(0:(K-1),K,replace = F)
  sampled_i =  sample(1:K,K,replace = F)
  sampled_j =  sample(1:K,K,replace = F)
  for(i in sampled_i){
    for(j in sampled_j){
      level_set <- abs(j-i)
      #proposing a new alpha
      alpha_prime <- sample_norm_trunc(1, alpha_current,s =sigma_prime,a = 0.01,b = 3)
      truncations_prime <- improper_prior5(K,beta_max,alpha = alpha_prime)
      
      #generating a proposal matrix
      p_prime = simulating_POMM_powerlaw1(K,alpha_prime,truncations_prime,beta_max)
      
      p_scanning = P_matrix
      p_scanning[i,j] = p_prime$matrix[i,j]
      p_scanning[j,i] = 1 - p_prime$matrix[i,j]
      diag(p_scanning) <- rep(0.5,K)
      #Updating p_ij_prime with the last membership
      
      matrix_z_p_scanning  = p_scanning%*%t(z_mat_current)
      p_n_scanning = z_mat_current%*%matrix_z_p_scanning
      p_ij_scanning =  p_n_scanning[upper.tri.non.zero]
      
      truncations_scanning = truncations_current
      truncations_scanning[c(level_set+1,level_set + 2)] = truncations_prime[c(level_set+1,level_set + 2)]
      
      A_prime = sum(dbinom(y_ij, n_ij, p_ij_scanning, log=T))
      C_prime =  l_like_p_ij(p_scanning,truncations_scanning) + dlnorm_param(alpha_prime)
      
      r = A_prime + C_prime - A_current - C_current
      
      alpha = min(1, exp(r))
      u = runif(1)
      if(u<alpha){
        #counting number of accepted proposals
        acc.count_p <- acc.count_p+1
        
        A_current <- A_prime
        C_current <- C_prime 
        
        #updating quantities
        P_matrix = p_scanning
        alpha_current <- alpha_prime
        truncations_current <- truncations_scanning
      }
    }
  }
  return(list(acc.count_p= acc.count_p, 
              alpha_current = alpha_current,
              p_current= P_matrix,
              A_current= A_current,
              C_current= C_current,
              truncations_current = truncations_current))
}











