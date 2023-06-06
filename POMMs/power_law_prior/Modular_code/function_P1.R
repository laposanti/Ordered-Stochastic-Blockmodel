
P_simple_update1 = function(z_current, P_matrix, K, n_ij,y_ij,A_current, C_current, upper.tri.non.zero){
  z_mat_current <- vec2mat(z_current)
  acc.count_p=0
  for(i in 1:K){
    for(j in 1:K){
      p_prime = p_proposal(p = P_matrix,sigma_p = .025,K = K)
      p_scanning = P_matrix
      p_scanning[i,j] = p_prime[i,j]
      p_scanning[j,i] = 1 - p_prime[i,j]
      #Updating p_ij_prime with the last membership
      
      matrix_z_p_scanning  = p_scanning%*%t(z_mat_current)
      p_n_scanning = z_mat_current%*%matrix_z_p_scanning
      p_ij_scanning =  p_n_scanning[upper.tri.non.zero]
      
      A_prime = sum(dbinom(y_ij, n_ij, p_ij_scanning, log=T))
      C_prime = get_B(p_scanning, 1)
                    
      r = A_prime + C_prime - A_current - C_current
      
      alpha = min(1, exp(r))
      u = runif(1)
      if(u<alpha){
        A_current <- A_prime
        C_current <- C_prime 
        #counting number of accepted proposals
        acc.count_p <- acc.count_p+1
        #updating quantities
        P_matrix = p_scanning
      }
    }
  }
  return(list(acc.count_p= acc.count_p, 
              p_current=P_matrix,
              A_current=A_current,
              C_current=C_current))
}
