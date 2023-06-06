P_matrix = synth$P_matrix
sigma0=2
non_negative_n_ij = upper.tri.non.zero

P_simple_update = function(z_mat_current, P_matrix,p_ij_current, K, n_ij,y_ij,non_negative_n_ij){
  acc.count_p=0
  p_prime = p_proposal(p = P_matrix,sigma_p = 1,K = K)
  
  
  for(i in 1:K){
    for(j in 1:K){
      p_scanning = P_matrix
      p_scanning[i,j] = p_prime[i,j]
      
      #Updating p_ij_prime with the last membership
      
      matrix_z_p_scanning  = p_scanning%*%t(z_mat_current)
      p_n_scanning = z_mat_current%*%matrix_z_p_scanning
      p_ij_scanning =  p_n_scanning[non_negative_n_ij]
      
      r = (sum(dbinom(y_ij, n_ij, p_ij_scanning, log=T)) +get_B(p_scanning, 1)) - 
        (sum(dbinom(y_ij, n_ij, p_ij_current, log = T)) +get_B(P_matrix, 1))
      
      alpha = min(1, exp(r))
      u = runif(1)
      if(u<alpha){
        #counting number of accepted proposals
        acc.count_p = acc.count_p+1
        #updating quantities
        p_ij_current = p_ij_scanning
        P_matrix = p_scanning
      }
    }
  }
  return(list(acc.count_p= acc.count_p, 
              p_current=P_matrix,
              p_ij_current=p_ij_current) )
}
  