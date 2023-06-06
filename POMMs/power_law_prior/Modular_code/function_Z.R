

#updating z function
z_update=function(z_current,p_nbyn_current, n_k_current, P_matrix,labels_available, K,gamma_vec,N, n_ij,y_ij,upper.tri.non.zero,aux){
  #N total number of players
  #n_ij_matrix contains the number of games between player i and j
  #y_ij_matrix contains the number of victories of player i vs j
  #upper.tri.non.zero= which(n_ij_matrix > 0 & upper.tri(n_ij_matrix))
  # n_ij = n_ij_matrix[upper.tri.non.zero]
  # y_ij = y_ij_matrix[upper.tri.non.zero]
  # n_k_current is the number of players in cluster k for every k=1,...K
  # p_nbyn_current contains the probabilities of player i winning against player j
  p_nbyn_prime=p_nbyn_current
  accepted=0

  #Complete sweeep of z vector
  #----
  sweeping_order = sample(x=1:N,size=N, replace=F)
  for(i in sweeping_order){
    #retrieving current label
    old_label = z_current[i]
    #proposing a new label
    new_label = adjacent_label_sampler(labels_available = labels_available, old_label)
    
    z_mat_prime <- z_mat_current
    
    # Remove the current entry value and updating
    z_mat_prime[i,old_label] <- 0 
    z_mat_prime[i, new_label] <- 1
    
    
    
    #updating labels
    # z_prime = z_current
    # z_prime[i] = new_label
    # z_mat_prime= vec2mat(z_prime)
    
    #computing the new victory probabilities
    while(TRUE){
      if(sum(colSums(z_mat_prime) >0)==K){
        break
      } else {
        ##if there is an error, undo the update
        z_mat_prime[i,new_label] <- 0
        #resample new_label and try again
        new_label = adjacent_label_sampler(labels_available=labels_available, z_current[i])
        z_mat_prime[i, new_label] <- 1
        print(sum(colSums(z_mat_prime)))
      }
    }
    
    z_prime=z_current
    z_prime[i] = new_label
    
    
    # Update the victory_probs matrix incrementally
    # z_mat_n <- z_mat_current[i, ]
    # vict_n_vs_all <- z_mat_n %*% aux
    # 
    #no need to update also the victory of j vs i, it will be disregarded
    # aux1 <- z_mat_current %*% P_matrix
    # vict_all_vs_n <- aux1 %*% z_mat_n
    # 
    # p_nbyn_prime[i, ] <- vict_n_vs_all
    # p_nbyn_prime[, i] <- vict_all_vs_n
    # 
    matrix_z_prime_p_current = P_matrix%*%t(z_mat_prime)
    p_n_z_prime = z_mat_prime%*%matrix_z_prime_p_current
    semi_similarity_checker(z_prime, p_n_z_prime)
    p_ij_z_prime =  p_n_z_prime[upper.tri.non.zero]
    
    n_k_prime = n_k_current
    n_k_prime[z_current[i]] = n_k_prime[z_current[i]] - 1
    n_k_prime[new_label] = n_k_prime[new_label] + 1
    
    #acceptance ratio
    r = (sum(dbinom(y_ij, n_ij, p_ij_z_prime, log = T)) + ddirichlet_multinomial(N,K,n_k = n_k_prime,my_alpha = gamma_vec)) - 
      (sum(dbinom(y_ij, n_ij, p_nbyn_current[upper.tri.non.zero], log = T)) + ddirichlet_multinomial(N,K,n_k = n_k_current ,my_alpha = gamma_vec))
    print(r)
    # r =ddirichlet_multinomial(N,K,n_k = n_k_prime,my_alpha = gamma_vec)-ddirichlet_multinomial(N,K,n_k = n_k_current ,my_alpha = gamma_vec)
    alpha_r = min(1, exp(r))
    u = runif(1)
    
    #if accepted
    if(u<alpha_r){
      accepted= accepted +1
      z_mat_current = z_mat_prime
      z_current[i]=new_label
      n_k_current = n_k_prime
      p_nbyn_current= p_ij_z_prime
      aux <- P_matrix %*% t(z_mat_current)
      #if not accepted
    }
  }
  return(list(acc = accepted,
              z_mat_cur=z_mat_current,
              z_cur= z_current,
              n_cur=n_k_current,
              p_nbyn_cur= p_nbyn_current,
              aux=aux))
}


