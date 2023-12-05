


z_update_adaptive = function(z_current, A_current,B_current,y_ij,n_ij,P_matrix,labels_available,upper.tri.non.zero,gamma_vec,K, acc.count_z, sigma_z){
  
  A_prime<- A_current
  B_prime<- B_current
  z_prime=z_current
  P_NbyN_prime<- calculate_victory_probabilities(vec2mat_0_P(z_prime,P_matrix),P_matrix)
  n_prime = table(z_prime)
  
  scanning_order = sample(1:N,N, replace=F)
  # full sweep
  for(ii in scanning_order){
    
    z_scanning=z_prime
    #save current label of z_ii
    k_prime<-z_prime[ii]
    
    # PROPOSE A Z_PRIME for z_prime
    
    # Calculate distances between labels
    label_distances <- abs(labels_available - k_cur)
    
    # Calculate probabilities based on label distances and sigma.z
    probabilities <- label_probability(label_distances, sigma_z[ii])
    
    # Adjust probabilities based on sigma.z
    probabilities <- probabilities[-k_cur] / sum(probabilities[-k_cur])
    
    # Sample a new label using the adjusted probabilities
    k_scanning <- sample(x = setdiff(labels_available, k_cur), size = 1, replace = F)
    
    z_scanning[ii] <- k_scanning
    
    while (TRUE) {
      if (length(unique(z_scanning)) == K) {
        print('ok')
        break
      } else {
        ## if there is an error, undo the update
        z_scanning[ii] <- k_prime
        
        # resample new_label and try again
        k_prime <- sample(labels_available,size = 1)
        z_scanning[ii] <- k_scanning
        
        print("!")
      }
    }
    
    #compute the likelihood of the data with the current assignment just for ii
    A_minus = sum(dbinom(y_ij[ii,], n_ij[ii,], P_NbyN_prime[ii,], log=T)) + sum(dbinom(y_ij[,ii], n_ij[,ii], P_NbyN_prime[,ii], log=T)) 
    
    #update P_NbyN
    P_NbyN_scanning = P_NbyN_prime
    for(nodes in 1:N){
      P_NbyN_scanning[ii,nodes]<- P_matrix[k_scanning,z_scanning[nodes]]
      P_NbyN_scanning[nodes,ii]<- P_matrix[z_scanning[nodes],k_scanning]
    }
    #compute the p_z'i,zj of the i-th row and the i-th column with the new assignment
    A_plus = sum(dbinom(y_ij[ii,], n_ij[ii,], P_NbyN_scanning[ii,], log=T)) + sum(dbinom(y_ij[,ii], n_ij[,ii], P_NbyN_scanning[,ii], log=T)) 
    
    #Updating the likelihood for the proposal
    A_scanning = A_prime - A_minus + A_plus
    
    
    n_k_scanning <- table(z_scanning)
    #n_k_scanning[c(k_cur, k_prime)] <- n_prime[c(k_cur, k_prime)] + c(-1, 1)
    
    B_scanning<- ddirichlet_multinomial(N,K,n_k = n_k_scanning,my_alpha = gamma_vec)
    
    log_r= A_scanning - A_prime + B_scanning - B_prime
    #create statements that check conditiond to accept move
    GS_condition= min(log_r,0)>=log(runif(1))
    if(GS_condition){
      acc.count_z[ii]=acc.count_z[ii]+1
      z_prime<-z_scanning
      A_prime<- A_scanning
      B_prime<- B_scanning
      P_NbyN_prime <- P_NbyN_scanning
      n_prime <- n_k_scanning
    }
    #labels_available are the same
    #else, z_prime[ii] stays equal to z_current[ii]
  }
  z_current<-z_prime
  A_current<- A_prime
  B_current<- B_prime
  
  
  
  return(list(acc.moves = acc.count_z,A_current = A_current,B_current = B_current, z_current= z_current))
} 

