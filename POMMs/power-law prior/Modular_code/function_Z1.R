

z_update_1 = function(z_current, A_current,B_current,y_ij,n_ij,P_matrix,labels_available,upper.tri.non.zero,gamma_vec,K){
  z_prime=z_current
  acc.value=0
  scanning_order = sample(1:N,N, replace=F)
  # full sweep
  for(ii in scanning_order){
    z_scanning=z_prime
    #save current label of z_ii
    k_cur<-z_current[ii]
    k_cur_position=which(labels_available==k_cur)
    # PROPOSE A Z_PRIME for z_prime
    # new clustering setup Z prime with only one node changed
    # reassign Z_[ii] randomly from a discrete uniform
    
    k_prime=sample(labels_available,size = 1)
    z_scanning[ii]=k_prime
    #computing p_ij_probabilities
    z_scanning_mat=vec2mat(z_scanning)
    while (TRUE) {
      if (sum(colSums(z_scanning_mat) > 0) == K) {
        break
      } else {
        ## if there is an error, undo the update
        z_scanning[ii] <- k_cur
        
        # resample new_label and try again
        k_prime <- sample(labels_available,size = 1)
        z_scanning[ii] <- k_prime
        z_scanning_mat <- vec2mat(z_scanning)
        print("!")
      }
    }
    p_ij_scanning= calculate_victory_probabilities(z_scanning_mat, P_matrix)
    p_ij_scanning = p_ij_scanning[upper.tri.non.zero]
    #compute new N and A
    # K_current unchanged 'cause we are in GS context
    n_k_scanning<-table(z_scanning)
    A_prime<- sum(dbinom(y_ij, n_ij, p_ij_scanning, log=T))
    B_prime<- ddirichlet_multinomial(N,K,n_k = n_k_scanning,my_alpha = gamma_vec)
    
    log_r= A_prime - A_current + B_prime - B_current
    #create statements that check conditiond to accept move
    GS_condition= min(log_r,0)>=log(runif(1))
    if(GS_condition){
      acc.value=acc.value+1
      z_prime<-z_scanning
      A_current<- A_prime
      B_current<- B_prime
    }
    #labels_available are the same
    #else, Z_prime[ii] stays equal to Z_current[ii]
  }
  z_current<-z_prime
  return(list(acc.moves = acc.value,A_current = A_current,B_current = B_current, z_current= z_current))
  
} 

