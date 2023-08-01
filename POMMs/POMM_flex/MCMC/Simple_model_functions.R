
P_simple_update = function(z_current, p_current, 
                                         A_current,C_current,y_ij,n_ij,labels_available,
                                         upper.tri.non.zero,K, diag0.5,acc.count_p,sigma_p){
  
  A_prime <- A_current
  C_prime <- sum(dbeta(p_current[upper.tri(p_current)],1,1,log=T))
  z_prime <- z_current
  p_prime <- p_current
  P_NbyN_prime <- calculate_victory_probabilities(vec2mat(z_prime),p_prime)
  
  j_start <- ifelse(diag0.5, yes = 1, no = 0)
  K_stop <- ifelse(diag0.5, yes = K-1, no = K)
  
  for(p_i in 1:K_stop){
    for(p_j in (p_i+j_start):K){
      
      #extracing just players in the updating clusters
      players_ii <- which(z_prime==p_i)
      players_jj <- which(z_prime==p_j)
      
      #the current likelihood and prior for those players
      A_minus = sum(dbinom(y_ij[players_ii,players_jj], n_ij[players_ii,players_jj], P_NbyN_prime[players_ii,players_jj], log=T)) + sum(dbinom(y_ij[players_jj,players_ii], n_ij[players_jj,players_ii], P_NbyN_prime[players_jj,players_ii], log=T))
      C_minus <- dbeta(p_prime[p_i,p_j],1,1,log = T)
      
      #proposing a new p_ij
      p_scanning = p_prime
      p_scanning[p_i,p_j] <- rtruncnorm(1, mean = p_prime[p_i,p_j],sd = sigma_p[p_i,p_j], a =   0.5, b  = beta_max)
      p_scanning[p_j,p_i] <- 1 - p_scanning[p_i,p_j]
      
      #updating P_NbyN_prime for all players in cluster ii and cluster jj
      P_NbyN_scanning = P_NbyN_prime
      P_NbyN_scanning[players_ii,players_jj] <- p_scanning[p_i,p_j]
      P_NbyN_scanning[players_jj,players_ii] <- p_scanning[p_j,p_i]
      
      #the new likelihood and prior values for those players
      A_plus = sum(dbinom(y_ij[players_ii,players_jj], n_ij[players_ii,players_jj], P_NbyN_scanning[players_ii,players_jj], log=T)) + sum(dbinom(y_ij[players_jj,players_ii], n_ij[players_jj,players_ii], P_NbyN_scanning[players_jj,players_ii], log=T))
      C_plus <- dbeta(p_scanning[p_i,p_j],1,1,log = T)
      
      A_scanning = A_prime - A_minus + A_plus
      C_scanning = C_prime - C_minus + C_plus
      
      
      #browser()
      log_r= A_scanning - A_prime + C_scanning - C_prime  
      
      #create statements that check conditiond to accept move
      MH_condition= min(log_r,0)>=log(runif(1))
      if(MH_condition){
        acc.count_p[p_i,p_j] =acc.count_p[p_i,p_j] +1
        acc.count_p[p_j,p_i] =acc.count_p[p_j,p_i] +1
        C_prime<- C_scanning
        p_prime<- p_scanning
        A_prime <- A_scanning
        P_NbyN_prime<-P_NbyN_scanning
      }
    }}
  
  A_current <- A_prime
  C_current <-C_prime
  z_current <-z_prime
  p_current <-p_prime
  
  return(list(acc.moves = acc.count_p,
              sigma_p=sigma_p,
              p_current=p_current,
              C_current = C_current, 
              A_current = A_current))
  
} 
