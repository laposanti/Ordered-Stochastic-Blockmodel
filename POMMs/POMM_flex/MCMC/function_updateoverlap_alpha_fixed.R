

P_POMM_overlap_update = function(z_current, p_current, 
                                 A_current,C_current,y_ij,n_ij,labels_available,
                                 upper.tri.non.zero,K,truncations_current, alpha_current,beta_max, overlap_current,acc.count_overlap, sigma_overlap){
  
  z_current_mat<- vec2mat(z_current)
  

  overlap_prime <- rtruncnorm(1,a = 0.1,b = 0.9,mean = overlap_current,sd = sigma_overlap)
  

  # #proposing a new P
  # p_prime <-  simulating_overlapping_POMM_powerlaw_norm(K, alpha_current, 
  #                                                       truncations_current, 
  #                                                       overlap = overlap_prime, 
  #                                                       beta_max = beta_max, diag0.5 = T)
  # 
  # #computing full probabilities
  # p_ij_prime_nbyn <- calculate_victory_probabilities(z_current_mat, P = p_prime)
  # p_ij_prime <- p_ij_prime_nbyn[upper.tri.non.zero]

  
  C_prime <- l_like_p_ij_normal_overlap(K = K, P_matrix = p_current,overlap = 
                                          overlap_prime, 
                                        truncations = truncations_current,
                                        diag0.5 = T) + dlnorm_param(alpha_current) + dlnormTruncAlt(x = overlap_prime,mean =0.5, cv = 1, min = 0.1,max=0.9)       


  #A_prime <- sum(dbinom(y_ij, n_ij, p_ij_prime, log = T))
  
  log_r=  C_prime -  C_current
  
  #create statements that check conditiond to accept move
  MH_condition_overlap= min(log_r,0)>=log(runif(1))
  if(MH_condition_overlap){
    acc.count_overlap=acc.count_overlap+1
    overlap_current <- overlap_prime
    C_current<- C_prime
    p_current<- p_current
    #A_current <- A_prime
  }
  
  
  return(list(acc.moves = acc.count_overlap,
              sigma_overlap=sigma_overlap,
              p_current=p_current,
              C_current = C_current, 
              A_current = A_current,
              overlap_current = overlap_current))
  
} 


