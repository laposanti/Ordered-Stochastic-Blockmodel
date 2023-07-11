
simulating_overlapping_POMM_powerlaw_norm = function(K,  alpha = 1, overlap=1, truncations, beta_max, diag0.5=T){
  
  #we map the proportions on a desired scale
  beta_0 <- truncations
  
  j_start = ifelse(diag0.5, yes = 1, no = 0)
  K_stop = ifelse(diag0.5, yes = K-1, no = K)
  N_levelset_i = K:1
  
  P_matrix = matrix(0, K,K)
  for( i in 1:K_stop){
    for(j in (i+j_start):K){
      level_set = abs(j-i) + abs(1-j_start)
      lb <- beta_0[level_set]
      ub <- beta_0[level_set+1]
      
      # Calculate the likelihood using a truncated distribution
      mu <- (lb + ub) / 2  # Mean of the truncated distribution
      #sigma <- (ub - lb) *overlap
      sigma <- overlap*(beta_max - 0.5)/K
      P_matrix[i,j] = rtruncnorm(1,0.5,beta_max,mu,sigma)
    }
  }
  P_matrix[lower.tri(P_matrix)] = 1 - t(P_matrix)[lower.tri(P_matrix)]
  if(diag0.5){
    diag(P_matrix) = rep(0.5,K)
  }
  
  return(P_matrix)}



#--------------------------#
# joint prior on the  #
#--------------------------#
l_like_p_ij_normal_overlap = function(K, P_matrix, overlap,truncations, diag0.5 = T) {
  #we map the proportions on a desired scale
  beta_0 <- truncations
  j_start = ifelse(diag0.5, yes = 1, no = 0)
  K_stop = ifelse(diag0.5, yes = K-1, no = K)
  N_levelset_i = K:1
  
  level_sets <- diag_split_matrix(P_matrix)
  
  # Consider or exclude the main diagonal
  lowest_level_set_index <- ifelse(diag0.5, 2, 1)
  lbindex <- ifelse(diag0.5, 1, 0)
  
  log_lik_matrix = matrix(0, K,K)
  for( ii in 1:K_stop){
    for(jj in (ii+j_start):K){
      level_set = abs(jj-ii) + abs(1-j_start)
      lb <- beta_0[level_set]
      ub <- beta_0[level_set+1]
      
      # Calculate the likelihood using a truncated distribution
      mu <- (lb + ub) / 2  # Mean of the truncated distribution
      sigma <- overlap *(beta_max - 0.5)/K
      
      log_lik_matrix[ii,jj] = log(dtruncnorm(P_matrix[ii,jj], a = 0.5, b = beta_max, mean = mu, sd = sigma))
    }
  }
  return(sum(log_lik_matrix[upper.tri(log_lik_matrix)]))
}

#--------------------------#
# SIMULATION OF TOURNAMENT #
#--------------------------#

simulating_tournament_new_overlap_norm<- function(N,alpha,overlap, beta_max, K,M, n_ij_max, gamma_vec, model, diag0.5){
  #simulating K from a truncated Poisson(1)
  
  K_true =  K
  labels_available = 1:K
  #expected number of players per block
  while(TRUE){
    expected_players = round(gamma_vec/sum(gamma_vec)*N,0)
    if(sum(expected_players)==N){
      break
    }else if(sum(expected_players) < N){
      expected_players[K] = expected_players[K] +1
      break
    }else{
      expected_players[K] = expected_players[K] -1
      break}
  }
  
  # #simulating z|K from dirichlet multinomial with gamma=1
  # while(TRUE){
  #   KtimesNmatrix = rdirichlet_multinomial(N,K,gamma_vec)
  #   if(sum((rowSums(KtimesNmatrix)>0)) == K){
  #     break
  #   }else{print(paste('Number of sampled blocks =', sum((rowSums(KtimesNmatrix)>0))))}}
  # 
  
  z = vector()
  # for(i in 1:ncol(KtimesNmatrix)){
  #   z[i]= which(KtimesNmatrix[,i] >0)}
  for(i in 1:length(labels_available)){
    z= append(z, rep(i, expected_players[i]))
  }
  
  
  
  if(model == 'POMM'){
    # simulating the matrix p according to the SST property
    truncations= improper_prior5(K,beta_max,alpha = alpha,diag0.5)
    p = simulating_overlapping_POMM_powerlaw_norm(K = K,alpha = alpha,overlap = overlap,truncations = truncations,beta_max = beta_max,diag0.5)
  }else if(model == 'Simple'){
    p = matrix(runif(K*K,1-beta_max,beta_max),K,K)
    diag(p) <- rep(0.5,K)
    p=semi_symmetric(p)
  }
  
  #creating the first dataframe
  z_players <- data.frame(id = 1:N, z = z)
  z_mat_true=vec2mat(z)
  p_n_true<- calculate_victory_probabilities(z_mat_true,p)
  
  similarity_plot(p_n_true,z,z)
  
  #stratfified sampling
  block_size = colSums(z_mat_true)
  print(paste("block size",c(1:K),"=",block_size))
  
  n_ij = matrix(0, N,N)
  i=1
  while(i<=M){
    #1) sample the blocks
    block1 = sample(x = labels_available,size = 1, prob = block_size/sum(block_size))
    block2 = sample(x = labels_available,size = 1, prob = block_size/sum(block_size))
    #2) sample a random player within that block
    # Get player IDs within each block
    players_block1 <- z_players[z_players$z == block1, "id"]
    players_block2 <- z_players[z_players$z == block2, "id"]
    
    # Avoid selecting the same player from both blocks
    repeat {
      pl_1_i <- sample(players_block1, size = 1)
      pl_2_i <- sample(players_block2, size = 1)
      
      if (pl_1_i != pl_2_i)  # Check if the players are different
        break
    }
    
    # Update n_ij matrix
    n_ij[cbind(pl_1_i, pl_2_i)] <- n_ij[cbind(pl_1_i, pl_2_i)] + 1
    n_ij[cbind(pl_2_i, pl_1_i)] <- n_ij[cbind(pl_2_i, pl_1_i)] + 1
    
    if(i %%(M/10) == 0){
      print(paste((i/M*100),'-th percent process complete'))
    }
    i =i+1
  }
  
  uppertri_nonzero_entries = which(upper.tri(n_ij) & n_ij > 0,arr.ind = T)
  
  y_ij = matrix(0, N,N)
  # Select upper triangular, non-zero entries in n_ij
  n_ij_upper <- n_ij[uppertri_nonzero_entries]
  p_ij <- p_n_true[uppertri_nonzero_entries]
  y_ij[uppertri_nonzero_entries] <- rbinom(length(uppertri_nonzero_entries)/2, n_ij_upper, p_ij)
  y_ij[lower.tri(y_ij)] <- n_ij[lower.tri(n_ij)] - t(y_ij)[lower.tri(y_ij)]
  
  
  
  
  
  print(paste("total number of matches:", sum(n_ij)/2 - sum(diag(n_ij))/2))
  print(paste("total number of victories:", sum(y_ij)))
  
  return(list(n_ij_true = n_ij, y_ij_true = y_ij, p_ij_true = p_n_true, P_matrix = p, z_true= z))
}

#-----------------------#
# POMM overlap P update #
#-----------------------#



P_POMM_update_fixed_alpha_S_z = function(z_current, p_current, 
                                        A_current,C_current,y_ij,n_ij,labels_available,
                                        upper.tri.non.zero,K,alpha_current,truncations_current,beta_max, overlap_current, diag0.5,acc.count_p,sigma_p){
  
  z_current_mat<- vec2mat(z_current)
  
  
  
  j_start = ifelse(diag0.5, yes = 1, no = 0)
  K_stop = ifelse(diag0.5, yes = K-1, no = K)
  
  
  for( ii in 1:K_stop){
    for(jj in (ii+j_start):K){
      
      
      p_prime = p_current
      p_prime[ii,jj] <- rnormTrunc(1, p_current[ii,jj],sd = sigma_p[ii,jj], min = 0.5, max = beta_max)
      p_prime[jj,ii] <- 1 - p_prime[ii,jj]
      
      #computing full probabilities
      p_ij_prime_nbyn <- calculate_victory_probabilities(z_current_mat, P  = p_prime)
      p_ij_prime <- p_ij_prime_nbyn[upper.tri.non.zero]
      
      C_prime <- l_like_p_ij_normal_overlap(K = K, P_matrix = p_prime,overlap = 
                                              overlap_current, 
                                            truncations = truncations_current,
                                            diag0.5 = T) + + dunif(alpha_current,0,4,log = T) + dunif(overlap_current,min = 0,max=1, log = T)
      
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






#------------------------#
# POMM overlap update (S)  #
#------------------------#




S_POMM_update_fixed_P_alpha_z = function(z_current, p_current, 
                                 A_current,C_current,y_ij,n_ij,labels_available,
                                 upper.tri.non.zero,K,truncations_current, alpha_current,beta_max, overlap_current,acc.count_overlap, sigma_overlap){
  
  
  overlap_prime <- rtruncnorm(1,a = 0.1,b = 0.9,mean = overlap_current,sd = sigma_overlap)
  
  
  
  C_prime <- l_like_p_ij_normal_overlap(K = K, P_matrix = p_current,overlap = 
                                          overlap_prime, 
                                        truncations = truncations_current,
                                        diag0.5 = T) + dunif(alpha_current,0,4,log = T) + dunif(overlap_prime,min = 0,max=1, log = T)
  
  
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





#---------
#alpha update
#---------



alpha_POMM_update_fixed_P_S_z = function(z_current, p_current, 
                               A_current,C_current,y_ij,n_ij,labels_available,
                               upper.tri.non.zero,K,alpha_current,beta_max, overlap_current, acc.count_alpha, sigma_alpha, truncations_current){
  
  z_current_mat<- vec2mat(z_current)
  
  #proposing a new alpha
  
  alpha_prime <- rtruncnorm(1,a = 0.1,b = 0.9,mean = alpha_current,sd = sigma_alpha)
  truncations_prime <- improper_prior5(K,beta_max,alpha = alpha_prime,diag0.5 = T)
  
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
                                          overlap_current, 
                                        truncations = truncations_prime,
                                        diag0.5 = T) + dunif(alpha_prime,0,4,log = T) + dunif(overlap_current,min = 0,max=1, log = T)
  
  
  
  
  #A_prime <- sum(dbinom(y_ij, n_ij, p_ij_prime, log = T))
  
  log_r= C_prime - C_current
  
  #create statements that check conditiond to accept move
  MH_condition_alpha= min(log_r,0)>=log(runif(1))
  if(MH_condition_alpha){
    acc.count_alpha= acc.count_alpha+1
    alpha_current <- alpha_prime
    truncations_current<- truncations_prime
    C_current<- C_prime
    #A_current <- A_prime
  }
  
  
  return(list(acc.moves = acc.count_alpha,
              sigma_alpha=sigma_alpha,
              C_current = C_current, 
              alpha_current = alpha_current,
              truncations_current=truncations_current))
  
} 





#-----------------
#z adaptive update
#-----------------




z_update_adaptive = function(z_current, A_current,B_current,y_ij,n_ij,P_matrix,labels_available,upper.tri.non.zero,gamma_vec,K, acc.count_z, sigma_z){
  z_prime=z_current
  scanning_order = sample(1:N,N, replace=F)
  # full sweep
  for(ii in scanning_order){
    z_scanning=z_prime
    #save current label of z_ii
    k_cur<-z_current[ii]
    
    # PROPOSE A Z_PRIME for z_prime
    
    # Calculate distances between labels
    label_distances <- abs(labels_available - k_cur)
    
    # Calculate probabilities based on label distances and sigma.z
    probabilities <- label_probability(label_distances, sigma_z[ii])
    
    # Adjust probabilities based on sigma.z
    probabilities <- probabilities[-k_cur] / sum(probabilities[-k_cur])
    
    # # Sample a new label using the adjusted probabilities
    k_prime <- sample(x = setdiff(labels_available, k_cur), size = 1, replace = F)
    
    z_scanning[ii] <- k_prime
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
      acc.count_z[ii]=acc.count_z[ii]+1
      z_prime<-z_scanning
      A_current<- A_prime
      B_current<- B_prime
    }
    #labels_available are the same
    #else, Z_prime[ii] stays equal to Z_current[ii]
  }
  z_current<-z_prime
  
  return(list(acc.moves = acc.count_z,A_current = A_current,B_current = B_current, z_current= z_current))
} 























#------------------
# Adaptive Proposal
#------------------


tuning_proposal<- function(iteration, acceptance_count, sigma, acceptanceTarget, min_sigma){
  acceptanceRate <- acceptance_count / iteration
  delta= min(0.01, iteration**(-1/2))
  lsi= log(sigma)
  if (acceptanceRate > acceptanceTarget) {
    lsi <- lsi + delta
  } else {
    lsi <- lsi - delta}
  
  # Update the proposal standard deviations
  sigma_updated <- max(min_sigma,exp(lsi))
  return(sigma_updated)
}



label_probability <- function(distance, sigma_z) {
  prob <- dnorm(distance, mean = 0, sd = sigma_z)
  return(prob)
}

# calculate_victory_probabilities_modified <- function(z_mat, P) {
#   # Convert z_mat to a sparse matrix if it's mostly zeros
#   if (sum(z_mat != 0) < 0.2 * length(z_mat)) {
#     z_mat <- sparseMatrix(i = row(z_mat), j = col(z_mat), x = z_mat)
#   }
#   
#   # Calculate intermediate result and cache it
#   aux <- P %*% t(z_mat)
#   
#   # Perform matrix multiplication using cached intermediate result
#   p_ij_scanning <- z_mat %*% aux
#   
#   return(p_ij_scanning)
# }
