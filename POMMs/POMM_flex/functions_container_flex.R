expected_min<- function(mu, sigma,n){
  return(mu + sigma*qnorm((1- pi/8)/(n-pi/4+1)))}
expected_max<- function(mu,sigma,n){
  return(mu + sigma*qnorm((n- pi/8)/(n-pi/4+1)))}

# Function to compute the parameter overlap given sigma
compute_S <- function(sigma, truncations,K,n_samples=1) {
  
  mid_points = vector()
  for(i in 1:(length(truncations))-1){
    mid_points <- append(mid_points, (truncations[i]+truncations[i+1])/2)
  }
  
  dim_level_set_i <- K:1
  expected_boundaries = matrix(0, nrow = length(mid_points), ncol=2)
  for(i in 1:length(mid_points)){
    expected_boundaries[i,1]<- expected_min(mu = mid_points[i],sigma = sigma,n = (dim_level_set_i[i])*n_samples)
    expected_boundaries[i,2]<- expected_max(mu = mid_points[i],sigma = sigma,n = (dim_level_set_i[i])*n_samples)
  }
  
  S <- 0
  possible_S <- 0
  num_expected_boundaries <- nrow(expected_boundaries)
  for (i in 1:(num_expected_boundaries - 1)) {
    for (j in (i + 1):num_expected_boundaries) {
      possible_S <- possible_S + max(expected_boundaries[j, 2],expected_boundaries[i, 2]) - min(expected_boundaries[i, 1], expected_boundaries[j, 1])
      if (expected_boundaries[i, 2] > expected_boundaries[j, 1] && expected_boundaries[j, 2] > expected_boundaries[i, 1]) {
        S <- S + (min(expected_boundaries[i, 2], expected_boundaries[j, 2]) - max(expected_boundaries[i, 1], expected_boundaries[j, 1]))
      }
    }
  }
  return(S / possible_S)
}

# Function to find the inverse of the S function
inverse_S <- function(S,truncations,K, beta_max,n_samples=1) {
  # Define the range of sigma values
  sigma_min <- 0  # Minimum possible sigma value
  sigma_max <- beta_max - 0.5  # Maximum possible sigma value
  
  # Define the desired precision
  epsilon <- 1e-6
  
  # Perform bisection method
  while (sigma_max - sigma_min > epsilon) {
    sigma_mid <- (sigma_min + sigma_max) / 2
    S_mid <- compute_S(sigma_mid, truncations,K,n_samples)
    
    if (S_mid < S) {
      sigma_min <- sigma_mid
    } else {
      sigma_max <- sigma_mid
    }
  }
  
  return(sigma_min)  # or sigma_max, they should be close at this point
}

simulating_overlapping_POMM_powerlaw_norm = function(K,  alpha = 1, S=1, truncations, beta_max, diag0.5=T, phi=0){
  
  #we map the proportions on a desired scale
  
  j_start = ifelse(diag0.5, yes = 1, no = 0)
  K_stop = ifelse(diag0.5, yes = K-1, no = K)
  N_levelset_i = K:1
  
  P_matrix = matrix(0, K,K)
  for( i in 1:K_stop){
    for(j in (i+j_start):K){
      level_set = abs(j-i) + abs(1-j_start)
      lb <- truncations[level_set]
      ub <- truncations[level_set+1]
      
      # Calculate the likelihood using a truncated distribution
      mu <- (lb + ub) / 2  # Mean of the truncated distribution
      #sigma <- (ub - lb) *overlap
      sigma <- S*((1-phi)+ phi*(j+i))
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
l_like_p_ij_normal_overlap <-function(K, P_matrix,S ,truncations, diag0.5 = T) {
  #we map the proportions on a desired scale
  beta_0 <- truncations
  j_start = ifelse(diag0.5, yes = 1, no = 0)
  K_stop = ifelse(diag0.5, yes = K-1, no = K)
  N_levelset_i = K:1
  
  level_sets <- diag_split_matrix(P_matrix)
  
  # Consider or exclude the main diagonal

  log_lik_matrix = matrix(0, K,K)
  sigma <- S
  for( ii in 1:K_stop){
    for(jj in (ii+j_start):K){
      
      
      level_set = abs(jj-ii) + abs(1-j_start)
      lb <- beta_0[level_set]
      ub <- beta_0[level_set+1]
      
      # Calculate the likelihood using a truncated distribution
      mu <- (lb + ub) / 2  # Mean of the truncated distribution
      
      log_lik_matrix[ii,jj] = log(dtruncnorm(P_matrix[ii,jj], a = 0.5, b = beta_max, mean = mu, sd = sigma))
    }
  }
  return(sum(log_lik_matrix[upper.tri(log_lik_matrix)]))
}

single_p_ij_normal_overlap = function(entry_i,entry_j, K, P_matrix,S ,truncations, diag0.5 = T) {
  #we map the proportions on a desired scale
  beta_0 <- truncations
  j_start = ifelse(diag0.5, yes = 1, no = 0)
  
  # Consider or exclude the main diagonal
  
  sigma <- S
  
  level_set = abs(entry_j -entry_i) + abs(1-j_start)
  lb <- beta_0[level_set]
  ub <- beta_0[level_set+1]
  
  # Calculate the likelihood using a truncated distribution
  mu <- (lb + ub) / 2  # Mean of the truncated distribution
  
  log_p = log(dtruncnorm(P_matrix[entry_i,entry_j], a = 0.5, b = beta_max, mean = mu, sd = sigma))

return(log_p)
}

#--------------------------#
# SIMULATION OF TOURNAMENT #
#--------------------------#

simulating_tournament_new_overlap_norm<- function(N,alpha,S, beta_max, K,M, n_ij_max, gamma_vec, model, diag0.5){
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
    p = simulating_overlapping_POMM_powerlaw_norm(K = K,alpha = alpha,S = S,truncations = truncations,beta_max = beta_max,diag0.5)
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
                                         upper.tri.non.zero,K,alpha_current,truncations_current,beta_max, S_current, diag0.5,acc.count_p,sigma_p){
  
  A_prime <- A_current
  C_prime <- l_like_p_ij_normal_overlap(K,p_current,S_current,truncations_current,diag0.5)
  z_prime <- z_current
  p_prime <- p_current
  P_NbyN_prime <- calculate_victory_probabilities(vec2mat_0_P(z_prime,p_prime),p_prime)
  
  j_start <- ifelse(diag0.5, yes = 1, no = 0)
  K_stop <- ifelse(diag0.5, yes = K-1, no = K)
  
  for(p_i in 1:K_stop){
    for(p_j in (p_i+j_start):K){
      
      #extracing just players in the updating clusters
      players_ii <- which(z_current==p_i)
      players_jj <- which(z_current==p_j)
      
      #the current likelihood and prior for those players
      A_minus = sum(dbinom(y_ij[players_ii,players_jj], n_ij[players_ii,players_jj], P_NbyN_prime[players_ii,players_jj], log=T)) + 
        sum(dbinom(y_ij[players_jj,players_ii], n_ij[players_jj,players_ii], P_NbyN_prime[players_jj,players_ii], log=T))
      C_minus <- single_p_ij_normal_overlap(entry_i = p_i,entry_j =p_j, K = K,P_matrix = p_prime,
                                            S = S_current,truncations = truncations_current,diag0.5 = diag0.5)
      
      
      #proposing a new p_ij
      p_scanning = p_prime
      p_scanning[p_i,p_j] <- rtruncnorm(1, mean = p_prime[p_i,p_j],sd = sigma_p[p_i,p_j], a =   0.5, b  = beta_max)
      p_scanning[p_j,p_i] <- 1 - p_scanning[p_i,p_j]
      
      #updating P_NbyN_prime for all players in cluster ii and cluster jj
      P_NbyN_scanning = P_NbyN_prime
      P_NbyN_scanning[players_ii,players_jj] <- p_scanning[p_i,p_j]
      P_NbyN_scanning[players_jj,players_ii] <- p_scanning[p_j,p_i]
      
      #the new likelihood and prior values for those players
      A_plus = sum(dbinom(y_ij[players_ii,players_jj], n_ij[players_ii,players_jj], P_NbyN_scanning[players_ii,players_jj], log=T)) + 
        sum(dbinom(y_ij[players_jj,players_ii], n_ij[players_jj,players_ii], P_NbyN_scanning[players_jj,players_ii], log=T))
      C_plus <- single_p_ij_normal_overlap(entry_i = p_i,entry_j =p_j, K = K, P_matrix = p_scanning,
                                           S = S_current,truncations = truncations_current,diag0.5 = diag0.5) 
      
      
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







#------------------------#
# POMM S update   #
#------------------------#




S_POMM_update_fixed_P_alpha_z = function(z_current, p_current, 
                                         A_current,C_current,y_ij,n_ij,labels_available,
                                         upper.tri.non.zero,K,truncations_current, alpha_current,beta_max, S_current,acc.count_S, sigma_S){
  

  S_prime <- rtruncnorm(1,a = 0.0000000001,b = 0.9,mean = S_current,sd = sigma_S)
  
  
  
  C_prime <- l_like_p_ij_normal_overlap(K = K, P_matrix = p_current,S = 
                                          S_prime, 
                                        truncations = truncations_current,
                                        diag0.5 = T) + dunif(alpha_current,0,4,log = T) + dunif(S_prime,min = 0,max=1, log = T)
  
  
  
  log_r=  C_prime -  C_current
  
  #create statements that check conditiond to accept move
  MH_condition_S= min(log_r,0)>=log(runif(1))
  if(MH_condition_S){
    acc.count_S=acc.count_S+1
    S_current <- S_prime
    C_current<- C_prime
    p_current<- p_current
    #A_current <- A_prime
  }
  
  
  return(list(acc.moves = acc.count_S,
              sigma_S=sigma_S,
              p_current=p_current,
              C_current = C_current, 
              A_current = A_current,
              S_current = S_current))
  
} 





#---------
#alpha update
#---------



alpha_POMM_update_fixed_P_S_z = function(z_current, p_current, 
                                         A_current,C_current,y_ij,n_ij,labels_available,
                                         upper.tri.non.zero,K,alpha_current,beta_max, S_current, acc.count_alpha, sigma_alpha, truncations_current){
  
  
  
  #proposing a new alpha
  
  alpha_prime <- rtruncnorm(1,a = 0.1,b = 0.9,mean = alpha_current,sd = sigma_alpha)
  truncations_prime <- improper_prior5(K,beta_max,alpha = alpha_prime,diag0.5 = T)
  
  
  C_prime <- l_like_p_ij_normal_overlap(K = K, P_matrix = p_current,S = 
                                          S_current, 
                                        truncations = truncations_prime,
                                        diag0.5 = T) + dunif(alpha_prime,0,4,log = T) + dunif(S_current,min = 0,max=1, log = T)
  
  
  
  
  #A_prime <- sum(dbinom(y_ij, n_ij, p_ij_prime, log = T))
  
  log_r= C_prime - C_current
  
  #create statements that check conditiond to accept move
  MH_condition_alpha= min(log_r,0)>=log(runif(1))
  if(MH_condition_alpha){
    acc.count_alpha= acc.count_alpha+1
    alpha_current <- alpha_prime
    truncations_current<- truncations_prime
    C_current<- C_prime
    
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






z_update_adaptive = function(z_current, A_current,B_current,y_ij,n_ij,P_matrix,labels_available,upper.tri.non.zero,gamma_vec,N,K, acc.count_z, sigma_z){
  A_prime<- A_current
  B_prime<- B_current
  z_prime=z_current
  P_NbyN_prime<- calculate_victory_probabilities(vec2mat_0_P(z_prime,P_matrix),P_matrix)
  n_prime = matrix(0,nrow(P_matrix),1)
  for(h in 1:nrow(P_matrix)){
    n_prime[h] = sum(which(z_prime==h))
  }
  
  scanning_order = sample(1:N,N, replace=F)
  # full sweep
  for(ii in scanning_order){
    
    z_scanning=z_prime
    #save current label of z_ii
    k_prime<-z_prime[ii]
    
    # PROPOSE A Z_PRIME for z_prime
    
    # Calculate distances between labels
    label_distances <- labels_available - k_prime
    
    # Calculate probabilities based on label distances and sigma.z
    probabilities <- label_probability(label_distances, sigma_z[ii])
    
    # Adjust probabilities based on sigma.z
    probabilities <- probabilities[-k_prime] / sum(probabilities[-k_prime])
    
    # Sample a new label using the adjusted probabilities
    k_scanning <- sample(x = setdiff(labels_available, k_prime), size = 1, replace = F)
    
    z_scanning[ii] <- k_scanning
    
    # while (TRUE) {
    #   if (length(unique(z_scanning)) == K) {
    #     break
    #   } else {
    #     # resample new_label and try again
    #     k_scanning <- sample(labels_available,size = 1)
    #     z_scanning[ii] <- k_scanning
    #     print("!")
    #   }
    # }
    
    #compute the likelihood of the data with the current assignment just for ii
    A_minus = sum(dbinom(y_ij[ii,], n_ij[ii,], P_NbyN_prime[ii,], log=T)) + sum(dbinom(y_ij[,ii], n_ij[,ii], P_NbyN_prime[,ii], log=T)) 
    
    #update P_NbyN
    P_NbyN_scanning = P_NbyN_prime
    for(nodes in 1:N){
      P_NbyN_scanning[ii,nodes]<- P_matrix[k_scanning,z_scanning[nodes]]
      P_NbyN_scanning[nodes,ii]<- P_matrix[z_scanning[nodes],k_scanning]
    }
    #compute the likelihood of the same points with the new assignment
    A_plus = sum(dbinom(y_ij[ii,], n_ij[ii,], P_NbyN_scanning[ii,], log=T)) + sum(dbinom(y_ij[,ii], n_ij[,ii], P_NbyN_scanning[,ii], log=T)) 
    
    #Updating the likelihood
    A_scanning = A_prime - A_minus + A_plus
    
    
    n_scanning<- n_prime
    n_scanning[c(k_prime, k_scanning)] <- n_prime[c(k_prime, k_scanning)] + c(-1, 1)
    
    B_scanning<- ddirichlet_multinomial(N,K,n_k = n_scanning,my_alpha = gamma_vec)
    
    log_r= A_scanning - A_prime + B_scanning - B_prime
    #create statements that check conditiond to accept move
    GS_condition= min(log_r,0)>=log(runif(1))
    if(GS_condition){
      acc.count_z[ii]=acc.count_z[ii]+1
      z_prime<-z_scanning
      A_prime<- A_scanning
      B_prime<- B_scanning
      P_NbyN_prime <- P_NbyN_scanning
      n_prime <- n_scanning
    }
    #labels_available are the same
    #else, z_prime[ii] stays equal to z_current[ii]
  }
  z_current<-z_prime
  A_current<- A_prime
  B_current<- B_prime
  
  
  
  return(list(acc.moves = acc.count_z,A_current = A_current,B_current = B_current, z_current= z_current))
} 
























#------------------
# Adaptive Proposal
#------------------


tuning_proposal<- function(iteration, acceptance_count, sigma, acceptanceTarget, min_sigma){
  #compute acceptance rate
  acceptanceRate <- acceptance_count / iteration
  #setting the change in the variance
  delta= min(0.01, iteration**(-1/2))
  #passing top the log scale, to have a finer scale
  lsi= log(sigma)
  #if we are accepting too much ==> increase variance, otherwise, reduce it
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
