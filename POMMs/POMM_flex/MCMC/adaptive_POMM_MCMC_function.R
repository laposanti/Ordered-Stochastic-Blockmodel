
adaptive_MCMC_POMM <- function(Yij_matrix, Nij_matrix,init , z_true, overlap, alpha,p_true, N,K, N_iter, targ_rate, beta_max, gamma_vec, diag0.5){
  
  upper.tri.non.zero = which(Nij_matrix > 0,arr.ind = T)
  #working with just with non-zero-values
  n_ij = Nij_matrix[upper.tri.non.zero]
  y_ij = Yij_matrix[upper.tri.non.zero]
  
  #setting containers
  z_container= matrix(0, nrow=N, ncol=N_iter)
  alpha_container = matrix(0, nrow=1, ncol=N_iter)
  overlap_container =  matrix(0, nrow=1, ncol=N_iter)
  p_container = array(0, dim=c(K,K,N_iter))
  A_container = matrix(0, nrow=1, ncol=N_iter)
  B_container = matrix(0, nrow=1, ncol=N_iter)
  C_container = matrix(0, nrow=1, ncol=N_iter)
  
  
  #initializing quantities
  z_current=  init$z0
  alpha_current<- init$alpha0
  truncations_current <- improper_prior5(K,beta_max,alpha = alpha_current,diag0.5 = diag0.5)
  overlap_current= init$overlap0
  p_current = simulating_POMM_powerlaw2(K,alpha_current,truncations_current,beta_max,diag0.5 = diag0.5)
  
  n_k_current = as.vector(table(z_current))
  z_mat_current = vec2mat(z_current)
  #p_ij_function = calculate_victory_probabilities(z_mat_current,P_true)
  aux = p_current%*%t(z_mat_current)
  p_nbyn_current = z_mat_current%*%aux
  p_ij_current = p_nbyn_current[upper.tri.non.zero]
  
  labels_available = 1:K
  
  A_current= sum(dbinom(y_ij, n_ij, p_ij_current, log = T))
  B_current=ddirichlet_multinomial(N,K,n_k = n_k_current ,my_alpha = gamma_vec)
  C_current =  l_like_p_ij_normal_overlap(K = K, P_matrix = p_current,overlap =overlap_current, truncations = truncations_current,diag0.5 = T) + dlnorm_param(alpha_current)
  
  z_container[,1] = z_current
  alpha_container[1] = alpha_current
  overlap_container[1] = overlap_current
  p_container[,,1] = p_current
  A_container[1]=A_current
  B_container[1]=B_current
  C_container[1]=C_current
  
  #containers for the counts of accepted proposals
  acc.count_z = rep(1,N)
  acc.count_overlap=1
  acc.count_p = matrix(1,K,K)
  acc.count_alpha = 1
  
  sigma_z <- rep(0.5,N)
  sigma_alpha = 0.2
  sigma_overlap =0.2 
  sigma_p= matrix(rep(0.2, K**2),K,K)
  
  sigma_z_container<- matrix(0, N, N_iter)
  sigma_alpha_container <- matrix(0, N_iter,1)
  sigma_overlap_container<- matrix(0, N_iter,1)
  sigma_p_container<- array(0, dim=c(K,K, N_iter))
  
  
  #setting time tracker
  pb=txtProgressBar(min=1,max=N_iter)
  j=2
  optimal_p = 0.3
  
  
  
  
  
  #READY TO BOMB!
  
  
  
  for(j in 1:N_iter){
    
    setTxtProgressBar(pb, j)
    # Create a vector of update indices in random order
    updateOrder <- sample(1:4, size = 4,replace = F)
    
    for (updateIndex in updateOrder) {
      if (updateIndex == 1) {
        #alpha_update
        alpha_update = P_POMM_alpha_update(z_current = z_current,
                                           p_current = p_current,
                                           K = K,n_ij = n_ij,
                                           y_ij = y_ij,
                                           A_current = A_current, C_current = C_current,
                                           upper.tri.non.zero = upper.tri.non.zero,
                                           alpha_current = alpha_current,
                                           beta_max = beta_max,overlap_current = overlap_current,acc.count_alpha = acc.count_alpha, 
                                           sigma_alpha = sigma_alpha, truncations_current = truncations_current,labels_available = labels_available)
        
        #updating quantities
        sigma_alpha = alpha_update$sigma_alpha 
        acc.count_alpha = alpha_update$acc.moves
        C_current = alpha_update$C_current
        alpha_current = alpha_update$alpha_current
        truncations_current = alpha_update$truncations_current
        #adaptive standard deviation
        if(j %% 50 == 0){
          sigma_alpha = tuning_proposal(iteration=j,acceptance_count = acc.count_alpha,sigma = sigma_alpha,acceptanceTarget = optimal_p, min_sigma = 0.01)
        }
      }else if (updateIndex == 2) {
        #z UPDATE----------------------------------------------------------------
        
        z_update = z_update_adaptive( z_current = z_current,
                                      P_matrix = p_current,
                                      K = K,n_ij = n_ij,
                                      y_ij = y_ij,
                                      A_current = A_current,B_current=B_current,
                                      upper.tri.non.zero = upper.tri.non.zero,labels_available = labels_available,
                                      gamma_vec = gamma_vec,acc.count_z = acc.count_z,sigma_z = sigma_z)
        
        
        acc.count_z = z_update$acc.moves
        if(j %% 50 == 0){
          for(t in 1:N){
            sigma_z[t] = tuning_proposal(iteration=j,acceptance_count = acc.count_z[t],sigma = sigma_z[t],acceptanceTarget = optimal_p, min_sigma = 0.2)
          }
          
        }
        #updating quantities
        z_current <- z_update$z_current
        B_current<- z_update$B_current
        A_current = z_update$A_current
        
        
        #overlap UPDATE----------------------------------------------------------------
        
      }else if (updateIndex == 3) {overlap_update= P_POMM_overlap_update(z_current = z_current,
                                                                         p_current = p_current,
                                                                         K = K,n_ij = n_ij,
                                                                         y_ij = y_ij,
                                                                         A_current = A_current, C_current = C_current,
                                                                         upper.tri.non.zero = upper.tri.non.zero,
                                                                         alpha_current = alpha_current,
                                                                         beta_max = beta_max,overlap_current = overlap_current,
                                                                         acc.count_overlap = acc.count_overlap,sigma_overlap = sigma_overlap,truncations_current = truncations_current)
      
      #updating quantities
      sigma_overlap = overlap_update$sigma_overlap
      acc.count_overlap = overlap_update$acc.moves
      C_current = overlap_update$C_current
      overlap_current = overlap_update$overlap_current
      if(j %% 50 == 0){
        sigma_overlap = tuning_proposal(iteration=j,acceptance_count = acc.count_overlap,sigma = sigma_overlap,acceptanceTarget = optimal_p, min_sigma = 0.01)
      }
      }else if (updateIndex == 4) {
        #P UPDATE----------------------------------------------------------------
        
        p_update= P_POMM_update_given_overlap1(z_current = z_current,
                                               p_current = p_current,
                                               K = K,n_ij = n_ij,
                                               y_ij = y_ij,
                                               A_current = A_current, C_current = C_current,
                                               upper.tri.non.zero = upper.tri.non.zero,
                                               alpha_current = alpha_current,
                                               beta_max = beta_max,overlap_current = overlap_current,diag0.5=T,
                                               labels_available = labels_available,
                                               acc.count_p =acc.count_p,sigma_p = sigma_p,truncations_current = truncations_current)
        
        #updating quantities
        sigma_p = p_update$sigma_p
        C_current = p_update$C_current
        A_current = p_update$A_current
        p_current = p_update$p_current
        acc.count_p = p_update$acc.moves
        if(j %% 50 == 0){
          j_start = ifelse(diag0.5, yes = 1, no = 0)
          K_stop = ifelse(diag0.5, yes = K-1, no = K)
          for( ii in 1:K_stop){
            for(jj in (ii+j_start):K){
              
              sigma_p[ii,jj] <- tuning_proposal(iteration=j,acceptance_count = acc.count_p[ii,jj],sigma = sigma_p[ii,jj],acceptanceTarget = optimal_p,min_sigma = 0.005)
            }
          }
        }
      }
    }
    if(j %% 1000 == 0){
      print( paste("iteration",j,"acceptance", acc.count_p))
    }
    
    #storing scales
    sigma_alpha_container[j] <- sigma_alpha
    sigma_overlap_container[j]<- sigma_overlap
    sigma_p_container[,,j]<- sigma_p
    sigma_z_container[,j] <- sigma_z
    
    #storing results for inference
    A_container[j] = A_current
    B_container[j] = B_current
    C_container[j]= C_current
    
    alpha_container[j]<- alpha_current
    z_container[,j] <- z_current
    overlap_container[j] = overlap_current
    p_container[,,j] = p_current
    
  }
  
  acceptance_rates <- list(acc.count_p = acc.count_p, acc.count_z = acc.count_z, acc.count_alpha = acc.count_alpha,acc.count_overlap=acc.count_overlap)
  st.deviations<- list(sd_p = sigma_p_container, sd_alpha = sigma_alpha_container, sd_overlap - sigma_overlap_container, sd_z = sigma_z_container )
  
  return(list(Yij_matrix=Yij_matrix, Nij_matrix=Nij_matrix,init = init, z_true=z_true, overlap=overlap, alpha=alpha,p_true=p_true,alpha_container= alpha_container,z_container= z_container,overlap_container = overlap_container,p_container= p_container, A_container= A_container, B_container= B_container, C_container= C_container, acceptance_rates= acceptance_rates, st.deviations=st.deviations))}
