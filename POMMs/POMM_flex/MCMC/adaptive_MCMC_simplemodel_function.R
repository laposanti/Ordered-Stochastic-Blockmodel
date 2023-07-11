

adaptive_MCMC_simple_model <- function(Yij_matrix, Nij_matrix,init, z_true,p_true, N,K, N_iter, gamma_vec, diag0.5, seed){
  
  
  upper.tri.non.zero = which(Nij_matrix > 0,arr.ind = T)
  #working with just with non-zero-values
  n_ij = Nij_matrix[upper.tri.non.zero]
  y_ij = Yij_matrix[upper.tri.non.zero]
  

  
  #------
  # setting containers
  
  z_container= matrix(0, nrow=N, ncol=N_iter)
  p_container = array(0, dim=c(K,K,N_iter))
  A_container = matrix(0, nrow=1, ncol=N_iter)
  B_container = matrix(0, nrow=1, ncol=N_iter)
  C_container = matrix(0, nrow=1, ncol=N_iter)
  
  #--------
  #initializing quantities
  
  p_current= matrix(runif(K**2,1-beta_max,beta_max),K,K)
  diag(p_current) <- rep(0.5,K)
  p_current = semi_symmetric(p_current)
  
  z_current= z_true
  n_k_current = as.vector(table(z_current))
  z_mat_current = vec2mat(z_current)
  
  aux = p_current%*%t(z_mat_current)
  p_nbyn_current = z_mat_current%*%aux
  p_ij_current = p_nbyn_current[upper.tri.non.zero]
  
  A_current= sum(dbinom(y_ij, n_ij, p_ij_current, log = T))
  B_current=ddirichlet_multinomial(N,K,n_k = n_k_current ,my_alpha = gamma_vec)
  C_current =  sum(dunif(p_current[upper.tri(p_current)],min = 1-beta_max, max = beta_max, log = T))
  
  
  
  labels_available = 1:K
  
  #---------
  acc.count_z = rep(1,N)
  acc.count_p = matrix(1,K,K)
  #updating containers
  z_container[,1] = z_current
  p_container[,,1] = p_current
  A_container[1]=A_current
  B_container[1]=B_current
  C_container[1]=C_current
  
  #containers for the counts of accepted proposals
  sigma_z <- rep(0.5,N)
  sigma_p= matrix(rep(0.2, K**2),K,K)
  
  sigma_z_container<- matrix(0, N, N_iter)
  sigma_p_container<- array(0, dim=c(K,K, N_iter))
  
  #-----------
  
  #setting time tracker
  pb=txtProgressBar(min=1,max=N_iter)
  j=2
  
  optimal_p =0.25
  #READY TO BOMB!
  while (j < N_iter + 1) {
    setTxtProgressBar(pb, j)
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
        sigma_z[t] = tuning_proposal(iteration=j,acceptance_count = acc.count_z[t],sigma = sigma_z[t],acceptanceTarget = optimal_p,min_sigma = 0.2)
      }
      
    }
    #updating quantities
    z_current <- z_update$z_current
    B_current<- z_update$B_current
    A_current = z_update$A_current
    
    #A_seq[j]= sum(dbinom(y_ij, n_ij, p_ij_current, log = T)) + ddirichlet_multinomial(N,K_true,n_k = n_k_current ,my_alpha = gamma_vec)
    z_container[, j] = z_current
    
    #P UPDATE----------------------------------------------------------------
    
    p_update= P_simple_update_adaptive(z_current = z_current,
                                       p_current  = p_current,K = K,n_ij = n_ij,y_ij = y_ij,diag0.5 = diag0.5,
                                       A_current = A_current,C_current = C_current,
                                       upper.tri.non.zero = upper.tri.non.zero,acc.count_p =  acc.count_p,sigma_p = sigma_p,beta_max=beta_max,labels_available = labels_available)
    
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
    
    #storing scales
    sigma_p_container[,,j]<- sigma_p
    sigma_z_container[,j] <- sigma_z
    
    #updating quantities
    p_current=  p_update$p_current
    A_current = p_update$A_current
    C_current = p_update$C_current
    
    #storing results for inference
    A_container[j] = A_current
    B_container[j] = B_current
    C_container[j]= C_current
    p_container[,,j] = p_current
    
    j=j+1
  }
  acceptance_rates <- list(acc.count_p = acc.count_p, acc.count_z = acc.count_z)
  st.deviations<- list(sd_p = sigma_p_container, sd_z = sigma_z_container )
  
  return(list(Yij_matrix=Yij_matrix, Nij_matrix=Nij_matrix,init=init,z_true=z_true,p_true=p_true,z_container= z_container,p_container= p_container, A_container= A_container, B_container= B_container, C_container= C_container,acceptance_rates=acceptance_rates,st.deviations=st.deviations, seed=seed))}
