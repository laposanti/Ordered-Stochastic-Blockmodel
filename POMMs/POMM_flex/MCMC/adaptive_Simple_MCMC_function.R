
#making inference on 
#z: the mixture membership vector
#alpha: the imbalance parameter in the preferences
#S: the transitivity of preferences
#P: the victory probabilities

#this code is modular
#estimation_control is a 4x1 vector, denoting which parameters must be estimated. 
#estimation_control[p]=1, means that the parameter p will be estimated. 
#estimation_control[p]=0, means that the parameter p will be fixed at the true value 

#ground_truth is a 4x1 vector, with the true value for the parameters
#if ground_truth is missing, then we are not in a simulation study

#hyper_params should contain a list of hyperparameters
#K: the number of blocks
#beta_max: the highest attainable probability in P
#gamma_vec: the hyperprior on the block sizes
#diag0.5: whether the main diagonal is set to 0.5 (diag0.5=T) or not (diag0.5=F)

adaptive_MCMC_simple <- function(Yij_matrix, Nij_matrix,init , estimation_control, ground_truth,N, N_iter, targ_rate, hyper_params, seed){
  
  #validation
  if(missing(ground_truth)){
    print("Estimation of real data: no ground truth available")
  }else if(sum(names(ground_truth) != c("z" ,     "P" ))!=0){
    print("please provide initial values with the following order and names: c(z, P)")
    stopifnot(FALSE)
  }
  if(sum(names(init) != c("z" ,    "P" ))!=0){
    print("please provide initial values with the following order and names: c(z, P)")
    stopifnot(FALSE)
  }
  
  if(sum(names(estimation_control) != c("z" ,    "P" ))!=0){
    print("please provide initial values with the following order and names: c(z, P)")
    stopifnot(FALSE)
  }
  
  print(paste0("Estimation of ", names(estimation_control)[which(estimation_control==1)]))
  upper.tri.non.zero = which(Nij_matrix > 0,arr.ind = T)
  #working with just with non-zero-values
  n_ij = Nij_matrix
  y_ij = Yij_matrix
  #setting hyperparams
  K = as.numeric(hyper_params$K)
  
  gamma_vec = as.vector(hyper_params$gamma_vec)
  diag0.5 = as.numeric(hyper_params$diag0.5)
  
  #setting containers
  z_container= matrix(0, nrow=N, ncol=N_iter)
  
  
  p_container = array(0, dim=c(K,K,N_iter))
  
  #if the parameters is fixed, setting it to the true value
  if(estimation_control$z==1){
    z_current= matrix(init$z,N,1)
  }else if(missing(ground_truth)){
    print('z must be estimated or fixed to a value')
    stopifnot(FALSE)
  }else{
    z_current=  matrix(ground_truth$z, N, 1)
  }
  
  if(estimation_control$P==1){
    p_current<- as.matrix(init$P)
  }else if(missing(ground_truth)){
    print('P must be estimated or fixed to a value')
    stopifnot(FALSE)
  }else{
    p_current=  as.matrix(ground_truth$P)
  }

  A_container = matrix(0, nrow=1, ncol=N_iter)
  B_container = matrix(0, nrow=1, ncol=N_iter)
  C_container = matrix(0, nrow=1, ncol=N_iter)
  
  
  #initializing quantities
  
  
  n_k_current = as.vector(table(z_current))
  z_mat_current = vec2mat(z_current)
  #p_ij_function = calculate_victory_probabilities(z_mat_current,P_true)
  aux = p_current%*%t(z_mat_current)
  p_nbyn_current = z_mat_current%*%aux

  
  labels_available = 1:K
  
  
  A_current= sum(dbinom(y_ij, n_ij, p_nbyn_current, log = T))
  B_current=ddirichlet_multinomial(N,K,n_k = n_k_current ,my_alpha = gamma_vec)
  C_current =  sum(dbeta(p_current[upper.tri(p_current)],1,1,log = T))
  
  z_container[,1] = z_current
  p_container[,,1] = p_current
  
  A_container[1]=A_current
  B_container[1]=B_current
  C_container[1]=C_current
  
  #containers for the counts of accepted proposals
  acc.count_z = rep(1,N)
  acc.count_alpha = 1
  acc.count_S=1
  acc.count_p = matrix(1,K,K)
  
  
  sigma_z <- rep(0.5,N)
  sigma_p= matrix(rep(0.2, K**2),K,K)
  
  sigma_z_container<- matrix(0, N, N_iter)
  sigma_p_container<- array(0, dim=c(K,K, N_iter))
  
  
  #setting time tracker
  pb=txtProgressBar(min=1,max=N_iter)
  j=2
  optimal_p = 0.3
  
  
  
  
  
  #READY TO BOMB!
  
  
  
  for(j in 1:N_iter){
    
    setTxtProgressBar(pb, j)
    # Create a vector of update indices in random order
    
    
    
    if (estimation_control$z == 1) {
      #z UPDATE----------------------------------------------------------------
      
      z_update = z_update_adaptive( z_current = z_current,
                                    P_matrix = p_current,
                                    K = K,N=N, n_ij = n_ij,
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
      
    }
    if (estimation_control$P== 1) {
      #P UPDATE----------------------------------------------------------------
      
      p_update= P_simple_update(z_current = z_current,
                                              p_current = p_current,
                                              K = K,n_ij = n_ij,
                                              y_ij = y_ij,
                                              A_current = A_current, C_current = C_current,diag0.5=diag0.5,
                                              labels_available = labels_available,
                                              acc.count_p =acc.count_p,sigma_p = sigma_p)
      
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
    
    
    
    #storing scales
    sigma_p_container[,,j]<- sigma_p
    sigma_z_container[,j] <- sigma_z
    
    #storing results for inference
    A_container[j] = A_current
    B_container[j] = B_current
    C_container[j]= C_current
    
    
    z_container[,j] <- z_current
    p_container[,,j] = p_current
    
    if(j %%1000 == 0){
      print(paste0("Iteration ",j))
    }
  }
  
  acceptance_rates <- list(acc.count_p = acc.count_p, acc.count_z = acc.count_z)
  st.deviations<- list(sd_p = sigma_p_container, sd_z = sigma_z_container )
  est_containers = list(z = z_container,P= p_container)
  control_containers = list(A = A_container, B = B_container,C = C_container)
  if(missing(ground_truth)){
    ground_truth= NA
  }
  return(list(Yij_matrix=Yij_matrix, Nij_matrix=Nij_matrix,init = init, ground_truth=ground_truth,est_containers=est_containers, control_containers=control_containers, acceptance_rates= acceptance_rates, st.deviations=st.deviations, seed=seed))}
