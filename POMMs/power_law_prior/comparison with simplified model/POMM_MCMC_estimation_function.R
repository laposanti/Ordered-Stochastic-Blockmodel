


POMM_MCMC_estimation = function(n_ij_matrix, y_ij_matrix, beta_max,K, N_iter, burn_in){
  K_true=K
  #selecting just those values such that there is inside a non-0 entry
  upper.tri_n_ij = upper.tri(n_ij_matrix)
  non_negative_n_ij = which(upper.tri_n_ij & n_ij_matrix > 0, arr.ind = T)
  
  
  
  #retrieving the non-0 values
  n_ij = n_ij_matrix[non_negative_n_ij]
  y_ij = y_ij_matrix[non_negative_n_ij]
  
  

  #Containers for p
  #---
  alpha.container = matrix(0,N_iter,1)
  p.container = array(0,c(K_true,K_true,N_iter))
  
  
  
  #Containers for z
  #---
  z.container = matrix(0,N,N_iter)
  A_seq = matrix(0,N_iter,1)
  
  
  #Containers for p
  #---
  alpha.container = matrix(0,N_iter,1)
  p.container = array(0,c(K_true,K_true,N_iter))
  
  
  #Initialization for the z vector
  #---
  init = kmeans(x = y_ij_matrix,centers =K)$cluster
  z_current= init
  z.container[,1] = z_current
  
  #Initialization for the p matrix
  #---
  alpha_current = 1
  p_current = simulating_POMM_powerlaw(K_true,alpha_current,beta_max = beta_max)
  
  alpha.container[1] = alpha_current
  p.container[,,1] = p_current$matrix
  
  #Setting up quantities needed within computations
  
  # here we are transforming a Nx1 vector, containg labels 1...K into
  # NXK matrix. z_mat_current[i,k] =1 if node i is in cluster k
  z_mat_current= vec2mat(z_current)
  n_k_current = colSums(z_mat_current)
 
  #creating an NxN matrix where p_n_current[i,j] is the probability that player i wins vs player j
  matrix_z_p_current = p_current$matrix%*%t(z_mat_current)
  p_n_current = z_mat_current%*%matrix_z_p_current
  p_ij_current =  p_n_current[non_negative_n_ij]
  
  #containers for the counts of accepted proposals
  acc.count_z = 0
  acc.count_p = 0
  
  #setting time tracker
  pb=txtProgressBar(min=1,max=N_iter)
  j=2
  
  
  #READY TO BOMB!
  while (j < N_iter + 1) {
    setTxtProgressBar(pb, j)
    
    #Complete sweeep of z vector
    #----
    sweeping_order = sample(x=c(1:N),size=N, replace=F)
    for(i in sweeping_order){
      #proposing a new label
      new_label = adjacent_label_sampler(labels_available = labels_available, z_current[i])
      
      #updating labels
      z_prime = z_current
      z_prime[i] = new_label
      
      #computing the new victory probabilities
      while(TRUE){
        z_mat_prime= vec2mat(z_prime)
        if(ncol(z_mat_prime) == K_true){
          break
        } else {
          #if there is an error, resample new_label and try again
          new_label = adjacent_label_sampler(labels_available=labels_available, z_current[i])
          z_prime[i] = new_label
        }
      }
      matrix_z_prime_p_current = p_current$matrix%*%t(z_mat_prime)
      p_n_z_prime = z_mat_prime%*%matrix_z_prime_p_current
      p_ij_z_prime =  p_n_z_prime[non_negative_n_ij]
      
      n_k_prime = colSums(z_mat_prime)
      #acceptance ratio
      r = (sum(dbinom(y_ij, n_ij, p_ij_z_prime, log=T)) + ddirichlet_multinomial(N,K_true,n_k = n_k_prime,my_alpha = gamma_vec)) - 
        (sum(dbinom(y_ij, n_ij, p_ij_current, log = T)) + ddirichlet_multinomial(N,K_true,n_k = n_k_current ,my_alpha = gamma_vec))
      
      alpha_r = min(1, exp(r))
      u = runif(1)
      
      #if accepted
      if(u<alpha_r){
        acc.count_z = acc.count_z + 1
        z_current=z_prime
        n_k_current = n_k_prime
        p_ij_current= p_ij_z_prime
        #if not accepted
      }
      z.container[, j] = z_current
    }
    
    #Update of the P matrix
    #----
    sigma_prime = 1
    
    #proposing a new alpha
    alpha_prime <- sample_norm_trunc(1,m = alpha_current,s =sigma_prime,a = 0.005,b=3)
    
    #generating a proposal matrix
    p_prime = simulating_POMM_powerlaw(K_true,alpha_prime,beta_max)
    
    
    
    #Updating p_ij_prime with the last membership
    z_mat_current = vec2mat(z_current)
    matrix_z_p_prime  = p_prime$matrix%*%t(z_mat_current)
    p_n_prime = z_mat_current%*%matrix_z_p_prime
    p_ij_prime =  p_n_prime[non_negative_n_ij]
    
    r= (sum(dbinom(y_ij, n_ij, p_ij_prime, log=T)+  l_like_p_ij(p_prime$matrix,p_prime$truncations)))  - 
      (sum(dbinom(y_ij, n_ij, p_ij_current, log = T)+  l_like_p_ij(p_current$matrix,p_current$truncations)))
    
    alpha_r = min(1, exp(r))
    u = runif(1)
    if(u<alpha_r){
      #counting number of accepted proposals
      acc.count_p = acc.count_p+1
      #updating quantities
      p_ij_current = p_ij_prime
      alpha_current = alpha_prime
      p_current = p_prime
    }else{
    }
    #storing results for diagnostics
    p.container[,,j] = p_current$matrix
    alpha.container[j]= alpha_current
    
    A_seq[j] = sum(dbinom(y_ij, n_ij, p_ij_current, log=T)) +  l_like_p_ij(p_current$matrix,p_current$truncations)
    j=j+1
  }
  
  z_burned = z.container[,-c(1:burn_in)]
  p_burned = p.container[,,-c(1:burn_in)]
  alpha_burned = alpha.container[-c(1:burn_in)]
  return(list(z=z_burned,p=p_burned,alpha=alpha_burned))
}

results = POMM_MCMC_estimation(synth$n_ij_true,synth$y_ij_true,beta_max = beta_max,K = K,N_iter = 10000,burn_in = 3000)


