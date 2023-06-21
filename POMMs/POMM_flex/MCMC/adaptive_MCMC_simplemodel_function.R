

MCMC_simple_model <- function(Yij_matrix, Nij_matrix,z_true,p_true, N,K, N_iter, gamma_vec, diag0.5){
  
  
  #upper.tri.non.zero= which(n_ij_matrix >= 0)
  #activate for simple model
  upper.tri.non.zero = which(Nij_matrix > 0,arr.ind = T)
  #working with just the upper triangular entries and with non-zero-values
  n_ij = Nij_matrix[upper.tri.non.zero]
  y_ij = Yij_matrix[upper.tri.non.zero]
  
  
  #setting containers
  z_container= matrix(0, nrow=N, ncol=N_iter)
  p_container = array(0, dim=c(K,K,N_iter))
  A_container = matrix(0, nrow=1, ncol=N_iter)
  B_container = matrix(0, nrow=1, ncol=N_iter)
  C_container = matrix(0, nrow=1, ncol=N_iter)
  
  
  #initializing quantities
  p_current= matrix(rbeta(K**2,1,1),K,K)
  p_current = semi_symmetric(p_current)
  
  z_current= kmeans(x = y_ij_matrix,centers = K)$cluster
  n_k_current = as.vector(table(z_current))
  z_mat_current = vec2mat(z_current)
  
  aux = p_current%*%t(z_mat_current)
  p_nbyn_current = z_mat_current%*%aux
  p_ij_current = p_nbyn_current[upper.tri.non.zero]
  
  labels_available = 1:K
  
  A_current= sum(dbinom(y_ij, n_ij, p_ij_current, log = T))
  B_current=ddirichlet_multinomial(N,K,n_k = n_k_current ,my_alpha = gamma_vec)
  C_current =  get_B(p_current,1)
  
  #initializing containers
  z_container[,1] = z_current
  p_container[,,1] = p_current
  A_container[1]=A_current
  B_container[1]=B_current
  C_container[1]=C_current
  
  #containers for the counts of accepted proposals
  acc.count_z = 0
  acc.count_p = 0
  
  
  #setting time tracker
  pb=txtProgressBar(min=1,max=N_iter)
  j=2
  
  
  #READY TO BOMB!
  
  
  while (j < N_iter + 1) {
    setTxtProgressBar(pb, j)
    
    
    z_sweep = z_update_1(z_current, A_current,B_current,y_ij,n_ij,p_current,labels_available = labels_available,upper.tri.non.zero = upper.tri.non.zero,gamma_vec = gamma_vec,K = K)
    
    acc.count_z = acc.count_z+ z_sweep$acc.moves
    z_current= z_sweep$z_current
    n_k_current =z_sweep$n_k_current
    A_current=z_sweep$A_current
    B_current=z_sweep$B_current
    
    #A_seq[j]= sum(dbinom(y_ij, n_ij, p_ij_current, log = T)) + ddirichlet_multinomial(N,K_true,n_k = n_k_current ,my_alpha = gamma_vec)
    z_container[, j] = z_current
    
    
    p_update= P_simple_update1(z_current = z_current,
                               P_matrix = p_current,
                               K = K,n_ij = n_ij,
                               y_ij = y_ij,
                               A_current = A_current,C_current = C_current,upper.tri.non.zero = upper.tri.non.zero)
    
    acc.count_p = acc.count_p + p_update$acc.count_p
    
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
  return(list(Yij_matrix=Yij_matrix, Nij_matrix=Nij_matrix,z_true=z_true,p_true=p_true,z_container= z_current,p_container= p_current, A_container= A_container, B_container= B_container, C_container= C_container))}
