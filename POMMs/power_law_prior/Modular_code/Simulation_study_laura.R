
  K_true =  K
  labels_available = 1:K
  #expected number of players per block
  while(TRUE){
    expected_players = round(gamma_vec/sum(gamma_vec)*N,0)
    if(sum(expected_players)==N){
      break
    }else if(sum(expected_players) < N){
      expected_players[3] = expected_players[3] +1
      break
    }else{
      expected_players[3] = expected_players[3] -1
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
  
  
  model='POMM'
  if(model == 'POMM'){
    # simulating the matrix p according to the SST property
    p = simulating_POMM_powerlaw(K,alpha,beta_max)$matrix
    diag(p) <- rep(0.5,K)
  }else if(model == 'Simple'){
    p = matrix(rbeta(K*K,1,1),K,K)
    diag(p) <- rep(0.5,K)
    p=semi_symmetric(p)
  }
  
  #creating the first dataframe
  z_players <- data.frame(id = 1:N, z = z)
  z_mat_true=vec2mat(z)
  matrix_z_p_true = p%*%t(z_mat_true)
  p_n_true = z_mat_true%*%matrix_z_p_true
  
  
  
  #-----------------------
  #Up to here, fixed part, and the data
  #------------------------
  #We have the true z_vector
  z_true = z
  #And, the true P matrix (KxK)
  P_true= p
  #We finally have the P_(z_i,z_j) (NxN)
  p_NbyN_true = p_n_true
  
  #Simulation study #
  test_size=1
  
  y_ij_container = array(0, dim=c(N,N,test_size))
  n_ij_container = array(0, dim=c(N,N,test_size))
  for(ii in 1:test_size){
  #Sampling the blocks, stratified sampling
  
  
  block_size = colSums(z_mat_true)

  
  n_ij = matrix(0, N,N)
  i=1
  while(i<=M){
    #1) sample the blocks
    block1 = sample(x = labels_available,size = 1,prob = 1/gamma_vec/1/sum(gamma_vec))
    block2 = sample(x = labels_available,size = 1,prob = 1/gamma_vec/1/sum(gamma_vec))
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
    
    
    i =i+1
  }
  
  uppertri_nonzero_entries = which(upper.tri(n_ij) & n_ij > 0,arr.ind = T)
  
  y_ij = matrix(0, N,N)
  
  # Select upper triangular, non-zero entries in n_ij
  n_ij_upper <- n_ij[uppertri_nonzero_entries]
  p_ij <- p_n_true[uppertri_nonzero_entries]

  
  y_ij[uppertri_nonzero_entries] <- rbinom(length(uppertri_nonzero_entries)/2, n_ij_upper, p_ij)
  y_ij[lower.tri(y_ij)] <- n_ij[lower.tri(n_ij)] - t(y_ij)[lower.tri(y_ij)]
  
  n_ij_container[,,ii]=n_ij
  y_ij_container[,,ii]=y_ij
  if(ii %%(test_size/10) == 0){
    print(paste((ii/test_size*100),'-th percent process complete'))
  }
  }
  
  #Generating random permutations of z_true

  random_permutations = matrix(0, nrow=N, ncol=test_size)
  for(i in 1:(test_size-1)){
    random_permutations[,i] = sample(x=z_true, size = N, replace=F)
  }
  #adding to the random permutations also the z_true partition
  random_permutations[,test_size] =z_true
  
  #for each permutation, we compute the P_(z_i,z_j) (NxN)
  P_zi_zj_test = array(0,dim=c(N,N, test_size))
  for( i in 1:test_size){
    z_mat_aux = vec2mat(random_permutations[,i])
    P_zi_zj_test[,,i]= calculate_victory_probabilities(z_mat = z_mat_aux,P =P_true )
  }
  
  MSE_P_zi_zj_test = matrix(0,nrow=test_size,ncol=1)
  for( i in 1:test_size){
    MSE_P_zi_zj_test[i]= mean((P_zi_zj_test[,,i]-p_NbyN_true)**2)
  }
  
  #computing the likelihood 
  A = matrix(0, nrow=test_size, ncol=1)
  for(j in 1:test_size){
    A_aux = matrix(0, nrow=test_size, ncol=1)
    for(i in 1:test_size){
      A_aux[i] = sum(dbinom(y_ij_container[,,i],n_ij_container[,,i],prob = P_zi_zj_test[,,j],log = T))
    }
    A[j] = sum(A_aux)
  
    if(j %%(test_size/10) == 0){
      print(paste((j/test_size*100),'-th percent process complete'))
    }
  }
  
  ggplot(simulation_study_df, aes(x = MSE, y = likelihood)) +
    geom_point(color = "black") +  # Set the default color to black
    geom_point(data = tail(simulation_study_df, 1), color = "red")  # Add a red point for the true z value
  
  
  
  # Sample corresponding entries in y_ij as binomial
  
  # y_ij = array(0, c(N,N, 4000))
  # for(k in 1:4000){
  # for(i in 1:N){
  #   for(j in 1:N){
  #     if(n_ij[i,j]>0 & (j>i)){
  #     y_ij[i,j,k] = rbinom(1, n_ij[i,j], p_n_true[i,j])
  #     y_ij[j,i,k] = n_ij[i,j] -y_ij[i,j,k]
  #   }}}
  # }
  
  # y_ij <- array(0, dim = c(N, N, 4000))
  # 
  # for (k in 1:4000) {
  #   y_ij_k <- matrix(0, nrow = N, ncol = N)
  #   
  #   y_ij_k[uppertri_nonzero_entries] <- rbinom(length(uppertri_nonzero_entries), n_ij_upper, p_ij)
  #   y_ij_k[lower.tri(y_ij_k)] <- n_ij[lower.tri(n_ij)] - t(y_ij_k)[lower.tri(y_ij_k)]
  #   
  #   y_ij[,,k] <- y_ij_k
  # }
  
  
  # y_ij_summary = matrix(0,N,N)
  # for(i in 1:N){
  #   for(j in 1:N){
  #     y_ij_summary[i,j] = mean(y_ij[i,j,])
  #   }
  # }
  # 
  # y_prox = y_ij_summary/n_ij
  # y_prox= na.replace(y_prox,0)
  # #y_ij[lower.tri(y_ij)] = n_ij[lower.tri(n_ij)] - t(y_ij)[lower.tri(y_ij)]
  # 
  # similarity_plot(y_prox,z,z)
  # 
  # check1 = data.frame(X = 1:length(y_ij[uppertri_nonzero_entries]),
  #                     Y = y_ij_summary[uppertri_nonzero_entries]/n_ij[uppertri_nonzero_entries],
  #                     P = p_ij)
  # ggplot(check1) +
  #   geom_line(aes(x = X, y =(Y)), stat = "identity", col = "red") +
  #   geom_line(aes(x = X, y = (P)), stat = "identity", col = "blue") +
  #   xlab("X") +
  #   ylab("Value") +
  #   ggtitle("Bar Plot of Y and P") +
  #   theme_minimal()
  
  # #---
  # a1 = matrix(y_ij, ncol =1)
  # b1 = matrix(p_n_true, ncol=1)
  # cor(a1,b1)

# Create sample data
n <- 100
mat <- matrix(runif(n*n,-1,1), nrow=n, ncol=n)

# Create index of non-negative cells
maty.n <- upper.tri(mat, diag = TRUE)
maty.n <- (upper.tri_n_ij & (mat > 0))

df = data.frame(var1 = 0, var2=0, y.n = 0)
for(i in 1:N){
  for(j in 1:N){
    df=  rbind(df, data.frame(var1 = i, var2=j, y.n= maty.n[i,j]))
  }
}



# Convert matrix to long format data.frame
df <- reshape2::melt(maty.n)
names(df) <- c("Var1", "Var2", "y.n")

# Create heatmap
ggplot(df, aes(x = var2, y = desc(var1), fill = factor(y.n))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("white", "red")) +
  labs(x = "", y = "") +
  theme_bw()

# Create heatmap
ggplot(df, aes(x = factor(var2,levels=my_levels.x), y = factor(var1,my_levels.y), fill = factor(y.n))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("white", "red")) +
  labs(x = "", y = "") +
  theme_bw()


ggplot(df, aes(x = var2, y = desc(var1), fill = factor(y.n))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("white", "red")) +
  labs(x = "", y = "") +
  theme_bw()


mean_mat = matrix(0,N,N)
for(i in 1:N){
  for(j in 1:N){
    if(n_ij[i,j]>0){
      mean_mat[i,j] = sum(y_ij_arr[i,j,])/(n_ij[i,j]*iters)
    }else{mean_mat[i,j] =0}
    
  }
}
