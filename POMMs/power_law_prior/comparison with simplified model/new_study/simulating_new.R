N=100
a=1
b=1
K=5 

M=8000
max_clust=3
min_clust= 3
gamma_vec = c(2/5, 3/5, 4/5)

simulating_tournament_powerlaw_new<- function(N,alpha, beta_max, K,M, n_ij_max, gamma_vec){
  #simulating K from a truncated Poisson(1)
  
  K_true =  K
  
  #simulating z|K from dirichlet multinomial with gamma=1
  while(TRUE){
    KtimesNmatrix = rdirichlet_multinomial(N,K,gamma_vec)
    if(sum((rowSums(KtimesNmatrix)>0)) == K){
      break
    }}
  
  print(rowSums(KtimesNmatrix)>0)
  
  z = matrix(0,N,1)
  for(i in 1:ncol(KtimesNmatrix)){
    z[i]= which(KtimesNmatrix[,i] >0)}
  
  
  # simulating the matrix p according to the SST property
  p = simulating_POMM_powerlaw(K,alpha,beta_max)$matrix

  
  #creating the first dataframe
  z_players <- data.frame(id = 1:N, z = z)
  z_mat_true=vec2mat(z)
  matrix_z_p_true = p%*%t(z_mat_true)
  p_n_true = z_mat_true%*%matrix_z_p_true

  
  cluster_proportions =(as.vector(table(z)/N))
  
  
  
  labels_available = seq(1,K,1)
  
  #the probability of playing is inversely related to the cluster size
  aux2=  data.frame(z = labels_available, p = 1/cluster_proportions)
  
  
  #probability of two players being extracted
  aux2= left_join( z_players,aux2,by="z", multiple='all')
  aux2$p = aux2$p/sum(aux2$p)
  
  n_ij = matrix(0, N,N)

  i=0
  while(i<M){
      pl_1_i <- sample(z_players$id, size=1, prob = aux2$p)
      pl_2_i <- sample(setdiff(z_players$id,pl_1_i), 
                      size=1, prob = aux2$p[setdiff(z_players$id,pl_1_i)])
      n_ij[cbind(pl_1_i,pl_2_i)] = n_ij[cbind(pl_1_i,pl_2_i)] +1 
      n_ij[cbind(pl_2_i,pl_1_i)] = n_ij[cbind(pl_2_i,pl_1_i)] +1 
      i =i+1
    }
  
  #selecting just non-0 entries in the upper.tri matrix
  upper.tri_n_ij = upper.tri(n_ij,diag = T)
  non_negative_n_ij = which(upper.tri_n_ij & (n_ij>0), arr.ind = T)


  
  # maty.n = upper.tri_n_ij & (n_ij>0)
  # 
  # df <- reshape2::melt(maty.n)
  # names(df) <- c("var1", "var2", "y.n")
  # 
  # z_ord = cbind(z, c(1:100))
  # 
  # my_levels.x = z_ord[order(z_ord[,1],decreasing = T),][,2]
  # my_levels.y = z_ord[order(z_ord[,1],decreasing = F),][,2]
  # 

  #----
  #iters = 10000
  
  

  # y_ij_arr = array(0,c(N,N,iters))
  # for(ii in 1:iters){
  #   for(i in 1:N){
  #     for(j in 1:N){
  #       if(n_ij[i,j]>0){
  #         y_ij_arr[i,j,ii] = rbinom(1, n_ij[i,j], p_n_true[i,j])
  #       }else{y_ij[i,j] = 0}
  #     }
  #   } }
  

  y_ij=matrix(0, N,N)
  for(i in 1:N){
    for(j in 1:N){
      if(n_ij[i,j]>0){
        y_ij[i,j] = rbinom(1, n_ij[i,j], p_n_true[i,j])
      }else{y_ij[i,j] = 0}
    }
  }
  y_ij[lower.tri(y_ij)] = n_ij[lower.tri(y_ij)] - t(y_ij)[lower.tri(y_ij)]

  
  
  similarity_plot(y_ij,z,z)


  #---
  print(paste("corr(y_ij,p_ij):",cor(y_ij[non_negative_n_ij],p_n_true[non_negative_n_ij])))
  print(paste("total number of matches:", sum(n_ij)/2))
  
  print(paste("total number of victories:", sum(y_ij)))
  

  #retrieving the non-0 values
  # 
  # n_ij = n_ij_matrix[non_negative_n_ij]
  #     temp = rbinom(1, 1, p_n_true[pl_1_i,pl_2_i ])
  #     if(temp ==1){
  #       y_ij[pl_1_i,pl_2_i ] = y_ij[pl_1_i,pl_2_i ] + 1
  #     }else{
  #       y_ij[pl_2_i,pl_1_i ] = y_ij[pl_2_i,pl_1_i ] + 1
  #     }
  #     i= i+1
  #   }

  return(list(n_ij_true = n_ij, y_ij_true = y_ij, p_ij_true = p_n_true, P_matrix = p, z_true= z))
}



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

