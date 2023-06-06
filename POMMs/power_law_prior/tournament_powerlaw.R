


simulating_tournament_powerlaw<- function(N,alpha, beta_max, min_clust,max_clust,M, n_ij_max){
  #simulating K from a truncated Poisson(1)
  K =  sample_truncated_poisson(N=1,lambda=1,a= min_clust,b = max_clust)
  
  #simulating z|K from dirichlet multinomial with gamma=1
  repeat {
    z<- dir_multinom_r(K,N,1)
    l = length(unique(z$y))
    #sample until number of sampled labels == number of clusters K
    if (l==K){
      break
    }}
  
  
  # simulating the matrix p according to the SST property
  p = simulating_POMM_powerlaw(K,alpha,beta_max)$matrix
  #creating the first dataframe
  z_players <- data.frame(id = 1:N, z = z$y)
  
  
  #the probability of playing is inversely related to the cluster size
  aux2=  data.frame(z = c(1:K), p = sort(z$p, decreasing = T))
  #probability of two players being extracted
  aux2= left_join( z_players,aux2,by="z", multiple='all')
  aux2$p = aux2$p/sum( aux2$p)
  
  df = data.frame(player_1=NA, player_2=NA, n_ij=NA, y_ij=NA, p_ij=NA)
  i=0
  while(i<M){
    pl_1_i <- sample(z_players$id, size=1, prob = aux2$p)
    pl_2_i <-sample(setdiff(z_players$id,pl_1_i), 
                    size=1, prob = aux2$p[setdiff(z_players$id,pl_1_i)])
    matches_i = data.frame(pl_1_i, pl_2_i)
    matches_i$z_1 <- z_players[which(z_players$id == pl_1_i),]$z
    matches_i$z_2 <- z_players[which(z_players$id == pl_2_i),]$z
    matches_i$n_ij = sample_truncated_poisson(1,p[cbind(matches_i$z_1,matches_i$z_2)]*n_ij_max ,a = 1,b = n_ij_max)
    
    matches_i$y_ij =  rbinom(n = 1, size = matches_i$n_ij,prob =  p[cbind(matches_i$z_1,matches_i$z_2)])  
    
    df = rbind(df, data.frame(player_1 = matches_i$pl_1_i, 
                              player_2= matches_i$pl_2_i,
                              n_ij = matches_i$n_ij,
                              y_ij = matches_i$y_ij,
                              p_ij = p[cbind(matches_i$z_1,matches_i$z_2)]))
    i=i+1
  }
  
  
  df = df[-1,]
  
  
  my_df = df
  rownames(my_df) = 1:nrow(my_df)
  head(my_df)
  data_clean = data.frame(player_1=NA, player_2=NA, n_ij=NA, y_ij=NA)
  n = nrow(my_df)
  while(n > 0){
    entries = c(which(my_df$player_1 == my_df$player_1[1] & my_df$player_2 == my_df$player_2[1]), which(my_df$player_1 == my_df$player_2[1] & my_df$player_2 == my_df$player_1[1]))
    df_support = my_df[entries,]
    n_ij = sum(df_support$n_ij)
    y_ij = 0
    for(i in 1:nrow(df_support)){
      if(df_support[i,]$player_1 == my_df$player_1[1]){
        y_ij = y_ij + df_support[i,]$y_ij
      }else if(df_support[i,]$player_1 == my_df$player_2[1]){
        y_ij = y_ij + df_support[i,]$n_ij - df_support[i,]$y_ij
      }
    }
    
    # if(n_victories == sum(which(df_support[,1] == edgelist[i,1] & df_support[,2] == edgelist[i,2]))){
    data_clean = rbind(data_clean, data.frame(player_1 = my_df$player_1[1], player_2 = my_df$player_2[1],n_ij = n_ij, y_ij = y_ij))
    # print("ok")
    # }else(print("error"))
    my_df = my_df[-entries,]
    n = nrow(my_df)
  }
  data_clean= data_clean[-1,]
  
  
  
  return(list(z_true = z_players, matches_results = data_clean, p_true = p, K_true = K))}
