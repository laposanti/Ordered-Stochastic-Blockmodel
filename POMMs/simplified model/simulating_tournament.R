


sample_truncated_poisson(10,10,1,n_ij_max)

K=4
N=10
min_clust=1
max_clust=5
alpha=1
beta=1

simulating_tournament<- function(N,alpha,beta, min_clust,max_clust,n_ij_max){
  #simulating K from a truncated Poisson(1)
  K =  sample_truncated_poisson(N=1,lambda=1,a= min_clust,b = max_clust)
  
  #simulating z|K from dirichlet multinomial with gamma=1
  repeat {
    z<- dir_multinom_r(K,N,1)
    l = length(unique(z))
    #sample until number of sampled labels == number of clusters K
    if (l==K){
      break
    }}
  
  # simulating the matrix p according to the SST property
  p = sampling_SST_matrix_beta(K,alpha, beta)
  p[is.na(p)] <- 0
  p = (1 -t(p*upper.tri(p))) * (lower.tri(p, diag = F)*1) + upper.tri(p, diag = T) *p
  
  #creating the first dataframe
  z_players <- data.frame(id = 1:N, z = z)
  
  
  player1=0
  player2=0
  same_set = 0
  same_K = 0 
  different_rows = 0
  #checking some conditions
  while((same_set+same_K+different_rows)!=3){
    player1= sample(x=c(1:N), size=N, replace = F)
    player2= sample(x=c(1:N), size=N, replace= T)
    same_set = setequal(z_players$id, unique(c(player1,player2)))
    same_K = ncol(p) == length(unique(z_players$z))
    different_rows = table(player1 != player2)["TRUE"] == N
  }
  
  z_player1 = left_join(data.frame(id = player1),z_players, by = "id")$z
  z_player2 = left_join(data.frame(id = player2),z_players, by = "id")$z
  
  #here I want the number of games to be inversely related to clusters' strength
  aux = matrix(0, K,K)
  aux = abs(row(aux) - col(aux))  + 1
  diag(aux) = (diag(row(aux))**(-1))/sum(diag(row(aux)**(-1)))
  aux = aux**(-1)/sum(aux**(-1))
  n_ij = sample_truncated_poisson(N,lambda= exp(aux[cbind(z_player1,z_player2)])*n_ij_max,a=1,b=n_ij_max)
  p_ij = p[cbind(z_player1,z_player2)]
  y_ij = rbinom(N, n_ij, p_ij)
  matches = data.frame(player1, player2, n_ij, y_ij)
  
  return(list(z_true = z_players, matches_results = matches, p_true = p, K_true = K))}
