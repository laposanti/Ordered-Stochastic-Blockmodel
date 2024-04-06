

#first function: everything done slowly
#inputs:
#1) current assignment
#2) data: Y_ij, N_ij
#3) other parameters: theta
#4) proposed assignment

#outputs
#1) the two likelihoods

naive_version = function(Y_ij, N_ij, z_current, theta, z_prime){
  z_mat = vec2mat_0_P(z_current, theta)
  P_ij = calculate_victory_probabilities(z_mat, theta)
  A_cur = dbinom(Y_ij[upper.tri(Y_ij)], N_ij[upper.tri(N_ij)], P_ij[upper.tri(P_ij)],log = T)

  z_mat_prime = vec2mat_0_P(z_prime, theta)
  P_ij_prime = calculate_victory_probabilities(z_mat_prime, theta)
  A_prime = dbinom(Y_ij[upper.tri(Y_ij)], N_ij[upper.tri(N_ij)], P_ij_prime[upper.tri(P_ij)],log = T)
  return(list(A_cur = sum(A_cur), A_prime = sum(A_prime)))
}

current_version = function(Y_ij, N_ij, z_current, theta,A, z_prime,item_changed,n){
  z_mat = vec2mat_0_P(z_current, theta)
  P_ij = calculate_victory_probabilities(z_mat, theta)
  
  A_cur = dbinom(Y_ij[upper.tri(Y_ij)], N_ij[upper.tri(N_ij)], P_ij[upper.tri(P_ij)],log = T)
  
  upper.tri_Y_ij = Y_ij*upper.tri(Y_ij)
  upper.tri_N_ij = N_ij*upper.tri(N_ij)
  upper.tri_P_ij = P_ij*upper.tri(P_ij)
  #compute the likelihood of the data with the current assignment just for i_th_turn
  A_minus = sum(dbinom(upper.tri_Y_ij[item_changed,], upper.tri_N_ij[item_changed,], upper.tri_P_ij[item_changed,], log=T)) + 
    sum(dbinom(upper.tri_Y_ij[,item_changed], upper.tri_N_ij[,item_changed], upper.tri_P_ij[,item_changed], log=T)) 
  
  #update P_NbyN
  P_ij_prime = upper.tri_P_ij
  for(nodes in 1:n){
    P_ij_prime[item_changed,nodes]<- theta[z_prime[item_changed],z_prime[nodes]]
    P_ij_prime[nodes,item_changed]<- theta[z_prime[nodes],z_prime[item_changed]]
  }
  
  upper.tri_P_ij_prime = P_ij_prime*upper.tri(P_ij_prime)
  #compute the likelihood of the same points with the new assignment
  A_plus = sum(dbinom(upper.tri_Y_ij[item_changed,], upper.tri_N_ij[item_changed,],
                      upper.tri_P_ij_prime[item_changed,], log=T)) + 
    sum(dbinom(upper.tri_Y_ij[,item_changed], upper.tri_N_ij[,item_changed], 
               upper.tri_P_ij_prime[,item_changed], log=T)) 
  
  #Updating the likelihood
  A_prime = sum(A_cur) - A_minus + A_plus
  
  return(list(A_cur = sum(A_cur), A_prime = sum(A_prime)))
}


new_version = function(Y_ij, ybar,mbar,lamdabar,n_minus_y1, coef1, 
                       theta,z_current,z_prime,K, item_changed, k_cur,k_prime){
 
  A_current = sum(lamdabar + ybar * log(theta)+ (mbar)*log(1 - theta))
  
  z_P_current = vec2mat_0_P(z_current,theta)
  z_P_prime =  vec2mat_0_P(z_prime,theta)
  upper.tri.Y_ij = Y_ij*upper.tri(Y_ij)
  #Kx1 vector
  i_th_victoriesvs_all_clusters = t(z_P_current)%*%(upper.tri.Y_ij)[item_changed,]
  all_clusters_victories_vs_ith = (upper.tri.Y_ij)[,item_changed]%*%z_P_current
  
  subtract_matrix = matrix(0,K,K)
  subtract_matrix[k_cur,]= i_th_victoriesvs_all_clusters
  subtract_matrix[,k_cur]= subtract_matrix[,k_cur] + all_clusters_victories_vs_ith
  
  
  i_th_games_vs_all_clusters = t(z_P_current)%*%(n_minus_y1)[item_changed,]
  all_clusters_games_vs_ith = (n_minus_y1)[,item_changed]%*%z_P_current
  
  
  subtract_matrix_N_ij = matrix(0,K,K)
  subtract_matrix_N_ij[k_cur,]= i_th_games_vs_all_clusters
  subtract_matrix_N_ij[,k_cur] = subtract_matrix_N_ij[,k_cur] + all_clusters_games_vs_ith
  
  i_th_possibilities_vs_all_clusters = t(z_P_current)%*%(coef1)[item_changed,]
  all_clusters_possibilities_vs_ith = (coef1)[,item_changed]%*%z_P_current
  
  subtract_matrix_lamda_ij = matrix(0,K,K)
  subtract_matrix_lamda_ij[k_cur,]= i_th_possibilities_vs_all_clusters
  subtract_matrix_lamda_ij[,k_cur] = subtract_matrix_lamda_ij[,k_cur] + all_clusters_possibilities_vs_ith
  
  A_minus = sum(subtract_matrix_lamda_ij+ subtract_matrix * log(theta)+(subtract_matrix_N_ij)*log(1 - theta))
  
  #updating ybar
  i_th_victoriesvs_all_clusters_scanning = t(z_P_prime)%*%(upper.tri.Y_ij)[item_changed,]
  all_clusters_victories_vs_ith_scanning = (upper.tri.Y_ij)[,item_changed]%*%z_P_prime
  

  to_add_matrix = matrix(0,K,K)
  to_add_matrix[k_prime,] = i_th_victoriesvs_all_clusters_scanning
  to_add_matrix[,k_prime] = to_add_matrix[,k_prime]+all_clusters_victories_vs_ith_scanning
  
  
  
  #updating mbar
  i_th_gamesvs_all_clusters_scanning = t(z_P_prime)%*%(n_minus_y1*upper.tri(n_minus_y1))[item_changed,]
  all_clusters_games_vs_ith_scanning = (n_minus_y1*upper.tri(n_minus_y1))[,item_changed]%*%z_P_prime
  
  to_add_matrix_N_ij = matrix(0,K,K)
  to_add_matrix_N_ij[k_prime,] = i_th_gamesvs_all_clusters_scanning
  to_add_matrix_N_ij[,k_prime] = to_add_matrix_N_ij[,k_prime] + all_clusters_games_vs_ith_scanning
  
  #updating lambdabar
  i_th_possibilities_all_clusters_scanning = t(z_P_prime)%*%coef1[item_changed,]
  all_clusters_possibilities_vs_ith_scanning = coef1[,item_changed]%*%z_P_prime
  
  to_add_matrix_lamda_ij = matrix(0,K,K)
  to_add_matrix_lamda_ij[k_prime,] = i_th_possibilities_all_clusters_scanning
  to_add_matrix_lamda_ij[,k_prime] = to_add_matrix_lamda_ij[,k_prime] + all_clusters_possibilities_vs_ith_scanning
  
  A_plus = sum(to_add_matrix_lamda_ij+ to_add_matrix * log(theta)+(to_add_matrix_N_ij)*log(1 - theta))

  A_prime = A_current - A_minus +A_plus

  return(list(A_cur = A_current, A_prime=A_prime))
}

my_df = data.frame(n=0, K=0,version=0, median=0)
n=100
K=3
for(n in seq(20,100,10)){
  for(K in 3:7){
    
Y_ij = matrix(rpois(n**2,4),n,n)
N_ij = matrix(max(Y_ij),n,n)
diag(N_ij) = 0 
Y_ij[lower.tri(Y_ij)] = N_ij[lower.tri(Y_ij)] - t(Y_ij)[lower.tri(Y_ij)]
diag(Y_ij) = 0
z_current = sample(x = c(1:K),size = n,replace = T)
z_prime = z_current
z_prime[1] <- 2
theta = matrix(runif(K**2),K,K)
theta[lower.tri(theta)] = 1 - t(theta)[lower.tri(theta)]
z_P = vec2mat_0_P(z_current, theta)

blocks_victories_vs_players = t(z_P)%*%(Y_ij*upper.tri(Y_ij))

players_victores_vs_blocks = (Y_ij*upper.tri(Y_ij))%*%z_P

ybar = t(z_P)%*%(Y_ij*upper.tri(Y_ij))%*%z_P
# number of missed victories between block p and block q
n_minus_y1 <- (N_ij-Y_ij)*upper.tri(N_ij)
# number of missed victories between block p and block q
mbar<- t(z_P)%*%n_minus_y1%*%z_P

coef1 = lchoose(N_ij, Y_ij)*upper.tri(N_ij)
lamdabar <- t(z_P)%*%(coef1)%*%z_P



naive_version(Y_ij = Y_ij, N_ij = N_ij, z_current = z_current, theta = theta,z_prime = z_prime)

current_version(Y_ij, N_ij = N_ij, z_current = z_current, theta = theta, z_prime = z_prime, item_changed = which(z_current != z_prime),n = length(z_current))

new_version(Y_ij = Y_ij, ybar = ybar,mbar = mbar,lamdabar = lamdabar,n_minus_y1 = n_minus_y1,coef1 = coef1,theta = theta, z_current = z_current,
            z_prime = z_prime,K = K,item_changed = 1,k_cur = z_current[1],k_prime = 2)


micro=microbenchmark::microbenchmark(
  naive_version(Y_ij = Y_ij, N_ij = N_ij, z_current = z_current, theta = theta,z_prime = z_prime),
  current_version(Y_ij, N_ij = N_ij, z_current = z_current, theta = theta, z_prime = z_prime, item_changed = which(z_current != z_prime),n = length(z_current)),
  new_version(Y_ij = Y_ij, ybar = ybar,mbar = mbar,lamdabar = lamdabar,n_minus_y1 = n_minus_y1,coef1 = coef1,theta = theta, z_current = z_current,
              z_prime = z_prime,K = K,item_changed = 1,k_cur = z_current[1],k_prime = 2),neval = 1000)

a = summary(micro)
a$median[c(1:3)]

my_df = rbind(my_df, data.frame(n= rep(n,3), K=rep(K,3), version=c('naive','current','new'), median = a$median[c(1:3)]))

}}
my_df =my_df[-1,]
my_df %>% ggplot(aes(n, median, color= version))+
  geom_point()
