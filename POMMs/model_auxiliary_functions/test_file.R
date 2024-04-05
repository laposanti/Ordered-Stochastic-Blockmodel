

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
  
  #compute the likelihood of the data with the current assignment just for i_th_turn
  A_minus = sum(dbinom(Y_ij[item_changed,], N_ij[item_changed,], P_ij[item_changed,], log=T)) + 
    sum(dbinom(Y_ij[,item_changed], N_ij[,item_changed], P_ij[,item_changed], log=T)) 
  
  #update P_NbyN
  P_ij_prime = P_ij
  for(nodes in 1:n){
    P_ij_prime[item_changed,nodes]<- theta[z_prime[item_changed],z_prime[nodes]]
    P_ij_prime[nodes,item_changed]<- theta[z_prime[nodes],z_prime[item_changed]]
  }
  #compute the likelihood of the same points with the new assignment
  A_plus = sum(dbinom(Y_ij[item_changed,], N_ij[item_changed,], P_ij_prime[item_changed,], log=T)) + 
    sum(dbinom(Y_ij[,item_changed], N_ij[,item_changed], P_ij_prime[,item_changed], log=T)) 
  
  #Updating the likelihood
  A_prime = sum(A_cur) - A_minus + A_plus
  
  return(list(A_cur = sum(A_cur), A_prime = sum(A_prime)))
}


new_version = function(Y_ij, ybar,mbar,lamdabar,n_minus_y1, coef1, 
                       theta,z_current,z_prime,K, item_changed,k_prime){
  A_current = sum(lamdabar+ ybar * log(theta)+(mbar)*log(1 - theta))
  
  z_P_current = vec2mat_0_P(z_current,theta)
  z_P_prime =  vec2mat_0_P(z_prime,theta)
  upper.tri.Y_ij = Y_ij*upper.tri(Y_ij)
  #Kx1 vector
  i_th_victoriesvs_all_clusters = t(z_P_current)%*%(upper.tri.Y_ij)[item_changed,]
  all_clusters_victories_vs_ith = (upper.tri.Y_ij)[,item_changed]%*%z_P_current
  
  subtract_matrix = matrix(0,K,K)
  subtract_matrix[item_changed,]= i_th_victoriesvs_all_clusters
  subtract_matrix[,item_changed]= subtract_matrix[,1]+ all_clusters_victories_vs_ith
  
  
  i_th_games_vs_all_clusters = t(z_P_current)%*%(n_minus_y1)[item_changed,]
  all_clusters_games_vs_ith = (n_minus_y1)[,item_changed]%*%z_P_current
  
  
  subtract_matrix_N_ij = matrix(0,K,K)
  subtract_matrix_N_ij[item_changed,]= i_th_games_vs_all_clusters
  subtract_matrix_N_ij[,item_changed] = subtract_matrix_N_ij[,1] + all_clusters_games_vs_ith
  
  i_th_possibilities_vs_all_clusters = t(z_P_current)%*%(coef1)[item_changed,]
  all_clusters_possibilities_vs_ith = (coef1)[,item_changed]%*%z_P_current
  
  subtract_matrix_lamda_ij = matrix(0,K,K)
  subtract_matrix_lamda_ij[item_changed,]= i_th_possibilities_vs_all_clusters
  subtract_matrix_lamda_ij[,item_changed] = subtract_matrix_lamda_ij[,item_changed] + all_clusters_possibilities_vs_ith
  
  A_minus = sum(subtract_matrix_lamda_ij+ subtract_matrix * log(theta)+(subtract_matrix_N_ij)*log(1 - theta))
  
  #updating ybar
  i_th_victoriesvs_all_clusters_scanning = t(z_P_prime)%*%(upper.tri.Y_ij)[1,]
  all_clusters_victories_vs_ith_scanning = (upper.tri.Y_ij)[,1]%*%z_P_prime
  

  to_add_matrix = matrix(0,K,K)
  to_add_matrix[k_prime,] = i_th_victoriesvs_all_clusters_scanning
  to_add_matrix[,k_prime] = to_add_matrix[,k_prime]+all_clusters_victories_vs_ith_scanning
  

  #updating mbar
  i_th_gamesvs_all_clusters_scanning = t(z_P_prime)%*%(n_minus_y1*upper.tri(n_minus_y1))[item_changed,]
  all_clusters_games_vs_ith_scanning = (n_minus_y1*upper.tri(n_minus_y1))[,item_changed]%*%z_P_prime
  
  to_add_matrix_N_ij = matrix(0,K,K)
  to_add_matrix_N_ij[k_prime,] = i_th_gamesvs_all_clusters_scanning
  to_add_matrix_N_ij[,k_prime] = to_add_matrix_N_ij[,k_prime] + all_clusters_games_vs_ith_scanning
  
  
  
  n_bar_scanning_prime = (mbar - subtract_matrix_N_ij + to_add_matrix_N_ij )
  
  
  i_th_possibilities_all_clusters_scanning = t(z_P_prime)%*%coef1[item_changed,]
  all_clusters_possibilities_vs_ith_scanning = coef1[,item_changed]%*%z_P_prime
  
  to_add_matrix_lamda_ij = matrix(0,K,K)
  to_add_matrix_lamda_ij[k_prime,] = i_th_possibilities_all_clusters_scanning
  to_add_matrix_lamda_ij[,k_prime] = to_add_matrix_lamda_ij[,k_prime] + all_clusters_possibilities_vs_ith_scanning
  
  A_plus = sum(to_add_matrix_lamda_ij+ to_add_matrix * log(theta)+(to_add_matrix_N_ij)*log(1 - theta))

  A_prime = A_current - A_minus +A_plus
  return(list(A_cur = A_current, A_prime=A_prime))
  }
n=5
Y_ij = matrix(rpois(36,4),6,6)
N_ij = matrix(max(Y_ij),6,6)
diag(N_ij) = 0 
Y_ij[lower.tri(Y_ij)] = N_ij[lower.tri(Y_ij)] - t(Y_ij)[lower.tri(Y_ij)]
diag(Y_ij) = 0
z_current = c(1,1,3,3,2,2) 
z_prime = c(2,1,3,3,2,2)
theta = matrix(runif(9),3,3)
theta = semi_symmetric(theta)
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
            z_prime = z_prime,K = 3,item_changed = 1,k_prime = 2)






