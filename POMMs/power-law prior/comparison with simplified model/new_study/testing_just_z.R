calculate_victory_probabilities <- function(z_mat, P) {
  aux <- P %*% t(z_mat)
  z_mat %*% aux
}

is.semi_symmetric=function(my_matrix){
  check = all(my_matrix*lower.tri(my_matrix) == (1 - t(my_matrix)*lower.tri(my_matrix)))
  if(check==TRUE){
    print("It is semi-simmetric")}
  else{
    diagnos_matrix= my_matrix*lower.tri(my_matrix) == (1 - t(my_matrix)*lower.tri(my_matrix))
    # diagnos_matrix = matrix(0, nrow(my_matrix),nrow(my_matrix))
    # 
    # for(i in 1:nrow(my_matrix)){
    #   for(j in i:nrow(my_matrix)){
    #     diagnos_matrix[j,i] = (my_matrix[j,i]==1-my_matrix[i,j])*1
    #   }}
  }
  return(diagnos_matrix)
}

semi_similarity_checker = function(z_current, p_nbyn_current){
  check_blox = matrix(c(1,1,2,2,3,3),ncol=2,nrow=3)
  for(p in 1:nrow(check_blox)){
    combi= matrix(0,2,nrow=1)
    first = which(z_current==check_blox[p,1])
    second = which(z_current==check_blox[p,2])
    for(i in 1:length(first)){
      for(j in 1:length(second)){
        combi = rbind(combi, c(first[i],second[j]))
      }
    }
    print(all( p_nbyn_current[cbind(combi)] == abs(1 - p_nbyn_current[cbind(combi[,2],combi[,1])])))
  }
}


spectral_clustering <- function(adj, K){
  #computing the degree vector
  d <- rowSums(adj)
  #computing the laplacian
  lapl <- diag(d) - adj
  #eigenvalues of the laplacian
  ev <- eigen(lapl, symmetric = TRUE)
  #extracting K  eigenvectors
  K_eigenvectors <- ev$vectors[ , (ncol(ev$vectors)-K+1):ncol(ev$vectors)] 
  #applying the kmeans algorithm
  ee <- kmeans(K_eigenvectors,K)$cluster
  return(ee)
}  

library(loo)
library(dbscan)
library(randnet)
library(fossil)
library(dplyr)
source("/Users/lapo_santi/Desktop/Nial/project/POMMs/power-law prior/Modular_code/function_Z1.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/Functions_priorSST.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/SaraWade.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/functions_byLapoSanti.R")

N_iter=10000
set.seed(34)

N=100
M= 10000
K=3
alpha=1






gamma_vec=c(1/5, 2/5, 3/5)
gamma_vec[1]/sum(gamma_vec)
gamma_vec[2]/sum(gamma_vec)
gamma_vec[3]/sum(gamma_vec)



beta_max=.75
synth = simulating_tournament_new(N = N, alpha = alpha,
                                           beta_max = beta_max,
                                           K=K, M = M,
                                           gamma_vec = gamma_vec,
                                           n_ij_max = 6,model = 'Simple')

n_ij_matrix=synth$n_ij_true
y_ij_matrix=synth$y_ij_true
p_ij_true=synth$p_ij_true
z_true=synth$z_true
P_true= synth$P_matrix


barplot(rowSums(y_ij_matrix)/rowSums(n_ij_matrix))
barplot(colSums(y_ij_matrix)/colSums(n_ij_matrix))



#upper.tri.non.zero= which(n_ij_matrix >= 0)
#activate for simple model
upper.tri.non.zero = which(n_ij_matrix > 0,arr.ind = T)
#working with just the upper triangular entries and with non-zero-values
n_ij = n_ij_matrix[upper.tri.non.zero]
y_ij = y_ij_matrix[upper.tri.non.zero]
p_ij = p_ij_true[upper.tri.non.zero]
# n_ij = n_ij_matrix
# y_ij = y_ij_matrix
# p_ij = p_ij_true


similarity_plot(y_ij_matrix,z_true,z_true)

#setting containers
z_container= matrix(0, nrow=N, ncol=N_iter)




A_container = matrix(0, nrow=1, ncol=N_iter)
B_container = matrix(0, nrow=1, ncol=N_iter)

#initializing_containers
init<-kmeans(x = y_ij_matrix,centers = K)$cluster
z_container[,1] = init





adj.rand.index(z_container[,1], z_true)
z_current= z_container[,1]


n_k_current = as.vector(table(z_current))
z_mat_current = vec2mat(z_current)
#p_ij_function = calculate_victory_probabilities(z_mat_current,P_true)
aux = P_true%*%t(z_mat_current)
p_nbyn_current = z_mat_current%*%aux
p_ij_current = p_nbyn_current[upper.tri.non.zero]

table(is.semi_symmetric(p_nbyn_current)[lower.tri(is.semi_symmetric(p_nbyn_current))])
semi_similarity_checker(z_current,p_nbyn_current)

similarity_plot(p_nbyn_current, z_current,z_current)



labels_available = 1:K
A_current= sum(dbinom(y_ij, n_ij, p_ij_current, log = T))
B_current=ddirichlet_multinomial(N,K,n_k = n_k_current ,my_alpha = gamma_vec)

A_container[1]=A_current
B_container[1]=B_current
#containers for the counts of accepted proposals
acc.count_z = 0
acc.count_p = 0


#setting time tracker
pb=txtProgressBar(min=1,max=N_iter)
j=2

diagnostic = matrix(0, nrow=N_iter,ncol=2)

diagnostic[1,1]=mean((p_nbyn_current-p_ij_true)**2)
diagnostic[1,2]=A_current+B_current

#READY TO BOMB!


while (j < N_iter + 1) {
  setTxtProgressBar(pb, j)
  
  
  z_sweep = z_update_1(z_current, A_current,B_current,y_ij,n_ij,P_true,labels_available = labels_available,upper.tri.non.zero = upper.tri.non.zero,gamma_vec = gamma_vec,K = K)
    
  acc.count_z = acc.count_z+ z_sweep$acc.moves
  z_current= z_sweep$z_current
  n_k_current =z_sweep$n_k_current
  A_current=z_sweep$A_current
  B_current=z_sweep$B_current

  #A_seq[j]= sum(dbinom(y_ij, n_ij, p_ij_current, log = T)) + ddirichlet_multinomial(N,K_true,n_k = n_k_current ,my_alpha = gamma_vec)
  z_container[, j] = z_current
  A_container[j] = A_current
  B_container[j] = B_current

  
  diag_mat = vec2mat(z_current)
  diag_nbyn=calculate_victory_probabilities(diag_mat,P = P_true)
  
  diagnostic[j,1]= mean((diag_nbyn-p_ij_true)**2)
  diagnostic[j,2]= A_current
  
  j=j+1
}




ts.plot(A_container[-c(1:20)])
acf((A_container[-c(1:10)]))


ts.plot(B_container[-c(1:100)])
acf((B_container[-c(1:100)]))


posterior_prop = A_container+B_container
ts.plot(posterior_prop[-c(1:100)])
MAP = z_container[, which(posterior_prop == max(posterior_prop))[1]]

acceptance_rate= acc.count_z/(N*N_iter)*100
print(acceptance_rate)

#mixing
# ts.plot(A_seq[-c(1:N_iter*0.5)], main="Traceplot of p(y_ij|z,P) p(z) p(P)",
#         xlab = "Iterations", ylab ="p(y_ij|z,P) p(z) p(P)")
# dev.off()
# acf(A_seq[-c(1:N_iter*0.5)],main="Autocorrelation plot",
#     xlab = "Lag", ylab ="ACF")
# dev.off()


#estimates
similarity_matrix = pr_cc(z_container[,-c(1:N_iter*0.25)])
point_est = minVI(similarity_matrix)$cl
adj.rand.index(point_est, z_true)


similarity_plot(y_ij_matrix,synth$z_true,synth$z_true)
similarity_plot(similarity_matrix,synth$z_true,synth$z_true)
diagnostic_df = data.frame(diagnostic[-2,])
colnames(diagnostic_df) =c("MSE_P_zi_zj", "Proportional_Posterior")

plot_df = diagnostic_df %>% filter(MSE_P_zi_zj < quantile(MSE_P_zi_zj,.95))


cor(diagnostic_df$MSE_P_zi_zj, diagnostic_df$Proportional_Posterior)

ggplot(plot_df,aes(x=MSE_P_zi_zj,y=Proportional_Posterior))+
  geom_point() +
  scale_x_discrete()

diagnostic = diagnostic[,-c(1:(N_iter*N*0.5))] 
diagnostic = t(diagnostic)
diagnostic= data.frame(diagnostic) %>% filter(X2>0)



# diagnostic[1,i+(j-2)*100]= z_sweep[[1]]
# diagnostic[2,i+(j-2)*100]= adj.rand.index(z_prime, z_sweep[[2]])
# diagnostic[3,i+(j-2)*100]= sum(dbinom(y_ij, n_ij, p_ij_z_prime, log=T)) + ddirichlet_multinomial(N,K_true,n_k = n_k_prime,my_alpha = gamma_vec)
# diagnostic[4,i+(j-2)*100]= sum((p_ij_z_prime- p_ij_true)**2)
# diagnostic[1,i+(j-2)*100]=0
# diagnostic[2,i+(j-2)*100]= adj.rand.index(z_prime, z.true)
# diagnostic[3,i+(j-2)*100]= sum(dbinom(y_ij, n_ij, p_ij_z_prime, log=T)) + ddirichlet_multinomial(N,K_true,n_k = n_k_prime,my_alpha = gamma_vec)
# diagnostic[4,i+(j-2)*100]= sum((p_ij_z_prime- p_ij_true)**2)
# 

# ------------------------------------
# LOUVAIN ALGORITHM
# ------------------------------------

# transform the adjacency matrix into an igraph object
net <- graph.adjacency(n_ij, mode=c("undirected"),diag = FALSE)
cluster_optimal(net)
# point estimate
Louv <- cluster_edge_betweenness(net)

# estimated H
length(table(Louv$membership))
# VI distance between estimated and true partition
mcclust::vi.dist(z_true,t(Louv))

# ------------------------------------
# Regularised SPECTRAL CLUSTERING
# ------------------------------------

reg_sp <- randnet::reg.SP(A = y_ij_matrix,K = 3,iter.max = 1000000000)$cluster
mcclust::vi.dist(z_true,reg_sp)

# ------------------------------------
# DBSCAN ALGORITHM
# ------------------------------------

dbscan::dbscan(y_ij)
# ------------------------------------
# K-MEANS ALGORITHM
# ------------------------------------



#measures
waic()
VI





