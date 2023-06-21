

library(loo)
library(dbscan)
library(randnet)
library(fossil)
library(dplyr)
#
#assegna a variabile /Users/lapo_santi/Desktop/Nial/project/
source("/Users/lapo_santi/Desktop/Nial/project/POMMs/power-law prior/Modular_code/function_Z1.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/Functions_priorSST.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/SaraWade.R")
source("/Users/lapo_santi/Desktop/Nial/project/POMMs/power-law prior/Modular_code/functionP_pomm.R")

N_iter=60000


N=100
M= 30000
K=3
alpha=.5

beta_max= .85



gamma_vec = vector()
for(i in 1:K){
  gamma_vec = append(gamma_vec, i/(K**2))
}

set.seed(34)
synth = simulating_tournament_new_norm(N = N, alpha = alpha,
                                       beta_max = beta_max,
                                       K=K, M = M,
                                       gamma_vec = gamma_vec,
                                       n_ij_max = 6,model = 'POMM',diag0.5 = T 
)

improper_prior5(K,beta_max = beta_max,alpha = alpha)



n_ij_matrix=synth$n_ij_true
y_ij_matrix=synth$y_ij_true
p_ij_true=synth$p_ij_true


z_true=synth$z_true
P_true= synth$P_matrix

similarity_plot(y_ij_matrix, z_true,z_true)




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


#initializing quantities
sigma_prime=.15
alpha_current = 1
truncations_current <- improper_prior5(K,beta_max,alpha = alpha_current)

#generating a proposal matrix
alpha_container = matrix(0, nrow=1, ncol=N_iter)
p_container = array(0, dim=c(K,K,N_iter))
p_current = simulating_POMM_powerlaw_norm(K,alpha_current,truncations_current,beta_max)



A_container = matrix(0, nrow=1, ncol=N_iter)
B_container = matrix(0, nrow=1, ncol=N_iter)
C_container = matrix(0, nrow=1, ncol=N_iter)
acc.count_p = 0

#initializing_containers
init<-kmeans(x = y_ij_matrix,centers = K)$cluster
z_current= init
adj.rand.index(init,z_true)




n_k_current = as.vector(table(z_current))
z_mat_current = vec2mat(z_current)
#p_ij_function = calculate_victory_probabilities(z_mat_current,P_true)
aux = p_current%*%t(z_mat_current)
p_nbyn_current = z_mat_current%*%aux
p_ij_current = p_nbyn_current[upper.tri.non.zero]


similarity_plot(p_nbyn_current, z_current,z_current)



labels_available = 1:K

A_current= sum(dbinom(y_ij, n_ij, p_ij_current, log = T))
B_current=ddirichlet_multinomial(N,K,n_k = n_k_current ,my_alpha = gamma_vec)
C_current =  l_like_p_ij_normal(K = K, P_matrix = p_current,truncations = truncations_current,diag0.5 = T) + dlnorm_param(alpha_current)

z_container[,1] = z_current
p_container[,,1] = p_current
alpha_container[1] = alpha_current
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
z_current=z_true


while (j < N_iter + 1) {
  setTxtProgressBar(pb, j)
  
  
  
  p_update= P_POMM_update2( z_current = z_current,
                               p_current = p_current,
                               K = K,n_ij = n_ij,
                               y_ij = y_ij,
                               A_current = A_current, C_current = C_current,
                               upper.tri.non.zero = upper.tri.non.zero,
                               alpha_current = alpha_current,
                               beta_max = beta_max)
  
  acc.count_p = acc.count_p + p_update$acc.moves
  #updating quantities
  p_current = p_update$p_current
  C_current = p_update$C_current
  A_current = p_update$A_current
  alpha_current = p_update$alpha_current
  if(j %% 1000 == 0){
    print( paste("iteration",j,"acceptance", acc.count_p))
  }
  
  #storing results for inference
  C_container[j]= C_current
  alpha_container[j] = alpha_current
  p_container[,,j] = p_current
  j=j+1
}



ts.plot(C_container[-c(1:20000)])

acf((C_container[-c(1:40000)]))






burnin_p = p_container[,,-c(1:N_iter*0.5)]
burnin_alpha = alpha_container[-c(1:N_iter*0.5)]
ts.plot(burnin_alpha)
mean(burnin_alpha)

alpha_container[which(C_container == max(C_container))]
#acceptance rate
acc.count_p/(K * N_iter)




plots = list()
for(i in 1:K) {
  for(j in 1:K) {
    y_try = data.frame(y = as.vector(burnin_p[i, j,]))
    p1 = ggplot(y_try, aes(y)) +
      geom_density(fill = "dodgerblue", alpha = 0.5) +
      scale_x_log10() +
      geom_vline(xintercept = synth$P_matrix[i, j], color = "red")+
      xlab("probability") +
      ylab("Density") +
      ggtitle(paste("Density plot of entry ", i, ",", j, sep = ""))
    
    plots[[length(plots) + 1]] <- p1
  }
}
p_combined = patchwork::wrap_plots(plots, ncol = K, nrow = K)
p_combined



mse_table = matrix(0,K,K)
for(i in 1:K){
  for(j in 1:K){
    mse_table[i,j]= (mean(burnin_p[i,j,]) - P_true[i,j])
  }
}
mse_table
