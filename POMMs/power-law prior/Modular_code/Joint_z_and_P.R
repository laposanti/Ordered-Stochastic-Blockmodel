library(loo)
library(dbscan)
library(randnet)
library(fossil)
library(dplyr)
source("/Users/lapo_santi/Desktop/Nial/project/POMMs/power-law prior/Modular_code/function_Z1.R")
source("/Users/lapo_santi/Desktop/Nial/project/POMMs/power-law prior/Modular_code/function_P1.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/Functions_priorSST.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/SaraWade.R")

N_iter= 10000
set.seed(32)
N=100
M= 30000
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

upper.tri.non.zero = which(n_ij_matrix > 0,arr.ind = T)
#working with just the upper triangular entries and with non-zero-values
n_ij = n_ij_matrix[upper.tri.non.zero]
y_ij = y_ij_matrix[upper.tri.non.zero]






#setting containers
z_container= matrix(0, nrow=N, ncol=N_iter)
p_container = array(0, dim=c(K,K,N_iter))
A_container = matrix(0, nrow=1, ncol=N_iter)
B_container = matrix(0, nrow=1, ncol=N_iter)
C_container = matrix(0,nrow=N_iter,ncol=1)

#initializing_containers
init<-kmeans(x = y_ij_matrix,centers = K)$cluster
z_container[,1] = init

adj.rand.index(z_container[,1], z_true)
z_current= z_container[,1]
n_k_current = as.vector(table(z_current))
z_mat_current = vec2mat(z_current)

#initializing quantities
p_current= matrix(rbeta(K**2,1,1),K,K)
p_container[,,1] = p_current

#p_ij_function = calculate_victory_probabilities(z_mat_current,P_true)
aux = p_current%*%t(z_mat_current)
p_ij_current = z_mat_current%*%aux
p_ij_current = p_ij_current[upper.tri.non.zero]
semi_similarity_checker(z_current,p_ij_current)

#similarity_plot(p_ij_current, z_current,z_current)

labels_available = 1:K

A_current= sum(dbinom(y_ij, n_ij, p_ij_current, log = T))
B_current=ddirichlet_multinomial(N,K,n_k = n_k_current ,my_alpha = gamma_vec)
C_current =  get_B(p_current,1)

A_container[1] = A_current
B_container[1]= B_current
C_container[1] = C_current



#containers for the counts of accepted proposals
acc.count_z = 0
acc.count_p = 0


#setting time tracker
pb=txtProgressBar(min=1,max=N_iter)
j=2


#READY TO BOMB!


while (j < N_iter + 1) {
  setTxtProgressBar(pb, j)
  
  #-------
  #Z update
  #-------
  z_sweep = z_update_1(z_current, A_current,B_current,y_ij,n_ij,p_current,labels_available,upper.tri.non.zero,gamma_vec,K)
  
  acc.count_z = acc.count_z+ z_sweep$acc.moves
  z_current= z_sweep$z_current
  n_k_current =z_sweep$n_k_current
  A_current=z_sweep$A_current
  B_current=z_sweep$B_current
  
  
  #-------
  #P update
  #-------
  
  p_update= P_simple_update1(z_current = z_current,
                             P_matrix = p_current,
                             K = K,n_ij = n_ij,
                             y_ij = y_ij,
                             A_current = A_current,C_current = C_current,
                             upper.tri.non.zero = upper.tri.non.zero)
  
  acc.count_p = acc.count_p + p_update$acc.count_p
  
  #updating quantities
  p_current=  p_update$p_current
  A_current = p_update$A_current
  C_current = p_update$C_current
  
  #storing results for inference
  z_container[, j] = z_current
  p_container[,,j] = p_current
  
  A_container[j] = A_current
  B_container[j] = B_current
  C_container[j]= C_current

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
similarity_matrix = pr_cc(z_container[,-c(1:N_iter*0.55)])
point_est = minVI(similarity_matrix)$cl
adj.rand.index(MAP, z_true)


similarity_plot(y_ij_matrix,synth$z_true,synth$z_true)
similarity_plot(similarity_matrix,synth$z_true,synth$z_true)

burnin_p = p_container[,,-(N_iter*0.5)]

#acceptance rate
acc.count_p/(K**2 * N_iter)

library(ggplot2)

plots = list()
for(i in 1:K) {
  for(j in 1:K) {
    y_try = data.frame(y = as.vector(burnin_p[i, j,]))
    p1 = ggplot(y_try, aes(y)) +
      geom_density(fill = "dodgerblue", alpha = 0.5) +
      scale_x_log10() +
      geom_vline(xintercept = p_true[i, j], color = "red")+
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
    mse_table[i,j]= (mean(burnin_p[i,j,]) - p_true[i,j])
  }
}

mse_table%>% pander::pander()



