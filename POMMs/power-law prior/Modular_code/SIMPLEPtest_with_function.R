
library(fossil)
library(dplyr)
library(ggplot2)
source("/Users/lapo_santi/Desktop/Nial/project/POMMs/power-law prior/Modular_code/function_Z1.R")
source("/Users/lapo_santi/Desktop/Nial/project/POMMs/power-law prior/Modular_code/function_P1.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/Functions_priorSST.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/SaraWade.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/functions_byLapoSanti.R")

N_iter=50000
set.seed(30)
N=100
M= 40000
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
z_true = synth$z_true
z_current = z_true
z_mat_current = vec2mat(z_true)
p_true = synth$P_matrix
labels_available = 1:K

upper.tri.non.zero= which(n_ij_matrix > 0 & upper.tri(n_ij_matrix))
n_ij = n_ij_matrix[upper.tri.non.zero]
y_ij = y_ij_matrix[upper.tri.non.zero]




#setting containers
p_container = array(0, dim=c(K,K,N_iter))
A_container = matrix(0,nrow=N_iter,ncol=1)
C_container = matrix(0,nrow=N_iter,ncol=1)



#initializing quantities
p_current= matrix(rbeta(K**2,1,1),K,K)
p_container[,,1] = p_current

z_aux = p_current%*% t(z_mat_current)
p_nbyn = z_mat_current%*%z_aux
p_ij_current = p_nbyn[upper.tri.non.zero]

A_current =  sum(dbinom(y_ij, n_ij, p_ij_current, log = T))
C_current =  get_B(p_current,1)

A_container[1] = A_current
C_container[1] = C_current
acc.count_p = 0

#setting time tracker
pb=txtProgressBar(min=1,max=N_iter)
j=2

#diagnostic = matrix(0, nrow=4,ncol=N_iter*N)


#READY TO BOMB!



while (j < N_iter + 1) {
  setTxtProgressBar(pb, j)
  
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
  C_container[j]= C_current
  p_container[,,j] = p_current
  
  j=j+1
}

burnin_p = p_container[,,-(N_iter*0.5)]
burnin_posterior = A_container[-c(1:5000)]

ts.plot(burnin_posterior)
summary(burnin_posterior)




#acceptance rate
acc.count_p/(K**2 * N_iter)


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


