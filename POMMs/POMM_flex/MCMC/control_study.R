

source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/functions_container_flex.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/Functions_priorSST.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/SaraWade.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/MCMC/adaptive_POMM_MCMC_function.R")
library(EnvStats)
library(ggplot2)
library(dplyr)
library(truncnorm)
M=10
N = 10
N_iter =10000
N_ij = matrix(M,N,N)
alpha=.5
S=.2
K=10
beta_max=0.75
gamma_vec = rep(1/K,K)
diag0.5=T
trunc = improper_prior5(K,beta_max ,alpha,diag0.5)
P=simulating_overlapping_POMM_powerlaw_norm(K,alpha,S,trunc,beta_max,diag0.5)
z = rep(1:K, N/K)
z_mat<- vec2mat(z)
P_NbyN<-calculate_victory_probabilities(z_mat,P)

Y_ij = matrix(0,N,N)
for(i in 1:N){
  for(j in 1:N){
    Y_ij[i,j] = rbinom(1,M, P_NbyN[i,j])
  }
}


z_k = which(z == 1)
P_NbyN[z_k,2]
mean(Y_ij[z_k,2]/M)



init =list(z = z,alpha=1,S=.3,P=P)
names(init)
estimation_control = list(z = 0,alpha=1,S=1,P=0)
ground_truth= list(z = z,alpha=alpha,S=S,P=P)
hyper_params = list(K = K,beta_max =beta_max,gamma_vec = gamma_vec,diag0.5=diag0.5)





TEST = adaptive_MCMC_POMM(Yij_matrix = Y_ij,Nij_matrix = N_ij,init = init,estimation_control = estimation_control,ground_truth = ground_truth,N = N,N_iter = N_iter,targ_rate = .22,hyper_params =hyper_params ,seed = 123)

alpha_container<- TEST$est_containers$alpha
p_container<- TEST$est_containers$P
S_container <- TEST$est_containers$S
plot(density(alpha_container))
mean(alpha_container)

plot(density(S_container))
mean(S_container)
plot_P(p_container,P,5000,K)

data <- data.frame(alpha = alpha_container[-c(1:6000)], S = S_container[-c(1:6000)])



# Create scatterplot with univariate distributions
ggplot(data, aes(x = alpha, y = S)) +
  geom_point(alpha=.8) +
  geom_point(x=alpha,y=S, col='red',size=2,shape=15)+
  geom_density_2d(size = .8, alpha=.7) +
  geom_xsidedensity(fill = 'red', alpha = 0.5) +
  geom_ysidedensity(fill = 'blue', alpha = 0.5) +
  labs(title = "Bivariate Scatterplot with Densities", x = "Alpha", y = "S") +
  theme_bw()


library(ggplot2)
library(ggside)



# Create scatterplot with univariate distributions
ggplot(data, aes(x = alpha, y = S)) +
  geom_point() +
  geom_xsidedensity(aes(fill = alpha), alpha = 0.5, position = "stack") +
  geom_ysidedensity(aes(fill = S), alpha = 0.5, position = "stack") +
  theme_minimal()
  
  

