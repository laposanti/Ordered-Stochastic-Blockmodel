# Load required libraries

library(doFuture)
library(progressr)
library(beepr)
library(foreach)
library(doParallel)
library(tidyverse)
library(EnvStats)
library(truncnorm)
library(dplyr)

source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/functions_container_flex.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/Functions_priorSST.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/SaraWade.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/MCMC/adaptive_POMM_MCMC_function.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/MCMC/adaptive_Simple_MCMC_function.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/MCMC/Inference_functions.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/MCMC/Simple_model_functions.R")

# Set up simulation parameters
N_values <- c(100)
K_values <- c(3,5,9)  # Range of K values to explore
S_values <- c(0.01,0.1,0.2)  # Range of overlap values to explore
alpha_values <- c(.5,1,1.5) # Range of alpha values to explore
switch_values <- c(1,0) # Range of models to explore: 1= POMM, 0 =Simple
M_values <- c(4000,10000,40000)
# Set other fixed parameters
N <- 100
targ_rate <-0.22
beta_max <- 0.85
n_chains<-4
diag0.5 <- T

# Set up parallel computing
cores <- 5

# Iterate over parameter combinations using foreach


setwd("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/MCMC/results_small_M/")

test_grid = expand_grid(K_values, switch_values, N_values)
test_grid = test_grid %>% filter(switch_values ==1)


for(iterazione in 1:nrow(test_grid)){
  #-----------------------------------------------------------------------------
  # Generating data
  #-----------------------------------------------------------------------------
  N=100
  M= 3
  N = test_grid$N_values[iterazione]
  N_ij = matrix(M,N,N)
  N_iter = 10000
  alpha=.5
  S=.01
  K= test_grid$K_values[iterazione]
  beta_max=0.8
  gamma_vec = rep(1/K,K)
  diag0.5=T
  trunc = improper_prior5(K,beta_max ,alpha,diag0.5)
  
  if(test_grid$switch_values[iterazione] == 1){
    
    #True Model selected: POMM
    model<- 'POMM'
    set.seed(123)
    P_true=simulating_overlapping_POMM_powerlaw_norm(K,alpha,S,trunc,beta_max,diag0.5)
    
    z_true = rep_len(1:K, N)
    
    z_mat=vec2mat(z_true)
    P_NbyN<-calculate_victory_probabilities(z_mat,P_true)
    
    Y_ij = matrix(0,N,N)
    for(i in 1:N){
      for(j in 1:N){
        Y_ij[i,j] = rbinom(1,M, P_NbyN[i,j])
      }
    }
  }else if(test_grid$switch_values[iterazione] == 0){
    
    #True Model Selected: Simple
    model<- 'Simple'
    set.seed(123)
    P_true= matrix(.5,K,K)
    P_true[upper.tri(P_true)]<- runif(K*(K-1)/2,0,beta_max)
    P_true[lower.tri(P_true)]<- 1- P_true[upper.tri(P_true)]
    
    z_true = rep_len(1:K, N)
    z_mat=vec2mat(z_true)
    
    P_NbyN<-calculate_victory_probabilities(z_mat,P_true)
    
    Y_ij = matrix(0,N,N)
    for(i in 1:N){
      for(j in 1:N){
        Y_ij[i,j] = rbinom(1,M, P_NbyN[i,j])
      }
    }
  }
  
  #-----------------------------------------------------------------------------
  # Estimation
  #-----------------------------------------------------------------------------
  
  
  
  #------
  #POMMM
  #------
  
  seed=123
  P_POMM = P_true
  z_POMM = z_true
  #initializing each chain
  print(paste0("Estimation of POMM model, K=",K))
  init_POMM = list()
  for(chain in 1:n_chains){
    alpha0=runif(1,0.1,3)
    
    S0=runif(1,0.1,.9)
    
    trunc=improper_prior5(K,beta_max,alpha = alpha0)
    P0_POMM= simulating_overlapping_POMM_powerlaw_norm(K,alpha0,S0,trunc,beta_max,diag0.5)
    
    z0=matrix(0,N,1)
    for(item in 1:N){
      z0[item]= sample(1:K,1)
    }
    init_POMM[[chain]] =list(z = z0,alpha=alpha0,S=S0,P=P0_POMM)
  }
  
  estimation_control = list(z = 1,alpha=1,S=1,P=1)
  ground_truth= list(z = z_true,alpha=alpha,S=S,P=P_true)
  hyper_params = list(K = K,beta_max =beta_max,gamma_vec = gamma_vec,diag0.5=diag0.5)
  chains_POMM = adaptive_MCMC_POMM(Yij_matrix = Y_ij,Nij_matrix = N_ij,init = init_POMM,
                                   estimation_control = estimation_control,
                                   ground_truth = ground_truth,N = N,n_chains = n_chains,
                                   N_iter = N_iter,targ_rate = .22,
                                   hyper_params =hyper_params ,seed = seed)
  my_names <- paste0("chain", 1:n_chains)
  names(chains_POMM)<-my_names 
  
  filename <- paste0("True_Model",model,"Est_model_POMM_","_N", N,"_K", K, "_S", S, "_alpha", alpha,"_M",M, "_seed", seed,".RDS")
  #saveRDS(chains_POMM, file = filename) #saving results
  
  #------
  #Simple 
  #------
  
  
  seed=123
  print(paste0("Estimation of Simple model, K=",K))
  init_Simple = list()
  for(chain in 1:n_chains){
    P0_Simple= matrix(.5,K,K)
    P0_Simple[upper.tri(P0_Simple)]<- runif(K*(K-1)/2,0.5,beta_max)
    P0_Simple[lower.tri(P0_Simple)]<- 1- P0_Simple[upper.tri(P0_Simple)]
    z0=vector()
    for(i in 1:N){
      z0=append(z0, sample(1:K,1))
    }
    init_Simple[[chain]]  =list(z = z0,P=P0_Simple)
  }

  
  estimation_control_Simple = list(z = 1,P=1)
  ground_truth_Simple= list(z = z_true,P=P_true)
  hyper_params_Simple = list(K = K,beta_max =beta_max,gamma_vec = gamma_vec,diag0.5=diag0.5)
  chains_Simple = adaptive_MCMC_simple(Yij_matrix = Y_ij,Nij_matrix = N_ij,
                                     init = init_Simple,estimation_control = estimation_control_Simple,
                                     ground_truth = ground_truth_Simple,N = N,N_iter = N_iter,n_chains = n_chains,
                                     targ_rate = .22,hyper_params =hyper_params_Simple, seed = seed)
  names(chains_Simple)<-my_names 
  filename_simple <- paste0("True_Model",model,"Est_model_Simple_","_N", N,"_K", K, "_S", S, "_alpha", alpha,"_M",M, "_seed", seed,".RDS")
  #saveRDS(chains_Simple, file = filename_simple) #saving results
  
}




