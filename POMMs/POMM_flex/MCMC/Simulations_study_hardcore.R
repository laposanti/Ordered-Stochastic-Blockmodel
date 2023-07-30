# Load required libraries
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
N_values <- c(30,100)
K_values <- c(3,5,9)  # Range of K values to explore
overlap_values <- c(0.2,0.4,0.6)  # Range of overlap values to explore
alpha_values <- c(.5,1,1.5) # Range of alpha values to explore
switch_values <- c(1,0) # Range of models to explore: 1= POMM, 0 =Simple
M_values <- c(4000,10000,40000)
# Set other fixed parameters
N <- 100

N_iter<- 1000 #number of iterations
targ_rate <-0.22
beta_max <- 0.85

diag0.5 <- T

# Set up parallel computing
cores <- 5

# Iterate over parameter combinations using foreach


setwd("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/MCMC/results_29_07/raw_results")

test_grid = expand_grid(K_values, switch_values, N_values)

cl <- makeCluster(cores)
registerDoParallel(cl)
foreach(iterazione = 9:nrow(test_grid)) %dopar% {
  library(foreach)
  library(doParallel)
  library(tidyverse)
  library(EnvStats)
  library(truncnorm)
  library(dplyr)
  #-----------------------------------------------------------------------------
  # Generating data
  #-----------------------------------------------------------------------------
  
  M= 13  
  N = test_grid$N_values[iterazione]
  N_iter = 40000
  N_ij = matrix(M,N,N)
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
    
    P=simulating_overlapping_POMM_powerlaw_norm(K,alpha,S,trunc,beta_max,diag0.5)
    
    z = rep_len(1:K, N)
    
    z_mat=vec2mat(z)
    P_NbyN<-calculate_victory_probabilities(z_mat,P)
    
    Y_ij = matrix(0,N,N)
    for(i in 1:N){
      for(j in 1:N){
        Y_ij[i,j] = rbinom(1,M, P_NbyN[i,j])
      }
    }
  }else if(test_grid$switch_values[iterazione] == 0){
    
    #True Model Selected: Simple
    model<- 'Simple'
    
    P= matrix(.5,K,K)
    P[upper.tri(P)]<- runif(K*(K-1)/2,0.5,beta_max)
    P[lower.tri(P)]<- 1- P[upper.tri(P)]
    
    z = rep_len(1:K, N)
    z_mat=vec2mat(z)
    
    P_NbyN<-calculate_victory_probabilities(z_mat,P)
    
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
  
  
  
  chains_POMM <- list()
  for(i in 1:4){
    seed=123
    alpha0=runif(1,0.1,3)
    trunc=improper_prior5(K,beta_max,alpha = alpha0)
    S0=runif(1,0.1,.9)
    P0_POMM= simulating_overlapping_POMM_powerlaw_norm(K,alpha0,S0,trunc,beta_max,diag0.5)
    init_POMM =list(z = rep_len(sample(1:K,K,F), N),alpha=alpha0,S=S0,P=P0_POMM)

    estimation_control = list(z = 1,alpha=0,S=1,P=1)
    ground_truth= list(z = z,alpha=alpha,S=S,P=P)
    hyper_params = list(K = K,beta_max =beta_max,gamma_vec = gamma_vec,diag0.5=diag0.5)
    TEST = adaptive_MCMC_POMM(Yij_matrix = Y_ij,Nij_matrix = N_ij,init = init_POMM,
                              estimation_control = estimation_control,
                              ground_truth = ground_truth,N = N,
                              N_iter = N_iter,targ_rate = .22,
                              hyper_params =hyper_params ,seed = seed)
    chains_POMM[[i]]= TEST
  }
  filename <- paste0("True_Model",model,"Est_model_POMM","N", N,"K", K, "S", S, "_alpha", alpha,"M",M, "_seed", seed,".RDS")
  saveRDS(chains_POMM, file = filename) #saving results
  
  #------
  #Simple 
  #------

  chains_Simple <- list()
  for(i in 1:4){
    seed=123
    P0_Simple= matrix(.5,K,K)
    P0_Simple[upper.tri(P0_Simple)]<- runif(K*(K-1)/2,0.5,beta_max)
    P0_Simple[lower.tri(P0_Simple)]<- 1- P0_Simple[upper.tri(P0_Simple)]
    init_Simple =list(z = rep_len(sample(1:K,K,F), N),P=P0_Simple)
    
    estimation_control_Simple = list(z = 1,P=1)
    ground_truth_Simple= list(z = z,P=P)
    hyper_params_Simple = list(K = K,beta_max =beta_max,gamma_vec = gamma_vec,diag0.5=diag0.5)
    TEST = adaptive_MCMC_simple(Yij_matrix = Y_ij,Nij_matrix = N_ij,
                                init = init_Simple,estimation_control = estimation_control_Simple,
                                ground_truth = ground_truth_Simple,N = N,N_iter = N_iter,
                                targ_rate = .22,hyper_params =hyper_params_Simple, seed = seed)
    chains_Simple[[i]]= TEST
  }
  filename_simple <- paste0("True_Model",model,"Est_model_Simple","N", N,"K", K, "S", S, "_alpha","M",M, alpha,"_seed", seed, ".RDS")
  saveRDS(chains_Simple, file = filename_simple) #saving results
  
}


stopCluster(cl)

