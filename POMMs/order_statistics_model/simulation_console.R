
library(doFuture)
library(progressr)
library(beepr)
library(foreach)
library(doParallel)
library(tidyverse)
library(EnvStats)
library(truncnorm)
library(dplyr)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(truncnorm)

source("/Users/lapo_santi/Desktop/Nial/project/simplified model/Functions_priorSST.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/order_statistics_model/MCMC_functions.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/order_statistics_model/MCMC_wrapper.R")


# Set up simulation parameters
K_max = 5

optimal_acceptance_rate=.22
N_values <- c(100)
K_values <- c(3,5,9)  # Range of K values to explore
sigma_squared_values <- c(0.01,0.1,0.2)  # Range of overlap values to explore

a <- c(.5,1,1.5) # Range of alpha values to explore
switch_values <- c(1,0) # Range of models to explore: 1= ordered, 0 =unordered

# Set other fixed parameters
n <- 100
targ_rate <-0.22

n_chains<-4
diag0.5 <- T

# Set up parallel computing
cores <- 5

# Iterate over parameter combinations using foreach


setwd("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/order_statistics_model/first_simulation_study(z,Kfixed)/")

test_grid = expand_grid(K_values, switch_values, N_values)
test_grid = test_grid %>% filter(switch_values ==1)

iterazione=1
for(iterazione in 1:nrow(test_grid)){
  #-----------------------------------------------------------------------------
  # Generating data
  #-----------------------------------------------------------------------------
  n=100
  M= 2
  N_iter = 40000
  K= test_grid$K_values[iterazione]
  
  
  if(test_grid$switch_values[iterazione] == 1){
    
    #True Model selected: POMM
    model<- 'POMM'
    set.seed(123)
    
    a = .70
    U = runif(K,0.5,a)
    sigma_squared = 0.001
    U_vec = sort(U)
    alpha_vec = rep(1/K,K)
    
    beta_params = beta_mean_var(U_vec,rep(sigma_squared,K) )
    a_k = beta_params$alpha
    b_k = beta_params$beta
    
    N_ij<- matrix(M,n,n)
    diag(N_ij)<-0
    
    P = matrix(0,K,K)
    for(k in 0:(K-1)){
      for(i in 1:(K-1)){
        for(j in (i+1):K){
          if((j-i)==k){
            P[i,j]= rbeta(1,a_k[k+1],b_k[k+1])
          }
        }
      }
    }
    
    P = P +  lower.tri(P)*(1-t(P))
    diag(P)<- rbeta(K,1,1)
    
    #simulating z
    z = matrix(0,n,1)
    z<- sample(1:K, n,replace=T)
    
    z_P<- vec2mat_0_P(z,P)
    
    P_nbyn<- calculate_victory_probabilities(z_P,P)
    
    #simulating Y_ij
    Y_ij <- matrix(0, n,n)
    for(i in 1:n){
      for(j in 1:n){
        Y_ij[i,j]<- rbinom(1,N_ij[i,j],P_nbyn[i,j])
      }
    }
    
    Y_ij[lower.tri(Y_ij)] = N_ij[lower.tri(N_ij)] - t(Y_ij)[lower.tri(Y_ij)]
    diag(Y_ij)<- 0
    
  }else if(test_grid$switch_values[iterazione] == 0){
    
    #True Model Selected: Simple
    model<- 'Simple'
    N_ij<- matrix(M,n,n)
    diag(N_ij)<-0
    
    P = matrix(0,K,K)
    for(k in 0:(K-1)){
      for(i in 1:(K-1)){
        for(j in (i+1):K){
          if((j-i)==k){
            P[i,j]= rbeta(1,1,1)
          }
        }
      }
    }
    
    P = P +  lower.tri(P)*(1-t(P))
    diag(P)<- rbeta(K,1,1)
    
    
    #simulating z
    z = matrix(0,N,1)
    z<- sample(1:K, N,replace=T)
    
    z_P<- vec2mat_0_P(z,P)
    
    P_nbyn<- calculate_victory_probabilities(z_P,P)
    
    #simulating Y_ij
    Y_ij <- matrix(0, n,n)
    for(i in 1:n){
      for(j in 1:n){
        Y_ij[i,j]<- rbinom(1,N_ij[i,j],P_nbyn[i,j])
      }
    }
    
    Y_ij[lower.tri(Y_ij)] = N_ij[lower.tri(N_ij)] - t(Y_ij)[lower.tri(Y_ij)]
    diag(Y_ij)<- 0
    
  }
  
  #-----------------------------------------------------------------------------
  # Estimation
  #-----------------------------------------------------------------------------
  
  
  
  #------
  #POMMM
  #------
  
  seed=123
  
  #initializing each chain
  print(paste0("Estimation of POMM model, K=",K))
  init = list()
  for(chain in 1:n_chains){
    
    
    a0=runif(1,0.5,1)

    
    U = runif(K,0.5,a0)
    U_vec0 = sort(U)
    
    P0 = matrix(0,K,K)
    for(k in 0:(K-1)){
      for(i in 1:(K-1)){
        for(j in (i+1):K){
          if((j-i)==k){
            P0[i,j]= rbeta(1,a_k[k+1],b_k[k+1])
          }
        }
      }
    }
    
    P0 = P0 +  lower.tri(P0)*(1-t(P0))
    diag(P0)<- rbeta(K,1,1)
    
    sigma_squared0= runif(1,0.001,0.25)
    
    z0=matrix(0,n,1)
    for(item in 1:n){
      z0[item]= sample(1:K,1)
    }
    
    init[[chain]] =list(z = z0,a=a0,sigma_squared=sigma_squared0, U_vec = U_vec0,K=K,P0=P0)
  }
  
  estimation_control = list(z = 1,a=1,sigma_squared=1, U_vec=1,K=0,P=1)
  ground_truth= list(z = z,a=a,sigma_squared=sigma_squared,U_vec = U_vec,K=K,P=P)
  hyper_params = list(K_max = K,alpha_vec =alpha_vec)
  chains_POMM = adaptive_MCMC_orderstats(Y_ij, N_ij,init , estimation_control, 
                                         ground_truth,n, N_iter,n_chains, 
                                         optimal_acceptance_rate=optimal_acceptance_rate, hyper_params, seed)
  my_names <- paste0("chain", 1:n_chains)
  names(chains_POMM)<-my_names 
  
  filename <- paste0("True_Model",model,"Est_model_POMM_","_N", n,"_K", K, "_S", S, "_alpha", alpha,"_M",M, "_seed", seed,".RDS")
  saveRDS(chains_POMM, file = filename) #saving results
  
  #------
  #Simple 
  #------
  
  
  seed=123
  print(paste0("Estimation of Simple model, K=",K))
  init_Simple = list()
  for(chain in 1:n_chains){
    P0_Simple= matrix(.5,K,K)
    P0_Simple[upper.tri(P0_Simple)]<- runif(K*(K-1)/2,0,beta_max)
    P0_Simple[lower.tri(P0_Simple)]<- 1- P0_Simple[upper.tri(P0_Simple)]
    z0=vector()
    for(i in 1:n){
      z0=append(z0, sample(1:K,1))
    }
    init_Simple[[chain]]  =list(z = z0,P=P0_Simple)
  }
  
  
  estimation_control_Simple = list(z = 1,P=1)
  ground_truth_Simple= list(z = z_true,P=P_true)
  hyper_params_Simple = list(K_max = K_max,beta_max =beta_max,gamma_vec = gamma_vec,diag0.5=diag0.5)
  chains_Simple = adaptive_MCMC_simple(Yij_matrix = Y_ij,Nij_matrix = N_ij,
                                       init = init_Simple,estimation_control = estimation_control_Simple,
                                       ground_truth = ground_truth_Simple,N = n,N_iter = N_iter,n_chains = n_chains,
                                       targ_rate = .22,hyper_params =hyper_params_Simple, seed = seed)
  names(chains_Simple)<-my_names 
  filename_simple <- paste0("True_Model",model,"Est_model_Simple_","_N", n,"_K", K, "_S", S, "_alpha", alpha,"_M",M, "_seed", seed,".RDS")
  saveRDS(chains_Simple, file = filename_simple) #saving results
  
}

