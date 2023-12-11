
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
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/order_statistics_model/MCMC_wrapper_ORDERED.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/order_statistics_model/MCMC_wrapper_UNORDERED.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/order_statistics_model/MCMC_functions.R")




################################################################################
#                        Set up parameters for the simulation study
################################################################################

#chosing where to save the files

setwd("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/order_statistics_model/first_simulation_study/")



K_values <- c(3)  # Range of K values to explore
sigma_squared_values <- c(0.001)  # Range of sigma_squared values to explore
model_selection <- c(1) # Range of models to explore: 1= SST, 0 =Simple/unordererd





test_grid = expand_grid(K_values, model_selection, sigma_squared_values)


for(iteration in 1:nrow(test_grid)){
  
  
  ###############################################################################
  # Generating data
  ###############################################################################
  
  if(test_grid$model_selection[iteration] == 1){
    model<- 'POMM'}
  else if(test_grid$model_selection[iteration] == 0){
    model<- 'Simple'
  }
  
  n=100
  M= 10
  N_iter = 40000
  K= test_grid$K_values[iteration]
  K_max = 5
  
  
  
  set.seed(123)
  
  a = .80
  U = runif((K-1),0.5,a)
  sigma_squared = test_grid$sigma_squared_values[iteration]
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
          if(test_grid$model_selection[iteration] == 1){
            #simulating ordered blocks
            P[i,j]= rbeta(1,a_k[k+1],b_k[k+1])}
          else if(test_grid$model_selection[iteration] == 0){
            #simulating unordered blocks
            P[i,j]= rbeta(1,1,1)}
        }
        
      }
    }
  }
  
  
  P = P +  lower.tri(P)*(1-t(P))
  diag(P)<- runif(K,0.4,0.6)
  
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
  
  
  
  ##############################################################################
  # Estimation
  ##############################################################################
  
  n_chains = 4
  optimal_acceptance_rate =.22
  
  #-----------------------------------------------------------------------------
  # SST MODEL
  #-----------------------------------------------------------------------------
  
  seed=123
  print(paste0("Estimation of the SST model, K=",K))
  
  #initializing each chain
  
  ground_truth= list()
  init = list()
  for(chain in 1:n_chains){
    
    a0=.8
    U = runif((K-1),0.5,a0)
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
    diag(P0)<- rep(0.5,K)
    
    sigma_squared0= runif(1,0.001,min(U*(1-U)))
    
    z0=matrix(0,n,1)
    for(item in 1:n){
      z0[item]= sample(1:K,1)
    }
    ground_truth[[chain]]=list(z = z,a=a,sigma_squared=sigma_squared,U_vec = U_vec,K=K,P=P)
    
    init[[chain]] =list(z = z0,a=a0,sigma_squared=sigma_squared0, U_vec = U_vec0,K=K,P0=P0)
  }
  
  estimation_control = list(z = 1,a=0,sigma_squared=0, U_vec=0,K=0,P=1)

  
  hyper_params = list(K_max = K,alpha_vec =alpha_vec)
  chains_POMM = adaptive_MCMC_orderstats(Y_ij = Y_ij, N_ij = N_ij,init = init , 
                                         estimation_control = estimation_control, 
                                         ground_truth = ground_truth, 
                                         N = n, N_iter = N_iter,n_chains = n_chains, 
                                         optimal_acceptance_rate=optimal_acceptance_rate, 
                                         hyper_params = hyper_params, seed = seed)
  my_names <- paste0("chain", 1:n_chains)
  names(chains_POMM)<-my_names 
  
  filename <- paste0("True_Model",model,"Est_model_POMM_","_N", n,"_K", K, "_sigma_squared", sigma_squared, "_a", a, "_seed", seed,".RDS")
  saveRDS(chains_POMM, file = filename) #saving results
  
  #-----------------------------------------------------------------------------
  # Simple model
  #-----------------------------------------------------------------------------
  
  
  seed=123
  print(paste0("Estimation of Simple model, K=",K))
  
  
  
  estimation_control = list(z = 1,P=1)
  ground_truth= list(z = z,P=P)
  hyper_params = list(K_max = K,alpha_vec =alpha_vec)
  chains_Simple = adaptive_MCMC_UNORDERED(Y_ij, N_ij,init , estimation_control, 
                                          ground_truth,n, N_iter,n_chains, 
                                          optimal_acceptance_rate=optimal_acceptance_rate, hyper_params, seed)
  names(chains_Simple)<-my_names 
  filename_simple <- paste0("True_Model",model,"Est_model_Simple_","_N", n,"_K", K, "_sigma_squared", sigma_squared, "_a", a, "_seed", seed,".RDS")
  saveRDS(chains_Simple, file = filename_simple) #saving results
  
}

