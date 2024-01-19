



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
#                         Set up parameters for the simulation study
################################################################################






K_values <- c(6)  # Range of K values to explore
sigma_squared_values <- c(0.001,0.01)  # Range of sigma_squared values to explore
model_selection <- c(1) # Range of models to explore: 1= SST, 0 =Simple/unordererd



true_model = 'Tennis_data'

test_grid = expand_grid(K_values)

###############################################################################
# uploading data
###############################################################################

Y_ij <- read.table("/Users/lapo_santi/Desktop/Nial/MCMC_results/applications_orderstats/tennis/rawdata/Y_ij.csv",header  = F,row.names = 1,sep = ",")
N_ij <- read.table("/Users/lapo_santi/Desktop/Nial/MCMC_results/applications_orderstats/tennis/rawdata/N_ij.csv",header  = F,row.names = 1,sep = ",")

Y_ij = as.matrix(Y_ij)
N_ij = as.matrix(N_ij)

#chosing where to save the files
setwd("/Users/lapo_santi/Desktop/Nial/MCMC_results/applications_orderstats/tennis/estimates_file/")


for(iteration in 1:nrow(test_grid)){
  
  
  n = nrow(Y_ij)
  diag(N_ij)<-0
  N_iter = 30000
  K = test_grid$K_values[iteration]
  K_max =  test_grid$K_values[iteration]
  alpha_vec = rep((1:K)/sum(1:K))
  
  
  set.seed(123)
  
  
  
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
    
    sigma_squared <- 0.001
    a0=.8
    U = runif((K-1),0.5,a0)
    U_vec0 = sort(U)
    

    
    beta_params = beta_mean_var(U_vec0,rep(sigma_squared,K-1) )
    a_k = beta_params$alpha
    b_k = beta_params$beta
    P0 = matrix(0,K,K)
    for(k in 1:(K-1)){
      for(i in 1:(K-1)){
        for(j in (i+1):K){
          if((j-i)==k){
            P0[i,j]= rbeta(1,a_k[k],b_k[k])
          }
        }
      }
    }
    
    P0 = P0 +  lower.tri(P0)*(1-t(P0))
    
    # blocks on the diagonal should have values close to 0.5
    diag(P0)= rep(0.5,K) + runif(K,-0.1,0.1) 
    

    
    z0=matrix(0,n,1)
    for(item in 1:n){
      z0[item]= sample(1:K,1)
    }

    ground_truth[[chain]]=list(z = NULL,a=.8,sigma_squared=sigma_squared,U_vec = NULL,K=K,P=NULL)
    
    init[[chain]] =list(z = z0,a=a0,sigma_squared=sigma_squared, U_vec = U_vec0,K=K,P0=P0)
  }
  
  estimation_control = list(z = 1,a=0,sigma_squared=0, U_vec=1,K=0,P=1)
  
  
  hyper_params = list(K_max = K,alpha_vec =alpha_vec)
  
  chains_SST = adaptive_MCMC_orderstats(Y_ij = Y_ij, N_ij = N_ij,init = init , 
                                        estimation_control = estimation_control, 
                                        ground_truth = ground_truth, 
                                        N = n, N_iter = N_iter,n_chains = n_chains, 
                                        optimal_acceptance_rate=optimal_acceptance_rate, 
                                        hyper_params = hyper_params, seed = seed)
  my_names <- paste0("chain", 1:n_chains)
  names(chains_SST)<-my_names 
  
  filename <- paste0("True_Model",true_model,"Est_model_SST","_N", n,"_K", K, "_seed", seed,".RDS")
  saveRDS(chains_SST, file = filename) #saving results
  
  #-----------------------------------------------------------------------------
  # WST MODEL
  #-----------------------------------------------------------------------------
 
  seed=123
  print(paste0("Estimation of the WST model, K=",K))
  
  #initializing each chain
  
  ground_truth= list()
  init = list()
  for(chain in 1:n_chains){
    
    sigma_squared <- 0.01
    a0=.8
    U = runif((K-1),0.5,a0)
    U_vec0 = sort(U)
    
    
    beta_params = beta_mean_var(U_vec0,rep(sigma_squared,K-1) )
    
    a_k = beta_params$alpha
    b_k = beta_params$beta
    
    P0 = matrix(0,K,K)
    for(k in 1:(K-1)){
      for(i in 1:(K-1)){
        for(j in (i+1):K){
          if((j-i)==k){
            P0[i,j]= rbeta(1,a_k[k],b_k[k])
          }
        }
      }
    }
    
    P0 = P0 +  lower.tri(P0)*(1-t(P0))
    
    # blocks on the diagonal should have values close to 0.5
    diag(P0)= rep(0.5,K) + runif(K,-0.1,0.1) 
    

    
    z0=matrix(0,n,1)
    for(item in 1:n){
      z0[item]= sample(1:K,1)
    }
    sigma_squared <- 0.01
    ground_truth[[chain]]=list(z = NULL,a=.8,sigma_squared=sigma_squared,U_vec = NULL,K=K,P=NULL)
    
    init[[chain]] =list(z = z0,a=a0,sigma_squared=sigma_squared, U_vec = U_vec0,K=K,P0=P0)
  }
  
  estimation_control = list(z = 1,a=0,sigma_squared=0, U_vec=1,K=0,P=1)
  
  
  hyper_params = list(K_max = K,alpha_vec =alpha_vec)
  
  chains_WST = adaptive_MCMC_orderstats(Y_ij = Y_ij, N_ij = N_ij,init = init , 
                                        estimation_control = estimation_control, 
                                        ground_truth = ground_truth, 
                                        N = n, N_iter = N_iter,n_chains = n_chains, 
                                        optimal_acceptance_rate=optimal_acceptance_rate, 
                                        hyper_params = hyper_params, seed = seed)
  my_names <- paste0("chain", 1:n_chains)
  names(chains_WST)<-my_names 
  
  filename <- paste0("True_Model",true_model,"Est_model_WST","_N", n,"_K", K,"_seed", seed,".RDS")
  saveRDS(chains_WST, file = filename) #saving results
  
  #-----------------------------------------------------------------------------
  # Simple model
  #-----------------------------------------------------------------------------
  
  
  seed=123
  print(paste0("Estimation of Simple model, K=",K))
  
  
  
  estimation_control = list(z = 1,P=1)
  ground_truth= list(z = NULL,P=NULL)
  hyper_params = list(K_max = K,alpha_vec =alpha_vec)
  chains_Simple = adaptive_MCMC_UNORDERED(Y_ij, N_ij, init , estimation_control, 
                                          ground_truth,n, N_iter,n_chains, 
                                          optimal_acceptance_rate=optimal_acceptance_rate, hyper_params, seed)
  names(chains_Simple)<-my_names 
  filename_simple <- paste0("True_Model",true_model,"Est_model_Simple_","_N", n,"_K", K,  "_seed", seed,".RDS")
  saveRDS(chains_Simple, file = filename_simple) #saving results
  
}

