



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

source("/Users/lapo_santi/Desktop/Nial/oldmaterial/project/simplified model/Functions_priorSST.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/Metropolis_within_Gibbs_code.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/model_auxiliary_functions/MCMC_functions.R")





################################################################################
#                         Set up parameters for the simulation study
################################################################################





K_values <- c(3)  # Range of K values to explore
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
setwd("/Users/lapo_santi/Desktop/Nial/MCMC_results/application_31Jan2024/raw/")


for(iteration in 1:nrow(test_grid)){
  
  
  n = nrow(Y_ij)
  diag(N_ij)<-0
  N_iter = 40000
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
  print(paste0("Begin cycle at:",date()))
  #initializing each chain
  
  ground_truth= list()
  init = list()
  for(chain in 1:n_chains){
    
    mu_SST = sort(rtruncnorm(K,a = 0, b = Inf, mean = 0,sd = 1))
    P_star_SST = matrix(NA, K, K)
    
    
    for(i in 0:(K-1)){
      P_star_SST[col(P_star_SST)-row(P_star_SST)==(i)] <- truncnorm::rtruncnorm(K-i,a = 0, b = 10,mean = mu_SST[i+1],sd = 0)
    }
    
    P_star_SST[lower.tri(P_star_SST)] = 1 - P_star_SST[upper.tri(P_star_SST)]
    
    sigma_i_SST=NA
    
    z0=matrix(0,n,1)
    for(item in 1:n){
      z0[item]= sample(1:K,1)
    }
    
    ground_truth[[chain]]=list(z = NA,sigma_squared=ifelse(true_model=='WST',sigma_squared,NA), mu_vec = NA,K=K,P=NA)
    
    init[[chain]] =list(z = z0,sigma_squared=sigma_i_SST, mu_vec = mu_SST,K=K,P0=P_star_SST)
  }
  
  estimation_control = list(z = 1,sigma_squared=0, mu_vec=1,K=0,P=1)
  
  
  hyper_params = list(K_max = K,alpha_vec =rep(1/K,K))
  
  chains_SST = adaptive_MCMC_orderstats(Y_ij = Y_ij, N_ij = N_ij,init = init , 
                                        estimation_control = estimation_control, 
                                        ground_truth = ground_truth, 
                                        N = n, N_iter = N_iter,n_chains = n_chains, 
                                        optimal_acceptance_rate=optimal_acceptance_rate, 
                                        hyper_params = hyper_params, seed = seed, model = 'SST')
  
  
  my_names <- paste0("chain", 1:n_chains)
  names(chains_SST)<-my_names 
  
  filename_SST <- paste0("True_Model",true_model,"Est_model_SST","_N", n,"_K", K, "_seed", seed,".RDS")
  saveRDS(chains_SST, file = filename_SST) #saving results
  beep("coin")
  
#-----------------------------------------------------------------------------
# WST MODEL
#-----------------------------------------------------------------------------


print(paste0("Estimation of the WST model, K=",K))
print(paste0("Begin cycle at:",date()))
#initializing each chain

init = list()
for(chain in 1:n_chains){
  
  mu_WST = sort(rtruncnorm(K,a = 0, b = Inf, mean = 0,sd = 1))
  P_star_WST = matrix(NA, K, K)
  
  sigma_i_WST= .5
  
  
  for(i in 0:(K-1)){
    P_star_WST[col(P_star_WST)-row(P_star_WST)==(i)] <- truncnorm::rtruncnorm(K-i,a = 0, b = 10,mean = mu_WST[i+1],sd = sigma_i_WST)
  }
  
  P_star_WST[lower.tri(P_star_WST)] = - t(P_star_WST)[lower.tri(P_star_WST)]
  
  z0=matrix(0,n,1)
  for(item in 1:n){
    z0[item]= sample(1:K,1)
  }
  
  ground_truth[[chain]]=list(z = NA,sigma_squared=NA, mu_vec = NA,K=K,P=NA)
  
  init[[chain]] =list(z = z0,sigma_squared=sigma_i_WST, mu_vec = mu_WST,K=K,P0=P_star_WST)
}

estimation_control = list(z = 1,sigma_squared=1, mu_vec=1,K=0,P=1)


hyper_params = list(K_max = K,alpha_vec =rep(1/K,K))


chains_WST = adaptive_MCMC_orderstats(Y_ij = Y_ij, N_ij = N_ij,init = init , 
                                      estimation_control = estimation_control, 
                                      ground_truth = ground_truth, 
                                      N = n, N_iter = N_iter,n_chains = n_chains, 
                                      optimal_acceptance_rate=optimal_acceptance_rate, 
                                      hyper_params = hyper_params, seed = seed, model = 'WST')
my_names <- paste0("chain", 1:n_chains)
names(chains_WST)<-my_names 

filename_WST <- paste0("True_Model",true_model,"Est_model_WST","_N", n,"_K", K, "_seed", seed,".RDS")
saveRDS(chains_WST, file = filename_WST) #saving results
beep("coin")
#-----------------------------------------------------------------------------
# Simple model
#-----------------------------------------------------------------------------


seed=123
print(paste0("Estimation of Simple model, K=",K))
print(paste0("Begin cycle at:",date()))


estimation_control = list(z = 1,P=1)


init = list()
for(chain in 1:n_chains){
  
  P_star_Simple = matrix(NA, K, K)
  
  for(i in 0:(K-1)){
    P_star_Simple[col(P_star_Simple)-row(P_star_Simple)==(i)] <- truncnorm::rtruncnorm(K-i,a = 0, b = 10, mean = 0,sd = .5)
  }
  
  P_star_Simple[lower.tri(P_star_Simple)] = - t(P_star_Simple)[lower.tri(P_star_Simple)]
  
  
  
  z0=matrix(0,n,1)
  for(item in 1:n){
    z0[item]= sample(1:K,1)
  }
  
  ground_truth[[chain]]=list(z = NA,sigma_squared=NA, mu_vec = NA,K=K,P=NA)
  
  init[[chain]] =list(z = z0,sigma_squared=NA, mu_vec = NA,K=K,P0=P_star_Simple)
}
estimation_control = list(z = 1,sigma_squared=0, mu_vec=0,K=0,P=1)


hyper_params = list(K_max = K,alpha_vec =rep(1/K,K))


chains_Simple = adaptive_MCMC_orderstats(Y_ij = Y_ij, N_ij = N_ij,init = init , 
                                         estimation_control = estimation_control, 
                                         ground_truth = ground_truth, 
                                         N = n, N_iter = N_iter,n_chains = n_chains, 
                                         optimal_acceptance_rate=optimal_acceptance_rate, 
                                         hyper_params = hyper_params, seed = seed, model = 'Simple')
names(chains_Simple)<-my_names 
filename_Simple <- paste0("True_Model",true_model,"Est_model_Simple","_N", n,"_K", K,  "_seed", seed,".RDS")
saveRDS(chains_Simple, file = filename_Simple) #saving results
beep("coin")
}

