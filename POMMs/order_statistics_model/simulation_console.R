
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



inverse_logit_f = function(x){
  y= exp(x)/(1+exp(x))
  return(y)
}

################################################################################
#                        Set up parameters for the simulation study
################################################################################

#chosing where to save the files

setwd("/Users/lapo_santi/Desktop/Nial/MCMC_results/simulation_study_orderstats/raw/")



K_values <- c(4,5,6)  # Range of K values to explore
sigma_squared_values <- c(0.001,0.01)  # Range of sigma_squared values to explore
model_selection <- c(1) # Range of models to explore: 1= SST, 0 =Simple/unordererd


is.simulation=T
true_model = 'WST'

test_grid = expand_grid(K_values)


for(iteration in 1:nrow(test_grid)){
  
  iteration=1
  ###############################################################################
  # Generating data
  ###############################################################################
  
  if(is.simulation == T)
    
    
  n = 100
  M = 12
  
  N_iter = 30000
  K = test_grid$K_values[iteration]
  K_max = test_grid$K_values[iteration]
  N_ij<- matrix(M,n,n)
  
  
  set.seed(1234)
  mu_vec = sort(rtruncnorm(K,a = 0, b = Inf, mean = 0,sd = 1))
  
  P_star = matrix(NA, K, K)
  
  if(true_model == 'SST'){
    sigma_i=0
  }else if(true_model =='WST'){
    sigma_i= .3
    print(paste0("True sigma value = ", sigma_i))
  }

  for(i in 0:(K-1)){
    P_star[col(P_star)-row(P_star)==(i)] <- truncnorm::rtruncnorm(K-i,a = 0, b = Inf,mean = mu_vec[i+1],sd = sigma_i)
  }
  
  P_star[lower.tri(P_star)] = 1 - P_star[upper.tri(P_star)]
  P<- inverse_logit_f(P_star)
 
 
  if(true_model =='Simple'){
    diag(P)<- runif(K,0.5, .9)
    P[col(P)- row(P) != 0]<- runif((K**2 - K),0.1, .3)
  }
  
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
  
  n_chains = 1
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
    
    mu_SST = sort(rtruncnorm(K,a = 0, b = Inf, mean = 0,sd = 1))
    P_star_SST = matrix(NA, K, K)
    
    sigma_i_SST=sigma_i

    
    for(i in 0:(K-1)){
      P_star_SST[col(P_star_SST)-row(P_star_SST)==(i)] <- truncnorm::rtruncnorm(K-i,a = 0, b = 10,mean = mu_SST[i+1],sd = sigma_i_SST)
    }
    
    P_star_SST[lower.tri(P_star_SST)] = 1 - P_star_SST[upper.tri(P_star_SST)]
    
    
    
    z0=matrix(0,n,1)
    for(item in 1:n){
      z0[item]= sample(1:K,1)
    }
    
    ground_truth[[chain]]=list(z = z,sigma_squared=sigma_i,mu_vec = mu_vec,K=K,P=P_star)
    
    init[[chain]] =list(z = z0,sigma_squared=sigma_i_SST, mu_vec = mu_SST,K=K,P0=P_star_SST)
  }
  
  estimation_control = list(z = 0,sigma_squared=0, mu_vec=0,K=0,P=1)
  
  
  hyper_params = list(K_max = K,alpha_vec =rep(1/K,K))
  
  chains_SST = adaptive_MCMC_orderstats(Y_ij = Y_ij, N_ij = N_ij,init = init , 
                                        estimation_control = estimation_control, 
                                        ground_truth = ground_truth, 
                                        N = n, N_iter = N_iter,n_chains = n_chains, 
                                        optimal_acceptance_rate=optimal_acceptance_rate, 
                                        hyper_params = hyper_params, seed = seed)
  my_names <- paste0("chain", 1:n_chains)
  names(chains_SST)<-my_names 
  
  filename <- paste0("True_Model",true_model,"Est_model_SST","_N", n,"_K", K, "true_sigma_squared", sigma_i, "_a", a, "_seed", seed,".RDS")
  saveRDS(chains_SST, file = filename) #saving results
  
  
  sigma_df = data.frame(iter = 1: N_iter,sigma= t(chains_SST$chain1$est_containers$sigma_squared))
  ggplot(sigma_df, aes(x=iter,y=sigma))+
    geom_line(group=1)+
    geom_hline(yintercept = chains_SST$chain1$ground_truth$sigma_squared, col='red', linetype = 2)
  
  mus = t(chains_SST$chain1$est_containers$mu_vec)
  mu_df = data.frame(iter = rep(1: N_iter,3) ,mu = c(mus[,1],mus[,2],mus[,3]), level_set = c(rep(1,N_iter), rep(2,N_iter), rep(3,N_iter)))
  
  ggplot(mu_df, aes(x=iter,y=mu, group= factor(level_set,levels = c(1,2,3)), color= factor(level_set)))+
    geom_line()+
    geom_hline(yintercept = inverse_logit_f(chains_SST$chain1$ground_truth$mu), linetype = 2)
  
  
  Ps <- as.numeric(chains_SST$chain1$est_containers$P[4,4,])
  

  
  apply(chains_SST$chain1$est_containers$P,1:2,mean)- P_star


  
  
  P_df = data.frame(iter = 1: N_iter,P= Ps)
  ggplot(P_df, aes(x=iter, y=P))+
    geom_line(group=1)+
    geom_hline(yintercept = P_star[4,4], linetype=2,col='red')
  
  #-----------------------------------------------------------------------------
  # WST MODEL
  #-----------------------------------------------------------------------------
  

  print(paste0("Estimation of the WST model, K=",K))
  
  #initializing each chain
  
  ground_truth= list()
  init = list()
  for(chain in 1:n_chains){
    
    sigma_squared= 0.01
    a0<- 0.8
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
    
    ground_truth[[chain]]=list(z = z,a=a,sigma_squared=sigma_squared,U_vec = U_vec,K=K,P=P)
    
    init[[chain]] =list(z = z0,a=a0,sigma_squared=sigma_squared, U_vec = U_vec0,K=K,P0=P0)
  }
  
  estimation_control = list(z = 1,a=0,sigma_squared=0, U_vec=1,K=0,P=1)
  
  
  hyper_params = list(K_max = K,alpha_vec = rep(1/K,K))
  
  chains_WST = adaptive_MCMC_orderstats(Y_ij = Y_ij, N_ij = N_ij,init = init , 
                                        estimation_control = estimation_control, 
                                        ground_truth = ground_truth, 
                                        N = n, N_iter = N_iter,n_chains = n_chains, 
                                        optimal_acceptance_rate=optimal_acceptance_rate, 
                                        hyper_params = hyper_params, seed = seed)
  my_names <- paste0("chain", 1:n_chains)
  names(chains_WST)<-my_names 
  
  filename <- paste0("True_Model",true_model,"Est_model_WST","_N", n,"_K", K, "true_sigma_squared", true_sigma_squared, "_a", a, "_seed", seed,".RDS")
  saveRDS(chains_WST, file = filename) #saving results
  
  #-----------------------------------------------------------------------------
  # Simple model
  #-----------------------------------------------------------------------------
  
  
  seed=123
  print(paste0("Estimation of Simple model, K=",K))
  
  
  
  estimation_control = list(z = 1,P=1)
  ground_truth= list(z = z,P=P)
  hyper_params = list(K_max = K,alpha_vec =alpha_vec)
  
  init = list()
  for(chain in 1:n_chains){
    
    P0 = matrix(0,K,K)
    for(k in 1:(K-1)){
      for(i in 1:(K-1)){
        for(j in (i+1):K){
          if((j-i)==k){
            P0[i,j]= rbeta(1,1,1)
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
    init[[chain]] =list(z = z0,K=K,P0=P0)
  }
  
  chains_Simple = adaptive_MCMC_UNORDERED(Y_ij, N_ij,init , estimation_control, 
                                          ground_truth,n, N_iter,n_chains, 
                                          optimal_acceptance_rate=optimal_acceptance_rate, hyper_params, seed)
  names(chains_Simple)<-my_names 
  filename_simple <- paste0("True_Model",true_model,"Est_model_Simple_","_N", n,"_K", K, "true_sigma_squared", true_sigma_squared, "_a", a, "_seed", seed,".RDS")
  saveRDS(chains_Simple, file = filename_simple) #saving results
  
}

