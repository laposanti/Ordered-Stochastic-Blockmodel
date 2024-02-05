
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
#                        Set up parameters for the simulation study
################################################################################

#chosing where to save the files

setwd("/Users/lapo_santi/Desktop/Nial/MCMC_results/simulation_31Jan2024/raw/")



K_values <- c(6)  # Range of K values to explore
sigma_values <- c(0.001,0.01)  # Range of sigma_squared values to explore
model_selection <- c(1) # Range of models to explore: 1= SST, 0 =Simple/unordererd


is.simulation=T
true_model = 'WST'

test_grid = expand_grid(K_values)


for(iteration in 1:nrow(test_grid)){
  
  ###############################################################################
  # Generating data
  ###############################################################################
  
  
  

  n = 100
  M = 12
  
  N_iter = 40000
  K = test_grid$K_values[iteration]
  K_max = test_grid$K_values[iteration]
  N_ij<- matrix(M,n,n)
  
  
  set.seed(1234)
  
  mu_vec = sort(rtruncnorm(K,a = 0, b = Inf, mean = 0,sd = 1))
  
  P_star = matrix(NA, K, K)
  if(true_model == 'SST'){
    print("Generating data from the SST model")
    P_star = matrix(NA, K, K)
    P_star[col(P_star)-row(P_star)==0] <- runif(K, min=-2, max = mu_vec[1])
    for(diag_i in 1:(K-1)){
      P_star[col(P_star)-row(P_star)==diag_i] <- runif(K-diag_i, min= mu_vec[diag_i],max = mu_vec[diag_i+1])
    }
    p_ij = inverse_logit_f(P_star)
    p_ij[lower.tri(p_ij)]<- 1-t(p_ij)[lower.tri(p_ij)]
    P_star = log(p_ij/(1-p_ij))
    P<- inverse_logit_f(P_star)
  }else if(true_model =='WST'){
    print("Generating data from the WST model")
    sigma_i= .3
    sigma_squared = sigma_i**2
    print(paste0("True sigma value = ", sigma_squared))
    P_star = matrix(NA, K, K)
    for(i in 0:(K-1)){
      P_star[col(P_star)-row(P_star)==(i)] <- truncnorm::rtruncnorm(K-i,a = 0, b = Inf,mean = mu_vec[i+1],sd = sqrt(sigma_squared))
    }
    p_ij = inverse_logit_f(P_star)
    p_ij[lower.tri(p_ij)]<- 1-t(p_ij)[lower.tri(p_ij)]
    P_star = log(p_ij/(1-p_ij))
    P<- inverse_logit_f(P_star)
  }else if(true_model == 'Simple'){
    print("Generating data from the Simple model")
    P = matrix(NA, K, K)
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
    
    ground_truth[[chain]]=list(z = z,sigma_squared=ifelse(true_model=='WST',sigma_squared,NA), mu_vec = mu_vec,K=K,P=P_star)
    
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

  # 
  # sigma_df = data.frame(iter = 1: N_iter,sigma= t(chains_SST$chain1$est_containers$sigma_squared))
  # 
  # ggplot(sigma_df, aes(x=iter,y=sigma))+
  #   geom_line(group=1)+
  #   geom_hline(yintercept = chains_SST$chain1$ground_truth$sigma_squared, col='red', linetype = 2)
  # 
  # plot(x=1:N_iter, y=chains_SST$chain1$st.deviations$tau_mu_vec)
  # 
  # mus = t(chains_SST$chain1$est_containers$mu_vec)
  # mu_df = data.frame(iter = rep(1: N_iter,K) ,
  #                    ground_truth = c(rep(chains_SST$chain1$ground_truth$mu[1],N_iter),
  #                                     rep(chains_SST$chain1$ground_truth$mu[2],N_iter),
  #                                     rep(chains_SST$chain1$ground_truth$mu[3],N_iter),
  #                                     rep(chains_SST$chain1$ground_truth$mu[4],N_iter),
  #                                     rep(chains_SST$chain1$ground_truth$mu[5],N_iter)),
  #                    mu = c(mus[,1],mus[,2],mus[,3],mus[,4],mus[,5]),
  #                    level_set = c(rep(1,N_iter), rep(2,N_iter), rep(3,N_iter), rep(4,N_iter),rep(5,N_iter)))
  # 
  # ggplot(mu_df, aes(x=iter,y=mu, group= factor(level_set,levels = c(1,2,3,4)), color= factor(level_set)))+
  #   geom_line()+
  #   geom_line(aes(y=ground_truth, x = iter), linetype=2, color='red')+
  #   facet_wrap(~level_set)
  # 
  # Ps <- as.numeric(chains_SST$chain1$est_containers$P[2,4,])  # 
  # P_df = data.frame(iter = 1: N_iter,P= Ps)
  # ggplot(P_df, aes(x=iter, y=P))+
  #   geom_line(group=1)+
  #   geom_hline(yintercept = P_star[2,4], linetype=2,col='red')
  
  
}
  
  #-----------------------------------------------------------------------------
  # WST MODEL
  #-----------------------------------------------------------------------------
  
  
  print(paste0("Estimation of the WST model, K=",K))
  print(paste0("Begin cycle at:",date()))
  #initializing each chain
  
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
    
    ground_truth[[chain]]=list(z = z,sigma_squared=sigma_squared, mu_vec = mu_vec,K=K,P=P_star)
    
    init[[chain]] =list(z = z0,sigma_squared=sigma_i_SST, mu_vec = mu_SST,K=K,P0=P_star_SST)
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
      P_star_SST[col(P_star_SST)-row(P_star_SST)==(i)] <- truncnorm::rtruncnorm(K-i,a = 0, b = 10,mean = mu_SST[i+1],sd = sigma_i_SST)
    }
    
    P_star_SST[lower.tri(P_star_SST)] = 1 - P_star_SST[upper.tri(P_star_SST)]
    
    
    
    z0=matrix(0,n,1)
    for(item in 1:n){
      z0[item]= sample(1:K,1)
    }
    
    ground_truth[[chain]]=list(z = z,sigma_squared=sigma_squared, mu_vec = mu_vec,K=K,P=P_star)
    
    init[[chain]] =list(z = z0,sigma_squared=sigma_i_SST, mu_vec = mu_SST,K=K,P0=P_star_SST)
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
  


