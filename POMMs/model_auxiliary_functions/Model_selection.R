
library(mcclust)
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
library(ggside)
library(parallel)
library(truncnorm)
library(label.switching)
library(doRNG)
source("/Users/lapo_santi/Desktop/Nial/oldmaterial/POMM_flex/functions_container_flex.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/model_auxiliary_functions/Functions_priorSST.R")
source("/Users/lapo_santi/Desktop/Nial/oldmaterial/project/simplified model/SaraWade.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/model_auxiliary_functions/Inference_orderstats.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/model_auxiliary_functions/MCMC_functions.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/Metropolis_within_Gibbs_code.R")





################################################################################
#                        step 1: upload the data
################################################################################

is.simulation=T

true_model = 'SST'
est_model = 'SST'
#where the data are stored
data_wd<- "/Users/lapo_santi/Desktop/Nial/MCMC_results/simulation_31Jan2024/SST_true/Simulated data/"


filenames <- list.files(pattern = paste0(true_model),path = data_wd)

#choose which data you want to use
K=3

file = K-2

data_to_be_estimated <- readRDS(paste0(data_wd,"/",filenames[file]))
data_to_be_estimated$data_plot

N_ij = data_to_be_estimated$N_ij
n = nrow(N_ij)
Y_ij = data_to_be_estimated$Y_ij
K = data_to_be_estimated$ground_truth$K
ground_truth =data_to_be_estimated$ground_truth

print(paste0("True data--->", filenames[file]))



################################################################################
#                 decide for how many Ks we want to compute the marginal posterior       
################################################################################

optimal_acceptance_rate =.30
N_iter= 30000


for(K_th in 3:7){

K_est = K_th
setwd(paste0("/Users/lapo_santi/Desktop/Nial/MCMC_results/simulation_model_selection/power posterior/K3true/raw/K",K_est))


n_chains = 1


chains_seeds = list(20)

################################################################################
#               step 2: run the chain for t=1
################################################################################
print(paste0("Estimation of the SST model, K=",K_est))
print(paste0("Begin cycle at:",date()))



estimation_control = list(z = 1,sigma_squared=0, mu_vec=1,K=0,P=1)

t=1
K=K_est
rm(custom_init)
chains = adaptive_MCMC_orderstats(Y_ij = Y_ij, N_ij = N_ij , 
                                  estimation_control = estimation_control, 
                                  ground_truth = ground_truth, 
                                  n = n, N_iter = N_iter,n_chains = n_chains,
                                  optimal_acceptance_rate=optimal_acceptance_rate, K = K_est,
                                  seed = chains_seeds, model = 'SST')


my_names <- paste0("chain", rep(1,n_chains))
names(chains)<-my_names 

################################################################################
# step 3: obtain point estimates for t=1 and use them as starting values for the other chains
################################################################################


#initialising z---------------
burnin = .4*N_iter
z_chain = chains$chain1$est_containers$z[,-c(1:burnin)]

#if we do not have the ground truth, use the MAP as pivotal partitioning for label switching
P_est <- apply(chains$chain1$est_containers$P[,,-c(1:burnin)], MARGIN = c(1,2), mean)
P_est <- inverse_logit_f(P_est)
P_true_upper <- upper.tri.extractor(chains$chain1$est_containers$P[,,1])
upper_tri_indices <- which(upper.tri(P_est, diag = T), arr.ind = TRUE)
P_chain = chains$chain1$est_containers$P[,,-c(1:burnin)]




#computing the likelihood of every partition
llik<- apply(z_chain, 2, function(z_chain) 
  log_lik_f_binom(N = N_ij,  
                  Y = Y_ij,
                  z = z_chain,
                  P = P_est))
llik = as.matrix(llik)
z_MAP <- z_chain[,which.max(llik)]
z_pivot = z_MAP

#if we have the ground truth, use it as pivotal partitioning for label switching
run_label_switch <- label.switching(method = "ECR" ,
                                    zpivot = z_pivot ,
                                    z = t(z_chain), 
                                    K = K_est,groundTruth = chains$chain1$ground_truth$z)


#permutations#permutationschains
permutations_z<-run_label_switch$permutations$ECR
z_init<- as.vector(run_label_switch$clusters)


#initialising P---------------
P_samples <- chains$chain1$est_containers$P[,,-c(1:burnin)]

P_permuted = array(NA, dim=c(K_est,K_est,nrow(permutations_z)))
for(i in 1: nrow(permutations_z)){
  # Permute the rows of matrix P
  P_permuted[,,i] <- P_samples[permutations_z[i,], permutations_z[i,],i]
}
P_samples <- P_permuted

P_init<- apply(P_samples, MARGIN = c(1,2), mean)



mu_vec_init <- apply(chains$chain1$est_containers$mu_vec,MARGIN = 1,FUN = mean)
P_init<- chains$chain1$est_containers$P[,,1]
mu_vec_init = chains$chain1$est_containers$mu_vec[,1]
z_init=chains$chain1$est_containers$z[,1]
################################################################################
# step 4: estimate all the other ts
################################################################################



n_ts = 50
i = seq(1:50)
t = (i/n_ts)^5


n_chains = 5

n = nrow(N_ij)

ground_truth =data_to_be_estimated$ground_truth

chains_seeds = list(20,09, 1997,24,25)
K_chain = rep(K_est, n_chains)
#initializing each chain
n_cycles = n_ts/n_chains
custom_init = list(z = z_init, P=P_init, mu_vec = mu_vec_init)
  
  
  
for(cycle in 1:n_cycles){
  t_i= t[c((1:5)+5*(cycle-1))]
  

  
  
  

  #-----------------------------------------------------------------------------
  # SST MODEL
  #-----------------------------------------------------------------------------
  

  print(paste0("Begin cycle at:",date()))
  print(paste0("Cycle ",cycle," out of ", n_cycles))
  
  
  estimation_control = list(z = 1,sigma_squared=0, mu_vec=1,K=0,P=1)
  
  
  K=K_chain

  chains = adaptive_MCMC_orderstats(Y_ij = Y_ij, N_ij = N_ij , 
                                    estimation_control = estimation_control, 
                                    ground_truth = ground_truth, 
                                    n = n, N_iter = N_iter,n_chains = n_chains,
                                    optimal_acceptance_rate=optimal_acceptance_rate, K = K_chain,
                                    seed = chains_seeds, model = 'SST',custom_init = custom_init, t = t_i)
  
  
  my_names <- paste0("chain", rep(1,n_chains))
  names(chains)<-my_names 
  for(i_to_be_saved in 1:n_chains){
  filename_SST <- paste0("True_Model",true_model,"Est_model_SST_True_K3","_N", n,"_t",t_i[i_to_be_saved],"_estK",K[i_to_be_saved],".RDS")
  saveRDS(object = chains[[i_to_be_saved]], file =filename_SST) #saving results
  }
  rm(K)
}


}







