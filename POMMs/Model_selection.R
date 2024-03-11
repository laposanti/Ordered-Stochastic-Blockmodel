
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
data_wd<- "/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/Data/Simulation_data/"


filenames <- list.files(pattern = paste0(true_model),path = data_wd)



#choose which data you want to use
# Choose which data you want to use


data_to_be_estimated <- readRDS(paste0(data_wd, "/", filenames[file]))
N_ij <- data_to_be_estimated$N_ij
n <- nrow(N_ij)
Y_ij <- data_to_be_estimated$Y_ij
K <- data_to_be_estimated$ground_truth$K
ground_truth <- data_to_be_estimated$ground_truth

cat(paste0("True data--->", filenames[file], "\n"))

################################################################################
# Decide for how many Ks we want to compute the marginal posterior
################################################################################

optimal_acceptance_rate <- 0.235


for(K_est in c(2:7)){
  
  setwd(paste0("/Users/lapo_santi/Desktop/Nial/MCMC_results/simulation_model_selection/power posterior/K", file + 2, "true/raw/K", K_est))
  
  N_iter <- 5000
  i <- seq(50, 0, -1.02)
  i = i[-c(22,24,26,28,30,32,34,36,38,41,43,44,45,46,48,49)]
  t <- (i / length(i)) ^ 5
  
  custom_init <- NA
  chains_seeds <- list(20)
  n_chains <- 1
  
  for (iteration in t) {
    
    cat(paste0("Estimation of the SST model, K=", K_est, "\n"))
    cat(paste0("Begin cycle at:", date(), "\n"))
    
    estimation_control <- list(z = 1, sigma_squared = 0, mu_vec = 1, K = 0, P = 1)
    
    chains <- adaptive_MCMC_orderstats(Y_ij = Y_ij, N_ij = N_ij,
                                       estimation_control = estimation_control,
                                       ground_truth = ground_truth,
                                       n = n, N_iter = N_iter, n_chains = n_chains,
                                       optimal_acceptance_rate = optimal_acceptance_rate, K = K_est,
                                       seed = chains_seeds, model = 'SST', t = iteration, custom_init = custom_init)
    
    my_names <- paste0("chain", rep(1, n_chains))
    names(chains) <- my_names 
    for (i_to_be_saved in 1:n_chains) {
      filename_SST <- paste0("True_Model", true_model, "Est_model_SST_True_K", file + 2, "_N", n, "iteration", which(t == iteration), "_estK", K_est, ".RDS")
      saveRDS(object = chains[[i_to_be_saved]], file = filename_SST) # saving results
    }
    
    custom_init <- list(z = chains$chain1$est_containers$z[, N_iter], 
                        P = chains$chain1$est_containers$P[, , N_iter],
                        mu_vec = chains$chain1$est_containers$mu_vec[, N_iter])
    
    cat(paste0(50 - which(t == iteration), " to go!!\n"))
  }
}

