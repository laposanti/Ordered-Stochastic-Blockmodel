# Load required libraries
library(foreach)
library(doParallel)

# Set up simulation parameters
K_values <- c(3,5,9)  # Range of K values to explore
overlap_values <- c(0.2,0.4,0.6)    # Range of overlap values to explore
alpha_values <- c(.5,1,1.5) 
switch_values <- c(1,0)
# Range of alpha values to explore
# Set other fixed parameters
N <- 100
N_iter <- 30000
targ_rate <-0.22
beta_max <- 0.85

diag0.5 <- T

# Set up parallel computing
cores <- 5
cl <- makeCluster(cores)
registerDoParallel(cl)


# Iterate over parameter combinations using foreach

combinations = expand_grid(K_values,alpha_values,overlap_values)
setwd("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/MCMC/results")

test_grid = expand_grid(K_values, switch_values)
samples <- 5

foreach(combn = 1:nrow(test_grid)) %dopar% {
  library(truncnorm)
  library(dplyr)
  library(EnvStats)
  
  overlap <- 0.2
  alpha <- 0.5
  K <- test_grid$K_values[combn]
  gamma_vec=vector()
  for(i in 1:K){
    gamma_vec = append(gamma_vec, i/(K**2))
  }
  seed = 79141
  set.seed(seed) 
  M <- 40000 #number of games in the tournament
  N_iter<-10000 #number of iterations
  switch_on_off <- test_grid$switch_values[combn] #selects from which model to generate
  model <- ifelse(switch_on_off, 'POMM', 'Simple')
  synth <- simulating_tournament_new_overlap_norm(N = N, alpha = alpha, overlap = overlap,
                                                  beta_max = beta_max,
                                                  K = K, M = M,
                                                  gamma_vec = gamma_vec,
                                                  n_ij_max = 6, model = model, diag0.5 = TRUE
  )
  
  Nij_matrix <- synth$n_ij_true
  Yij_matrix <- synth$y_ij_true
  z_true <- synth$z_true
  p_true <- synth$P_matrix
  
  #initialising the relevant quantities
  overlap0 <- runif(1,0.2,.9)
  alpha0 <- runif(1,0.2,1.5)
  z0 <- sample(x = c(1:K),size = N, replace = T)
  init = list(overlap0 = overlap0, alpha0 = alpha0, z0 = z0)
  # Call adaptive_MCMC_POMM function and store results
  samples <- adaptive_MCMC_POMM(Yij_matrix, Nij_matrix, init, z_true, overlap, alpha, p_true, N, K, N_iter, targ_rate, beta_max, gamma_vec, diag0.5)

  filename <- paste0("True_Model",model,"Est_model_POMM","K", K, "_overlap", overlap, "_alpha", alpha, "_seed", seed,".RDS")
  #filename <- paste0("ModelPOMMPOMM","K", K, "_overlap", overlap, "_alpha", alpha, ".RDS")
  saveRDS(samples, file = filename)

  N_iter=30000
  samples_simple <- adaptive_MCMC_simple_model(Yij_matrix, Nij_matrix,init,z_true,p_true, N,K, N_iter, gamma_vec, diag0.5)
  filename_simple <- paste0("True_Model",model,"Est_model_Simple","K", K, "_overlap", overlap, "_alpha", alpha, ".RDS")
  # saveRDS(samples, file = filename_simple)



}
