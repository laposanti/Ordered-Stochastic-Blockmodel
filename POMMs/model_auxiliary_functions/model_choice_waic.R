

library(googledrive)
library(ggplot2)
library(dplyr)
library(label.switching)
library(foreach)

setwd("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/")
source("./model_auxiliary_functions/MCMC_functions.R")
source("./model_auxiliary_functions/Functions_priorSST.R")

subject = "lapo.santi@ucdconnect.ie"
service_account_key = "./sonic-426715-75af23aca274.json"
googledrive::drive_deauth()
googledrive::drive_auth_configure(path = "./client_secret_573831164304-jqqj3i5mhvubbkkuifvtgkfsut8lse3g.apps.googleusercontent.com.json")
googledrive::drive_auth(email = subject)
# 
# filenames <- list.files(pattern = paste0('Data_from',true_model),path = data_wd)

folder_url <- "https://drive.google.com/drive/u/1/folders/1T0mrKQgDnn2QyW-NpBp2o_apFuraC9Uy"




folder <- drive_get(as_id(folder_url))

# List all files in the folder
files_in_folder <- drive_ls(path = folder)




for(est_model in c('SST','WST', 'Simple')){
  true_model = 'Citation_data'
  files_in_folder <- list.files(pattern = paste0(true_model),path = data_wd)
  
  
  # est_model_files = grep(pattern = paste0('est_model',est_model),
  #                        x = filenames,value = T,ignore.case = F)
  #
  # print(est_model_files)
  # Define the expression or pattern to match file names
  
  # Example: "report" to match all files with "report" in the name
  
  # Filter files based on the pattern?
  
  
  # pattern <- est_model
  # matching_files <- files_in_folder[grep(pattern, files_in_folder$name), ]
  # 
  
  
  pattern <- est_model
  files_given_model = grep(pattern, files_in_folder)
  
  
  #uploaded_results<- readRDS(paste0(data_wd,"/",est_model_files[file]))
  
  data_wd = './results/MCMC_output/model_choice/WAIC_method/diag_free/citations_true///'
  
  df_model_choice = data.frame(K_true = numeric(),
                               true_model = character(),
                               est_model = character(),
                               K_est = numeric(),
                               looic = numeric(),
                               num_iter = numeric(),
                               seed = numeric(),
                               recovery_level = numeric())
  
  df_diagnostics = data.frame(K_true = numeric(),
                              true_model = character(),
                              est_model = character(),
                              K_est = numeric(),
                              NA_count = numeric(),
                              Inf_count = numeric(),
                              bad_values_count = numeric(),
                              r_effs_less100 = numeric(),
                              num_iter = numeric())
  
  for(file in files_given_model){
    
    # save_path = paste0(data_wd,matching_files$name[file])
    # drive_download(file = matching_files$id[file], path = save_path,overwrite = T)
    uploaded_results = readRDS(paste0(data_wd, files_in_folder[file]))
    
    for(chain in 1:(length(uploaded_results))){
      print(paste0("Estimating now:", uploaded_results[[chain]]$control_containers$est_model, " K=",dim(uploaded_results[[chain]]$est_containers$theta)[1]))
      
      Y_ij = uploaded_results[[chain]]$Y_ij
      N_ij = uploaded_results[[chain]]$N_ij
      n = nrow(N_ij)
      N_iter = uploaded_results[[chain]]$control_containers$N_iter
      num_samples = uploaded_results[[chain]]$control_containers$N_iter_eff
      
      z_chain = uploaded_results[[chain]]$est_containers$z
      theta_chain = uploaded_results[[chain]]$est_containers$theta
      
      
      upper_tri_indices <- matrix(upper.tri(N_ij),n,n,byrow = F)
      
      # Get the indices where the matrix elements are greater than zero
      non_zero_indices <- matrix(N_ij > 0,n,n)
      
      # Find the common indices (both upper.tri and strictly positive)
      common_indices <- upper_tri_indices*non_zero_indices
      
      # Convert back to a logical matrix
      filtering_obs <- matrix(as.logical(common_indices), nrow = n,ncol = n)
      
      
      upper.tri.Y_ij = Y_ij[filtering_obs]
      upper.tri.N_ij = N_ij[filtering_obs]
      
      Y_pred = matrix(NA, nrow = num_samples,ncol = length(upper.tri.Y_ij))
      
      LL_list <- foreach(i=1, .packages='foreach')%do%{
        
        K = nrow(theta_chain[,,1])
        llik = matrix(NA,  nrow = num_samples, ncol = length(upper.tri.Y_ij))
        
        for(t in 1:num_samples){
          
          z_chain_mat = vec2mat_0_P(z_chain[,t], K)
          
          P_entry = inverse_logit_f(theta_chain[,,t])
          
          P_ij = calculate_victory_probabilities(z_mat =z_chain_mat, P = P_entry)
          
          llik[t,] =  dbinom(x = upper.tri.Y_ij, 
                             size = upper.tri.N_ij,
                             prob =  P_ij[filtering_obs],log = T) 
          
          Y_pred[t+(i-1)*num_samples,] <- rbinom(length(upper.tri.Y_ij), upper.tri.N_ij, P_ij[filtering_obs])
        }
        print(i)
        return(llik)
      }
      
      LLik_sum <- lapply(LL_list,FUN = rowSums)
      
      c_log_lik = rbind(LL_list[[1]])
      
      r_effs = loo::relative_eff(c_log_lik,c(rep(1,num_samples)))
      
      
      loo_model_fit = loo::loo(c_log_lik,cores = 3,save_psis = T,r_eff = r_effs,is_method = 'psis')
      
      # saveRDS(loo_model_fit,file = paste0(data_wd,'loo_fit',
      #                                     uploaded_results[[i]]$control_containers$est_model,
      #                                     nrow(uploaded_results[[i]]$est_containers$theta[,,1]),
      #                                     ".rds"))
      # # 
      influence_values = loo::pareto_k_influence_values(loo_model_fit)
      NA_count = sum(is.na(influence_values))
      NA_position = which(is.na(influence_values))
      INF_count = sum(influence_values==Inf,na.rm = T)
      INF_position = which(influence_values==Inf)
      # 
      r_eff_count = sum(loo_model_fit$diagnostics$n_eff<100)
      bad_values_count = sum(influence_values>0.7)
      
      
      
      if(is.na(uploaded_results[[chain]]$ground_truth[1]) == T){
        recovery_level = NA
        K_true = NA
        true_model = pattern
      }else{
        recovery_level = uploaded_results$recovery_level
        K_true = uploaded_results[[i]]$ground_truth$K
        true_model = uploaded_results[[i]]$ground_truth$model
      }
      
      df_model_choice = rbind(df_model_choice, data.frame(K_tru = K_true,
                                                          true_model = true_model,
                                                          est_model = uploaded_results[[chain]]$control_containers$est_model,
                                                          K_est = K,
                                                          looic =loo_model_fit$looic,
                                                          num_iter = num_samples,
                                                          recovery_capability = recovery_level,
                                                          seed = uploaded_results[[chain]]$seed,
                                                          simulation = file))
      
      df_diagnostics= rbind(df_diagnostics, data.frame(K_tru = K_true,
                                                       true_model = true_model,
                                                       est_model = uploaded_results[[chain]]$control_containers$est_model,
                                                       K_est = K,
                                                       bad_values_count = bad_values_count,
                                                       NA_count = NA_count,
                                                       Inf_count = INF_count,
                                                       r_effs_less100 = r_eff_count,
                                                       num_iter = num_samples,
                                                       recovery_capability= recovery_level,
                                                       seed  = uploaded_results[[chain]]$seed,
                                                       simulation = file))
    }
  }
}

# Get the seed for each model when K = 3
df_model_choice <- read.csv("./results/MCMC_output/model_choice/WAIC_method/diag_free/citations_true//model_choice1.csv")
df_model_choice = df_model_choice %>%
  filter(num_iter == 42000)%>%
  mutate(seedK3 = rep(df_model_choice$seed[seq(1,130,13)],13))

df_model_choice%>%
  ggplot(aes(K_est,looic,color = est_model,shape = factor(seedK3)))+
  geom_point(alpha=0.5)+
  labs(title = "Looic for different models",
       subtitle = 'Lower values are better',
       shape = 'Seed',
       x = 'K',
       color = 'Fitted Model',
       caption = paste0("True data ~ ",df_model_choice$true_model[1]))+
  theme_minimal()

df_model_choice %>%
  group_by(simulation) %>%
  slice(which.min(looic))

print(df_model_choice[which.min(df_model_choice$looic),])

write.csv(df_model_choice, "./results/MCMC_output/model_choice/WAIC_method/diag_free/citations_true//model_choice1.csv")
write.csv(df_diagnostics, "./results/MCMC_output/model_choice/WAIC_method/diag_free/citations_true//model_diagnostics1.csv")




marginal_lik_for_bridge = function(Y_ij, N_ij, K, theta, z){
  N_iter =ncol(z)
  n = nrow(Y_ij)
  A=rep(1,K)
  log_container = matrix(NA, N_iter, length(Y_ij[upper.tri(Y_ij)]))
  for(t in 1:N_iter){
    log_matrix = matrix(0, n,n)
    for(i in 1:n){
      for(j in 1:n){
        if(N_ij[i,j] !=0){
          theta_t = theta[,,t]
          P_t = inverse_logit_f(theta_t)
          dbin_i = 0
          
          for(jj in 1:K){
            for(ii in 1:jj){
              
              vec_0 = rep(0,K)
              vec_0[ii]=1
              vec_1 = rep(0,K)
              vec_1[jj]=1
              dbin = dbinom(Y_ij[i,j],N_ij[i,j],P_t[ii,jj]) *ddirichlet_multinomial(1,K,vec_0,A)*ddirichlet_multinomial(1,K,vec_1,A)
              dbin_i = dbin_i + (dbin)
            }
          }
          
          log_matrix[i,j] <- log(dbin_i)
        }
      }
    }
    log_container[t,] = (log_matrix)[upper.tri(log_matrix)]
  }
  return(log_container)
}

marginal_lik_for_bridge1 <- function(Y_ij, N_ij, K, theta, z) {
  N_iter <- ncol(z)
  n <- nrow(Y_ij)
  A <- rep(1, K)
  
  # Precompute log container dimensions
  upper_tri_indices <- upper.tri(Y_ij)
  upper_tri_len <- sum(upper_tri_indices)
  log_container <- matrix(NA, N_iter, upper_tri_len)
  
  for (t in 1:N_iter) {
    theta_t <- theta[,,t]
    P_t <- inverse_logit_f(theta_t)
    
    # Compute all binomials for this iteration
    dbin_matrix <- matrix(0, n, n)
    
    non_zero_idx <- which(N_ij != 0, arr.ind = TRUE)
    for (idx in seq_len(nrow(non_zero_idx))) {
      i <- non_zero_idx[idx, 1]
      j <- non_zero_idx[idx, 2]
      
      # Vectorized inner operations over K
      ii_vals <- rep(1:K, times = K:1)
      jj_vals <- unlist(lapply(1:K, function(j) rep(j, j)))
      
      # Create matrices for vec_0 and vec_1
      vec_0 <- matrix(0, nrow = length(ii_vals), ncol = K)
      vec_1 <- matrix(0, nrow = length(jj_vals), ncol = K)
      vec_0[cbind(1:length(ii_vals), ii_vals)] <- 1
      vec_1[cbind(1:length(jj_vals), jj_vals)] <- 1
      
      # Precompute Dirichlet values for all pairs
      dirich_0 <- apply(vec_0, 1, function(v) ddirichlet_multinomial(1, K, v, A))
      dirich_1 <- apply(vec_1, 1, function(v) ddirichlet_multinomial(1, K, v, A))
      
      # Compute binomial likelihoods
      binom_vals <- dbinom(Y_ij[i, j], N_ij[i, j], P_t[ii_vals, jj_vals])
      
      # Combine them together
      dbin_i <- sum(binom_vals * dirich_0 * dirich_1)
      dbin_matrix[i, j] <- log(dbin_i)
    }
    
    # Extract upper triangular part
    log_container[t, ] <- dbin_matrix[upper_tri_indices]
  }
  
  return(log_container)
}

