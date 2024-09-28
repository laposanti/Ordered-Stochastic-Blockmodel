

library(googledrive)
library(ggplot2)
library(dplyr)
library(label.switching)
library(foreach)
library(loo)
# setwd("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/")
source("./model_auxiliary_functions/MCMC_functions.R")
source("./model_auxiliary_functions/Functions_priorSST.R")

subject = "lapo.santi@ucdconnect.ie"
service_account_key = "./sonic-426715-75af23aca274.json"
googledrive::drive_deauth()
googledrive::drive_auth_configure(path = "./client_secret_573831164304-jqqj3i5mhvubbkkuifvtgkfsut8lse3g.apps.googleusercontent.com.json")
googledrive::drive_auth(email = subject)
# 
# filenames <- list.files(pattern = paste0('Data_from',true_model),path = data_wd)

folder_df = data.frame(Ks = c('K3_true','K5_true','K6_true','K7_true'),
                       url = c('https://drive.google.com/drive/u/1/folders/1QLqA5DfE1LSfqq7GJqDy3u5O2K7vwc8N',
                               'https://drive.google.com/drive/u/1/folders/1WMb0GZuW1Je4crjh3Fv5TtjwhQF-k24G',
                               'https://drive.google.com/drive/u/1/folders/1XKZh_eKciM3XqEc2bQX-e_H9iqOQ5dcZ',
                               'https://drive.google.com/drive/u/1/folders/1zFH46_l_4BUmqE_2G0-ac6GtQLJ3B8lV')
)
folder_df = folder_df%>%
  mutate(diag_type = rep('diag_free',nrow(folder_df)))

folder_url <- folder_df$url[2]




folder <- drive_get(as_id(folder_url))

# List all files in the folder
files_in_folder <- drive_ls(path = folder)




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



#uploaded_results<- readRDS(paste0(da ta_wd,"/",est_model_files[file]))

data_wd = paste0('./results/MCMC_output/model_choice/WAIC_method/diag_free/',folder_df$Ks[2])

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


for(est_model in c('SST','WST', 'Simple')){
# 
#   files_in_folder <- list.files(pattern = paste0(true_model),path = data_wd)
#   
  pattern <- paste0("est_model",est_model)
  files_given_model = grep(pattern, files_in_folder$name)
  
  for(file in files_given_model){
    
    # save_path = paste0(data_wd,matching_files$name[file])
    save_path = paste0(data_wd,'/temp1.rds')
    drive_download(file = files_in_folder$id[file], path = save_path,overwrite = T)
    uploaded_results = readRDS(save_path)

    for(chain in 1:(length(uploaded_results)-1)){
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
    
      
      loo_model_fit = loo::loo(c_log_lik,cores = 1,
                               save_psis = T,
                               r_eff = r_effs,is_method = 'psis')
      
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
        recovery_level = NAx
        K_true = NA
        true_model = true_model
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
                                                          simulation =     uploaded_results[[chain]]$ground_truth$data_ref))
      
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
                                                       simulation = uploaded_results[[chain]]$ground_truth$data_ref))
    }
  }
}

df_model_choice_path = paste0(data_wd,"/model_choice1.csv")
df_diagnostic_path = paste0(data_wd,"/model_diagnostics1.csv")
write.csv(df_model_choice,df_model_choice_path )
write.csv(df_diagnostics, df_diagnostic_path)
drive_put(media = df_model_choice_path, path = folder)
drive_put(media = df_diagnostic_path, path = folder)
#df_model_choice_m = read.csv("./results/MCMC_output/model_choice/WAIC_method/diag_free/K6_true/model_choice1.csv")
# Get the seed for each model when K = 3
#df_model_choice <- read.csv("./results/MCMC_output/model_choice/WAIC_method/diag_free/citations_true//model_choice1.csv")
# df_model_choice_m = df_model_choice_m %>%
#   mutate(seedK3_neutro = rep(seq(0,8,1), nrow(df_model_choice_m)/9))%>%
#   mutate(seedK3 = seed - seedK3_neutro)
# 
# df_model_choice_m%>%
#   ggplot(aes(K_est,looic,color = est_model, group=simulation))+
#   geom_point()+
#   labs(title = "Looic for different models",
#        subtitle = 'Lower values are better',
#        shape = 'Seed',
#        x = 'K',
#        color = 'Fitted Model',
#        caption = paste0("True data ~ ",df_model_choice$true_model[1]))+
#   theme_minimal()+
#   facet_wrap(~recovery_capability)
# 
# df_model_choice_m %>%
#   group_by(simulation) %>%
#   slice_min(looic, with_ties = FALSE) %>%
#   View()
# 
# 
# 
# 
# 
# 
# df_model_choice %>%
#   group_by(simulation) %>%
#   slice(which.min(looic))
# 
# print(df_model_choice[which.min(df_model_choice$looic),])
# 

# 
# marginal_lik_for_bridge = function(Y_ij, N_ij, K, theta, z){
#   N_iter =ncol(z)
#   n = nrow(Y_ij)
#   A=rep(1,K)
#   log_container = matrix(NA, N_iter, length(Y_ij[upper.tri(Y_ij)]))
#   for(t in 1:N_iter){
#     log_matrix = matrix(0, n,n)
#     for(i in 1:n){
#       for(j in 1:n){
#         if(N_ij[i,j] !=0){
#           theta_t = theta[,,t]
#           P_t = inverse_logit_f(theta_t)
#           dbin_i = 0
#           
#           for(jj in 1:K){
#             for(ii in 1:jj){
#               
#               vec_0 = rep(0,K)
#               vec_0[ii]=1
#               vec_1 = rep(0,K)
#               vec_1[jj]=1
#               dbin = dbinom(Y_ij[i,j],N_ij[i,j],P_t[ii,jj]) *ddirichlet_multinomial(1,K,vec_0,A)*ddirichlet_multinomial(1,K,vec_1,A)
#               dbin_i = dbin_i + (dbin)
#             }
#           }
#           
#           log_matrix[i,j] <- log(dbin_i)
#         }
#       }
#     }
#     log_container[t,] = (log_matrix)[upper.tri(log_matrix)]
#   }
#   return(log_container)
# }
# 
# marginal_lik_for_bridge1 <- function(Y_ij, N_ij, K, theta, z) {
#   N_iter <- ncol(z)
#   n <- nrow(Y_ij)
#   A <- rep(1, K)
#   
#   # Precompute log container dimensions
#   upper_tri_indices <- upper.tri(Y_ij)
#   upper_tri_len <- sum(upper_tri_indices)
#   log_container <- matrix(NA, N_iter, upper_tri_len)
#   
#   for (t in 1:N_iter) {
#     theta_t <- theta[,,t]
#     P_t <- inverse_logit_f(theta_t)
#     
#     # Compute all binomials for this iteration
#     dbin_matrix <- matrix(0, n, n)
#     
#     non_zero_idx <- which(N_ij != 0, arr.ind = TRUE)
#     for (idx in seq_len(nrow(non_zero_idx))) {
#       i <- non_zero_idx[idx, 1]
#       j <- non_zero_idx[idx, 2]
#       
#       # Vectorized inner operations over K
#       ii_vals <- rep(1:K, times = K:1)
#       jj_vals <- unlist(lapply(1:K, function(j) rep(j, j)))
#       
#       # Create matrices for vec_0 and vec_1
#       vec_0 <- matrix(0, nrow = length(ii_vals), ncol = K)
#       vec_1 <- matrix(0, nrow = length(jj_vals), ncol = K)
#       vec_0[cbind(1:length(ii_vals), ii_vals)] <- 1
#       vec_1[cbind(1:length(jj_vals), jj_vals)] <- 1
#       
#       # Precompute Dirichlet values for all pairs
#       dirich_0 <- apply(vec_0, 1, function(v) ddirichlet_multinomial(1, K, v, A))
#       dirich_1 <- apply(vec_1, 1, function(v) ddirichlet_multinomial(1, K, v, A))
#       
#       # Compute binomial likelihoods
#       binom_vals <- dbinom(Y_ij[i, j], N_ij[i, j], P_t[ii_vals, jj_vals])
#       
#       # Combine them together
#       dbin_i <- sum(binom_vals * dirich_0 * dirich_1)
#       dbin_matrix[i, j] <- log(dbin_i)
#     }
#     
#     # Extract upper triangular part
#     log_container[t, ] <- dbin_matrix[upper_tri_indices]
#   }
#   
#   return(log_container)
# }
# 
