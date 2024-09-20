

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
 
folder_url <- "https://drive.google.com/drive/u/1/folders/1sBqyRFO1xSP5wKzFWrIBBnJFRzcmtWYz"




folder <- drive_get(as_id(folder_url))

# List all files in the folder
files_in_folder <- drive_ls(path = folder)




# for(est_model in c('SST','WST', 'Simple')){
# filenames <- list.files(pattern = paste0('True_Model',true_model,'Est_model_', est_model),path = data_wd)


# est_model_files = grep(pattern = paste0('est_model',est_model), 
#                        x = filenames,value = T,ignore.case = F)
# 
# print(est_model_files)
# Define the expression or pattern to match file names

# Example: "report" to match all files with "report" in the name

# Filter files based on the pattern?


pattern <- "SST" 
matching_files <- files_in_folder[grep(pattern, files_in_folder$name), ]


files_10 <- grep("Kest130", matching_files$name)

total_files = 1:length(matching_files$name)
file_to_analyse = setdiff(total_files,files_10)

#uploaded_results<- readRDS(paste0(data_wd,"/",est_model_files[file]))

data_wd = './results/MCMC_output/model_choice/WAIC_method/diag_fixed//K3_true//'

df_model_choice = data.frame(K_true = numeric(),
                             true_model = character(),
                             est_model = character(),
                             K_est = numeric(),
                             looic = numeric(),
                             num_iter = numeric())

df_diagnostics = data.frame(K_true = numeric(),
                            true_model = character(),
                            est_model = character(),
                            K_est = numeric(),
                            NA_count = numeric(),
                            Inf_count = numeric(),
                            r_effs_less100 = numeric(),
                            num_iter = numeric())

for(file in c(1:length(matching_files$name))){
  
  save_path = paste0(data_wd,matching_files$name[file])
  drive_download(file = matching_files$id[file], path = save_path,overwrite = T)
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
    
    
    loo_model_fit = loo::loo(c_log_lik,cores = 3,save_psis = T,r_eff = r_effs,is_method = 'psis')
    plot(loo_model_fit)
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
    
    
    
    
    if(is.na(uploaded_results[[chain]]$ground_truth[1]) == T){
      K_true = NA
      true_model = pattern
    }else{
      K_true = uploaded_results[[i]]$ground_truth$K
      true_model = uploaded_results[[i]]$ground_truth$model
    }
    
    df_model_choice = rbind(df_model_choice, data.frame(K_tru = K_true,
                                                        true_model = true_model,
                                                        est_model = uploaded_results[[chain]]$control_containers$est_model,
                                                        K_est = K,
                                                        looic =loo_model_fit$looic,
                                                        num_iter = num_samples,
                                                        recovery_capability = uploaded_results$recovery_level,
                                                        simulation = file))
    
    df_diagnostics= rbind(df_diagnostics, data.frame(K_tru = K_true,
                                                     true_model = true_model,
                                                     est_model = uploaded_results[[chain]]$control_containers$est_model,
                                                     K_est = K,
                                                     NA_count = NA_count,
                                                     Inf_count = INF_count,
                                                     r_effs_less100 = r_eff_count,
                                                     num_iter = num_samples,
                                                     recovery_capability= uploaded_results$recovery_level,
                                                     simulation = file))
  }
}


df_model_choice%>%
  ggplot(aes(K_est,looic,color = est_model,shape = factor(recovery_capability)))+
  geom_point(alpha=0.5)+
  labs(title = "Looic for different models",
       subtitle = 'Lower values are better',
       caption = paste0("True data ~ ",df_model_choice$true_model[1], " model ,K = ",df_model_choice$K_true[1]))+
  theme_minimal()

print(df_model_choice[which.min(df_model_choice$looic),])

write.csv(df_model_choice, "./results/MCMC_output/model_choice/WAIC_method/diag_fixed//K3_true/model_choice1.csv")
write.csv(df_diagnostics, "./results/MCMC_output/model_choice/WAIC_method/diag_fixed/K3_true//model_diagnostics1.csv")





