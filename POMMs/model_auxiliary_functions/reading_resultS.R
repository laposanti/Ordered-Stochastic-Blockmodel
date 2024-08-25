
library(ggside)
library(ggrepel)
library(igraph)
library(ggplot2)
library(abind)
library(dplyr)
library(label.switching)
library(collpcm)
library(loo)
library(gt)
library(doParallel)
library(coda)
library(mcclust)
library(LaplacesDemon)
library(cowplot)
library(tidyr)
library(googledrive)
source("/Users/lapo_santi/Desktop/Nial/oldmaterial/project/simplified model/SaraWade.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/model_auxiliary_functions/Functions_priorSST.R")

source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/model_auxiliary_functions/Inference_orderstats.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/model_auxiliary_functions/MCMC_functions.R")




#FLAG is.simulation=T IF YOU ARE READING THE RESULTS FOR A SIMULATION STUDY
is.simulation = F

if(is.simulation==F){
  
  true_model = "Tennis_data"
  
  if(true_model == "Tennis_data"){
    #where the data are stored
    setwd('/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/')
    Y_ij <-readRDS("./Data/Tennis application/Y_new.rds")
    N_ij <-readRDS("./Data/Tennis application/N_new.rds")
    median_rank  <-readRDS("./Data/Tennis application/median_rank105.rds")
    Y_ij = as.matrix(Y_ij)
    N_ij = as.matrix(N_ij)
    
    data_wd<- "/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/results/MCMC_output/Fixed_K/Application/"
    #where the data are saved
    processed_wd <- "/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/results/MCMC_output/Fixed_K/Application/Tennis_processed//"
    
    players_df = data.frame(Id = rownames(Y_ij)) 
    #now, for each game I want to filter just those players in the top one-hundred
    
    
    colnames(Y_ij) <- rownames(Y_ij)
    rownames(N_ij) <- rownames(Y_ij)
    colnames(N_ij) <- rownames(N_ij)
    
    indices <- expand.grid(row = rownames(Y_ij), col = colnames(Y_ij))
    # Convert the matrix to a data frame
    z_df_complete <- data.frame(
      row = as.character(indices$row),
      col = as.character(indices$col),
      Y = NA
    )
    
    for (i in seq_len(nrow(z_df_complete))) {
      z_df_complete$Y[i] <- Y_ij[z_df_complete$row[i], z_df_complete$col[i]]
    }
    for (i in seq_len(nrow(z_df_complete))) {
      z_df_complete$degree_pl_row[i] <- sum(Y_ij[z_df_complete$row[i],],na.rm = T)/sum(Y_ij[,z_df_complete$row[i]],na.rm = T)
    }
    for (i in seq_len(nrow(z_df_complete))) {
      z_df_complete$N[i] <- N_ij[z_df_complete$row[i], z_df_complete$col[i]]
    }
    
    
    z_df_complete$L = (z_df_complete$N -  z_df_complete$Y)
    Y_mat = matrix(z_df_complete$Y,ncol=1)
    L_mat = matrix(z_df_complete$L,ncol=1)
    
    
    BT_mat = cbind(Y_mat,L_mat)
    player2 = factor(z_df_complete$col, levels = unique(c(z_df_complete$row, z_df_complete$col)))
    player1 = factor(z_df_complete$row, levels = unique(c(z_df_complete$row, z_df_complete$col)))
    
    library(BradleyTerry2)
    
    BTM = BTm(outcome = BT_mat,player1 =player1 ,player2 =player2)
    
    abilities = BTabilities(BTM)
    bt_rank_df = data.frame(Id = rownames(abilities)[order(abilities[,1])], BTrank = nrow(Y_ij):1)
    players_df = players_df %>% inner_join(bt_rank_df, by = 'Id') %>%
      mutate(degree_pl = rowSums(Y_ij)/colSums(Y_ij))%>%
      inner_join(median_rank, by = c('Id' = 'players'))
    
    
    #     
    #     my_df = data.frame(rank = players_df$BTrank, p_top = as.numeric(log(p)), cluster = factor(est_df$est_cl))
    #     
    #     my_df %>%
    #       ggplot(aes(x = rank, y = p_top, color = cluster))+
    #       geom_point()+
    #       labs(title = 'The posterior allocation log-probability to the top block',
    #            y = 'Log p (top-block)',
    #            x = 'Bradley - Terry Ranking')
    #       theme_light()
    #     
    #     
    
    comparison = plot(players_df$last_rank, players_df$BTrank)
    
  }else if( true_model == 'Citation_data'){
    
    
    data_wd<- "/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/results/MCMC_output/Fixed_K/Application/"
    #where the data are saved
    processed_wd <- "/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/results/MCMC_output/Fixed_K/Application/processed/Citations/"
    
    scores=  read.csv("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/Data/Citations_application/journal-scores.csv")
    
    
    
    # Creating the data frame
    journal_data <- data.frame(
      Rank = 1:47,
      Journal = c("JRSS-B", "AoS", "Bka", "JASA", "Bcs", "JRSS-A", "Bern", "SJS", "Biost", "JCGS",
                  "Tech", "AmS", "JTSA", "ISR", "AISM", "CJS", "StSin", "StSci", "LDA", "JRSS-C",
                  "StMed", "ANZS", "StCmp", "StataJ", "SPL", "StNee", "Envr", "JABES", "Mtka",
                  "StMod", "JSPI", "SMMR", "BioJ", "JMA", "EES", "CSDA", "JNS", "CmpSt", "Stats",
                  "Test", "CSTM", "JSS", "JBS", "JSCS", "CSSC", "StPap", "JAS"),
      SM = c(2.09, 1.38, 1.29, 1.26, 0.85, 0.70, 0.69, 0.66, 0.66, 0.64, 0.53, 0.40, 0.37, 0.33,
             0.32, 0.30, 0.29, 0.11, 0.10, 0.09, 0.06, 0.06, 0.04, 0.02, -0.09, -0.10, -0.11,
             -0.16, -0.18, -0.22, -0.33, -0.35, -0.40, -0.45, -0.48, -0.52, -0.53, -0.64, -0.65,
             -0.70, -0.74, -0.80, -0.83, -0.92, -1.26, -1.35, -1.41),
      QSE = c(0.11, 0.07, 0.08, 0.06, 0.07, 0.19, 0.15, 0.12, 0.11, 0.12, 0.15, 0.18, 0.20, 0.25,
              0.16, 0.14, 0.09, 0.11, 0.17, 0.15, 0.07, 0.21, 0.15, 0.33, 0.09, 0.25, 0.18, 0.23,
              0.17, 0.21, 0.07, 0.16, 0.12, 0.08, 0.25, 0.07, 0.15, 0.22, 0.18, 0.15, 0.10, 0.19,
              0.16, 0.15, 0.14, 0.20, 0.15),
      SMgrouped = c(1.87, 1.17, 1.11, 1.11, 0.65, 0.31, 0.31, 0.31, 0.31, 0.31, 0.31, 0.04, 0.04,
                    0.04, 0.04, 0.04, 0.04, -0.04, -0.04, -0.04, -0.04, -0.04, -0.04, -0.04, -0.04,
                    -0.04, -0.04, -0.04, -0.04, -0.04, -0.31, -0.31, -0.31, -0.36, -0.36, -0.36,
                    -0.36, -0.36, -0.36, -0.36, -0.36, -0.36, -0.36,
                    -0.36, -0.88, -0.88, -0.88)
    )
    
    # Displaying the data frame
    print(journal_data)
    
    
    Y_ij=read.csv("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/Data/Citations_application/cross-citation-matrix.csv",header = T,row.names = 1)
    diag(Y_ij) = 0
    
    
    N_ij= matrix(0,47,47) +Y_ij*upper.tri(Y_ij)+
      t(Y_ij)*upper.tri(Y_ij)+Y_ij*lower.tri(Y_ij)+t(Y_ij)*lower.tri(Y_ij)
    rownames(Y_ij) %in% scores[,1]
    
    # Define the mapping between the acronyms in Y_ij and scores
    acronym_mapping <- c("AmS" = "AmerStatist", "AISM" = "AnnInstStatMath", "AoS" = "AnnStat", "ANZS" = "ANZJStat", 
                         "Bern" = "Bernoulli", "BioJ" = "BiometJ", "Bcs" = "Biometrics", "Bka" = "Biometrika", 
                         "Biost" = "Biostatistics", "CJS" = "CanJStat", "CSSC" = "CommStatTM", "CSTM" = "CommStatTM", 
                         "CmpSt" = "ComputatStat", "CSDA" = "CSDA", "EES" = "Econometrica", "Envr" = "Environmetrics", 
                         "ISR" = "ISR", "JABES" = "JABES", "JASA" = "JASA", "JAS" = "JAP", "JBS" = "JbiopharmStat", 
                         "JCGS" = "JCGS", "JMA" = "JMVA", "JNS" = "JRSSA", "JRSS-A" = "JRSSA", "JRSS-B" = "JRSSB", 
                         "JRSS-C" = "JRSSC", "JSCS" = "JSCS", "JSPI" = "JSPI", "JSS" = "JStatSoft", "JTSA" = "JTSA", 
                         "LDA" = "LDA", "Mtka" = "MCAP", "SJS" = "SJS", "StataJ" = "StatComput", "StCmp" = "StatMed", 
                         "Stats" = "Statistics", "StMed" = "StatMed", "SMMR" = "SMMR", "StMod" = "StatMod", "StNee" = "StatNeerlandica", 
                         "StPap" = "StatProbLett", "SPL" = "SPA", "StSci" = "StatSci", "StSin" = "StatSinica", "Tech" = "Technometrics", 
                         "Test" = "TPA")
    
    
    # Use mapvalues function to replace the entries
    scores$shortName <- plyr::mapvalues(scores$shortName, from = acronym_mapping, to = names(acronym_mapping))
    
    
    players_df = data.frame(Id = journal_data$Journal, points= journal_data$SM, median_rank = journal_data$Rank)
    
    
    #now, for each game I want to filter just those players in the top one-hundred
    
    
    colnames(Y_ij) <- rownames(Y_ij)
    rownames(N_ij) <- rownames(Y_ij)
    colnames(N_ij) <- rownames(N_ij)
    
    indices <- expand.grid(row = rownames(Y_ij), col = colnames(Y_ij))
    # Convert the matrix to a data frame
    z_df_complete <- data.frame(
      row = as.character(indices$row),
      col = as.character(indices$col),
      Y = NA
    )
    
    for (i in seq_len(nrow(z_df_complete))) {
      z_df_complete$Y[i] <- Y_ij[z_df_complete$row[i], z_df_complete$col[i]]
    }
    for (i in seq_len(nrow(z_df_complete))) {
      z_df_complete$degree_pl_row[i] <- sum(Y_ij[z_df_complete$row[i],])/sum(Y_ij[,z_df_complete$row[i]])
    }
    for (i in seq_len(nrow(z_df_complete))) {
      z_df_complete$N[i] <- N_ij[z_df_complete$row[i], z_df_complete$col[i]]
    }
    
    
    z_df_complete$L = (z_df_complete$N -  z_df_complete$Y)
    Y_mat = matrix(z_df_complete$Y,ncol=1)
    L_mat = matrix(z_df_complete$L,ncol=1)
    
    
    BT_mat = cbind(Y_mat,L_mat)
    player2 = factor(z_df_complete$col, levels = unique(c(z_df_complete$row, z_df_complete$col)))
    player1 = factor(z_df_complete$row, levels = unique(c(z_df_complete$row, z_df_complete$col)))
    library(BradleyTerry2)
    BTM = BTm(outcome = BT_mat,player1 =player1 ,player2 =player2)
    
    abilities = BTabilities(BTM)
    bt_rank_df = data.frame(Id = rownames(abilities)[order(abilities[,1])], BTrank = nrow(Y_ij):1)
    players_df = players_df %>% inner_join(bt_rank_df, by = 'Id') %>%
      mutate(degree_pl = rowSums(Y_ij)/colSums(Y_ij))
    
    
    
    
    
    
  }
  
  
  
}else if(is.simulation == T){
  true_model = "SST"
  data_wd = "/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/results/MCMC_output/Fixed_K/Simulation/"
  processed_wd <- "/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/results/MCMC_output/Fixed_K/Simulation/processed/"
}

setwd("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/")
subject = "lapo.santi@ucdconnect.ie"
service_account_key = "./sonic-426715-75af23aca274.json"
googledrive::drive_deauth()
googledrive::drive_auth_configure(path = "./client_secret_573831164304-jqqj3i5mhvubbkkuifvtgkfsut8lse3g.apps.googleusercontent.com.json")
googledrive::drive_auth(email = subject)
# 
# filenames <- list.files(pattern = paste0('Data_from',true_model),path = data_wd)
# print(filenames)
folder_url <- "https://drive.google.com/drive/u/1/folders/1V-lQDh0DCWSx57YJ1hHf7ebwd6UinY6Z"




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

true_model = 'Tennis_data'
pattern <- "Tennis_data"  # Example: "report" to match all files with "report" in the name

# Filter files based on the pattern
matching_files <- files_in_folder[grep(pattern, files_in_folder$name), ]


files_10 <- grep("Kest10", matching_files$name)
total_files = 1:length(matching_files$name)
file_to_analyse = setdiff(total_files,files_10)

check_P_posterior <- data.frame(model = character(), K= numeric(), t = numeric(), is.SST=logical(), is.WST=logical(),N_iter_eff = numeric())

for(file in file_to_analyse){
  
  #uploaded_results<- readRDS(paste0(data_wd,"/",est_model_files[file]))


  
  save_path = paste0(data_wd,matching_files$name[file])
  drive_download(file = matching_files $id[file], path = save_path,overwrite = T)
  uploaded_results = readRDS(save_path)
  

  burnin=0
  print(paste0('Now estimating ', matching_files$name[file]))
  # The file name
  filename <-  matching_files$name[file]
  
  # Extract the desired part using regular expressions
  match <- sub(".*est_model([^_]+)_.*", "\\1", filename)
  
  # Print the result
  est_model =  match
  
  print(paste0(length(matching_files$name)-file+1,' within the same class left '))
  
  N= nrow(uploaded_results$chain1$Y_ij)
  n=N
  N_iter = dim(uploaded_results$chain1$est_containers$z)[[2]]
  K = dim(uploaded_results$chain1$est_containers$theta)[[1]]
  Y_ij <- uploaded_results$chain1$Y_ij
  N_ij <- uploaded_results$chain1$N_ij
  
  
  
  #-------------------------------------------------------------------------------
  # P temporary estimate
  #-------------------------------------------------------------------------------
  theta_burned = uploaded_results$chain1$est_containers$theta[,,1:N_iter]
  z_burned =  uploaded_results$chain1$est_containers$z[,1:N_iter]
  if(est_model != 'Simple'){
    m_vec_burned = uploaded_results$chain1$est_containers$mu_vec[,1:N_iter]
  }
  if(est_model == 'WST'){
    sigma_squared_burned = uploaded_results$chain1$est_containers$sigma_squared[1:N_iter]
  }
  
  K0 <- apply(X = z_burned, MARGIN = 2, FUN = function(col) count_nonempty_clusters(col, K))
  
  K0_hat <- mean(K0)
  
  
  
  theta = apply(theta_burned, c(1,2), mean)
  P_est = inverse_logit_f(theta)
  

  
  
  my_z_est<- z_plot(z_burned = z_burned,  A = uploaded_results$chain1$control_containers$A[1:N_iter],
                    Y_ij = Y_ij, N_ij = N_ij, true_model= true_model,P_est = P_est,
                    est_model = est_model, true_value =is.simulation, 
                    diag0.5 =diag0.5 , K=  max(z_burned), N=nrow(uploaded_results$chain1$Y_ij),
                    z_true = uploaded_results$chain1$ground_truth$z ,
                    burnin =  burnin ,label_switch = T,tap= processed_wd)
  
  
  point_est_z<- as.vector(my_z_est$point_est)
  
  if(est_model == 'SST'&K==6){
    write.csv(point_est_z,paste0(processed_wd,"z_est_K_6modelSST.csv"))
  }
  
  
  table(point_est_z)
  
  
  
  
  K_est<- length(unique(point_est_z))
  permutations_z<- my_z_est$permutations
  z_chain_permuted<- my_z_est$relabeled_chain
  
  K0 <- apply(X = z_chain_permuted, MARGIN = 2, FUN = function(col) count_nonempty_clusters(col, K))
  
  K0_hat_switched <- mean(K0)
  
  
  # 
  # if(is.simulation==F){
  #   #conmputing the a-posteriori probability for each number of clusters
  #   p_clusters <- apply(z_burned, 1, function(row) {
  #     table(factor(row, levels = 1:K))
  #   })
  #   #matrix NxK where the entry n,k contains the number of times the individual n was assigned to cluster k
  #   p_clusters = t(p_clusters)
  #   
  #   z_prob = p_clusters/rowSums(p_clusters,na.rm = T)
  #   
  #   top_block_df = data.frame(items = rownames(Y_ij), p_top_block = z_prob[,1], z = point_est_z)%>%
  #     mutate(model = est_model)%>%
  #     mutate(n_clust = K)
  #   
  #   
  #   
  #   
  #   if(!exists('top_block_df_container')){
  #     top_block_df_container = top_block_df
  #   }else{
  #     top_block_df_container = rbind(top_block_df_container,top_block_df)
  #   }
  # }
  # 
  #---------------------------------------------------------------------------
  # P parameter estimate
  #-------------------------------------------------------------------------------
  
  if(est_model == 'Simple'&&true_model=='SST'){

    permutation_matrix <- permutations(n = K, r = K, v = 1:K)
    N_iter_eff = ncol(z_burned)
    for(t in 1:N_iter_eff){
      
      is.SST <- FALSE
      is.WST <- FALSE
      
      for(j in 1:nrow(permutation_matrix)){
        
        mat_permuted <- theta_burned[permutation_matrix[j,], permutation_matrix[j,], t]
        mat_permuted <- inverse_logit_f(mat_permuted)
        if(check_SST(mat_permuted)){
          is.SST <- TRUE
        }
        if(check_WST(mat_permuted)){
          is.WST <- TRUE
        }
      }
      
      # Add the results to the dataframe
      check_P_posterior <- rbind(check_P_posterior, 
                                 data.frame(model = 'Simple', 
                                            K = K, 
                                            t = t, 
                                            is.SST = is.SST, 
                                            is.WST = is.WST,
                                            N_iter_eff = N_iter_eff))
    }
  }
 
  
  
  
  
  P_s_table <- P_summary_table(P_burned = theta_burned,
                               true_value = is.simulation,
                               permutations_z = permutations_z,
                               diag0.5 = TRUE,
                               K = K, P_true = uploaded_results$chain1$ground_truth$theta,
                               burnin = burnin,
                               label_switch = T)
  
  
  P_s_table_save <-P_s_table$table
  theta_chain_permuted <- P_s_table$P_permuted
  P_est_relabeled<- inverse_logit_f(P_s_table$P_hat)
  
  
  
  #some traceplots
  upper_tri_indices= which(upper.tri(P_est_relabeled, diag=T),arr.ind = T)
  
  if(is.simulation ==T){
    theta_trace_df_post_switch <- do.call(rbind, lapply(1:(N_iter-burnin), function(j) {
      data.frame(iteration = j,
                 theta = upper.tri.extractor(theta_chain_permuted[,,j]), 
                 theta_true = upper.tri.extractor(uploaded_results$chain1$ground_truth$theta), 
                 entry = paste0(upper_tri_indices[,1], upper_tri_indices[,2]))
    }))
    theta_trace_df_post_switch= theta_trace_df_post_switch%>% 
      mutate(P = inverse_logit_f(theta))%>%
      mutate(P_true = inverse_logit_f(theta_true))
    
    traceplot_P = ggplot(theta_trace_df_post_switch, aes(x = iteration, color = entry, group=entry))+
      geom_line(aes(y=P), alpha=.3)+
      geom_line(aes(y=P_true), linetype=2, color='red')+
      facet_wrap(~entry)+
      theme_bw()
    
  }else if(is.simulation ==F){
    P_trace_df_post_switch <- do.call(rbind, lapply(1:(N_iter-burnin), function(j) {
      data.frame(iteration = j+burnin,
                 theta = upper.tri.extractor(theta_chain_permuted[,,j]), 
                 entry = paste0(upper_tri_indices[,1], upper_tri_indices[,2]))
    }))
    P_trace_df_post_switch=P_trace_df_post_switch%>% mutate(P = inverse_logit_f(theta))
    
    traceplot_P = ggplot(P_trace_df_post_switch, aes(x = iteration, color = entry, group=entry))+
      geom_line(aes(y=P), alpha=.3)+
      facet_wrap(~entry)+
      theme_bw()
  }
  
  # if(est_model == 'SST'){
  #   mu_vec_df <- do.call(rbind, lapply(1:(N_iter-burnin), function(j) {
  #     data.frame(iteration = j+burnin,
  #                mu = uploaded_results$chain1$est_containers$mu_vec[,j], 
  #                entry = factor(1:(K)))
  #   }
  #   ))
  
  # ggplot(mu_vec_df, aes(x = iteration, color = entry, group=entry))+
  #   geom_line(aes(y=mu), alpha=.3)+
  #   theme_bw()
  
  
  # mu_vec_df %>% group_by(entry) %>%
  #   summarise(mean = mean(mu),
  #             quantile5 = quantile(probs = 0.05, x = mu),
  #             quantile95 = quantile(probs = 0.95, x = mu))
  # 
  # 
  # }
  
  
  
  
  
  
  
  
  plot_name_traceplot_P <- paste0(processed_wd,"//P_traceplot",true_model,est_model,"K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
  png(plot_name_traceplot_P,width = 800, height = 800)
  par(mar = c(1.5, 1.5,1.5,1.5))
  print(traceplot_P)
  
  dev.off()
  # 
  # P_variance_df <- do.call(rbind, lapply(1:(N_iter), function(j) {
  #   data.frame(iteration = j,
  #              P = upper.tri.extractor(uploaded_results$chain1$st.deviations$tau_theta[,,j]),
  #              P_ij = paste0(upper_tri_indices[,1], upper_tri_indices[,2]))
  # }))
  # 
  # 
  # traceplot_proposal_P = ggplot(P_variance_df, aes(x = iteration, color = P_ij, group=P_ij))+
  #   geom_line(aes(y=P), alpha=.3)+
  #   facet_wrap(~P_ij)+
  #   theme_bw()+
  #   labs(title = "Adaptive variances proposals for P", x = 'Iterations')
  # 
  # plot_name_traceplot_proposal_P <- paste0(processed_wd,"//P_traceplot_proposal",true_model,est_model,"K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
  # png(plot_name_traceplot_proposal_P,width = 800, height = 800)
  # par(mar = c(1.5, 1.5,1.5,1.5))
  # print(traceplot_proposal_P)
  # dev.off()
  
  
  if(is.simulation==T){
    P_s_table_sum <- P_s_table_save%>%
      summarise(
        average_credible_length = mean(abs(quantile95 - quantile05)),
        MAE = mean(MAE)
      ) %>% round(3) 
  }else{
    P_s_table_sum <- P_s_table_save%>%
      summarise(
        average_credible_length = mean(abs(quantile95 - quantile05))
      ) %>% round(3) 
  }
  P_s_table_sum = P_s_table_sum %>% mutate(model = est_model)%>% mutate(n_clust = K)
  #adjusting colnames for the current number of clusters K
  if(!exists('Pcontainer')){
    Pcontainer = P_s_table_sum
  }else{
    Pcontainer =  rbind(Pcontainer,P_s_table_sum)
  }
  
  
  #-------------------------------------------------------------------------------
  # relabeling the chains to correct for label switching
  #-------------------------------------------------------------------------------
  # Relabel remaining chains
  z_burned_1 = uploaded_results$chain1$est_containers$z[,1:N_iter]
  z_burned_2 = uploaded_results$chain2$est_containers$z[,1:N_iter]
  z_burned_3 = uploaded_results$chain3$est_containers$z[,1:N_iter]
  z_burned_4 = uploaded_results$chain4$est_containers$z[,1:N_iter]
  
  theta_burned_1 = uploaded_results$chain1$est_containers$theta[,,1:N_iter]
  theta_burned_2 = uploaded_results$chain2$est_containers$theta[,,1:N_iter]
  theta_burned_3 = uploaded_results$chain3$est_containers$theta[,,1:N_iter]
  theta_burned_4 = uploaded_results$chain4$est_containers$theta[,,1:N_iter]
  
  z_burned_list = list(z_burned_1,
                       z_burned_2,
                       z_burned_3,
                       z_burned_4)
  
  theta_burned_list = list(theta_burned_1,
                           theta_burned_2,
                           theta_burned_3,
                           theta_burned_4)
  # chain_relabeled2 <- relabel_chain(2, permutations_z =  permutations_z, z_chain = z_burned_2, 
  #                                   ncol_iter = N_iter - burnin,n=n)
  # chain_relabeled3 <- relabel_chain(3,  permutations_z = permutations_z, z_chain = z_burned_3, 
  #                                   ncol_iter =  N_iter - burnin,n=n)
  # chain_relabeled4 <- relabel_chain(4, permutations_z = permutations_z, z_chain = z_burned_4,
  #                                   ncol_iter = N_iter - burnin,n=n)
  # 
  # 
  # # Permute P matrices for the remaining chains
  # theta_permuted2 <- permute_P(chain_index = 2,  permutations_z = permutations_z, P_chain = theta_burned_2,K=K)
  # theta_permuted3 <- permute_P(chain_index = 3, permutations_z = permutations_z, P_chain = theta_burned_3,K=K)
  # theta_permuted4 <- permute_P(4, permutations_z = permutations_z, P_chain = theta_burned_4,K=K)
  # 
  
  #-------------------------------------------------------------------------------
  # computing the estimated loglikelihood for each chain
  #-------------------------------------------------------------------------------
  
  
  num_samples = 1332
  filtering_obs = which(upper.tri(N_ij) & N_ij!= 0,arr.ind = T)
  upper.tri.Y_ij = Y_ij[filtering_obs]
  upper.tri.N_ij = N_ij[filtering_obs]
  
  Y_pred = matrix(NA, nrow = num_samples,ncol = length(upper.tri.Y_ij))
  
  LL_list <- foreach(i=1:4, .packages='foreach')%do%{
    
    z_chain <- z_burned_list[[i]]
    theta_chain <- theta_burned_list[[i]]
    
    llik = matrix(NA,  nrow = num_samples, ncol = length(upper.tri.Y_ij))
    
    for(t in 1:num_samples){
      
      z_chain_mat = vec2mat_0_P(z_chain[,t], theta_burned[,,1])
      
      P_entry = inverse_logit_f(theta_chain[,,t])
      
      P_ij = calculate_victory_probabilities(z_mat =z_chain_mat, P = P_entry)
      
      llik[t,] =  dbinom(upper.tri.Y_ij, upper.tri.N_ij, P_ij[filtering_obs],log = T)
      
      Y_pred[t,] <- rbinom(length(upper.tri.Y_ij), upper.tri.N_ij,
                           P_ij[filtering_obs])
    }
    print(i)
    return(llik)
  }
  
  
  LLik_sum <- lapply(LL_list,FUN = rowSums)
  
  # 
  # saveRDS(LLik_sum,file = paste0(processed_wd,"//loglik",true_model,est_model,K,".RDS"))
  
  #-------------------------------------------------------------------------------
  # printing traceplots of the likelihood
  #-------------------------------------------------------------------------------
  df_traceplot = data.frame(chain = c(rep(1,num_samples), rep(2,num_samples), rep(3,num_samples),rep(4,num_samples)),
                            log_likelihood = c(LLik_sum[[1]],LLik_sum[[2]],LLik_sum[[3]],LLik_sum[[4]]),
                            iterations = rep(1:num_samples,4))
  
  df_traceplot = df_traceplot %>% mutate(chain = factor(chain, 
                                                        levels = 1:4))
  
  
  
  # Plot the density lines
  stationary_density =  df_traceplot %>%
    ggplot(aes(y = log_likelihood, color = chain), alpha=0.5) + # Set 'x' to log_likelihood
    geom_density() +
    labs(title = "Stationary density for the 4 chains",
         subtitle = paste0("Number of iterations: ", 
                           N_iter," || Burnin: ", burnin),
         x = 'Density',
         y = "log-likelihood sum",
         color = "Chain",
         caption = paste0("True data: ", true_model, 
                          ", Fitted model: ", est_model, ", K = ", K))+ # Add a label for the legend
    theme_minimal() + # Use a minimal theme for better appearance
    theme(legend.position = "right") # Position the legend to the right
  
  
  stationary_name <- paste0(processed_wd,"//stationary_density",est_model, "_K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
  png(stationary_name,width = 500, height = 250)
  print(stationary_density)
  dev.off()
  # Convert 'chain' to a factor
  
  
  my_sexy_traceplot<- ggplot(df_traceplot, 
                             aes(x = iterations, 
                                 y = log_likelihood,
                                 color = factor(chain), 
                                 group=chain))+
    geom_line(alpha = .45)+
    labs(title = "Log likelihood for the 4 chains",
         subtitle = paste0("Number of iterations: ", 
                           N_iter," || Burnin: ", burnin),
         x = "Iterations",
         y = "Log likelihood",
         color = "Chain",
         caption = paste0("True data: ", true_model, 
                          ", Fitted model: ", est_model, ", K = ", K))+
    theme_bw()
  
  traceplot_name <- paste0(processed_wd,"//traceplot",est_model, "_K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
  png(traceplot_name,width = 500, height = 250)
  print(my_sexy_traceplot)
  dev.off()
  
  
  
  
  
  library(bayesplot)
  
  r_effs = loo::relative_eff(LL_list[[1]],rep(1,num_samples))
  
  loo_model_fit = loo(LL_list[[1]],cores = 3,save_psis = T,r_eff = r_effs)
  saveRDS(loo_model_fit, paste0(processed_wd,"/modelcheck",est_model,K,".RDS"))
  plot(loo_model_fit)
  waic_model_fit = waic(LL_list[[1]])
  
  
  
  problematic_values = pareto_k_ids(loo_model_fit)
  
  
  
  
  prob_df = data.frame(row = rownames(Y_ij)[row(Y_ij)[problematic_values]], 
                       col = colnames(Y_ij)[col(Y_ij)[problematic_values]]) 
  
  for(i in 1:length(prob_df)){
    prob_df$z_row[i] =  point_est_z[which(rownames(Y_ij)==prob_df$row[i])]
    prob_df$z_col[i] =  point_est_z[which(rownames(Y_ij)==prob_df$col[i])]
  }
  
  for(i in 1:length(prob_df)){
    prob_df$P_MCMC[i] =  P_est[prob_df$z_row[i],prob_df$z_col[i]]
    prob_df$P_hat[i] =  Y_ij[prob_df$row[i],prob_df$col[i]]/N_ij[prob_df$row[i],prob_df$col[i]]
    
  }
  
  
  pareto_k_values(loo_model_fit)
  
  # 
  # 
  # bad_values = loo::pareto_k_ids(loo_model_fit)
  # bad_values_df = data.frame(item = 0, value=0, bad=T,iteration=0)
  # 
  # good_values = setdiff(1:ncol(LL_list[[1]]), bad_values)
  # for(t in 1:25){
  #   bad_df = data.frame(item= bad_values, value = LL_list[[1]][t,bad_values], bad=T, iteration=t)
  #   good_df= data.frame(item= good_values, value = LL_list[[1]][t,good_values], bad=F, iteration=t) 
  #   bad_values_df = rbind(bad_values_df, bad_df, good_df)
  # }
  # bad_values_df%>%
  #   ggplot(aes(x=item, y=value, color=bad))+
  #   geom_point(alpha=0.4)
  # 
  # 
  # 
  # upper.tri.Y_ij[which(bad_values == Inf)] > upper.tri.N_ij[which(bad_values == Inf)]
  # upper.tri.Y_ij[is.na(bad_values)] > upper.tri.N_ij[is.na(bad_values)]
  # 
  # dbinom(upper.tri.Y_ij[which(bad_values == Inf)] , upper.tri.N_ij[which(bad_values == Inf)], P_ij[which(bad_values == Inf)], log = T)
  # dbinom(upper.tri.Y_ij[is.na(bad_values)] , upper.tri.N_ij[is.na(bad_values)], P_ij[is.na(bad_values)], log = T)
  # 
  # dbinom(upper.tri.Y_ij[961] , upper.tri.N_ij[961], P_ij[961], log = T)
  # 
  diagnostic1<- ppc_loo_pit_qq(
    y = upper.tri.Y_ij,
    yrep = Y_pred,
    lw = weights(loo_model_fit$psis_object)
  )+
    labs(x = 'uniform',title = paste0('Posterior Predictive LOO-PIT'),
         subtitle = paste0('Model =',est_model,", K=",K),
         caption = paste0("data source:",true_model))
  ggsave(plot=diagnostic1,filename = paste0(processed_wd,'diagnostic1',est_model,K,".png"))
  
  
  diagnostic2 = ppc_loo_pit_overlay( y = upper.tri.Y_ij,
                                     yrep = Y_pred,
                                     lw = weights(loo_model_fit$psis_object))+
    labs(x = 'uniform',title = paste0('Posterior Predictive LOO-PIT'),
         subtitle = paste0('Model =',est_model,", K=",K),
         caption = paste0("data source:",true_model))
  ggsave(plot=diagnostic2,filename = paste0(processed_wd,'diagnostic2',est_model,K,".png"))
  
  # 
  diagnostic3 = ppc_bars(y = upper.tri.Y_ij,
                         yrep = Y_pred)+
    labs(x = 'Y_ij values',title = paste0('Posterior Predictive Check'),
         subtitle = paste0('Model =',est_model,", K=",K),
         caption = paste0("data source:",true_model))
  ggsave(plot=diagnostic3,filename = paste0(processed_wd,'diagnostic3',est_model,K,".png"))
  
  
  if(is.simulation==F){
    z_s_table = data.frame(lone_out = loo_model_fit$estimates[1],lone_out_se = loo_model_fit$estimates[4],
                           waic = waic_model_fit$estimates[1], waic_se = waic_model_fit$estimates[4],
                           percent_bad_values=    round(pareto_k_table(loo_model_fit)[8],4)*100,
                           K0_hat =     K0_hat,
                           perc_label_switch = my_z_est$label_switch_count/N*100,
                           K0_hat_switched = K0_hat_switched)
  }else if(is.simulation==T){
    z_s_table = data.frame(lone_out = loo_model_fit$estimates[1],lone_out_se = loo_model_fit$estimates[4],
                           waic = waic_model_fit$estimates[1], waic_se = waic_model_fit$estimates[4],
                           percent_bad_values=    round(pareto_k_table(loo_model_fit)[8],4)*100,
                           vi_dist = vi.dist(uploaded_results$chain1$ground_truth$z, point_est_z),   
                           K0_hat_nolabel_sw=    K0_hat,
                           perc_label_switch = my_z_est$label_switch_count/N*100,
                           K0_hat_switched = K0_hat_switched) 
  }
  
  
  z_s_table = z_s_table %>% mutate(model=est_model)%>% mutate(n_clust = K)
  if(!exists('z_container')){
    z_container = z_s_table
  }else{
    z_container =  rbind(z_container,z_s_table)
  }
  
  if(!exists('long_df')){
    # Convert to long format
    long_df <- data.frame(loo_model_fit$estimates) %>%
      tibble::rownames_to_column(var = "Metric") %>%  # Convert row names to a column
      mutate(model = est_model) %>%
      mutate(num_clust = K)
    long_df1 <- data.frame(Metric = names(WAIC(LL_list[[1]])),Estimate = as.numeric(WAIC(LL_list[[1]])))%>%
      mutate(SE = NA)%>%
      mutate(model = est_model) %>%
      mutate(num_clust = K)
    long_df = rbind(long_df,long_df1)
  }else{
    long_1 = data.frame(loo_model_fit$estimates) %>%
      tibble::rownames_to_column(var = "Metric") %>%  # Convert row names to a column
      mutate(model = est_model) %>%
      mutate(num_clust = K)
    long_2 <- data.frame(Metric = names(WAIC(LL_list[[1]])),Estimate = as.numeric(WAIC(LL_list[[1]])))%>%
      mutate(SE = NA)%>%
      mutate(model = est_model) %>%
      mutate(num_clust = K)
    long_df = rbind(long_df,long_1,long_2)
  }
  
  
  # 
  # if(est_model ==  'WST'){
  #   
  #   #-------------------------------------------------------------------------------
  #   # sigma^2 parameter estimate
  #   #-------------------------------------------------------------------------------
  #   
  #   sigma_squared_s_table<- sigma_squared_summary_table(sigma_burned = sigma_squared_burned, sigma_true = uploaded_results$chain1$ground_truth$sigma_squared,
  #                                                       true_value = is.simulation*(true_model=='WST') , 
  #                                                       diag0.5 = TRUE, K = K, burnin = burnin)
  #   
  #   
  #   sigma_squared_s_table = sigma_squared_s_table %>% mutate(model=est_model)%>% mutate(n_clust = K)
  #   if(!exists('sigma_squared_container')){
  #     sigma_squared_container = sigma_squared_s_table}
  #   else{
  #     sigma_squared_container =  rbind(sigma_squared_container,sigma_squared_s_table)
  #   }
  # }
  if(est_model == 'SST'){
    #-------------------------------------------------------------------------------
    # mu parameter estimate
    #-------------------------------------------------------------------------------
    #
    
    mu_chain<- uploaded_results$chain1$est_containers$mu_vec
    mu_df = data.frame(mu_chain = 0, n_iter = 0, level_set =0)
    for(i in 1:(K)){
      mu_df = rbind(mu_df, 
                    data.frame(mu_chain = mu_chain[i,], 
                               n_iter = 1:ncol(mu_chain), 
                               level_set = i))
    }
    mu_df = mu_df[-1,]
    if(is.simulation ==F){
      mu_plot =  mu_df %>%ggplot(aes(n_iter, mu_chain, 
                                     group = level_set,color= as.factor(level_set)))+
        geom_line() 
      ggsave(mu_plot, filename = paste0(processed_wd,"mu_trace",true_model, est_model,K,".png"))
    }else if(is.simulation == T){
      
      mu_plot = mu_df$ground_truth = rep(uploaded_results$chain1$ground_truth$mu_vec_star,nrow(mu_df)/K)
      
      
      mu_plot = mu_df %>%ggplot(aes(n_iter, mu_chain, 
                                    group = level_set,color= as.factor(level_set)))+
        geom_line()+
        geom_hline(aes(yintercept=ground_truth))
      
      ggsave(mu_plot, filename = paste0(processed_wd,"mu_trace",true_model, est_model,K,".png"))
      
    }
    mu_vec_s_table<- mu_vec_summary_table(chains = uploaded_results, true_value = is.simulation*(true_model!='Simple'),
                                          diag0.5 = TRUE, K = K, burnin = burnin)
    
    mu_vec_s_table = mu_vec_s_table %>% mutate(model=rep(est_model,nrow(mu_vec_s_table))) %>% mutate(n_clust = rep(K,nrow(mu_vec_s_table)))
    if(is.simulation == T){
      mu_vec_s_table %>% 
        mutate(true_value = inverse_logit_f(uploaded_results$chain1$ground_truth$mu_vec_star))
    }
    if(!exists('mu_vec_container')){
      mu_vec_container = mu_vec_s_table
    }else{
      mu_vec_container =  rbind(mu_vec_container,mu_vec_s_table)
    }
    
  }
  
  
  
  
  #-------------------------------------------------------------------------------
  # P diagnostics 
  #-------------------------------------------------------------------------------
  # 
  # 
  # P_d_table<- P_diagnostic_table(chains = theta_list_relab, true_value = is.simulation, 
  #                                permutations_z = permutations_z,
  #                                diag0.5 = TRUE,K = K,
  #                                P = uploaded_results$chain1$ground_truth$theta,
  #                                burnin = burnin, N_iter = N_iter, label_switch =T)
  # 
  # 
  # P_d_table_save <- P_d_table$results 
  # 
  # P_d_table_save = P_d_table_save %>% mutate(model= est_model)%>% mutate(n_clust = K)
  # if(est_model=='SST'&file==1){
  #   P_d_container = P_d_table_save
  # }else{
  #   P_d_container =  rbind(P_d_container,P_d_table_save)
  # }
  # 
  

  # -------------------------------------------------------------------------------
  # z diagnostics
  # -------------------------------------------------------------------------------
  # 
  z_d_table <- z_diagnostic_table(chains = uploaded_results, true_value = is.simulation, diag0.5 = TRUE,
                                  K = K, burnin = N_iter*0.25, N_iter=N_iter,label_switch=F)
  
  z_d_table = z_d_table %>% mutate(model= est_model) %>% mutate(n_clust = K)
  
  
  if(!exists('z_d_container')){
    z_d_container = z_d_table
  }else{
    z_d_container =  rbind(z_d_container,z_d_table)
  }
  
  count_labels <- function(my_row, labels) {
    counts_list = table(factor(my_row,levels = labels))
    counts_list
  }
  
  # Assuming chains_SST$chain1$est_containers$z is a matrix or data frame
  mixing_labels <- apply(uploaded_results$chain1$est_containers$z, 1, count_labels, labels = 1:K)
  
  if(!is.null(rownames(Y_ij))){
    label_df = data.frame(item = rownames(Y_ij), t(mixing_labels)) %>%
      pivot_longer(-item)
  }else{
    label_df = data.frame(item = 1:n, t(mixing_labels)) %>%
      pivot_longer(-item)
  }
  order_df = label_df%>%
    filter(name=='X1')%>%
    arrange(value)%>%
    mutate(order = 1:n)%>%
    select(item, order)
  
  labels_plot = label_df %>% left_join(order_df , by =c("item")) %>%
    ggplot(aes(x = reorder(item,order), y = value, fill= name))+
    geom_col(alpha=0.7)+
    theme_bw()+
    labs(x = 'Items', y = 'Number of allocations for each label')+
    theme(axis.text.x = element_text(angle=90))
  
  # label_df <- label_df %>%
  #   group_by(item) %>%
  #   mutate(proportion = value / sum(value))
  # library(scatterpie)
  # 
  # label_df = data.frame(item = rownames(Y_ij), t(mixing_labels), win_prop = rowSums(Y_ij)/colSums(N_ij))%>%
  #   mutate(item = factor(item[order(win_prop,decreasing = T)]))
  # 
  # label_df =label_df%>%
  #   mutate(lab = order(label_df$win_prop,decreasing = T))%>%
  #   mutate_at(vars(starts_with("X")), ~ . / rowSums(label_df[,paste0('X',1:K)]))
  # 
  
  
  # label_df %>% 
  #   left_join(order_df , by =c("item")) %>%
  #   ggplot(aes(x = item, y = proportion, fill = name, group=name)) +
  #   geom_scatterpie(data = label_df,aes(x=item, y=lat, group=region, r=radius))+
  #   theme_void() +
  #   theme(legend.position = "right") +
  #   labs(fill = "Item")
  
  
  
  
  
  plot_name_auto <- paste0(processed_wd,"//labels_plot",true_model,est_model,"K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
  png(plot_name_auto,width = 800, height = 800)
  print(labels_plot)
  dev.off()
  
  
  
  
  
  # #
  # 
  
  # if(est_model == 'WST'){
  # 
  #   #-------------------------------------------------------------------------------
  #   # sigma^2 diagnostics
  #   #-------------------------------------------------------------------------------
  # 
  #   sigma_df = data.frame(iterations = 1:num_samples, 
  #                       sigma = uploaded_results$chain1$est_containers$sigma_squared[1:N_iter]) %>%
  #   ggplot(aes(x = iterations, y = sigma))+
  #   geom_line()
  #   
  #   ggsave(plot = sigma_df, filename = paste0(processed_wd,"sigma_trace",true_model, est_model,K,".png"))
  # 
  # 
  # 
  # 
  #   sigma_squared_d_table <- sigma_squared_diagnostic_table(chains = uploaded_results,
  #                                                           true_value = is.simulation*(true_model=='WST'), diag0.5 = TRUE, K = K,
  #                                                           burnin = burnin, N_iter = N_iter)
  # 
  #   sigma_squared_d_table = sigma_squared_d_table %>% mutate(model= est_model) %>% mutate(n_clust = K)
  #   if(!exists('sigma_squared_d_container')){
  #     sigma_squared_d_container = sigma_squared_d_table
  #   }else{
  #     sigma_squared_d_container =  rbind(sigma_squared_d_container,sigma_squared_d_table)
  #   }
  # #   
  # }
  # # if(est_model!="Simple"){
  # #   
  # #   #-------------------------------------------------------------------------
  # #   # mu diagnostics
  # #   #-------------------------------------------------------------------------
  # #   
  # #   mu_vec_d_table <- mu_vec_diagnostic_table(chains = uploaded_results, true_value = is.simulation*(true_model!='Simple'), 
  # #                                             diag0.5 = TRUE,
  # #                                             K = K, burnin = burnin,N_iter = N_iter)
  # #   
  # #   
  # #   mu_vec_d_table_save = mu_vec_d_table$results %>% mutate(model= est_model)%>% mutate(n_clust = K)
  # #   
  # #   if(est_model=='SST'&file==1){
  # #     mu_vec_d_container = mu_vec_d_table_save
  # #   }else{
  # #     mu_vec_d_container =  rbind(mu_vec_d_container,mu_vec_d_table_save)
  # #   }
  # # }
  # # 
  # 
  
  # #---------------------------------------------------------------------------
  # # Saving Plots and matrices
  # #---------------------------------------------------------------------------
  P_est_title <- paste0(processed_wd,'/P_est_matrix',true_model,est_model,K, '.csv')
  P_est <- round(P_est_relabeled,3) %>% data.frame()
  P_est %>% write.csv(file = P_est_title)

  if(is.simulation==T){
    P_true_title <- paste0(processed_wd,'/P_true_matrix',true_model,K, '.csv')
    P_true <- round(uploaded_results$chain1$ground_truth$theta,3) %>% inverse_logit_f() %>% data.frame()
    P_true %>% write.csv(file = P_true_title)
  }
  
  
  # 
  # # convergence diagnostics plot -----------------------------------------------
  # 
  # P_list <- P_d_table$plots_list
  # mu_list <- mu_vec_d_table$plots_list
  # 
  #---------------------------------------------------------------------------
  #DIAGNOSTICS PLOTS FOR P 
  #---------------------------------------------------------------------------
  
  # plot_name_gelman <- paste0(processed_wd,"//P_gelman_rubin_plot%03d",true_model,est_model,"K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
  # # Save the plot with the constructed file name
  # png(plot_name_gelman,width = 800, height = 800)
  # par(mar = c(1.5, 1.5,1.5,1.5))
  # gelman.plot(P_list)
  # dev.off()
  # 
  
  #Crosscorrelation
  # plot_name_cross <- paste0(processed_wd,"//P_crosscorr_plot",true_model,est_model,"K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
  # png(plot_name_cross,width = 800, height = 800)
  # par(mfrow = c(1,1))
  # crosscorr.plot(P_list)
  # dev.off()
  # #
  # #Autocorrelation
  # auto_plot <- acfplot(P_list, type = 'l')
  # plot_name_auto <- paste0(processed_wd,"//P_autocorr_plot%03d",true_model,est_model,"K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
  # png(plot_name_auto,width = 800, height = 800)
  # print(auto_plot)
  # dev.off()
  
  
  
  #---------------------------------------------------------------------------
  #DIAGNOSTICS PLOTS FOR mu
  #---------------------------------------------------------------------------
  
  # if(est_model != 'Simple'){
  #   tryCatch({
  #     #Gelman Rubin
  #     plot_name_gelman_U <- paste0(processed_wd,"//U_gelman_rubin_plot%03d",true_model,est_model,"K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
  #     # Save the plot with the constructed file name
  #     png(plot_name_gelman_U,width = 800, height = 800)
  #     par(mar = c(1.5, 1.5,1.5,1.5))
  #     gelman.plot(mu_list)
  #     dev.off()
  #     
  #     
  #     
  #     #Crosscorrelation
  #     plot_name_cross_U <- paste0(processed_wd,"//U_crosscorr_plot",true_model,est_model,"K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
  #     png(plot_name_cross_U,width = 800, height = 800)
  #     par(mfrow = c(1,1))
  #     crosscorr.plot(mu_list)
  #     dev.off()
  #     
  #     #Autocorrelation
  #     auto_plot_U = acfplot(mu_list, type = 'l')
  #     plot_name_auto_U <- paste0(processed_wd,"//U_autocorr_plot%03d",true_model,est_model,"K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
  #     png(plot_name_auto_U,width = 800, height = 800)
  #     print(auto_plot_U)
  #     dev.off()
  #   }, error = function(e){
  #     message("An error occurred:\n", e)}
  #   )
  #   
  # }
  # 
  #---------------------------------------------------------------------------
  # COMPARISON WITH THE EXTERNAL RANKING
  #---------------------------------------------------------------------------
  
  if(is.simulation==F){
    
    
    plot_name <- paste0(processed_wd,"//RankvsClust_Est_model",est_model, "_K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
    # Save the plot with the constructed file name
    
    est_df<- data.frame(Id = rownames(Y_ij), 
                        marginal_victories = rowSums(Y_ij),
                        marginal_losses = colSums(Y_ij),
                        est_cl = point_est_z, 
                        unique_identifier = runif(N, -0.001,0.001))%>%
      mutate(unique_identifier = unique_identifier+est_cl)%>%
      mutate(relative_victories = marginal_victories/(marginal_losses+marginal_victories))
    
    
    if(est_model =='Simple'){
      new_order = est_df %>%
        group_by(est_cl)%>%
        summarise(mean_win = mean(relative_victories))%>%
        arrange(-mean_win)
      
      # Example new order (a permutation of 1:K)
      new_z = vector()
      for(item in 1:n){
        new_z = append(new_z,which(new_order$est_cl == point_est_z[item]))
      }
      
      est_df$est_cl = new_z
      
      
      write.csv(new_z, paste0(processed_wd,"/z_est_K_",K,"model",est_model,".csv"))
    }
    
    combined_df = players_df %>% inner_join(est_df, by = 'Id')%>% dplyr::arrange(unique_identifier)
    
    if(true_model=='Tennis_data'){
      degree_plot <- ggplot(combined_df, aes(x = reorder(Id, -relative_victories), y = relative_victories, fill =factor(est_cl) )) +
        geom_bar(stat = "identity") +
        labs(x = "Player Name", y = "Percentage Victories", fill='Cluster', title = "Percentage of victories for each player",
             subtitle = 'Players are sorted according to their victory %') +
        scale_color_discrete(name = "Cluster") +
        theme_bw() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
          text = element_text(size = 10, family = "Arial"),
          plot.title = element_text(size = 13, face = "bold", margin = margin(r = 10)),
          plot.subtitle = element_text(size = 10, margin = margin(t = 10, r = 10, b = 10)),  # Adjust the top margin
          plot.margin = margin(20, 20, 20, 20)
        )
    }else if(true_model == 'Citation_data'){
      degree_plot <- ggplot(combined_df, aes(x = reorder(Id, median_rank), y = relative_victories, fill =factor(est_cl) )) +
        geom_bar(stat = "identity") +
        labs(x = "Journal Name", y = "Citations received / Citations made", fill='Cluster', 
             title = "Citations received / citations made for each journal",
             subtitle = 'Players are ranked by the Stigler model of Varin et al (2016)') +
        scale_color_discrete(name = "Cluster") +
        theme_bw() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
          text = element_text(size = 10, family = "Arial"),
          plot.title = element_text(size = 13, face = "bold", margin = margin(r = 10)),
          plot.subtitle = element_text(size = 10, margin = margin(t = 10, r = 10, b = 10)),  # Adjust the top margin
          plot.margin = margin(20, 20, 20, 20)
        )
      
      
    }
    png(plot_name,width = 800, height = 627)
    print(degree_plot)
    dev.off()
    # 
    #---------------------------------------------------------------------------
    # CHECKING THE HOMOGENEITY OF THE CLUSTERS: HEATMAP
    #---------------------------------------------------------------------------
    
    ggplot(combined_df, aes(x=relative_victories, y=last_rank)) +
      # Add points
      geom_point(color = "blue", size = 2) +
      # Add a red dotted line at 45 degrees passing through the origin
      geom_smooth(method = "lm", color = "green", se = FALSE, linetype = "solid", size = 1) +
      # Add title and labels
      ggtitle("Scatter Plot of Relative Victories vs Median Rank") +
      xlab("Relative Victories %") +
      ylab("Median ATP Rank") +
      scale_x_continuous(labels = scales::percent) +
      # Improve the overall appearance
      theme_minimal() +
      scale_y_reverse()+
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        panel.grid.major = element_line(color = "grey", linetype = "dotted"),
        panel.grid.minor = element_blank()
      )
    
    
    
    plot_name <- paste0(processed_wd,"//RankvsClust_Est_model",est_model, "_K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
    
    colnames(Y_ij) <- rownames(Y_ij)
    
    
    for(i in seq_len(nrow(z_df_complete))){
      z_df_complete$marginal_victories_row[i] <- sum(Y_ij[z_df_complete$row[i],]/est_df$unique_identifier,na.rm = T)/sum(Y_ij[,z_df_complete$row[i]]*est_df$unique_identifier,na.rm = T)
    }
    for(i in seq_len(nrow(z_df_complete))){
      z_df_complete$marginal_victories_col[i] <- sum(Y_ij[z_df_complete$col[i],]/est_df$unique_identifier,na.rm = T)/sum(Y_ij[,z_df_complete$col[i]]*est_df$unique_identifier,na.rm = T)
    }
    
    
    
    
    plot_df = z_df_complete%>%
      inner_join(est_df, by = c("row" = "Id")) %>%
      dplyr::rename(row_z = est_cl) %>%
      inner_join(est_df, by = c("col" = "Id")) %>%
      dplyr::rename(col_z = est_cl) %>%
      mutate(row = factor(row, levels = unique(row[order(row_z, -marginal_victories_row, decreasing = TRUE)])),
             col = factor(col, levels = unique(col[order(col_z, -marginal_victories_col, decreasing = F)]))) %>%
      mutate(perc_success = Y/N)
    
    
    
    v_lines_list = list()
    for(k in 1:(K-1)){
      lines_df = plot_df %>% filter(col_z == k, row_z == k)
      v_lines_list[[k]] <- lines_df$row[which.min(lines_df$marginal_victories_row)]
    }
    
    
    est_df_disordered = est_df %>% mutate(est_cl =order(Id))
    plot_df_disordered = z_df_complete%>%
      inner_join(est_df_disordered, by = c("col" = "Id")) %>%
      inner_join(est_df_disordered, by = c("row" = "Id")) %>%
      mutate(perc_success = Y/N) %>%
      mutate(col = factor(col, levels = unique(col[order(est_cl.x, col)])),
             row = factor(row, levels = unique(row[order(est_cl.y , row, decreasing = TRUE)])))
    
    
    
    rowSums(Y_ij)
    if(true_model == 'Tennis_data'){
      adjacency_m   <- ggplot(plot_df, aes(x = col, y = row)) +
        geom_tile(aes(fill = perc_success), color = "gray", show.legend = T) +
        scale_fill_gradient(low = "white", high = "red") +
        geom_ysidecol(aes(x = degree_pl_row, color=factor(row_z))) +
        geom_vline(xintercept = unlist(v_lines_list), color='black')+
        geom_hline(yintercept = unlist(v_lines_list), color='black')+
        labs(x = paste0("Players ordered by blocks"),
             y = paste0("Playersordered by blocks"),
             fill = "% victories",
             color = "Block")+
        theme(legend.position = 'bottom', legend.direction = 'horizontal')+
        theme_bw() +
        theme(axis.text.y = element_blank(),
              axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank())
      
      adjacency_m_disordered<- ggplot(plot_df_disordered, aes(x = col, y = row)) +
        geom_tile(aes(fill = perc_success), color = "gray", show.legend = F) +
        scale_fill_gradient(low = "white", high = "red") +
        labs(x = paste0("Players ordered by blocks"),
             y = paste0("Players ordered by blocks"),
             fill = "% relative citations",
             color = "Block")+
        theme(legend.position = 'bottom', legend.direction = 'horizontal')+
        theme_bw() +
        theme(axis.text.y = element_blank(),
              axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.title.y = element_blank())
      
      
      
      
      # Assuming adjacency_m and adjacency_m_disordered are your two ggplot objects
      
      # Combine the plots
      # Extract the legend from the first plot
      legend <- get_legend(adjacency_m)
      
      # Combine the plot without the legend and the extracted legend
      combined_plot_with_legend <- plot_grid(adjacency_m + theme(legend.position = "none"), 
                                             legend, 
                                             rel_widths = c(1, 0.2))
      
      # Now combine with the second plot
      final_combined_plot <- plot_grid(adjacency_m_disordered, combined_plot_with_legend, 
                                       ncol = 2, labels = 'auto',
                                       rel_widths = c(0.75, 1))
      #saving the heatmap
      plot_name1<- paste0(processed_wd,"//OrderedVsUnordered",est_model, "_K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
      png(plot_name1,width = 900, height = 450)
      print(final_combined_plot)
      dev.off()
      
      
    }else if(true_model == 'Citation_data'){
      adjacency_m   <- ggplot(plot_df, aes(x = col, y = row)) +
        geom_tile(aes(fill = perc_success), color = "gray", show.legend = T) +
        scale_fill_gradient(low = "white", high = "red") +
        geom_ysidecol(aes(x = degree_pl_row, color=factor(row_z))) +
        geom_vline(xintercept = unlist(v_lines_list), color='black')+
        geom_hline(yintercept = unlist(v_lines_list), color='black')+
        labs(x = paste0("Journals ordered by blocks"),
             y = paste0("Journals ordered by blocks"),
             fill = "% cit. received/made",
             color = "Block")+
        theme(legend.position = 'bottom', legend.direction = 'horizontal')+
        theme_bw()+
        theme(axis.text.y = element_blank(),
              axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank())
      if(est_model == 'SST'&file==4){
        adjacency_m_disordered<- ggplot(plot_df_disordered, aes(x = col, y = row)) +
          geom_tile(aes(fill = perc_success), color = "gray", show.legend = F) +
          scale_fill_gradient(low = "white", high = "red") +
          labs(x = paste0("Journals ordered by blocks"),
               y = paste0("Journals ordered by blocks"),
               fill = "% citations received/made",
               color = "Block")+
          theme(legend.position = 'bottom', legend.direction = 'horizontal')+
          theme_bw() +
          theme(axis.text.y = element_blank(),
                axis.text.x = element_blank(),
                axis.title.x = element_blank(),
                axis.title.y = element_blank())
        
        
        
        # Assuming adjacency_m and adjacency_m_disordered are your two ggplot objects
        
        # Combine the plots
        # Extract the legend from the first plot
        legend <- get_legend(adjacency_m)
        
        # Now combine with the second plot
        final_combined_plot <- plot_grid(adjacency_m_disordered, adjacency_m + guides(color = FALSE, fill = FALSE),
                                         legend, 
                                         ncol = 3,rel_widths = c(1,1.2,0.5))
        #saving the heatmap
        plot_name1<- paste0(processed_wd,"//cowplot_plot",est_model, "_K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
        png(plot_name1,width = 1100, height = 428)
        print(final_combined_plot)
        dev.off()
      }
    }
    
    
    #saving the heatmap
    plot_name1<- paste0(processed_wd,"//Combined_plot",est_model, "_K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
    png(plot_name1,width = 800, height = 594)
    print(adjacency_m)
    dev.off()
    
    
    
    #---------------------------------------------------------------------------
    # CHECKING THE RANKING and THE CLUSTERING
    #---------------------------------------------------------------------------
    
    percentage_to_display <- 10
    set.seed(23)
    # Randomly sample a subset of labels to display
    
    
    
    
    if(true_model == 'Citation_data'){
      rank_boxplot<- ggplot(combined_df, aes(x = factor(est_cl), y = median_rank,color = factor(est_cl))) +
        geom_boxplot(aes(fill=factor(est_cl)),alpha=.3) +
        geom_label_repel(
          data = combined_df,  # Use the sampled labels for display
          aes(label = Id),
          size = 3,
          hjust = .5,
          vjust = 0,
          show.legend = F,
          alpha=.8
        ) +
        labs(title= "Rank of the players divided into blocks",
             subtitle = "Not all names are displayed to avoid overlapping",
             x = "Clusters",
             y = "Varin et al rank",
             color = "Cluster",
             fill = "Cluster")+
        theme_classic()+
        theme(legend.position = "bottom",
              plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
              plot.subtitle = element_text(face = "bold.italic", hjust = 0.5),
              plot.caption = element_text(face = "italic"))
    }else if(true_model == 'Tennis_data'){
      
      for(k in 1:K){
        
        k_members = combined_df %>% filter(est_cl == k) 
        proportion_wanted = max(round(percentage_to_display / 100 * nrow(k_members)),1)
        members_to_sample = proportion_wanted
        sampled_row = sample(1:nrow(k_members), size = members_to_sample)
        if(!exists('sampled_labels')){
          sampled_labels <- k_members[sampled_row, ]
        }else{
          sampled_labels = rbind(sampled_labels,k_members[sampled_row, ])
        }
      }
      
      rank_boxplot<- ggplot(combined_df, aes(x = factor(est_cl), y = relative_victories, color = factor(est_cl))) +
        geom_boxplot(alpha=.3) +
        geom_label_repel(
          data = sampled_labels,  # Use the sampled labels for display
          aes(label = Id),
          size = 3,
          hjust = .5,
          vjust = 0,
          show.legend = F,
          alpha=.8
        ) +
        labs(title= "Rank of the players divided into blocks",
             subtitle = "Not all names are displayed to avoid overlapping",
             x = "Clusters",
             y = "Win proportions",
             color = "Cluster",
             fill = "Cluster")+
        theme_classic()+
        theme(legend.position = "bottom",
              plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
              plot.subtitle = element_text(face = "bold.italic", hjust = 0.5),
              plot.caption = element_text(face = "italic"))
      rm(sampled_labels)
      
      
      rank_label<- ggplot(combined_df, aes(x = factor(est_cl), y = relative_victories, color = factor(est_cl))) +
        geom_point()+
        geom_label_repel(aes(label = Id),
                         size = 3,
                         hjust = .5,
                         vjust = 0,
                         show.legend = F,
                         alpha=.8
        ) +
        labs(title= "Rank of the players divided into blocks",
             subtitle = "Not all names are displayed to avoid overlapping",
             x = "Clusters",
             y = "Win proportions",
             color = "Cluster",
             fill = "Cluster")+
        theme_classic()+
        theme(legend.position = "bottom",
              plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
              plot.subtitle = element_text(face = "bold.italic", hjust = 0.5),
              plot.caption = element_text(face = "italic"))
      rm(sampled_labels)
    }
    
    plot_name2<- paste0(processed_wd,"//Boxplot",est_model, "_K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
    png(plot_name2,width = 800, height = 594)
    print(rank_boxplot)
    dev.off()
    
    plot_name3<- paste0(processed_wd,"//Rank_Labels",est_model, "_K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
    png(plot_name3,width = 1200, height = 1000)
    print(rank_label)
    dev.off()
  }
  beepr::beep('coin')
}


if(is.simulation == F){
  top_block_df_container = top_block_df_container %>% arrange(n_clust) %>% 
    relocate (n_clust) %>% relocate(where(is.character)) %>% rename(K = n_clust)
  top_block_df_container %>%  write.csv(file = paste0(processed_wd,"/top_block_df.csv"))
}

z_container = z_container %>% arrange(n_clust) %>% 
  relocate (n_clust) %>% relocate(where(is.character)) %>% rename(K = n_clust)

z_container %>%  write.csv(file = paste0(processed_wd,"/z_container.csv"))


Pcontainer = Pcontainer %>% arrange(n_clust ) %>% relocate (n_clust) %>% 
  relocate(where(is.character)) %>% rename(K = n_clust) 

Pcontainer %>% write.csv(file = paste0(processed_wd,"/Pcontainer.csv"))

# sigma_squared_container = sigma_squared_container %>% arrange(n_clust ) %>% relocate (n_clust) %>% 
#   relocate(where(is.character)) %>% rename(K = n_clust) 
# 
# sigma_squared_container%>% write.csv(file = paste0(processed_wd,"/sigma_squared_container.csv"))
# 

write.csv(long_df,paste0(processed_wd,"/long_df_model_choice.csv"))

check_P_posterior_sum = check_P_posterior%>%
  group_by(model,K)%>%
  summarise(obs = n(),
            is.SST_sum = sum(is.SST),
            is.WST_sum = sum(is.WST))

write.csv(check_P_posterior_sum,paste0(processed_wd,"/check_P_transitivity.csv"))

mu_vec_container <- mu_vec_container %>% arrange(n_clust) %>% arrange(n_clust ) %>% relocate (n_clust) %>% 
  relocate(where(is.character)) %>% rename(K = n_clust) 

mu_vec_container%>% write.csv(file = paste0(processed_wd,"/mu_vec_container.csv"))


model_selection_plot = long_df %>%
  filter(Metric == 'looic')%>%
  mutate(K = as.factor(num_clust))%>%
  ggplot(aes(x = K, y = Estimate, color = model))+
  geom_point(shape = 18, size=3)+
  labs(title = paste0('Leave-one-out for different K and models'), y = 'Looic',
       subtitle = 'Lower values are better')


plot_name2<- paste0(processed_wd,"//Model_choice",true_model,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
png(plot_name2,width = 500, height = 400)
print(model_selection_plot)
dev.off()


uploaded_results$chain1$est_containers$theta[1,1,1] ==uploaded_results$chain1$est_containers$theta[1,1,5] 
# 
# 
# z_d_container = z_d_container %>% arrange(n_clust ) %>% relocate (n_clust) %>%
#   relocate(where(is.character))
# 
# z_d_container %>% write.csv(file = paste0(processed_wd,"/z_d_container.csv"))
# 
# P_d_container = P_d_container %>% arrange(n_clust ) %>% relocate (n_clust) %>%
#   relocate(where(is.character)) %>% rename(K = n_clust)
# 
# P_d_container%>% write.csv(file = paste0(processed_wd,"/P_d_container.csv"))
# 
# sigma_squared_d_container = sigma_squared_d_container %>% arrange(n_clust ) %>%
#   relocate (n_clust) %>%
#   relocate(where(is.character)) %>% rename(K = n_clust)
# 
# sigma_squared_d_container %>% write.csv(file = paste0(processed_wd,"/sigma_squared_d_container.csv"))
# 
# 
# mu_vec_d_container = mu_vec_d_container %>% arrange(n_clust ) %>% relocate (n_clust) %>%
#   relocate(where(is.character)) %>% rename(K = n_clust)
# 
# mu_vec_d_container %>%  write.csv(file = paste0(processed_wd,"/mu_vec_d_container.csv"))

if(is.simulation==F){
  #BLOCKWISE HEATMAP ---------------------------------------------------------
  
  
  chosen_model <- z_container[which(z_container$lone_out ==  max(z_container$lone_out)),]
  
  
  
  uploaded_results<- readRDS(paste0(data_wd, 'Data_from', true_model,'_est_model',chosen_model$model[1],'_Kest',chosen_model$K[1],'.RDS'))
  K<- chosen_model$K[1]
  theta_est <- apply(uploaded_results$chain1$est_containers$theta[,,1:N_iter], MARGIN = c(1,2), mean)
  P_est <- inverse_logit_f(theta_est)
  z_burned <- uploaded_results$chain1$est_containers$z[,1:N_iter]
  
  my_z_est<- z_plot(z_burned = z_burned ,
                    A =uploaded_results$chain1$control_containers$A[,1:N_iter], 
                    Y_ij = Y_ij, N_ij = N_ij, 
                    true_model= true_model,P_est = P_est,est_model = chosen_model$model[1]
                    , true_value =is.simulation, 
                    diag0.5 =diag0.5 , K=K, N=nrow(uploaded_results$chain1$Y_ij), z = uploaded_results$chain1$ground_truth$z ,
                    burnin =  burnin ,label_switch = T,tap= processed_wd)
  
  
  point_est<- as.vector(my_z_est$point_est)
  table(point_est)
  
  
  est_df<- data.frame(player_slug = rownames(Y_ij), est_cl = point_est)
  
  
  
  
  my_beauti_plottini = list()
  K<- length(unique( point_est))
  
  for(block_i in 1:K){
    for(block_j in 1:K){
      
      players12 <- est_df %>% mutate(player_in_block_i = est_cl %in% block_i)%>%
        mutate(player_in_block_j = est_cl %in% block_j)
      
      A_subset<- Y_ij[players12$player_in_block_i, players12$player_in_block_j]
      N_subset<- N_ij[players12$player_in_block_i, players12$player_in_block_j]
      
      
      new_df<- as.data.frame(as.table(as.matrix(A_subset)))%>%
        rename(player1 = Var1)%>%
        rename(player2 = Var2)%>%
        rename(n_victories = Freq)
      
      new_df1<- as.data.frame(as.table(as.matrix(N_subset)))%>%
        rename(player1 = Var1)%>%
        rename(player2 = Var2)%>%
        rename(n_games = Freq)
      
      subset_df = inner_join(new_df, new_df1, by = c("player1","player2")) %>%
        mutate(perc_victories= n_victories/n_games)
      #counting the number of victories and games played by each player in player2 column vs all the others
      total_victories_and_games_df_player2 <- subset_df %>%
        group_by(player2) %>%
        summarize(total_victories_player2 = sum(n_victories), total_games_player2 = sum(n_games)) %>%
        mutate(total_victories_player2 = total_games_player2- total_victories_player2)
      
      #counting the number of victories and games played by each player in player1 column vs all the others
      total_victories_and_games_df_player1 <- subset_df %>%
        group_by(player1) %>%
        summarize(total_victories_player1 = sum(n_victories), total_games_player1 = sum(n_games))
      
      # Merge the total_victories and total_games back to the original data frame
      final_df = subset_df %>% inner_join(total_victories_and_games_df_player1, by=c("player1")) %>%
        inner_join(total_victories_and_games_df_player2, by=c("player2")) %>%
        mutate(marginal_player1 = total_victories_player1/total_games_player1)%>%
        mutate(marginal_player2 = total_victories_player2/total_games_player2) %>%
        mutate(color_marginal_1 = if_else(as.numeric(marginal_player1) >.5,"light","dark")) %>%
        mutate(color_marginal_2 = if_else(as.numeric(marginal_player2) > .5,"light","dark"))
      
      
      
      combined_plot12<- ggplot(final_df)+
        geom_tile(aes(x=player2, y= player1, fill= perc_victories), color="grey9")+
        geom_ysidetile(aes(y= player1,x = "marginal_y", fill = `marginal_player1`)) +
        geom_ysidetext(aes(y= player1,x = "marginal_y",
                           label = round(`marginal_player1`,1))) +
        geom_xsidetile(aes(x= player2,y = "marginal_x", fill = `marginal_player2`)) +
        geom_xsidetext(aes(x= player2,y = "marginal_x",
                           label = round(`marginal_player2`,1)),
                       angle = 90) +
        scale_fill_gradient(low = "white", high = "red") +
        theme_bw()+
        theme(legend.direction = "vertical", 
              axis.text.x = element_text(angle=45, hjust = 1),
              plot.title =  element_text(face = "bold", hjust = 0.5))+
        labs( x = paste0("Players in block ", block_j),
              y = paste0("Players in block ", block_i),
              fill = "% victories",
              color = "Block")+
        guides(fill = 'none', color='none')
      
      my_beauti_plottini[[paste0("blocks",block_i, block_j)]] <- combined_plot12
    }
  }
  
  for(i in 1:K){
    j=i-1
    png(paste0(processed_wd,"/decomposed_heatmap",i,".png"),width = 2000, height = 400)
    print(cowplot::plot_grid(plotlist = my_beauti_plottini[c(1:K)+(K*j)], nrow=1))
    dev.off()
  }
  
}










  