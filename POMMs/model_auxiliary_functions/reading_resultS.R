
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

source("/Users/lapo_santi/Desktop/Nial/oldmaterial/POMM_flex/functions_container_flex.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/model_auxiliary_functions/Functions_priorSST.R")
source("/Users/lapo_santi/Desktop/Nial/oldmaterial/project/simplified model/SaraWade.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/model_auxiliary_functions/Inference_orderstats.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/model_auxiliary_functions/MCMC_functions.R")



#FLAG is.simulation=T IF YOU ARE READING THE RESULTS FOR A SIMULATION STUDY
is.simulation = F

if(is.simulation==F){
  
  true_model = "Tennis_data"
  if(true_model == "Tennis_data"){
    #where the data are stored
    setwd('/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/')
    Y_ij <- read.table("./Data/Tennis application/Y_ij.csv",header  = F,row.names = 1,sep = ",")
    N_ij <- read.table("./Data/Tennis application/N_ij.csv",header  = F,row.names = 1,sep = ",")
    
    Y_ij = as.matrix(Y_ij)
    N_ij = as.matrix(N_ij)
    
    data_wd<- "/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/results/application/Tennis_data/"
    #where the data are saved
    processed_wd <- "/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/results/application/Tennis_data/processed/"

    df_rank <- readRDS("/Users/lapo_santi/Desktop/Nial/weekly material/Tennis application/data/df_rank.RData")
    

    
    ranks= df_rank   %>%
      filter(week_year==2017)
    
    
    df_match <- readRDS("/Users/lapo_santi/Desktop/Nial/weekly material/Tennis application/data/df_match.RData")
    ranks= df_rank   %>%
      filter(week_year==2017)  %>% group_by(player_slug) %>% summarise(median_rank = median(rank_number),
                                                                       max_r = max(rank_number),
                                                                       min_r = min(rank_number))
    
    top100players = ranks %>% filter(median_rank <= 100) %>% arrange(median_rank)

    
    players_df = data.frame(Id = rownames(Y_ij)) %>% inner_join(top100players,by= c('Id'='player_slug'))
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
    bt_rank_df = data.frame(Id = rownames(abilities)[order(abilities[,1])], BTrank = 95:1)
    players_df = players_df %>% inner_join(bt_rank_df, by = 'Id') %>%
      mutate(degree_pl = rowSums(Y_ij)/colSums(Y_ij))

    
    
    comparison = plot(players_df$median_rank,players_df$BTrank)
    
    
    
    
    
    
    
    
    
  }else if( true_model == 'Citation_data'){
    data_wd<- "/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/results/application/Citation_data/"
    #where the data are saved
    processed_wd <- "/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/results/application/Citation_data/processed/"
    
    
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
    
    
    
    comparison = plot(players_df$median_rank,players_df$BTrank)
    
    
  }
  
  
  
}else if(is.simulation == T){
  true_model = "Simple"
  data_wd = "/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/results/simulation/simulation/Simple_true//"
  processed_wd <- "/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/results/simulation/simulation/Simple_true/processed/"
  
  
  
}


for(est_model in c('SST','WST','Simple')){
  est_model = 'SST'
  # filenames <- list.files(pattern = paste0('True_Model',true_model,'Est_model_', est_model),path = data_wd)
  filenames <- list.files(pattern = paste0("Est_model_",est_model),path = data_wd)
  print(filenames)
  
  for(file in 1:length(filenames)){

    uploaded_results<- readRDS(paste0(data_wd,"/",filenames[file]))
    uploaded_results <- chains_SST
    print(paste0('Now estimating ', filenames[file]))
    print(paste0(length(filenames)-file+1,' within the same class left '))
    N= nrow(uploaded_results$chain1$Y_ij)
    n=N
    N_iter = dim(uploaded_results$chain1$est_containers$z)[[2]]
    K = dim(uploaded_results$chain1$est_containers$theta)[[1]]
    burnin = max(N_iter - 40000,1)
    Y_ij <- uploaded_results$chain1$Y_ij
    N_ij <- uploaded_results$chain1$N_ij
    
    
    
    #-------------------------------------------------------------------------------
    # P temporary estimate
    #-------------------------------------------------------------------------------
    P_burned = uploaded_results$chain1$est_containers$theta[,,-c(1:burnin)]
    z_burned =  uploaded_results$chain1$est_containers$z[,-c(1:burnin)]
    if(est_model != 'Simple'){
      m_vec_burned = uploaded_results$chain1$est_containers$mu_vec[,-c(1:burnin)]
    }
    if(est_model == 'WST'){
      sigma_squared_burned = uploaded_results$chain1$est_containers$sigma_squared[-c(1:burnin)]
    }
    
    
    
    
    theta = apply(P_burned, c(1,2), mean)
    P_est = inverse_logit_f(theta)
    
    my_z_est<- z_plot(z_burned = z_burned,  A = uploaded_results$chain1$control_containers$A[-c(1:burnin)],Y_ij = Y_ij, N_ij = N_ij, true_model= true_model,P_est = P_est,
                      est_model = est_model, true_value =is.simulation, 
                      diag0.5 =diag0.5 , K=K, N=nrow(uploaded_results$chain1$Y_ij), z_true = uploaded_results$chain1$ground_truth$z ,
                      burnin =  burnin ,label_switch = T,tap= processed_wd)
    
    
    point_est_z<- as.vector(my_z_est$point_est)
    
    table(point_est_z)
    
    
    K_est<- length(unique(point_est_z))
    permutations_z<- my_z_est$permutations
    z_chain_permuted<- my_z_est$relabeled_chain
    
    
    
    
    #---------------------------------------------------------------------------
    # P parameter estimate
    #-------------------------------------------------------------------------------
    
    P_s_table <- P_summary_table(P_burned = P_burned,
                                 true_value = is.simulation,
                                 permutations_z = permutations_z,
                                 diag0.5 = TRUE,
                                 K = K, P_true = uploaded_results$chain1$ground_truth$theta,
                                 burnin = burnin,
                                 label_switch = T)
    
    P_s_table_save <-P_s_table$table
    theta_chain_permuted <- P_s_table$P_permuted
    P_est_relabeled<- P_s_table$P_hat
    
    inverse_logit_f(P_est_relabeled)
    
    #some traceplots
    upper_tri_indices= which(upper.tri(P_est_relabeled, diag=T),arr.ind = T)
    if(is.simulation ==T){
      P_trace_df_post_switch <- do.call(rbind, lapply(1:(N_iter-burnin), function(j) {
        data.frame(iteration = j,
                   theta = upper.tri.extractor(theta_chain_permuted[,,j]), 
                   theta_true = upper.tri.extractor(uploaded_results$chain1$ground_truth$theta), 
                   entry = paste0(upper_tri_indices[,1], upper_tri_indices[,2]))
      }))
      P_trace_df_post_switch=P_trace_df_post_switch%>% mutate(P = inverse_logit_f(theta))%>%
        mutate(P_true = inverse_logit_f(theta_true))
      
      traceplot_P = ggplot(P_trace_df_post_switch, aes(x = iteration, color = entry, group=entry))+
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
    
    
    
    plot_name_traceplot_P <- paste0(processed_wd,"//P_traceplot",true_model,est_model,"K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
    png(plot_name_traceplot_P,width = 800, height = 800)
    par(mar = c(1.5, 1.5,1.5,1.5))
    print(traceplot_P)
    dev.off()
    
    P_variance_df <- do.call(rbind, lapply(1:(N_iter), function(j) {
      data.frame(iteration = j,
                 P = upper.tri.extractor(uploaded_results$chain1$st.deviations$tau_theta[,,j]),
                 P_ij = paste0(upper_tri_indices[,1], upper_tri_indices[,2]))
    }))


    traceplot_proposal_P = ggplot(P_variance_df, aes(x = iteration, color = P_ij, group=P_ij))+
      geom_line(aes(y=P), alpha=.3)+
      facet_wrap(~P_ij)+
      theme_bw()+
      labs(title = "Adaptive variances proposals for P", x = 'Iterations')

    plot_name_traceplot_proposal_P <- paste0(processed_wd,"//P_traceplot_proposal",true_model,est_model,"K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
    png(plot_name_traceplot_proposal_P,width = 800, height = 800)
    par(mar = c(1.5, 1.5,1.5,1.5))
    print(traceplot_proposal_P)
    dev.off()
    
    
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
    if(est_model== 'SST' & file ==1){
      Pcontainer = P_s_table_sum
    }else{
      Pcontainer =  rbind(Pcontainer,P_s_table_sum)
    }
    
    
    
    #-------------------------------------------------------------------------------
    # relabeling the chains to correct for label switching
    #-------------------------------------------------------------------------------
    # Relabel remaining chains
    
    z_burned_2 = uploaded_results$chain2$est_containers$z[,-c(1:burnin)]
    z_burned_3 = uploaded_results$chain3$est_containers$z[,-c(1:burnin)]
    z_burned_4 = uploaded_results$chain3$est_containers$z[,-c(1:burnin)]
    
    theta_burned_2 = uploaded_results$chain2$est_containers$theta[,,-c(1:burnin)]
    theta_burned_3 = uploaded_results$chain3$est_containers$theta[,,-c(1:burnin)]
    theta_burned_4 = uploaded_results$chain3$est_containers$theta[,,-c(1:burnin)]
    
    
    
    chain_relabeled2 <- relabel_chain(2, permutations_z =  permutations_z, z_chain = z_burned_2, 
                                      ncol_iter = N_iter - burnin,n=n)
    chain_relabeled3 <- relabel_chain(3,  permutations_z = permutations_z, z_chain = z_burned_3, 
                                      ncol_iter =  N_iter - burnin,n=n)
    chain_relabeled4 <- relabel_chain(4, permutations_z = permutations_z, z_chain = z_burned_4,
                                      ncol_iter = N_iter - burnin,n=n)
    
    
    # Permute P matrices for the remaining chains
    theta_permuted2 <- permute_P(chain_index = 2,  permutations_z = permutations_z, P_chain = theta_burned_2,K=K)
    theta_permuted3 <- permute_P(chain_index = 3, permutations_z = permutations_z, P_chain = theta_burned_3,K=K)
    theta_permuted4 <- permute_P(4, permutations_z = permutations_z, P_chain = theta_burned_4,K=K)
    
    z_list_relab = list(z1 = z_chain_permuted,z2=chain_relabeled2,z3=chain_relabeled3,z4=chain_relabeled4)
    theta_list_relab = list(P1 = theta_chain_permuted,P2=theta_permuted2,P3=theta_permuted3,P4=theta_permuted4)
    
    #-------------------------------------------------------------------------------
    # computing the estimated loglikelihood for each chain
    #-------------------------------------------------------------------------------


    num_samples = N_iter - burnin
    
    LL_list <- foreach(i=1:4, .packages='foreach') %do% {
      z_chain <- z_list_relab[[i]]
      theta_chain <- theta_list_relab[[i]]
      
      u_y = Y_ij[upper.tri(N_ij)&(N_ij!=0)]
      u_n = N_ij[upper.tri(N_ij)&(N_ij!=0)]

      
      
      P_chain = array(apply(X = theta_chain, MARGIN = 3, FUN = inverse_logit_f),dim=c(K,K,N_iter-burnin))
      z_mat_array = array(apply(X = z_chain, MARGIN = 2, FUN = vec2mat_0_P,P=P_chain[,,1]), dim=c(n, K, N_iter-burnin))
      
      llik = matrix(NA, length(u_n), N_iter-burnin)
      
      for(t in 1:(N_iter- burnin)){
        P_ij= calculate_victory_probabilities(z_mat_array[,,t], P_chain[,,t])
        llik[,t] = dbinom(x = u_y, size = u_n, P_ij[upper.tri(N_ij)&(N_ij!=0)], log=T)
      }
      
      print(i)
      return(llik)
    }

    # Stop the cluster after the loop
        
    
    LLik_sum <- lapply(LL_list,FUN = colSums)
    
    
    
    
    saveRDS(LLik_sum,file = paste0(processed_wd,"//loglik",true_model,est_model,K,".RDS"))
    LL = LL_list[[2]]
    #-------------------------------------------------------------------------------
    # printing traceplots of the likelihood
    #-------------------------------------------------------------------------------
    df_traceplot = data.frame(chain = c(rep(1,ncol(LL)), rep(2,ncol(LL)), rep(3,ncol(LL)),rep(4,ncol(LL))),
                              log_likelihood = c(LLik_sum[[1]],LLik_sum[[2]],LLik_sum[[3]],LLik_sum[[4]]),
                              iterations = rep(((burnin+1):N_iter),4))
    df_traceplot = df_traceplot %>% mutate(chain = factor(chain, levels = 1:4))
    
    
    my_sexy_traceplot<- ggplot(df_traceplot, aes(x = iterations, y = log_likelihood, 
                                                 color = factor(chain), group=chain))+
      geom_line(alpha = .5)+
      labs(title = "Log likelihood for the 4 chains",
           subtitle = paste0("Number of iterations: ", N_iter," || Burnin: ", burnin), 
           x = "Iterations",
           y = "Log likelihood",
           color = "Chain",
           caption = paste0("True data: ", true_model, ", Fitted model: ", est_model, ", K = ", K))+
      theme_bw()
    traceplot_name <- paste0(processed_wd,"//traceplot",est_model, "_K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
    png(traceplot_name,width = 500, height = 250)
    print(my_sexy_traceplot)
    dev.off()
    
    loo_model_fit = loo(t(LL))
    
    
    plot(loo_model_fit)
    waic_model_fit = waic(t(LL))
    
    
    if(is.simulation==F){
    z_s_table = data.frame(lone_out = loo_model_fit$estimates[3],lone_out_se = loo_model_fit$estimates[6],
                           waic = waic_model_fit$estimates[3], waic_se = waic_model_fit$estimates[6],
                           percent_bad_values=    round(pareto_k_table(loo_model_fit)[8],4)*100 )
    }else if(is.simulation==T){
      z_s_table = data.frame(lone_out = loo_model_fit$estimates[3],lone_out_se = loo_model_fit$estimates[6],
                             waic = waic_model_fit$estimates[3], waic_se = waic_model_fit$estimates[6],
                             percent_bad_values=    round(pareto_k_table(loo_model_fit)[8],4)*100,
                             vi_dist = vi.dist(uploaded_results$chain1$ground_truth$z, point_est_z)) 
    }
    
    
    z_s_table = z_s_table %>% mutate(model=est_model)%>% mutate(n_clust = K)
    if(est_model=='SST'&file==1){
      z_container = z_s_table
    }else{
      z_container =  rbind(z_container,z_s_table)
    }
    
    
    
    
    if(est_model ==  'WST'){
      
      #-------------------------------------------------------------------------------
      # sigma^2 parameter estimate
      #-------------------------------------------------------------------------------
      
      sigma_squared_s_table<- sigma_squared_summary_table(sigma_burned = sigma_squared_burned, sigma_true = uploaded_results$chain1$ground_truth$sigma_squared,
                                                          true_value = is.simulation*(true_model=='WST') , 
                                                          diag0.5 = TRUE, K = K, burnin = burnin)
      
      
      sigma_squared_s_table = sigma_squared_s_table %>% mutate(model=est_model)%>% mutate(n_clust = K)
      if(est_model=='WST' & file==1){
        sigma_squared_container = sigma_squared_s_table}
      else{
        sigma_squared_container =  rbind(sigma_squared_container,sigma_squared_s_table)
      }
    }
    if(est_model != 'Simple'){
      #-------------------------------------------------------------------------------
      # mu parameter estimate
      #-------------------------------------------------------------------------------
      #
      
      mu_chain<- uploaded_results$chain1$est_containers$mu_vec
      mu_df = data.frame(mu_chain = 0, n_iter = 0, level_set =0)
      for(i in 1:(K+1)){
      mu_df = rbind(mu_df, 
                    data.frame(mu_chain = mu_chain[i,], n_iter = 1:ncol(mu_chain), 
                               level_set = i))
      }
      
      
      mu_plot = mu_df %>%ggplot(aes(n_iter, mu_chain, group = level_set,color= as.factor(level_set)))+
        geom_line() 
      ggsave(filename = paste0(processed_wd,"mu_trace",true_model, est_model,K,".png"))
      
      
      
      
      mu_vec_s_table<- mu_vec_summary_table(chains = uploaded_results, true_value = is.simulation*(true_model!='Simple'),
                                            diag0.5 = TRUE, K = K, burnin = burnin)
      
      mu_vec_s_table = mu_vec_s_table %>% mutate(model=rep(est_model,nrow(mu_vec_s_table))) %>% mutate(n_clust = rep(K,nrow(mu_vec_s_table)))
      if(is.simulation == T){
        mu_vec_s_table %>% mutate(true_value = inverse_logit_f(uploaded_results$chain1$ground_truth$mu_vec_star))
      }
      if(est_model=='SST'&file==1){
        mu_vec_container = mu_vec_s_table
      }else{
        mu_vec_container =  rbind(mu_vec_container,mu_vec_s_table)
      }
      
    }
    
    
    
    
    #-------------------------------------------------------------------------------
    # P diagnostics 
    #-------------------------------------------------------------------------------
    # 
    
    P_d_table<- P_diagnostic_table(chains = theta_list_relab, true_value = is.simulation, 
                                   permutations_z = permutations_z,
                                   diag0.5 = TRUE,K = K,
                                   P = uploaded_results$chain1$ground_truth$theta,
                                   burnin = burnin, N_iter = N_iter, label_switch =T)
    
    
    P_d_table_save <- P_d_table$results 
    
    P_d_table_save = P_d_table_save %>% mutate(model= est_model)%>% mutate(n_clust = K)
    if(est_model=='SST'&file==1){
      P_d_container = P_d_table_save
    }else{
      P_d_container =  rbind(P_d_container,P_d_table_save)
    }
    
    
    
    # -------------------------------------------------------------------------------
    # z diagnostics
    # -------------------------------------------------------------------------------
    
    z_d_table <- z_diagnostic_table(chains = uploaded_results, true_value = is.simulation, diag0.5 = TRUE,
                                    K = K, burnin = N_iter*0.25, N_iter=N_iter,label_switch=F)
    
    z_d_table = z_d_table %>% mutate(model= est_model) %>% mutate(n_clust = K)
    
    if(est_model=='SST'&file==1){
      z_d_container = z_d_table
    }else{
      z_d_container =  rbind(z_d_container,z_d_table)
    }
    
    #
    
    
    if(est_model == 'WST'){
      
      #-------------------------------------------------------------------------------
      # sigma^2 diagnostics 
      #-------------------------------------------------------------------------------
      
      sigma_squared_d_table <- sigma_squared_diagnostic_table(chains = uploaded_results, 
                                                              true_value = is.simulation*(true_model=='WST'), diag0.5 = TRUE, K = K, 
                                                              burnin = burnin, N_iter = N_iter)
      
      sigma_squared_d_table = sigma_squared_d_table %>% mutate(model= est_model) %>% mutate(n_clust = K)
      if(est_model=='WST'&file==1){
        sigma_squared_d_container = sigma_squared_d_table
      }else{
        sigma_squared_d_container =  rbind(sigma_squared_d_container,sigma_squared_d_table)
      }
      
    }
    if(est_model!="Simple"){

      #-------------------------------------------------------------------------
      # mu diagnostics
      #-------------------------------------------------------------------------

      mu_vec_d_table <- mu_vec_diagnostic_table(chains = uploaded_results, true_value = is.simulation*(true_model!='Simple'), 
                                                diag0.5 = TRUE,
                                                K = K, burnin = burnin,N_iter = N_iter)


      mu_vec_d_table_save = mu_vec_d_table$results %>% mutate(model= est_model)%>% mutate(n_clust = K)

      if(est_model=='SST'&file==1){
        mu_vec_d_container = mu_vec_d_table_save
      }else{
        mu_vec_d_container =  rbind(mu_vec_d_container,mu_vec_d_table_save)
      }
    }

    
    
    #---------------------------------------------------------------------------
    # Saving Plots and matrices
    #---------------------------------------------------------------------------
    P_est_title <- paste0(processed_wd,'/P_est_matrix',true_model,est_model,K, '.csv')
    P_est <- round(P_est_relabeled,3) %>% data.frame() 
    P_est %>% write.csv(file = P_est_title)
    
    if(is.simulation==T){
      P_true_title <- paste0(processed_wd,'/P_true_matrix',true_model,K, '.csv')
      P_true <- round(uploaded_results$chain1$ground_truth$theta,3) %>% inverse_logit_f() %>% data.frame() 
      P_true %>% write.csv(file = P_true_title)
    }
    
    # convergence diagnostics plot -----------------------------------------------
    
    P_list <- P_d_table$plots_list
    mu_list <- mu_vec_d_table$plots_list
    # 
    #---------------------------------------------------------------------------
    #DIAGNOSTICS PLOTS FOR P 
    #---------------------------------------------------------------------------
    
    plot_name_gelman <- paste0(processed_wd,"//P_gelman_rubin_plot%03d",true_model,est_model,"K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
    # Save the plot with the constructed file name
    png(plot_name_gelman,width = 800, height = 800)
    par(mar = c(1.5, 1.5,1.5,1.5))
    gelman.plot(P_list)
    dev.off()
    
    
    #Crosscorrelation
    plot_name_cross <- paste0(processed_wd,"//P_crosscorr_plot",true_model,est_model,"K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
    png(plot_name_cross,width = 800, height = 800)
    par(mfrow = c(1,1))
    crosscorr.plot(P_list)
    dev.off()
    #
    #Autocorrelation
    auto_plot <- acfplot(P_list, type = 'l')
    plot_name_auto <- paste0(processed_wd,"//P_autocorr_plot%03d",true_model,est_model,"K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
    png(plot_name_auto,width = 800, height = 800)
    print(auto_plot)
    dev.off()
    
    
    
    #---------------------------------------------------------------------------
    #DIAGNOSTICS PLOTS FOR mu
    #---------------------------------------------------------------------------
    
    if(est_model != 'Simple'){
      tryCatch({
      #Gelman Rubin
      plot_name_gelman_U <- paste0(processed_wd,"//U_gelman_rubin_plot%03d",true_model,est_model,"K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
      # Save the plot with the constructed file name
      png(plot_name_gelman_U,width = 800, height = 800)
      par(mar = c(1.5, 1.5,1.5,1.5))
      gelman.plot(mu_list)
      dev.off()
      
      
      
      #Crosscorrelation
      plot_name_cross_U <- paste0(processed_wd,"//U_crosscorr_plot",true_model,est_model,"K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
      png(plot_name_cross_U,width = 800, height = 800)
      par(mfrow = c(1,1))
      crosscorr.plot(mu_list)
      dev.off()
      
      #Autocorrelation
      auto_plot_U = acfplot(mu_list, type = 'l')
      plot_name_auto_U <- paste0(processed_wd,"//U_autocorr_plot%03d",true_model,est_model,"K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
      png(plot_name_auto_U,width = 800, height = 800)
      print(auto_plot_U)
      dev.off()
      }, error = function(e){
        message("An error occurred:\n", e)}
    )
    
    }
    
    #---------------------------------------------------------------------------
    # COMPARISON WITH THE EXTERNAL RANKING
    #---------------------------------------------------------------------------
    if(is.simulation==F){
      
      
      plot_name <- paste0(processed_wd,"//RankvsClust_Est_model",est_model, "_K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
      # Save the plot with the constructed file name
      
  

      est_df<- data.frame(Id = rownames(Y_ij), est_cl = point_est_z) 

      combined_df = players_df %>% inner_join(est_df, by = 'Id')%>% dplyr::arrange(median_rank)
      
      if(true_model=='Tennis_data'){
        degree_plot <- ggplot(combined_df, aes(x = reorder(Id, BTrank), y = degree_pl, fill =factor(est_cl) )) +
          geom_bar(stat = "identity") +
          labs(x = "Player Name", y = "Percentage Victories", fill='Cluster', title = "Percentage of victories for each player",
               subtitle = 'Players are sorted according to their Bradley Terry model') +
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
        degree_plot <- ggplot(combined_df, aes(x = reorder(Id, median_rank), y = degree_pl, fill =factor(est_cl) )) +
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
      
      
      plot_name <- paste0(processed_wd,"//RankvsClust_Est_model",est_model, "_K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
      
      
      
      
      plot_df = z_df_complete%>%
        inner_join(est_df, by = c("col" = "Id")) %>%
        inner_join(est_df, by = c("row" = "Id")) %>%
        mutate(col = factor(col, levels = unique(col[order(est_cl.x, col)])),
               row = factor(row, levels = unique(row[order(est_cl.y , row, decreasing = TRUE)])))%>%
        mutate(perc_success = Y/N)
      
      # est_df_disordered = est_df %>% mutate(est_cl = sample(1:5, 47,T))
      # plot_df_disordered = z_df_complete%>%
      #   inner_join(est_df_disordered, by = c("col" = "Id")) %>%
      #   inner_join(est_df_disordered, by = c("row" = "Id")) %>%
      #   mutate(perc_success = Y/N) %>%
      #   mutate(col = factor(col, levels = unique(col[order(est_cl.x, col)])),
      #          row = factor(row, levels = unique(row[order(est_cl.y , row, decreasing = TRUE)])))
      #   
      
      
      
      if(true_model == 'Tennis_data'){
      adjacency_m   <- ggplot(plot_df, aes(x = col, y = row)) +
          geom_tile(aes(fill = perc_success), color = "gray", show.legend = T) +
          scale_fill_gradient(low = "white", high = "red") +
          geom_ysidecol(aes(x = degree_pl_row, color=factor(est_cl.y))) +
          labs(title = 'Heatmap filled with victory percentages',
               x = paste0("Players ordered by blocks"),
               y = paste0("Playersordered by blocks"),
               fill = "% victories",
               color = "Block")+
          theme(legend.position = 'bottom', legend.direction = 'horizontal')+
          theme_minimal()+
          theme(axis.text.y = element_blank(),
                axis.text.x = element_blank(),
                axis.title.x = element_blank(),
                axis.title.y = element_blank())
      }else if(true_model == 'Citation_data'){
        adjacency_m<- ggplot(plot_df, aes(x = col, y = row)) +
          geom_tile(aes(fill = perc_success), color = "gray", show.legend = T) +
          scale_fill_gradient(low = "white", high = "red") +
          geom_ysidecol(aes(x = degree_pl_row, color=factor(est_cl.y))) +
          labs(title = 'Heatmap filled with percetage citation received over citations made',
               x = paste0("Players ordered by blocks"),
               y = paste0("Playersordered by blocks"),
               fill = "% relative citations",
               color = "Block")+
          theme(legend.position = 'bottom', legend.direction = 'horizontal')+
          theme_minimal() +
          theme(axis.text.y = element_text(angle=30),
                axis.text.x = element_blank(),
                axis.title.x = element_blank(),
                axis.title.y = element_blank())
        
        adjacency_m_disordered<- ggplot(plot_df_disordered, aes(x = col, y = row)) +
          geom_tile(aes(fill = perc_success), color = "gray", show.legend = T) +
          scale_fill_gradient(low = "white", high = "red") +
          geom_ysidecol(aes(x = degree_pl_row, color=factor(est_cl.y))) +
          labs(title = 'Heatmap filled with percetage citation received over citations made',
               x = paste0("Players ordered by blocks"),
               y = paste0("Playersordered by blocks"),
               fill = "% relative citations",
               color = "Block")+
          theme(legend.position = 'bottom', legend.direction = 'horizontal')+
          theme_minimal() +
          theme(axis.text.y = element_text(angle=30),
                axis.text.x = element_blank(),
                axis.title.x = element_blank(),
                axis.title.y = element_blank())
      }
      
      
      #saving the heatmap
      plot_name1<- paste0(processed_wd,"//Combined_plot",est_model, "_K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
      png(plot_name1,width = 800, height = 594)
      print(adjacency_m)
      dev.off()
      
      
      
      #---------------------------------------------------------------------------
      # CHECKING THE RANKING and THE CLUSTERING
      #---------------------------------------------------------------------------
      
      percentage_to_display <- 30
      set.seed(23)
      # Randomly sample a subset of labels to display
      sampled_labels <- combined_df[sample(nrow(combined_df), size = round(percentage_to_display / 100 * nrow(combined_df))), ]
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
        rank_boxplot<- ggplot(combined_df, aes(x = factor(est_cl), y = BTrank,color = factor(est_cl))) +
          geom_boxplot(aes(fill=factor(est_cl)),alpha=.3) +
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
               y = "Bradley Terry rank",
               color = "Cluster",
               fill = "Cluster")+
          theme_classic()+
          theme(legend.position = "bottom",
                plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
                plot.subtitle = element_text(face = "bold.italic", hjust = 0.5),
                plot.caption = element_text(face = "italic"))
      }
      plot_name2<- paste0(processed_wd,"//Boxplot",est_model, "_K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
      png(plot_name2,width = 800, height = 594)
      print(rank_boxplot)
      dev.off()
    }
    beepr::beep('coin')
  }
}


z_container = z_container %>% arrange(n_clust) %>% 
  relocate (n_clust) %>% relocate(where(is.character)) %>% rename(K = n_clust)

z_container %>%  write.csv(file = paste0(processed_wd,"/z_container.csv"))

Pcontainer = Pcontainer %>% arrange(n_clust ) %>% relocate (n_clust) %>% 
  relocate(where(is.character)) %>% rename(K = n_clust) 

Pcontainer %>% write.csv(file = paste0(processed_wd,"/Pcontainer.csv"))

sigma_squared_container = sigma_squared_container %>% arrange(n_clust ) %>% relocate (n_clust) %>% 
  relocate(where(is.character)) %>% rename(K = n_clust) 

sigma_squared_container%>% write.csv(file = paste0(processed_wd,"/sigma_squared_container.csv"))



mu_vec_container <- mu_vec_container %>% arrange(n_clust) %>% arrange(n_clust ) %>% relocate (n_clust) %>% 
  relocate(where(is.character)) %>% rename(K = n_clust) 

mu_vec_container%>% write.csv(file = paste0(processed_wd,"/mu_vec_container.csv"))



z_d_container = z_d_container %>% arrange(n_clust ) %>% relocate (n_clust) %>%
  relocate(where(is.character))

z_d_container %>% write.csv(file = paste0(processed_wd,"/z_d_container.csv"))

P_d_container = P_d_container %>% arrange(n_clust ) %>% relocate (n_clust) %>%
  relocate(where(is.character)) %>% rename(K = n_clust)

P_d_container%>% write.csv(file = paste0(processed_wd,"/P_d_container.csv"))

sigma_squared_d_container = sigma_squared_d_container %>% arrange(n_clust ) %>%
  relocate (n_clust) %>%
  relocate(where(is.character)) %>% rename(K = n_clust)

sigma_squared_d_container %>% write.csv(file = paste0(processed_wd,"/sigma_squared_d_container.csv"))


mu_vec_d_container = mu_vec_d_container %>% arrange(n_clust ) %>% relocate (n_clust) %>%
  relocate(where(is.character)) %>% rename(K = n_clust)

mu_vec_d_container %>%  write.csv(file = paste0(processed_wd,"/mu_vec_d_container.csv"))

if(is.simulation==F){
  #BLOCKWISE HEATMAP ---------------------------------------------------------
  
  chosen_model <- z_container[which(z_container$waic ==  min(z_container$waic)),]
  
  
  
  uploaded_results<- readRDS(paste0(data_wd, 'True_Model', true_model,'Est_model_',chosen_model$model,'_N',n,'_K',chosen_model$K,'.RDS'))
  K<- chosen_model$K
  theta_est <- apply(uploaded_results$chain1$est_containers$theta[,,-c(1:burnin)], MARGIN = c(1,2), mean)
  P_est <- inverse_logit_f(theta_est)
  z_burned <- uploaded_results$chain1$est_containers$z[,-c(1:burnin)]
  
  my_z_est<- z_plot(z_burned = z_burned ,A =uploaded_results$chain1$control_containers$A[,-c(1:burnin)], Y_ij = Y_ij, N_ij = N_ij, 
                    true_model= true_model,P_est = P_est,est_model = chosen_model$model
                    , true_value =is.simulation, 
                    diag0.5 =diag0.5 , K=chosen_model$K, N=nrow(uploaded_results$chain1$Y_ij), z = uploaded_results$chain1$ground_truth$z ,
                    burnin =  burnin ,label_switch = T,tap= processed_wd)
  
  
  point_est<- as.vector(my_z_est$point_est)
  table(point_est)
  
  g_df =  data.frame(vertex_attr(g)) %>%
    rename(player_slug= name) %>%
    left_join(players_df, by="player_slug") %>%
    mutate(degree_pl = degree(g,mode = 'out')/degree(g,mode = 'all')) %>%
    arrange()
  est_df<- data.frame(player_slug = rownames(Y_ij), est_cl = point_est)
  combined_df<- inner_join(g_df,est_df,by = 'player_slug')
  
  
  
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
    png(paste0(processed_wd,"/decomposed_heatmap",i,".png"),width = 1800, height = 400)
    print(cowplot::plot_grid(plotlist = my_beauti_plottini[c(1:K)+(K*j)], nrow=1))
    dev.off()
  }
  
}
