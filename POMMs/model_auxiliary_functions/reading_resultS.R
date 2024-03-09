

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


#where the data are stored
data_wd<- "/Users/lapo_santi/Desktop/Nial/MCMC_results/application_13Feb2024/raw/"
#where the data are saved
processed_wd <- "/Users/lapo_santi/Desktop/Nial/MCMC_results/application_13Feb2024/tennis/estimates/"

#FLAG is.simulation=T IF YOU ARE READING THE RESULTS FOR A SIMULATION STUDY
is.simulation = F

if(is.simulation==F){
  df_rank <- readRDS("/Users/lapo_santi/Desktop/Nial/weekly material/Tennis application/data/df_rank.RData")
  df_match <- readRDS("/Users/lapo_santi/Desktop/Nial/weekly material/Tennis application/data/df_match.RData")
  ranks= df_rank   %>%
    filter(week_year==2017)  %>% group_by(player_slug) %>% summarise(median_rank = median(rank_number),max_r = max(rank_number),min_r = min(rank_number))
  
  top100players = ranks %>% filter(median_rank <= 100) %>% arrange(median_rank)
  
  #adding one extra column with the player id
  df_r =  inner_join(ranks,df_rank%>% select(player_slug,player_id), by='player_slug')
  
  #now, for each game I want to filter just those players in the top one-hundred
  df_match = df_match %>% filter(winner_slug %in% top100players$player_slug) %>% filter(loser_slug %in% top100players$player_slug)
  
  
  my_edges = df_match %>% select(winner_slug, loser_slug)
  g =graph_from_edgelist(as.matrix(my_edges),directed = T)
  my_name = data.frame(player_slug=vertex_attr(g)$name)
  players_df = inner_join(my_name,top100players, by="player_slug")
  
  A = as_adjacency_matrix(g)
  true_model = "Tennis_data"
}else if(is.simulation == T){
  true_model = "SST"
}


for(est_model in c('SST','WST','Simple')){
  
  filenames <- list.files(pattern = paste0('True_Model',true_model,'Est_model_', est_model),path = data_wd)
  print(filenames)
  
  for(file in 1:length(filenames)){
    
    uploaded_results<- readRDS(paste0(data_wd,"/",filenames[file]))
    
    
    print(paste0('Now estimating ', filenames[file]))
    print(paste0(length(filenames)-file+1,' within the same class left '))
    N= nrow(uploaded_results$chain1$Y_ij)
    n=N
    N_iter = dim(uploaded_results$chain1$est_containers$z)[[2]]
    K = dim(uploaded_results$chain1$est_containers$P)[[1]]
    burnin = N_iter-20000
    Y_ij <- uploaded_results$chain1$Y_ij
    N_ij <- uploaded_results$chain1$N_ij
    
    #-------------------------------------------------------------------------------
    # P temporary estimate
    #-------------------------------------------------------------------------------
    
    P_est <- apply(uploaded_results$chain1$est_containers$P[,,-c(1:burnin)], MARGIN = c(1,2), mean)
    P_est <- inverse_logit_f(P_est)
    
    P_true_upper <- upper.tri.extractor(uploaded_results$chain1$est_containers$P[,,1])
    upper_tri_indices <- which(upper.tri(P_est, diag = T), arr.ind = TRUE)
    P_chain = uploaded_results$chain1$est_containers$P[,,-c(1:burnin)]
    
    
    
    #---------------------------------------------------------------------------
    # SUMMARY PLOT PLOTS FOR z
    #---------------------------------------------------------------------------
    
    # setwd(plots_dir)
    my_z_est<- z_plot(chains = uploaded_results , true_model= true_model,P_est = P_est,
                      est_model = est_model, true_value =is.simulation, 
                      diag0.5 =diag0.5 , K=K, N=nrow(uploaded_results$chain1$Y_ij), z = uploaded_results$chain1$ground_truth$z ,
                      burnin =  burnin ,label_switch = T,tap= processed_wd)
    
    
    point_est_z<- as.vector(my_z_est$point_est)
    
    K_est<- length(unique(point_est_z))
    permutations_z<- my_z_est$permutations
    z_chain_permuted<- my_z_est$relabeled_chain
    
    
    #---------------------------------------------------------------------------
    # P parameter estimate
    #-------------------------------------------------------------------------------
    
    P_s_table <- P_summary_table(chains = uploaded_results,
                                 true_value = is.simulation,
                                 permutations_z = permutations_z,
                                 diag0.5 = TRUE,
                                 K = K, P = uploaded_results$chain1$ground_truth$P,
                                 burnin = burnin,
                                 label_switch = T)
    
    P_s_table_save <-P_s_table$table
    P_chain_permuted <- P_s_table$P_permuted
    P_est_relabeled<- P_s_table$P_hat
    
    
    
    P_trace_df_post_switch <- do.call(rbind, lapply(1:10000, function(j) {
      data.frame(iteration = j,
                 P = upper.tri.extractor(P_chain_permuted[,,j]), 
                 P_true = P_true_upper, 
                 P_ij = paste0(upper_tri_indices[,1], upper_tri_indices[,2]))
    }))
    
    
    traceplot_P = ggplot(P_trace_df_post_switch, aes(x = iteration, color = P_ij, group=P_ij))+
      geom_line(aes(y=P), alpha=.3)+
      geom_line(aes(y=P_true), linetype=2, color='red')+
      facet_wrap(~P_ij)+
      theme_bw()
    
    
    
    plot_name_traceplot_P <- paste0(processed_wd,"//P_traceplot",true_model,est_model,"K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
    png(plot_name_traceplot_P,width = 800, height = 800)
    par(mar = c(1.5, 1.5,1.5,1.5))
    print(traceplot_P)
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
    chain_relabeled2 <- relabel_chain(2, permutations_z, uploaded_results, N_iter - burnin,n=n)
    chain_relabeled3 <- relabel_chain(3, permutations_z, uploaded_results, N_iter - burnin,n=n)
    chain_relabeled4 <- relabel_chain(4, permutations_z, uploaded_results, N_iter - burnin,n=n)
    
    
    # Permute P matrices for the remaining chains
    P_permuted2 <- permute_P(2, permutations_z, uploaded_results,K=K)
    P_permuted3 <- permute_P(3, permutations_z, uploaded_results,K=K)
    P_permuted4 <- permute_P(4, permutations_z, uploaded_results,K=K)
    
    z_list_relab = list(z1 = z_chain_permuted,z2=chain_relabeled2,z3=chain_relabeled3,z4=chain_relabeled4)
    P_list_relab = list(P1 = P_chain_permuted,P2=P_permuted2,P3=P_permuted3,P4=P_permuted4)
    
    #-------------------------------------------------------------------------------
    # computing the estimated loglikelihood for each chain
    #-------------------------------------------------------------------------------
    
    # Set up parallel backend
    cl <- makeCluster(4)  # Adjust the number of cores accordingly
    registerDoParallel(cl)
    
    # Export necessary variables to the workers
    clusterExport(cl, list("P_list_relab", "compute_likelihood_foreach","z_list_relab", "Y_ij", "N_ij", "inverse_logit_f", "vec2mat_0_P", "calculate_victory_probabilities", "dbinom", "P_chain"), envir = .GlobalEnv)
    
    num_samples = N_iter - burnin
    # Perform parallel computation using foreach
    LL_list <- foreach(i = 1:4 ) %do% {
      z_chain <- z_list_relab[[i]]
      P_chain <- P_list_relab[[i]]
      
      # Define the number of chunks in which to split the likelihood
      num_chunks <- 5
      # Split the columns into chunks for parallel processing
      chunk_size <- ceiling(num_samples / num_chunks)
      chunks <- split(1:num_samples, cut(1:num_samples, breaks = num_chunks, labels = FALSE))
      
      # Apply function to each chunk
      LL <- foreach(chunk_idx = 1:num_chunks, .combine = "cbind") %dopar% {
        chunk <- chunks[[chunk_idx]]
        sapply(chunk, compute_likelihood_foreach,)
      }
      LL
    }
    # Stop the cluster after the loop
    stopCluster(cl)
    
    
    LLik_sum <- lapply(LL_list,FUN = colSums)
    saveRDS(LLik_sum,file = paste0(processed_wd,"//loglik",true_model,est_model,K))
    #-------------------------------------------------------------------------------
    # printing traceplots of the likelihood
    #-------------------------------------------------------------------------------
    df_traceplot = data.frame(chain = c(rep(1,ncol(LL)), rep(2,ncol(LL)), rep(3,ncol(LL)),rep(4,ncol(LL))),
                              log_likelihood = c(LLik_sum[[1]],LLik_sum[[2]],LLik_sum[[3]],LLik_sum[[4]]),
                              iterations = rep(1:(20000),4))
    df_traceplot = df_traceplot %>% mutate(chain = factor(chain, levels = 1:4))
    
    my_sexy_traceplot<- ggplot(df_traceplot, aes(x = iterations, y = log_likelihood, color = factor(chain), group=chain))+
      geom_line(alpha = .5)+
      labs(title = "Log likelihood for the 4 chains",
           subtitle = paste0("Number of iterations: ", N_iter," || Burnin: ", burnin), 
           x = "Iterations",
           y = "Log likelihood",
           color = "Chain")+
      theme_bw()
    traceplot_name <- paste0(processed_wd,"//traceplot",est_model, "_K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
    png(traceplot_name,width = 500, height = 250)
    print(my_sexy_traceplot)
    dev.off()
    
    #computing the 
    WAIC_est_1 = waic(t(LL_list[[1]]))
    LOO<-loo(t(LL_list[[1]]))
    
    z_s_table = data.frame(WAIC_est =WAIC_est_1$estimates[3],
                           elpd_waic = WAIC_est_1$estimates[1],
                           p_waic =  WAIC_est_1$estimates[2],
                           WAIC_SE =  WAIC_est_1$estimates[6],
                           looic =LOO$estimates[3],
                           elpd_loo = LOO$estimates[1],
                           p_loo = LOO$estimates[2],
                           loiic_SE = LOO$estimates[6])
    
    
    if(is.simulation == T){
      z_true = uploaded_results$chain1$ground_truth$z
      similarity_matrix <- pr_cc(z_list_relab[[1]])
      point_est_minVI = minVI(similarity_matrix)$cl
      
      z_s_table$VIdist_minVI <- vi.dist(point_est_minVI, z_true)
      z_s_table$VIdist_MAP <- vi.dist(point_est_z, z_true)
      
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
      
      sigma_squared_s_table<- sigma_squared_summary_table(chains = uploaded_results, 
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
      
      mu_vec_s_table<- mu_vec_summary_table(chains = uploaded_results, true_value = is.simulation*(true_model!='Simple'),
                                            diag0.5 = TRUE, K = K, burnin = burnin)
      
      mu_vec_s_table = mu_vec_s_table %>% mutate(model=rep(est_model,nrow(mu_vec_s_table))) %>% mutate(n_clust = rep(K,nrow(mu_vec_s_table)))
      if(est_model=='SST'&file==1){
        mu_vec_container = mu_vec_s_table
      }
      mu_vec_container =  rbind(mu_vec_container,mu_vec_s_table)
      
      
    }
    
    
    
    
    #-------------------------------------------------------------------------------
    # P diagnostics 
    #-------------------------------------------------------------------------------
    # 
    
    P_d_table<- P_diagnostic_table(chains = uploaded_results, true_value = is.simulation, 
                                   permutations_z = permutations_z,
                                   diag0.5 = TRUE,K = K,
                                   P = uploaded_results$chain1$ground_truth$P,
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
      
      mu_vec_d_table <- mu_vec_diagnostic_table(chains = uploaded_results, true_value = is.simulation*(true_model!='Simple'), diag0.5 = TRUE,
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
      P_true <- round(uploaded_results$chain1$ground_truth$P,3) %>% data.frame()
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
    }
    
    
    
    #---------------------------------------------------------------------------
    # COMPARISON WITH THE EXTERNAL RANKING
    #---------------------------------------------------------------------------
    if(is.simulation==F){
      my_names <- read.csv("/Users/lapo_santi/Desktop/Nial/MCMC_results/applications_orderstats/tennis/rawdata/players_df.csv")
      
      plot_name <- paste0(processed_wd,"//RankvsClust_Est_model",est_model, "_K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
      # Save the plot with the constructed file name
      
      g_df =  data.frame(vertex_attr(g)) %>%
        rename(player_slug= name) %>%
        left_join(players_df, by="player_slug") %>%
        mutate(degree_pl = degree(g,mode = 'out')/degree(g,mode = 'all')) %>%
        arrange()
      est_df<- data.frame(player_slug = my_names$player_slug, est_cl = point_est_z)
      combined_df<- inner_join(g_df,est_df,by = 'player_slug')
      png(plot_name,width = 800, height = 627)
      print(rank_vs_cluster(combined_df, combined_df$est_cl,est_model = est_model))
      dev.off()
      # 
      #---------------------------------------------------------------------------
      # CHECKING THE HOMOGENEITY OF THE CLUSTERS: HEATMAP
      #---------------------------------------------------------------------------
      
      my_names <- read.csv("/Users/lapo_santi/Desktop/Nial/MCMC_results/applications_orderstats/tennis/rawdata/players_df.csv")
      
      plot_name <- paste0(processed_wd,"//RankvsClust_Est_model",est_model, "_K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
      
      A = uploaded_results$chain1$Y_ij
      N = uploaded_results$chain1$N_ij
      
      
      
      
      players_considered=rownames(uploaded_results$chain1$Y_ij)
      #creating a dataframe with all the possible combinations among players
      all_possible_combinations=expand.grid(players_considered,players_considered)
      #counting how many times a player has won against every other player
      for(i in 1:nrow(all_possible_combinations)){
        all_possible_combinations$victories[i]<- A[all_possible_combinations$Var1[i],all_possible_combinations$Var2[i]]
      }
      #counting how many times a player has played against every other player
      for(i in 1:nrow(all_possible_combinations)){
        all_possible_combinations$games_played[i]<- N[all_possible_combinations$Var1[i],all_possible_combinations$Var2[i]]
      }
      #computing the percentage of victories as games_won/games_played
      all_possible_combinations = all_possible_combinations %>%
        mutate(percentage_won = victories/games_played)
      #computing the percentage of defeats as (games_played-games_won)/games_played
      all_possible_combinations = all_possible_combinations %>%
        mutate(percentage_lost = (games_played-victories)/games_played)
      #renaming player1 and player2 columns
      all_possible_combinations = all_possible_combinations %>%
        rename(player1= Var1, player2= Var2 )
      
      n_col_df<- ncol(all_possible_combinations)
      #joining the informations about player 1
      all_possible_combinations = all_possible_combinations %>%
        inner_join(combined_df, by = c("player1" = "player_slug")) %>% #renaming varibales with suffix
        rename_at(vars((n_col_df+1):(n_col_df+1+ncol(combined_df)-2)),function(x) paste0(x,"_player1"))
      #joining infos about player 2
      all_possible_combinations = all_possible_combinations %>%
        inner_join(combined_df, by = c("player2" = "player_slug"))
      
      #FULL HEATMAP --------------------------------------------------------------
      combined_plot <- ggplot(all_possible_combinations, aes(x=reorder(player2,est_cl,decreasing = F), y= reorder(player1,est_cl_player1,decreasing = T)))+
        geom_tile(aes(fill= percentage_won), color="grey9")+
        geom_ysidecol(aes(x = degree_pl_player1, color=factor(est_cl_player1))) +
        geom_xsidecol(aes(y = 1-degree_pl, color=factor(est_cl))) +
        scale_fill_gradient(low = "white", high = "red") +
        theme_bw()+ theme(legend.direction = "vertical",
                          axis.text.x = element_blank(),
                          axis.text.y = element_blank())+
        labs(title = 'Heatmap filled with victory percentages',
             x = paste0("Players ordered by blocks"),
             y = paste0("Playersordered by blocks"),
             fill = "% victories",
             color = "Block")
      
      #saving the heatmap
      plot_name1<- paste0(processed_wd,"//Combined_plot",est_model, "_K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
      png(plot_name1,width = 800, height = 594)
      print(combined_plot)
      dev.off()
      
      
      
      #---------------------------------------------------------------------------
      # CHECKING THE RANKING and THE CLUSTERING
      #---------------------------------------------------------------------------
      
      percentage_to_display <- 30
      set.seed(23)
      # Randomly sample a subset of labels to display
      sampled_labels <- combined_df[sample(nrow(combined_df), size = round(percentage_to_display / 100 * nrow(combined_df))), ]
      
      rank_boxplot<- ggplot(combined_df, aes(x = factor(est_cl), y = median_rank,color = factor(est_cl))) +
        geom_boxplot(aes(fill=factor(est_cl)),alpha=.3) +
        geom_label_repel(
          data = sampled_labels,  # Use the sampled labels for display
          aes(label = player_slug),
          size = 3,
          hjust = .5,
          vjust = 0,
          show.legend = F,
          alpha=.8
        ) +
        labs(title= "Rank of the players divided into blocks",
             subtitle = "Not all names are displayed to avoid overlapping",
             x = "Clusters",
             y = "Median Rank 2017",
             color = "Cluster",
             fill = "Cluster")+
        theme_classic()+
        theme(legend.position = "bottom",
              plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
              plot.subtitle = element_text(face = "bold.italic", hjust = 0.5),
              plot.caption = element_text(face = "italic"))
      
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
  
  chosen_model <- z_container[which(z_container$WAIC_est ==  min(z_container$WAIC_est)),]
  
  
  
  uploaded_results<- readRDS(paste0(data_wd, 'True_ModelTennis_dataEst_model_',chosen_model$model,'_N95_K',chosen_model$K,'.RDS'))
  K<- chosen_model$K
  P_est <- apply(uploaded_results$chain1$est_containers$P[,,-c(1:burnin)], MARGIN = c(1,2), mean)
  P_est <- inverse_logit_f(P_est)
  my_z_est<- z_plot(chains = uploaded_results , true_model= true_model,P_est = P_est,
                    est_model = est_model, true_value =is.simulation, 
                    diag0.5 =diag0.5 , K=chosen_model$K, N=nrow(uploaded_results$chain1$Y_ij), z = uploaded_results$chain1$ground_truth$z ,
                    burnin =  burnin ,label_switch = T,tap= processed_wd)
  
  
  point_est<- as.vector(my_z_est$point_est)
  table(point_est)
  
  g_df =  data.frame(vertex_attr(g)) %>%
    rename(player_slug= name) %>%
    left_join(players_df, by="player_slug") %>%
    mutate(degree_pl = degree(g,mode = 'out')/degree(g,mode = 'all')) %>%
    arrange()
  est_df<- data.frame(player_slug = my_names$player_slug, est_cl = point_est)
  combined_df<- inner_join(g_df,est_df,by = 'player_slug')
  
  colnames(A)<- rownames(A)
  colnames(N)<- colnames(A)
  rownames(N)<- rownames(A)
  
  
  my_beauti_plottini = list()
  K<- length(unique( point_est))
  
  for(block_i in 1:K){
    for(block_j in 1:K){
      
      players12 <- est_df %>% mutate(player_in_block_i = est_cl %in% block_i)%>%
        mutate(player_in_block_j = est_cl %in% block_j)
      
      A_subset<- A[players12$player_in_block_i, players12$player_in_block_j]
      N_subset<- N[players12$player_in_block_i, players12$player_in_block_j]
      
      
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
