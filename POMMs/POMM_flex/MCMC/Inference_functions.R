rank_vs_cluster<- function(data_with_statistics, clustering, est_model){
  

  # Create the ggplot plot with error bars and modifications
  main_plot<-ggplot(data_with_statistics, aes(x = reorder(player_slug, median_rank), y = median_rank, color = factor(clustering))) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = min_r, ymax = max_r), size = 1) +
    labs(x = "Player Name", y = "Ranking", title = paste0(est_model," Estimated Block Membership and Ranking"),
         subtitle = 'Players are sorted in ascending order relative to their median ranking in 2017') +
    scale_color_discrete(name = "Cluster") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
      text = element_text(size = 10, family = "Arial"),
      plot.title = element_text(size = 13, face = "bold", margin = margin(r = 10)),
      plot.subtitle = element_text(size = 10, margin = margin(t = 10, r = 10, b = 10)),  # Adjust the top margin
      legend.text = element_text(size = 12),
      plot.margin = margin(20, 20, 20, 20)
    )
  degree_plot <- ggplot(data_with_statistics, aes(x = reorder(player_slug, median_rank), y = degree_pl, fill =factor(clustering) )) +
    geom_bar(stat = "identity") +
    labs(x = "Player Name", y = "Percentage Victories", fill='Cluster', title = "Percentage of victories for each player",
         subtitle = 'Players are sorted in descending order relative to their percentage of victories') +
    scale_color_discrete(name = "Cluster") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
      text = element_text(size = 10, family = "Arial"),
      plot.title = element_text(size = 13, face = "bold", margin = margin(r = 10)),
      plot.subtitle = element_text(size = 10, margin = margin(t = 10, r = 10, b = 10)),  # Adjust the top margin
      plot.margin = margin(20, 20, 20, 20)
    )
  return(plot_grid(main_plot, degree_plot, ncol = 1, align = "v"))}


plot_P = function(p_container, p_true, burnin,K){
  
  burnin_p <- p_container[,,-c(1:burnin)]
  plots = list()
  for(i in 1:K) {
    for(j in 1:K) {
      y_try = data.frame(y = as.vector(burnin_p[i, j,]))
      p1 = ggplot(y_try, aes(y)) +
        geom_density(fill = "dodgerblue", alpha = 0.5) +
        scale_x_log10() +
        geom_vline(xintercept = p_true[i, j], color = "red")+
        xlab("probability") +
        ylab("Density") +
        ggtitle(paste("Density plot of entry ", i, ",", j, sep = ""))
      
      plots[[length(plots) + 1]] <- p1
    }
  }
  p_combined = patchwork::wrap_plots(plots, ncol = K, nrow = K)
  return(p_combined)}



z_permute<-function (z_container, permutations,K) {
  
  m <- dim(permutations)[1]
  K <- dim(permutations)[2]
  J <- nrow(z_container)
  
  t_z<- t(z_container)
  
  z_array<- array(0,dim=c(m,K,J))  
  for(iter in 1:m){
    z_array[iter,,]<- t(vec2mat(t_z[iter,]))
  }
  
  
  mcmc.permuted <- z_array
  for (iter in 1:m) {
    for (j in 1:J) {
      mcmc.permuted[iter, , j] <- z_array[iter, permutations[iter, 
      ], j]
    }
  }
  
  output<- matrix(0,J,m)
  for(iter in 1:m){
    for(j in 1:J){
      output[j,iter]<- which(mcmc.permuted[iter,,j]>0)
    }
  }
  return(output)
}
permute_array <-function(array_samples, perm_matrix) {
  
  N_iter <- dim(array_samples)[3]  # Number of iterations
  K <- dim(array_samples)[1]       # Dimension of the array (K by K)
  
  permuted_array <- array(dim = c(K, K, N_iter))  # Initialize permuted array
  
  for (i in 1:N_iter) {
    perm_indices <- perm_matrix[i, ]  # Permutation indices for the current iteration
    permuted_array[, , i] <- array_samples[perm_indices, perm_indices, i]
  }
  
  return(permuted_array)
}
uncertainty_labels<- function(z_container,K)
{
  N <- nrow(z_container)
  m <- ncol(z_container)
  label_proportions <- matrix(0, N, K)
  colnames(label_proportions)<-1:K
  # Calculate label proportions with zero counts
  for (i in 1:N) {
    label_counts <- table(z_container[i, ])
    
    # Fill in the label proportions matrix
    label_proportions[i, names(label_counts)] <- round(label_counts / m,2)
  }
  return(label_proportions)
}


estimator_P <- function(P_container){
  K<-nrow(P_container)
  P_hat<- matrix(0,K,K )
  for(p in 1:K){
    for(q in 1:K){
      P_hat[p,q]<- mean(P_container[p,q,])
    }
  }
  return(P_hat)
}

calculate_misclassification_rate <- function(N_new, z_new, N, p_true, z_true, sampled_games, labels_available, P_est, z_est) {
  # create empty matrix of edges between the N_new nodes and those in the original network
  Yij_new <- matrix(0, N_new, N)
  
  # simulate the new edges
  for (i in 1:N_new){
    for (j in 1:sampled_games){
      Yij_new[i, j] <- rbinom(1, 1, prob = p_true[z_new[i], z_true[j]])
    }
  }
  
  z_proportion <- table(z_est) / N
  
  K <- length(labels_available)
  z_predicted_prob <- matrix(0, N_new, K)
  
  for (k in labels_available){
    for (i in 1:N_new){
      lik_i <- 0 
      for (j in sampled_games){
        lik_i <- lik_i + dbinom(Yij_new[i, j], 1, P_est[k, z_est[j]], log = TRUE)
      }
      lik_i <- lik_i + log(z_proportion[k])
      z_predicted_prob[i, k] <- lik_i 
    }
  }
  
  misclassification_rate <- (N_new - sum(diag(table(z_new, apply(z_predicted_prob, 1, which.max))))) / N_new
  
  return(misclassification_rate)
}


calculate_waic_matrix <- function(n_ij_matrix, z_container, N_iter, p_container, y_ij_matrix) {
  
  upper.tri.non.zero.waic <- which(n_ij_matrix > 0, arr.ind = TRUE)
  waic_matrix_container <- matrix(0, nrow = N_iter, ncol = length(upper.tri.non.zero.waic))
  
  for (ii in 1:N_iter) {
    z_mat_waic <- vec2mat(z_container[, ii])
    p_ij_waic <- calculate_victory_probabilities(z_mat_waic, p_container[,, ii])
    waic_matrix_container[ii, ] <- dbinom(y_ij_matrix[upper.tri.non.zero.waic], 
                                          size = n_ij_matrix[upper.tri.non.zero.waic], 
                                          p_ij_waic[upper.tri.non.zero.waic])
  }
  
  waic.matrix(waic_matrix_container)
}

MSE_p_matrix = function(burnin, p_container, p_true){
  K<- nrow(p_true)
  burned_p <- p_container[,,-c(1:burnin)]
  mse_table = matrix(0,K,K)
  for(i in 1:K){
    for(j in 1:K){
      mse_table[i,j]= (mean(burned_p[i,j,]) - p_true[i,j])
    }
  }
  return(mse_table)
}


Est_p_matrix= function(burnin, p_container, p_true){
  K<- nrow(p_container[,,1])
  burned_p <- p_container[,,-c(1:burnin)]
  est_table = matrix(0,K,K)
  for(i in 1:K){
    for(j in 1:K){
      est_table[i,j]= mean(burned_p[i,j,])
    }
  }
  return(est_table)
}


save_table_to_file <- function(table_code, filename, title = NULL, subtitle = NULL) {
  table_code = data.frame(table_code) 
  if(all(is.numeric(table_code)) != F){
    table_code %>% select(where(is.numeric)) %>% round(digits = 3)}
  # Create a gt table
  table <- table_code %>%
    gt() %>%
    tab_header(title = title, subtitle = subtitle) %>%
    # Apply a custom LaTeX theme
    tab_options(table.align = 'left') %>%
    as_latex() %>%
    as.character()
  
  # Write the LaTeX table to the file
  writeLines(con = filename, text = table)
}

z_plot<- function(test_output, true_model, est_model, true_value, diag0.5 , K, N, z , burn_in,label_switch ){
  
  z_container_POMM <- assembling_chains(test_output,burn_in,'z')
  
  Y_ij<-test_output$chain1$Yij_matrix
  #extracting similarity matrix
  
  
  if(true_value == T){
    if(label_switch==T){
      runPOMM<- label.switching(method = 'ECR' ,zpivot = z ,z = t(z_container_POMM), K = K)
      z_container_POMM<- z_permute(z_container_POMM, permutations = runPOMM$permutations$ECR)
    }
    similarity_matrixPOMMM = pr_cc(z_container_POMM)
    #plotting it
    plot_name <- paste0("adjacency_",true_model,est_model, "_K",K,"_N",N,".png")
    # Save the plot with the constructed file name
    png(plot_name,width = 800, height = 800)
    similarity_plot(Y_ij, z, z) #checking mixing
    # Close the device to save the plot
    dev.off()
    
    #plotting it
    plot_name <- paste0("similarity_",true_model,est_model, "_K",K,"_N",N,".png")
    # Save the plot with the constructed file name
    png(plot_name,width = 800, height = 800)
    similarity_plot(similarity_matrixPOMMM, z, z) #checking mixing
    # Close the device to save the plot
    dev.off()
  }else{
    if(label_switch==T){
      runPOMM<- label.switching(method = 'DATA-BASED',z = t(z_container_POMM), K = K,data = rowSums(Yij_matrix)/colSums(Yij_matrix))
      z_container_POMM<- z_permute(z_container_POMM, permutations = runPOMM$permutations$`DATA-BASED`)
    }
    similarity_matrixPOMMM = pr_cc(z_container_POMM)
    A_container_POMM <- cbind(test_output$chain1$control_containers$A[-c(1:burnin)],
                              test_output$chain2$control_containers$A[-c(1:burnin)],
                              test_output$chain3$control_containers$A[-c(1:burnin)],
                              test_output$chain4$control_containers$A[-c(1:burnin)]) 
    
    z_MAP_POMM= z_container_POMM[,which(A_container_POMM == max(A_container_POMM))[1]]
    #plotting it
    point_est_POMM = minVI(similarity_matrixPOMMM)$cl
    plot_name <- paste0("adjacency_",true_model,est_model, "_K",K,"_N",N,".png")
    # Save the plot with the constructed file name
    png(plot_name,width = 800, height = 800)
    similarity_plot(Y_ij, z_MAP_POMM, z_MAP_POMM) #checking mixing
    # Close the device to save the plot
    dev.off()
    
    #plotting it
    plot_name <- paste0("similarity_",true_model,est_model,"K",K,"_N",N,".png")
    # Save the plot with the constructed file name
    png(plot_name,width = 800, height = 800)
    
    similarity_plot(similarity_matrixPOMMM, z_MAP_POMM, z_MAP_POMM) #checking mixing
    # Close the device to save the plot
    dev.off()
    return(point_est_POMM)
  }
  
}


z_summary_table<- function(test_output , true_value, diag0.5 , K, burn_in, label_switch = F){
  
  
  z_container_POMM <- assembling_chains(test_output,burnin = burn_in,parameter = 'z')
  
  #Data used to generate the data -----
  if(true_value==T){
    
    if(label_switch==T){
      runPOMM<- label.switching(method = 'ECR' ,zpivot = z ,z = t(z_container_POMM), K = K)
      z_container_POMM<- z_permute(z_container_POMM, permutations = runPOMM$permutations$ECR)
    }
    Yij_matrix<- test_output$chain1$Yij_matrix
    Nij_matrix<- test_output$chain1$Nij_matrix
    
    
    
    
    # Create a data frame to store the results
    results <- data.frame(
      MAP_vi_dist = 0,
      MINVI_vi_dist = 0,
      WAIC_est = 0,
      WAIC_se = 0
    )
    
    A_container_POMM <- cbind(test_output$chain1$control_containers$A[-c(1:burnin)],
                              test_output$chain2$control_containers$A[-c(1:burnin)],
                              test_output$chain3$control_containers$A[-c(1:burnin)],
                              test_output$chain4$control_containers$A[-c(1:burnin)]) 
    z_truePOMM <- test_output$chain1$ground_truth$z #true underlying value
    
    
    similarity_matrixPOMMM<- pr_cc(z_container_POMM)
    #point est 1
    point_est_POMM = minVI(similarity_matrixPOMMM)$cl
    #point est 2
    z_MAP_POMM= z_container_POMM[,which(A_container_POMM == max(A_container_POMM))[1]]
    
    
    results$MAP_vi_dist <- vi.dist(z_MAP_POMM, z_truePOMM)
    results$MINVI_vi_dist<- vi.dist(point_est_POMM, z_truePOMM)
    
    
    #computing WAIC
    
    Pcontainer<- assembling_chains(test_output,burnin = burn_in,parameter = 'P')
    if(label_switch==T){
      Pcontainer<- permute_array(array_samples = Pcontainer, perm_matrix = runPOMM$permutations$ECR)
    }
    WAIC<- calculate_waic_matrix(n_ij_matrix = Nij_matrix,z_container = z_container_POMM,N_iter = ncol(z_container_POMM),p_container = Pcontainer,y_ij_matrix = Yij_matrix )
    results$WAIC_est <- WAIC$estimates[3,1]
    results$WAIC_se <- WAIC$estimates[3,2]
    
    
    
  }else{
    Yij_matrix<- test_output$chain1$Yij_matrix
    Nij_matrix<- test_output$chain1$Nij_matrix
    # apply the permutations returned by typing:
    if(label_switch == T){
      runPOMM<- label.switching(method = 'DATA-BASED',z = t(z_container_POMM), K = K,data = rowSums(Yij_matrix)/colSums(Yij_matrix))
      z_container_POMM<- z_permute(z_container_POMM, permutations = runPOMM$permutations$`DATA-BASED`)
    }
    results <- data.frame(
      WAIC_est = 0,
      WAIC_se = 0
    )
    
    similarity_matrixPOMMM<- pr_cc(z_container_POMM)
    #point est 1
    point_est_POMM = minVI(similarity_matrixPOMMM)$cl
    #point est 2
    
    
    A_container_POMM <- cbind(test_output$chain1$control_containers$A[-c(1:burnin)],
                              test_output$chain2$control_containers$A[-c(1:burnin)],
                              test_output$chain3$control_containers$A[-c(1:burnin)],
                              test_output$chain4$control_containers$A[-c(1:burnin)]) 
    
    z_MAP_POMM= z_container_POMM[,which(A_container_POMM == max(A_container_POMM))[1]]
    
    #computing WAIC
    Yij_matrix<- test_output$chain1$Yij_matrix
    Nij_matrix<- test_output$chain1$Nij_matrix
    
    Pcontainer<- assembling_chains(test_output,burnin = burn_in,parameter = 'P')
    if(label_switch==T){
      Pcontainer<- permute_array(array_samples = Pcontainer, perm_matrix = runPOMM$permutations$`DATA-BASED`)
    }    
    
    #computing WAIC
    WAIC<- calculate_waic_matrix(n_ij_matrix = Nij_matrix,
                                 z_container = z_container_POMM,
                                 N_iter = ncol(z_container_POMM),
                                 p_container = Pcontainer,y_ij_matrix = Yij_matrix )
    
    results$WAIC_est <- WAIC$estimates[3,1]
    results$WAIC_se <- WAIC$estimates[3,2]
    
  }
  return(list(table=results, memb = z_MAP_POMM))
}

z_diagnostic_table<- function(chains, true_value, diag0.5,K,burn_in,N_iter){
  stopifnot(length(chains)==4)
  
  test1<-chains$chain1
  test2<-chains$chain2
  test3<-chains$chain3
  test4<-chains$chain4
  
  N<- nrow(chains$chain1$Yij_matrix)
  if(true_value == F){
    
    results = data.frame(ESS = 0, LAG_30=0, acceptance_rate=0,Gelman_rubin=0)
    
    mm<-mcmc.list(chains_list = mcmc.list(mcmc(t(test1$est_containers$z[,-c(1:burn_in)])),
                                          mcmc(t(test2$est_containers$z[,-c(1:burn_in)])),
                                          mcmc(t(test3$est_containers$z[,-c(1:burn_in)])),
                                          mcmc(t(test4$est_containers$z[,-c(1:burn_in)]))))
    
    mm_A<-mcmc.list(chains_list = mcmc.list(mcmc(test1$control_containers$A[-c(1:burn_in)]),
                                            mcmc(test2$control_containers$A[-c(1:burn_in)]),
                                            mcmc(test3$control_containers$A[-c(1:burn_in)]),
                                            mcmc(test4$control_containers$A[-c(1:burn_in)])))
    #ESS
    results$ESS <- round(sum(simplify2array(lapply(mm, effectiveSize))),0)/N
    #Gelman Rubin
    gelman_vector =vector()
    for(i in 1:N){
      individual_i_chain<- mcmc.list(chains_list = mcmc.list(mcmc((test1$est_containers$z[i,-c(1:burn_in)])),
                                                   mcmc((test2$est_containers$z[i,-c(1:burn_in)])),
                                                   mcmc((test3$est_containers$z[i,-c(1:burn_in)])),
                                                   mcmc((test4$est_containers$z[i,-c(1:burn_in)]))))
      gelman_d<-gelman.diag(individual_i_chain)[[1]]
      gelman_vector<-append(gelman_vector, gelman_d[1])
    }
    results$Gelman_rubin <- median(gelman_vector,na.rm = T)
    #Autocorrelation at lag=30
    results$LAG_30 <- round(mean(simplify2array(lapply(mm_A,autocorr.diag,lag=30))),3)
    
    mm_acc<-list(median(test1$acceptance_rates$acc.count_z),
                 median(test2$acceptance_rates$acc.count_z),
                 median(test3$acceptance_rates$acc.count_z),
                 median(test4$acceptance_rates$acc.count_z))
    
    results$acceptance_rate<- mean(unlist(mm_acc))/N_iter*100
  }else if(true_value == T){
    
    results = data.frame(ESS = 0, LAG_30=0, acceptance_rate=0, Gelman_rubin=0,MAP=0)
    
    mm<-mcmc.list(chains_list = mcmc.list(mcmc(t(test1$est_containers$z[,-c(1:burn_in)])),
                                          mcmc(t(test2$est_containers$z[,-c(1:burn_in)])),
                                          mcmc(t(test3$est_containers$z[,-c(1:burn_in)])),
                                          mcmc(t(test4$est_containers$z[,-c(1:burn_in)]))))
    
    mm_A<-mcmc.list(chains_list = mcmc.list(mcmc(test1$control_containers$A[-c(1:burn_in)]),
                                            mcmc(test2$control_containers$A[-c(1:burn_in)]),
                                            mcmc(test3$control_containers$A[-c(1:burn_in)]),
                                            mcmc(test4$control_containers$A[-c(1:burn_in)])))
    #ESS
    results$ESS <- round(sum(simplify2array(lapply(mm, effectiveSize))),0)/N
    
    gelman_vector =vector()
    for(i in 1:N){
      individual_i_chain<- mcmc.list(chains_list = mcmc.list(mcmc((test1$est_containers$z[i,-c(1:burn_in)])),
                                                             mcmc((test2$est_containers$z[i,-c(1:burn_in)])),
                                                             mcmc((test3$est_containers$z[i,-c(1:burn_in)])),
                                                             mcmc((test4$est_containers$z[i,-c(1:burn_in)]))))
      gelman_d<-gelman.diag(individual_i_chain)[[1]]
      gelman_vector<-append(gelman_vector, gelman_d[1])
    }
    results$Gelman_rubin <- median(gelman_vector,na.rm = T)
    
    #Autocorrelation at lag=30
    results$LAG_30 <- round(mean(simplify2array(lapply(mm_A,autocorr.diag,lag=30))),3)
    
    mm_acc<-list(median(test1$acceptance_rates$acc.count_z),
                 median(test2$acceptance_rates$acc.count_z),
                 median(test3$acceptance_rates$acc.count_z),
                 median(test4$acceptance_rates$acc.count_z))
    
    results$acceptance_rate<- mean(unlist(mm_acc))/N_iter*100
    A_bind = unlist(mm_A)
    brrr <-simplify2array(mm)
    z<- test1$ground_truth$z
    z_bind= rbind(mm[[1]],mm[[2]],mm[[3]],mm[[4]])
    results$MAP <- vi.dist(z_bind[which(A_bind==max(A_bind))[1],],z)
  }
  return(results)}

###-----------------------------------------------------------------------------
# P summary and diagnostics
###


P_summary_table <- function(test_output, true_value, diag0.5, P, K, burn_in, label_switch) {
  
  Yij_matrix <- test_output$chain1$Yij_matrix
  
  j_start <- ifelse(diag0.5, yes = 1, no = 0)
  K_stop <- ifelse(diag0.5, yes = K - 1, no = K)
  
  entries_df <- data.frame(entry_i = 0, entry_j = 0)
  for (ii in 1:K_stop) {
    for (jj in (ii + j_start):K) {
      entries_df <- rbind(entries_df, data.frame(entry_i = ii, entry_j = jj))
    }
  }
  entries_df <- entries_df[-1,]
  
  if (true_value == FALSE) {
    
    Pcontainer <- assembling_chains(test_output, burnin = burn_in, parameter = 'P')
    if (label_switch == TRUE) {
      z_container_POMM <- assembling_chains(test_output, burnin = burn_in, parameter = 'z')
      runPOMM <- label.switching(method = 'DATA-BASED', z = t(z_container_POMM), K = K, data = rowSums(Yij_matrix) / colSums(Yij_matrix))
      Pcontainer <- permute_array(array_samples = Pcontainer, perm_matrix = runPOMM$permutations$`DATA-BASED`)
    }
    
    results <- cbind(entries_df, data.frame(mean_est = rep(0, nrow(entries_df)),
                                            credible_interval_05 = rep(0, nrow(entries_df)),
                                            credible_interval_95 = rep(0, nrow(entries_df))))
    for (i in 1:nrow(results)) {
      m <- mcmc(Pcontainer[results$entry_i[i], results$entry_j[i], ])
      results$mean_est[i] <- mean(m)
      HPD <- round(cbind(coda::HPDinterval(m)), 2)
      results$credible_interval_95[i] <- HPD[2]
      results$credible_interval_05[i] <- HPD[1]
    }
    P_hat <- estimator_P(P_container = Pcontainer)
  } else if (true_value == TRUE) {
    z_true <- test_output$chain1$ground_truth$z
    
    Pcontainer <- assembling_chains(test_output, burnin = burn_in, parameter = 'P')
    
    if (label_switch == TRUE) {
      z_container_POMM <- assembling_chains(test_output, burnin = burn_in, parameter = 'z')
      runPOMM <- label.switching(method = 'ECR', z = t(z_container_POMM), K = K, zpivot = z_true)
      Pcontainer <- permute_array(array_samples = Pcontainer, perm_matrix = runPOMM$permutations$ECR)
    }
    results <- cbind(entries_df, data.frame(mean_est = rep(0, nrow(entries_df)),
                                            credible_interval_05 = rep(0, nrow(entries_df)),
                                            credible_interval_95 = rep(0, nrow(entries_df)),
                                            true_value = rep(0, nrow(entries_df))))
    for (i in 1:nrow(results)) {
      m <- mcmc(Pcontainer[results$entry_i[i], results$entry_j[i], ])
      results$mean_est[i] <- round(mean(m), 4)
      HPD <- round(cbind(coda::HPDinterval(m)), 4)
      results$credible_interval_95[i] <- HPD[2]
      results$credible_interval_05[i] <- HPD[1]
      results$true_value[i] <- round(P[results$entry_i[i], results$entry_j[i]], 4)
    }
    P_hat <- estimator_P(P_container = Pcontainer)
  }
  
  return(list(table = results, P_hat = P_hat))
}



#Diagnostics for P#


P_diagnostic_table<- function(chains, true_value, diag0.5,P,K,burn_in,N_iter){
  stopifnot(length(chains)==4)
  
  test1<-chains$chain1
  test2<-chains$chain2
  test3<-chains$chain3
  test4<-chains$chain4
  
  
  
  j_start = ifelse(diag0.5, yes = 1, no = 0)
  K_stop = ifelse(diag0.5, yes = K-1, no = K)
  
  entries_df <- data.frame(entry_i = 0 ,entry_j =0 )
  for( ii in 1:K_stop){
    for(jj in (ii+j_start):K){
      entries_df <- rbind(entries_df, data.frame(entry_i= ii, entry_j = jj))
    }
  }
  entries_df=entries_df[-1,]   
  
  if(true_value == F){
    results = cbind(entries_df, data.frame(ESS = rep(0,nrow(entries_df)),
                                           LAG_30=rep(0,nrow(entries_df)),
                                           acceptance_rate=rep(0,nrow(entries_df)),
                                           Gelman_rubin=rep(0,nrow(entries_df))))
    for(i in 1:nrow(results)){
      
      mm<-mcmc.list(chains_list = mcmc.list(mcmc(test1$est_containers$P[results$entry_i[i],results$entry_j[i],-c(1:burn_in)]),
                                            mcmc(test2$est_containers$P[results$entry_i[i],results$entry_j[i],-c(1:burn_in)]),
                                            mcmc(test3$est_containers$P[results$entry_i[i],results$entry_j[i],-c(1:burn_in)]),
                                            mcmc(test4$est_containers$P[results$entry_i[i],results$entry_j[i],-c(1:burn_in)])))
      #ESS
      results$ESS[i] <- round(mean(simplify2array(lapply(mm, effectiveSize))),0)
      #Gelman Rubin
      results$Gelman_rubin[i]<-  round(gelman.diag(mm)[1]$psrf[1],3)
      #Autocorrelation at lag=30
      results$LAG_30[i] <- round(mean(simplify2array(lapply(mm,autocorr.diag,lag=30))),3)
      
      mm_acc<-list(test1$acceptance_rates$acc.count_p[results$entry_i[i],results$entry_j[i]],
                   test2$acceptance_rates$acc.count_p[results$entry_i[i],results$entry_j[i]],
                   test3$acceptance_rates$acc.count_p[results$entry_i[i],results$entry_j[i]],
                   test4$acceptance_rates$acc.count_p[results$entry_i[i],results$entry_j[i]])
      
      results$acceptance_rate[i]<- mean(unlist(mm_acc))/N_iter*100
    }
  }else if(true_value == T){
    results = cbind(entries_df, data.frame(ESS = rep(0,nrow(entries_df)),
                                           LAG_30=rep(0,nrow(entries_df)),
                                           acceptance_rate=rep(0,nrow(entries_df)),
                                           Gelman_rubin=rep(0,nrow(entries_df)),
                                           MAE = rep(0,nrow(entries_df))))
    for(i in 1:nrow(results)){
      
      mm<-mcmc.list(chains_list = mcmc.list(mcmc(test1$est_containers$P[results$entry_i[i],results$entry_j[i],-c(1:burn_in)]),
                                            mcmc(test2$est_containers$P[results$entry_i[i],results$entry_j[i],-c(1:burn_in)]),
                                            mcmc(test3$est_containers$P[results$entry_i[i],results$entry_j[i],-c(1:burn_in)]),
                                            mcmc(test4$est_containers$P[results$entry_i[i],results$entry_j[i],-c(1:burn_in)])))
      #ESS
      results$ESS[i] <- round(mean(simplify2array(lapply(mm, effectiveSize))),0)
      #Gelman Rubin
      results$Gelman_rubin[i]<-  round(gelman.diag(mm)[1]$psrf[1],3)
      #Autocorrelation at lag=30
      results$LAG_30[i] <- round(mean(simplify2array(lapply(mm,autocorr.diag,lag=30))),3)
      
      mm_acc<-list(test1$acceptance_rates$acc.count_p[results$entry_i[i],results$entry_j[i]],
                   test2$acceptance_rates$acc.count_p[results$entry_i[i],results$entry_j[i]],
                   test3$acceptance_rates$acc.count_p[results$entry_i[i],results$entry_j[i]],
                   test4$acceptance_rates$acc.count_p[results$entry_i[i],results$entry_j[i]])
      
      results$acceptance_rate[i]<- mean(unlist(mm_acc))/N_iter*100
      results$MAE[i]=round(abs(mean(simplify2array(lapply(mm, mean))) - P[results$entry_i[i],results$entry_j[i]]),4)
    }
  }
  return(results)}

#S inference and diagnostics
S_summary_table<- function(test_output, true_value, diag0.5,S,K,burn_in){
  
  MCMC_samples<- c(test_output$chain1$est_containers$S[-c(1:burn_in)],
                   test_output$chain2$est_containers$S[-c(1:burn_in)],
                   test_output$chain3$est_containers$S[-c(1:burn_in)],
                   test_output$chain4$est_containers$S[-c(1:burn_in)])
  
  m<-mcmc(MCMC_samples)
  if(true_value == F){
    results = data.frame(mean_est = 0, credible_interval_95 =0)
    
    results$mean_est <- mean(m)
    HPD <- round(cbind(coda::HPDinterval(m)),2)
    results$credible_interval_95<- paste0("[",HPD[1],",",HPD[2],"]")
  }else if(true_value == T){
    results = data.frame(mean_est = 0, credible_interval_95 =0, true_value =0)
    
    results$mean_est<- round(mean(m),4)
    HPD <- round(cbind(coda::HPDinterval(m)),4)
    results$credible_interval_95<- paste0("[",HPD[1],",",HPD[2],"]")
    results$true_value<- ifelse(S=='na','na',round(S,4))
    
  }
  return(results)}



S_diagnostic_table<- function(chains, true_value, diag0.5,S,K,burn_in,N_iter){
  stopifnot(length(chains)==4)
  
  test1<-chains$chain1
  test2<-chains$chain2
  test3<-chains$chain3
  test4<-chains$chain4
  
  mm<-mcmc.list(chains_list = mcmc.list(mcmc(test1$est_containers$S[-c(1:burn_in)]),
                                        mcmc(test2$est_containers$S[-c(1:burn_in)]),
                                        mcmc(test3$est_containers$S[-c(1:burn_in)]),
                                        mcmc(test4$est_containers$S[-c(1:burn_in)])))
  
  if(true_value == F){
    results = data.frame(ESS = 0, LAG_30=0, acceptance_rate=0,Gelman_rubin=0)
    
    
    #ESS
    results$ESS <- round(mean(simplify2array(lapply(mm, effectiveSize))),0)
    #Gelman Rubin
    results$Gelman_rubin<-  round(gelman.diag(mm)[1]$psrf[1],3)
    #Autocorrelation at lag=30
    results$LAG_30 <- round(mean(simplify2array(lapply(mm,autocorr.diag,lag=30))),3)
    
    mm_acc<-list(test1$acceptance_rates$acc.count_S,
                 test2$acceptance_rates$acc.count_S,
                 test3$acceptance_rates$acc.count_S,
                 test4$acceptance_rates$acc.count_S)
    
    results$acceptance_rate<- mean(unlist(mm_acc))/N_iter*100
  }else if(true_value == T){
    results = data.frame(ESS = 0,
                         LAG_30=0,
                         acceptance_rate=0,
                         Gelman_rubin=0,
                         MAE = 0)
    
    
    
    #ESS
    results$ESS <- round(mean(simplify2array(lapply(mm, effectiveSize))),0)
    #Gelman Rubin
    results$Gelman_rubin<-  round(gelman.diag(mm)[1]$psrf[1],3)
    #Autocorrelation at lag=30
    results$LAG_30 <- round(mean(simplify2array(lapply(mm,autocorr.diag,lag=30))),3)
    
    mm_acc<-list(test1$acceptance_rates$acc.count_S,
                 test2$acceptance_rates$acc.count_S,
                 test3$acceptance_rates$acc.count_S,
                 test4$acceptance_rates$acc.count_S)
    
    results$acceptance_rate<- mean(unlist(mm_acc))/N_iter*100
    S<- test1$ground_truth$S
    results$MAE=round(abs(mean(simplify2array(lapply(mm, mean))) - S),4)
    
  }
  return(results)}
#alpha inference and diagnosics
alpha_summary_table<- function(test_output, true_value, diag0.5,alpha,K,burn_in){
  alpha_all_chain <- c(test_output$chain1$est_containers$alpha[-c(1:burn_in)],
                       test_output$chain2$est_containers$alpha[-c(1:burn_in)],
                       test_output$chain3$est_containers$alpha[-c(1:burn_in)],
                       test_output$chain4$est_containers$alpha[-c(1:burn_in)])
  
  MCMC_samples = alpha_all_chain
  m<-mcmc(MCMC_samples)
  if(true_value == F){
    results = data.frame(mean_est = 0, credible_interval_95 =0)
    results$mean_est <- mean(m)
    HPD <- round(cbind(coda::HPDinterval(m)),2)
    results$credible_interval_95<- paste0("[",HPD[1],",",HPD[2],"]")
  }else if(true_value == T){
    results = data.frame(mean_est = 0, credible_interval_95 =0, true_value =0)
    results$mean_est<- round(mean(m),4)
    HPD <- round(cbind(coda::HPDinterval(m)),4)
    results$credible_interval_95<- paste0("[",HPD[1],",",HPD[2],"]")
    results$true_value<- round(alpha,4)
    
  }
  return(results)}

alpha_diagnostic_table<- function(chains, true_value, diag0.5,K,burn_in,N_iter){
  stopifnot(length(chains)==4)
  
  test1<-chains$chain1
  test2<-chains$chain2
  test3<-chains$chain3
  test4<-chains$chain4
  mm<-mcmc.list(chains_list = mcmc.list(mcmc(test1$est_containers$alpha[-c(1:burn_in)]),
                                        mcmc(test2$est_containers$alpha[-c(1:burn_in)]),
                                        mcmc(test3$est_containers$alpha[-c(1:burn_in)]),
                                        mcmc(test4$est_containers$alpha[-c(1:burn_in)])))
  
  if(true_value == F){
    results = data.frame(ESS = 0, LAG_30=0,acceptance_rate=0, Gelman_rubin=0 )
    
    
    #ESS
    results$ESS <- round(mean(simplify2array(lapply(mm, effectiveSize))),0)
    #Gelman Rubin
    results$Gelman_rubin<-  round(gelman.diag(mm)[1]$psrf[1],3)
    #Autocorrelation at lag=30
    results$LAG_30 <- round(mean(simplify2array(lapply(mm,autocorr.diag,lag=30))),3)
    
    mm_acc<-list(test1$acceptance_rates$acc.count_alpha,
                 test2$acceptance_rates$acc.count_alpha,
                 test3$acceptance_rates$acc.count_alpha,
                 test4$acceptance_rates$acc.count_alpha)
    
    results$acceptance_rate<- mean(unlist(mm_acc))/N_iter*100
  }else if(true_value == T){
    results = data.frame(ESS = 0,
                         LAG_30=0,
                         acceptance_rate=0,
                         Gelman_rubin=0,
                         MAE = 0)
    
    
    #ESS
    results$ESS <- round(mean(simplify2array(lapply(mm, effectiveSize))),0)
    #Gelman Rubin
    results$Gelman_rubin<-  round(gelman.diag(mm)[1]$psrf[1],3)
    #Autocorrelation at lag=30
    results$LAG_30 <- round(mean(simplify2array(lapply(mm,autocorr.diag,lag=30))),3)
    
    mm_acc<-list(test1$acceptance_rates$acc.count_alpha,
                 test2$acceptance_rates$acc.count_alpha,
                 test3$acceptance_rates$acc.count_alpha,
                 test4$acceptance_rates$acc.count_alpha)
    
    results$acceptance_rate<- mean(unlist(mm_acc))/N_iter*100
    alpha<- test1$ground_truth$alpha
    results$MAE=round(abs(mean(simplify2array(lapply(mm, mean))) - alpha),4)
    
  }
  return(results)}

assembling_chains <- function(chains, burnin, parameter){
  test1<-chains[[1]]
  test2<-chains[[2]]
  test3<-chains[[3]]
  test4<-chains[[4]]
  if(parameter == 'S'){
    assembled<- c(test1$est_containers$S[-c(1:burnin)],test2$est_containers$S[-c(1:burnin)],test3$est_containers$S[-c(1:burnin)],test4$est_containers$S[-c(1:burnin)])
    return(assembled)
  }else if(parameter == 'alpha'){
    assembled<- c(test1$est_containers$alpha[-c(1:burnin)],test2$est_containers$alpha[-c(1:burnin)],test3$est_containers$alpha[-c(1:burnin)],test4$est_containers$alpha[-c(1:burnin)])
    return(assembled)
  }else if(parameter == 'P'){
    assembled<- abind::abind(test1$est_containers$P[,,-c(1:burnin)],test2$est_containers$P[,,-c(1:burnin)],test3$est_containers$P[,,-c(1:burnin)],test4$est_containers$P[,,-c(1:burnin)],along = 3)
    return(assembled)
  }else if(parameter == 'z'){
    assembled<- cbind(test1$est_containers$z[,-c(1:burnin)],test2$est_containers$z[,-c(1:burnin)],test3$est_containers$z[,-c(1:burnin)],test4$est_containers$z[,-c(1:burnin)])
    return(assembled)
  }else{
    print('Please provide a valide name for the parameters: alpha,S,P,Z')
  }
}

