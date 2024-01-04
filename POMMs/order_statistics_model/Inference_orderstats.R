rank_vs_cluster<- function(data_with_statistics, clustering, est_model){
  
  
  # Create the ggplot plot with error bars and modifications
  # main_plot<-ggplot(data_with_statistics, aes(x = reorder(player_slug, median_rank), y = median_rank, color = factor(clustering))) +
  #   geom_point(size = 3) +
  #   geom_errorbar(aes(ymin = min_r, ymax = max_r), size = 1) +
  #   labs(x = "Player Name", y = "Ranking", title = paste0(est_model," Estimated Block Membership and Ranking"),
  #        subtitle = 'Players are sorted in ascending order relative to their median ranking in 2017') +
  #   scale_color_discrete(name = "Cluster") +
  #   theme_bw() +
  #   theme(
  #     axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
  #     text = element_text(size = 10, family = "Arial"),
  #     plot.title = element_text(size = 13, face = "bold", margin = margin(r = 10)),
  #     plot.subtitle = element_text(size = 10, margin = margin(t = 10, r = 10, b = 10)),  # Adjust the top margin
  #     legend.text = element_text(size = 12),
  #     plot.margin = margin(20, 20, 20, 20)
  #   )
  degree_plot <- ggplot(data_with_statistics, aes(x = reorder(player_slug, degree_pl,decreasing = T), y = degree_pl, fill =factor(clustering) )) +
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
  # return(plot_grid(main_plot, degree_plot, ncol = 1, align = "v"))}
  return(degree_plot)}

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
  P=matrix(0,K,K)
  for(iter in 1:m){
    z_array[iter,,]<- t(vec2mat_0_P(t_z[iter,],P = P))
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


estimator_P <- function(Pcontainer){
  K<-nrow(Pcontainer)
  P_hat<- matrix(0,K,K )
  for(p in 1:K){
    for(q in 1:K){
      P_hat[p,q]<- mean(Pcontainer[p,q,])
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


LL_edges <- function(N_ij, Y_ij, z, P){
  z_mat = vec2mat_0_P(z,P)
  P_ij<- calculate_victory_probabilities(z_mat,P)
  ll_lik = dbinom(x = Y_ij[upper.tri(Y_ij,diag=F)], size =  N_ij[upper.tri(N_ij,diag=F)], prob = P_ij[upper.tri(P_ij,diag=F)], log=T)
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


z_plot<- function(chains, true_model, est_model, true_value, diag0.5 , K, N, z , burnin,label_switch,tap){
  
  z_container_POMM <- assembling_chains(chains,burnin,'z')
  z<- chains$chain1$ground_truth$z
  Y_ij<-chains$chain1$Y_ij
  
  
  if(true_value == T){
    if(label_switch==T){
      runPOMM<- label.switching(method = 'ECR' ,zpivot = z ,z = t(z_container_POMM), K = K)
      z_container_POMM<- z_permute(z_container_POMM, permutations = runPOMM$permutations$ECR)
    }
    similarity_matrixPOMM = pr_cc(z_container_POMM)
    #plotting it
    plot_name <- paste0(tap,"//adjacency_",true_model,est_model, "_K",K,"_N",N,".png")
    # Save the plot with the constructed file name
    png(plot_name,width = 800, height = 800)
    similarity_plot(Y_ij, z, z) #checking mixing
    # Close the device to save the plot
    dev.off()
    
    #plotting it
    plot_name <- paste0(tap,"//similarity_",true_model,est_model, "_K",K,"_N",N,".png")
    # Save the plot with the constructed file name
    png(plot_name,width = 800, height = 800)
    similarity_plot(similarity_matrixPOMM, z, z) #checking mixing
    # Close the device to save the plot
    dev.off()
    
  }else{
    if(label_switch==T){
      runPOMM<- label.switching(method = 'DATA-BASED',z = t(z_container_POMM), K = K,data = rowSums(Y_ij)/colSums(Y_ij))
      z_container_POMM<- z_permute(z_container_POMM, permutations = runPOMM$permutations$`DATA-BASED`)
    }
    similarity_matrixPOMM = pr_cc(z_container_POMM)
    A_container_POMM <- cbind(chains$chain1$control_containers$A[-c(1:burnin)],
                              chains$chain2$control_containers$A[-c(1:burnin)],
                              chains$chain3$control_containers$A[-c(1:burnin)],
                              chains$chain4$control_containers$A[-c(1:burnin)]) 
    
    z_MAP_POMM= z_container_POMM[,which(A_container_POMM == max(A_container_POMM))[1]]
    #plotting it
    point_est_POMM = minVI(similarity_matrixPOMM)$cl
    plot_name <- paste0(tap,"//adjacency_",true_model,est_model, "_K",K,"_N",N,".png")
    # Save the plot with the constructed file name
    png(plot_name,width = 800, height = 800)
    similarity_plot(Y_ij, z_MAP_POMM, z_MAP_POMM) #checking mixing
    # Close the device to save the plot
    dev.off()
    
    #plotting it
    plot_name <- paste0(tap,"//similarity_",true_model,est_model,"K",K,"_N",N,".png")
    # Save the plot with the constructed file name
    png(plot_name,width = 800, height = 800)
    
    similarity_plot(similarity_matrixPOMM, z_MAP_POMM, z_MAP_POMM) #checking mixing
    # Close the device to save the plot
    dev.off()
    return(point_est_POMM)
  }
  
}

z_summary_table<- function(chains , true_value, diag0.5 , K, burnin, label_switch = F,tap){
 
  z_container_POMM <- assembling_chains(chains,burnin = burnin,parameter = 'z')
  
  Y_ij<- chains$chain1$Y_ij
  N_ij<- chains$chain1$N_ij
  n =nrow(Y_ij)
  #Data used to generate the data -----
  
  Pcontainer <- assembling_chains(chains, burnin = burnin, parameter = 'P')
  if(label_switch==T){
    runPOMM<- label.switching(method = 'ECR' ,zpivot = z ,z = t(z_container_POMM), K = K)
    z_container_POMM<- z_permute(z_container_POMM, permutations = runPOMM$permutations$`ECR-1`, K = K)
    Pcontainer<- permute_array(array_samples = Pcontainer, perm_matrix = runPOMM$permutations$`ECR-1`)
  }
  
  A_container_POMM <- c(chains$chain1$control_containers$A[-c(1:burnin)],
                        chains$chain2$control_containers$A[-c(1:burnin)],
                        chains$chain3$control_containers$A[-c(1:burnin)],
                        chains$chain4$control_containers$A[-c(1:burnin)])
  
  z_truePOMM <- chains$chain1$ground_truth$z #true underlying value
  similarity_matrixPOMM<- pr_cc(z_container_POMM)
  #point est 1
  point_est_POMM = minVI(similarity_matrixPOMM)$cl
  #point est 2
  z_MAP_POMM= z_container_POMM[,which(A_container_POMM == max(A_container_POMM))[1]]
  
  if(true_value==T){ 
    results = data.frame(MAP_vi_dist = as.numeric(vi.dist(z_MAP_POMM, z_truePOMM)))
    results = results %>% mutate(MINVI_vi_dist= vi.dist(point_est_POMM, z_truePOMM))
  }else{
    results = data.frame(WAIC_est= 0)
  }
  #computing WAIC
  
  Pcontainer<- assembling_chains(chains,burnin = burnin,parameter = 'P')
  
  
  if(true_value==F & label_switch == T){
    # apply the permutations returned by typing:
    
    runPOMM<- label.switching(method = 'DATA-BASED',z = t(z_container_POMM), K = K,data = rowSums(Y_ij)/colSums(Y_ij))
    z_container_POMM<- z_permute(z_container_POMM, permutations = runPOMM$permutations$`DATA-BASED`)
    Pcontainer<- permute_array(array_samples = Pcontainer, perm_matrix = runPOMM$permutations$`DATA-BASED`)
  }   
  
  
  #computing WAIC
  
  LL <- matrix(nrow=n*(n-1)/2,ncol=N_iter-burnin)
  
  for(t in 1:ncol(LL)){
    LL[,t]<- LL_edges(N_ij, Y_ij, z_container_POMM[,t], Pcontainer[,,t])
    if(t%%3000==0){print(paste0("iteration number----->",t))}
  }
  
  
  results = results %>% mutate(WAIC_est = WAIC(LL)$WAIC)
  index_traceplot <- sample(c(1:(n*(n-1)/2)),1)
  
  plot_name <- paste0(tap,"//traceplot_",true_model,est_model, "_K",K,"_N",n,".png")
  # Save the plot with the constructed file name
  png(plot_name,width = 600, height = 326)
  my_plot = plot(ts(LL[index_traceplot,]),xlab="",ylab="")
  dev.off()
  return(list(table=results, memb = z_MAP_POMM,my_plot = my_plot ))
}

z_diagnostic_table<- function(chains, true_value, diag0.5,K,burnin,N_iter,label_switch){
  stopifnot(length(chains)==4)
  
  
  N<- nrow(chains$chain1$Y_ij)
  
  
  results = data.frame(ESS = 0, LAG_30=0, acceptance_rate=0,Gelman_rubin=0)
  
  mm<-mcmc.list(chains_list = mcmc.list(mcmc(t(chains$chain1$est_containers$z[,-c(1:burnin)])),
                                        mcmc(t(chains$chain2$est_containers$z[,-c(1:burnin)])),
                                        mcmc(t(chains$chain3$est_containers$z[,-c(1:burnin)])),
                                        mcmc(t(chains$chain4$est_containers$z[,-c(1:burnin)]))))
  
  mm_A<-mcmc.list(chains_list = mcmc.list(mcmc(chains$chain1$control_containers$A[-c(1:burnin)]),
                                          mcmc(chains$chain2$control_containers$A[-c(1:burnin)]),
                                          mcmc(chains$chain3$control_containers$A[-c(1:burnin)]),
                                          mcmc(chains$chain4$control_containers$A[-c(1:burnin)])))
  #ESS
  results$ESS <- round(sum(simplify2array(lapply(mm, effectiveSize))),0)/N
  #Gelman Rubin
  gelman_vector =vector()
  for(i in 1:N){
    individual_i_chain<- mcmc.list(chains_list = mcmc.list(mcmc((chains$chain1$est_containers$z[i,-c(1:burnin)])),
                                                           mcmc((chains$chain2$est_containers$z[i,-c(1:burnin)])),
                                                           mcmc((chains$chain3$est_containers$z[i,-c(1:burnin)])),
                                                           mcmc((chains$chain4$est_containers$z[i,-c(1:burnin)]))))
    gelman_d<-gelman.diag(individual_i_chain)[[1]]
    gelman_vector<-append(gelman_vector, gelman_d[1])
  }
  results$Gelman_rubin <- median(gelman_vector,na.rm = T)
  #Autocorrelation at lag=30
  results$LAG_30 <- round(mean(simplify2array(lapply(mm_A,autocorr.diag,lag=30))),3)
  
  mm_acc<-list(median(chains$chain1$acceptance_rates$acc.count_z),
               median(chains$chain2$acceptance_rates$acc.count_z),
               median(chains$chain3$acceptance_rates$acc.count_z),
               median(chains$chain4$acceptance_rates$acc.count_z))
  
  results$acceptance_rate<- mean(unlist(mm_acc))/N_iter*100
  if(true_value == T){
    
    if(label_switch==T){
      for(i in 1:length(mm)){
        
        z<- chains$chain1$ground_truth$z
        z_container_POMM <- t(mm[[i]])
        runPOMM<- label.switching(method = 'ECR' ,zpivot = z ,z = t(z_container_POMM), K = K)
        mm[[i]]<- t(z_permute(z_container_POMM, permutations = runPOMM$permutations$ECR))
      }
    }
    
    A_bind = unlist(mm_A)
    z<- chains$chain1$ground_truth$z
    z_bind= rbind(mm[[1]],mm[[2]],mm[[3]],mm[[4]])
    results$MAP <- vi.dist(z_bind[which(A_bind==max(A_bind))[1],],z)
  }
  return(results)
}

###-----------------------------------------------------------------------------
# P summary and diagnostics
###


P_summary_table <- function(chains, true_value, diag0.5, P, K, burnin, label_switch){
  
  Y_ij <- chains$chain1$Y_ij
  
  entries_df <- data.frame(entry_i = 0 ,entry_j =0 )
  for( ii in 1:K){
    for(jj in ii:K){
      entries_df <- rbind(entries_df, data.frame(entry_i= ii, entry_j = jj))
    }
  }
  entries_df=entries_df[-1,]   
  
  
  if (label_switch == TRUE){
    P_true_switched = list()
    P_switched = list()
    for(i in 1:length(chains)){
      z_container <- chains[[i]]$est_containers$z
      Pcontainer<- chains[[i]]$est_containers$P
      if(true_value == T){
        
        P_true_rep <- array(chains[[i]]$ground_truth$P,c(dim(Pcontainer)[1],
                                                         dim(Pcontainer)[2],
                                                         dim(Pcontainer)[3]))
        runPOMM <- label.switching(method = 'ECR', z = t(z_container), K = K, zpivot = chains[[i]]$ground_truth$z)
        P_permuted_container <- permute_array(array_samples = Pcontainer, perm_matrix = runPOMM$permutations$ECR)
        P_true_switched_i<- permute_array(array_samples = P_true_rep, perm_matrix = runPOMM$permutations$ECR)
        
        P_switched[[i]]<- P_permuted_container
        P_true_switched[[i]]<- P_true_switched_i
        
        
        
      }else if(true_value == F){
        runPOMM <- label.switching(method = 'DATA-BASED',  z = t(z_container), K = K, 
                                   data = rowSums(chains$chain1$Y_ij)/rowSums(chains$chain1$N_ij))
        P_permuted_container <- permute_array(array_samples = Pcontainer, perm_matrix = runPOMM$permutations$`DATA-BASED`)
        
        P_switched[[i]]<- P_permuted_container
      }
    }
    
    Pcontainer = abind(P_switched[[1]],P_switched[[2]],P_switched[[3]],P_switched[[4]], along=3)
    P_true_switched =  abind(P_true_switched[[1]],P_true_switched[[2]],P_true_switched[[3]],P_true_switched[[4]], along=3)
  }else{
    Pcontainer <- assembling_chains(chains, burnin = burnin, parameter = 'P')
  }
  
  
  results <- cbind(entries_df, data.frame(mean_est = rep(0, nrow(entries_df)),
                                          credible_interval_05 = rep(0, nrow(entries_df)),
                                          credible_interval_95 = rep(0, nrow(entries_df))))
  
  for (i in 1:nrow(results)){

    m <- mcmc(Pcontainer[results$entry_i[i], results$entry_j[i], ])
    results$mean_est[i] <- mean(m)
    HPD <- round(cbind(coda::HPDinterval(m)), 2)
    results$credible_interval_95[i] <- HPD[2]
    results$credible_interval_05[i] <- HPD[1]
    if (true_value == TRUE) {
      if(label_switch == T){
        results$MAE[i] <- round(mean(Pcontainer[results$entry_i[i], results$entry_j[i],] - P_true_switched[results$entry_i[i], results$entry_j[i],]), 4)
      }else{
        
        results$MAE[i] <- round(mean(Pcontainer[results$entry_i[i], results$entry_j[i],] - P[results$entry_i[i], results$entry_j[i]]), 4)
      }
    }
  }
  P_hat <- estimator_P(Pcontainer = Pcontainer)
  
  
  return(list(table = results, P_hat = P_hat))
}



#Diagnostics for P#

P_diagnostic_table<- function(chains, true_value, diag0.5,P,K,burnin,N_iter, label_switch){
  stopifnot(length(chains)==4)
  
  
  entries_df <- data.frame(entry_i = 0 ,entry_j =0 )
  for( ii in 1:K){
    for(jj in ii:K){
      entries_df <- rbind(entries_df, data.frame(entry_i= ii, entry_j = jj))
    }
  }
  entries_df=entries_df[-1,]   
  
  
  results = cbind(entries_df, data.frame(ESS = rep(0,nrow(entries_df)),
                                         LAG_30=rep(0,nrow(entries_df)),
                                         acceptance_rate=rep(0,nrow(entries_df)),
                                         Gelman_rubin=rep(0,nrow(entries_df))))
  
  if (label_switch == TRUE) {
    
    P_true_switched = list()
    P_switched = list()
    for(i in 1:length(chains)){
      z_container <- chains[[i]]$est_containers$z
      Pcontainer<- chains[[i]]$est_containers$P
      
      if(true_value == TRUE){
        P_true_rep <- array(chains[[i]]$ground_truth$P,c(dim(Pcontainer)[1],
                                                         dim(Pcontainer)[2],
                                                         dim(Pcontainer)[3]))
        runPOMM <- label.switching(method = 'ECR', z = t(z_container), K = K, zpivot = chains[[i]]$ground_truth$z)
        Pcontainer <- permute_array(array_samples = Pcontainer, perm_matrix = runPOMM$permutations$ECR)
        P_true_switched_i<- permute_array(array_samples = P_true_rep, perm_matrix = runPOMM$permutations$ECR)
        P_switched[[i]]<-Pcontainer
        P_true_switched[[i]]<- P_true_switched_i
      }else if(true_value == F){
        
        runPOMM <- label.switching(method = 'DATA-BASED',  z = t(z_container), K = K, 
                                   data = rowSums(chains$chain1$Y_ij)/rowSums(chains$chain1$N_ij))
        Pcontainer <- permute_array(array_samples = Pcontainer, perm_matrix = runPOMM$permutations$`DATA-BASED`)
        P_switched[[i]]<-Pcontainer
      }
    }
  }
  
  P_par_mcmc = array(NA, dim=c((N_iter*0.75), 4, (K*(K-1))/2 +K ))
  my_var_names = vector()
  for(i in 1:nrow(results)){
    
    if(label_switch == F){
      mm<-mcmc.list(mcmc(chains$chain1$est_containers$P[results$entry_i[i],results$entry_j[i],-c(1:burnin)]),
                    mcmc(chains$chain2$est_containers$P[results$entry_i[i],results$entry_j[i],-c(1:burnin)]),
                    mcmc(chains$chain3$est_containers$P[results$entry_i[i],results$entry_j[i],-c(1:burnin)]),
                    mcmc(chains$chain4$est_containers$P[results$entry_i[i],results$entry_j[i],-c(1:burnin)]))
    }else if(label_switch == T){
      mm<-mcmc.list(mcmc(P_switched[[1]][results$entry_i[i],results$entry_j[i],-c(1:burnin)]),
                    mcmc(P_switched[[2]][results$entry_i[i],results$entry_j[i],-c(1:burnin)]),
                    mcmc(P_switched[[3]][results$entry_i[i],results$entry_j[i],-c(1:burnin)]),
                    mcmc(P_switched[[4]][results$entry_i[i],results$entry_j[i],-c(1:burnin)]))
    }
    
    P_par_mcmc[,1,i]<- mcmc(mm[[1]])
    P_par_mcmc[,2,i]<- mcmc( mm[[2]])
    P_par_mcmc[,3,i]<- mcmc(mm[[3]])
    P_par_mcmc[,4,i]<- mcmc(mm[[4]])
    
    my_var_names <- append(my_var_names, paste0("P_",results$entry_i[i],results$entry_j[i]))
  }
  
  P_list = mcmc.list( mcmc(P_par_mcmc[,1,]), 
                      mcmc(P_par_mcmc[,2,]),
                      mcmc(P_par_mcmc[,3,]),
                      mcmc(P_par_mcmc[,4,]))
  varnames(P_list)<- my_var_names
  
  
  # convergence diagnostics statistics -----------------------------------------
  
  #ESS
  results$ESS[i] <- round(mean(simplify2array(lapply(mm, effectiveSize))),0)
  #Gelman Rubin
  results$Gelman_rubin[i]<-  round(gelman.diag(mm)[1]$psrf[1],3)
  #Autocorrelation at lag=30
  results$LAG_30[i] <- round(mean(simplify2array(lapply(mm,autocorr.diag,lag=30))),3)
  
  mm_acc<-list(chains$chain1$acceptance_rates$acc.count_P[results$entry_i[i],results$entry_j[i]],
               chains$chain2$acceptance_rates$acc.count_P[results$entry_i[i],results$entry_j[i]],
               chains$chain3$acceptance_rates$acc.count_P[results$entry_i[i],results$entry_j[i]],
               chains$chain4$acceptance_rates$acc.count_P[results$entry_i[i],results$entry_j[i]])
  
  results$acceptance_rate[i]<- mean(unlist(mm_acc))/N_iter*100
  
  if(true_value == T){
    if(label_switch == T){
      results$MAE[i]=round(abs(mean(simplify2array(lapply(mm, mean))) - P_true_switched[[1]][results$entry_i[i],results$entry_j[i],]),4)
    }
    results$MAE[i]=round(abs(mean(simplify2array(lapply(mm, mean))) - P[results$entry_i[i],results$entry_j[i]]),4)
  }
  
  
  return(list(results = results, 
              plots_list = P_list))
}

#_squared inference and diagnostics
sigma_squared_summary_table<- function(chains, true_value, diag0.5,K,burnin){
  
  MCMC_samples<- c(chains$chain1$est_containers$sigma_squared[-c(1:burnin)],
                   chains$chain2$est_containers$sigma_squared[-c(1:burnin)],
                   chains$chain3$est_containers$sigma_squared[-c(1:burnin)],
                   chains$chain4$est_containers$sigma_squared[-c(1:burnin)])
  
  m<-mcmc(MCMC_samples)
  
  results = data.frame(mean_est = 0, credible_interval_95 =0)
  
  results$mean_est <- mean(m)
  HPD <- round(cbind(coda::HPDinterval(m)),2)
  results$credible_interval_95<- paste0("[",HPD[1],",",HPD[2],"]")
  if(true_value == T){
    
    results$true_value<- chains$chain1$ground_truth$sigma_squared
    
  }
  return(results)}



sigma_squared_diagnostic_table<- function(chains, true_value, diag0.5,K,burnin,N_iter){
  stopifnot(length(chains)==4)
  
  
  mm<-mcmc.list(chains_list = mcmc.list(mcmc(chains$chain1$est_containers$sigma_squared[1,-c(1:burnin)]),
                                        mcmc(chains$chain2$est_containers$sigma_squared[1,-c(1:burnin)]),
                                        mcmc(chains$chain3$est_containers$sigma_squared[1,-c(1:burnin)]),
                                        mcmc(chains$chain4$est_containers$sigma_squared[1,-c(1:burnin)])))
  
  
  results = data.frame(ESS = 0, LAG_30=0, acceptance_rate=0,Gelman_rubin=0)
  
  
  #ESS
  results$ESS <- round(mean(simplify2array(lapply(mm, effectiveSize))),0)
  #Gelman Rubin
  results$Gelman_rubin<-  round(gelman.diag(mm)[1]$psrf[1],3)
  #Autocorrelation at lag=30
  results$LAG_30 <- round(mean(simplify2array(lapply(mm,autocorr.diag,lag=30))),3)
  
  mm_acc<-list(chains$chain1$acceptance_rates$acc.count_sigma_squared,
               chains$chain1$acceptance_rates$acc.count_sigma_squared,
               chains$chain1$acceptance_rates$acc.count_sigma_squared,
               chains$chain1$acceptance_rates$acc.count_sigma_squared)
  
  results$acceptance_rate<- mean(unlist(mm_acc))/N_iter*100
  if(true_value == T){
    
    sigma_squared<- chains$chain1$ground_truth$sigma_squared
    results$MAE=round(abs(mean(simplify2array(lapply(mm, mean))) - sigma_squared),4)
    
  }
  return(results)}
#alpha inference and diagnosics

a_summary_table<- function(chains, true_value, diag0.5,alpha,K,burnin){
  MCMC_samples <- c(chains$chain1$est_containers$a[-c(1:burnin)],
                    chains$chain2$est_containers$a[-c(1:burnin)],
                    chains$chain3$est_containers$a[-c(1:burnin)],
                    chains$chain4$est_containers$a[-c(1:burnin)])
  
  
  m<-mcmc(MCMC_samples)
  
  results = data.frame(mean_est = 0, credible_interval_95 =0)
  results$mean_est <- mean(m)
  HPD <- round(cbind(coda::HPDinterval(m)),2)
  results$credible_interval_95<- paste0("[",HPD[1],",",HPD[2],"]")
  if(true_value == T){
    results$true_value<- round(a,4)
    
  }
  return(results)}

a_diagnostic_table<- function(chains, true_value, diag0.5,K,burnin,N_iter){
  stopifnot(length(chains)==4)
  
  mm<-mcmc.list(chains_list = mcmc.list(mcmc(chains$chain1$est_containers$a[1,-c(1:burnin)]),
                                        mcmc(chains$chain1$est_containers$a[1,-c(1:burnin)]),
                                        mcmc(chains$chain1$est_containers$a[1,-c(1:burnin)]),
                                        mcmc(chains$chain1$est_containers$a[1,-c(1:burnin)])))
  
  
  results = data.frame(ESS = 0, LAG_30=0,acceptance_rate=0, Gelman_rubin=0 )
  
  
  #ESS
  results$ESS <- round(mean(simplify2array(lapply(mm, effectiveSize))),0)
  #Gelman Rubin
  results$Gelman_rubin<-  round(gelman.diag(mm)[1]$psrf[1],3)
  #Autocorrelation at lag=30
  results$LAG_30 <- round(mean(simplify2array(lapply(mm,autocorr.diag,lag=30))),3)
  
  mm_acc<-list(chains$chain1$acceptance_rates$acc.count_a,
               chains$chain2$acceptance_rates$acc.count_a,
               chains$chain3$acceptance_rates$acc.count_a,
               chains$chain4$acceptance_rates$acc.count_a)
  
  results$acceptance_rate<- mean(unlist(mm_acc))/N_iter*100
  if(true_value == T){
    a<- chains$chain1$ground_truth$a
    results$MAE=round(abs(mean(simplify2array(lapply(mm, mean))) - a),4)
    
  }
  return(results)
}

assembling_chains <- function(chains, burnin, parameter){
  
  if(parameter == 'sigma_squared'){
    assembled<- c(chains$chain1$est_containers$sigma_squared[1,-c(1:burnin)],
                  chains$chain2$est_containers$sigma_squared[1,-c(1:burnin)],
                  chains$chain3$est_containers$sigma_squared[1,-c(1:burnin)],
                  chains$chain4$est_containers$sigma_squared[1,-c(1:burnin)])
    return(assembled)
  }else if(parameter == 'a'){
    assembled<- c(chains$chain1$est_containers$a[1,-c(1:burnin)],
                  chains$chain2$est_containers$a[1,-c(1:burnin)],
                  chains$chain3$est_containers$a[1,-c(1:burnin)],
                  chains$chain4$est_containers$a[1,-c(1:burnin)])
    return(assembled)
  }else if(parameter == 'P'){
    assembled<- abind::abind(chains$chain1$est_containers$P[,,-c(1:burnin)],
                             chains$chain2$est_containers$P[,,-c(1:burnin)],
                             chains$chain3$est_containers$P[,,-c(1:burnin)],
                             chains$chain4$est_containers$P[,,-c(1:burnin)],
                             along = 3)
    return(assembled)
  }else if(parameter == 'z'){
    assembled<- cbind(chains$chain1$est_containers$z[,-c(1:burnin)],
                      chains$chain2$est_containers$z[,-c(1:burnin)],
                      chains$chain3$est_containers$z[,-c(1:burnin)],
                      chains$chain4$est_containers$z[,-c(1:burnin)])
    return(assembled)
  }else{
    print('Please provide a valide name for the parameters: a,sigma_squared,P,Z')
  }
}









U_vec_summary_table<- function(chains, true_value, diag0.5,K,burnin){
  
  MCMC_samples<- cbind(chains$chain1$est_containers$U_vec[,-c(1:burnin)],
                       chains$chain2$est_containers$U_vec[,-c(1:burnin)],
                       chains$chain3$est_containers$U_vec[,-c(1:burnin)],
                       chains$chain4$est_containers$U_vec[,-c(1:burnin)])
  
  m<- mcmc(MCMC_samples)
  
  results = data.frame(mean_est = apply(m,1,mean))
  
  HPD <- round(cbind(coda::HPDinterval(t(m))),2)
  results$credible_interval_95<- paste0("[",HPD[,1],",",HPD[,2],"]")
  if(true_value == T){
    
    results$true_value<- chains$chain1$ground_truth$U_vec
    
  }
  
  return(results)}




U_vec_diagnostic_table<- function(chains, true_value, diag0.5,K,burnin,N_iter){
  stopifnot(length(chains)==4)
  
  
  mm<-mcmc.list(chains_list = mcmc.list(mcmc(t(chains$chain1$est_containers$U_vec[,-c(1:burnin)])),
                                        mcmc(t(chains$chain2$est_containers$U_vec[,-c(1:burnin)])),
                                        mcmc(t(chains$chain3$est_containers$U_vec[,-c(1:burnin)])),
                                        mcmc(t(chains$chain4$est_containers$U_vec[,-c(1:burnin)]))))
  
  
  results = data.frame(ESS = effectiveSize(mm))
  
  
  #Gelman Rubin
  results$Gelman_rubin<-  round(gelman.diag(mm)[1]$psrf[,1],3)
  #Autocorrelation at lag=30
  U_autocorr_chain_i<-simplify2array(lapply(mm,autocorr.diag,lag=30))
  results$LAG_30 <- round(apply(U_autocorr_chain_i,2,mean),3)
  
  mm_acc<-list(chains$chain1$acceptance_rates$acc.count_U[1],
               chains$chain1$acceptance_rates$acc.count_U[1],
               chains$chain1$acceptance_rates$acc.count_U[1],
               chains$chain1$acceptance_rates$acc.count_U[1])
  
  results$acceptance_rate<- mean(unlist(mm_acc))/N_iter*100
  if(true_value == T){
    
    U_vec<- chains$chain1$ground_truth$U_vec
    results$MAE=round(abs(mean(simplify2array(lapply(mm, mean))) - U_vec),4)
    
  }
  return(results)}
