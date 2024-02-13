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
  P<- inverse_logit_f(P)
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


z_plot<- function(chains, true_model, est_model, true_value, P_est, diag0.5 , K, N, z , burnin,label_switch,tap){

  Y_ij <-chains$chain1$Y_ij
  N_ij <- chains$chain1$N_ij
  z_chain = chains$chain1$est_containers$z[,-c(1:burnin)]
  
  psm<- comp.psm(t(z_chain))


  
  if(label_switch==T){
    if(true_value == F){
      #if we do not have the ground truth, use the MAP as pivotal partitioning for label switching
      
      #computing the MAP --------------------------
      #computing the likelihood of every partition
      llik<- apply(z_chain, 2, function(z_chain) 
        log_lik_f_binom(N = N_ij,  
                        Y = Y_ij,
                        z = z_chain,
                        P = P_est))
      llik = as.matrix(llik)
      z_MAP <- z_chain[,which.max(llik)]
      z_pivot = z_MAP
    }else if(true_value==T){
      #if we have the ground truth, use it as pivotal partitioning for label switching
      z_pivot = chains$chain1$ground_truth$z
    }
    
    run_label_switch <- label.switching(method = 'ECR' ,
                                        zpivot = z_pivot ,
                                        z = t(z_chain), 
                                        K = K)
    #permutations
    permutations_z<-run_label_switch$permutations$ECR
    point_est_z<- as.vector(run_label_switch$clusters)
    
    #relabeled chain
    chain_relabeled = matrix(NA, nrow=N, ncol = N_iter-burnin)
    for(i in 1:ncol(chain_relabeled)){
      chain_relabeled[,i] <- permutations_z[i,][z_chain[,i]]
    }
  }else if(label_switch==F){
    point_est_z <- minVI(psm = psm)$cl
  }
  
  # Create row and column indices
  indices <- expand.grid(row = 1:N, col = 1:N)
  
  if(true_value == T){
    z_df <- data.frame(items = 1:N, 
                       z = chains$chain1$ground_truth$z)
  }else if(true_value ==F){
    z_df <- data.frame(items = 1:N, 
                       z = as.vector(point_est_z))
  }
  
  # Convert the matrix to a data frame
  z_df_complete <- data.frame(
    row = indices$row,
    col = indices$col,
    similarity_value = NA,
    Y = NA
  )
  
  for (i in seq_len(nrow(z_df_complete))) {
    z_df_complete$Y[i] <- Y_ij[z_df_complete$col[i], z_df_complete$row[i]]
  }
  for (i in seq_len(nrow(z_df_complete))) {
    z_df_complete$similarity_value[i] <- psm[z_df_complete$col[i], z_df_complete$row[i]]
  }
  
  plot_df = z_df_complete%>%
    inner_join(z_df, by = c("row" = "items")) %>%
    rename(row_z = z) %>%
    inner_join(z_df, by = c("col" = "items")) %>%
    rename(col_z = z) %>%
    mutate(row = factor(row, levels = unique(row[order(row_z, row)])),
           col = factor(col, levels = unique(col[order(col_z, col, decreasing = TRUE)])))
  
  
  
  similarity_m <- ggplot(plot_df, aes(x = reorder(row, row_z), y = reorder(col, col_z, decreasing=T))) +
    geom_tile(aes(fill = similarity_value), color = "gray",show.legend = F) +
    scale_fill_gradient(low = "white", high = "black") +
    geom_ysidetile(aes(color=factor(col_z)), show.legend = F, width=.5)+
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  adjacency_m<- ggplot(plot_df, aes(x = row, y = col)) +
    geom_tile(aes(fill = Y), color = "gray", show.legend = FALSE) +
    scale_fill_gradient(low = "white", high = "black") +
    geom_ysidetile(aes(color = factor(col_z)), show.legend = FALSE, width = 0.5) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  
  #plotting it
  plot_name <- paste0(tap,"//adjacency_",true_model,est_model, "_K",K,"_N",N,".png")
  # Save the plot with the constructed file name
  png(plot_name,width = 800, height = 800)
  #similarity_plot(Y_ij, z, z) #checking mixing
  print(adjacency_m)
  # Close the device to save the plot
  dev.off()
  
  #plotting it
  plot_name <- paste0(tap,"//similarity_",true_model,est_model, "_K",K,"_N",N,".png")
  # Save the plot with the constructed file name
  png(plot_name,width = 800, height = 800)
  #similarity_plot(similarity_matrixPOMM, z, z) #checking mixing
  print(similarity_m)
  # Close the device to save the plot
  dev.off()
  if(label_switch == T){
    return(list(point_est= point_est_z, relabeled_chain=chain_relabeled, permutations = permutations_z ))
  }else{
    return(point_est= point_est_z)
  }
}

z_summary_table<- function(chains , true_value, z_chain_permuted, P_chain_permuted, z_est,  diag0.5 , K, burnin, label_switch = F,tap){
  
  
  Y_ij<- chains$chain1$Y_ij
  N_ij<- chains$chain1$N_ij
  n =nrow(Y_ij)
  
  if (label_switch == T){
    z_chain = z_chain_permuted
    P_chain = P_chain_permuted
  }else if(label_switch==F){
    z_chain = chains$chain1$est_containers$z
    P_chain = chains$chain1$est_containers$P
  }
  

  
  Y_ij_vector = Y_ij[upper.tri(Y_ij,diag=F)]
  N_ij_vector = N_ij[upper.tri(N_ij,diag=F)]
  N_ij_vector_clean = N_ij_vector[N_ij_vector!=0]
  Y_ij_vector_clean= Y_ij_vector[N_ij_vector!=0]
  LL <- matrix(NA,length(N_ij_vector_clean),ncol=N_iter-burnin)
  for(t in 1:ncol(LL)){

    P<- inverse_logit_f(P_chain[,,t])
    z_mat = vec2mat_0_P(z_chain[,t],P)
    P_ij<- calculate_victory_probabilities(z_mat,P)
    P_ij_vector<- P_ij[upper.tri(P_ij,diag=F)]
    
    LL[,t] = dbinom(x = Y_ij_vector_clean, 
                    size =  N_ij_vector_clean, 
                    prob = P_ij_vector[-N_ij_vector!=0], log=T)
    if(t%%3000==0){print(paste0("iteration number----->",t))}
  }
  

  
  WAIC_est_1 = waic(t(LL))
  LOO<-loo(t(LL))
  
  results = data.frame(WAIC_est =WAIC_est_1$estimates[3],
                       WAIC_SE =  WAIC_est_1$estimates[6],
                       looic =LOO$estimates[3],
                       loiic_SE = LOO$estimates[6])
  
  
  if(true_value == T){
    z_true = chains$chain1$ground_truth$z
    similarity_matrix <- pr_cc(z_chain)
    point_est_minVI = minVI(similarity_matrix)$cl
    
    results$VIdist_minVI <- vi.dist(point_est_minVI, z_true)
    results$VIdist_MAP <- vi.dist(z_est, z_true)
    return(list(table=results, VIdist_minVI = results$VIdist_minVI,LL=LL))
  }
  
  return(list(table=results, LL=LL))
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


P_summary_table <- function(chains, permutations_z, true_value, diag0.5, P, K, burnin, label_switch){

  
  P_samples <- chains$chain1$est_containers$P[,,-c(1:burnin)]
  if (label_switch == TRUE){
    P_permuted = array(NA, dim=c(K,K,nrow(permutations_z)))
    for(i in 1: nrow(permutations_z)){
      # Permute the rows of matrix P
      P_permuted[,,i] <- P_samples[permutations_z[i,], permutations_z[i,],i]
    }
    P_samples <- P_permuted
  }
  P_est<- apply(P_samples, MARGIN = c(1,2), mean)
  P_05<- apply(P_samples, MARGIN = c(1,2), function(P_samples) quantile(P_samples, probs = 0.05))
  P_95<- apply(P_samples, MARGIN = c(1,2), function(P_samples) quantile(P_samples, probs = 0.95))
  
  upper_tri_indices <- which(upper.tri(P_samples[,,1], diag = T), arr.ind = TRUE)
  
  
  results = data.frame(entry_i = upper_tri_indices[,1],
                       entry_j = upper_tri_indices[,2],
                       mean_est = P_est[upper.tri(P_est,diag = T)],
                       quantile05 = P_05[upper.tri(P_05,diag = T)],
                       quantile95 = P_95[upper.tri(P_95,diag = T)]) %>%
    mutate(mean_est = inverse_logit_f(mean_est)) %>%
    mutate(quantile05 = inverse_logit_f(quantile05))%>%
    mutate(quantile95 = inverse_logit_f(quantile95))
  
  if(true_value ==T){
    P_true = chains$chain1$ground_truth$P
    P_true = inverse_logit_f(P_true)
    results=results %>% mutate(P_true = P_true[upper.tri(P_true,diag = T)]) %>%
      mutate(P_true = inverse_logit_f(P_true))%>%
      mutate(MAE = abs(P_true - mean_est))
  }
  
  if(label_switch == T){
    return(list(table=results, P_hat = P_est, P_permuted = P_permuted))
  }else if(label_switch ==F){
    return(list(table = results, P_hat = P_est))
  }
}



#Diagnostics for P#

P_diagnostic_table<- function(chains, true_value, permutations_z, diag0.5,P,K,burnin,N_iter, label_switch){
  stopifnot(length(chains)==4)
  Y_ij <- chains$chain1$Y_ij
  
  
  P_samples_list = list(P_1 = chains$chain1$est_containers$P[,,-c(1:burnin)],
                        P_2 = chains$chain2$est_containers$P[,,-c(1:burnin)],
                        P_3 = chains$chain3$est_containers$P[,,-c(1:burnin)],
                        P_4 = chains$chain4$est_containers$P[,,-c(1:burnin)])
  
  if (label_switch == TRUE){
    for(chain in 1:4){
      P_permuted = array(NA, dim=c(K,K,nrow(permutations_z)))
      P_chain_i = P_samples_list[[chain]]
      for(i in 1: nrow(permutations_z)){
        # Permute the rows of matrix P
        P_permuted[,,i] <- P_chain_i[permutations_z[i,], permutations_z[i,],i]
      }
      P_samples_list[[chain]]<- P_permuted
    }
  }
  upper.tri.extractor = function(x){ x[upper.tri(x,diag = T)]}
  

  mm<-mcmc.list(mcmc(t(apply(P_samples_list[[1]], MARGIN = c(3), FUN = upper.tri.extractor))),
                mcmc(t(apply(P_samples_list[[2]], MARGIN = c(3), FUN = upper.tri.extractor))),
                mcmc(t(apply(P_samples_list[[3]], MARGIN = c(3), FUN = upper.tri.extractor))),
                mcmc(t(apply(P_samples_list[[4]], MARGIN = c(3), FUN = upper.tri.extractor))))
  
  upper_tri_indices <- which(upper.tri(P_samples_list$P_1[,,1], diag = T), arr.ind = TRUE)
  
  results = data.frame(entry_i = upper_tri_indices[,1],
                       entry_j = upper_tri_indices[,2],
                       ESS=effectiveSize(mm),
                       gelman_rubin = gelman.diag(mm)$psrf[,1],
                       acceptance_rate =  AcceptanceRate(mm)*100)
  
  convergence_within_chain = matrix(NA, nrow = nrow(results), ncol=4)
  #using geweke convergence diagnostics to compute the z_score for the difference in means from the first part of the chain and the last part of the chain
  for(chain in 1:4){
  #finding chains that haven't converged: ---------------
  #value 1 for z_score outside the 95percent interval,
  #value 0 for z_score within the 95percent interval
  convergence_within_chain[,chain]<- (unlist(coda::geweke.diag(mm)[[chain]]$z) < qnorm(.05))+(unlist(coda::geweke.diag(mm)[[chain]]$z) >qnorm(.95))
  }
  
  # results$n_chain_converged = rep(4,nrow(results)) - rowSums(convergence_within_chain)
  # best_chain = which.min(colSums(convergence_within_chain))
  # results$best_chain = rep(which.min(colSums(convergence_within_chain)), nrow(results))
  # 
  return(list(results = results, 
              plots_list = mm))
}

#_squared inference and diagnostics
sigma_squared_summary_table<- function(chains, true_value, diag0.5,K,burnin){
  
  sigma_chain<- c(chains$chain1$est_containers$sigma_squared[-c(1:burnin)])
  results = data.frame(est_mean = mean(sigma_chain),
                       quantile05 = quantile(sigma_chain, probs = 0.05),
                       quantile95 = quantile(sigma_chain, probs = 0.95))
  if(true_value == T){
    results = mutate(true_sigma = chains$chain1$ground_truth$sigma_squared)
  }
 plot_df = data.frame(y = sigma_chain, x = 1:length(sigma_chain))
 

  return(results)
}



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









mu_vec_summary_table<- function(chains, true_value, diag0.5,K,burnin){
  
  
  mu_vec_chain = chains$chain1$est_containers$mu_vec
  mu_vec_chain = inverse_logit_f(mu_vec_chain)
  results=data.frame(est_mean = apply(mu_vec_chain,MARGIN = 1,FUN = mean),
                     quantile05 = apply(mu_vec_chain,MARGIN = 1,
                                        FUN = function(mu_vec_chain) quantile(x = mu_vec_chain, probs=0.05)),
                     quantile95 = apply(mu_vec_chain, MARGIN = 1, 
                                        FUN = function(mu_vec_chain) quantile(x = mu_vec_chain,probs=0.95)))
  
  if(true_value == T){
    
    results$true_value<- inverse_logit_f(chains$chain1$ground_truth$mu_vec)
    
  }
  
  return(results)
}



mu_vec_diagnostic_table<- function(chains, true_value, diag0.5,K,burnin,N_iter){
  stopifnot(length(chains)==4)
  
 mu_samples_list = list(mu_1 = chains$chain1$est_containers$mu_vec[,-c(1:burnin)],
                        mu_2 = chains$chain2$est_containers$mu_vec[,-c(1:burnin)],
                        mu_3 = chains$chain3$est_containers$mu_vec[,-c(1:burnin)],
                        mu_4 = chains$chain4$est_containers$mu_vec[,-c(1:burnin)])
  
  if (label_switch == TRUE){
    for(chain in 1:4){
      mu_permuted = matrix(NA,K,nrow(permutations_z))
      mu_i = mu_samples_list[[chain]]
      for(i in 1: nrow(permutations_z)){
        # Permute the rows of matrix P
        mu_permuted[,i] <- mu_i[permutations_z[i,],i]
      }
      mu_samples_list[[chain]]<- mu_permuted
    }
  }
 
  
  
  mm<-mcmc.list(mcmc(t(mu_samples_list[[1]])),
                mcmc(t(mu_samples_list[[2]])),
                mcmc(t(mu_samples_list[[3]])),
                mcmc(t(mu_samples_list[[4]])))
  

  
  results = data.frame(mu = 1:K,
                       gelman.diag(mm)$psrf[,1],
                       ESS=effectiveSize(mm),
                       acceptance_rate =  AcceptanceRate(mm))
  
  convergence_within_chain = matrix(NA, nrow = nrow(results), ncol=4)
  #using geweke convergence diagnostics to compute the z_score for the difference in means from the first part of the chain and the last part of the chain
  for(chain in 1:4){
    #finding chains that haven't converged: ---------------
    #value 1 for z_score outside the 95percent interval,
    #value 0 for z_score within the 95percent interval
    convergence_within_chain[,chain]<- (unlist(coda::geweke.diag(mm)[[chain]]$z) < qnorm(.05))+
      (unlist(coda::geweke.diag(mm)[[chain]]$z) >qnorm(.95))
  }
  
  results$n_chain_converged = rep(4,nrow(results)) - rowSums(convergence_within_chain)
  results$best_chain = rep(which.min(colSums(convergence_within_chain)), nrow(results))
  

  return(list(results=results, plots_list=mm))
}









