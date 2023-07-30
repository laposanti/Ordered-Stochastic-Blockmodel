

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



permute_array <- function(array_samples, perm_matrix) {
  N_iter <- dim(array_samples)[3]  # Number of iterations
  K <- dim(array_samples)[1]       # Dimension of the array (K by K)
  
  permuted_array <- array(dim = c(K, K, N_iter))  # Initialize permuted array
  
  for (i in 1:N_iter) {
    perm_indices <- perm_matrix[i, ]  # Permutation indices for the current iteration
    permuted_array[, , i] <- array_samples[perm_indices, perm_indices, i]
  }
  
  return(permuted_array)
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
  K<- nrow(p_true)
  burned_p <- p_container[,,-c(1:burnin)]
  est_table = matrix(0,K,K)
  for(i in 1:K){
    for(j in 1:K){
      est_table[i,j]= mean(burned_p[i,j,])
    }
  }
  return(est_table)
}

z_plot <- function(test1, model,burn_in,directory,K,M){
  obj_POMM<- test1
  z_container_POMM <- test1$est_containers$z #container matrix
  z_truePOMM <- obj_POMM$ground_truth$z #true underlying value
  #extracting similarity matrix
  similarity_matrixPOMMM = pr_cc(z_container_POMM[,-c(1:10000)])
  setwd(directory)
  #plotting it
  plot_name <- paste0("adjacency_",model,"K",K,"_M",M,".png")
  # Save the plot with the constructed file name
  png(plot_name,width = 800, height = 800)
  similarity_plot(obj_POMM$Yij_matrix, z_truePOMM, z_truePOMM) #checking mixing
  # Close the device to save the plot
  dev.off()
  
  #plotting it
  plot_name <- paste0("similarity_",model,"K",K,"_M",M,".png")
  # Save the plot with the constructed file name
  png(plot_name,width = 800, height = 800)
  similarity_plot(similarity_matrixPOMMM, z_truePOMM, z_truePOMM) #checking mixing
  # Close the device to save the plot
  dev.off()
  
}

z_summary_table<- function(test1, model,burn_in){
  
  obj_POMM<- test1
  
  #Data used to generate the data -----
  K = nrow(obj_POMM$init$P)
  M = sum(obj_POMM$Yij_matrix)
  # Create a data frame to store the results
  results <- data.frame(
    MAP_vi_dist = 0,
    MINVI_vi_dist = 0,
    WAIC_est = 0,
    WAIC_se = 0
  )
  
  getwd()
  A_container_POMM <- obj_POMM$control_containers$A #likelihood across iterations
  z_container_POMM <- test1$est_containers$z #container matrix
  z_truePOMM <- obj_POMM$ground_truth$z #true underlying value

 
  
  #point est 1
  point_est_POMM = minVI(similarity_matrixPOMMM)$cl
  #point est 2
  z_MAP_POMM= z_container_POMM[,which(A_container_POMM == max(A_container_POMM))[1]]
  
  #computing VI distance
  print(paste("MAP",vi.dist(z_MAP_POMM, z_truePOMM)))
  print(paste("MINVI",vi.dist(point_est_POMM, z_truePOMM)))
  
  results$MAP[1] <- vi.dist(z_MAP_POMM, z_truePOMM)
  results$MINVI<- vi.dist(point_est_POMM, z_truePOMM)
  
  #computing WAIC
  WAIC<- calculate_waic_matrix(n_ij_matrix = obj_POMM$Nij_matrix,z_container = z_container_POMM,N_iter = ncol(z_container_POMM),p_container = obj_POMM$est_containers$P,y_ij_matrix = obj_POMM$Yij_matrix )
  results$WAIC_est <- WAIC$estimates[3,1]
  results$WAIC_se <- WAIC$estimates[3,2]
  # #computing MISCLASS
  # N_new = 60
  # z_new_init = sample(x=c(1:K),size = N_new,replace = T)
  # sampled_games <- 40
  # 
  #-------we need the point estimates for p---
  #pomm
  # runPOMM<- label.switching(method = 'ECR' ,zpivot = obj_POMM$z_true,z = t(obj_POMM$z_container), K = K)
  # # apply the permutations returned by typing:
  # perm.POMM<- permute_array(array_samples = obj_POMM$p_container, perm_matrix = runPOMM$permutations$ECR)
  # #obtaining the point estimate
  # p_est_POMM<- Est_p_matrix(10000,p_container = perm.POMM,p_true = obj_POMM$p_true)
  # 
  #------ here is the misclass
  #new games
  # misss<- calculate_misclassification_rate(N_new = N_new,z_new = z_new_init,N = nrow(obj_POMM$Nij_matrix),
  #                                          p_true =obj_POMM$p_true,z_true =  obj_POMM$z_true,sampled_games = sampled_games,
  #                                          labels_available = c(1:K),P_est =p_est_POMM ,z_est = z_MAP_POMM)
  # 
  # 
  # print(paste0("MISCLASSERROR = ",misss))
  # results$MISCLASSERROR<- misss
  
  
  
  # Save the table to a file in LaTeX format
  rownames(results) <- paste0(model,"K", K, "_M", M)
  return(results)
}

z_diagnostic_table<- function(chains, true_value, diag0.5,z,K,burn_in){
  stopifnot(length(chains)==4)
  
  test1<-chains[[1]]
  test2<-chains[[2]]
  test3<-chains[[3]]
  test4<-chains[[4]]
  
  
  if(true_value == F){
    results = data.frame(ESS = 0, LAG_30=0, Gelman_rubin=0, acceptance_rate=0)
    
    mm<-mcmc.list(chains_list = mcmc.list(mcmc(test1$est_containers$z[,-c(1:burn_in)]),
                                          mcmc(test2$est_containers$z[,-c(1:burn_in)]),
                                          mcmc(test3$est_containers$z[,-c(1:burn_in)]),
                                          mcmc(test4$est_containers$z[,-c(1:burn_in)])))
    mm_A<-mcmc.list(chains_list = mcmc.list(mcmc(t(test1$control_containers$A[-c(1:burn_in)])),
                                            mcmc(t(test2$control_containers$A[-c(1:burn_in)])),
                                            mcmc(t(test3$control_containers$A[-c(1:burn_in)])),
                                            mcmc(t(test4$control_containers$A[-c(1:burn_in)]))))
    #ESS
    results$ESS <- round(mean(simplify2array(lapply(mm, effectiveSize))),0)
    #Gelman Rubin
    results$Gelman_rubin<-  round(gelman.diag(mm_A)[1]$psrf[1],3)
    #Autocorrelation at lag=30
    results$LAG_30 <- round(mean(simplify2array(lapply(mm_A,autocorr.diag,lag=30))),3)
    
    mm_acc<-list(median(test1$acceptance_rates$acc.count_z),
                 median(test2$acceptance_rates$acc.count_z),
                 median(test3$acceptance_rates$acc.count_z),
                 median(test4$acceptance_rates$acc.count_z))
    
    results$acceptance_rate<- mean(unlist(mm_acc))/N_iter*100
  }else if(true_value == T){
    results = data.frame(ESS = 0, LAG_30=0, Gelman_rubin=0, acceptance_rate=0, MAP=0)
    
    mm<-mcmc.list(chains_list = mcmc.list(mcmc(test1$est_containers$z[,-c(1:burn_in)]),
                                          mcmc(test2$est_containers$z[,-c(1:burn_in)]),
                                          mcmc(test3$est_containers$z[,-c(1:burn_in)]),
                                          mcmc(test4$est_containers$z[,-c(1:burn_in)])))
    mm_A<-mcmc.list(chains_list = mcmc.list(mcmc(t(test1$control_containers$A[-c(1:burn_in)])),
                                            mcmc(t(test2$control_containers$A[-c(1:burn_in)])),
                                            mcmc(t(test3$control_containers$A[-c(1:burn_in)])),
                                            mcmc(t(test4$control_containers$A[-c(1:burn_in)]))))
    #ESS
    results$ESS <- round(mean(simplify2array(lapply(mm, effectiveSize))),0)
    #Gelman Rubin
    results$Gelman_rubin<-  round(gelman.diag(mm_A)[1]$psrf[1],3)
    #Autocorrelation at lag=30
    results$LAG_30 <- round(mean(simplify2array(lapply(mm_A,autocorr.diag,lag=30))),3)
    
    mm_acc<-list(median(test1$acceptance_rates$acc.count_z),
                 median(test2$acceptance_rates$acc.count_z),
                 median(test3$acceptance_rates$acc.count_z),
                 median(test4$acceptance_rates$acc.count_z))
    
    results$acceptance_rate<- mean(unlist(mm_acc))/N_iter*100
    A_bind = unlist(mm_A)
    brrr <-simplify2array(mm)
    z_bind = brrr[,,1]
    z_bind= cbind(z_bind,brrr[,,2],brrr[,,3],brrr[,,4])
    results$MAP<- vi.dist(z_bind[,which(A_bind==max(A_bind))],z)
  }
  return(results)}

###-----------------------------------------------------------------------------
# P summary and diagnostics
###


P_summary_table<- function(test_output, true_value, diag0.5,P,K, burn_in){
  
  MCMC_samples = test_output$est_containers$P[,,-c(1:burn_in)]
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
    results = cbind(entries_df, data.frame(mean_est = rep(0,nrow(entries_df)),
                                           credible_interval_95 =rep(0,nrow(entries_df))))
    for(i in 1:nrow(results)){
      m<-mcmc(MCMC_samples[results$entry_i[i],results$entry_j[i],])
      results$mean_est[i] <- mean(m)
      HPD <- round(cbind(coda::HPDinterval(m)),2)
      results$credible_interval_95[i]<- paste0("[",HPD[1],",",HPD[2],"]")
    }
  }else if(true_value == T){
    results = cbind(entries_df, data.frame(mean_est = rep(0,nrow(entries_df)),
                                           credible_interval_95 =rep(0,nrow(entries_df)), true_value =rep(0,nrow(entries_df))))
    for(i in 1:nrow(results)){
      m<-mcmc(MCMC_samples[results$entry_i[i],results$entry_j[i],])
      results$mean_est[i] <- round(mean(m),4)
      HPD <- round(cbind(coda::HPDinterval(m)),4)
      results$credible_interval_95[i]<- paste0("[",HPD[1],",",HPD[2],"]")
      results$true_value[i]<- round(P[results$entry_i[i],results$entry_j[i]],4)
    }
  }
  return(results)}


#Diagnostics for P#


P_diagnostic_table<- function(chains, true_value, diag0.5,P,K,burn_in){
  stopifnot(length(chains)==4)
  
  test1<-chains[[1]]
  test2<-chains[[2]]
  test3<-chains[[3]]
  test4<-chains[[4]]
  
  
  
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
                                           Gelman_rubin=rep(0,nrow(entries_df)),
                                           acceptance_rate=rep(0,nrow(entries_df))))
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
                                           Gelman_rubin=rep(0,nrow(entries_df)),
                                           acceptance_rate=rep(0,nrow(entries_df)),
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
  MCMC_samples = test_output$est_containers$S[-c(1:burn_in)]
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
    results$true_value<- round(S,4)
    
  }
  return(results)}



S_diagnostic_table<- function(chains, true_value, diag0.5,S,K,burn_in){
  stopifnot(length(chains)==4)
  
  test1<-chains[[1]]
  test2<-chains[[2]]
  test3<-chains[[3]]
  test4<-chains[[4]]
  
  mm<-mcmc.list(chains_list = mcmc.list(mcmc(test1$est_containers$S[-c(1:burn_in)]),
                                        mcmc(test2$est_containers$S[-c(1:burn_in)]),
                                        mcmc(test3$est_containers$S[-c(1:burn_in)]),
                                        mcmc(test4$est_containers$S[-c(1:burn_in)])))
  
  if(true_value == F){
    results = data.frame(ESS = 0, LAG_30=0, Gelman_rubin=0, acceptance_rate=0)
    

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
    results = data.frame(ESS = rep(0,nrow(entries_df)),
                         LAG_30=rep(0,nrow(entries_df)),
                         Gelman_rubin=rep(0,nrow(entries_df)),
                         acceptance_rate=rep(0,nrow(entries_df)),
                         MAE = rep(0,nrow(entries_df)))
    
    

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
    results$MAE=round(abs(mean(simplify2array(lapply(mm, mean))) - S),4)
    
  }
  return(results)}
#alpha inference and diagnosics
alpha_summary_table<- function(test_output, true_value, diag0.5,alpha,K,burn_in){
  MCMC_samples = test_output$est_containers$alpha[-c(1:burn_in)]
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

alpha_diagnostic_table<- function(chains, true_value, diag0.5,alpha,K,burn_in){
  stopifnot(length(chains)==4)
  
  test1<-chains[[1]]
  test2<-chains[[2]]
  test3<-chains[[3]]
  test4<-chains[[4]]
  
  mm<-mcmc.list(chains_list = mcmc.list(mcmc(test1$est_containers$alpha[-c(1:burn_in)]),
                                        mcmc(test2$est_containers$alpha[-c(1:burn_in)]),
                                        mcmc(test3$est_containers$alpha[-c(1:burn_in)]),
                                        mcmc(test4$est_containers$alpha[-c(1:burn_in)])))
  
  if(true_value == F){
    results = data.frame(ESS = 0, LAG_30=0, Gelman_rubin=0, acceptance_rate=0)
    

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
    results = data.frame(ESS = rep(0,nrow(entries_df)),
                         LAG_30=rep(0,nrow(entries_df)),
                         Gelman_rubin=rep(0,nrow(entries_df)),
                         acceptance_rate=rep(0,nrow(entries_df)),
                         MAE = rep(0,nrow(entries_df)))
    
    
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
    results$MAE=round(abs(mean(simplify2array(lapply(mm, mean))) - alpha),4)
    
  }
  return(results)}



