
library(label.switching)
library(collpcm)
library(loo)

#where the data are stored
data_wd<- "/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/MCMC/results/old_studies"

setwd(data_wd)

getwd()

model<- "POMM"
filename <- list.files(pattern = paste0("True_ModelSimpleEst_model_", model))
print(filename)

#uploading the data
results_gen <- data.frame(
  MAP = 0,
  MINVI = 0,
  MISCLASSERROR = 0,
  WAIC_est = 0,
  WAIC_se = 0
)

#where to save the results
setwd("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/MCMC/results/old_studies/plots")
for(i in 1:length(filename)){
  a<- obtaining_resultsI(filename[i],model = model)
  results_gen<- rbind(results_gen,a)
}

results_gen = results_gen[-1,]
results_gen = round(results_gen,2)


setwd("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/MCMC/results/old_studies/plots")
resultsII_gen = data.frame(alphaest = 0, alpha0.05 = 0 , alpha0.95 =0 , overlapest =0, overlap0.05 = 0 , overlap0.95 =0)
resultsII_gen<- resultsII_gen[-1,]
for(i in 1:length(filename)){
  a<- obtaining_resultsII(filename[i])
  resultsII_gen<- rbind(resultsII_gen,a)
}
resultsII_gen<- round(resultsII_gen,2)
resultsII_gen = results_genII[-1,]

setwd("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/MCMC/results/old_studies/plots")
resultsIII_gen = data.frame(p_est = 0, p_est0.05 = 0 , p_est0.95 =0 , mean_MSE = 0)
for(i in 1:length(filename)){
  a<- obtaining_resultsII(filename[i])
  resultsII_gen<- rbind(resultsII_gen,a)
}
resultsII_gen<- resultsII_gen[-1,]


resultsIII_gen = data.frame(counter5perc = 0, counter1perc =0 , mean_MSE = 0)

for(i in 1:length(filename)){
  a<- obtaining_resultsIII(filename[i],"Simple")
  resultsIII_gen<- rbind(resultsIII_gen,a)
}
resultsIII_gen<- resultsIII_gen[-1,]
resultsIII_gen <- round(resultsIII_gen,2)


obj_POMM<- readRDS(paste0(data_wd,"/",filename[1]))



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



  results = data.frame(counter5perc = 0, counter1perc =0 , mean_MSE = 0)
  obj_POMM<- readRDS(paste0(data_wd,"/",filename[1]))
  print(filename)
  
  K = nrow(obj_POMM$p_true)
  M = sum(obj_POMM$Yij_matrix)
  
  entries<- which(upper.tri(obj_POMM$p_true),arr.ind = T)
 
  #pomm
  
  
  if(model=="POMM"){
    mean_MSE<- mean(abs(MSE_p_matrix(10000,obj_POMM$p_container,p_true = obj_POMM$p_true)[upper.tri(matrix(0,K,K))]))
    results$mean_MSE <-  mean_MSE
    
    p_true_POMM <- obj_POMM$p_true
    burnin_p <- obj_POMM$p_container[,,-c(1:10000)]
    plotsPOMM = list()
    for(i in 1:K) {
      for(j in 1:K) {
        y_try = data.frame(y = as.vector(burnin_p[i, j,]))
        p1 = ggplot(y_try, aes(y)) +
          geom_density(fill = "dodgerblue", alpha = 0.5) +
          scale_x_log10() +
          geom_vline(xintercept = obj_POMM$p_true[i, j], color = "red")+
          xlab("probability") +
          ylab("Density") +
          ggtitle(paste("Density plot of entry ", i, ",", j, sep = ""))
        
        plotsPOMM[[length(plotsPOMM) + 1]] <- p1
      }
    }
    p_combinedPOMM = patchwork::wrap_plots(plotsPOMM, ncol = K, nrow = K)
    # Construct the file name using paste() or paste0()
    plot_name <- paste0("P_est_",model,"K",K,"_M",M,".png")
    # Save the plot with the constructed file name
    png(plot_name,width = 800, height = 518)
    print(p_combinedPOMM)
    dev.off()
    
    counter5perc=0
    counter1perc=0
    for(i in 1:nrow(entries)){
      lb <- quantile(obj_POMM$p_container[entries[i,1],entries[i,2],-c(1:10000)],probs = 0.05)
      ub <- quantile(obj_POMM$p_container[entries[i,1],entries[i,2],-c(1:10000)], probs = 0.95)
      contained5<- obj_POMM$p_true[entries[i,1],entries[i,2]]<ub & obj_POMM$p_true[entries[i,1],entries[i,2]]>lb
      counter5perc = counter5perc + contained5
      
      lb <- quantile(obj_POMM$p_container[entries[i,1],entries[i,2],-c(1:10000)],probs = 0.01)
      ub <- quantile(obj_POMM$p_container[entries[i,1],entries[i,2],-c(1:10000)], probs = 0.99)
      contained1<- obj_POMM$p_true[entries[i,1],entries[i,2]]<ub & obj_POMM$p_true[entries[i,1],entries[i,2]]>lb
      counter1perc = counter1perc + contained1
    }
    results$counter5perc <- counter5perc/nrow(entries)
    results$counter1perc<- counter1perc/nrow(entries)
  }else{
    runPOMM<- label.switching(method = 'ECR' ,zpivot = obj_POMM$z_true,z = t(obj_POMM$z_container), K = K)
    # apply the permutations returned by typing:
    perm.POMM<- permute_array(array_samples = obj_POMM$p_container, perm_matrix = runPOMM$permutations$ECR)
    mean_MSE<- mean(abs(MSE_p_matrix(10000,perm.POMM,p_true = obj_POMM$p_true)[upper.tri(matrix(0,K,K))]))
    results$mean_MSE <-  mean_MSE
    burnin_p <- perm.POMM[,,-c(1:10000)]
    plotsPOMM = list()
    for(i in 1:K) {
      for(j in 1:K) {
        y_try = data.frame(y = as.vector(burnin_p[i, j,]))
        p1 = ggplot(y_try, aes(y)) +
          geom_density(fill = "dodgerblue", alpha = 0.5) +
          scale_x_log10() +
          geom_vline(xintercept = obj_POMM$p_true[i, j], color = "red")+
          xlab("probability") +
          ylab("Density") +
          ggtitle(paste("Density plot of entry ", i, ",", j, sep = ""))
        
        plotsPOMM[[length(plotsPOMM) + 1]] <- p1
      }
    }
    p_combinedPOMM = patchwork::wrap_plots(plotsPOMM, ncol = K, nrow = K)
    # Construct the file name using paste() or paste0()
    plot_name <- paste0("P_est_",model,"K",K,"_M",M,".png")
    # Save the plot with the constructed file name
    png(plot_name,width = 800, height = 518)
    print(p_combinedPOMM)
    dev.off()
    
    counter5perc=0
    counter1perc=0
    for(i in 1:nrow(entries)){
      lb <- quantile(perm.POMM[entries[i,1],entries[i,2],-c(1:10000)],probs = 0.05)
      ub <- quantile(perm.POMM[entries[i,1],entries[i,2],-c(1:10000)], probs = 0.95)
      contained5<- obj_POMM$p_true[entries[i,1],entries[i,2]]<ub && obj_POMM$p_true[entries[i,1],entries[i,2]]>lb
      counter5perc = counter5perc + contained5
      
      lb <- quantile(perm.POMM[entries[i,1],entries[i,2],-c(1:10000)],probs = 0.01)
      ub <- quantile(perm.POMM[entries[i,1],entries[i,2],-c(1:10000)], probs = 0.99)
      contained1<- obj_POMM$p_true[entries[i,1],entries[i,2]]<ub && obj_POMM$p_true[entries[i,1],entries[i,2]]>lb
      counter1perc = counter1perc + contained1
    }
    results$counter5perc <- counter5perc/nrow(entries)
    results$counter1perc<- counter1perc/nrow(entries)
  }
  rownames(results) <- paste0(model,"K", K, "_M", M)
  return(results)
}


obtaining_resultsII<- function(filename){
  
  obj_POMM<- readRDS(paste0(data_wd,"/",filename))
  print(filename)
  
  K = nrow(obj_POMM$p_true)
  M = sum(obj_POMM$Yij_matrix)
  print(obj_POMM$alpha)
  print(obj_POMM$overlap)
  resultsII = data.frame(alphaest = 0, alpha0.05 = 0 , alpha0.95 =0 , overlapest =0, overlap0.05 = 0 , overlap0.95 =0)
  
  
  resultsII$alphaest<- mean(obj_POMM$alpha_container[-c(1:10000)])
  resultsII$alpha0.05<- quantile(obj_POMM$alpha_container[-c(1:10000)],probs = 0.05)
  resultsII$alpha0.95<- quantile(obj_POMM$alpha_container[-c(1:10000)], probs = 0.95)
  
  resultsII$overlapest<- mean(obj_POMM$overlap_container[-c(1:10000)])
  resultsII$overlap0.05<- quantile(obj_POMM$overlap_container[-c(1:10000)],probs = 0.05)
  resultsII$overlap0.95<- quantile(obj_POMM$overlap_container[-c(1:10000)], probs = 0.95)
  print(MSE_p_matrix(10000,obj_POMM$p_container,p_true = obj_POMM$p_true))
  MSE_sum<- mean(abs(MSE_p_matrix(10000,obj_POMM$p_container,p_true = obj_POMM$p_true)[upper.tri(matrix(0,K,K))]))
  results$MSE_sum <-  MSE_sum
  
  p_true_POMM <- obj_POMM$p_true
  
  burnin_p <- obj_POMM$p_container[,,-c(1:10000)]
  plotsPOMM = list()
  for(i in 1:K) {
    for(j in 1:K) {
      y_try = data.frame(y = as.vector(burnin_p[i, j,]))
      p1 = ggplot(y_try, aes(y)) +
        geom_density(fill = "dodgerblue", alpha = 0.5) +
        scale_x_log10() +
        geom_vline(xintercept = p_true_POMM[i, j], color = "red")+
        xlab("probability") +
        ylab("Density") +
        ggtitle(paste("Density plot of entry ", i, ",", j, sep = ""))
      
      plotsPOMM[[length(plotsPOMM) + 1]] <- p1
    }
  }
  p_combinedPOMM = patchwork::wrap_plots(plotsPOMM, ncol = K, nrow = K)
  # Construct the file name using paste() or paste0()
  plot_name <- paste0("P_est_","POMM_","K",K,"_M",M,".png")
  # Save the plot with the constructed file name
  png(plot_name,width = 800, height = 518)
  print(p_combinedPOMM)
  dev.off()
  
  rownames(resultsII) <- paste0("K", K, "_M", M)
  return(resultsII)
}






obtaining_resultsI<- function(filename, model){
  
  obj_POMM<- readRDS(paste0(data_wd,"/",filename))
  print(filename)
  #Data used to generate the data -----
  K = nrow(obj_POMM$p_true)
  M = sum(obj_POMM$Yij_matrix)
  # Create a data frame to store the results
  results <- data.frame(
    MAP = 0,
    MINVI = 0,
    MISCLASSERROR = 0,
    WAIC_est = 0,
    WAIC_se = 0
  )
  A_container_POMM <- obj_POMM$A_container #likelihood across iterations
  z_container_POMM <- obj_POMM$z_container #similarity matrix
  z_truePOMM <- obj_POMM$z_true #true underlying value
  
  # Construct the file name using paste() or paste0()
  plot_name <- paste0("traceplot_",model,"K",K,"_M",M,".png")
  
  
  # Save the plot with the constructed file name
  png(plot_name,width = 800, height = 518)
  plot(ts(A_container_POMM[-c(1:N_iter*0.50)])) #checking mixing
  # Close the device to save the plot
  dev.off()
  
  
  # Construct the file name using paste() or paste0()
  plot_name <- paste0("autocorrplot_",model,"K",K,"_M",M,".png")
  # Save the plot with the constructed file name
  png(plot_name,width = 800, height = 518)
  
  acf(A_container_POMM[-c(1:N_iter*0.25)])
  dev.off()
  
  #extracting similarity matrix
  similarity_matrixPOMMM = pr_cc(z_container_POMM[,-c(1:10000)])
  
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
  
  
  #point est 1
  point_est_POMM = minVI(similarity_matrixPOMMM)$cl
  #point est 2
  z_MAP_POMM= obj_POMM$z_container[,which(obj_POMM$A_container == max(obj_POMM$A_container))[1]]
  
  #computing VI distance
  print(paste("MAP",vi.dist(z_MAP_POMM, z_truePOMM)))
  print(paste("MINVI",vi.dist(point_est_POMM, z_truePOMM)))
  
  results$MAP[1] <- vi.dist(z_MAP_POMM, z_truePOMM)
  results$MINVI<- vi.dist(point_est_POMM, z_truePOMM)
  
  #computing WAIC
  WAIC<- calculate_waic_matrix(n_ij_matrix = obj_POMM$Nij_matrix,z_container = obj_POMM$z_container,N_iter = 40000,p_container = obj_POMM$p_container,y_ij_matrix = obj_POMM$Yij_matrix )
  results$WAIC_est <- WAIC$estimates[3,1]
  results$WAIC_se <- WAIC$estimates[3,2]
  #computing MISCLASS
  N_new = 60
  z_new_init = sample(x=c(1:K),size = N_new,replace = T)
  sampled_games <- 40
  #-------we need the point estimates for p---
  #pomm
  runPOMM<- label.switching(method = 'ECR' ,zpivot = obj_POMM$z_true,z = t(obj_POMM$z_container), K = K)
  # apply the permutations returned by typing:
  perm.POMM<- permute_array(array_samples = obj_POMM$p_container, perm_matrix = runPOMM$permutations$ECR)
  #obtaining the point estimate
  p_est_POMM<- Est_p_matrix(10000,p_container = perm.POMM,p_true = obj_POMM$p_true)
  
  #------ here is the misclass
  #new games
  misss<- calculate_misclassification_rate(N_new = N_new,z_new = z_new_init,N = nrow(obj_POMM$Nij_matrix),
                                           p_true =obj_POMM$p_true,z_true =  obj_POMM$z_true,sampled_games = sampled_games,
                                           labels_available = c(1:K),P_est =p_est_POMM ,z_est = z_MAP_POMM)
  
  
  print(paste0("MISCLASSERROR = ",misss))
  results$MISCLASSERROR<- misss
  
  
  
  # Save the table to a file in LaTeX format
  rownames(results) <- paste0(model,"K", K, "_M", M)
  return(results)
}

P_summary_table<- function(MCMC_samples, true_value, diag0.5,P,K){
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
      results$true_value[i]<- P[results$entry_i[i],results$entry_j[i]]
    }
  }
  return(results)}