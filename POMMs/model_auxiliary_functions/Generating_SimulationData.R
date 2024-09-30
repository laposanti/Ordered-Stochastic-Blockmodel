library(doFuture)
library(progressr)
library(beepr)
library(foreach)
library(doParallel)
library(tidyverse)
library(EnvStats)
library(truncnorm)
library(dplyr)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggside)
library(truncnorm)
library(doRNG)
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/gitignore/oldmaterial/project/simplified model/Functions_priorSST.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/Metropolis_within_Gibbs_code.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/model_auxiliary_functions/MCMC_functions.R")




generate_theta_from_theta_prior = function(K, model){
  
  print(paste0('You are simulating theta according to the ', model, ' model prior'))
  if(model == 'SST'){
    
    mu_vec_01_sort = seq(from = 0.5,to = 0.9, length.out = K)
    mu_vec_sort = log(mu_vec_01_sort/(1-mu_vec_01_sort))
    
    ut <- upper.tri(matrix(0,K,K),diag = T) # get the logical matrix for upper triangular elements
    Pcombn = which(ut, arr.ind = TRUE) # get the indices of the upper triangular elements
    
    uo<- data.frame(Pcombn[sample(nrow(Pcombn)), ])# permuting the order of the rows
    n_P_entries<- nrow(uo)
    
    theta_prime<-matrix(0,K,K)
    for(i_th in 1:n_P_entries){
      
      i_star<- uo$row[i_th]
      j_star<- uo$col[i_th]
      
      theta_prime[i_star,j_star]<- mu_vec_sort[j_star - i_star + 1]
      
    }
    
    theta_prime[lower.tri(theta_prime)] = - t(theta_prime)[lower.tri(theta_prime)]
    
    P_prime = inverse_logit_f(theta_prime)
    
  }else if(model == 'WST'){
    
    
    P_prime = matrix(NA, K, K)
    P_prime[col(P_prime)-row(P_prime)==0] <-runif(K, 0,1)
    for(diag_i in 1:(K-1)){
      P_prime[col(P_prime)-row(P_prime)==diag_i] <- runif( K-diag_i,min = 0.5,max = 1)
    }
    
    diag(P_prime) = 0.5
    P_prime[lower.tri(P_prime)] = 1-t(P_prime)[lower.tri(P_prime)]
    
    
  }else if(model == 'Simple'){
    #upper triangular entries should not be greater than 0.5 (WST axiom) 
    #nor increasing in the columns and decreasing in the rows (SST axiom)
    
    P_prime = matrix(NA, K, K)
    P_prime[col(P_prime)-row(P_prime)==0] <-runif(K, 0,1)
    for(diag_i in 1:(K-1)){
      P_prime[col(P_prime)-row(P_prime)==diag_i] <- runif( K-diag_i,min = 0,max = 1)
    }
    diag(P_prime) = 0.5
    P_prime[lower.tri(P_prime)] = 1-t(P_prime)[lower.tri(P_prime)]
    #check for P
    violating_WST_percent = 1- sum(P_prime >=0.5 & upper.tri(P_prime,diag = T))/sum(upper.tri(P_prime,diag = T))
    
    while(violating_WST_percent < .4){
      P_prime = matrix(NA, K, K)
      P_prime[col(P_prime)-row(P_prime)==0] <-runif(K, 0,1)
      for(diag_i in 1:(K-1)){
        P_prime[col(P_prime)-row(P_prime)==diag_i] <- runif( K-diag_i,min = 0,max = 1)
      }
      violating_WST_percent = 1- sum(P_prime >=0.5 & upper.tri(P_prime,diag = T))/sum(upper.tri(P_prime,diag = T))
    }
    
  }
  
  if(model == 'SST'){
    to_be_returned = list(P = P_prime, mu= mu_vec_sort)
  }else{
    to_be_returned = list(P = P_prime)
  }
  
  
  return(to_be_returned)
}

###############################################################################
# Generating data from the SST
###############################################################################

true_model = 'SST'
saving_directory="/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/Data/Sim_2_data////"
#from 1, 'very difficult' to 5, 'very easy'
sparsity_level = c(0.3,0.5,0.7)

for(sparsity_i in sparsity_level){
  for(k in 4:7){
    for(seed_i in 1:3){
      
      n = 100
      K=k
      seed =2021+seed_i
      set.seed(seed)
      if(true_model =='SST'){
        prior_SST = generate_theta_from_theta_prior(K=k, 
                                                    model = 'SST')
        P <- prior_SST$P
        mu_vec =  prior_SST$mu
        
        theta = log(P/(1-P))
        
      }else if(true_model == 'WST'){
        
        prior_WST = generate_theta_from_theta_prior(K, 
                                                    model = 'WST')
        P <- prior_WST$P
        mu_vec =  prior_WST$mu
        
        theta = log(P/(1-P))
        
      }else if( true_model == 'Simple'){
        
        prior_Simple =  generate_theta_from_theta_prior(K,
                                                        model = 'Simple')
        P = prior_Simple$P
        P[lower.tri(P)] = 1- t(P)[lower.tri(P)]
        print(P)
        theta = log(P/(1-P))
        
      }
      
      
      
      K = nrow(P)
      z <- sample(1:K, n,replace=T)
      z_P<- vec2mat_0_P(clust_lab = z,K=K)
      P_nbyn<- calculate_victory_probabilities(z_mat = z_P,P = P)
      
      #Here we compute a KxK matrix cotaining the average number of comparisons between pairs of blocks
      # num_comparisons = recovery_capability:(recovery_capability+2)
      # N_blocks = matrix(sample(x = num_comparisons,
      #                          size = K**2,
      #                          replace = T), 
      #                   nrow = K,ncol = K)
      # N_blocks = make_symmetric(N_blocks)
      # 
      # N_ij = matrix(NA, n, n)
      # for(i in 1:n){
      #   for(j in 1:n){
      #     N_ij[i,j] = rpois(1,N_blocks[z[i], z[j]])
      #   }
      # }
      # N_ij[lower.tri(N_ij)] = t(N_ij)[lower.tri(N_ij)]
      # diag(N_ij) = 0
      
      #simulating Y_ij
      
      
      P_interact_i = rnorm(n,1/z,sd = 0.2)
      P_interact_n_by_n = P_interact_i%*%t(P_interact_i)
      inverse_logit_f(P_interact_n_by_n)
      
      n_trials = 3
      Y_ij <- matrix(0, n,n)
      N_ij = matrix(0, n,n)
      for(j in 1:n){
        for(i in 1:j){
          interact_not_interact = rbinom(1,1, 1-sparsity_i) #here you control the sparsity
          if(interact_not_interact==1){
            N_ij[i,j] <- n_trials
            N_ij[j,i] <- n_trials
            for(trial in 1:n_trials){
              win_loose <-rbinom(1,1, P_nbyn[i,j])
              if(win_loose ==1){
                Y_ij[i,j]<- Y_ij[i,j]+ 1
                
              }else{
                Y_ij[j,i]<- Y_ij[j,i]+ 1
              }
            }
          }
        }
      }
      
      
      
      
      
      indices <- expand.grid(row = 1:n, col = 1:n)
      
      
      z_df <- data.frame(items = 1:n, 
                         z = z)
      
      
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
        z_df_complete$N[i] <- N_ij[z_df_complete$col[i], z_df_complete$row[i]]
      }
      z_df_complete$y_prop = z_df_complete$Y/z_df_complete$N
      
      
      plot_df = z_df_complete%>%
        inner_join(z_df, by = c("row" = "items")) %>%
        rename(row_z = z) %>%
        inner_join(z_df, by = c("col" = "items")) %>%
        rename(col_z = z) %>%
        mutate(row = factor(row, levels = unique(row[order(row_z, row)])),
               col = factor(col, levels = unique(col[order(col_z, col, decreasing = TRUE)])))
      
      
      
      adjacency_m <- ggplot(plot_df, aes(x = row, y = col)) +
        geom_tile(aes(fill = y_prop), color = 'gray30') +
        scale_fill_gradient(low = 'white', high = 'red') +
        geom_ysidetile(aes(color = factor(col_z)), show.legend = FALSE, width = 0.5) +
        theme_minimal() +
        theme(axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              legend.position = "bottom", # Position legend below the graph
              legend.direction = "horizontal", # Extend legend horizontally,
              legend.key.width = unit(1.1, 'cm'),
              
        ) +
        guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))+
        labs(fill = 'Success %',
             caption = paste0("Data ~ ",true_model," SBM, K= ",k,", seed=",seed, ', ',sparsity_i*100, "% zeros"))# Center the title horizontally
      
      
      print(adjacency_m)
      
      
      ground_truth = list(z = z,
                          mu_vec_star = mu_vec,
                          K=K,
                          theta=theta,
                          model = true_model,
                          identification = paste0(true_model,"K",k, "-",seed,"-",sparsity_i))
      
      
      
      to_be_saved = list(Y_ij=Y_ij, N_ij =N_ij, ground_truth = ground_truth, 
                         data_plot = adjacency_m, 
                         sparsity = sparsity_i,
                         seed=seed)
      
      saveRDS(to_be_saved, paste0(saving_directory, true_model, K,"_sparsity",sparsity_i,"seed",seed,".RDS"))
    }
  }
}
