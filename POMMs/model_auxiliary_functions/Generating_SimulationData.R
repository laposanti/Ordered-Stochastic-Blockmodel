

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
source("/Users/lapo_santi/Desktop/Nial/oldmaterial/project/simplified model/Functions_priorSST.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/Metropolis_within_Gibbs_code.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/model_auxiliary_functions/MCMC_functions.R")




generate_theta_from_theta_prior = function(K, model='WST', sigma=0){
  
  print(paste0('You are simulating theta according to the ', model, ' model prior'))
  if(model == 'SST'){
    sigma=0
    print('If model is SST, sigma should be zero')
    
    
    mu_vec_sort = seq(0.4,0.9, (0.9-0.4)/(K))
    ut <- upper.tri(matrix(0,K,K),diag = T) # get the logical matrix for upper triangular elements
    Pcombn = which(ut, arr.ind = TRUE) # get the indices of the upper triangular elements
    
    uo<- data.frame(Pcombn[sample(nrow(Pcombn)), ])# permuting the order of the rows
    n_P_entries<- nrow(uo)
    
    P_prime<-matrix(0,K,K)
    for(i_th in 1:n_P_entries){
      
      i_star<- uo$row[i_th]
      j_star<- uo$col[i_th]
      
      lower.bound = mu_vec_sort[j_star - i_star + 1] 
      upper.bound = mu_vec_sort[j_star - i_star + 2] 
      
      
      P_prime[i_star,j_star]<- runif(1, min  = lower.bound , max = upper.bound)
      
    }
    
    P_prime[lower.tri(P_prime)] = 1- t(P_prime)[lower.tri(P_prime)]
    
    
  }else if(model == 'WST'&sigma!=0){

    mu_vec_sort = seq(0.4,0.9, (0.9-0.4)/(K))
    
    ut <- upper.tri(matrix(0,K,K),diag = T) # get the logical matrix for upper triangular elements
    Pcombn = which(ut, arr.ind = TRUE) # get the indices of the upper triangular elements
    
    uo<- data.frame(Pcombn[sample(nrow(Pcombn)), ])# permuting the order of the rows
    n_P_entries<- nrow(uo)
    
    P_prime<-matrix(0,K,K)
    for(i_th in 1:n_P_entries){
      
      i_star<- uo$row[i_th]
      j_star<- uo$col[i_th]
      
      lower.bound = mu_vec_sort[j_star - i_star + 1] - sigma
      upper.bound = mu_vec_sort[j_star - i_star + 2] + sigma
      
      
      P_prime[i_star,j_star]<- runif(1, min  = lower.bound , max = upper.bound)
      
    }
    
    P_prime[lower.tri(P_prime)] = 1- t(P_prime)[lower.tri(P_prime)]
    
  }else if(model == 'Simple'){
    #upper triangular entries should not be greater than 0.5 (WST axiom) 
    #nor increasing in the columns and decreasing in the rows (SST axiom)
    
    P = matrix(NA, K, K)
    P[col(P)-row(P)==0] <-runif(K, 0,1)
    for(diag_i in 1:(K-1)){
      P[col(P)-row(P)==diag_i] <- runif( K-diag_i,0,1)
    }
    #check for P
    violating_WST_percent = 1- sum(P >=0.5 & upper.tri(P,diag = T))/sum(upper.tri(P,diag = T))
    while(violating_WST_percent < .4){
      scrambler_matrix = matrix(0,K,K)
      scrambler_matrix[upper.tri(scrambler_matrix,diag = T)] = sample(x = c(1,-1),size = (K*(K-1))/2 + K,replace = T)
      P = scrambler_matrix*P
      violating_WST_percent = 1- sum(P >=0.5 & upper.tri(P,diag = T))/sum(upper.tri(P,diag = T))
      }
    
  }
  
  if(model != 'Simple'){
  to_be_returned = list(P = P_prime, mu= mu_vec_sort)
  }else{
    to_be_returned = list(P = P_prime)
  }
  
  
  return(to_be_returned)
}
###############################################################################
# Generating data from the SST
###############################################################################

true_model = 'Simple'
saving_directory="/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/Data/Sim1_data///"
simulations = 1
k=3
for(k in 3:6){
  
  for(n_simul in 1:simulations){
    n = 80
    K=k
    seed =2021+n_simul-1
    set.seed(seed)
    if(true_model =='SST'){
      prior_SST = generate_theta_from_SST_prior(K, model = 'SST',sigma = 0)
      theta<- prior_SST$P
      mu_vec =  prior_SST$mu
    }else if(true_model == 'WST'){
      sigma_squared = 0.3
      prior_WST = generate_theta_from_SST_prior(K, model = 'WST',sigma = sigma_squared)
      theta<- prior_WST$P
      mu_vec =  prior_WST$mu
    }else if( true_model == 'Simple'){
      generate_theta_from_theta_prior(K,model = 'Simple',sigma = NA)
      P[lower.tri(P)] = 1-t(P)[lower.tri(P)]
      theta = log(P/(1-P))
    }
    
    P = inverse_logit_f(theta)
    P
    z <- sample(1:K, n,replace=T)
    z_P<- vec2mat_0_P(clust_lab = z,P = theta)
    P_nbyn<- calculate_victory_probabilities(z_mat = z_P,P = P)
    
    N_blocks = matrix(sample(x = (6:12),size = K**2,replace = T), nrow = K,ncol = K)
    N_blocks = make_symmetric(N_blocks)
    
    N_ij= calculate_victory_probabilities(z_mat = z_P,P = N_blocks)
    diag(N_ij) = 0
    #simulating Y_ij
    
    Y_ij <- matrix(0, n,n)
    for(i in 1:n){
      for(j in 1:n){
        Y_ij[i,j]<-rbinom(1,N_ij[i,j], P_nbyn[i,j])
      }
    }
    
    Y_ij[lower.tri(Y_ij)] = N_ij[lower.tri(N_ij)] - t(Y_ij)[lower.tri(Y_ij)]
    diag(Y_ij)<- 0
    
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
    
    plot_df = z_df_complete%>%
      inner_join(z_df, by = c("row" = "items")) %>%
      rename(row_z = z) %>%
      inner_join(z_df, by = c("col" = "items")) %>%
      rename(col_z = z) %>%
      mutate(row = factor(row, levels = unique(row[order(row_z, row)])),
             col = factor(col, levels = unique(col[order(col_z, col, decreasing = TRUE)])))
    
    
    
    adjacency_m<- ggplot(plot_df, aes(x = row, y = col)) +
      geom_tile(aes(fill = Y), color = "gray", show.legend = FALSE) +
      scale_fill_gradient(low = "white", high = "black") +
      geom_ysidetile(aes(color = factor(col_z)), show.legend = FALSE, width = 0.5) +
      theme_minimal() +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank())
    
    print(adjacency_m)
    
    
    if(true_model == 'SST'){
      ground_truth = list(z = z,
                          sigma_squared=NA, 
                          mu_vec_star = mu_vec,
                          K=K,theta=theta,
                          model = true_model) 
    }else if(true_model == 'WST'){
      ground_truth = list(z = z,
                          sigma_squared=sigma_squared, 
                          mu_vec_star = mu_vec,
                          K=K,theta=theta,
                          model = true_model) 
    }else if(true_model == 'Simple'){
      ground_truth = list(z = z,
                          sigma_squared=NA, 
                          mu_vec_star = NA,
                          K=K,
                          theta=theta,
                          model = true_model)
    }
    
    
    to_be_saved = list(Y_ij=Y_ij, N_ij =N_ij, ground_truth = ground_truth, 
                       data_plot = adjacency_m, 
                       seed=seed)
    
    saveRDS(to_be_saved, paste0(saving_directory, true_model, K,"_data",n_simul,".RDS"))
  }
}
