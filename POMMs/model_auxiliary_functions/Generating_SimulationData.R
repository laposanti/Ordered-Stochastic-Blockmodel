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





###############################################################################
# Generating data from the SST
###############################################################################

true_model = 'Simple'
saving_directory="/Users/lapo_santi/Desktop/Nial/MCMC_results/simulation_31Jan2024/Simple_true/Simulated data/"


n = 100
M = 12
K=3
N_ij<- matrix(M,n,n)
seed =1234
set.seed(seed)
if(true_model =='SST'){
  mu_vec = seq(0.40,.85, by= (.85-.40)/(K))
  #mapping mu_vec in the R space
  mu_vec_star = log(mu_vec/(1-mu_vec))
  P_star = matrix(NA, K, K)
  for(diag_i in 0:(K-1)){
    P_star[col(P_star)-row(P_star)==diag_i] <- rep((mu_vec_star[diag_i+1]+ mu_vec_star[diag_i+2])/2, K-diag_i)
  }
  P_star[lower.tri(P_star)] = -t(P_star)[lower.tri(P_star)]
}else if( true_model == 'WST'){
  mu_vec = seq(0.40,.85, by= (.85-.40)/(K-1))
  #mapping mu_vec in the R space
  mu_vec_star = log(mu_vec/(1-mu_vec))
  sigma_i = .4
  sigma_squared = sigma_i**2
  P_star = matrix(NA, K, K)
  for(diag_i in 0:(K-1)){
    P_star[col(P_star)-row(P_star)==diag_i] <- rep((mu_vec_star[diag_i+1]+ mu_vec_star[diag_i+2])/2, K-diag_i)
  }
  P_star[upper.tri(P_star,diag = T)] =  P_star[upper.tri(P_star,diag=T)] + 
    rnorm(n = length(P_star[upper.tri(P_star,diag = T)]), mean = 0, sd = sigma_i)
  P_star[lower.tri(P_star)] = -t(P_star)[lower.tri(P_star)]
}else if( true_model == 'Simple'){
  
  P = matrix(NA, K, K)
  P[col(P)-row(P)==0] <-runif(K, .6,.85)
  for(diag_i in 1:(K-1)){
    P[col(P)-row(P)==diag_i] <- runif( K-diag_i,0.01,.9)
  }
  P[lower.tri(P)] = 1-t(P)[lower.tri(P)]
  P_star = log(P/(1-P))
}

P = inverse_logit_f(P_star)
z <- sample(1:K, n,replace=T)
z_P<- vec2mat_0_P(z,P)

P_nbyn<- calculate_victory_probabilities(z_P,P)
#simulating Y_ij

Y_ij <- matrix(0, n,n)
for(i in 1:n){
  for(j in 1:n){
    Y_ij[i,j]<- rbinom(1, N_ij[i,j],P_nbyn[i,j])
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

adjacency_m


if(true_model == 'SST'){
  ground_truth = list(z = z,
                      sigma_squared=NA, 
                      mu_vec_star = mu_vec_star,
                      K=K,P=P_star) 
}else if(true_model == 'WST'){
  ground_truth = list(z = z,
                      sigma_squared=sigma_squared, 
                      mu_vec_star = mu_vec_star,
                      K=K,P=P_star) 
}else if(true_model == 'Simple'){
  ground_truth = list(z = z,
                      sigma_squared=NA, 
                      mu_vec_star = NA,
                      K=K,
                      P=P_star) 
}


to_be_saved = list(Y_ij=Y_ij, N_ij =N_ij, ground_truth = ground_truth, 
                   data_plot = adjacency_m, 
                   seed=seed)

saveRDS(to_be_saved, paste0(saving_directory, true_model, K,".RDS"))


