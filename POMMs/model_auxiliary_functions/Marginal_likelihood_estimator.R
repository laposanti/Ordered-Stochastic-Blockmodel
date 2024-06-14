


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



#estimating marginal likelihood
est_model = 'WST'
est_marg_lik = function(directory, est_model, is.simulation){
  
  filenames <- list.files(pattern = paste0(est_model),path = directory)
  if(is.simulation ==T){
    log_z_t = data.frame(t= 0, expected_evidence = 0, k=0, sd = 0,k_true=0 )
  }else{
    log_z_t = data.frame(t= 0, expected_evidence = 0, k=0, sd = 0)
  }
  
  n_temperatures = length(filenames)
  
  print(paste0('Computing power posterior with ',n_temperatures, 'temperatures'))
  for(i in 1:n_temperatures){
    if(i %% 10 == 0){
      print(paste0("estimating K=", k_est, " ---- ", n_temperatures - i, ' to go'))
    }
    #---------------------------------------------------------------------------
    #---------------------------------------------------------------------------
    #importing the MCMC function output
    MCMC_output <- readRDS(file = paste0(directory,"/",filenames[i]))
    Y_ij = MCMC_output$Y_ij
    N_ij = MCMC_output$N_ij
    n = nrow(Y_ij)
    #number of clusters
    k_est = dim(MCMC_output$est_containers$theta)[1]
    # parameters sampled
    z_burned = MCMC_output$est_containers$z
    theta_burned = MCMC_output$est_containers$theta
    n_samples = dim(z_burned)[2]
    
    # z_pivot <- z_burned[,which.max(MCMC_output$control_containers$A[-c(1:burnin)])]
    # 
    # 
    # # Apply function to each chunk
    # run_label_switch <- label.switching(method = "ECR" ,
    #                                     zpivot = z_pivot ,
    #                                     z = t(z_burned),
    #                                     K = k_est)
    # if(MCMC_output$t == 1){
    #   z_est = as.numeric(run_label_switch$clusters)
    #   z_container = rbind(z_container, data.frame(z = z_est, k_est = rep(k_est,nrow(Y_ij))))
    # }
    # 
    # permutation_z = run_label_switch$permutations$ECR
    # 
    # 
    # z_permuted = matrix(NA, nrow=nrow(Y_ij), ncol = num_samples)
    # for(iii in 1:ncol(z_permuted)){
    #   z_permuted[,iii] <- permutation_z[iii,][z_burned[,iii]]
    # }
    # 
    # theta_permuted = array(NA, dim=c(K,K,nrow(permutation_z)))
    # for(uuu in 1:num_samples){
    #   # Permute the rows of matrix P
    #   theta_permuted[,,uuu] <- theta_burned[permutation_z[uuu,], permutation_z[uuu,],uuu]
    # }
    
    #----------------
    # computing the likelihood of a given interaction outcome between each pair of items
    
    #likelihood container
    
    u_y = Y_ij[upper.tri(N_ij)&(N_ij!=0)]
    u_n = N_ij[upper.tri(N_ij)&(N_ij!=0)]
    
    
    LL = matrix(NA, length(u_n), n_samples)
    
    P_chain = array(apply(X = theta_burned, MARGIN = 3, FUN = inverse_logit_f),dim=c(k_est,k_est,n_samples))
    z_mat_array = array(apply(X = z_burned, MARGIN = 2, FUN = vec2mat_0_P,P=theta_burned[,,1]), dim=c(n, k_est, n_samples))
    

    for(iter in 1:n_samples){
      P_ij= calculate_victory_probabilities(z_mat_array[,,iter], P_chain[,,iter])
      LL[,iter] = dbinom(x = u_y, size = u_n, prob = P_ij[upper.tri(N_ij)&(N_ij!=0)], log=T)
    }
    if(MCMC_output$t == 1){
      WAIC_t51 = waic(t(LL))$estimates[3]
    }
    #obtaining the likelihood for a given iteration
    LLik_sum_item = colSums(LL)

    #mean likelihood across iterations (expected deviance) of a given temperature t
    expected_evidence = mean(LLik_sum_item)
    
    
    if(is.simulation ==T){
      log_z_t = rbind(log_z_t, data.frame(t = MCMC_output$t , expected_evidence= expected_evidence,
                                          sd= sd(LLik_sum_item),
                                          k=k_est, k_true = MCMC_output$ground_truth$K))
    }else{
      log_z_t = rbind(log_z_t, data.frame(t = MCMC_output$t , expected_evidence= expected_evidence, 
                                          sd= sd(LLik_sum_item),
                                          k=k_est))
    }
    
  }
  
  log_z_t = log_z_t[-1,]
  
  log_z_t = log_z_t %>% arrange(t)
  log_z_t$riemann = rep(NA, nrow(log_z_t))
  for(row_i in 1:(nrow(log_z_t)-1)){
    ith_sum = (log_z_t$t[row_i+1] - log_z_t$t[row_i])*(log_z_t$expected_evidence[row_i]+log_z_t$expected_evidence[row_i+1])/2
    log_z_t$riemann[row_i] = ith_sum
  }
  
  
  marginal_likelihood = sum(log_z_t$riemann,na.rm = T)
  
  return(list(df = log_z_t, marginal_likelihood = marginal_likelihood, 
              WAIC = WAIC_t51, z_container = z_container, Y_ij= Y_ij, N_ij = N_ij))
}

is.simulation = T
if(is.simulation ==F){
  sav_dir = 'application'
  k_true = NA
}else{
  sav_dir = 'model_choice'
  k_true = 3
}

general_df = data.frame(t =NA, expected_evidence=NA , sd=NA , k=NA, riemann=NA,k_true = NA)
marginal_likelohood_df = data.frame( k_est=NA, marginal_likelihood=NA,WAIC_est=NA,k_true = NA)
z_container = data.frame(z = 0, k_est = 0,k_true =0)

true_model = 'SST'
est_model = 'SST'



for(k_est in 2:7){
  
  # directory = paste0('/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/results/',
  #                    sav_dir,'/',true_model,'/',est_model,'/K',k_est,'/')
  # 
  directory =  paste0('/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/results/model_selection/K3true/raw/K',k_est,'/')
  estimation = est_marg_lik(directory = directory,est_model = est_model,is.simulation = T)
  general_df = rbind(general_df, estimation$df)
  marginal_likelohood_df = rbind(marginal_likelohood_df, 
                                 data.frame(k_est=k_est, 
                                            marginal_likelihood=estimation$marginal_likelihood,
                                            WAIC_est = estimation$WAIC,
                                            k_true = k_true)) 
  
  
  print(marginal_likelohood_df)
  
}

marginal_likelohood_df = marginal_likelohood_df[-1,]

# saveRDS(marginal_likelohood_df, paste0('./results/',sav_dir,'/',true_model,'/study_on_model_selection_cite_4June.RDS'))
# saveRDS(general_df, paste0('./results/',sav_dir,'/',true_model,'/general_df_4June.RDS'))
# 


marginal_likelihood_df1 = readRDS("./results/application/Tennis_data/study_on_model_selection_cite_4June.RDS")
marginal_likelihood_df6 = readRDS("./Desktop/Nial/POMM_pairwise/POMMs/results/model_selection/study_on_model_selection6.RDS")

marginal_likelihood_df1 = marginal_likelihood_df1[-1]
marginal_likelihood_df1
top_block_df_container = read.csv('./results/application/Tennis_data/processed/top_block_df.csv')
top_block_df_container = top_block_df_container[-c((381+95):nrow(top_block_df_container)),]

top_block_df_container = top_block_df_container %>% filter(items != 'juan-monaco')

a_posteriori_df = data.frame(items = top_block_df_container$items[which(top_block_df_container$K == 6)])
for(k in 3:7){
  partial_thing =  top_block_df_container[which(top_block_df_container$K == k),c('items', "p_top_block","z")]
  colnames(partial_thing)[2] <- paste0('p_top_n_', k)
  colnames(partial_thing)[3] <- paste0('z', k)
  
  partial_thing[,2]  = log(partial_thing[,2]) + 
    marginal_likelihood_df1$marginal_likelihood[k]
  a_posteriori_df = inner_join(a_posteriori_df, partial_thing, by = 'items')
  
  
}

a_posteriori_df = a_posteriori_df %>% mutate(exp_prob = rowSums(exp(a_posteriori_df[, c(2,4,6,8,10)]+300),na.rm = T))
a_posteriori_df['normalised']<- as.numeric(log(mpfr(a_posteriori_df$exp_prob,20)/exp(300)))
a_posteriori_df

# 
# marginal_likelihood_mod_df <- marginal_likelohood_df[-1,] %>%
#   mutate(max_log_lik = ifelse(marginal_likelihood == max(marginal_likelihood), "Max Value", "")) %>%
#   mutate(marginal_likelihood = round(marginal_likelihood,1))%>%
#   mutate(n_clust = paste0("K = ",k_est)) 
# 
# marginal_likelihood_mod_df = marginal_likelihood_mod_df%>%
#   filter(k_est >3 & k_est < 6)%>%
#   mutate(exp_lik = exp(marginal_likelihood + 1000))%>%
#   mutate(normalised = exp_lik/sum(exp_lik))%>%
#   mutate(log_k_givenY = marginal_likelihood + dpois(k_est, 1,T) - 
#            log(sum(marginal_likelihood * dpois(k_est, 1,T))))%>%
#     mutate(lprop_post_K = marginal_likelihood + dpois(k_est, 1,T))%>%
#     mutate(reg_2 =  exp(log_k_givenY - (max(log_k_givenY) + log(sum(exp(log_k_givenY - max(log_k_givenY))))) ))
#   
#   


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
  z_df_complete$degree_pl_row[i] <- sum(Y_ij[z_df_complete$row[i],],na.rm = T)/sum(Y_ij[,z_df_complete$row[i]],na.rm = T)
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

bt_rank_df = data.frame(Id = rownames(abilities), 
                        BTrank = abilities[,"ability"],
                        SE = abilities[,"s.e."])

a_posteriori_plot = a_posteriori_df %>% 
  inner_join(bt_rank_df, by = c('items'='Id')) %>%
  inner_join(est_df, by = c('items'='Id'))


a_posteriori_plot %>%
  ggplot(aes(x = BTrank, y = normalised, color = factor(z5)))+
  geom_point()+
  geom_errorbarh(aes(xmin =BTrank-1.96*SE , xmax =  BTrank+1.96*SE))+
  labs(title = 'The posterior allocation log-probability to the top block',
       y = 'Log p (top-block)',
       x = 'Bradley - Terry Abilities',
       color = 'Clustering K=5')+
  theme_light()


a_posteriori_df %>%
  ggplot(aes(x = abilities, y = log(degree_pl), color = factor(est_cl)))+
  geom_point()+
  labs(title = 'The posterior allocation log-probability to the top block',
       y = 'Log p (top-block)',
       x = 'Bradley - Terry Ranking')+
  theme_light()



a_posteriori_df%>%filter(est_cl == 7)















# Define a vector with your desired colors
my_colors <- c("lightblue", "skyblue", "deepskyblue", "dodgerblue", "blue", "darkblue", "midnightblue")
marginal_likelihood_df %>%
  ggplot(aes(k_est, prob)) +
  geom_point(aes(color = prob)) +
  geom_segment(aes(yend = yend, xend = k_est, color = prob)) +
  labs(fill = 'Marginal Likelihood') +
  theme_minimal() +
  labs(y = 'Posterior probability',
       x = 'Number of clusters')+
  scale_color_gradient(low = "pink1", high = "red2") +
  theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )# Set the panel background to white

z_container = z_container[-c(1:2),]
z_weight = matrix(0, 47)
for(k in 2:7){
  z_weight_df = z_container %>% filter(k_est == 2) 
  z_k = as.numeric(z_weight_df$z) * p[k-1] 
  z_weight[,1] = z_weight[,1] + z_k
}






