


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
est_model = 'SST'
directory = '/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/results/model_selection/K3true/raw/K2/'
est_marg_lik = function(directory, est_model, burnin, is.simulation){
  
  filenames <- list.files(pattern = paste0(est_model),path = directory)
  if(is.simulation ==T){
    log_z_t = data.frame(t= 0, expected_evidence = 0, k=0, sd = 0,k_true=0 )
  }else{
    log_z_t = data.frame(t= 0, expected_evidence = 0, k=0, sd = 0)
    z_container = data.frame(z = 0, k_est = 0)
  }
  
  for(i in 1:length(filenames)){
    
    #---------------------------------------------------------------------------
    #---------------------------------------------------------------------------
    
    item <- readRDS(file = paste0(directory,"/",filenames[i]))
    Y_ij = item$Y_ij
    N_ij = item$N_ij
    K = nrow(item$est_containers$theta)
    N_iter = ncol(item$est_containers$z)
    num_samples = burnin
    burnin= N_iter-num_samples
    
    k_est = dim(item$est_containers$theta)[1]
    z_burned = item$est_containers$z[,-(1:burnin)]
    theta_burned = item$est_containers$theta[,,-c(1:burnin)]
    
    z_pivot <- z_burned[,which.max(item$control_containers$A[-c(1:burnin)])]


    # Apply function to each chunk
    run_label_switch <- label.switching(method = "ECR" ,
                                        zpivot = z_pivot ,
                                        z = t(z_burned),
                                        K = k_est)
    if(item$t == 1){
      z_est = as.numeric(run_label_switch$clusters)
      z_container = rbind(z_container, data.frame(z = z_est, k_est = rep(k_est,nrow(Y_ij))))
    }

    permutation_z = run_label_switch$permutations$ECR


    z_permuted = matrix(NA, nrow=nrow(Y_ij), ncol = num_samples)
    for(iii in 1:ncol(z_permuted)){
      z_permuted[,iii] <- permutation_z[iii,][z_burned[,iii]]
    }

    theta_permuted = array(NA, dim=c(K,K,nrow(permutation_z)))
    for(uuu in 1:num_samples){
      # Permute the rows of matrix P
      theta_permuted[,,uuu] <- theta_burned[permutation_z[uuu,], permutation_z[uuu,],uuu]
    }



    LL <- matrix(0, ((nrow(Y_ij)*(nrow(Y_ij)-1))/2), num_samples)
    p=0
    for(ii in 1:num_samples){
      p=p+1
      z_ii<- z_permuted[,ii]
      theta_ii <- theta_permuted[,,ii]
      P_ii<- inverse_logit_f(theta_ii)
      P_ij<- calculate_victory_probabilities(vec2mat_0_P(z_ii,P = P_ii), P_ii)
      LL[,p]<- dbinom(x = Y_ij[upper.tri(Y_ij)], size = N_ij[upper.tri(N_ij)],prob = P_ij[upper.tri(P_ij)], log = T)
    }
    
    LLik_sum_item = colSums(LL)
    expected_evidence = sum(LLik_sum_item)*(1/num_samples)
    if(is.simulation ==T){
      log_z_t = rbind(log_z_t, data.frame(t = item$t , expected_evidence= expected_evidence, sd= sd(LLik_sum_item),
                                          k=k_est, k_true = item$ground_truth$K))
    }else{
      log_z_t = rbind(log_z_t, data.frame(t = item$t , expected_evidence= expected_evidence, sd= sd(LLik_sum_item),
                                          k=k_est))
    }
    
    print(paste0("estimating K=", k_est, "  ^_^ ", length(filenames) - i, ' to go *_*'))
  }
  
  log_z_t = log_z_t[-1,]
  
  
  log_z_t = log_z_t %>% arrange(t)
  log_z_t$riemann = rep(NA, nrow(log_z_t))
  for(row_i in 1:(nrow(log_z_t)-1)){
    ith_sum = (log_z_t$t[row_i+1] - log_z_t$t[row_i])*(log_z_t$expected_evidence[row_i]+log_z_t$expected_evidence[row_i+1])/2
    log_z_t$riemann[row_i] = ith_sum
  }
  
  
  marginal_likelihood = sum(log_z_t$riemann,na.rm = T)
  
  return(list(df = log_z_t, marginal_likelihood = marginal_likelihood, z_container = z_container))
}


est_model = 'SST'
directories = list()
general_df = data.frame(t =0, expected_evidence=0 , sd=0 , k=0, riemann=0)
marginal_likelohood_df = data.frame( k_est=0, marginal_likelihood=0)
z_container = data.frame(z = 0, k_est = 0)
for(k_est in 2:7){
  
  directory = paste0('/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/results/application/Citation_data/Citation_data/SST/K',k_est,'/')
  
  estimation = est_marg_lik(directory = directory,est_model = est_model,
                            burnin = 30000,is.simulation = F)
  general_df = rbind(general_df, estimation$df)
  marginal_likelohood_df = rbind(marginal_likelohood_df, data.frame(k_est=k_est, 
                                                                    marginal_likelihood=estimation$marginal_likelihood)) 
  z_container = rbind(z_container, estimation$z_container)
}


saveRDS(marginal_likelohood_df, "./results/application/Citation_data/study_on_model_selection_cite.RDS")
saveRDS(general_df, "./Desktop/Nial/POMM_pairwise/POMMs/results/application/Tennis_data/general_df_tennis.RDS")
saveRDS(z_container, "./Desktop/Nial/POMM_pairwise/POMMs/results/application/Tennis_data/z_container_citations.RDS")

z_container = readRDS('./results/application/Tennis_data/z_container_tennis.RDS')
marginal_likelihood_df1 = readRDS("./results/application/Citation_data/study_on_model_selection_cite.RDS")
marginal_likelihood_df6 = readRDS("./Desktop/Nial/POMM_pairwise/POMMs/results/model_selection/study_on_model_selection6.RDS")

marginal_likelihood_df <- marginal_likelihood_df1[-1,] %>%
  mutate(max_log_lik = ifelse(marginal_likelihood == max(marginal_likelihood), "Max Value", "")) %>%
  mutate(marginal_likelihood = round(marginal_likelihood,1))%>%
  mutate(n_clust = paste0("K = ",k_est)) %>%
  mutate(exp_lik = exp(marginal_likelihood + 1300))

marginal_likelihood_df %>%
  mutate(normalised = exp_lik/sum(exp_lik))
library(Rmpfr)

lprobabilities = marginal_likelihood_df$marginal_likelihood

x <- mpfr(lprobabilities,4)
p =exp(x)
k_given_y = p*dpois(x = 2:7, lambda = 1)/(ppois(7,1)-ppois(0,1))
a_posteriori = k_given_y/sum(k_given_y)
a_posteriori
marginal_likelihood_df = marginal_likelihood_df %>%
  mutate(prob= as.numeric(a_posteriori))%>%
  mutate(yend =rep(0,6))


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






