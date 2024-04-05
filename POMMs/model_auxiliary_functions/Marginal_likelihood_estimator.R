




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
  }
  for(i in 1:length(filenames)){
    
    #---------------------------------------------------------------------------
    #---------------------------------------------------------------------------
    
    item <- readRDS(file = paste0(directory,"/",filenames[i]))
    Y_ij = item$Y_ij
    N_ij = item$N_ij
    K = nrow(item$est_containers$P)
    N_iter = ncol(item$est_containers$z)
    num_samples = burnin
    burnin= N_iter-num_samples
    
    k_est = dim(item$est_containers$P)[1]
    z_burned = item$est_containers$z[,-(1:burnin)]
    P_burned = item$est_containers$P[,,-c(1:burnin)]
    
    z_pivot <- z_burned[,which.max(item$control_containers$A[-c(1:burnin)])]
    
    
    # Apply function to each chunk
    run_label_switch <- label.switching(method = "ECR" ,
                                        zpivot = z_pivot ,
                                        z = t(z_burned), 
                                        K = k_est)
    
    
    permutation_z = run_label_switch$permutations$ECR
    
    z_permuted = matrix(NA, nrow=nrow(Y_ij), ncol = num_samples)
    for(iii in 1:ncol(z_permuted)){
      z_permuted[,iii] <- permutation_z[iii,][z_burned[,iii]]
    }
    
    P_permuted = array(NA, dim=c(K,K,nrow(permutation_z)))
    for(uuu in 1:num_samples){
      # Permute the rows of matrix P
      P_permuted[,,uuu] <- P_burned[permutation_z[uuu,], permutation_z[uuu,],uuu]
    }
    
    
    
    LL <- matrix(0, ((nrow(Y_ij)*(nrow(Y_ij)-1))/2), num_samples)
    p=0
    for(ii in 1:num_samples){
      p=p+1
      z_ii<- z_permuted[,ii]
      P_ii <- P_permuted[,,ii]
      theta_ii<- inverse_logit_f(P_ii)
      P_ij<- calculate_victory_probabilities(vec2mat_0_P(z_ii,P = theta_ii), theta_ii)
      LL[,p]<- dbinom(x = Y_ij[upper.tri(Y_ij)], size = N_ij[upper.tri(N_ij)],prob = P_ij[upper.tri(P_ij)], log = T)
    }
    
    LLik_sum_item = colSums(LL)
    expected_evidence = sum(LLik_sum_item)*(1/num_samples)
    if(is.simulation ==T){
    log_z_t = rbind(log_z_t, data.frame(t = item$t , expected_evidence= expected_evidence, sd= sd(LLik_sum_item),
                                        k=k, k_true = item$ground_truth$K))
    }else{
      log_z_t = rbind(log_z_t, data.frame(t = item$t , expected_evidence= expected_evidence, sd= sd(LLik_sum_item),
                                          k=k))
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
  
  return(list(df = log_z_t, marginal_likelihood = marginal_likelihood))
}


est_model = 'SST'
directories = list()
general_df = data.frame(t =0, expected_evidence=0 , sd=0 , k=0, k_true =0, riemann=0)
marginal_likelohood_df = data.frame(k_true=0, k_est=0, marginal_likelihood=0)
for(k_true in 3:5){
  for(k_est in 2:5){
    directory = paste0('/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/results/model_selection/K',k_true,'true/raw/K',k_est,'/')
    estimation = est_marg_lik(directory = directory,est_model = est_model,burnin = 30000,is.simulation = T)
    general_df = rbind(general_df, estimation$df)
    marginal_likelohood_df = rbind(marginal_likelohood_df, data.frame(k_true=k_true, k_est=k_est, marginal_likelihood=estimation$marginal_likelihood)) 
    }
}

saveRDS(marginal_likelohood_df, "./Desktop/Nial/POMM_pairwise/POMMs/results/model_selection/study_on_model_selection.RDS")
saveRDS(general_df, "./Desktop/Nial/POMM_pairwise/POMMs/results/model_selection/general_df.RDS")

marginal_likelihood_df = readRDS("./Desktop/Nial/POMM_pairwise/POMMs/results/model_selection/study_on_model_selection.RDS")
marginal_likelihood_df = marginal_likelihood_df[-1,] %>% 
 mutate(marginal_likelihood= round(marginal_likelihood,0))


marginal_likelihood_df <- marginal_likelihood_df %>%
  mutate(k_true = paste0("K true = ", k_true)) %>%
  group_by(k_true) %>%
  mutate(max_log_lik = ifelse(marginal_likelihood == max(marginal_likelihood), "Max Value", "")) 



marginal_likelihood_df %>% ggplot(aes(k_est,marginal_likelihood))+
  geom_point(aes(color = max_log_lik),show.legend = F)+
  facet_wrap(~k_true)+
  labs(title = 'Power Posterior estimates of marginal likelihood',
       subtitle = 'The highest marginal is highlighted in green',
       x = 'Estimated K', y='Marginal Likelihood',
       caption = 'Data simulated from the SST model \n Each model has been fitted with a fixed value for K')+
  theme_bw()














