
#new functions
library(loo)
library(dbscan)
library(randnet)
library(fossil)
library(dplyr)
library(truncnorm)
source("/Users/lapo_santi/Desktop/Nial/project/POMMs/power-law prior/Modular_code/function_Z1.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/Functions_priorSST.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/SaraWade.R")
source("/Users/lapo_santi/Desktop/Nial/project/POMMs/power-law prior/Modular_code/functionP_POMM2.R")


#generating the overlapping sets




simulating_overlapping_POMM_powerlaw_norm = function(K,  alpha = 1, overlap=1, truncations, beta_max, diag0.5=T){
  
  #we map the proportions on a desired scale
  beta_0 <- truncations
  
  j_start = ifelse(diag0.5, yes = 1, no = 0)
  K_stop = ifelse(diag0.5, yes = K-1, no = K)
  N_levelset_i = K:1
  
  P_matrix = matrix(0, K,K)
  for( i in 1:K_stop){
    for(j in (i+j_start):K){
      level_set = abs(j-i) + abs(1-j_start)
      lb <- beta_0[level_set]
      ub <- beta_0[level_set+1]
      
      # Calculate the likelihood using a truncated distribution
      mu <- (lb + ub) / 2  # Mean of the truncated distribution
      sigma <- (ub - lb) *overlap
      
      P_matrix[i,j] = rtruncnorm(1,0.5,beta_max,mu,sigma)
    }
  }
  P_matrix[lower.tri(P_matrix)] = 1 - t(P_matrix)[lower.tri(P_matrix)]
  if(diag0.5){
    diag(P_matrix) = rep(0.5,K)
  }
  
  return(P_matrix)}

l_like_p_ij_normal_overlap = function(K, P_matrix, overlap,truncations, diag0.5 = T) {
  #we map the proportions on a desired scale
  beta_0 <- truncations
  j_start = ifelse(diag0.5, yes = 1, no = 0)
  K_stop = ifelse(diag0.5, yes = K-1, no = K)
  N_levelset_i = K:1
  
  level_sets <- diag_split_matrix(P_matrix)
  
  # Consider or exclude the main diagonal
  lowest_level_set_index <- ifelse(diag0.5, 2, 1)
  lbindex <- ifelse(diag0.5, 1, 0)
  
  log_lik_matrix = matrix(0, K,K)
  for( ii in 1:K_stop){
    for(jj in (ii+j_start):K){
      level_set = abs(jj-ii) + abs(1-j_start)
      lb <- beta_0[level_set]
      ub <- beta_0[level_set+1]
      
      # Calculate the likelihood using a truncated distribution
      mu <- (lb + ub) / 2  # Mean of the truncated distribution
      sigma <- (ub - lb) *overlap
      
      log_lik_matrix[ii,jj] = log(dtruncnorm(P_matrix[ii,jj], a = 0.5, b = beta_max, mean = mu, sd = sigma))
    }
  }
  return(sum(log_lik_matrix[upper.tri(log_lik_matrix)]))
}

truncations_try = improper_prior5(4,.8,1,T)

simulating_overlapping_POMM_powerlaw_norm(4,1,.3,truncations_try,T)

K=3
n_samples=1000
overlap = .5
beta_max = .8
alpha=1
diag0.5=T
true_alpha<-alpha


#creating a sample of P matrices
p_container = array(0, dim=c(K,K,n_samples))
for(i in 1:n_samples){
  trunc = improper_prior5(K,alpha = alpha,diag0.5 = diag0.5, beta_max = beta_max )
  p_container[,,i] = simulating_overlapping_POMM_powerlaw_norm(K, alpha=alpha, overlap  = overlap, beta_max = beta_max,truncations = trunc,diag0.5 = diag0.5)
}


# Combine the four levels into a list
level_list_p_container <- generalized_levels(p_container,K,N = n_samples,diag0.5 = diag0.5)




blue_purple <-generate_color_gradient(K)

ggplot() +
  # Add a layer for each level
  lapply(seq_along(level_list_p_container), function(i) {
    geom_density(data = data.frame(x = level_list_p_container[[i]]), aes(x = x, y = ..density.., fill = paste0("Level ", i)), alpha = .5)
  }) +
  # Set the x-axis limits
  scale_fill_manual(values=blue_purple)+
  # Set the legend title
  labs(fill = "Levels", x = "Points", y = "Density", title = paste("Density Plot of the ", K, " Level sets [alpha=",alpha,", overlap=",overlap,", diag=",0.5,"]", sep = ""))+
  theme_bw()


print(paste("Densityplot",K, "Level_sets_alpha_",alpha,"_overlap_",overlap,"_diag_",0.5, sep = "_"))

alpha_test = seq(0.1,3,0.1)

#set the containers
likelihood_est = matrix(0, nrow=n_samples, ncol=length(alpha_test))


for(j in 1: n_samples){
  for(i in 1:length(alpha_test)){
    trunc_i = improper_prior5(K,beta_max,alpha = alpha_test[i],diag0.5 = diag0.5)
    likelihood_est[j,i]= l_like_p_ij_normal_overlap(K,p_container[,,j],overlap = overlap,trunc_i,diag0.5 = diag0.5) + dlnorm_param( alpha_test[i])
  }
}

# Define the function to maximize
likelihood_function <- function(alpha) {
  likelihood_est = matrix(0, nrow=n_samples, ncol=1)
  for(j in 1: n_samples){
    trunc_i = improper_prior5(K,beta_max,alpha,diag0.5 = diag0.5)
    likelihood_est[j]= l_like_p_ij_normal_overlap(K,p_container[,,j],overlap = overlap,trunc_i,diag0.5 = diag0.5) + dlnorm_param( alpha)
  }
  return(sum(likelihood_est))
}

# Find the alpha value that maximizes the likelihood_est
max_likelihood <- optimize(likelihood_function, interval = c(0.1, 3), maximum = TRUE)
max_alpha <- max_likelihood$max

max_alpha

# Define the blue shades for the color gradient

df_diagnostic = data.frame(alpha = alpha_test, likelihood= colSums(likelihood_est)) %>% filter(likelihood != - Inf)

ggplot(df_diagnostic, aes(x = alpha, y = scale(likelihood))) +
  geom_line(colour = blue_purple[3]) +
  geom_vline(aes(xintercept = true_alpha, colour = blue_purple[1]), alpha = 0.8, show.legend = TRUE) +
  geom_vline(aes(xintercept = max_alpha, colour = "red" ), alpha = 0.8, linetype = "dashed", show.legend = TRUE) +
  theme_bw() +
  labs(x = "Alpha", y = "Likelihood", title = paste("Likelihood of alphas = [0.1,3] | ", K, " Level sets | alpha =", alpha, " | overlap =", overlap, " | diag =", 0.5, ".", sep = "")) +
  scale_color_manual(values = c("red", blue_purple[4]), labels = c("True Value", "Estimate")) +
  guides(colour = guide_legend(title = "Legend"))

print(paste("Likelihood",K, "Level_sets_alpha_",alpha,"_overlap_",overlap,"_diag_",0.5, sep = "_"))


abs(max_alpha - true_alpha)


### -----------------
# no plots
# univariate
#Simulation study without plots: recovering overlap, fixing alpha
###------------------

K=3
n_samples=1000
beta_max = .8
alpha=1
diag0.5=T
true_alpha<-alpha
overlap <- 0.35
true_overlap<-overlap


K_value_test = seq(3,11,2)
alpha_value_test = seq(0.5,3, 0.5)
results=matrix(0,length(K_value_test), length(alpha_value_test))
for(p in 1:length(K_value_test)){
  for(q in 1:length(alpha_value_test)){
    #creating a sample of P matrices
    K = K_value_test[p]
    alpha = alpha_value_test[q]
    p_container = array(0, dim=c(K,K,n_samples))
    for(i in 1:n_samples){
      trunc = improper_prior5(K,alpha = alpha,diag0.5 = diag0.5, beta_max = beta_max )
      p_container[,,i] = simulating_overlapping_POMM_powerlaw_norm(K, alpha=alpha, overlap  = overlap, beta_max = beta_max,truncations = trunc,diag0.5 = diag0.5)
    }
    # Define the function to maximize
    likelihood_function <- function(overlap) {
      likelihood_est = matrix(0, nrow=n_samples, ncol=1)
      for(j in 1: n_samples){
        trunc_i = improper_prior5(K,beta_max,alpha,diag0.5 = diag0.5)
        likelihood_est[j]= l_like_p_ij_normal_overlap(K,p_container[,,j],overlap = overlap,trunc_i,diag0.5 = diag0.5) + dlnorm_param( alpha)
      }
      return(sum(likelihood_est))
    }
    # Find the alpha value that maximizes the likelihood_est
    max_likelihood <- optimize(likelihood_function, interval = c(0.1, 1), maximum = TRUE)
    max_overlap <- max_likelihood$max
    
     results[p,q] <- abs(max_overlap - overlap)
    
  }
}

### -----------------
# no plots
# univariate
#Simulation study without plots: recovering alpha, fixing overlap
###

K=3
n_samples=1000
overlap = .5
beta_max = .8
alpha=1
diag0.5=T
true_alpha<-alpha

K_value_test = seq(3,11,2)
overlap_value_test = seq(0.1,0.9, 0.1)
results=matrix(0,length(K_value_test), length(overlap_value_test))
for(p in 1:length(K_value_test)){
  for(q in 1:length(overlap_value_test)){
    #creating a sample of P matrices
    K = K_value_test[p]
    overlap = overlap_value_test[q]
    p_container = array(0, dim=c(K,K,n_samples))
    for(i in 1:n_samples){
      trunc = improper_prior5(K,alpha = alpha,diag0.5 = diag0.5, beta_max = beta_max )
      p_container[,,i] = simulating_overlapping_POMM_powerlaw_norm(K, alpha=alpha, overlap  = overlap, beta_max = beta_max,truncations = trunc,diag0.5 = diag0.5)
    }
    # Define the function to maximize
    likelihood_function <- function(alpha) {
      likelihood_est = matrix(0, nrow=n_samples, ncol=1)
      for(j in 1: n_samples){
        trunc_i = improper_prior5(K,beta_max,alpha,diag0.5 = diag0.5)
        likelihood_est[j]= l_like_p_ij_normal_overlap(K,p_container[,,j],overlap = overlap,trunc_i,diag0.5 = diag0.5) + dlnorm_param( alpha)
      }
      return(sum(likelihood_est))
    }
    # Find the alpha value that maximizes the likelihood_est
    max_likelihood <- optimize(likelihood_function, interval = c(0.1, 3), maximum = TRUE)
    max_alpha <- max_likelihood$max
    
    results[p,q] <- abs(max_alpha - true_alpha)
    
  }
}


### -----------------
K = 3
n_samples = 10000
beta_max = 0.8
diag0.5 = TRUE
alpha=1
overlap=0.2
K_value_test = seq(3, 11, 2)
results = matrix(0, length(K_value_test), 2)

for (p in 1:length(K_value_test)) {
  K = K_value_test[p]
  p_container = array(0, dim = c(K, K, n_samples))
  
  for (i in 1:n_samples) {
    trunc = improper_prior5(K, alpha = alpha, diag0.5 = diag0.5, beta_max = beta_max)
    p_container[,,i] = simulating_overlapping_POMM_powerlaw_norm(K, alpha = alpha, overlap = overlap, beta_max = beta_max, truncations = trunc, diag0.5 = diag0.5)
  }
  
  # Define the likelihood function to maximize
  likelihood_function <- function(params) {
    alpha = params[1]
    overlap = params[2]
    
    likelihood_est = matrix(0, nrow = n_samples, ncol = 1)
    
    for (jjj in 1:n_samples) {
      trunc_i = improper_prior5(K, beta_max, alpha, diag0.5 = diag0.5)
      likelihood_est[jjj] = l_like_p_ij_normal_overlap(K, p_container[,,jjj], overlap = overlap, trunc_i, diag0.5 = diag0.5) + dlnorm_param(alpha)
    }
    
    likelihood_est[likelihood_est==-Inf] <- -2**31
   
    
    return(sum(likelihood_est))
  }
  
  # lik_explore= matrix(0,nrow = nrow(param_test), ncol=1)
  # param_test <- as.matrix(expand_grid(overlap_value_test,alpha_value_test))
  # param_test = matrix(c(0.31,0.1),nrow = 1,ncol=2)
  # for(k in 1:nrow(param_test)){
  #   lik_explore[k]<- likelihood_function(param_test[k,])
  # }
  
  
  
  # Find the alpha and overlap values that maximize the likelihood_est
  max_likelihood <- optim(c(1, 0.35), likelihood_function, method = "L-BFGS-B", lower = c(0.3, 0.1), upper = c(2.7, .9))
  max_alpha <- max_likelihood$par[1]
  max_overlap <- max_likelihood$par[2]
  
  results[p,1] <- abs(max_alpha - alpha) 
  results[p,2] <- abs(max_overlap - overlap)
}

#
#
#
K = 3
n_samples = 1000
beta_max = 0.8
diag0.5 = TRUE
alpha=
overlap=0.2
K_value_test = seq(3, 11, 2)
results = matrix(0, length(K_value_test), 2)

params <- expand.grid(alpha_value_test, overlap_value_test)

results = matrix(0, nrow=nrow(params),ncol=1)

K = 5
p_container = array(0, dim = c(K, K, n_samples))

for (i in 1:n_samples) {
  trunc = improper_prior5(K, alpha = alpha, diag0.5 = diag0.5, beta_max = beta_max)
  p_container[,,i] = simulating_overlapping_POMM_powerlaw_norm(K, alpha = alpha, overlap = overlap, beta_max = beta_max, truncations = trunc, diag0.5 = diag0.5)
}


for(rows in 1:nrow(params)){
  
  
  alpha = params[rows,1]
  overlap = params[rows,2]
  
  likelihood_est = matrix(0, nrow = n_samples, ncol = 1)
  
  for (jjj in 1:n_samples) {
    trunc_i = improper_prior5(K, beta_max, alpha, diag0.5 = diag0.5)
    likelihood_est[jjj] = l_like_p_ij_normal_overlap(K, p_container[,,jjj], overlap = overlap, trunc_i, diag0.5 = diag0.5) + dlnorm_param(alpha)
  }
  
  #likelihood_est[likelihood_est==-Inf] <- -2**31
  
  
  results[rows,1]<- sum(likelihood_est)
  
}

# Reshape results for plotting
results_df <- data.frame(alpha<- params[,1], overlap<-params[,2], likelihood = results) %>% filter(likelihood>0) %>% mutate(likelihood= likelihood**3) %>%  mutate(likelihood= scale(likelihood))


colnames(results_df) <- c("alpha", "overlap", "likelihood")
# Create the bivariate plot
ggplot(results_df, aes(x = alpha, y = overlap, fill = likelihood)) +
  geom_tile() +
  labs(x = "Alpha", y = "overlap", fill = "Likelihood") +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal()


K=3
n_samples=1000
overlap = .5
beta_max = .8
alpha=1
diag0.5=T
true_alpha<-alpha


#creating a sample of P matrices
p_container = array(0, dim=c(K,K,n_samples))
for(i in 1:n_samples){
  trunc = improper_prior5(K,alpha = alpha,diag0.5 = diag0.5, beta_max = beta_max )
  p_container[,,i] = simulating_overlapping_POMM_powerlaw_norm(K, alpha=alpha, overlap  = overlap, beta_max = beta_max,truncations = trunc,diag0.5 = diag0.5)
}


# Combine the four levels into a list
level_list_p_container <- generalized_levels(p_container,K,N = n_samples,diag0.5 = diag0.5)




blue_purple <-generate_color_gradient(K)

ggplot() +
  # Add a layer for each level
  lapply(seq_along(level_list_p_container), function(i) {
    geom_density(data = data.frame(x = level_list_p_container[[i]]), aes(x = x, y = ..density.., fill = paste0("Level ", i)), alpha = .5)
  }) +
  # Set the x-axis limits
  scale_fill_manual(values=blue_purple)+
  # Set the legend title
  labs(fill = "Levels", x = "Points", y = "Density", title = paste("Density Plot of the ", K, " Level sets [alpha=",alpha,", overlap=",overlap,", diag=",0.5,"]", sep = ""))+
  theme_bw()


print(paste("Densityplot",K, "Level_sets_alpha_",alpha,"_overlap_",overlap,"_diag_",0.5, sep = "_"))

alpha_test = seq(0.1,3,0.1)

#set the containers
likelihood_est = matrix(0, nrow=n_samples, ncol=length(alpha_test))


for(j in 1: n_samples){
  for(i in 1:length(alpha_test)){
    trunc_i = improper_prior5(K,beta_max,alpha = alpha_test[i],diag0.5 = diag0.5)
    likelihood_est[j,i]= l_like_p_ij_normal_overlap(K,p_container[,,j],overlap = overlap,trunc_i,diag0.5 = diag0.5) + dlnorm_param( alpha_test[i])
  }
}

# Define the function to maximize
likelihood_function <- function(alpha) {
  likelihood_est = matrix(0, nrow=n_samples, ncol=1)
  for(j in 1: n_samples){
    trunc_i = improper_prior5(K,beta_max,alpha,diag0.5 = diag0.5)
    likelihood_est[j]= l_like_p_ij_normal_overlap(K,p_container[,,j],overlap = overlap,trunc_i,diag0.5 = diag0.5) + dlnorm_param( alpha)
  }
  return(sum(likelihood_est))
}

# Find the alpha value that maximizes the likelihood_est
max_likelihood <- optimize(likelihood_function, interval = c(0.1, 3), maximum = TRUE)
max_alpha <- max_likelihood$max

max_alpha

# Define the blue shades for the color gradient

df_diagnostic = data.frame(alpha = alpha_test, likelihood= colSums(likelihood_est)) %>% filter(likelihood != - Inf)

ggplot(df_diagnostic, aes(x = alpha, y = scale(likelihood))) +
  geom_line(colour = blue_purple[3]) +
  geom_vline(aes(xintercept = true_alpha, colour = blue_purple[1]), alpha = 0.8, show.legend = TRUE) +
  geom_vline(aes(xintercept = max_alpha, colour = "red" ), alpha = 0.8, linetype = "dashed", show.legend = TRUE) +
  theme_bw() +
  labs(x = "Alpha", y = "Likelihood", title = paste("Likelihood of alphas = [0.1,3] | ", K, " Level sets | alpha =", alpha, " | overlap =", overlap, " | diag =", 0.5, ".", sep = "")) +
  scale_color_manual(values = c("red", blue_purple[4]), labels = c("True Value", "Estimate")) +
  guides(colour = guide_legend(title = "Legend"))

print(paste("Likelihood",K, "Level_sets_alpha_",alpha,"_overlap_",overlap,"_diag_",0.5, sep = "_"))


abs(max_alpha - true_alpha)


###








