

source("/Users/lapo_santi/Desktop/Nial/project/POMMs/power-law prior/Modular_code/function_Z1.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/Functions_priorSST.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/SaraWade.R")
source("/Users/lapo_santi/Desktop/Nial/project/POMMs/power-law prior/Modular_code/functionP_pomm.R")

library(pander)


N_iter=40000
set.seed(34)

N=100
M= 10000
K=
alpha=1
beta_max= .85
gamma_vec=c(1/5, 2/5, 3/5)

synth = simulating_tournament_new(N = N, alpha = alpha,
                                  beta_max = beta_max,
                                  K=K, M = M,
                                  gamma_vec = gamma_vec,
                                  n_ij_max = 6,model = 'POMM',diag0.5 = T 
)

n_ij_matrix=synth$n_ij_true
y_ij_matrix=synth$y_ij_true
p_ij_true=synth$p_ij_true


z_true=synth$z_true
P_true= synth$P_matrix

similarity_plot(y_ij_matrix, z_true,z_true)




barplot(rowSums(y_ij_matrix)/rowSums(n_ij_matrix))
barplot(colSums(y_ij_matrix)/colSums(n_ij_matrix))



#upper.tri.non.zero= which(n_ij_matrix >= 0)
#activate for simple model
upper.tri.non.zero = which(n_ij_matrix > 0,arr.ind = T)
#working with just the upper triangular entries and with non-zero-values
n_ij = n_ij_matrix[upper.tri.non.zero]
y_ij = y_ij_matrix[upper.tri.non.zero]
p_ij = p_ij_true[upper.tri.non.zero]





alpha_prime = sample_norm_trunc(1, alpha_current,s =sigma_prime,a = 0.01,b = 3)
truncations_prime =  improper_prior5(K,beta_max,alpha = alpha_prime,diag0.5 = T)
p_prime = simulating_POMM_powerlaw2(K = K,alpha = alpha_prime,truncations = truncations_current,beta_max = beta_max)


A_prime<- sum(dbinom(y_ij, n_ij, p_ij_scanning, log=T))
C_prime<- l_like_p_ij1(p_scanning,truncations_scanning,diag0.5 = T) + dlnorm_param(alpha_scanning)

#------------#
# Test 1:
# fix P, simulate different alpha and truncations values, check where the likelihood is maximised
#------------#

K=4
#fix beta_max
beta_max = 0.75
#fix alpha
alpha =1
#set the number of samples
n_samples=20000

#creating a sample of P matrices
p_container = array(0, dim=c(K,K,n_samples))
for(i in 1:n_samples){
  trunc = improper_prior5(K,alpha = alpha,diag0.5 = T, beta_max = beta_max )
  p_container[,,i] = simulating_POMM_powerlaw_norm(K, alpha=alpha, beta_max = beta_max,truncations = trunc,diag0.5 = T)
}

p_container[,,4]


# Combine the four levels into a list
level_list_p_container <- generalized_levels(p_container,K=K+1,N = n_samples)


# Create a density plot of points for each level
ggplot() +
  # Add a layer for each level
  lapply(seq_along(level_list_p_container), function(i) {
    geom_density(data = data.frame(x = level_list_p_container[[i]]), aes(x = x, y = ..density.., fill = paste0("Level ", i)), alpha = alpha)
  }) +
  # Set the x-axis limits
  scale_x_continuous(limits = c(0.5, beta_max)) +
  # Set the legend title
  labs(fill = "Levels", x = "Points", y = "Density", title = "Density Plot of Points for Each Level")

heat_map_blue(p_container[,,2],title = "pipoz")






#alpha_values to be tested
alpha_test = seq(0.1,3,by=0.1)


#----new method


diag_split_matrix = function(matrix_0){
  K = nrow(matrix_0)
  diag_indicator <- (row(matrix_0)-col(matrix_0))
  diag_indicator<- upper.tri(diag_indicator,diag = T)*diag_indicator + matrix(-1,nrow(diag_indicator),ncol(diag_indicator))*upper.tri(diag_indicator,diag = T)
  L_k = split(matrix_0, diag_indicator)
  Level_sets = list()
  for(i in (K):1){
    Level_sets <- append(Level_sets, L_k[i])
  }
  return(Level_sets)
}


simulating_POMM_powerlaw_norm = function(K,  alpha = 1,truncations, beta_max, diag0.5=T){
  
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
      sigma <- (ub - lb) / 6  # Standard deviation of the truncated distribution
      
      P_matrix[i,j] = rtruncnorm(1,lb,ub,mu,sigma)
    }
  }
  P_matrix[lower.tri(P_matrix)] = 1 - t(P_matrix)[lower.tri(P_matrix)]
  if(diag0.5){
    diag(P_matrix) = rep(0.5,K)
  }
  
  return(P_matrix)}


l_like_p_ij_normal = function(K, P_matrix, truncations, diag0.5 = T) {
  K <- ifelse(diag0.5, K - 1, K)
  level_sets <- diag_split_matrix(P_matrix)
  
  # Consider or exclude the main diagonal
  lowest_level_set_index <- ifelse(diag0.5, 2, 1)
  lbindex <- ifelse(diag0.5, 1, 0)
  
  log_lik <- 0
  for (i in (1 + lbindex):(K + 1)) {
    level_set_i <- level_sets[i]
    lb <- truncations[i - lbindex]
    ub <- truncations[i + 1 - lbindex]
    
    # Calculate the likelihood using a truncated distribution
    mu <- (lb + ub) / 2  # Mean of the truncated distribution
    sigma <- (ub - lb) / 6  # Standard deviation of the truncated distribution
    
    # Calculate the log-likelihood using the truncated distribution
    log_lik <- log_lik + sum(log(dtruncnorm(unlist(level_set_i), a = lb, b = ub, mean = mu, sd = sigma)))/(mu*2)
  }
  
  return(log_lik)
}




#set the containers
likelihood_est = matrix(0, nrow=n_samples, ncol=length(alpha_test))

for(j in 1: n_samples){
  for(i in 1:length(alpha_test)){
    trunc_i = improper_prior5(K,beta_max,alpha = alpha_test[i],diag0.5 = T)
    likelihood_est[j,i]= l_like_p_ij_normal(K,p_container[,,j],trunc_i,diag0.5 = T) + dlnorm_param( alpha_test[i])
  }
}




non_inf_counts <- apply(likelihood_est, 2, function(x) sum(!is.infinite(x)))
likelihood_est[is.infinite(likelihood_est)] <- 0

plot(x= alpha_test, y = colSums(likelihood_est), "l")
abline(v = alpha)
# Count non-Inf values in each column

alpha_test[which(non_inf_counts == max(non_inf_counts))]

# Find the minimum count of non-Inf values

min_non_inf_count <- min(non_inf_counts)

non_inf_likelihood_est <- matrix(NA, nrow = min_non_inf_count, ncol = sum(non_inf_counts >= min_non_inf_count))

col_idx <- 1
for (i in 1:length(non_inf_counts)) {
  col_values <- likelihood_est[, i]
  non_inf_indices <- which(!is.infinite(col_values))
  
  if (length(non_inf_indices) >= min_non_inf_count) {
    selected_indices <- sample(non_inf_indices, min_non_inf_count)
    non_inf_likelihood_est[, col_idx] <- col_values[selected_indices]
    col_idx <- col_idx + 1
  }
}

# Count non-Inf values in each column

# Calculate the column sums of the new container
likelihood_est_sum <- colSums(non_inf_likelihood_est)


plot(x=alpha_test, y=likelihood_est_sum, "l")
abline(v=alpha)

























#----old method

#set the containers
likelihood_est = matrix(0, nrow=n_samples, ncol=length(alpha_test))

for(j in 1: n_samples){
  for(i in 1:length(alpha_test)){
    trunc_i = improper_prior5(K,beta_max,alpha = alpha_test[i],diag0.5 = F)
    likelihood_est[j,i]= l_like_p_ij(p_container[,,j],trunc_i) + dlnorm_param( alpha_test[i])
  }
}

# Count non-Inf values in each column
non_inf_counts <- apply(likelihood_est, 2, function(x) sum(!is.infinite(x)))

# Find the minimum count of non-Inf values
min_non_inf_count <- min(non_inf_counts)



non_inf_likelihood_est <- matrix(NA, nrow = min_non_inf_count, ncol = sum(non_inf_counts >= min_non_inf_count))

col_idx <- 1
for (i in 1:length(non_inf_counts)) {
  col_values <- likelihood_est[, i]
  non_inf_indices <- which(!is.infinite(col_values))
  
  if (length(non_inf_indices) >= min_non_inf_count) {
    selected_indices <- sample(non_inf_indices, min_non_inf_count)
    non_inf_likelihood_est[, col_idx] <- col_values[selected_indices]
    col_idx <- col_idx + 1
  }
}

# Count non-Inf values in each column

# Calculate the column sums of the new container
likelihood_est_sum <- colSums(likelihood_est)


plot(x=alpha_test, y=likelihood_est_sum)
abline(v=alpha)

