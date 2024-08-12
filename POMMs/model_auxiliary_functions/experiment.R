#experiment

source("/Users/lapo_santi/Desktop/Nial/oldmaterial/project/simplified model/Functions_priorSST.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/Metropolis_within_Gibbs_code.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/model_auxiliary_functions/MCMC_functions.R")

library(truncnorm)
library(tidyr)
library(ggplot2)
library(dplyr)
# Define a function to check the matrix
check_SST <- function(mat) {
  n <- nrow(mat)
  if (n != ncol(mat)) {
    stop("The matrix is not square.")
  }
  
  # Extract upper triangular part (excluding diagonal)
  # Extract upper triangular part (excluding diagonal)
  upper_tri <- mat
  upper_tri[lower.tri(upper_tri, diag = TRUE)] <- NA
  
  
  
  # Check if entries in each column are decreasing
  col_increasing <- sapply(2:n, function(j) all(na.omit(upper_tri[, j]) == sort(na.omit(upper_tri[, j]), decreasing = T)))
  
  # Check if entries in each row are icreasing
  row_decreasing <- sapply(1:(n-1), function(i) all(na.omit(upper_tri[i, ]) == sort(na.omit(upper_tri[i, ]), decreasing = F)))
  
  # Combine the results
  boolean=  all(col_increasing) && all(row_decreasing) && all(mat[upper.tri(mat,diag =F)] >= 0.5)
  return(boolean)
}

check_WST <- function(mat) {
  n <- nrow(mat)
  if (n != ncol(mat)) {
    stop("The matrix is not square.")
  }
  
  
  
  
  # Check if entries in each column are increasing
  greater_than_0.5 = all(mat[upper.tri(mat,diag = F)] >= 0.5)
  
  # Combine the results
  return(greater_than_0.5)
}



experiment_df = data.frame(K = 2:7, check_across_samplesSST = NA, check_across_samplesWST = NA)
for(K in 2:7){
  
  check_across_samplesWST = 0
  check_across_samplesSST = 0
  
  for(i in 1:1000){
    check_across_permutation_WST = 0
    check_across_permutation_SST = 0
    P_sample = matrix(0,K,K)
    P_sample[upper.tri(P_sample,diag =T)]<-runif((K*(K-1))/2+K)
    P_sample[lower.tri(P_sample)] = 1 - P_sample[upper.tri(P_sample)]
    
    
    permutation_matrix = permutations( n = K, r = K, v = 1:K)
    for(j in 1:nrow(permutation_matrix)){
      mat_permuted <- P_sample[permutation_matrix[j,], permutation_matrix[j,]]
      
      
      check_across_permutation_SST=  check_across_permutation_SST+ check_SST(mat_permuted)
      
      
      
      check_across_samplesWST =check_across_samplesWST+ check_WST(mat_permuted)
      
      
      check_across_samplesSST=check_across_samplesSST+check_SST(mat_permuted)
      
    }
  }
  experiment_df[which(experiment_df$K == K),2] = check_across_samplesSST
  experiment_df[which(experiment_df$K == K),3] = check_across_samplesWST
  
}
experiment_df


#experiment 2
#check how the SST condition is modified when individual covariates come into play


n=100
w = rtruncnorm(n,0,Inf)

theta_f = generate_theta_from_theta_prior(6,'SST')
theta = log(theta_f$P/(1-theta_f$P))

library(tidyr)

combn = data.frame(col = tidyr::expand_grid(w,w))
colnames(combn) = c('col_1', 'col_2')
for(i in 1:nrow(combn)){
   p = inverse_logit_f(theta + combn$col_1[i] + combn$col_2[i])
   check = check_SST(p)
   combn$checkSST[i] = check
}

combn %>%
  ggplot(aes(x = col_1, y = col_2,color= checkSST))+
  geom_point()

combn =combn%>% mutate(diff_v = abs(col_1 - col_2))%>%
  mutate(sum_v = abs(col_1 +col_2))

combn%>%
  ggplot(aes(x = sum_v, y = diff_v, color=checkSST))+
  geom_point(alpha=.5)

