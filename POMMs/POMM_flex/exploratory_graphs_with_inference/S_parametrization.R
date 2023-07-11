

alpha=1
beta_max = 0.8
K=5
truncations <- improper_prior5(K,beta_max,alpha,T)

means = vector()
for(i in 1:(length(truncations))-1){
  means <- append(means, (truncations[i]+truncations[i+1])/2)
}

sigma <- (beta_max - 0.5)/K
dim_level_set_i <- K:1
expected_boundaries = matrix(0, nrow = length(means), ncol=2)
for(i in 1:length(means)){
  expected_boundaries[i,1]<- expected_min(mu = means[i],sigma = sigma,n = dim_level_set_i[i])
  expected_boundaries[i,2]<- expected_max(mu = means[i],sigma = sigma,n = dim_level_set_i[i])
}

overlap_between_sets = vector()
for(i in 1:((nrow(expected_boundaries))-1)){
  overlap_i= max(0,(expected_boundaries[i,2]-expected_boundaries[i+1,1]))
  overlap_between_sets <- append(overlap_between_sets, overlap_i)
}
  
overlap_ = sum(overlap_between_sets)/length(overlap_between_sets)
overlap_

# Function to compute the parameter overlap given sigma
compute_overlap <- function(sigma, truncations,K,nsamples=1) {

  
  means = vector()
  for(i in 1:(length(truncations))-1){
    means <- append(means, (truncations[i]+truncations[i+1])/2)
  }
  

  dim_level_set_i <- K:1
  expected_boundaries = matrix(0, nrow = length(means), ncol=2)
  for(i in 1:length(means)){
    expected_boundaries[i,1]<- expected_min(mu = means[i],sigma = sigma,n = (dim_level_set_i[i])*nsamples)
    expected_boundaries[i,2]<- expected_max(mu = means[i],sigma = sigma,n = (dim_level_set_i[i])*nsamples)
  }
  
  overlap_between_sets = vector()
  for(i in 1:((nrow(expected_boundaries))-1)){
    overlap_i= max(0,(expected_boundaries[i,2]-expected_boundaries[i+1,1]))
    overlap_between_sets <- append(overlap_between_sets, overlap_i)
  }
  
  overlap_ = sum(overlap_between_sets)/length(overlap_between_sets)
  return(overlap_)
}

# Function to find the inverse of the overlap function
inverse_overlap <- function(overlap,truncations,K, beta_max) {
  # Define the range of sigma values
  sigma_min <- 0  # Minimum possible sigma value
  sigma_max <- beta_max - 0.5  # Maximum possible sigma value
  
  # Define the desired precision
  epsilon <- 1e-6
  
  # Perform bisection method
  while (sigma_max - sigma_min > epsilon) {
    sigma_mid <- (sigma_min + sigma_max) / 2
    overlap_mid <- compute_overlap(sigma_mid, truncations,K)
    
    if (overlap_mid < overlap) {
      sigma_min <- sigma_mid
    } else {
      sigma_max <- sigma_mid
    }
  }
  
  return(sigma_min)  # or sigma_max, they should be close at this point
}

# Example usage
desired_overlap <- 0.03861841  # Replace with the desired overlap value
inverse_sigma <- inverse_overlap(desired_overlap)
inverse_sigma
