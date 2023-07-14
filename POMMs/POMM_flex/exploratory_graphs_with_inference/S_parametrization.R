

expected_min<- function(mu, sigma,n){
  return(mu + sigma*qnorm((1- pi/8)/(n-pi/4+1)))}
expected_max<- function(mu,sigma,n){
  return(mu + sigma*qnorm((n- pi/8)/(n-pi/4+1)))}

# Function to compute the parameter overlap given sigma
compute_S <- function(sigma, truncations,K,nsamples=1) {

  mid_points = vector()
  for(i in 1:(length(truncations))-1){
    mid_points <- append(mid_points, (truncations[i]+truncations[i+1])/2)
  }
  
  dim_level_set_i <- K:1
  expected_boundaries = matrix(0, nrow = length(mid_points), ncol=2)
  for(i in 1:length(mid_points)){
    expected_boundaries[i,1]<- expected_min(mu = mid_points[i],sigma = sigma,n = (dim_level_set_i[i])*nsamples)
    expected_boundaries[i,2]<- expected_max(mu = mid_points[i],sigma = sigma,n = (dim_level_set_i[i])*nsamples)
  }
  
  S <- 0
  possible_S <- 0
  num_expected_boundaries <- nrow(expected_boundaries)
  for (i in 1:(num_expected_boundaries - 1)) {
    for (j in (i + 1):num_expected_boundaries) {
      possible_S <- possible_S + max(expected_boundaries[j, 2],expected_boundaries[i, 2]) - min(expected_boundaries[i, 1], expected_boundaries[j, 1])
      if (expected_boundaries[i, 2] > expected_boundaries[j, 1] && expected_boundaries[j, 2] > expected_boundaries[i, 1]) {
        S <- S + (min(expected_boundaries[i, 2], expected_boundaries[j, 2]) - max(expected_boundaries[i, 1], expected_boundaries[j, 1]))
      }
    }
  }
  return(S / possible_S)
}

# Function to find the inverse of the S function
inverse_S <- function(S,truncations,K, beta_max,nsamples=1) {
  # Define the range of sigma values
  sigma_min <- 0  # Minimum possible sigma value
  sigma_max <- beta_max - 0.5  # Maximum possible sigma value
  
  # Define the desired precision
  epsilon <- 1e-6
  
  # Perform bisection method
  while (sigma_max - sigma_min > epsilon) {
    sigma_mid <- (sigma_min + sigma_max) / 2
    S_mid <- compute_S(sigma_mid, truncations,K,n_samples)
    
    if (S_mid < S) {
      sigma_min <- sigma_mid
    } else {
      sigma_max <- sigma_mid
    }
  }
  
  return(sigma_min)  # or sigma_max, they should be close at this point
}

# Example usage
desired_S <- 1.7 # Replace with the desired S value
sigma <- inverse_S(desired_S,truncations = truncations,K = K,beta_max = beta_max)
sigma



desired_S <- 1.7 # Replace with the desired overlap value
sigma <- inverse_overlap1(desired_overlap,truncations = truncations,K = K,beta_max = beta_max)
sigma

