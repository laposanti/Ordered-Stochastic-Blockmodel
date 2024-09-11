
#-------------------------------------------------------------------------------
############################ MCMC functions ####################################
#-------------------------------------------------------------------------------

tuning_proposal<- function(iteration, acceptance_count, sigma, acceptanceTarget, min_sigma){
  #compute acceptance rate
  acceptanceRate <- acceptance_count / iteration
  #setting the change in the variance
  
  #passing top the log scale, to have a finer scale
  lsi= log(sigma) + (iteration**(-0.8))*(acceptanceRate - acceptanceTarget)
  
  # Update the proposal standard deviations
  sigma_updated <- min(max(min_sigma,exp(lsi)),0.3)
  return(sigma_updated)
}


label_probability <- function(distance, sigma_z) {
  prob <- dnorm(distance, mean = 0, sd = sigma_z)
  return(prob)
}

upper.tri.extractor = function(x){ x[upper.tri(x,diag = T)]}

inverse_logit_f = function(x){
  y= exp(x)/(1+exp(x))
  return(y)
}

logit_f = function(x){
  y = log(x/(1-x))
}



vec2mat_0_P <- function(clust_lab,K){
  # in: vector clust_lab of length V s.t. clust_lab[v]=h if node v is in cluster h
  # out: binary VxH matrix M s.t. M[v,h]=1{node v is in cluster h}
  n <- length(clust_lab)
  z_mat <- matrix(0,n,K)
  for (N in 1:n){
    z_mat[N,clust_lab[N]] <- 1
  }
  return(z_mat)
}


