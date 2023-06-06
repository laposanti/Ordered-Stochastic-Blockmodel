
#generating from the model


adjl_creation <- function(M){
  # create an indicator for all diagonals in the matrix
  d <- (row(M) - col(M))
  d<- upper.tri(d,diag = T)*d + matrix(-1,nrow(d),ncol(d))*upper.tri(d,diag = T)
  
  # use split to group on these values
  
  adjl= data.frame(p_ij =NA, adjl_1=NA,adjl_2=NA )
  L_k = split(M, d)
  for(i in 1:(length(L_k)-2)){
    n = length(as.vector(L_k[[i]]))
    for(j in 1:n){
      adjl = rbind(adjl, data.frame(p_ij =as.vector(L_k[[i]])[j], 
                                    adjl_1 = as.vector(L_k[[i+1]])[j],
                                    adjl_2 = as.vector(L_k[[i+1]])[j+1]))
    }
  }
  adjl = adjl[-1,]
  return(adjl)
}

estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}
]



simulating_POMM = function(K){

#simulating from the prior
gamma_prior = matrix(runif(K**2,0,1),K,K)

#matrix needed to start building the probabilities
aux_matrix= matrix(NA,K,K)
diag(aux_matrix)<- 0.5

# create an indicator for all diagonals in the matrix
diag_indicator <- (row(aux_matrix) - col(aux_matrix))
diag_indicator<- upper.tri(diag_indicator,diag = T)*diag_indicator + matrix(-1,nrow(diag_indicator),ncol(diag_indicator))*upper.tri(diag_indicator,diag = T)

#L_k is a list containing each diagonal of the upper triangular matrix
L_k = split(aux_matrix, diag_indicator)

#scaling the gamma parameter
gamma_prior =gamma_prior/100
#creating the splits also for the gamma matrix
L_gamma = split(gamma_prior, diag_indicator)
#the variance is increasing in the probabilities
my_var = split(aux_matrix, diag_indicator)
# Compute beta values
for(i in (K-1):1){
  n_items = length(L_k[[i]])
  max_adj = matrix(NA, nrow=n_items,ncol= 1)
  for(j in 1:n_items){
    #obtaining the maximum of the 2 adjl neighborhoods
    max_adj[j] = max(L_k[[i+1]][j], L_k[[i+1]][j+1])
  }
  #max of the adjil + the gamma advantage factor
  mu = max_adj + L_gamma[[i]]
  #the variance is increasing in i
  var= -my_var[[i]]*0.001
  #setting alpha and beta s.t. E(X) = M_adjl + gamma_ij
  alpha = estBetaParams(mu, var)$alpha
  beta = estBetaParams(mu, var)$beta
  #generating the beta values
  L_k[[i]] = rbeta(n_items,alpha,beta)
}
#putting back everything together in the matrix
for(i in 1:K){
  j = (K+1)-i
  d[which(d==-i)]<-L_k[[j]]
}

return(d)}



simulating_POMM(6)

# Define a function to compute the posterior mean of the matrix
posterior_mean <- function(samples) {
  apply(samples, c(1, 2), mean)
}


# Generate posterior samples using the simulating_POMM function
n_samples <- 1000
samples <- replicate(n_samples, simulating_POMM(K))

# Compute the posterior mean
posterior_mean <- posterior_mean(samples)

# Compute the posterior standard deviation
posterior_sd <- apply(samples, c(1, 2), sd)

# Print the results
print(posterior_mean)
print(posterior_sd)

