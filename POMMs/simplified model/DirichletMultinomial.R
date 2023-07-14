


# Define the parameters
K <- 3
n <- 10
my_gamma <- rep(1, K)

#computing the prior on z|K
dir_multinom_d<- function(N,K,n_k,my_log =T){
  my_gamma = rep(1,K)
  alpha_post <- my_gamma + n_k
  A <- sum(my_gamma)
  N_post <- sum(alpha_post)
  marginal <- gamma(A) / gamma(N_post + A) * prod(gamma(n_k + my_gamma) / gamma(my_gamma))
  if(my_log ==T){
    return(log(marginal))
  }else{return(marginal)}}

dir_multinom_d(N, K_true, n_current)
my_gamma = rep(1,K_true)
n_k =n_current
alpha_post <- my_gamma + n_current
A <- sum(alpha_post)


#Sampling from the prior on z|K
dir_multinom_r<- function(K,N,gamma=1){
  alpha=rep(1,K) # value of alpha, here chosen at random
  p=rgamma(K,alpha) # pre-simulation of the Dirichlet
  y=sample(1:K,N,prob=p/sum(p),rep=TRUE) # Multinomial
  return(y)}
