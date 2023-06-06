

dirichlet_multinomial_sample <- function(n_samples, n, K, gamma) {
  # Sample from Dirichlet distribution with parameter vector gamma
  theta <- rgamma(K, gamma, 1)
  # Normalize to get probabilities
  theta <- theta / sum(theta)
  # Sample from multinomial distribution with probabilities theta
  rmultinom(n_samples, n, theta)
}

# Function for computing the probability mass function of Dirichlet Multinomial distribution
# Input: z - vector of length K representing the counts for each value k, gamma - parameter vector
# Output: probability mass function of Dirichlet Multinomial distribution evaluated at z
dirichlet_multinomial_pmf <- function(z, gamma) {
  B_gamma <- multibeta(gamma)
  sum(lgamma(z + gamma) - lgamma(gamma)) - B_gamma
}

# Function for computing the multivariate Beta function
# Input: alpha - parameter vector
# Output: value of the multivariate Beta function
multibeta <- function(alpha) {
  sum(lgamma(alpha)) - lgamma(sum(alpha))
}

N_iter = 1000
n=100
K=3
gamma=1
gammavec = rep(1,K)



set.seed(12)
rdirichlet_multinomial <- function(N, K, alpha_vec) {
  # Step 1: Draw theta from Dirichlet distribution
  theta <- MCMCpack::rdirichlet(1, alpha_vec)
  
  # Step 2: Draw n_k from Multinomial distribution for each k
  n <- rmultinom(N, 1, theta)
  
  #Step 3: 
  return(n)
}
set.seed(12)
rdirichlet_multinomial(10,3,rep(1,3))
rowSums(rdirichlet_multinomial(10,3,rep(1,3)))


ddirichlet_multinomial<- function(N, K, n_k, my_alpha) {
    A <- sum(my_alpha)
    constant_term <- lgamma(A) + lgamma(N+1) - lgamma(N+A)
    variable_term <- sum(lgamma(n_k + my_alpha) - lgamma(my_alpha) - lgamma(n_k+1))
    return(constant_term + variable_term)
  }
  
  
  
#setup
set.seed(14)
#number of samples
N=1000
#Number of items
n=10
#containers
K=6
#true gamma value
true_gamma=1
gamma_vec = rep(true_gamma,K)
#sampling

sample_container = matrix(0, K,N)
for(i in 1:N){
  shot1 = rdirichlet_multinomial(n,K,gamma_vec)
  sample_container[,i] = rowSums(shot1)}

z= matrix(0,n,1)
shot1
for(i in 1:ncol(shot1)){
  z[i]= which(shot1[,i] >0)
}

table(z) == rowSums(shot1)

my_log_prob(n,K,n_k, as.vector(gamma_vec))

gammas = seq(.3,5,.1)

lik_container = matrix(0, nrow=length(gammas),ncol=1)
for( i in 1:length(gammas)){
  l_lik = matrix(0,nrow=N,ncol=1)
  for(j in 1:N){
  gamma_vec = rep(gammas[i],K)
  l_lik[j] = my_log_prob(n, K, sample_container[,j], gamma_vec)
  lik_container[i] = sum(l_lik)
}
  }
plot(x=gammas,y = lik_container, type = "l")
abline(v=true_gamma)

abline(v= gammas[which(l_lik == max(l_lik))], col="red")



