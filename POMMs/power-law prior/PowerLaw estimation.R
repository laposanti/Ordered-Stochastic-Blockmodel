
K=80
alpha=1
test = simulating_POMM_powerlaw(K,alpha,beta_max)
test_m = test$matrix
try = test$truncations

#function to sample the truncations
improper_prior5 = function(K,beta_max,alpha){
  x_min =  0
  x_max = (beta_max-0.5)^(1/alpha)
  delta = (x_max - x_min)/(K-1)
  delta_sum = rep(delta,K-1)
  points = cumsum(c(x_min,delta_sum))
  truncations = points^alpha + 0.5
  return(truncations)}

#given a value of alpha, this function generates a matrix accordingly
simulating_POMM_powerlaw = function(K,  alpha = 1, beta_max){
  
  #we map the proportions on a desired scale
  beta_0 <-  improper_prior5(K,beta_max,alpha)
  
  #L_k is a list containing all the diagonals of the upper triangular matrix
  L_k = diag_split_NA(K)
  #the diagonal is made by 0.5 values
  L_k[[K]] = rep(0.5,K)
  # Compute beta values
  for(i in (K-1):1){
    n_items = length(L_k[[i]])
    #generating the beta values using the proportions as truncations
    L_k[[i]] = sample_beta_trunc(n_items,1,1,beta_0[K-i], beta_0[K-i+1])
  }
  
  #putting back everything together in the matrix
  result = diag_matrix(K)
  for(i in 0:(K-1)){
    result[which(result==-i)]<-L_k[[K-i]]
  }
  matrix = semi_symmetric(result)
  returning = list("truncations"= beta_0, "matrix"=matrix)
  return(returning)}




l_like_trunc = function(truncations, alpha, beta_max, K){
  container = matrix(0, 1, (length(truncations)-1))
  for(i in 2:(length(truncations)-1)){
    container[1,i-1] = log(pdf_power_law(truncations[i], i, K, beta_max,alpha))
  }
  
  return(rowSums(container))
}

pdf_power_law <- function(y, k, K, beta_max, alpha) {
  abs(2/(log((k-1)/(K-1))*(2*y-1))) *
    (1/sqrt(2*pi)) *
    exp(-0.5*( (log(y-0.5) - log(beta_max-0.5)/(k-1))/(log((k-1)/(K-1))) - alpha)^2)
}



l_like_p_ij = function(matrix, truncations){
  K = nrow(matrix)
  split_K=  diag_split_matrix(matrix)
  pdf_container = matrix(0,1,K-1)
  for(i in (K-1):1){
    pdf_container[1,i] = log_lik_dtruncbeta(split_K[[i]],1,1,a = truncations[-i+K],b=truncations[-i+1+K])}
  joint_p = rowSums(pdf_container)
  return(joint_p)
}



get_B_power_law = function(matrix, truncations, alpha, beta_max, K){
  return(l_like_p_ij(matrix,truncations)+ l_like_trunc(truncations, alpha, beta_max, K))
}

proposal_distribution = function(alpha0,sigma0){sample_norm_trunc(1,alpha0,sigma0,0.01,Inf)}

#sample a new value of alpha
#obtain the truncations
#generate a POMM matrix

#evaluate the ratio according to
#r = get_A () + joint_truncations



############
#metropolis hastings
############

#hyperparameters
N_iter = 30000
sigma0=1
K=80
alpha=1.5
beta_max=.8
#data
data = simulating_POMM_powerlaw(K,alpha,beta_max)
data_m = test$matrix

#containers
alpha_seq= matrix(0,1,N_iter)
A_seq = matrix(0,1,N_iter)
alpha_current <- proposal_distribution(1,sigma0)
alpha_seq[1] = alpha_current
trunc_current = improper_prior5(K,beta_max,alpha_current)
A_seq[1]= joint_truncations(data_m,alpha_current) + joint_p_ij_level_sets(data_m, trunc_current)

acc.count=0
pb=txtProgressBar(min=1,max=N_iter)
j=2
while(j<N_iter+1){
  setTxtProgressBar(pb, j)

  #proposing a new alpha
  alpha_prime <- proposal_distribution(alpha_current,sigma0)
  
  #generating a proposal matrix
  p_prime = simulating_POMM_powerlaw(K,alpha_prime,beta_max)
  
  
  #mapping it into the aux dataframe
  data_prime = df_aux_fast(synth_matches,synth_players,p_prime)
  

  #evaluating the likelihood
  r = (get_A(data_prime$n_ij, data_prime$y_ij,data_prime$p_ij)+
        l_like_p_ij(p_prime$matrix,p_prime$truncations) + 
         l_like_trunc(truncations = p_prime$truncations,alpha = alpha_prime,beta_max = beta_max,K = K)) - 
    (get_A(data_prime$n_ij, data_prime$y_ij,data_current$p_ij)+
       l_like_p_ij(p_current$matrix,p_current$truncations) + 
       l_like_trunc(truncations = p_current$truncations,alpha = alpha_current,beta_max = beta_max,K = K))
  
  u=runif(1)
  if(log(u)<r){
    acc.count=acc.count+1
    p_current = p_prime
    alpha_current = alpha_prime
    alpha_seq[j] = alpha_prime
    p_seq[,,j] = p_prime
    A_seq[j] = get_A(data_prime$n_ij, data_prime$y_ij,data_prime$p_ij)+get_B(p_prime, beta_0)
  }else{
    alpha_seq[j] = alpha_seq[j-1]
    p_seq[,,j] = p_current
    A_seq[j] = A_seq[j-1]
  }
  j=j+1
}
close(pb)
cat("accepted ", signif(100*acc.count/(N_iter),2), "%\n",sep="")

#checking the mixing of the distribution
ts.plot(A_seq)

#looking at the MSE between the MAP and the true one
mse_seq = matrix(0,N_iter,1)
for(i in 20000:N_iter){
  mse_seq[i] = sum((p_seq[,,i] - synth_p)**2)
}

plot(as.vector(A_seq), mse_seq)

burnin_p = p_seq[,,20001:30000]

########
# SImulation Study Nial
######  

mu=2
alpha = sample_norm_trunc(10000,mu,mu,0,Inf)

K=10
container_L_k = matrix(0,K+1,length(alpha))
for(i in 1:length(alpha)){
  container_L_k[,i] = improper_prior5(10,.75,alpha[i])
}

mean = rowSums(container_L_k)/length(alpha)

plot(x = 1:9, y = diff(mean[-1]), ylim=  c(0,.1))

mean(diff(mean[-1]))




#########
#
#####
K=10
alpha=.5
beta_max=.8
test = simulating_POMM_powerlaw(K,alpha = alpha, beta_max = beta_max)

ground_truth = test$truncations





x = seq(.01,10,.1)
y=matrix(0,nrow= length(x),1)
for(i in 1:length(x)){
y[i] = l_like_trunc(ground_truth,alpha = x[i],beta_max = beta_max,K = K)
}
plot(x,y, type='l')
abline(v = x[which(y==max(y))])
abline(v = alpha, col="red")
x[which(y==max(y))]



inverse_alpha_k = function(y,k,beta_max,K){
return( (log(y - 0.5) -  log( (beta_max - 0.5))) / log((k-1)/(K-1)))}

inverse = matrix(0,K, 1)
for(i in 1:(K)){
inverse[i] = inverse_alpha_k(ground_truth[i],i,beta_max = beta_max ,K = K)
}

