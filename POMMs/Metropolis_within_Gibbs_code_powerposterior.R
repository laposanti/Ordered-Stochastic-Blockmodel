
#making inference on dd
#z: the mixture membership vector
#a: the max attainable probability
#sigma_squared: the transitivity of preferences
#mu: the means of the level sets

#this code is modular
#estimation_control is a 4x1 vector, denoting which parameters must be estimated. 
#estimation_control[p]=1, means that the parameter p will be estimated. 
#estimation_control[p]=0, means that the parameter p will be fixed at the true value 

#ground_truth is a 4x1 vector, with the true value for the parameters
#if ground_truth is missing, then we are not in a simulation study

#hyper_params should contain a list of hyperparameters
#K: the number of blocks
#beta_max: the highest attainable probability in P
#gamma_vec: the hyperprior on the block sizes
#diag0.5: whether the main diagonal is set to 0.5 (diag0.5=T) or not (diag0.5=F)

adaptive_MCMC_orderstats_powerposterior <- function(Y_ij, N_ij , estimation_control,
                                                    ground_truth=NA,n, N_iter, burnin, data_description,
                                                    K_est, seed, model, saving_directory, 
                                                    custom_init=NA,power_posterior_apprach = T, thin=1,
                                                    diag0.5){
  #setting for each chain a different seed
  #if the given seed is 20, the chains' seeds will be 21 for chain 1, 22 for chain 2 and so on...
  
  
  
  
  
  
  
  
  
  # where_to_save =list()
  # for(i in 1:length(K_est)){
  #   if(power_posterior_apprach==T){
  #     where_to_save[[i]] =  file.path(saving_directory, paste0("/MCMC_output/powerposterior/Data_",data_description,"/Est_",model,"/K", K_est[[i]],"/"))
  #     dir.create(where_to_save[[i]], showWarnings = F, recursive = T)
  #   }
  #   
  #   else{
  #     where_to_save[[i]] =  NA
  #   }
  #   
  # }
  # 
  # 
  
  n_chains = length(K_est)
  if(detectCores() < n_chains){
    cat('Warning: Number of cores exceeds the ')
  }
  
  
  variables_to_add = c('Y_ij', 'N_ij' , 'estimation_control', 
                       'ground_truth','n', 'N_iter','n_chains', 
                       'optimal_acceptance_rate_theta', 'optimal_acceptance_rate_mu', 'K_est','thin', 'burnin', 'seed','model','data_description',
                       'power_posterior_apprach' ,'true_model', 'custom_init','p','n_temperatures','where_to_save','diag0.5')
  
  registerDoFuture()
  reprex <- local({
    handlers(global = TRUE)
    p <- progressor(steps = n_chains * N_iter/5000)
    plan(multisession, workers= n_chains, gc = TRUE)
    y <- foreach(chain = 1:n_chains,.options.future = list(globals = structure(F, 
                                                                               add=variables_to_add),options(future.globals.onReference = "error"),
                                                           seed=TRUE)) %dofuture%{ 
                                                             
                                                             library(doFuture,quietly = T)
                                                             library(progressr,quietly = T)
                                                             library(foreach,quietly = T)
                                                             library(doParallel,quietly = T)
                                                             library(EnvStats,quietly = T)
                                                             library(truncnorm,quietly = T)
                                                             library(dplyr,quietly = T)
                                                             library(parallel,quietly = T)
                                                             library(truncnorm,quietly = T)
                                                             library(doRNG,quietly = T)
                                                             
                                                             source("./model_auxiliary_functions/Functions_priorSST.R")
                                                             source("./model_auxiliary_functions/MCMC_functions.R")
                                                             
                                                             ############################################
                                                             #Define the prior model of choice for theta
                                                             ############################################
                                                             
                                                             # Prior and proposal on theta for the SST model
                                                             #--------------------------------
                                                             if(model =='SST'){
                                                               #if we impose that the main diagonal should be fixed to 0.5
                                                               if(diag0.5 == T){
                                                                 #Prior distribution on theta
                                                                 theta_prior_probability <-function(theta,K){
                                                                   fac <- sum(lfactorial((K-1):1))
                                                                   joint_density<- log(dtruncnorm(theta[1,2:K],a = 0,mean = 0,sd = 1)) - 
                                                                     log(0.5)*((K-1):1)
                                                                   
                                                                   log_p_sum = joint_density + fac 
                                                                   return(log_p_sum)
                                                                 }
                                                               }else if(diag0.5 == F){
                                                                 
                                                                 
                                                                 theta_prior_probability <-function(theta,K){
                                                                   fac <- sum(lfactorial((K-1):1))
                                                                   P_diag = inverse_logit_f(diag(theta))
                                                                   
                                                                   
                                                                   P_diag_density = dbeta(P_diag,shape1 = 2,shape2 = 2)*
                                                                     (exp(P_diag)/((1+exp(P_diag))**2))
                                                                   
                                                                   joint_density<- log(dtruncnorm(theta[1,2:K],a = 0,mean = 0,sd = 1)) - 
                                                                     log(0.5)*((K-1):1)
                                                                   
                                                                   log_p_sum = joint_density + fac + sum(log(P_diag_density))
                                                                   return(log_p_sum)
                                                                 }
                                                               }
                                                                 #Proposal distribution on theta
                                                                 r_d_theta_proposal <- function(theta_prime, sd_proposal, i_star, j_star){
                                                                   mu <- j_star - i_star
                                                                   K<- nrow(theta_prime)
                                                                  
                                                                   #if it's the main diagonal, we have no boundary below it. We fix arbitrarily -9.21
                                                                   
                                                                   bounds = data.frame(mu = 0:(K-1), 
                                                                              lb = c(-9.21,0,theta_prime[1,2:(K-1)]), 
                                                                              ub = c(9.21, theta_prime[1,3:(K)], 9.21))
                                                              
                                                                   theta_scanning_ij = rtruncnorm(n = 1, 
                                                                                                  a = bounds$lb[mu+1], 
                                                                                                  b = bounds$ub[mu+1], 
                                                                                                  mean = theta_prime[i_star,j_star],
                                                                                                  sd =sd_proposal[i_star,j_star])
                                                                   
                                                                   p_prime_given_scanning = dtruncnorm(x = theta_prime[i_star,j_star], 
                                                                                                       a = bounds$lb[mu+1], 
                                                                                                       b = bounds$ub[mu+1], 
                                                                                                       mean = theta_scanning_ij,
                                                                                                       sd =sd_proposal[i_star,j_star])
                                                                   
                                                                   p_scanning_given_prime = dtruncnorm(x = theta_scanning_ij, 
                                                                                                       a = bounds$lb[mu+1], 
                                                                                                       b = bounds$ub[mu+1],  
                                                                                                       mean = theta_prime[i_star,j_star],
                                                                                                       sd =sd_proposal[i_star,j_star])
                                                                   
                                                                   return(list(theta_scanning_ij= theta_scanning_ij,
                                                                               p_prime_given_scanning = p_prime_given_scanning,
                                                                               p_scanning_given_prime = p_scanning_given_prime))
                                                                 
                                                                 
                                                                 
                                                                 
                                                                 
                                                                 
                                                               }
                                                             }
                                                             
                                                             #Prior on theta for the WST model
                                                             #--------------------------------
                                                             if(model =='WST'){
                                                               #double check that the diagonal is fixed to 0.5
                                                               if(diag0.5 == T){
                                                                 #Prior distribution for the WST model
                                                                 
                                                                 theta_prior_probability <- function(theta,K){
                                                                   P = inverse_logit_f(theta)
                                                                   P_upper.tri = P[upper.tri(P,diag = F)]
                                                                   
                                                                   p_prod = dunif(P_upper.tri,0.5,0.9999)*
                                                                     (exp(P_upper.tri)/((1+exp(P_upper.tri))**2))
                                                                   
                                                                   
                                                                   log_p_sum = sum(log(p_prod))
                                                                   return(log_p_sum)
                                                                   
                                                                 }
                                                               }else if(diag0.5 == F){
                                                                 #Prior distribution for the WST model
                                                                 
                                                                 theta_prior_probability <- function(theta,K){
                                                                   P = inverse_logit_f(theta)
                                                                   
                                                                   
                                                                   P_diag = diag(theta)
                                                                   P_diag_density = dbeta(P_diag,shape1 = 2,shape2 = 2)*
                                                                     (P_diag/((1+P_diag)**2))
                                                                   
                                                                   P_upper.tri = P[upper.tri(P,diag = F)]
                                                                   
                                                                   p_prod = dunif(P_upper.tri,0.5,0.9999)*
                                                                     (exp(P_upper.tri)/((1+exp(P_upper.tri))**2))
                                                                   
                                                                   
                                                                   log_p_sum = sum(log(p_prod))
                                                                   return(log_p_sum)
                                                                   
                                                                 }
                                                               }
                                                               #Proposal distribution on theta
                                                               r_d_theta_proposal <- function(theta_prime, sd_proposal, i_star, j_star){
                                                                 
                                                                 if(mu == 0){
                                                                   ub = 9.21
                                                                   lb = -9.21
                                                                 }else{
                                                                   lb = 0 
                                                                   ub = 9.21
                                                                 }
                                                                 theta_scanning_ij = rtruncnorm(n = 1, 
                                                                                                a = lb, b = ub, 
                                                                                                mean = theta_prime[i_star,j_star],
                                                                                                sd =sd_proposal[i_star,j_star])
                                                                 
                                                                 p_prime_given_scanning = dtruncnorm(x = theta_prime[i_star,j_star], 
                                                                                                     a = lb, 
                                                                                                     b = ub, 
                                                                                                     mean = theta_scanning_ij,
                                                                                                     sd =sd_proposal[i_star,j_star])
                                                                 
                                                                 p_scanning_given_prime = dtruncnorm(x = theta_scanning_ij, 
                                                                                                     a = lb, 
                                                                                                     b = ub, 
                                                                                                     mean = theta_prime[i_star,j_star],
                                                                                                     sd =sd_proposal[i_star,j_star])
                                                                 
                                                                 return(list(theta_scanning_ij= theta_scanning_ij,
                                                                             p_prime_given_scanning = p_prime_given_scanning,
                                                                             p_scanning_given_prime = p_scanning_given_prime))
                                                                 
                                                               }
                                                               
                                                               
                                                             }
                                                             
                                                             #Prior on theta for the Simple model
                                                             #--------------------------------
                                                             if(model =='Simple'){
                                                               
                                                               #double check that the diagonal is not fixed to 0.5
                                                               diag0.5 = F
                                                               
                                                               #Prior distribution for the Simple model 
                                                               theta_prior_probability <- function(theta,K){
                                                                 P = inverse_logit_f(theta)
                                                                 P_upper.tri = P[upper.tri(P,diag = F)]
                                                                 
                                                                 p_prod = dunif(P_upper.tri,0.0001,0.9999)*
                                                                   (exp(P_upper.tri)/((1+exp(P_upper.tri))**2))
                                                                 
                                                                 log_p_sum = sum(log(p_prod))
                                                                 
                                                                 return(log_p_sum)
                                                               }
                                                               #Proposal distribution for the Simple model 
                                                               r_d_theta_proposal <- function(theta_prime, sd_proposal, i_star, j_star){
                                                                 lb = -9.21
                                                                 ub = 9.21
                                                                 
                                                                 theta_scanning_ij = rtruncnorm(n = 1, a = lb, b = ub, 
                                                                                                mean = theta_prime[i_star,j_star],
                                                                                                sd = sd_proposal[i_star,j_star])
                                                                 
                                                                 p_prime_given_scanning = dtruncnorm(x = theta_prime[i_star,j_star], 
                                                                                                     a = lb, 
                                                                                                     b = ub, 
                                                                                                     mean = theta_scanning_ij,
                                                                                                     sd = sd_proposal[i_star,j_star])
                                                                 
                                                                 p_scanning_given_prime = dtruncnorm(x = theta_scanning_ij, 
                                                                                                     a = lb, 
                                                                                                     b = ub, 
                                                                                                     mean = theta_prime[i_star,j_star],
                                                                                                     sd =sd_proposal[i_star,j_star])
                                                                 
                                                                 return(list(theta_scanning_ij= theta_scanning_ij,
                                                                             p_prime_given_scanning = p_prime_given_scanning,
                                                                             p_scanning_given_prime = p_scanning_given_prime))
                                                               }
                                                               
                                                             }
                                                             
                                                             #Binomial Log-likelihood function
                                                             #relevant_indices is a logical nxn matrix that makes the code faster:
                                                             #the likelihood is computed just for the entries that
                                                             #1) are in the upper triangular adjacency matrix
                                                             #2) have N_ij > 0 
                                                             
                                                             ll_computation <- function(Y_ij, N_ij, P_nbyn, relevant_indices){
                                                               log_likelihood = dbinom(x = Y_ij[relevant_indices], 
                                                                                       size = N_ij[relevant_indices], 
                                                                                       prob = P_nbyn[relevant_indices],
                                                                                       log = T)
                                                               return(sum(log_likelihood))
                                                             }
                                                             
                                                             #---------------------------
                                                             # Preparatory Steps for the MCMC to run
                                                             #--------------------------------------------------------------
                                                             #setting the seed for reproducibility
                                                             
                                                             seeds = seed + (1:n_chains)*10
                                                             set.seed(seeds[[chain]])
                                                             
                                                             #specific location of saving for the directory
                                                             save_dir = where_to_save[[chain]]
                                                             
                                                             #fixing the number of clusters
                                                             K <- as.numeric(K_est[[chain]])
                                                             #number of samples considering burnin and thinning
                                                             num_samples = (burnin+1):N_iter
                                                             N_iter_eff = sum(num_samples %% thin == 0) 
                                                             
                                                             
                                                             #----------------------------------
                                                             #Computation of the relevant indices
                                                             #---------------------------------
                                                             #Relevant indices are
                                                             # 1) i,j : i<j (upper triangular entries)
                                                             # 2) i,j : N_ij > 0 (non-zero number of interactions)
                                                             
                                                             # Get the indices of the upper triangular part of the matrix
                                                             upper_tri_indices <- matrix(upper.tri(N_ij),n,n,byrow = F)
                                                             
                                                             # Get the indices where the matrix elements are greater than zero
                                                             non_zero_indices <- matrix(N_ij > 0,n,n)
                                                             
                                                             # Find the common indices (both upper.tri and strictly positive)
                                                             relevant_indices <- upper_tri_indices*non_zero_indices
                                                             
                                                             # Convert back to a logical matrix
                                                             relevant_indices <- matrix(as.logical(relevant_indices), nrow = n,ncol = n)
                                                             
                                                             #--------------------------
                                                             # Self -initilization if custom init is not provided
                                                             # --> if you do not provide custom initial values, the MH self-initializes
                                                             if(all(is.na(custom_init))){
                                                               
                                                               # Self-initializing z parameter
                                                               #-------------------------------------------------------------------------
                                                               
                                                               if(estimation_control$z==1){
                                                                 
                                                                 z_current=  sample(x = c(1:K), size = n,replace = T)
                                                               }else{
                                                                 #if the parameters is not to be estimated, initialize it to the "true" value
                                                                 z_current=  matrix(ground_truth$z, n, 1)
                                                               }
                                                               
                                                               
                                                               # Self-initializing theta parameter
                                                               #-------------------------------------------------------------------------
                                                               
                                                               if(estimation_control$theta==1){
                                                                 
                                                                 mu_vec_current <- NA
                                                                 theta_current = matrix(0,K,K)
                                                                 
                                                                 if(model == 'SST'){
                                                                   #initialise mu hyperparameter
                                                                   mu_vec_current<- sort(rtruncnorm(K-1,a = 0, b = 9.20, 
                                                                                                    mean = 0,sd = 1))
                                                                   mu_vec_current = c(0,mu_vec_current)
                                                                   
                                                                   
                                                                   #initialize theta current
                                                                   for(d in 0:(K-1)){
                                                                     theta_current[col(theta_current)-row(theta_current)==d]<- mu_vec_current[d+1]
                                                                   }
                                                                   
                                                                 }else if(model =='WST'){
                                                                   #-------------------------------------------------------------------------
                                                                   # IF ESTIMATING WST OR SIMPLE: Initializing the theta parameter
                                                                   #-------------------------------------------------------------------------
                                                                   
                                                                   theta_current[upper.tri(theta_current)] <- runif(n = K*(K-1)/2,
                                                                                                                    min = 0, 
                                                                                                                    max = 9.20)
                                                                   
                                                                 }else if(model =='Simple'){
                                                                   theta_current[upper.tri(theta_current,diag = T)] <- runif(n = K*(K-1)/2+K,
                                                                                                                             min = -9.20, 
                                                                                                                             max = 9.20)
                                                                   
                                                                 }
                                                                 
                                                                 theta_current[lower.tri(theta_current)] = - t(theta_current)[lower.tri(theta_current)]
                                                               }else{
                                                                 theta_current = matrix(ground_truth$theta, K, K)
                                                                 mu_vec_current = theta_current[1,]
                                                               }
                                                             }else{ 
                                                               #-------------------------------------
                                                               # If custom initialization is chosen
                                                               #-------------------------------------
                                                               z_current = custom_init$z
                                                               theta_current = custom_init$theta
                                                             }
                                                             
                                                             # e_0 <- length(theta_current[upper.tri(theta_current,diag = T)])
                                                             # if(model != 'Simple'){
                                                             #   e_0 = e_0 + length(mu_vec_current)
                                                             # }
                                                             # e_0 = e_0/2+1
                                                             # e_0 = 1
                                                             
                                                             
                                                             if(power_posterior_apprach==T){
                                                               n_temperatures <- 50
                                                               i <- 0:n_temperatures
                                                               t_list <- (i / n_temperatures)^ 5
                                                               
                                                               evidence_df = data.frame(t = t_list, evidence = rep(NA, n_temperatures+1), K_est = rep(K, n_temperatures+1))
                                                             }else{
                                                               #estimating just one chain
                                                               t_list = 1
                                                               t=1
                                                             }
                                                             
                                                             for(t in t_list){
                                                               
                                                               #---------------------------------
                                                               #initializing auxiliary quantities
                                                               
                                                               labels_available<- 1:K
                                                               label_counts <- table(factor(z_current, levels = labels_available))
                                                               n_k = as.numeric(label_counts) #number of nodes in each cluster
                                                               alpha_vec = as.vector(rep(1,K)) #Dirichlet Multinomial weights
                                                               
                                                               #checking that we have exactly K labels
                                                               #if not, sample again z_current until all clusters are non-empty
                                                               while(any(n_k==0)){
                                                                 k_missing = which(n_k == 0)
                                                                 for(i in 1:length(k_missing)){
                                                                   z_current[sample(n, size = n*1/K, replace = F)] <- k_missing[i]
                                                                   
                                                                   label_counts <- table(factor(z_current, levels = labels_available))
                                                                   n_k = as.numeric(label_counts)
                                                                 }
                                                               }
                                                               #---------------------------------
                                                               #initialising priors and likelihood
                                                               #---------------------------------
                                                               z_mat_current = vec2mat(z_current)
                                                               P_current = inverse_logit_f(theta_current)
                                                               P_nbyn_current = calculate_victory_probabilities(z_mat_current, P_current)
                                                               ll_current =  ll_computation(Y_ij = Y_ij, N_ij= N_ij, 
                                                                                            P_nbyn = P_nbyn_current,
                                                                                            relevant_indices = relevant_indices)*t
                                                               
                                                               #log prior on P
                                                               prior_theta<- theta_prior_probability(theta = theta_current, K = K)
                                                               
                                                               #log prior on z
                                                               prior_z <- ddirichlet_multinomial(N = n, K = K, n_k = n_k, my_alpha = alpha_vec)
                                                               
                                                               #computing the whole log proportional posterior
                                                               check <- ll_current + prior_theta + prior_z
                                                               if(is.numeric(check)==T){
                                                                 print("Check completed")
                                                               }else{
                                                                 break
                                                               }
                                                               #-------------------------------
                                                               
                                                               #--------------------------------------------------------------------------
                                                               # PREPARATORY STEPS FOR THE MCMC
                                                               #--------------------------------------------------------------------------
                                                               
                                                               #defining the containers to store results
                                                               ll_container <- matrix(0, ncol = sum(relevant_indices) , nrow = N_iter_eff)
                                                               z_container <- matrix(0, nrow = n, ncol = N_iter_eff)
                                                               mu_vec_container <- matrix(0, nrow = K, ncol = N_iter_eff)
                                                               theta_container <- array(0, dim = c(K,K,N_iter_eff))
                                                               
                                                               #defining the containers to store acceptance counts
                                                               acc.count_z <- rep(1,n)
                                                               acc.count_mu_vec <- rep(1, K)
                                                               acc.count_theta<- matrix(1,K,K)
                                                               
                                                               #defining the proposals' variances
                                                               tau_mu_vec= rep(0.3,K)
                                                               tau_theta = matrix(0.25,K,K)
                                                               
                                                               #-----------------------------------------
                                                               # dataframe of the theta entries that needs to be updated
                                                               #-----------------------------------------
                                                               #Updating each entry of P, one at the time
                                                               ut <- upper.tri(theta_current, diag= !diag0.5) #if diag0.5 ==F, also the main diagonal is updated
                                                               theta_combn = which(ut, arr.ind = TRUE) # get the indices of the upper triangular elements
                                                               
                                                               uo<- data.frame(row = theta_combn[,1], col = theta_combn[,2] ) %>%
                                                                 mutate(diff = col-row) # computing the level sets
                                                               if(model =='SST'){
                                                                 uo = uo%>% 
                                                                   group_by(diff) %>%   # Group by diff to ensure each group has unique diff values
                                                                   slice(1)  # just one scan for each level set
                                                               }
                                                               #-----------------------------------------
                                                               #READY TO GO!
                                                               iteration_time= vector()
                                                               save_count = 0
                                                               for(j in 2:N_iter){
                                                                 
                                                                 start_time <- Sys.time()
                                                                 
                                                                 #z UPDATE-------------------------------------------------------------
                                                                 if (estimation_control$z == 1) {
                                                                   # ------------------------------------------
                                                                   #storing current values
                                                                   # ------------------------------------------
                                                                   P_prime<- inverse_logit_f(theta_current)
                                                                   z_prime= z_current
                                                                   z_mat_prime = vec2mat_0_P(clust_lab = z_prime,K = K)
                                                                   P_nbyn_prime<- calculate_victory_probabilities(z_mat_prime,P_prime)
                                                                   
                                                                   
                                                                   ll_prime<- ll_computation(Y_ij = Y_ij,
                                                                                             N_ij = N_ij,
                                                                                             P_nbyn = P_nbyn_prime,
                                                                                             relevant_indices = relevant_indices)
                                                                   
                                                                   B_prime<- ddirichlet_multinomial(N = n,
                                                                                                    K = K,
                                                                                                    n_k = n_k, 
                                                                                                    my_alpha =  alpha_vec)
                                                                   
                                                                   
                                                                   n_prime = table(factor(z_prime, levels = labels_available))
                                                                   # ------------------------------------------
                                                                   
                                                                   # ------------------------------------------
                                                                   # full sweep of z entries
                                                                   # ------------------------------------------
                                                                   for(i_th_turn in 1:n){
                                                                     
                                                                     z_scanning = z_prime
                                                                     #save current label of z_ii
                                                                     k_prime <- z_prime[i_th_turn]
                                                                     
                                                                     # Sample a new label using the adjusted probabilities
                                                                     k_scanning <- sample(x = setdiff(x = labels_available, 
                                                                                                      y = k_prime), 
                                                                                          size = 1, replace = F)
                                                                     
                                                                     #new proposed partition
                                                                     z_scanning[i_th_turn] <- k_scanning
                                                                     
                                                                     #getting the nodes that interact with node i_th
                                                                     logical_matrix1 = matrix(FALSE, n,n)
                                                                     logical_matrix1[i_th_turn,]<-TRUE #outgoing edges
                                                                     logical_matrix1[,i_th_turn]<-TRUE #incoming edges
                                                                     
                                                                     # Create a matrix that is TRUE only at the positions that:
                                                                     # 1) are in the upper triangle matrix
                                                                     # 2) have N_ij > 0 
                                                                     # 3) Interacted with the i_th node 
                                                                     filtering_matrix = (logical_matrix1)*relevant_indices == T
                                                                     
                                                                     #compute the likelihood of the data with the current assignment just for i_th_turn
                                                                     ll_minus = ll_computation(Y_ij = Y_ij,
                                                                                               N_ij = N_ij,
                                                                                               P_nbyn = P_nbyn_prime,
                                                                                               relevant_indices = filtering_matrix)
                                                                     
                                                                     #update P_nbyn with the new partition
                                                                     P_nbyn_scanning = P_nbyn_prime
                                                                     for(nodes in 1:n){
                                                                       P_nbyn_scanning[i_th_turn,nodes]<- P_prime[k_scanning,z_scanning[nodes]]
                                                                       P_nbyn_scanning[nodes,i_th_turn]<- P_prime[z_scanning[nodes],k_scanning]
                                                                     }
                                                                     
                                                                     #compute the likelihood of the same positions with the new assignment
                                                                     ll_plus = ll_computation(Y_ij = Y_ij,
                                                                                              N_ij = N_ij,
                                                                                              P_nbyn = P_nbyn_scanning,
                                                                                              relevant_indices = filtering_matrix)
                                                                     
                                                                     #Updating the likelihood
                                                                     ll_scanning = (ll_prime - ll_minus + ll_plus)
                                                                     
                                                                     #updating the prior on z
                                                                     n_scanning<- n_prime
                                                                     n_scanning[c(k_prime, k_scanning)] <- n_prime[c(k_prime, k_scanning)] + c(-1, 1)
                                                                     
                                                                     B_scanning<- ddirichlet_multinomial(N = n,
                                                                                                         K = K,
                                                                                                         n_k = n_scanning,
                                                                                                         my_alpha = alpha_vec)
                                                                     # the prior on theta cancels out, therefore is omitted
                                                                     # the proposal is symmetric and proportional to 1/(K-1): it cancels out too
                                                                     
                                                                     # log ratio
                                                                     log_r= t*(ll_scanning - ll_prime) + B_scanning - B_prime
                                                                     
                                                                     #ACCEPT/REJECT CONDITION
                                                                     MH_condition_Z= min(log_r,0)>=log(runif(1))
                                                                     
                                                                     if(MH_condition_Z){
                                                                       acc.count_z[i_th_turn]=acc.count_z[i_th_turn]+1
                                                                       z_prime<-z_scanning
                                                                       ll_prime<- ll_scanning
                                                                       B_prime<- B_scanning
                                                                       P_nbyn_prime <- P_nbyn_scanning
                                                                       n_prime <- n_scanning
                                                                     }
                                                                     
                                                                   }
                                                                   
                                                                   ll_current <- ll_prime
                                                                   P_nbyn_current <- P_nbyn_prime
                                                                   z_current = z_prime
                                                                   label_counts <- table(factor(z_current, levels = labels_available))
                                                                   n_k = as.numeric(label_counts)  
                                                                   
                                                                 }
                                                                 if (estimation_control$theta == 1) {
                                                                   #----------------------------
                                                                   #theta UPDATE----------------
                                                                   #----------------------------
                                                                   
                                                                   #storing current values
                                                                   theta_prime <- theta_current
                                                                   P_prime<- inverse_logit_f(theta_prime)
                                                                   z_mat_prime = vec2mat_0_P(clust_lab = z_current,K = K)
                                                                   P_nbyn_prime<- P_nbyn_current
                                                                   ll_prime<- ll_current
                                                                   theta_prior_prime <- theta_prior_probability(theta = theta_prime, K=K)
                                                                   
                                                                   #update of each upper.triangular entry of theta
                                                                   for(i_th in 1:nrow(uo)){
                                                                     
                                                                     theta_scanning <- theta_prime
                                                                     P_nbyn_scanning = P_nbyn_prime
                                                                     #indices of the entry
                                                                     i_star<- uo$row[i_th]
                                                                     j_star<- uo$col[i_th]
                                                                     #level set
                                                                     mu = uo$diff[i_th]
                                                                     
                                                                     
                                                                     #proposing a new value for theta p_q
                                                                     theta_ij_proposal <- r_d_theta_proposal(theta_prime = theta_prime,
                                                                                                             sd_proposal = tau_theta,
                                                                                                             i_star = i_star, j_star = j_star)
                                                                     #storing new proposed value
                                                                     theta_ij_scanning = theta_ij_proposal$theta_scanning_ij
                                                                     #theta_scanning is the newly proposed theta matrix
                                                                     if(model =='SST'){
                                                                       theta_scanning[col(theta_prime)- row(theta_prime) == mu] <- theta_ij_scanning
                                                                       if(mu != 0){
                                                                         theta_scanning[col(theta_prime)- row(theta_prime) == -mu] <- -theta_ij_scanning
                                                                       }
                                                                     }else{
                                                                       theta_scanning[i_star,j_star] <- theta_ij_scanning
                                                                       if(mu != 0){
                                                                       theta_scanning[j_star,i_star] <- -theta_ij_scanning
                                                                       }
                                                                     }
                                                                     
                                                                     #wi
                                                                     if(model == 'SST'){
                                                                       Z_i_star = which(z_current %in% (1:(K-mu)))
                                                                       Z_j_star = which(z_current %in% ((mu+1):K))
                                                                     }else{
                                                                       Z_i_star = which(z_current %in% i_star)
                                                                       Z_j_star =  which(z_current %in% j_star)
                                                                     }
                                                                     
                                                                     
                                                                     logical_matrix1 = matrix(FALSE, n,n)
                                                                     logical_matrix1[Z_i_star,]<-TRUE
                                                                     logical_matrix1[,Z_j_star]<-TRUE
                                                                     
                                                                     logical_matrix2 = matrix(FALSE, n,n)
                                                                     logical_matrix2[Z_j_star,]<-TRUE
                                                                     logical_matrix2[,Z_i_star]<-TRUE
                                                                     
                                                                     # Create a matrix that is TRUE only at the positions that are both in the upper triangle and in the selected clusters
                                                                     filtering_matrix = (logical_matrix1|logical_matrix2)*relevant_indices == T
                                                                     # Filter the relevant entries
                                                                     
                                                                     #First, subtract the contribution to the likelihood of those items
                                                                     ll_minus = ll_computation(Y_ij = Y_ij, 
                                                                                               N_ij = N_ij, 
                                                                                               P_nbyn = P_nbyn_prime,
                                                                                               relevant_indices = filtering_matrix)
                                                                     
                                                                     #Second, subtract also the contribution of the theta[i_star,j_star] entry to the prior
                                                                     theta_prior_prime = theta_prior_probability(theta_prime, K = K)
                                                                     
                                                                     #recompute the interaction success probabilities, just those that were affected
                                                                     
                                                                     
                                                                     P_scanning = inverse_logit_f(theta_scanning)
                                                                     P_nbyn_scanning = calculate_victory_probabilities(z_mat_prime, P_scanning) 
                                                                     
                                                                     #Third,add the contribution to the likelihood of those items with new interaction probabilities
                                                                     ll_plus = ll_computation(Y_ij = Y_ij, 
                                                                                              N_ij = N_ij, 
                                                                                              P_nbyn = P_nbyn_scanning, 
                                                                                              relevant_indices = filtering_matrix)
                                                                     
                                                                     #Fourth,add the contribution to the prior
                                                                     theta_prior_scanning = theta_prior_probability(theta_scanning, K = K)
                                                                     
                                                                     #Updating the likelihood
                                                                     ll_scanning = ll_prime - ll_minus + ll_plus
                                                                     
                                                                     
                                                                     
                                                                     log_r= t*(ll_scanning - ll_prime) + 
                                                                       theta_prior_scanning - theta_prior_prime + 
                                                                       log(theta_ij_proposal$p_prime_given_scanning) - log(theta_ij_proposal$p_scanning_given_prime)
                                                                     
                                                                     
                                                                     
                                                                     #create statements that check conditiond to accept move
                                                                     MH_update_theta = min(log_r,0)>=log(runif(1))
                                                                     
                                                                     if(MH_update_theta){
                                                                       
                                                                       if(model =='SST'){
                                                                         acc.count_theta[col(theta_prime)- row(theta_prime) == mu] <- acc.count_theta[col(theta_prime)- row(theta_prime) == mu] +1
                                                                       }else{
                                                                         acc.count_theta[i_star,j_star] <- acc.count_theta[i_star,j_star] +1
                                                                       }
                                                                       ll_prime = ll_scanning
                                                                       theta_prime <- theta_scanning
                                                                       P_nbyn_prime = P_nbyn_scanning
                                                                     }
                                                                     
                                                                   }
                                                                   ll_current <- ll_prime
                                                                   theta_current = theta_prime
                                                                   P_nbyn_current = P_nbyn_prime
                                                                 }
                                                                 
                                                                 #---------------------------
                                                                 #storing results
                                                                 #---------------------------
                                                                 if(j > burnin & j%%thin==0){
                                                                   
                                                                   save_count = save_count +1 
                                                                   z_container[,save_count] <- z_current
                                                                   theta_container[,,save_count] <- theta_current
                                                                   if(model=='SST'){
                                                                     mu_vec_container[,save_count] <- theta_current[1,]
                                                                   }
                                                                   
                                                                   ll_container[save_count,] = ll_current
                                                                   
                                                                 }
                                                                 
                                                                 
                                                                 
                                                                 
                                                                 end_time <- Sys.time()
                                                                 
                                                                 iteration_time<-append(iteration_time,as.numeric(difftime(end_time, 
                                                                                                                           start_time, 
                                                                                                                           units = "secs")))
                                                                 
                                                                 if(j%%5000==0){
                                                                   avg_iteration<- mean(iteration_time)
                                                                   current_time <- Sys.time() # Get the current time
                                                                   
                                                                   # Calculate the expected finishing time
                                                                   expected_finishing_time <- current_time + (avg_iteration * (N_iter - j) )
                                                                   formatted_expected_finishing_time <- format(expected_finishing_time, "%H:%M:%S")
                                                                   
                                                                   p(sprintf("Chain %d: mean acc.rate z %.3f%%,
                    mean acc.rate theta %.3f%%,
                    mean acc.rate mu_vec %.3f%% single_iter_time '%*.3f' seconds, will_finish_at '%s'",
                    chain,
                    100 * mean(acc.count_z/j),
                    100 * mean(acc.count_theta/j),
                    100 * mean(acc.count_mu_vec/j),
                    4, avg_iteration, formatted_expected_finishing_time), class = "sticky")
                                                                   
                                                                   
                                                                 }
                                                               }
                                                               
                                                               
                                                               
                                                               if(power_posterior_apprach ==T){
                                                            
                                                                 #obtaining the likelihood for a given iteration
                                                                 LL = colSums(ll_container)
                                                                 
                                                                 #mean likelihood across iterations (expected deviance) of a given temperature t
                                                                 evidence = mean(LL)
                                                                 evidence_df[which(evidence_df$t == t), evidence]
                                                               }         
                                                             }
                                                             
                                                             
                                                             acceptance_rates <- list(acc.count_theta =acc.count_theta, 
                                                                                      acc.count_z = acc.count_z, 
                                                                                      acc.count_mu_vec= acc.count_mu_vec)
                                                             
                                                             est_containers = list(z = z_container,
                                                                                   theta = theta_container, 
                                                                                   mu_vec = mu_vec_container)
                                                             
                                                             control_containers = list(est_model = model,
                                                                                       N_iter = N_iter,
                                                                                       thin = thin,
                                                                                       N_iter_eff = N_iter_eff)
                                                             
                                                             
                                                             if(power_posterior_apprach==F){
                                                             return(list(Y_ij= Y_ij, N_ij = N_ij, 
                                                                         ground_truth=ground_truth,
                                                                         est_containers=est_containers, 
                                                                         control_containers=control_containers, 
                                                                         acceptance_rates= acceptance_rates, 
                                                                         t=t, seed = seed + chain))
                                                             }else{
                                                               evidence_df = evidence_df %>% arrange(t)
                                                               
                                                               for(row_i in 1:(nrow(evidence_df)-1)){
                                                                 ith_sum = (evidence_df$t[row_i+1] - evidence_df$t[row_i])*(evidence_df$evidence[row_i]+complete_df$evidence[row_i+1])/2
                                                                 evidence_df$riemann[row_i] = ith_sum
                                                               }
                                                               
                                                               
                                                               marginal_likelihood = sum(complete_df$riemann)
                                                               
                                                               return(list(Y_ij= Y_ij, N_ij = N_ij, 
                                                                           ground_truth=ground_truth,
                                                                           est_containers=est_containers, 
                                                                           control_containers=control_containers, 
                                                                           acceptance_rates= acceptance_rates, 
                                                                           evidence_df = evidence_df,
                                                                           marginal_likelihood = marginal_likelihood,
                                                                           t=t, seed = seed + chain))
                                                             }
                                                             
                                                           }
    
    
    
  })
  
  
  
  
  
  return(reprex)
} 




