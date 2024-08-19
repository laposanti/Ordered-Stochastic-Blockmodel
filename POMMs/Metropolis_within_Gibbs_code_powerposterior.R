
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
                                                    custom_init=NA,power_posterior_apprach = T, thin=1){
  #setting for each chain a different seed
  #if the given seed is 20, the chains' seeds will be 21 for chain 1, 22 for chain 2 and so on...
  

  
  if(power_posterior_apprach == T){
    n_temperatures=50
  }
  
  where_to_save =list()
  for(i in 1:length(K_est)){
    if(power_posterior_apprach==T){
      where_to_save[[i]] =  file.path(saving_directory, paste0("/MCMC_output/powerposterior/Data_",data_description,"/Est_",model,"/K", K_est[[i]],"/"))
      dir.create(where_to_save[[i]], showWarnings = F, recursive = T)
      }

    else{
        where_to_save[[i]] =  NA
    }
    
  }


  
  n_chains = length(K_est)
  if(detectCores() < n_chains){
    cat('Warning: Number of cores exceeds the ')
  }
  
  variables_to_add = c('Y_ij', 'N_ij' , 'estimation_control', 
                       'ground_truth','n', 'N_iter','n_chains', 
                       'optimal_acceptance_rate_theta', 'optimal_acceptance_rate_mu', 'K_est','thin', 'burnin', 'seed','model','data_description',
                       'power_posterior_apprach' ,'true_model', 'custom_init','p','n_temperatures','where_to_save')
  
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
                                                             
                                                             
                                                             
                                                             seeds = seed + 1:n_chains
                                                             
                                                             set.seed(seeds[[chain]])
                                                             
                                                             save_dir = where_to_save[[chain]]
                                                             #setting hyperparams
                                                             K <- as.numeric(K_est[[chain]])
                                                             
                                                             N_iter_eff = length(seq(1,N_iter-burnin, thin))
                                                             
                                                             
                                                             #if you do not provide custom initial values, the MH auto initialises starting from the seed
                                                             if(all(is.na(custom_init))){
                                                               #-------------------------------------------------------------------------
                                                               #Initializing z allocation vector 
                                                               #-------------------------------------------------------------------------
                                                               
                                                               #if the parameters is fixed, setting it to the true value
                                                               if(estimation_control$z==1){
                                                                 
                                                                 z_current=  sample(x = c(1:K), size = n,replace = T)
                                                                 
                                                               }else{
                                                                 z_current=  matrix(ground_truth$z, n, 1)
                                                               }
                                                      
                                                               #-------------------------------------------------------------------------
                                                               #Initializing mu parameter for different models
                                                               #-------------------------------------------------------------------------
                                                               
                                                               # if model WST ore SST is selected: initialise the means of the level sets mu_vec
                                                               if(model == 'SST'){
                                                                 if(estimation_control$mu_vec==1){
                                                                   
                                                                   mu_vec_current<- sort(rtruncnorm(K-1,a = 0, b = 10, mean = 0,sd = 1.5))
                                                                   mu_vec_current = c(0,mu_vec_current)
                                                                   theta_current = matrix(NA,K,K)
                                                                   if(model =='SST'){
                                                                     for(d in 0:(K-1)){
                                                                       theta_current[col(theta_current)-row(theta_current)==d]<- mu_vec_current[d+1]
                                                                     }
                                                                   }
                                                                   theta_current[lower.tri(theta_current)] = - t(theta_current)[lower.tri(theta_current)]
                                                                 }else{
                                                                   mu_vec_current=  as.numeric(ground_truth$mu_vec)
                                                                 }
                                                               }else if(model != 'SST'){
                                                                 mu_vec_current =NA
                                                               }
                                                               #-------------------------------------------------------------------------
                                                               #Initializing theta matrix for different models
                                                               #-------------------------------------------------------------------------
                                                               if(estimation_control$theta==1 & model != 'SST'){
                                                                 theta_current = matrix(NA,K,K)
                                                                 diag(theta_current) =0
                                                               if(model =='WST'){
                                                                   for(d in 1:(K-1)){
                                                                     theta_current[col(theta_current)-row(theta_current)==d]<- runif(K-d,0, 9.21023)
                                                                   }
                                                                 }else if( model =='Simple'){
                                                                   for(d in 1:(K-1)){
                                                                     theta_current[col(theta_current)-row(theta_current)==d]<- runif(K-d, min  = -9.21023, 
                                                                                                                                     max = 9.21023)
                                                                   }
                                                                 }
                                                                 theta_current[lower.tri(theta_current)] = - t(theta_current)[lower.tri(theta_current)]
                                                               }
                                                               
                                                             }else{ #the custom initialisation, if provided
                                                               z_current = custom_init$z
                                                               theta_current = custom_init$theta
                                                               if(model == 'SST'|| model =='WST'){
                                                                 mu_vec_current = custom_init$mu_vec
                                                               }
                                                               if(model == 'WST'){
                                                                 sigma_current = custom_init$sigma
                                                               }
                                                             }
                                                             
                                                             e_0 <- length(theta_current[upper.tri(theta_current,diag = T)])
                                                             if(model != 'Simple'){
                                                               e_0 = e_0 + length(mu_vec_current)
                                                             }
                                                             e_0 = e_0/2+1
                                                             # e_0 = 1
                                                             alpha_vec = as.vector(rep(e_0,K))
                                                             
                                                             if(power_posterior_apprach==T){
                                                               i <- 0:n_temperatures
                                                               t_list <- (i / n_temperatures) ^ 5
                                                             }else{
                                                               #estimating just one chain
                                                               t_list = 1
                                                               t=1
                                                             }
                                                             
                                                             for(t in t_list){
                                                               
                                                               
                                                               #initializing quantities
                                                               
                                                               labels_available<- 1:K
                                                               #checking that we have exactly K labels
                                                               label_counts <- table(factor(z_current, levels = labels_available))
                                                               n_k = as.numeric(label_counts)
                                                               
                                                               while(any(n_k==0)){
                                                                 k_missing = which(n_k == 0)
                                                                 for(i in 1:length(k_missing)){
                                                                   z_current[sample(n, size = n*1/K, replace = F)] <- k_missing[i]
                                                                   
                                                                   label_counts <- table(factor(z_current, levels = labels_available))
                                                                   n_k = as.numeric(label_counts)
                                                                 }
                                                               }
                                                               
                                                               z_P = vec2mat(z_current)
                                                               
                                                               labels_available<- 1:K
                                                               #checking that we have exactly K labels
                                                               label_counts <- table(factor(z_current, levels = labels_available))
                                                               n_k = as.numeric(label_counts)
                                                               while(any(n_k==0)){
                                                                 k_missing = which(n_k == 0)
                                                                 for(i in 1:length(k_missing)){
                                                                   z_current[sample(n, size = n*1/K, replace = F)] <- k_missing[i]
                                                                   
                                                                   label_counts <- table(factor(z_current, levels = labels_available))
                                                                   n_k = as.numeric(label_counts)
                                                                 }
                                                               }
                                                               
                                                               z_P = vec2mat(z_current)
                                                               # number of victories between block p and block q
                                                               # number of victories between block p and block q
                                                               # ybar = t(z_P)%*%(Y_ij*upper.tri(Y_ij))%*%z_P
                                                               # # number of missed victories between block p and block q
                                                               # n_minus_y1 <- (N_ij-Y_ij)*upper.tri(N_ij)
                                                               # # number of missed victories between block p and block q
                                                               # mbar<- t(z_P)%*%n_minus_y1%*%z_P
                                                               # 
                                                               # coef1 = lchoose(N_ij, Y_ij)*upper.tri(N_ij)
                                                               # lamdabar <- t(z_P)%*%(coef1)%*%z_P
                                                               # 
                                                               
                                                               
                                                               
                                                               A_current =  ll_naive(z = z_current, theta = theta_current, 
                                                                                     Y_ij = Y_ij, N_ij= N_ij)*t
                                                               
                                                               check = lprop_posterior(z = z_current, 
                                                                                       Y_ij = Y_ij,
                                                                                       N_ij = N_ij,
                                                                                       theta =theta_current,
                                                                                       alpha_vec = alpha_vec,
                                                                                       mu_vec = mu_vec_current,
                                                                                       n_k = n_k,
                                                                                       K = K, 
                                                                                       model = model,  
                                                                                       t=t)
                                                               if(is.numeric(check)==T){
                                                                 print("Check completed")
                                                               }else{
                                                                 break
                                                               }
                                                               #--------------------------------------------------------------------------
                                                               #setting and initialising containers
                                                               #--------------------------------------------------------------------------
                                                               
                                                               #initialising the chain
                                                               A_container <- matrix(0, nrow = 1, ncol = N_iter_eff)
                                                               
                                                               A_container[1] <- A_current
                                                               #initialising the chain
                                                               z_container <- matrix(0, nrow = n, ncol = N_iter_eff)
                                                               z_container[,1] <- z_current
                                                              
                                                               if(model == 'WST'||model == 'SST'){
                                                                 #initialising the chain
                                                                 mu_vec_container = matrix(0, nrow = K, ncol = N_iter_eff)
                                                                 mu_vec_container[,1] <- mu_vec_current
                                                                 #initialising the adaptive variance
                                                                 tau_mu_vec <- rep(0.3, K)
                                                                 tau_mu_vec_container = matrix(0,K, N_iter)
                                                                 tau_mu_vec_container[,1] <- tau_mu_vec
                                                               }else{
                                                                 mu_vec_container <- NA
                                                                 tau_mu_vec <- NA
                                                                 tau_mu_vec_container <- NA
                                                               }
                                                               #initialing theta and its adaptive variance container
                                                               theta_container = array(0, dim = c(K,K,N_iter_eff))
                                                               theta_container[,,1] <- theta_current
                                                               
                                                               tau_theta_container = array(0,dim=c(K,K,N_iter_eff))
                                                               tau_theta =matrix(0.25,K,K)
                                                               tau_theta_container[,,1] = tau_theta
                                                               
                                                               
                                                               
                                                               #containers for the counts of accepted proposals
                                                               acc.count_z = rep(1,n)
                                                               acc.count_sigma_squared=1
                                                               acc.count_mu_vec = rep(1, K)
                                                               acc.count_theta<- matrix(1,K,K)
                                                               
                                                               #READY TO BOMB!
                                                               iteration_time= vector()
                                                               save_count = 1
                                                               for(j in 2:N_iter){
                                                                 
                                                                 start_time <- Sys.time()
                                                                 
                                                                 
                                                                 if (estimation_control$z == 1) {
                                                                   #z UPDATE-------------------------------------------------------------
                                                                   
                                                                   
                                                                   z_update = z_update_f(z = z_current, N_ij = N_ij, 
                                                                                         Y_ij = Y_ij,theta =theta_current,
                                                                                         lamdabar = lamdabar,ybar=ybar,
                                                                                         mbar=mbar,alpha_vec = alpha_vec,
                                                                                         n_k = n_k,K = K,
                                                                                         acc.count_z = acc.count_z,
                                                                                         labels_available = labels_available,
                                                                                         model, t=t)
                                                                   
                                                                   
                                                                   llik = z_update$A_prime
                                                                   z_current = z_update$z
                                                                   acc.count_z = z_update$acc.moves
                                                                   
                                                                   label_counts <- table(factor(z_current, levels = labels_available))
                                                                   n_k = as.numeric(label_counts)  
                                                                   
                                                                 }
                                                                 if (estimation_control$theta == 1) {
                                                                   #theta UPDATE-------------------------------------------------------------
                                                                   
                                                                   theta_update = theta_update_f(z = z_current, N_ij = N_ij, 
                                                                                                 Y_ij = Y_ij,
                                                                                                 theta = theta_current,alpha_vec = alpha_vec,
                                                                                                 n_k = n_k,
                                                                                                 mu_vec = mu_vec_current,
                                                                                                 K = K,tau_theta =tau_theta,
                                                                                                 acc.count_theta =acc.count_theta,model=model,t=t)
                                                                   llik = theta_update$llik
                                                                   theta_current = theta_update$theta
                                                                   acc.count_theta =theta_update$acc.moves
                                                                   
                                                                   
                                                                   # if(j %% 50 == 0){
                                                                   #   
                                                                   #   for(my_p in 1:K){
                                                                   #     for(my_q in my_p:K){
                                                                   #       
                                                                   #       tau_theta[my_p,my_q] = tuning_proposal(iteration = j,
                                                                   #                                              acceptance_count = acc.count_theta[my_p,my_q],
                                                                   #                                              sigma = tau_theta[my_p,my_q],
                                                                   #                                              acceptanceTarget = optimal_acceptance_rate_theta,
                                                                   #                                              min_sigma = 0.00002)
                                                                   #       
                                                                   #     }
                                                                   #   }
                                                                   # }
                                                                 }
                                                                 
                                                                 if(estimation_control$mu== 1) {
                                                                   #mu UPDATE----------------------------------------------------------------
                                                                  
                                                                   mu_update=  mu_update_f(z = z_current, N_ij = N_ij, llik=llik,
                                                                                           Y_ij = Y_ij,  theta =theta_current,
                                                                                           alpha_vec =  alpha_vec, n_k = n_k,
                                                                                           mu_vec = mu_vec_current,K = K, 
                                                                                           tau_mu_vec = tau_mu_vec,
                                                                                           acc.count_mu_vec,model,t=t)
                                                                   
                                                                   #updating quantities
                                                                   mu_vec_current = mu_update$mu_vec
                                                                   acc.count_mu_vec = mu_update$acc.moves
                                                                   theta_current = mu_update$theta
                                                                   
                                                                   
                                                                   # 
                                                                   # if(j %% 50 == 0){
                                                                   #   tau_mu_vec <- tuning_proposal(iteration=j,acceptance_count = acc.count_mu_vec,
                                                                   #                                 sigma = tau_mu_vec,
                                                                   #                                 acceptanceTarget = optimal_acceptance_rate_mu,
                                                                   #                                 min_sigma = 0.002)
                                                                   # }
                                                                   # 
                                                                   
                                                                   
                                                                 }
                                                                 
                                                                 
                                                                 #storing scales
                                                                 # tau_sigma_squared_container[j]<- tau_sigma_squared
                                                                 # tau_mu_vec_container[j]<- tau_mu_vec
                                                                 # tau_theta_container[,,j]<- tau_theta
                                                                 # 
                                                                 #storing results for inference
                                                                 
                                                                 if(j > burnin & N_iter%%thin==0){
                                                                   
                                                                   z_container[,save_count] <- z_current
                                                                   theta_container[,,save_count] <- theta_current
                                                                   if(model=='SST'){
                                                                     mu_vec_container[,save_count] <- mu_vec_current
                                                                   }
                                                                   save_count = save_count +1 
                                                                 }
                                                                 
                                                                 # #storing scales
                                                                 # tau_sigma_squared_container[j]<- tau_sigma_squared
                                                                 # tau_mu_vec_container[j]<- tau_mu_vec
                                                                 # tau_theta_container[,,j]<- tau_theta
                                                                 # 
                                                                 # 
                                                                 
                                                                 
                                                                 
                                                                 end_time <- Sys.time()
                                                                 
                                                                 iteration_time<-append(iteration_time,as.numeric(difftime(end_time, start_time, units = "secs")))
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
                                                                 acceptance_rates <- list(acc.count_theta = acc.count_theta, 
                                                                                          acc.count_z = acc.count_z,
                                                                                          acc.count_mu_vec= acc.count_mu_vec)
                                                                 
                                                                 
                                                                 
                                                                 est_containers = list(z = z_container,theta = theta_container,
                                                                                       mu_vec = mu_vec_container)
                                                                 
                                                                 control_containers = list(A = A_container,
                                                                                           N_iter = N_iter,
                                                                                           thin = thin,
                                                                                           N_iter_eff = N_iter_eff)
                                                                 
                                                                 st.deviations<- list(tau_theta =tau_theta_container,
                                                                                      tau_mu_vec= tau_mu_vec_container)
                                                                 
                                                                 chains = list(Y_ij= Y_ij, N_ij = N_ij, ground_truth=ground_truth,est_containers=est_containers,
                                                                               control_containers=control_containers, acceptance_rates= acceptance_rates,
                                                                               st.deviations=st.deviations, t=t, seed= seed + chain)
                                                                 
                                                                 #storing the results of each chain
                                                                 my_names <- paste0("chain")
                                                                 
                                                                 my_filename <- paste0(save_dir,data_description, "Est_model",model,"_estK_",
                                                                                       K,"_N", n, "iteration", which(t == t_list),".RDS")
                                                                 saveRDS(object = chains, file = my_filename) # saving results
                                                                 
                                                               }         
                                                             }
                                                             
                                                             
                                                             acceptance_rates <- list(acc.count_theta =acc.count_theta, acc.count_z = acc.count_z, 
                                                                                      acc.count_mu_vec= acc.count_mu_vec)
                                                             
                                                             st.deviations<- list(tau_theta =tau_theta_container,
                                                                                  tau_mu_vec= tau_mu_vec_container)
                                                             
                                                             est_containers = list(z = z_container,theta = theta_container, mu_vec = mu_vec_container)
                                                             
                                                             control_containers = list(A = A_container, alpha_vec = alpha_vec)
                                                             
                                                             return(list(Y_ij= Y_ij, N_ij = N_ij, ground_truth=ground_truth,est_containers=est_containers, 
                                                                         control_containers=control_containers, acceptance_rates= acceptance_rates, 
                                                                         st.deviations=st.deviations, t=t, seed = seed + chain))
                                                             
                                                             
                                                           }
    
    
    
  })
  
  
  
  
  
  return(reprex)
} 

