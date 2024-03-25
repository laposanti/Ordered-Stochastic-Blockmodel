
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
                                     ground_truth,n, N_iter,n_chains, true_model,
                                     optimal_acceptance_rate_P, optimal_acceptance_rate_mu, K, seed,est_model,n_temperatures, saving_directory, custom_init=NA){
  
  variables_to_add = c('Y_ij', 'N_ij' , 'estimation_control', 
                       'ground_truth','n', 'N_iter','n_chains', 
                       'optimal_acceptance_rate_P', 'optimal_acceptance_rate_mu', 'K', 'seed','est_model', 'true_model', 'custom_init','p','n_temperatures','saving_directory')
  
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
                                                             
                                                             
                                                             
                                                             
                                                             
                                                             
                                                             set.seed(seed[[chain]])
                                                             
                                                             #setting hyperparams
                                                             K <- as.numeric(K[[chain]])
                                                             alpha_vec = as.vector(rep(1/K,K))
                                                             
                                                             #if you do not provide custom initial values, the MH auto initialises starting from the seed
                                                             if(all(is.na(custom_init))){
                                                               #-------------------------------------------------------------------------
                                                               #Initializing z allocation vector 
                                                               #-------------------------------------------------------------------------
                                                               
                                                               #if the parameters is fixed, setting it to the true value
                                                               if(estimation_control$z==1){
                                                                 
                                                                 z_current=  sample(1:K, n, replace=T)
                                                                 while(length(unique(z_current))!= K){
                                                                   z_current=  sample(1:K, n, replace=T)
                                                                 }
                                                               }else{
                                                                 z_current=  matrix(ground_truth$z, n, 1)
                                                               }
                                                               #-------------------------------------------------------------------------
                                                               #Initializing sigma_squared boundaries 
                                                               #-------------------------------------------------------------------------
                                                               # if model WST is selected: initialise sigma_squared
                                                               if(est_model == 'WST'){
                                                                 if(estimation_control$sigma_squared==1){
                                                                   sigma_squared_current <- runif(1,0.01,0.4)
                                                                 }else{
                                                                   sigma_squared_current <- as.numeric(ground_truth$sigma_squared)
                                                                 }
                                                               }else{
                                                                 sigma_squared_current <-  NA
                                                               }
                                                               #-------------------------------------------------------------------------
                                                               #Initializing mu parameter for different models
                                                               #-------------------------------------------------------------------------
                                                               # if model WST ore SST is selected: initialise the means of the level sets mu_vec
                                                               if(est_model == 'WST' || est_model == 'SST'){
                                                                 if(estimation_control$mu_vec==1){
                                                                   mu_vec_1_K_1<- sort(rtruncnorm(K,a = 0, b = 10, mean = 0,sd = 1))
                                                                   mu_vec0 = rtruncnorm(1,a = -Inf, b = min(mu_vec_1_K_1), mean = 0,sd = 1)
                                                                   mu_vec_current = c(mu_vec0,mu_vec_1_K_1)
                                                                 }else{
                                                                   mu_vec_current=  as.numeric(ground_truth$mu_vec)
                                                                 }
                                                               }else if(est_model == 'Simple'){
                                                                 mu_vec_current =NA
                                                               }
                                                               
                                                               #-------------------------------------------------------------------------
                                                               #Initializing P matrix for different models
                                                               #-------------------------------------------------------------------------
                                                               if(estimation_control$P==1){
                                                                 P_current = matrix(NA,K,K)
                                                                 if(est_model =='SST'){
                                                                   for(d in 0:(K-1)){
                                                                     P_current[col(P_current)-row(P_current)==d]<- runif(K-d, 
                                                                                                                         min  = mu_vec_current[d+1] , 
                                                                                                                         max = mu_vec_current[d+2])
                                                                   }
                                                                 }else if(est_model =='WST'){
                                                                   sigma=sqrt(sigma_squared_current)
                                                                   for(d in 0:(K-1)){
                                                                     P_current[col(P_current)-row(P_current)==d]<- runif(K-d, min  = mu_vec_current[d+1]-sigma , 
                                                                                                                         max = mu_vec_current[d+2]+sigma)
                                                                   }
                                                                 }else if( est_model =='Simple'){
                                                                   for(d in 0:(K-1)){
                                                                     P_current[col(P_current)-row(P_current)==d]<- runif(K-d, min  = -5 , 
                                                                                                                         max = +10)
                                                                   }
                                                                 }
                                                                 P_current[lower.tri(P_current)] = - t(P_current)[lower.tri(P_current)]
                                                               }else{
                                                                 P_current=  as.matrix(ground_truth$P)
                                                               }
                                                               
                                                             }else{#the custom initialisation, if provided
                                                               z_current = custom_init$z
                                                               P_current = custom_init$P
                                                               if(est_model == 'SST'|| est_model =='WST'){
                                                                 mu_vec_current = custom_init$mu_vec
                                                               }
                                                               if(est_model == 'WST'){
                                                                 sigma_current = custom_init$sigma
                                                               }
                                                             }
                                                             
                                                             
                                                             
                                                            
                                                             
                                                             
                                                             
                                                             i <- 0:n_temperatures
                                                             t_list <- (i / n_temperatures) ^ 5
                                                             tau_sigma_squared <- 0.2 
                                                             tau_mu_vec <- 0.05
                                                             tau_P = matrix(0.2,K,K)
                                                             
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
                                                               # number of victories between block p and block q
                                                               # number of victories between block p and block q
                                                               ybar = t(z_P)%*%(Y_ij*upper.tri(Y_ij))%*%z_P
                                                               # number of missed victories between block p and block q
                                                               n_minus_y1 <- (N_ij-Y_ij)*upper.tri(N_ij)
                                                               # number of missed victories between block p and block q
                                                               mbar<- t(z_P)%*%n_minus_y1%*%z_P
                                                               
                                                               coef1 = lchoose(N_ij, Y_ij)*upper.tri(N_ij)
                                                               lamdabar <- t(z_P)%*%(coef1)%*%z_P
                                                               
                                                               
                                                               
                                                               
                                                               A_current = llik_over_blocks_f_binomial(lamdabar = lamdabar, ybar = ybar,
                                                                                                       mbar = mbar,P = P_current,K = K, t=t)
                                                               
                                                               check = lprop_posterior_withP(lamdabar = lamdabar,
                                                                                             ybar = ybar,
                                                                                             mbar = mbar,
                                                                                             P = P_current,
                                                                                             alpha_vec = alpha_vec,
                                                                                             n_k = n_k,
                                                                                             sigma_squared = sigma_squared_current,
                                                                                             mu_vec = mu_vec_current,
                                                                                             K = K, 
                                                                                             model = est_model,  t=t)
                                                               #--------------------------------------------------------------------------
                                                               #setting and initialising containers
                                                               #--------------------------------------------------------------------------
                                                               #initialising the chain
                                                               z_container= matrix(0, nrow = n, ncol = N_iter)
                                                               z_container[,1] = z_current
                                                               if(est_model == 'WST'){
                                                                 #initialising the chain
                                                                 sigma_squared_container =  matrix(0, nrow = 1, ncol = N_iter)
                                                                 sigma_squared_container[1] <- sigma_squared_current
                                                                 #initialising the adaptive variance
                                                                 
                                                                 tau_sigma_squared_container = matrix(0,1, N_iter)
                                                                 tau_sigma_squared_container[1] <- tau_sigma_squared
                                                               }else{
                                                                 tau_sigma_squared <- NA
                                                                 sigma_squared_container <- NA
                                                                 tau_sigma_squared_container <- NA
                                                               }
                                                               if(est_model == 'WST'||est_model == 'SST'){
                                                                 #initialising the chain
                                                                 mu_vec_container = matrix(0, nrow = K+1, ncol = N_iter)
                                                                 mu_vec_container[,1] <- mu_vec_current
                                                                 #initialising the adaptive variance
                                                                 
                                                                 tau_mu_vec_container = matrix(0,1, N_iter)
                                                                 tau_mu_vec_container[1] <- tau_mu_vec
                                                               }else{
                                                                 mu_vec_container <- NA
                                                                 tau_mu_vec <- NA
                                                                 tau_mu_vec_container <- NA
                                                               }
                                                               #initialing P and its adaptive variance container
                                                               P_container = array(0, dim = c(K,K,N_iter))
                                                               P_container[,,1] <- P_current
                                                               
                                                               tau_P_container = array(0,dim=c(K,K,N_iter))
                                                               
                                                               tau_P_container[,,1] = tau_P
                                                               
                                                               
                                                               A_container= matrix(0, nrow=1, ncol=N_iter)
                                                               A_container[1] <- A_current
                                                               
                                                               #containers for the counts of accepted proposals
                                                               acc.count_z = rep(1,n)
                                                               acc.count_sigma_squared=1
                                                               acc.count_mu_vec = 1
                                                               acc.count_P<- matrix(1,K,K)
                                                               
                                                               #READY TO BOMB!
                                                               iteration_time= vector()
                                                               for(j in 2:N_iter){
                                                                 start_time <- Sys.time()
                                                                 
                                                                 
                                                                 if (estimation_control$z == 1) {
                                                                   #z UPDATE-------------------------------------------------------------
                                                                   
                                                                   z_update = z_update_f_withP(z = z_current, N_ij = N_ij, Y_ij = Y_ij,P = P_current,
                                                                                               lamdabar = lamdabar,ybar=ybar,mbar=mbar,alpha_vec = alpha_vec,
                                                                                               n_k = n_k,K = K,
                                                                                               acc.count_z = acc.count_z,labels_available = labels_available, model = est_model, t=t)
                                                                   
                                                                   lamdabar = z_update$lamdabar
                                                                   ybar = z_update$ybar
                                                                   mbar = z_update$mbar
                                                                   
                                                                   z_current = z_update$z
                                                                   acc.count_z = z_update$acc.moves
                                                                   
                                                                   
                                                                   
                                                                   
                                                                 }
                                                                 if (estimation_control$P == 1) {
                                                                   #P UPDATE-------------------------------------------------------------
                                                                   
                                                                   P_update = P_update_f(lamdabar = lamdabar,ybar = ybar,mbar = mbar,
                                                                                         P = P_current,alpha_vec = alpha_vec,
                                                                                         n_k = n_k,sigma_squared = sigma_squared_current,
                                                                                         mu_vec = mu_vec_current,K = K,tau_P = tau_P,
                                                                                         acc.count_P = acc.count_P,model = est_model,t=t)
                                                                   
                                                                   P_current = P_update$P
                                                                   acc.count_P = P_update$acc.moves
                                                                   
                                                                   if(j %% 50 == 0){
                                                                     
                                                                     for(my_p in 1:K){
                                                                       for(my_q in my_p:K){
                                                                         
                                                                         tau_P[my_p,my_q] = tuning_proposal(iteration=j,acceptance_count = acc.count_P[my_p,my_q],
                                                                                                            sigma = tau_P[my_p,my_q],
                                                                                                            acceptanceTarget = optimal_acceptance_rate_P,
                                                                                                            min_sigma = 0.00002)
                                                                         
                                                                       }
                                                                     }
                                                                   }
                                                                 }
                                                                 
                                                                 if (estimation_control$sigma_squared == 1) {
                                                                   #sigma_squared UPDATE----------------------------------------------------------------
                                                                   
                                                                   sigma_squared_update <- sigma_squared_update_f_withP(lamdabar = lamdabar, P = P_current,
                                                                                                                        ybar = ybar,mbar=mbar,
                                                                                                                        alpha_vec = alpha_vec, n_k = n_k,
                                                                                                                        sigma_squared = sigma_squared_current, 
                                                                                                                        mu_vec = mu_vec_current,
                                                                                                                        K = K, tau_sigma_squared = tau_sigma_squared,
                                                                                                                        acc.count_sigma_squared = acc.count_sigma_squared,model = est_model,t=t)
                                                                   #updating quantities
                                                                   
                                                                   acc.count_sigma_squared = sigma_squared_update$acc.moves
                                                                   sigma_squared_current = sigma_squared_update$sigma_squared
                                                                   if(j %% 50 == 0){
                                                                     tau_sigma_squared <- tuning_proposal(iteration=j,
                                                                                                          acceptance_count = acc.count_sigma_squared,
                                                                                                          sigma = tau_sigma_squared,
                                                                                                          acceptanceTarget = optimal_acceptance_rate_P,
                                                                                                          min_sigma = 0.02)
                                                                   }
                                                                 }
                                                                 if (estimation_control$mu== 1) {
                                                                   #P UPDATE----------------------------------------------------------------
                                                                   
                                                                   mu_update=  mu_update_f_withP(lamdabar = lamdabar, ybar = ybar,mbar = mbar,  P = P_current,
                                                                                                 alpha_vec =  alpha_vec, n_k = n_k,
                                                                                                 sigma_squared = sigma_squared_current, 
                                                                                                 mu_vec = mu_vec_current,K = K, tau_mu_vec = tau_mu_vec,
                                                                                                 acc.count_mu_vec,model = est_model,t=t)
                                                                   #updating quantities
                                                                   mu_vec_current = mu_update$mu_vec
                                                                   acc.count_mu_vec = mu_update$acc.moves
                                                                   
                                                                   
                                                                   
                                                                   if(j %% 50 == 0){
                                                                     tau_mu_vec <- tuning_proposal(iteration=j,acceptance_count = acc.count_mu_vec,
                                                                                                   sigma = tau_mu_vec,
                                                                                                   acceptanceTarget = optimal_acceptance_rate_mu,
                                                                                                   min_sigma = 0.002)
                                                                   }
                                                                   
                                                                   
                                                                   
                                                                 }
                                                                 
                                                                 #storing scales
                                                                 tau_sigma_squared_container[j]<- tau_sigma_squared
                                                                 tau_mu_vec_container[j]<- tau_mu_vec
                                                                 tau_P_container[,,j]<- tau_P
                                                                 
                                                                 #storing results for inference
                                                                 A_container[j] = llik_over_blocks_f_binomial( lamdabar = lamdabar, ybar = ybar,mbar = mbar, P = P_current, K=K, t=1)
                                                                 z_container[,j] <- z_current
                                                                 P_container[,,j] <- P_current
                                                                 if(est_model == 'WST'){
                                                                   sigma_squared_container[1,j] = sigma_squared_current
                                                                 }
                                                                 if(est_model == 'WST'|| est_model=='SST'){
                                                                   mu_vec_container[,j] <- mu_vec_current
                                                                 }
                                                                 
                                                                 
                                                                 
                                                                 
                                                                 
                                                                 end_time <- Sys.time()
                                                                 
                                                                 iteration_time<-append(iteration_time,as.numeric(difftime(end_time, start_time, units = "secs")))
                                                                 
                                                                 if(j%%5000==0){
                                                                   avg_iteration<- mean(iteration_time)
                                                                   current_time <- Sys.time() # Get the current time
                                                                   
                                                                   # Calculate the expected finishing time
                                                                   expected_finishing_time <- current_time + (avg_iteration * (N_iter - j) )
                                                                   formatted_expected_finishing_time <- format(expected_finishing_time, "%H:%M:%S")
                                                                   
                                                                   p(sprintf("Chain %d, iteration %d out of %d, mean acc.rate z %.3f%%, 
                    mean acc.rate P %.3f%%,
                    mean acc.rate sigma_squared %.3f%%,
                    mean acc.rate mu_vec %.3f%% single_iter_time '%*.3f' seconds, will_finish_at '%s'",
                    chain, which(t==t_list), n_temperatures,
                    100 * mean(acc.count_z/j),
                    100 * mean(acc.count_P/j),
                    100 * mean(acc.count_sigma_squared/j),
                    100 * mean(acc.count_mu_vec/j),
                    4, avg_iteration, formatted_expected_finishing_time), class = "sticky")
                                                                   
                                                                   
                                                                 }
                                                               }
                                                               
                                                               acceptance_rates <- list(acc.count_P = acc.count_P, acc.count_z = acc.count_z,
                                                                                        acc.count_sigma_squared=acc.count_sigma_squared, acc.count_mu_vec= acc.count_mu_vec)
                                                               
                                                               st.deviations<- list(tau_P = tau_P_container,tau_sigma_squared = tau_sigma_squared_container, 
                                                                                    tau_mu_vec= tau_mu_vec_container)
                                                               
                                                               est_containers = list(z = z_container,P = P_container,
                                                                                     sigma_squared= sigma_squared_container, mu_vec = mu_vec_container)
                                                               
                                                               control_containers = list(A = A_container)
                                                               
                                                               chains = list(Y_ij= Y_ij, N_ij = N_ij,
                                                                             seed = seed, ground_truth=ground_truth,est_containers=est_containers, 
                                                                             control_containers=control_containers, acceptance_rates= acceptance_rates, 
                                                                             st.deviations=st.deviations, t=t, seed=seed)
                                                               
                                                               my_names <- paste0("chain")
                                                               where_to_save = saving_directory[[chain]]
                                                               my_filename <- paste0(where_to_save,true_model, "Est_model",est_model,"_estK_", K,"_N", n, "iteration", which(t == t_list),".RDS")
                                                               saveRDS(object = chains, file = my_filename) # saving results
                                                               
                                                               
                                                               
                                                             }
                                                             
                                                             
                                                             return(list(Y_ij= Y_ij, N_ij = N_ij,
                                                                         seed = seed, ground_truth=ground_truth,est_containers=est_containers, 
                                                                         control_containers=control_containers, acceptance_rates= acceptance_rates, 
                                                                         st.deviations=st.deviations, t=t, seed=seed))
                                                           }
    
    
    
    
  })

  return(reprex)
} 

