
#making inference on 
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

adaptive_MCMC_orderstats <- function(Y_ij, N_ij,init , estimation_control, 
                                     ground_truth,N, N_iter,n_chains, 
                                     optimal_acceptance_rate, hyper_params, seed){
  
  handlers(global = TRUE)
  registerDoFuture()

  reprex <- local({
    p <- progressor(steps = n_chains * N_iter/5000)
    plan(multisession, workers= min(n_chains,(detectCores()-2)))
    y <- foreach(chain = 1:n_chains,.options.future = list(globals = structure(TRUE, add=c('init','ground_truth')), seed=TRUE)) %dofuture% { 

      init <- init[[chain]]
      ground_truth<- ground_truth[[chain]]
      N=n
      #setting hyperparams
      K_max = as.numeric(hyper_params$K_max)
      alpha_vec = as.vector(hyper_params$alpha_vec)
      A_container= matrix(0, nrow=1, ncol=N_iter)
      K = nrow(init$P)
      #setting containers
      z_container= matrix(0, nrow = N, ncol = N_iter)
      sigma_squared_container =  matrix(0, nrow = 1, ncol = N_iter)
      mu_vec_container = matrix(0, nrow = K, ncol = N_iter)
      K_container = matrix(0, nrow = 1, ncol = N_iter)
      P_container = array(0, dim = c(K,K,N_iter))
      
      #if the parameters is fixed, setting it to the true value
      if(estimation_control$z==1){
        z_current= matrix(init$z,N,1)
      }else{
        z_current=  matrix(ground_truth$z, N, 1)
      }
      

      if(estimation_control$sigma_squared==1){
        sigma_squared_current<- as.numeric(init$sigma_squared)
      }else{
        sigma_squared_current=  as.numeric(ground_truth$sigma_squared)
      }
      
      if(estimation_control$mu_vec==1){
        mu_vec_current<- as.numeric(init$mu_vec)
      }else{
        mu_vec_current=  as.numeric(ground_truth$mu_vec)
      }
      
      if(estimation_control$K==1){
        K_current<- as.numeric(init$K)
      }else{
        K_current=  as.numeric(ground_truth$K)
      }
      if(estimation_control$P==1){
        P_current<- as.matrix(init$P0)
      }else{
        P_current=  as.matrix(ground_truth$P)
      }
      
      labels_available<- 1:K_current
      
      #initializing quantities
      n_k = as.vector(table(z_current))
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
      
      
      
      
      A_current = llik_over_blocks_f_binomial( lamdabar = lamdabar, ybar = ybar,mbar = mbar,P = P_current)
      
      z_container[,1] = z_current

      sigma_squared_container[1] = sigma_squared_current
      mu_vec_container[,1] = mu_vec_current
      P_container[,,1]<-P_current
      
      A_container[1]=A_current
      
      #containers for the counts of accepted proposals
      acc.count_z = rep(1,N)

      acc.count_sigma_squared=1
      acc.count_mu_vec = 1
      acc.count_P<- matrix(1,K,K)
      
      #adaptive variances initialisation
      tau_sigma_squared <- 0.2 
      tau_mu_vec <- rep(0.5,(K))
      tau_P <- matrix(0.2,K,K)


      tau_sigma_squared_container<- matrix(0,1, N_iter)
      tau_mu_vec_container<- 1
      tau_P_container <- array(0,dim=c(K,K,N_iter))
      
      
      
      
      
      
      
      
      #READY TO BOMB!
      
      set.seed(seed)
      iteration_time= vector()
      for(j in 2:N_iter){
        start_time <- Sys.time()
        
        # Create a vector of update indices in random order
        
        
        
        if (estimation_control$z == 1) {
          #z UPDATE-------------------------------------------------------------
          
          z_update = z_update_f_withP(z = z_current, N_ij = N_ij, Y_ij = Y_ij,P = P_current,
                                      lamdabar = lamdabar,ybar=ybar,mbar=mbar,alpha_vec = alpha_vec,
                                      n_k = n_k,K = K,
                                      acc.count_z = acc.count_z,labels_available = labels_available)
          
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
                                mu_vec = mu_vec_current,K = K,tau_P = tau_P,acc.count_P = acc.count_P)
          
          P_current = P_update$P
          acc.count_P = P_update$acc.moves
          
          if(j %% 50 == 0 && j< N_iter*0.3){
            
            for(my_p in 1:K){
              for(my_q in my_p:K){
                
                tau_P[my_p,my_q] = tuning_proposal(iteration=j,acceptance_count = acc.count_P[my_p,my_q],
                                                   sigma = tau_P[my_p,my_q],
                                                   acceptanceTarget = optimal_acceptance_rate,
                                                   min_sigma = 0.00002)
                tau_P[K,K]<-.2
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
                                                               acc.count_sigma_squared = acc.count_sigma_squared)
          #updating quantities
          
          acc.count_sigma_squared = sigma_squared_update$acc.moves
          sigma_squared_current = sigma_squared_update$sigma_squared
          if(j %% 50 == 0 && j< N_iter*0.3){
            tau_sigma_squared <- tuning_proposal(iteration=j,
                                                 acceptance_count = acc.count_sigma_squared,
                                                 sigma = tau_sigma_squared,
                                                 acceptanceTarget = optimal_acceptance_rate,
                                                 min_sigma = 0.02)
          }
        }
        if (estimation_control$mu== 1) {
          #P UPDATE----------------------------------------------------------------
          
          mu_update=  mu_update_f_withP(lamdabar = lamdabar, ybar = ybar,mbar = mbar,  P = P_current,
                                        alpha_vec =  alpha_vec, n_k = n_k,
                                      sigma_squared = sigma_squared_current, 
                                      mu_vec = mu_vec_current,K = K, tau_mu_vec = tau_mu_vec,
                                      acc.count_mu_vec)
          #updating quantities
          mu_vec_current = mu_update$mu_vec
          acc.count_mu_vec = mu_update$acc.moves
          
          if(j %% 50 == 0 && j< N_iter*0.3){
              tau_mu_vec <- tuning_proposal(iteration=j,acceptance_count = acc.count_mu_vec,
                                               sigma = tau_mu_vec,
                                               acceptanceTarget = optimal_acceptance_rate,
                                               min_sigma = 0.002)
          }
          
          
        
        }
        
        #storing scales
        tau_sigma_squared_container[j]<- tau_sigma_squared
        tau_mu_vec_container<- tau_mu_vec
        tau_P_container[,,j]<- tau_P
        
        #storing results for inference
        A_container[j] = llik_over_blocks_f_binomial( lamdabar = lamdabar, ybar = ybar,mbar = mbar, P = P_current)
        z_container[,j] <- z_current
        sigma_squared_container[1,j] = sigma_squared_current
        mu_vec_container[,j] <- mu_vec_current
        P_container[,,j] <- P_current
        
        
        
        
        end_time <- Sys.time()
        
        iteration_time<-append(iteration_time,as.numeric(difftime(end_time, start_time, units = "secs")))
        
        if(j%%5000==0){
          avg_iteration<- mean(iteration_time)
          current_time <- Sys.time() # Get the current time

          # Calculate the expected finishing time
          expected_finishing_time <- current_time + (avg_iteration * (N_iter - j) )
          formatted_expected_finishing_time <- format(expected_finishing_time, "%H:%M:%S")

          p(sprintf("Chain %d: mean acc.rate z %.3f%%, mean acc.rate P %.3f%%,mean acc.rate sigma_squared %.3f%%,mean acc.rate mu_vec %.3f%% single_iter_time '%*.3f' seconds, will_finish_at '%s'",
                    chain,
                    100 * mean(acc.count_z/j),
                    100 * mean(acc.count_P/j),
                    100 * mean(acc.count_sigma_squared/j),
                    100 * mean(acc.count_mu_vec/j),
                    4, avg_iteration, formatted_expected_finishing_time), class = "sticky")


        }
      }
      
      acceptance_rates <- list(acc.count_P = acc.count_P, acc.count_z = acc.count_z,
                               acc.count_sigma_squared=acc.count_sigma_squared, acc.count_mu_vec= acc.count_mu_vec)
      
      st.deviations<- list(tau_P = tau_P,tau_sigma_squared = tau_sigma_squared, 
                           tau_mu_vec= tau_mu_vec)
      
      est_containers = list(z = z_container,P = P_container,
                            sigma_squared= sigma_squared_container, mu_vec = mu_vec_container)
      
      control_containers = list(A = A_container)
      
      return(list(Y_ij= Y_ij, N_ij = N_ij,
                  init = init, ground_truth=ground_truth,est_containers=est_containers, 
                  control_containers=control_containers, acceptance_rates= acceptance_rates, 
                  st.deviations=st.deviations, seed=seed))
    }
  })
  
  return(reprex)
} 

