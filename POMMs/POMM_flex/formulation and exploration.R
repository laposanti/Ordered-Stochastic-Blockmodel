
#new functions
library(loo)
library(dbscan)
library(randnet)
library(fossil)
library(dplyr)
library(truncnorm)

source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/functions_container_flex.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/Functions_priorSST.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/SaraWade.R")




K=5
n_samples=1000
S = 0.01
beta_max = .8
alpha=1
diag0.5=T
true_alpha<-alpha



#creating a sample of P matrices
p_container = array(0, dim=c(K,K,n_samples))
for(i in 1:n_samples){
  trunc = improper_prior5(K,alpha = alpha,diag0.5 = diag0.5, beta_max = beta_max )
  p_container[,,i] = simulating_overlapping_POMM_powerlaw_norm(K, alpha=alpha, S  = S, beta_max = beta_max,truncations = trunc,diag0.5 = diag0.5)
}

level_list_p_container<- generalized_levels(p_container,K,n_samples, diag0.5 = diag0.5)
# Combine the four levels into a list

ggplot(combined_data, aes(x = levels, y = values)) +
  geom_point() +
  geom_line(data = mean_values, aes(x = levels, y = values, group=1), color = "red", linewidth = 1.2) +
  labs(x = "Levels", y = "Values") +
  scale_x_discrete(labels = c("1", "2", "3", "4")) +
  stat_summary(fun = mean, geom = "errorbar", color = "red", width = 0.2) +
  stat_summary(fun = mean, geom = "point", color = "red", size = 3, shape = 18)+
  labs(x = "Level sets", y = "P_ij", title = paste0("Distribution of ",n_samples,' simulated P matrices with K=', K-1 , " ,alpha=",alpha,",overlap=",overlap))+
  theme_bw()

paste0('alpha_and_overlap_investigation_alpha',alpha,'overlap',overlap,'K',K)

#interpretation of overlap

overlap_space = vector()
for(i in 2:length(level_list_p_container)){
  overlap_space = append(overlap_space,abs(max(level_list_p_container[[i]]) - min(level_list_p_container[[i-1]])))
}
mean(overlap_space)/(beta_max- 0.5)




blue_purple <- generate_color_gradient_K(K)

#goood plot here
ggplot() +
  # Add a layer for each level
  lapply(seq_along(level_list_p_container), function(i) {
    geom_density(data = data.frame(x = level_list_p_container[[i]]), aes(x = x, y = ..density.., fill = paste0("Level ", i)), alpha = .5)
  }) +
  # Set the x-axis limits
  scale_fill_manual(values=blue_purple)+
  # Set the legend title
  labs(fill = "Levels", x = "Points", y = "Density", title = paste("Density Plot of the ", K, " Level sets [alpha=",alpha,", overlap=",S,", diag=",0.5,"]", sep = ""))+
  theme_bw()


print(paste("Densityplot",K, "Level_sets_alpha_",alpha,"_overlap_",overlap,"_diag_",0.5, sep = "_"))


#########
# Test on the likelihood
#######



alphas_to_be_tested <- c(0.1,0.5,1,1.5,3)
overlap_to_be_tested<- c(0.2,0.4,0.7)
combinations<- expand.grid(alphas_to_be_tested, overlap_to_be_tested)
results<- data.frame( alpha = combinations$Var1, S =combinations$Var2, MAE = rep(0, nrow(combinations)))

for (nrow in 1: nrow(results)){
  K=5
  n_samples=1000
  overlap = results$S[nrow]
  beta_max = .8
  alpha=results$alpha[nrow]
  diag0.5=T
  true_alpha<-alpha
  
  
  #creating a sample of P matrices
  p_container = array(0, dim=c(K,K,n_samples))
  for(i in 1:n_samples){
    trunc = improper_prior5(K,alpha = alpha,diag0.5 = diag0.5, beta_max = beta_max )
    p_container[,,i] = simulating_overlapping_POMM_powerlaw_norm(K, alpha=alpha, overlap  = overlap, beta_max = beta_max,truncations = trunc,diag0.5 = diag0.5)
  }
  
  alpha_test = seq(0.1,4,0.1)
  overlap_test = seq(0.1,.9,0.05)
  #set the containers
  likelihood_est = matrix(0, nrow=n_samples, ncol=length(overlap_test))
  
  
  for(j in 1: n_samples){
    for(i in 1:length(overlap_test)){
      trunc_i = improper_prior5(K,beta_max,alpha = alpha,diag0.5 = diag0.5)
      likelihood_est[j,i]= l_like_p_ij_normal_overlap(K,p_container[,,j],overlap = overlap_test[i],trunc_i,diag0.5 = diag0.5) 
    }
  }
  
  # Define the function to maximize
  likelihood_function <- function(overlap) {
    likelihood_est = matrix(0, nrow=n_samples, ncol=1)
    for(j in 1: n_samples){
      trunc_i = improper_prior5(K,beta_max,alpha,diag0.5 = diag0.5)
      likelihood_est[j]= l_like_p_ij_normal_overlap(K,p_container[,,j],overlap = overlap,trunc_i,diag0.5 = diag0.5)
    }
    return(sum(likelihood_est))
  }
  
  # Find the alpha value that maximizes the likelihood_est
  max_likelihood <- optimize(likelihood_function, interval = c(0.1, 1), maximum = TRUE)
  max_value <- max_likelihood$max
  overlap_true <- overlap
  
  results$MAE[nrow] = abs(max_value - overlap)
  
  df_diagnostic = data.frame(overlap = overlap_test, likelihood= colSums(likelihood_est)) %>% filter(likelihood != - Inf)
  
  print(paste("Overlap_Likelihood",K, "Level_sets_alpha_",alpha,"_overlap_",overlap,"_diag_",0.5, sep = "_"))
  setwd("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/exploratory_graphs_with_inference")
  plot_name <-paste0("Overlap_LikelihoodK",K,"_alpha",alpha,"_overlap",overlap,"_diag_",".png")
  
  lik_plot<- ggplot(df_diagnostic, aes(x = overlap, y = scale(likelihood))) +
    geom_line(colour = blue_purple[3]) +
    geom_vline(aes(xintercept = overlap_true, colour = blue_purple[1]), alpha = 0.8, show.legend = TRUE) +
    geom_vline(aes(xintercept = max_value, colour = "red" ), alpha = 0.8, linetype = "dashed", show.legend = TRUE) +
    theme_bw() +
    labs(x = "Alpha", y = "Likelihood", title = paste0("Likelihood of overlap = [0.1,3], K= ", K, ", true alpha = ", alpha, ", S = ", overlap)) +
    scale_color_manual(values = c("red", blue_purple[4]), labels = c("True Value", "Estimate")) +
    guides(colour = guide_legend(title = "Legend"))
  
  png(plot_name,width = 800, height = 518)
  print(lik_plot)
  dev.off()
}

library(gt)
library(tidyr)
results_summary <- results %>%
  pivot_wider(names_from = alpha, values_from = MAE)%>%
  as.data.frame()

results_summary<-results_summary %>%gt(rowname_col= 'S')%>%
  tab_header(title = md('Mean absolute error for different S values'),
             subtitle = md('Each row is a different S value. Sample size = 1000')) %>%
  cols_label(
    '0.1' =  md("alpha = 0.1"),
    '0.5' =  md("alpha = 0.5"),
    '1' =  md("alpha = 1"),
    '1.5' =  md("alpha = 1.5"),
    '3' =  md("alpha = 3"),
  )%>%
  opt_align_table_header(align = 'left') %>% 
  opt_table_outline()%>%
  fmt_auto()

results_summary

as.character(as_latex(results_summary)[1])%>%writeLines(con = "latex_out.tex")



# Define the blue shades for the color gradient
blue_purple<- generate_color_gradient_K(K)
df_diagnostic = data.frame(alpha = alpha_test, likelihood= colSums(likelihood_est)) %>% filter(likelihood != - Inf)

print(paste("Overlap_Likelihood",K, "Level_sets_alpha_",alpha,"_overlap_",overlap,"_diag_",0.5, sep = "_"))
setwd("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/exploratory_graphs_with_inference")
plot_name <-paste0("Overlap_LikelihoodK",K,"_alpha",alpha,"_overlap",overlap,"_diag_",".png")

lik_plot<- ggplot(df_diagnostic, aes(x = alpha, y = scale(likelihood))) +
  geom_line(colour = blue_purple[3]) +
  geom_vline(aes(xintercept = true_alpha, colour = blue_purple[1]), alpha = 0.8, show.legend = TRUE) +
  geom_vline(aes(xintercept = max_alpha, colour = "red" ), alpha = 0.8, linetype = "dashed", show.legend = TRUE) +
  theme_bw() +
  labs(x = "Alpha", y = "Likelihood", title = paste0("Likelihood of alphas = [0.1,3], K= ", K, ", true alpha = ", alpha, ", S = ", overlap)) +
  scale_color_manual(values = c("red", blue_purple[4]), labels = c("True Value", "Estimate")) +
  guides(colour = guide_legend(title = "Legend"))

png(plot_name,width = 800, height = 518)
print(lik_plot)
dev.off()

print(lik_plot)

abs(max_alpha - true_alpha)


### -----------------

# no plots
# univariate
# Simulation study without plots: recovering overlap, fixing alpha

### ------------------

K=3
n_samples=1000
beta_max = .8
alpha=1
diag0.5=T
true_alpha<-alpha
overlap <- 0.35
true_overlap<-overlap


K_value_test = seq(3,11,2)
alpha_value_test = seq(0.5,3, 0.5)
results=matrix(0,length(K_value_test), length(alpha_value_test))
for(p in 1:length(K_value_test)){
  for(q in 1:length(alpha_value_test)){
    #creating a sample of P matrices
    K = K_value_test[p]
    alpha = alpha_value_test[q]
    p_container = array(0, dim=c(K,K,n_samples))
    for(i in 1:n_samples){
      trunc = improper_prior5(K,alpha = alpha,diag0.5 = diag0.5, beta_max = beta_max )
      p_container[,,i] = simulating_overlapping_POMM_powerlaw_norm(K, alpha=alpha, overlap  = overlap, beta_max = beta_max,truncations = trunc,diag0.5 = diag0.5)
    }
    # Define the function to maximize
    likelihood_function <- function(overlap) {
      likelihood_est = matrix(0, nrow=n_samples, ncol=1)
      for(j in 1: n_samples){
        trunc_i = improper_prior5(K,beta_max,alpha,diag0.5 = diag0.5)
        likelihood_est[j]= l_like_p_ij_normal_overlap(K,p_container[,,j],overlap = overlap,trunc_i,diag0.5 = diag0.5) 
      }
      return(sum(likelihood_est))
    }
    # Find the alpha value that maximizes the likelihood_est
    max_likelihood <- optimize(likelihood_function, interval = c(0.1, 1), maximum = TRUE)
    max_overlap <- max_likelihood$max
    
    results[p,q] <- abs(max_overlap - overlap)
    
  }
}


### -----------------
# no plots
# univariate
# Simulation study without plots: recovering alpha, fixing overlap
### -----------------

K=3
n_samples=1000
overlap = .5
beta_max = .8
alpha=1
diag0.5=T
true_alpha<-alpha

K_value_test = seq(3,11,2)
overlap_value_test = seq(0.1,0.9, 0.1)
results=matrix(0,length(K_value_test), length(overlap_value_test))
for(p in 1:length(K_value_test)){
  for(q in 1:length(overlap_value_test)){
    #creating a sample of P matrices
    K = K_value_test[p]
    overlap = overlap_value_test[q]
    p_container = array(0, dim=c(K,K,n_samples))
    for(i in 1:n_samples){
      trunc = improper_prior5(K,alpha = alpha,diag0.5 = diag0.5, beta_max = beta_max )
      p_container[,,i] = simulating_overlapping_POMM_powerlaw_norm(K, alpha=alpha, overlap  = overlap, beta_max = beta_max,truncations = trunc,diag0.5 = diag0.5)
    }
    # Define the function to maximize
    likelihood_function <- function(alpha) {
      likelihood_est = matrix(0, nrow=n_samples, ncol=1)
      for(j in 1: n_samples){
        trunc_i = improper_prior5(K,beta_max,alpha,diag0.5 = diag0.5)
        likelihood_est[j]= l_like_p_ij_normal_overlap(K,p_container[,,j],overlap = overlap,trunc_i,diag0.5 = diag0.5) + dlnorm_param( alpha)
      }
      return(sum(likelihood_est))
    }
    # Find the alpha value that maximizes the likelihood_est
    max_likelihood <- optimize(likelihood_function, interval = c(0.1, 3), maximum = TRUE)
    max_alpha <- max_likelihood$max
    
    results[p,q] <- abs(max_alpha - true_alpha)
    
  }
}


### -----------------
K = 3
n_samples = 10000
beta_max = 0.8
diag0.5 = TRUE
alpha=1
overlap=0.2
K_value_test = seq(3, 11, 2)
results = matrix(0, length(K_value_test), 2)

for (p in 1:length(K_value_test)) {
  K = K_value_test[p]
  p_container = array(0, dim = c(K, K, n_samples))
  
  for (i in 1:n_samples) {
    trunc = improper_prior5(K, alpha = alpha, diag0.5 = diag0.5, beta_max = beta_max)
    p_container[,,i] = simulating_overlapping_POMM_powerlaw_norm(K, alpha = alpha, overlap = overlap, beta_max = beta_max, truncations = trunc, diag0.5 = diag0.5)
  }
  
  # Define the likelihood function to maximize
  likelihood_function <- function(params) {
    alpha = params[1]
    overlap = params[2]
    
    likelihood_est = matrix(0, nrow = n_samples, ncol = 1)
    
    for (jjj in 1:n_samples) {
      trunc_i = improper_prior5(K, beta_max, alpha, diag0.5 = diag0.5)
      likelihood_est[jjj] = l_like_p_ij_normal_overlap(K, p_container[,,jjj], overlap = overlap, trunc_i, diag0.5 = diag0.5) + dlnorm_param(alpha)
    }
    
    likelihood_est[likelihood_est==-Inf] <- -2**31
    
    
    return(sum(likelihood_est))
  }
  
  # lik_explore= matrix(0,nrow = nrow(param_test), ncol=1)
  # param_test <- as.matrix(expand_grid(overlap_value_test,alpha_value_test))
  # param_test = matrix(c(0.31,0.1),nrow = 1,ncol=2)
  # for(k in 1:nrow(param_test)){
  #   lik_explore[k]<- likelihood_function(param_test[k,])
  # }
  
  
  
  # Find the alpha and overlap values that maximize the likelihood_est
  max_likelihood <- optim(c(1, 0.35), likelihood_function, method = "L-BFGS-B", lower = c(0.3, 0.1), upper = c(2.7, .9))
  max_alpha <- max_likelihood$par[1]
  max_overlap <- max_likelihood$par[2]
  
  results[p,1] <- abs(max_alpha - alpha) 
  results[p,2] <- abs(max_overlap - overlap)
}

#
# MCMC on P_ij | alpha and overlap
#

N_iter = 10000
overlap = 0.2
alpha=0.5
beta_max=0.8
K=5
truncations= improper_prior5(K=K,beta_max = beta_max,alpha = alpha,diag0.5 = T)
P_data = simulating_overlapping_POMM_powerlaw_norm(K=5,alpha = alpha,overlap = overlap,truncations = truncations,beta_max = beta_max,diag0.5 = T)

p_container = array(0, dim=c(K,K,N_iter))
acc.count_p = matrix(0, K,K)
truncations_current = truncations
overlap_current=overlap
alpha_current=alpha
p_current = simulating_overlapping_POMM_powerlaw_norm(K=5,alpha = 0.8,overlap = 0.5,truncations = truncations,beta_max = beta_max,diag0.5 = T)
C_current <- l_like_p_ij_normal_overlap(K = K, P_matrix = p_current,overlap = 
                                          overlap_current, 
                                        truncations = truncations_current,
                                        diag0.5 = T) +  dlnorm_mu_sigma(alpha_current) + dlnorm_mu_sigma(overlap_current,mu = 0.3,sigma   = 0.5)         



j_start = ifelse(diag0.5, yes = 1, no = 0)
K_stop = ifelse(diag0.5, yes = K-1, no = K)

for( iter in 1:N_iter){
  for( ii in 1:K_stop){
    for(jj in (ii+j_start):K){
      
      
      p_prime = p_current
      p_prime[ii,jj] <- rnormTrunc(1, p_current[ii,jj],sd = 0.02, min = 0.5, max = beta_max)
      p_prime[jj,ii] <- 1 - p_prime[ii,jj]
      
      C_prime <- l_like_p_ij_normal_overlap(K = K, P_matrix = p_prime,overlap = 
                                              overlap_current, 
                                            truncations = truncations_current,
                                            diag0.5 = T) + + dlnorm_mu_sigma(alpha_current) + dlnorm_mu_sigma(overlap_current,mu = 0.3,sigma   = 0.5)         
      
      log_r=  C_prime - C_current
      
      #create statements that check conditiond to accept move
      MH_condition= min(log_r,0)>=log(runif(1))
      if(MH_condition){
        acc.count_p[ii,jj] =acc.count_p[ii,jj] +1
        acc.count_p[jj,ii] =acc.count_p[jj,ii] +1
        C_current<- C_prime
        p_current<- p_prime
      }
    }}
  
  p_container[,,iter]<- p_current
  
}

burnin_p = p_container[,,-(N_iter*0.25)]
plots = list()
for(i in 1:K) {
  for(j in 1:K) {
    y_try = data.frame(y = as.vector(burnin_p[i, j,]))
    p1 = ggplot(y_try, aes(y)) +
      geom_density(fill = "dodgerblue", alpha = 0.5) +
      scale_x_log10() +
      geom_vline(xintercept = P_data[i, j], color = "red")+
      xlab("probability") +
      ylab("Density") +
      ggtitle(paste("Density plot of entry ", i, ",", j, sep = ""))
    
    plots[[length(plots) + 1]] <- p1
  }
}
p_combined = patchwork::wrap_plots(plots, ncol = K, nrow = K)
p_combined


#-----
# MCMC on alpha | P_ij and overlap
#-----

N_iter = 10000
overlap = 0.2
alpha=0.1
beta_max=0.8
K=5
truncations= improper_prior5(K=K,beta_max = beta_max,alpha = alpha,diag0.5 = T)
P_data = simulating_overlapping_POMM_powerlaw_norm(K=K,alpha = alpha,overlap = overlap,truncations = truncations,beta_max = beta_max,diag0.5 = T)

alpha_container = matrix(0, N_iter, 1)
acc.count_p = matrix(0, K,K)

overlap_current=overlap
alpha_current=0.8

truncations_current = improper_prior5(K=K,beta_max = beta_max,alpha = alpha_current,diag0.5 = T)
p_current = P_data
C_current <- l_like_p_ij_normal_overlap(K = K, P_matrix = p_current,overlap = 
                                          overlap_current, 
                                        truncations = truncations_current,
                                        diag0.5 = T) +  dlnorm_mu_sigma(alpha_current) + dlnorm_mu_sigma(overlap_current,mu = 0.3,sigma   = 0.5)         


acc.count_alpha = 0
j_start = ifelse(diag0.5, yes = 1, no = 0)
K_stop = ifelse(diag0.5, yes = K-1, no = K)

for( iter in 1:N_iter){
  
  #proposing a new overlap
  
  alpha_prime <- rtruncnorm(1,a = 0.1,b = 0.9,mean = alpha_current,sd = 0.02)
  truncations_prime <- improper_prior5(K = K,beta_max = beta_max,alpha = alpha_prime,diag0.5 = T)
  
  
  C_prime <- l_like_p_ij_normal_overlap(K = K, P_matrix = p_current,overlap = overlap_current, truncations = truncations_prime,diag0.5 = T) +  
    dlnorm_mu_sigma(alpha_prime,mu = 0.5, sigma = 1,log = T) + dlnorm_mu_sigma(overlap_current,mu = 0.3,sigma  = 0.5,log = T)      
  
  
  
  
  #A_prime <- sum(dbinom(y_ij, n_ij, p_ij_prime, log = T))
  
  log_r= C_prime - C_current
  
  #create statements that check conditiond to accept move
  MH_condition_alpha= min(log_r,0)>=log(runif(1))
  if(MH_condition_alpha){
    acc.count_alpha= acc.count_alpha+1
    alpha_current <- alpha_prime
    truncations_current<- truncations_prime
    C_current<- C_prime
    #A_current <- A_prime
    
    
  }
  
  alpha_container[iter]<- alpha_current
  
}

ts.plot(alpha_container)

mean(alpha_container[which(alpha_container!= 0)])



#
# MCMC on overlap | P_ij and alpha
#

N_iter = 30000
overlap = 0.7
alpha=.5
beta_max=0.8
K=5

truncations= improper_prior5(K=K,beta_max = beta_max,alpha = alpha,diag0.5 = T)
P_data = simulating_overlapping_POMM_powerlaw_norm(K=K,alpha = alpha,overlap = overlap,truncations = truncations,beta_max = beta_max,diag0.5 = T)

overlap_container = matrix(0, N_iter, 1)
acc.count_overlap = 0
alpha_current=alpha

truncations_current = improper_prior5(K=K,beta_max = beta_max,alpha = alpha_current,diag0.5 = T)

p_current = P_data
C_current <- l_like_p_ij_normal_overlap(K = K, P_matrix = p_current,overlap = 
                                          overlap_current, 
                                        truncations = truncations_current,
                                        diag0.5 = T) 


for( iter in 1:N_iter){
  
  #proposing a new overlap
  
  overlap_prime <- rtruncnorm(1,a = 0.1,b = 0.95,mean = overlap_current,sd = 0.02)
  
  
  
  C_prime <- l_like_p_ij_normal_overlap(K = K, P_matrix = p_current,overlap = 
                                          overlap_prime, 
                                        truncations = truncations_current,
                                        diag0.5 = T) 
  
  
  
  log_r=  C_prime -  C_current
  
  #create statements that check conditiond to accept move
  MH_condition_overlap= min(log_r,0)>=log(runif(1))
  if(MH_condition_overlap){
    acc.count_overlap=acc.count_overlap+1
    overlap_current <- overlap_prime
    C_current<- C_prime
  }
  
  
  overlap_container[iter]<- overlap_current
  
}

ts.plot(overlap_container[-c(1:N_iter*0.25)])
acf(overlap_container[-c(1:N_iter*0.25)])
mean(overlap_container[-c(1:N_iter*0.25)])


###
#bivariate analysis
###


K = 3
n_samples = 1000
beta_max = 0.8
diag0.5 = TRUE
#true values
alpha= 0.5
overlap=0.2
#values combinations
alpha_value_test = seq(0.1,3,0.1)
overlap_value_test= seq(0.1,.8,0.1)

results <- expand.grid(alpha_value_test, overlap_value_test) %>% mutate(likelihood = rep(0,nrow(params)))




p_container = array(0, dim = c(K, K, n_samples))

for (i in 1:n_samples) {
  trunc = improper_prior5(K, alpha = alpha, diag0.5 = diag0.5, beta_max = beta_max)
  p_container[,,i] = simulating_overlapping_POMM_powerlaw_norm(K, alpha = alpha, overlap = overlap, beta_max = beta_max, truncations = trunc, diag0.5 = diag0.5)
}


for(rows in 1:nrow(params)){
  
  
  alpha = results$Var1[rows]
  overlap = results$Var2[rows]
  
  likelihood_est = matrix(0, nrow = n_samples, ncol = 1)
  
  for (jjj in 1:n_samples) {
    trunc_i = improper_prior5(K, beta_max, alpha, diag0.5 = diag0.5)
    likelihood_est[jjj] = l_like_p_ij_normal_overlap(K, p_container[,,jjj], overlap = overlap, trunc_i, diag0.5 = diag0.5) }
  
  #likelihood_est[likelihood_est==-Inf] <- -2**31
  
  
  results$likelihood[rows]<- sum(likelihood_est)
  
}


head(results)


# Reshape results for plotting
results_df <- data.frame(alpha<- params[,1], overlap<-params[,2], likelihood = results) %>% filter(likelihood>0) %>% mutate(likelihood= likelihood**3) %>%  mutate(likelihood= scale(likelihood))


colnames(results_df) <- c("alpha", "overlap", "likelihood")
# Create the bivariate plot
ggplot(results_df, aes(x = alpha, y = overlap, fill = likelihood)) +
  geom_tile() +
  labs(x = "Alpha", y = "overlap", fill = "Likelihood") +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal()


K=3
n_samples=1000
overlap = .5
beta_max = .8
alpha=1
diag0.5=T
true_alpha<-alpha


#creating a sample of P matrices
p_container = array(0, dim=c(K,K,n_samples))
for(i in 1:n_samples){
  trunc = improper_prior5(K,alpha = alpha,diag0.5 = diag0.5, beta_max = beta_max )
  p_container[,,i] = simulating_overlapping_POMM_powerlaw_norm(K, alpha=alpha, overlap  = overlap, beta_max = beta_max,truncations = trunc,diag0.5 = diag0.5)
}


# Combine the four levels into a list
level_list_p_container <- generalized_levels(p_container,K,N = n_samples,diag0.5 = diag0.5)




blue_purple <-generate_color_gradient(K)

ggplot() +
  # Add a layer for each level
  lapply(seq_along(level_list_p_container), function(i) {
    geom_density(data = data.frame(x = level_list_p_container[[i]]), aes(x = x, y = ..density.., fill = paste0("Level ", i)), alpha = .5)
  }) +
  # Set the x-axis limits
  scale_fill_manual(values=blue_purple)+
  # Set the legend title
  labs(fill = "Levels", x = "Points", y = "Density", title = paste("Density Plot of the ", K, " Level sets [alpha=",alpha,", overlap=",overlap,", diag=",0.5,"]", sep = ""))+
  theme_bw()


print(paste("Densityplot",K, "Level_sets_alpha_",alpha,"_overlap_",overlap,"_diag_",0.5, sep = "_"))

alpha_test = seq(0.1,3,0.1)

#set the containers
likelihood_est = matrix(0, nrow=n_samples, ncol=length(alpha_test))


for(j in 1: n_samples){
  for(i in 1:length(alpha_test)){
    trunc_i = improper_prior5(K,beta_max,alpha = alpha_test[i],diag0.5 = diag0.5)
    likelihood_est[j,i]= l_like_p_ij_normal_overlap(K,p_container[,,j],overlap = overlap,trunc_i,diag0.5 = diag0.5) + dlnorm_param( alpha_test[i])
  }
}

# Define the function to maximize
likelihood_function <- function(alpha) {
  likelihood_est = matrix(0, nrow=n_samples, ncol=1)
  for(j in 1: n_samples){
    trunc_i = improper_prior5(K,beta_max,alpha,diag0.5 = diag0.5)
    likelihood_est[j]= l_like_p_ij_normal_overlap(K,p_container[,,j],overlap = overlap,trunc_i,diag0.5 = diag0.5) + dlnorm_param( alpha)
  }
  return(sum(likelihood_est))
}

# Find the alpha value that maximizes the likelihood_est
max_likelihood <- optimize(likelihood_function, interval = c(0.1, 3), maximum = TRUE)
max_alpha <- max_likelihood$max

max_alpha

# Define the blue shades for the color gradient

df_diagnostic = data.frame(alpha = alpha_test, likelihood= colSums(likelihood_est)) %>% filter(likelihood != - Inf)

ggplot(df_diagnostic, aes(x = alpha, y = scale(likelihood))) +
  geom_line(colour = blue_purple[3]) +
  geom_vline(aes(xintercept = true_alpha, colour = blue_purple[1]), alpha = 0.8, show.legend = TRUE) +
  geom_vline(aes(xintercept = max_alpha, colour = "red" ), alpha = 0.8, linetype = "dashed", show.legend = TRUE) +
  theme_bw() +
  labs(x = "Alpha", y = "Likelihood", title = paste("Likelihood of alphas = [0.1,3] | ", K, " Level sets | alpha =", alpha, " | overlap =", overlap, " | diag =", 0.5, ".", sep = "")) +
  scale_color_manual(values = c("red", blue_purple[4]), labels = c("True Value", "Estimate")) +
  guides(colour = guide_legend(title = "Legend"))

print(paste("Likelihood",K, "Level_sets_alpha_",alpha,"_overlap_",overlap,"_diag_",0.5, sep = "_"))


abs(max_alpha - true_alpha)


###








