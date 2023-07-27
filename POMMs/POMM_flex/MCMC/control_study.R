

source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/functions_container_flex.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/Functions_priorSST.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/SaraWade.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/MCMC/adaptive_POMM_MCMC_function.R")

library(EnvStats)
library(ggplot2)
library(dplyr)
library(truncnorm)
library(fossil)
M= 13
N = 30
N_iter =5000
N_ij = matrix(M,N,N)
alpha=.5
S=.05
K=3
beta_max=0.75
gamma_vec = rep(1/K,K)
diag0.5=T
trunc = improper_prior5(K,beta_max ,alpha,diag0.5)
P=simulating_overlapping_POMM_powerlaw_norm(K,alpha,S,trunc,beta_max,diag0.5)

z = rep_len(1:K, N)

z_mat=vec2mat(z)
P_NbyN<-calculate_victory_probabilities(z_mat,P)

Y_ij = matrix(0,N,N)
for(i in 1:N){
  for(j in 1:N){
    Y_ij[i,j] = rbinom(1,M, P_NbyN[i,j])
  }
}


z_k = which(z == 1)
P_NbyN[z_k,2]
mean(Y_ij[z_k,2]/M)


#chains= list()
#for(i in 1:4){
alpha=runif(1,0.1,3)
trunc=improper_prior5(K,beta_max,alpha = alpha)
S=runif(1,0.1,.9)
P_start= simulating_overlapping_POMM_powerlaw_norm(K,alpha,S,trunc,beta_max,diag0.5)
init =list(z = rep_len(sample(1:K,K,F), N),alpha=alpha,S=S,P=P_start)
names(init)
estimation_control = list(z = 0,alpha=0,S=0,P=1)
ground_truth= list(z = z,alpha=alpha,S=S,P=P)
hyper_params = list(K = K,beta_max =beta_max,gamma_vec = gamma_vec,diag0.5=diag0.5)


TEST = adaptive_MCMC_POMM(Yij_matrix = Y_ij,Nij_matrix = N_ij,init = init,estimation_control = estimation_control,ground_truth = ground_truth,N = N,N_iter = N_iter,targ_rate = .22,hyper_params =hyper_params ,seed = 123)

TEST$acceptance_rates$acc.count_p

mean(TEST$st.deviations$sd_p[1,2,])
mean(TEST$st.deviations$sd_p[1,3,])
mean(TEST$st.deviations$sd_p[2,3,])
for(i in 1:K){
  for(j in 1:K)
    print(mean(TEST$est_containers$P[i,j,]) - P[i,j])
}

library(gt)
library(coda)

#creating a table with entry, mean, HPD 95%, if true value is provided also the true value 

P_summary_table<- function(MCMC_samples, true_value, diag0.5,P,K){
  j_start = ifelse(diag0.5, yes = 1, no = 0)
  K_stop = ifelse(diag0.5, yes = K-1, no = K)
  
  entries_df <- data.frame(entry_i = 0 ,entry_j =0 )
  for( ii in 1:K_stop){
    for(jj in (ii+j_start):K){
      entries_df <- rbind(entries_df, data.frame(entry_i= ii, entry_j = jj))
    }
  }
  entries_df=entries_df[-1,]   
  
  if(true_value == F){
    results = cbind(entries_df, data.frame(mean_est = rep(0,nrow(entries_df)),
                                           credible_interval_95 =rep(0,nrow(entries_df))))
    for(i in 1:nrow(results)){
      m<-mcmc(MCMC_samples[results$entry_i[i],results$entry_j[i],])
      results$mean_est[i] <- mean(m)
      HPD <- round(cbind(coda::HPDinterval(m)),2)
      results$credible_interval_95[i]<- paste0("[",HPD[1],",",HPD[2],"]")
    }
  }else if(true_value == T){
    results = cbind(entries_df, data.frame(mean_est = rep(0,nrow(entries_df)),
                                           credible_interval_95 =rep(0,nrow(entries_df)), true_value =rep(0,nrow(entries_df))))
    for(i in 1:nrow(results)){
      m<-mcmc(MCMC_samples[results$entry_i[i],results$entry_j[i],])
      results$mean_est[i] <- round(mean(m),4)
      HPD <- round(cbind(coda::HPDinterval(m)),4)
      results$credible_interval_95[i]<- paste0("[",HPD[1],",",HPD[2],"]")
      results$true_value[i]<- P[results$entry_i[i],results$entry_j[i]]
    }
  }
  return(results)}
x=10
while(min(which(abs(coda::autocorr.diag(m,1:x)) <0.01))==Inf){
  x=x+1
  print(x)
  if(x>1000)
    stop
}
min(which(abs(coda::autocorr.diag(m,1:10)) <0.01))


P_summary_table(MCMC_samples = TEST$est_containers$P,true_value = T,diag0.5 = T,P = P,K=K)

P_summary_table<- function(MCMC_samples, true_value, diag0.5,P,K){
  j_start = ifelse(diag0.5, yes = 1, no = 0)
  K_stop = ifelse(diag0.5, yes = K-1, no = K)
  
  entries_df <- data.frame(entry_i = 0 ,entry_j =0 )
  for( ii in 1:K_stop){
    for(jj in (ii+j_start):K){
      entries_df <- rbind(entries_df, data.frame(entry_i= ii, entry_j = jj))
    }
  }
  entries_df=entries_df[-1,]   
  
  if(true_value == F){
    results = cbind(entries_df, data.frame(mean_est = rep(0,nrow(entries_df)),
                                           credible_interval_95 =rep(0,nrow(entries_df))))
    for(i in 1:nrow(results)){
      m<-mcmc(MCMC_samples[results$entry_i[i],results$entry_j[i],])
      results$mean_est[i] <- mean(m)
      HPD <- round(cbind(coda::HPDinterval(m)),2)
      results$credible_interval_95[i]<- paste0("[",HPD[1],",",HPD[2],"]")
    }
  }else if(true_value == T){
    results = cbind(entries_df, data.frame(mean_est = rep(0,nrow(entries_df)),
                                           credible_interval_95 =rep(0,nrow(entries_df)), true_value =rep(0,nrow(entries_df))))
    for(i in 1:nrow(results)){
      m<-mcmc(MCMC_samples[results$entry_i[i],results$entry_j[i],])
      results$mean_est[i] <- round(mean(m),4)
      HPD <- round(cbind(coda::HPDinterval(m)),4)
      results$credible_interval_95[i]<- paste0("[",HPD[1],",",HPD[2],"]")
      results$true_value[i]<- P[results$entry_i[i],results$entry_j[i]]
    }
  }
  return(results)}




for(i in 1:nrow(results)){
  m<-mcmc(TEST$est_containers$P[results$entry_i[i],results$entry_j[i],])
  results$mean_est[i] <- mean(m)
  HPD <- round(cbind(coda::HPDinterval(m)),2)
  results$credible_interval_95[i]<- paste0("[",HPD[1],",",HPD[2],"]")
}

m<-mcmc(TEST$est_containers$P[1,2,])
coda::HPDinterval(m)

mean(m)






#chains[[i]]= TEST
#}

test1<-chains[[1]]
test2<-chains[[2]]
test3<-chains[[3]]
test4<-chains[[4]]

mcmc_list<- mcmc.list(mcmc(t(test1$est_containers$S)), mcmc(t(test2$est_containers$S)),
                      mcmc(t(test3$est_containers$S)),mcmc(t(test4$est_containers$S)))

#computing gelman rubin diagnostics
unlist(gelman.diag(mcmc_list))[1]

# mean effective sample size
mean_acceptance_rate<- mean(test1$acceptance_rates$acc.count_S,test2$acceptance_rates$acc.count_S,
                            test3$acceptance_rates$acc.count_S,test4$acceptance_rates$acc.count_S)/N_iter

mean_effective_size<-mean(unlist(lapply(mcmc_list, effectiveSize)))

#highest posterior density interval
lapply(mcmc_list, HPDinterval)[[1]][1]
coda::HPDinterval(mcmc_list)

S_diagnostic<- data.frame(gelman_diagnostic = unlist(gelman.diag(mcmc_list))[1], 
                          effective_sample_size = mean_effective_size, mean_acceptance_rate=mean_acceptance_rate*100,
                          HPD_region_0.05 = lapply(mcmc_list, HPDinterval)[[1]][1],
                          HPD_region_0.95 = lapply(mcmc_list, HPDinterval)[[1]][2])



TEST$est_containers$z[,which(TEST$control_containers$A == max(TEST$control_containers$A))[1]]






plot(TEST$st.deviations$sd_alpha)
plot(density(TEST$est_containers$alpha[-c(1:2000)]))

#extracting similarity matrix
similarity_matrixPOMMM = pr_cc(TEST$est_containers$z[,-c(1:1000)])
similarity_plot(similarity_matrixPOMMM, z, z) #checking mixing
library(coda)

#point est 1
point_est_POMM = minVI(similarity_matrixPOMMM,method = 'avg')$cl



data <- data.frame(alpha = alpha_container[-c(1:6000)], S = S_container[-c(1:6000)])
alpha_true=alpha
S_true<- S
# Create scatterplot with univariate distributions
ggplot(data, aes(x = alpha, y = S)) +
  geom_point(alpha=.8) +
  geom_point(x=alpha,y=S, col='red',size=2,shape=15)+
  geom_density_2d(size = .8, alpha=.7) +
  geom_xsidedensity(fill = 'red', alpha = 0.5) +
  geom_ysidedensity(fill = 'blue', alpha = 0.5) +
  labs(title = "Bivariate Scatterplot with Densities", x = "Alpha", y = "S") +
  theme_bw()


library(ggplot2)
library(ggside)



# Create scatterplot with univariate distributions
ggplot(data, aes(x = alpha, y = S)) +
  geom_point() +
  geom_xsidedensity(aes(fill = alpha), alpha = 0.5, position = "stack") +
  geom_ysidedensity(aes(fill = S), alpha = 0.5, position = "stack") +
  theme_minimal()



