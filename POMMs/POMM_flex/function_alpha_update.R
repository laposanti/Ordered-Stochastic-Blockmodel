

P_POMM_alpha_update = function(z_current, p_current, 
                               A_current,C_current,y_ij,n_ij,labels_available,
                               upper.tri.non.zero,K,alpha_current,beta_max, overlap_current, acc.count_alpha, sigma_alpha, truncations_current){
  
  z_current_mat<- vec2mat(z_current)
  
#proposing a new overlap

alpha_prime <- rtruncnorm(1,a = 0.1,b = 0.9,mean = overlap_current,sd = sigma_alpha)
truncations_prime <- improper_prior5(K,beta_max,alpha = alpha_prime,diag0.5 = T)

# #proposing a new P
# p_prime <-  simulating_overlapping_POMM_powerlaw_norm(K, alpha_current, 
#                                                       truncations_current, 
#                                                       overlap = overlap_prime, 
#                                                       beta_max = beta_max, diag0.5 = T)
# 
# #computing full probabilities
# p_ij_prime_nbyn <- calculate_victory_probabilities(z_current_mat, P = p_prime)
# p_ij_prime <- p_ij_prime_nbyn[upper.tri.non.zero]


C_prime <- l_like_p_ij_normal_overlap(K = K, P_matrix = p_current,overlap = 
                                        overlap_current, 
                                      truncations = truncations_prime,
                                      diag0.5 = T) + dlnorm_param(alpha_prime) + dlnormTruncAlt(x = overlap_current,mean =0.5, cv = 1, min = 0.1,max=0.9)       




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


return(list(acc.moves = acc.count_alpha,
            sigma_alpha=sigma_alpha,
            C_current = C_current, 
            alpha_current = alpha_current,
            truncations_current=truncations_current))

} 


