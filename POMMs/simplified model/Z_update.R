
library(dplyr)
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/Functions_priorSST.R")
setwd("/Users/lapo_santi/Desktop/Nial/project/simplified model")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/SaraWade.R")


###########
###INPUT
#N: number of players
#K:number of blocks
#z: label controlling the assignment of each player to one block
#p_ij: probabilities of victory against different blocks
#n_ij:number of matches between player ij
#y_ij:number of victories of player i against player j
#consider that p_ji ==1 - p_ij
#consider also that y_ji = n_ij - y_ij



N=100
a=1
b=3
K_max=4
M=300
max_clust=3
min_clust=3
max_number_games = 200
N_iter = 100000
############

############
#SIMULATING THE TOURNAMENT
synth = simulating_tournament_test(N=N,alpha = 1,beta = 3,min_clust = min_clust,max_clust = max_clust , M = M, n_ij_max =max_number_games )
synth_matches = synth$matches_results
synth_players = synth$z_true
synth_p = synth$p_true
K_true = synth$K_true

test = df_aux_fast(synth_matches,synth_players,synth_p)
n_test = table(synth_players$z)
n_test
cor(test$y_ij/test$n_ij, test$p_ij)
###################

#########
#INITIALIZING THE MCMC
##########


labels_available = 1:K_true
z_proposal = sample(labels_available, 1)


z_current = data.frame(id = synth_players$id, z = synth_players$z)
df_current = df_aux_fast(synth_matches, z_current, synth_p)
n_current = table(z_current$z)
A_seq = matrix(0, N_iter, 1)
z_seq = matrix(0, N, N_iter)
z_seq[, 1] =  z_current$z
A_seq[1] = get_A(df_current$n_ij, df_current$y_ij, df_current$p_ij) + dir_multinom_d(N, K_true, n_current)
acc.count = 0
pb = txtProgressBar(min = 1, max = N_iter)
j = 2
max_reassign = 3  # maximum number of nodes to reassign
min_acc_ratio = 0.3  # minimum acceptance ratio
reassign_nodes = 2  # initial number of nodes to reassign

#######
##METROPOLIS-HASTINGS

while (j < N_iter + 1) {
  setTxtProgressBar(pb, j)
  z_prime = z_current
  df_prime = df_current
  n_prime = n_current
  for (ii in 1:N) {
    reassign_nodes <- sample(size = 1, x=c(1:max_reassign), p=c(0.5,seq(from = 0.3,to= 0.01,by= ((0.01 - 0.3 )/(max_reassign - 2)))))
    if (reassign_nodes == 1) {
      z_scanning = get_proposal1(z_prime, labels_available)
    } else {
      z_scanning = get_proposal2(z_prime, labels_available, reassign_nodes)
    }
    df_scanning = df_aux_fast(synth_matches, z_scanning, synth_p)
    n_scanning = table(factor(z_scanning$z, levels = labels_available))
    r = (get_A(df_scanning$n_ij, df_scanning$y_ij, df_scanning$p_ij) + dir_multinom_d(N, labels_available, n_scanning)) -
      (get_A(df_prime$n_ij, df_prime$y_ij, df_prime$p_ij) + dir_multinom_d(N, labels_available, n_prime))
    u = runif(1)
    if (log(u) < r) {
      acc.count = acc.count + 1
      z_prime = z_scanning
      df_prime = df_scanning
      n_prime = n_scanning
      A_seq[j] = get_A(df_prime$n_ij, df_prime$y_ij, df_prime$p_ij) + dir_multinom_d(N, labels_available, n_prime)
    } else {
      A_seq[j] = A_seq[j - 1]
      # if (acc.count / j < min_acc_ratio && reassign_nodes < max_reassign) {
      #   reassign_nodes = reassign_nodes + 1
      # }
    }
    z_current = z_prime
    df_current = df_prime
    n_current = n_prime
    z_seq[, j - 1] = z_prime$z
    j = j + 1
  }
}

##########


#######
#RESULTS
close(pb)
cat("accepted ", signif(100*acc.count/(N_iter*N),2), "%\n",sep="")
MAP = which(A_seq[1:(N_iter-1)] == max(A_seq[1:(N_iter-1)]))


ts.plot(A_seq[50000:N_iter])



similarity<- pr_cc(z_seq[,c(50000:N_iter)])

memb_Z_DP_VI <- minVI(similarity,method="avg",max.k=5)

#POINT ESTIMATES
memb_Z_DP_VI$cl
z_seq[,MAP[1]]


similarity_plot(similarity,z_0 = synth_players$z, z_est =z_seq[,MAP[1]])






############