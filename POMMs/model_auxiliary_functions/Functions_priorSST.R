

########
##DIRICHLET-MULTINOMIAL
#########


rdirichlet_multinomial <- function(N, K, alpha_vec) {
  # Step 1: Draw theta from Dirichlet distribution
  theta <- MCMCpack::rdirichlet(1, alpha_vec)
  
  # Step 2: Draw n_k from Multinomial distribution for each k
  n <- rmultinom(N, 1, theta)
  
  #Step 3: 
  return(n)
}



ddirichlet_multinomial<- function(N, K, n_k, my_alpha) {
  A <- sum(my_alpha)
  constant_term <- lgamma(A) + lgamma(N+1) - lgamma(N+A)
  variable_term <- sum(lgamma(n_k + my_alpha) - lgamma(my_alpha) - lgamma(n_k+1))
  return(constant_term + variable_term)
}

########
##POISSON
#########


sample_truncated_poisson = function(N,lambda, a ,b){
  return(qpois(((ppois(b, lambda) -ppois(a, lambda))*(runif(N))+ppois(a, lambda)),lambda))}

########
##BETA
#########
sample_beta_trunc = function(N,alpha,beta,a,b){
  u = runif(N)
  return( qbeta(((pbeta(b, alpha, beta) -pbeta(a, alpha, beta))*u+pbeta(a, alpha, beta)),alpha,beta))}



pdf_beta_trunc = function(y,alpha,beta,a,b){
  pdf= (dbeta(y,alpha,beta)*(y>a)*(y<=b))/(pbeta(b,alpha,beta)-pbeta(a,alpha,beta))
  return(pdf)}



cdf_beta_trunc = function(y,alpha,beta,a,b){
  pdf= (pbeta(y,alpha,beta)*(y>a)*(y<=b))/(pbeta(b,alpha,beta)-pbeta(a,alpha,beta))
  return(pdf)}

log_lik_dtruncbeta = function(y,alpha,beta,a,b, my_log=T){
  p = pdf_beta_trunc(y,alpha,beta,a,b)
  if(my_log==T)
  {return(sum(log(p)))}
  else{return(prod(p))}
}


log_lik_dtruncbeta_prime = function(y,alpha,beta,a,b, my_log=T){
  y = as.matrix(y)
  y = y[upper.tri(y)]
  p = pdf_beta_trunc(y,alpha,beta,a,b)
  if(my_log==T)
  {return(sum(log(p)))}
  else{return(prod(p))}
}


sampling_SST_matrix_beta = function(k,alpha,beta){
  Y = matrix(NA, k, k)
  diag(Y) = 0.5
  
  for(i in 2:k) {
    Y[1, i] = sample_beta_trunc(N=1,  alpha, beta, a = Y[1, i-1], b = Inf)
  }
  
  for(i in 2:(k-1)) {
    for(j in (i+1):k) {
      Y[i, j] = sample_beta_trunc(N=1, alpha, beta, a = Y[i, j-1], b = Y[i-1, j])
    }
  }
  return(Y)
}

########
##EXP
#########
exp_trunc = function(N, lambda, li, ui) {
  rtruncated(N=N, lo = li, hi = ui, pf = pexp, qf = qexp, rate = lambda) 
}


p_exp_truncation <- function(x,lam,a,b,my_log=F){
  return((dexp(x,rate=lam,log = my_log)*1*(x>a)*(x<b))/(pexp(b,rate=lam, log.p = my_log)-pexp(a,rate=lam, log.p = my_log)))
}

log_lik_dtruncexp = function(y, eta, lb = 0, ub = Inf){
  p = dtruncexp(y, eta, lb , ub )
  return(prod((p)))
}


sampling_SST_matrix_exp = function(k,lam){
  Y = matrix(NA, k, k)
  diag(Y) = 0
  
  for(i in 2:k) {
    Y[1, i] = exp_trunc(N=1, lambda = lam, li = Y[1, i-1], ui = Inf)
  }
  
  for(i in 2:(k-1)) {
    for(j in (i+1):k) {
      Y[i, j] = exp_trunc(N=1, lambda = lam, li = Y[i, j-1], ui = Y[i-1, j])
    }
  }
  return(Y)
}
########
##NORMAL
#########

pdf_norm_trunc = function(y,mu,sigma,a,b){
  pdf= (dnorm(y,mu,sigma)*(y>a)*(y<=b))/(pnorm(b,mu,sigma)-pnrom(a,mu,sigma))
  return(pdf)}


#working with eta
log_lik_dtruncnorm = function(y, my_eta, lb = 0, ub = Inf){
  p = log(dtruncnorm(y, eta = c(0, my_eta), lb , ub ))
  return(sum(p))
}

sampling_SST_matrix_norm = function(k,sigma){
  Y = matrix(NA, k, k)
  diag(Y) = 0
  
  for(i in 2:k) {
    Y[1, i] = rtruncnorm(N=1, 0, sigma, a = Y[1, i-1], b = Inf)
  }
  
  for(i in 2:(k-1)) {
    for(j in (i+1):k) {
      Y[i, j] = rtruncnorm(N=1,  0, sigma, a = Y[i, j-1], b = Y[i-1, j])
    }
  }
  return(Y)
}


sample_norm_trunc=function(N, m=0, s=1, a=0, b=1) {
  u=runif(N)
  x=qnorm(((pnorm(b, mean= m, sd =  s) -pnorm(a, mean=m,sd =   s))*u + pnorm(a, mean=m, sd =  s)),mean = m,sd = s)
  return(x)
}

########
##logNORMAL: pdf,cdf,rsample, truncation
#########

dlnorm_mu_sigma = function(y,mu=1, sigma=1, log=T){
  lmu = log(mu**2/(sqrt(mu**2+sigma**2)))
  lsigma = sqrt(log(1+ sigma**2/mu**2))
  return(dlnorm(y,lmu, lsigma,log = log))
}


calculate_victory_probabilities_modified <- function(z_mat, P) {
  # Convert z_mat to a sparse matrix if it's mostly zeros
  z_mat <- sparseMatrix(i = row(z_mat), j = col(z_mat), x = z_mat)
  
  
  # Calculate intermediate result and cache it
  aux <- P %*% t(z_mat)
  
  # Perform matrix multiplication using cached intermediate result
  p_ij_scanning <- z_mat %*% aux
  
  return(p_ij_scanning)
}
plnorm_mu_sigma = function(y,mu=1, sigma=1,log=T){
  lmu = log(mu**2/(sqrt(mu**2+sigma**2)))
  lsigma = sqrt(log(1+ sigma**2/mu**2))
  return(plnorm(y,lmu, lsigma,log = log))
}
rlnorm_mu_sigma = function(N,mu=1,sigma=1){
  lmu = log(mu**2/(sqrt(mu**2+sigma**2)))
  lsigma = sqrt(log(1+ sigma**2/mu**2))
  rlnorm(N,lmu, lsigma)
}

dlnorm_mu_sigma_trunc = function(y,mu,sigma,a,b,log=T){
  pdf= (dlnorm_mu_sigma(y = y,mu = mu,sigma = sigma,log = log)*(y>a)*(y<=b))/(plnorm_mu_sigma(y = b,mu = mu,sigma = sigma,log = log)-plnorm_mu_sigma(y = a,mu = mu,sigma = sigma,log = log))
  return(pdf)
}



###############
##AUXILIARIES
###############
generate_color_gradient_K <- function(K) {
  # Define the starting and ending colors
  start_color <- "deepskyblue"
  end_color <- "darkorchid3"
  
  # Generate the gradient of colors
  gradient_colors <- colorRampPalette(c(start_color, end_color))(K)
  
  return(gradient_colors)
}


semi_symmetric = function(Y){
  Y = (1 -t(Y*upper.tri(Y))) * (lower.tri(Y, diag = F)*1) + upper.tri(Y, diag = T) *Y
  return(Y)}

make_symmetric = function(a_given_matrix){
  a_given_matrix[lower.tri(a_given_matrix)] <- t(a_given_matrix)[lower.tri(a_given_matrix)]
  if(isSymmetric.matrix(a_given_matrix)){
    print("congrats, you have been symmetrified!")
    return(a_given_matrix)}
  else{print("something went wrong")}
}

inv_transformation = function(X)(log(X/(1-X)))


matrix_to_dataframe = function(M){
  my_k= ncol(M)
  df = data.frame(y=0 ,lb=0, ub=0)
  
  for(i in 2:my_k) {
    df= rbind(df, data.frame(y=M[1, i], lb= M[1, i-1], ub= Inf))
  }
  
  for(i in 2:(my_k-1)) {
    for(j in (i+1):my_k) {
      df= rbind(df, data.frame(y=M[i, j], lb=  M[i, j-1], ub=M[i-1, j]))
    }
  }
  df[,2] = df[,2]
  return(df[-1,])
}

#creating the auxiliary dataframe
df_aux = function(matches, z_players, p_matrix){
  ok1 = (ncol(z_players) == 2)
  ok2 = setequal(unique(z_players[,1]),unique(c(matches[,1]),matches[,2]))
  ok3 = setequal((1:ncol(p_matrix)), unique(z_players[,2]))
  ok4 = (ncol(matches) == 4)
  
  if(ok1 * ok2 * ok3 * ok4){
    colnames(z_players) <- c("id", "z_current")
    colnames(matches) <- c("player1", "player2","n_ij", "y_ij") 
    N= nrow(z_players)
    matches = cbind(matches, z_player1 = rep(0,N),z_player2 =rep(0,N), p_ij =rep(0,N) )
    for(i in 1:N){
      matches[i,]$z_player1 = z_players[which(z_players$id == matches[i,]$player1),]$z_current
      matches[i,]$z_player2 = z_players[which(z_players$id == matches[i,]$player2),]$z_current
      matches$p_ij[i] = p_matrix[matches$z_player1[i],matches$z_player2[i]]
    }
    return(matches)
  }else{
    violated_conditions =which(c(ok1,ok2,ok3,ok4) != 1)
    print(paste("Check conditions", violated_conditions))}
}


df_aux_fast <- function(matches, z_players, p_matrix) {
  # if(ncol(z_players) == 2 &&
  #       setequal(unique(z_players[, 1]), unique(c(matches[, 1], matches[, 2]))) &&
  #       setequal(1:ncol(p_matrix), unique(z_players[, 2])) &&
  #       ncol(matches) == 4){print("ok")}else{
  #   stop("Check conditions")
  # }
  # 
  colnames(z_players) <- c("id", "z_current")
  colnames(matches) <- c("player1", "player2", "n_ij", "y_ij")
  
  matches$z_player1 <- left_join(data.frame(id = matches$player1),z_players, by = "id")$z
  matches$z_player2 <- left_join(data.frame(id = matches$player2),z_players, by = "id")$z
  # matches$z_player1 <- z_players$z_current[match(matches$player1, z_players$id)]
  # matches$z_player2 <- z_players$z_current[match(matches$player2, z_players$id)]
  matches$p_ij <- p_matrix[cbind(matches$z_player1, matches$z_player2)]
  
  return(matches)
}

#Function to have the vector of labels available to assign
labels_add<-function(LABELS){
  CLUSTERS<-length(LABELS)
  #first scenario is the first is missing
  first_one_missing<- LABELS[1]!=1
  if(first_one_missing){
    #scenario where the first is missing
    NEW_LABELS_first_missing<-c(1,LABELS)
    return(list(NEW_LABEL_VECTOR=NEW_LABELS_first_missing, LABEL_ADDED=1))
  }else{
    diff_lab<-diff(LABELS)
    all_ones<- all(diff_lab==1)
    if(all_ones){
      #I need to add the last label
      NEW_LABELS_last_missing<-c(LABELS, CLUSTERS+1)
      return(list(NEW_LABEL_VECTOR=NEW_LABELS_last_missing, LABEL_ADDED=CLUSTERS+1))
    }else{
      label_to_insert<-which(diff_lab!=1)[1]+1
      NEW_LABELS_middle<-append(x=LABELS,values=label_to_insert, after=label_to_insert-1 )
      return(list(NEW_LABEL_VECTOR=NEW_LABELS_middle, LABEL_ADDED=label_to_insert))
    }
  }
}


#Simulating a synthethic tournament

generating_synthetic_data = function(alpha, beta, K, N, max_number_games ){
  # alpha: first parameter of the beta(alpha,beta) distribution
  # beta:second parameter of the beta(alpha,beta) distribution
  # K: number of clusters in the generated network
  # N: number of players 
  # max_number_games: maximum number of times 2 players are allowed to play against each other
  
  
  #P matrix generating the SST P matrix ---
  # the code is explained in the "prior" file
  p_fake = sampling_SST_matrix_beta(K,alpha, beta)
  p_fake[is.na(p_fake)] <- 0
  p_fake = (1 -t(p_fake*upper.tri(p_fake))) * (lower.tri(p_fake, diag = F)*1) + upper.tri(p_fake, diag = T) *p_fake
  
  #Z vector: assigning to players a random block membership ----
  players_fake = data.frame(id = c(1:N), z_true = sample(1:K, N, replace = T) )
  #----
  
  ## Simulating the tournament: the two first columns contain the pair of players 
  #who have playerd against each other ------------
  
  
  data_clean_fake = data.frame(player1 = rep(1,N), player2 = rep(1,N), z_player1=rep(3,N),z_player2=rep(4,N))
  #check not to have two identical rows
  check1 <- sum(data_clean_fake$player1 == data_clean_fake$player2)<1
  #check whether at least every player has played once
  check2 <- (setequal(unique(players_fake[,1]),unique(c(data_clean_fake[,1],data_clean_fake[,2]))))
  while(check1 + check2 !=2){
    data_clean_fake = data.frame(player1 = sample(c(1:N),N,replace = F), player2 = sample(c(1:N),N,replace = F), z_player1=rep(0,N),z_player2=rep(0,))
    check1 <- sum(data_clean_fake$player1 == data_clean_fake$player2)<1
    check2 <- (setequal(unique(players_fake[,1]),unique(c(data_clean_fake[,1],data_clean_fake[,2]))))
  }
  
  
  #including the block membership for player1 and player 2 into the main dataframe
  for(i in 1:nrow(data_clean_fake)){
    data_clean_fake$z_player1[i] = players_fake[which(players_fake$id == data_clean_fake$player1[i]),]$z_true
    data_clean_fake$z_player2[i] =players_fake[which(players_fake$id == data_clean_fake$player2[i]),]$z_true
  }
  
  #sampling n_ij for each pair of players
  data_clean_fake = cbind(data_clean_fake, n_ij = sample((1:max_number_games),N, replace = T))
  #sampling y_ij for each pair of players according to y ~ Binom(n_ij, p_{zi,zj})
  
  data_clean_fake = cbind(data_clean_fake, y_ij = rep(0,N))
  for(i in 1:N){
    data_clean_fake$y_ij[i] = rbinom(1, size = data_clean_fake$n_ij[i], p= p_fake[data_clean_fake$z_player1[i],data_clean_fake$z_player2[i]])
  }
  
  #removing the block membership from the dataframe (this info appears into player fake df)
  data_clean_fake = data_clean_fake[,-c(3:4)]
  
  #------------
  
  # The function returns:--------
  # players_fake: a Nx2 dataframe: players ID in the first column, players block membership in the second one
  # matches_result: sa Nx4 dataframe: player1 ID in the first column,player2 ID in the second column players 
  # n_ij in the third column, y_ij in the fourth colomun
  # p_true: KxK matrix containing the interblock connection probabilities
  #-------------
  return(list(z_true = players_fake, matches_results = data_clean_fake, p_true = p_fake))}


#compute the prior on p|alpha
# l_like_p_ij = function(my_matrix, truncations){
#   K = nrow(my_matrix)
#   split_K=  diag_split_matrix(my_matrix)
#   pdf_container = matrix(0,1,K)
#   for(i in (K):1){
#     pdf_container[1,i] = log_lik_dtruncbeta(split_K[[i]],1,1,a = truncations[-i+K+1],b=truncations[-i+2+K])}
#   joint_p = rowSums(pdf_container)
#   return(joint_p)
# }

# p_prime = simulating_POMM_powerlaw(6,1,.6)

l_like_p_ij = function(P_matrix, truncations, diag0.5 =T){
  lowest_level_set_index = ifelse(diag0.5,2,1)
  K = nrow(P_matrix)
  #splitting the matrix into its level sets
  split_K=  diag_split_matrix(P_matrix)
  #container for probabilities
  dlog=matrix(0,K,1)
  for(i in K:lowest_level_set_index){
    #number of elements in level set i
    n_elements = length(unlist(split_K[[i]]))
    #probability
    dlog[i]= log(truncations[K+2-i]-truncations[K+1-i])*n_elements 
    indicator = (sum(split_K[[i]]>truncations[-i+K+1]&split_K[[i]]<truncations[-i+K+2]) == n_elements)*1
  }
  return(sum(-(dlog)))
  # if(indicator == 1){
  # return(sum(-(dlog)))}else{
  #   print("values outside truncations")
  # }
}


l_like_p_ij1 = function(P_matrix, truncations, diag0.5 =T){
  K <- nrow(P_matrix)
  level_sets = diag_split_matrix(P_matrix)
  #considering or excluding the main diagonal
  lowest_level_set_index = ifelse(diag0.5,2,1)
  lbindex = ifelse(diag0.5,1,0)
  log_lik = 0 
  for(i in (1+lbindex):(K)){
    level_set_i <- level_sets[i]
    lb = truncations[i-lbindex]
    ub = truncations[i+1-lbindex]
    log_lik = log_lik + sum(dunif(unlist(level_set_i),min = lb,max = ub,log = T))
  }
  return(log_lik)
}





# computing the likelihood
get_A = function(nij,yij,pij){
  return(sum(lchoose(nij,yij)+(yij*log(pij))+((nij- yij)*log(1-pij))))}
#computing the prior on p|beta matrix
get_B = function(p,b){
  sum((b-1)*log(1- p[upper.tri(p,diag = F)]))
}

#computing the prior on z|K
get_C = function(N,K,n){
  lgamma(K) - lgamma(N+K) - log(factorial(K)) + sum(lgamma(n +1))
}






p_proposal = function(p, sigma_p,K){
  p_new <- matrix(sample_norm_trunc(K*K, p, sigma_p,0,1),K,K)
  return(p_new)
}


simulating_tournament<- function(N,alpha,beta, min_clust,max_clust,n_ij_max){
  #simulating K from a truncated Poisson(1)
  K =  sample_truncated_poisson(N=1,lambda=1,a= min_clust,b = max_clust)
  
  #simulating z|K from dirichlet multinomial with gamma=1
  repeat {
    z<- dir_multinom_r(K,N,1)
    l = length(unique(z))
    #sample until number of sampled labels == number of clusters K
    if (l==K){
      break
    }}
  
  # simulating the matrix p according to the SST property
  p = sampling_SST_matrix_beta(K,alpha, beta)
  p[is.na(p)] <- 0
  p = (1 -t(p*upper.tri(p))) * (lower.tri(p, diag = F)*1) + upper.tri(p, diag = T) *p
  
  #creating the first dataframe
  z_players <- data.frame(id = 1:N, z = z)
  
  
  
  player1=0
  player2=0
  same_set = 0
  same_K = 0 
  different_rows = 0
  #checking some conditions
  while((same_set+same_K+different_rows)!=3){
    player1= sample(x=c(1:N), size=N, replace = F)
    player2= sample(x=c(1:N), size=N, replace= T)
    same_set = setequal(z_players$id, unique(c(player1,player2)))
    same_K = ncol(p) == length(unique(z_players$z))
    different_rows = table(player1 != player2)["TRUE"] == N
  }
  
  z_player1 = left_join(data.frame(id = player1),z_players, by = "id")$z
  z_player2 = left_join(data.frame(id = player2),z_players, by = "id")$z
  
  #here I want the number of games to be inversely related to clusters' strength
  aux = matrix(0, K,K)
  aux = abs(row(aux) - col(aux))  + 1
  diag(aux) = (diag(row(aux))**(-1))/sum(diag(row(aux)**(-1)))
  aux = aux**(-1)/sum(aux**(-1))
  n_ij = sample_truncated_poisson(N,lambda= exp(aux[cbind(z_player1,z_player2)])*n_ij_max,a=1,b=n_ij_max)
  p_ij = p[cbind(z_player1,z_player2)]
  y_ij = rbinom(N, n_ij, p_ij)
  matches = data.frame(player1, player2, n_ij, y_ij)
  
  return(list(z_true = z_players, matches_results = matches, p_true = p, K_true = K))}


# simulating_tournament<- function(N,alpha,beta, min_clust,n_ij_max){
#   require(gtools)
#   #simulating K from a truncated Poisson(1)
#   K =  qpois(((ppois(Inf, 1) -ppois(min_clust, 1))*(runif(1))+ppois(min_clust, 1)),1)
#   
#   #simulating theta
#   theta = rdirichlet(1,c(rep(1/K,K)))
#   
#   # simulating from the multinomial              
#   
#   z = rmultinom(N,1,theta)
#   
#   # simulating the matrix p according to the SST property
#   p = sampling_SST_matrix_beta(K,alpha, beta)
#   p[is.na(p)] <- 0
#   p = (1 -t(p*upper.tri(p))) * (lower.tri(p, diag = F)*1) + upper.tri(p, diag = T) *p
#   
#   z_players <- data.frame(id = 1:N, z = apply(z, 2, which.max))
#   
#   
#   player1=0
#   player2=0
#   same_set = 0
#   same_K = 0 
#   different_rows = 0
#   while((same_set+same_K+different_rows)!=3){
#     player1= sample(x=c(1:N), size=N, replace = F)
#     player2= sample(x=c(1:N), size=N, replace= T)
#     same_set = setequal(z_players$id, unique(c(player1,player2)))
#     same_K = ncol(p) == length(unique(z_players$z))
#     different_rows = table(player1 != player2)["TRUE"] == N
#   }
#   
#   z_player1 = left_join(data.frame(id = player1),z_players, by = "id")$z
#   z_player2 = left_join(data.frame(id = player2),z_players, by = "id")$z
#   
#   
#   
#   n_ij = sample(size = N, x = c(1:n_ij_max),replace = T)
#   p_ij = p[cbind(z_player1,z_player2)]
#   y_ij = rbinom(N, n_ij, p_ij)
#   matches = data.frame(player1, player2, n_ij, y_ij)
#   
#   return(list(z_true = z_players, matches_results = matches, p_true = p))}
vec2mat_0_P <- function(clust_lab,P){
  # in: vector clust_lab of length V s.t. clust_lab[v]=h if node v is in cluster h
  # out: binary VxH matrix M s.t. M[v,h]=1{node v is in cluster h}
  V <- length(clust_lab)
  H <- nrow(P)
  M <- matrix(0,V,H)
  for (v in 1:V){
    M[v,clust_lab[v]] <- 1
  }
  return(M)
}

vec2mat <- function(clust_lab){
  # in: vector clust_lab of length V s.t. clust_lab[v]=h if node v is in cluster h
  # out: binary VxH matrix M s.t. M[v,h]=1{node v is in cluster h}
  V <- length(clust_lab)
  H <- max(clust_lab)
  M <- matrix(0,V,H)
  for (v in 1:V){
    M[v,clust_lab[v]] <- 1
  }
  return(M)
}

# vec2mat <- function(clust_lab) {
#   V <- length(clust_lab)
#   H <- max(clust_lab)
#   M <- matrix(0, V, H)
#   idx <- seq_len(V)
#   M[cbind(idx, clust_lab)] <- 1
#   return(M)
# }



pr_cc <- function(z_post, P){
  # in: posterior sample of assignments (VxN_iter matrix)
  # out: VxV matrix c with elements c[vu]=fraction of iterations in which v and u are in the same cluster
  V <- nrow(z_post)    
  N_iter <- ncol(z_post)
  c <- matrix(0,V,V)
  for (t in 1:N_iter){
    Z <- vec2mat_0_P(z_post[,t], P)
    c <- c + Z%*%t(Z)
  }
  return(c/N_iter)
}

similarity_plot <- function(adj, z_0, z_est){
  ID <- c(1:dim(adj)[1])
  df_a <- data.frame(adj)
  #sorting columns accordin to z_0 
  #(descending order to have the clusters along the main diagonal)
  z_0_a <- data.frame(z_0, ID )
  df_a <- data.frame(adj, ID)
  df_a <- left_join(df_a, z_0_a, by="ID")
  df_a <- df_a[order(df_a$z_0,ID, decreasing = T),]
  df_a$z_0 <- NULL
  df_a$ID <- NULL
  #sorting rows according to z_0
  
  df_a <- t(df_a)
  df_a <- data.frame(df_a, ID)
  df_a <- left_join(df_a, z_0_a, by="ID")
  df_a <- df_a[order(df_a$z_0,ID, decreasing = F),]
  df_a$z_0 <- NULL
  df_a$ID <- NULL
  
  #ordering the nodes of z_est for consistency with z_0
  z_aux <- (z_0_a[order(z_0),])
  df_est_aux <- data.frame(z_est, ID)
  # to plot the membership on the columns
  z_est_ord <- vector()
  for(i in 1:nrow(z_aux)){
    z_est_ord <- append(z_est_ord, z_est[z_aux$ID[i]])
  }
  
  
  #to plot the membership on the rows (descending order)
  z_est_desc <- vector()
  for(i in nrow(z_aux):1){
    z_est_desc <- append(z_est_desc, z_est[z_aux$ID[i]])
  }
  
  
  a_aux <- as.matrix(df_a)
  col_fun <- colorRampPalette(c("white", "gray", "black"))
  col_range <- seq(min(a_aux), max(a_aux), length = 100)
  
  side_col_fun <- colorRampPalette(c("gray", "blue","pink","purple"))
  side_col_range <- seq(1, nlevels(factor(z_est_ord)), length = nlevels(factor(z_est_ord)))
  
  
  
  h_aux <- heatmap(a_aux , col=col_fun(100) , symm = TRUE, Rowv = NA, Colv = NA, 
                   ColSideColors = side_col_fun(nlevels(factor(z_est_ord)))[z_est_ord], 
                   RowSideColors = side_col_fun(nlevels(factor(z_est_desc)))[z_est_desc])
  return(h_aux)
}

get_proposal1 = function(z_current, labels_available) {
  ii = sample(seq_len(nrow(z_current)), 1)
  k_proposal = sample(setdiff(labels_available,z_current[ii,]$z),size=1)
  z_proposal = z_current
  z_proposal[ii, "z"] = k_proposal
  return(z_proposal)
}


get_proposal2 = function(z_current, labels_available, num_nodes) {
  ii = sample(seq_len(nrow(z_current)), num_nodes)
  k_proposal = sample(labels_available,size=1)
  z_proposal = z_current
  z_proposal[ii, "z"] = k_proposal
  return(z_proposal)}


simulating_tournament_test<- function(N,alpha,beta, min_clust,max_clust,M, n_ij_max){
  #simulating K from a truncated Poisson(1)
  K =  sample_truncated_poisson(N=1,lambda=1,a= min_clust,b = max_clust)
  
  #simulating z|K from dirichlet multinomial with gamma=1
  repeat {
    z<- dir_multinom_r(K,N,1)
    l = length(unique(z$y))
    #sample until number of sampled labels == number of clusters K
    if (l==K){
      break
    }}
  
  # simulating the matrix p according to the SST property
  p = sampling_SST_matrix_beta(K,alpha, beta)
  p[is.na(p)] <- 0
  p = (1 -t(p*upper.tri(p))) * (lower.tri(p, diag = F)*1) + upper.tri(p, diag = T) *p
  
  
  #creating the first dataframe
  z_players <- data.frame(id = 1:N, z = z$y)
  
  
  #the probability of playing is inversely related to the cluster size
  aux2=  data.frame(z = c(1:K), p = sort(z$p, decreasing = T))
  #probability of two players being extracted
  aux2= left_join( z_players,aux2,by="z", multiple='all')
  aux2$p = aux2$p/sum( aux2$p)
  
  df = data.frame(player_1=NA, player_2=NA, n_ij=NA, y_ij=NA, p_ij=NA)
  i=0
  while(i<M){
    pl_1_i <- sample(z_players$id, size=1, prob = aux2$p)
    pl_2_i <-sample(setdiff(z_players$id,pl_1_i), 
                    size=1, prob = aux2$p[setdiff(z_players$id,pl_1_i)])
    matches_i = data.frame(pl_1_i, pl_2_i)
    matches_i$z_1 <- z_players[which(z_players$id == pl_1_i),]$z
    matches_i$z_2 <- z_players[which(z_players$id == pl_2_i),]$z
    matches_i$n_ij = sample_truncated_poisson(1,p[cbind(matches_i$z_1,matches_i$z_2)]*n_ij_max ,a = 1,b = n_ij_max)
    
    matches_i$y_ij =  rbinom(n = 1, size = matches_i$n_ij,prob =  p[cbind(matches_i$z_1,matches_i$z_2)])  
    
    df = rbind(df, data.frame(player_1 = matches_i$pl_1_i, 
                              player_2= matches_i$pl_2_i,
                              n_ij = matches_i$n_ij,
                              y_ij = matches_i$y_ij,
                              p_ij = p[cbind(matches_i$z_1,matches_i$z_2)]))
    i=i+1
  }
  
  
  df = df[-1,]
  
  
  my_df = df
  rownames(my_df) = 1:nrow(my_df)
  head(my_df)
  data_clean = data.frame(player_1=NA, player_2=NA, n_ij=NA, y_ij=NA)
  n = nrow(my_df)
  while(n > 0){
    entries = c(which(my_df$player_1 == my_df$player_1[1] & my_df$player_2 == my_df$player_2[1]), which(my_df$player_1 == my_df$player_2[1] & my_df$player_2 == my_df$player_1[1]))
    df_support = my_df[entries,]
    n_ij = sum(df_support$n_ij)
    y_ij = 0
    for(i in 1:nrow(df_support)){
      if(df_support[i,]$player_1 == my_df$player_1[1]){
        y_ij = y_ij + df_support[i,]$y_ij
      }else if(df_support[i,]$player_1 == my_df$player_2[1]){
        y_ij = y_ij + df_support[i,]$n_ij - df_support[i,]$y_ij
      }
    }
    
    # if(n_victories == sum(which(df_support[,1] == edgelist[i,1] & df_support[,2] == edgelist[i,2]))){
    data_clean = rbind(data_clean, data.frame(player_1 = my_df$player_1[1], player_2 = my_df$player_2[1],n_ij = n_ij, y_ij = y_ij))
    # print("ok")
    # }else(print("error"))
    my_df = my_df[-entries,]
    n = nrow(my_df)
  }
  data_clean= data_clean[-1,]
  
  
  
  return(list(z_true = z_players, matches_results = data_clean, p_true = p, K_true = K))}
#############
# Diagonals
###########

semi_symmetric = function(Y){
  Y = (1 -t(Y*upper.tri(Y))) * (lower.tri(Y, diag = F)*1) + upper.tri(Y, diag = T) *Y
  return(Y)}



plot_matrix <- function(my_matrix, cilower, ciupper, model_name, colorscale, width = 400, height = 300) {
  n <- nrow(my_matrix)
  plot_ly(x = 1:n, y = 1:n, z = my_matrix, type = "heatmap", colorscale = colorscale, width = width, height = height) %>%
    #add_surface(x = 1:n, y = 1:n, z = cilower, opacity = 0.4, showscale = FALSE) %>%
    #add_surface(x = 1:n, y = 1:n, z = ciupper, opacity = 0.4, showscale = FALSE) %>%
    colorbar(title = "Mean") %>%
    layout(title = model_name,
           xaxis = list(title = "Column index"),
           yaxis = list(title = "Row index"))
}

beta_mean_var <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu^2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

diag_matrix = function(K){
  my_diagonal_matrix  = matrix(NA,K,K)
  my_diagonal_matrix = row(my_diagonal_matrix) - col(my_diagonal_matrix)
  return(my_diagonal_matrix)
}

generalized_levels <- function(concave_matrix, K, N,diag0.5) {
  excluding_level0 = ifelse(diag0.5,1 ,0)
  level_list <- vector("list", K-excluding_level0)
  
  for (i in 1:N) {
    split_matrices <- diag_split_matrix(concave_matrix[,,i])
    
    for (j in 1:(K-excluding_level0)) {
      level_list[[j]] <- append(level_list[[j]], split_matrices[[j+excluding_level0]])
    }
  }
  
  return(level_list)
}

diag_split_NA = function(K){
  aux_matrix= matrix(NA,K,K)
  diag_indicator <- (row(aux_matrix) - col(aux_matrix))
  diag_indicator<- upper.tri(diag_indicator,diag = T)*diag_indicator + matrix(-1,nrow(diag_indicator),ncol(diag_indicator))*upper.tri(diag_indicator,diag = T)
  L_k = split(aux_matrix, diag_indicator)
  return(L_k)
}

diag_split_matrix = function(matrix_0){
  K = nrow(matrix_0)
  diag_indicator <- (row(matrix_0)-col(matrix_0))
  diag_indicator<- upper.tri(diag_indicator,diag = T)*diag_indicator + matrix(-1,nrow(diag_indicator),ncol(diag_indicator))*upper.tri(diag_indicator,diag = T)
  L_k = split(matrix_0, diag_indicator)
  Level_sets = list()
  for(i in (K):1){
    Level_sets <- append(Level_sets, L_k[i])
  }
  return(Level_sets)
}

#function to sample the truncations
improper_prior5 = function(K,beta_max,alpha, diag0.5 = T){
  K = ifelse(diag0.5, K-1, K)
  x_min =  0
  x_max = (beta_max-0.5)^(1/alpha)
  delta = (x_max - x_min)/(K)
  delta_sum = rep(delta,K)
  points = cumsum(c(x_min,delta_sum))
  truncations = points^alpha + 0.5
  return(truncations)}



simulating_POMM_powerlaw = function(K,  alpha = 1, beta_max){
  K=3
  #we map the proportions on a desired scale
  beta_0 <-  improper_prior5(K,beta_max,alpha)
  
  #L_k is a list containing all the diagonals of the upper triangular matrix
  L_k = diag_split_NA(K)
  
  # Compute beta values
  for(i in (K):1){
    n_items = length(L_k[[i]])
    #generating the beta values using the proportions as truncations
    L_k[[i]] = sample_beta_trunc(n_items,1,1,beta_0[-i+K], beta_0[-i+K+1])
  }
  
  #putting back everything together in the matrix
  result = diag_matrix(K)
  for(i in 0:(K-1)){
    result[which(result==-i)]<-L_k[[K-i]]
  }
  matrix = semi_symmetric(result)
  returning = list("truncations"= beta_0, "matrix"=matrix)
  return(returning)}

simulating_POMM_powerlaw1 = function(K,  alpha = 1,truncations, beta_max){
  
  #we map the proportions on a desired scale
  beta_0 <- truncations
  
  #L_k is a list containing all the diagonals of the upper triangular matrix
  L_k = diag_split_NA(K)
  # Compute beta values
  for(i in (K):1){
    n_items = length(L_k[[i]])
    #generating the beta values using the proportions as truncations
    L_k[[i]] = sample_beta_trunc(n_items,1,1,beta_0[-i+K+1], beta_0[-i+K+2])
  }
  
  #putting back everything together in the matrix
  result = diag_matrix(K)
  for(i in 0:(K-1)){
    result[which(result==-i)]<-L_k[[K-i]]
  }
  matrix = semi_symmetric(result)
  returning = list("truncations"= beta_0, "matrix"=matrix)
  return(returning)}



simulating_POMM_powerlaw2 = function(K,  alpha = 1,truncations, beta_max, diag0.5=T){
  
  #we map the proportions on a desired scale
  beta_0 <- truncations
  
  j_start = ifelse(diag0.5, yes = 1, no = 0)
  K_stop = ifelse(diag0.5, yes = K-1, no = K)
  N_levelset_i = K:1
  
  P_matrix = matrix(0, K,K)
  for( i in 1:K_stop){
    for(j in (i+j_start):K){
      
      level_set = abs(j-i) + abs(1-j_start)
      P_matrix[i,j] = runif(1,beta_0[level_set],beta_0[level_set+1])
    }
  }
  P_matrix[lower.tri(P_matrix)] = 1 - t(P_matrix)[lower.tri(P_matrix)]
  if(diag0.5){
    diag(P_matrix) = rep(0.5,K)
  }
  
  return(P_matrix)}



simulating_tournament_powerlaw<- function(N,alpha, beta_max, min_clust,max_clust,M, n_ij_max, gamma_vec){
  #simulating K from a truncated Poisson(1)
  K =  min_clust
  
  #simulating z|K from dirichlet multinomial with gamma=1
  while(TRUE){
    KtimesNmatrix = rdirichlet_multinomial(N,K,gamma_vec)
    if(sum((rowSums(KtimesNmatrix)>0)) == K){
      break
    }}
  print(rowSums(KtimesNmatrix)>0)
  z = matrix(0,N,1)
  for(i in 1:ncol(KtimesNmatrix)){
    z[i]= which(KtimesNmatrix[,i] >0)}
  
  
  
  # simulating the matrix p according to the SST property
  p = simulating_POMM_powerlaw(K,alpha,beta_max)$matrix
  #creating the first dataframe
  z_players <- data.frame(id = 1:N, z = z)
  
  cluster_proportions = sort(as.vector(table(z)/N), decreasing = T)
  labels_available = seq(1,K,1)
  #the probability of playing is inversely related to the cluster size
  aux2=  data.frame(z = labels_available, p = cluster_proportions)
  #probability of two players being extracted
  aux2= left_join( z_players,aux2,by="z", multiple='all')
  aux2$p = aux2$p/sum( aux2$p)
  
  df = data.frame(player_1=NA, player_2=NA, n_ij=NA, y_ij=NA, p_ij=NA)
  i=0
  while(i<M){
    pl_1_i <- sample(z_players$id, size=1, prob = aux2$p)
    pl_2_i <-sample(setdiff(z_players$id,pl_1_i), 
                    size=1, prob = aux2$p[setdiff(z_players$id,pl_1_i)])
    matches_i = data.frame(pl_1_i, pl_2_i)
    matches_i$z_1 <- z_players[which(z_players$id == pl_1_i),]$z
    matches_i$z_2 <- z_players[which(z_players$id == pl_2_i),]$z
    matches_i$n_ij = sample_truncated_poisson(1,p[cbind(matches_i$z_1,matches_i$z_2)]*n_ij_max ,a = 1,b = n_ij_max)
    
    matches_i$y_ij =  rbinom(n = 1, size = matches_i$n_ij,prob =  p[cbind(matches_i$z_1,matches_i$z_2)])  
    
    df = rbind(df, data.frame(player_1 = matches_i$pl_1_i, 
                              player_2= matches_i$pl_2_i,
                              n_ij = matches_i$n_ij,
                              y_ij = matches_i$y_ij,
                              p_ij = p[cbind(matches_i$z_1,matches_i$z_2)]))
    i=i+1
  }
  
  
  df = df[-1,]
  
  
  my_df = df
  rownames(my_df) = 1:nrow(my_df)
  head(my_df)
  data_clean = data.frame(player_1=NA, player_2=NA, n_ij=NA, y_ij=NA)
  n = nrow(my_df)
  while(n > 0){
    entries = c(which(my_df$player_1 == my_df$player_1[1] & my_df$player_2 == my_df$player_2[1]), which(my_df$player_1 == my_df$player_2[1] & my_df$player_2 == my_df$player_1[1]))
    df_support = my_df[entries,]
    n_ij = sum(df_support$n_ij)
    y_ij = 0
    for(i in 1:nrow(df_support)){
      if(df_support[i,]$player_1 == my_df$player_1[1]){
        y_ij = y_ij + df_support[i,]$y_ij
      }else if(df_support[i,]$player_1 == my_df$player_2[1]){
        y_ij = y_ij + df_support[i,]$n_ij - df_support[i,]$y_ij
      }
    }
    
    # if(n_victories == sum(which(df_support[,1] == edgelist[i,1] & df_support[,2] == edgelist[i,2]))){
    data_clean = rbind(data_clean, data.frame(player_1 = my_df$player_1[1], player_2 = my_df$player_2[1],n_ij = n_ij, y_ij = y_ij))
    # print("ok")
    # }else(print("error"))
    my_df = my_df[-entries,]
    n = nrow(my_df)
  }
  data_clean= data_clean[-1,]
  
  
  
  return(list(z_true = z_players, matches_results = data_clean, p_true = p, K_true = K))}





heat_map_blue = function(matrix, title){
  K = nrow(matrix)
  # convert matrix to data frame
  df <- as.data.frame(t(matrix))
  x = c(1:K)
  y= c(K:1)
  data <- expand.grid(X=x, Y=y)
  df_long <- reshape2::melt(df)
  
  df_plot = data.frame(X =data$X, Y= data$Y, Z = df_long$value)
  # Create heatmap
  ggplot(df_plot, aes(x=X, y=Y, fill=Z)) +
    geom_tile() +
    scale_fill_gradient(low="white", high="blue") +
    scale_x_discrete(limits=as.character(x), name="Column")+
    scale_y_discrete(limits=(as.character(y)), name="Row")+
    labs(title=title, x="Column", y="Row", fill="Prob")+
    geom_text(aes(label = round(Z, 2)), size =2)
}


simulating_tournament_simple_model<- function(N,alpha, beta_max, min_clust,max_clust,M, n_ij_max, gamma_vec){
  #simulating K from a truncated Poisson(1)
  K =  min_clust
  
  #simulating z|K from dirichlet multinomial with gamma=1
  while(TRUE){
    KtimesNmatrix = rdirichlet_multinomial(N,K,gamma_vec)
    if(sum((rowSums(KtimesNmatrix)>0)) == K){
      break
    }}
  print(rowSums(KtimesNmatrix)>0)
  z = matrix(0,N,1)
  for(i in 1:ncol(KtimesNmatrix)){
    z[i]= which(KtimesNmatrix[,i] >0)}
  
  
  
  # simulating the matrix p according to the SST property
  values= runif((K*(K-1))/2 + K,0.5,beta_max)
  p = matrix(0,K,K)
  p[upper.tri(p, diag = T)] = values
  p[is.na(p)] <- 0
  p = (1 -t(p*upper.tri(p))) * (lower.tri(p, diag = F)*1) + upper.tri(p, diag = T) *p
  
  
  #creating the first dataframe
  z_players <- data.frame(id = 1:N, z = z)
  
  cluster_proportions = sort(as.vector(table(z)/N), decreasing = T)
  labels_available = seq(1,K,1)
  #the probability of playing is inversely related to the cluster size
  aux2=  data.frame(z = labels_available, p = cluster_proportions)
  #probability of two players being extracted
  aux2= left_join( z_players,aux2,by="z", multiple='all')
  aux2$p = aux2$p/sum( aux2$p)
  
  df = data.frame(player_1=NA, player_2=NA, n_ij=NA, y_ij=NA, p_ij=NA)
  i=0
  while(i<M){
    pl_1_i <- sample(z_players$id, size=1, prob = aux2$p)
    pl_2_i <-sample(setdiff(z_players$id,pl_1_i), 
                    size=1, prob = aux2$p[setdiff(z_players$id,pl_1_i)])
    matches_i = data.frame(pl_1_i, pl_2_i)
    matches_i$z_1 <- z_players[which(z_players$id == pl_1_i),]$z
    matches_i$z_2 <- z_players[which(z_players$id == pl_2_i),]$z
    matches_i$n_ij = sample_truncated_poisson(1,p[cbind(matches_i$z_1,matches_i$z_2)]*n_ij_max ,a = 1,b = n_ij_max)
    
    matches_i$y_ij =  rbinom(n = 1, size = matches_i$n_ij,prob =  p[cbind(matches_i$z_1,matches_i$z_2)])  
    
    df = rbind(df, data.frame(player_1 = matches_i$pl_1_i, 
                              player_2= matches_i$pl_2_i,
                              n_ij = matches_i$n_ij,
                              y_ij = matches_i$y_ij,
                              p_ij = p[cbind(matches_i$z_1,matches_i$z_2)]))
    i=i+1
  }
  
  
  df = df[-1,]
  
  
  my_df = df
  rownames(my_df) = 1:nrow(my_df)
  head(my_df)
  data_clean = data.frame(player_1=NA, player_2=NA, n_ij=NA, y_ij=NA)
  n = nrow(my_df)
  while(n > 0){
    entries = c(which(my_df$player_1 == my_df$player_1[1] & my_df$player_2 == my_df$player_2[1]), which(my_df$player_1 == my_df$player_2[1] & my_df$player_2 == my_df$player_1[1]))
    df_support = my_df[entries,]
    n_ij = sum(df_support$n_ij)
    y_ij = 0
    for(i in 1:nrow(df_support)){
      if(df_support[i,]$player_1 == my_df$player_1[1]){
        y_ij = y_ij + df_support[i,]$y_ij
      }else if(df_support[i,]$player_1 == my_df$player_2[1]){
        y_ij = y_ij + df_support[i,]$n_ij - df_support[i,]$y_ij
      }
    }
    
    # if(n_victories == sum(which(df_support[,1] == edgelist[i,1] & df_support[,2] == edgelist[i,2]))){
    data_clean = rbind(data_clean, data.frame(player_1 = my_df$player_1[1], player_2 = my_df$player_2[1],n_ij = n_ij, y_ij = y_ij))
    # print("ok")
    # }else(print("error"))
    my_df = my_df[-entries,]
    n = nrow(my_df)
  }
  data_clean= data_clean[-1,]
  
  
  
  return(list(z_true = z_players, matches_results = data_clean, p_true = p, K_true = K))}

adjacent_label_sampler = function(labels_available, z_ith) {
  first = min(labels_available)
  last = max(labels_available)
  
  if (z_ith == first) {
    new_label = first + 1
  } else if (z_ith == last) {
    new_label = last - 1
  } else {
    new_label = ifelse(runif(1) > 0.5, z_ith + 1, z_ith - 1)
  }
  
  return(new_label)
}









calculate_victory_probabilities <- function(z_mat, P) {
  aux <- P %*% t(z_mat)
  z_mat %*% aux
}

is.semi_symmetric=function(my_matrix){
  check = sum(my_matrix[lower.tri(my_matrix)] == 1 - my_matrix[upper.tri(my_matrix)]) == (nrow(my_matrix)*(nrow(my_matrix)-1))/2
  if(check==TRUE){
    print("It is semi-simmetric")}
  else{
    diagnos_matrix= my_matrix*lower.tri(my_matrix) == (1 - t(my_matrix)*lower.tri(my_matrix))
    # diagnos_matrix = matrix(0, nrow(my_matrix),nrow(my_matrix))
    # 
    # for(i in 1:nrow(my_matrix)){
    #   for(j in i:nrow(my_matrix)){
    #     diagnos_matrix[j,i] = (my_matrix[j,i]==1-my_matrix[i,j])*1
    #   }}
  }
  return(diagnos_matrix)
}

semi_similarity_checker = function(z_current, p_nbyn_current){
  check_blox = matrix(c(1,1,2,2,3,3),ncol=2,nrow=3)
  for(p in 1:nrow(check_blox)){
    combi= matrix(0,2,nrow=1)
    first = which(z_current==check_blox[p,1])
    second = which(z_current==check_blox[p,2])
    for(i in 1:length(first)){
      for(j in 1:length(second)){
        combi = rbind(combi, c(first[i],second[j]))
      }
    }
    print(all( p_nbyn_current[cbind(combi)] == 1- p_nbyn_current[cbind(combi[,2],combi[,1])]))
  }
}


simulating_tournament_new<- function(N,alpha, beta_max, K,M, n_ij_max, gamma_vec, model, diag0.5){
  #simulating K from a truncated Poisson(1)
  
  K_true =  K
  labels_available = 1:K
  #expected number of players per block
  while(TRUE){
    expected_players = round(gamma_vec/sum(gamma_vec)*N,0)
    if(sum(expected_players)==N){
      break
    }else if(sum(expected_players) < N){
      expected_players[K] = expected_players[K] +1
      break
    }else{
      expected_players[K] = expected_players[K] -1
      break}
  }
  
  # #simulating z|K from dirichlet multinomial with gamma=1
  # while(TRUE){
  #   KtimesNmatrix = rdirichlet_multinomial(N,K,gamma_vec)
  #   if(sum((rowSums(KtimesNmatrix)>0)) == K){
  #     break
  #   }else{print(paste('Number of sampled blocks =', sum((rowSums(KtimesNmatrix)>0))))}}
  # 
  
  z = vector()
  # for(i in 1:ncol(KtimesNmatrix)){
  #   z[i]= which(KtimesNmatrix[,i] >0)}
  for(i in 1:length(labels_available)){
    z= append(z, rep(i, expected_players[i]))
  }
  
  
  
  if(model == 'POMM'){
    # simulating the matrix p according to the SST property
    truncations= improper_prior5(K,beta_max,alpha = alpha,diag0.5)
    p = simulating_POMM_powerlaw2(K = K,alpha = alpha,truncations = truncations,beta_max = beta_max,diag0.5)
  }else if(model == 'Simple'){
    p = matrix(rbeta(K*K,1,1),K,K)
    diag(p) <- rep(0.5,K)
    p=semi_symmetric(p)
  }
  
  #creating the first dataframe
  z_players <- data.frame(id = 1:N, z = z)
  z_mat_true=vec2mat(z)
  matrix_z_p_true = p%*%t(z_mat_true)
  p_n_true = z_mat_true%*%matrix_z_p_true
  
  similarity_plot(p_n_true,z,z)
  
  #stratfified sampling
  block_size = colSums(z_mat_true)
  print(paste("block size",c(1:K),"=",block_size))
  
  n_ij = matrix(0, N,N)
  i=1
  while(i<=M){
    #1) sample the blocks
    block1 = sample(x = labels_available,size = 1,prob = 1/block_size/1/sum(block_size))
    block2 = sample(x = labels_available,size = 1,prob = 1/block_size/1/sum(block_size))
    #2) sample a random player within that block
    # Get player IDs within each block
    players_block1 <- z_players[z_players$z == block1, "id"]
    players_block2 <- z_players[z_players$z == block2, "id"]
    
    # Avoid selecting the same player from both blocks
    repeat {
      pl_1_i <- sample(players_block1, size = 1)
      pl_2_i <- sample(players_block2, size = 1)
      
      if (pl_1_i != pl_2_i)  # Check if the players are different
        break
    }
    
    # Update n_ij matrix
    n_ij[cbind(pl_1_i, pl_2_i)] <- n_ij[cbind(pl_1_i, pl_2_i)] + 1
    n_ij[cbind(pl_2_i, pl_1_i)] <- n_ij[cbind(pl_2_i, pl_1_i)] + 1
    
    if(i %%(M/10) == 0){
      print(paste((i/M*100),'-th percent process complete'))
    }
    i =i+1
  }
  
  uppertri_nonzero_entries = which(upper.tri(n_ij) & n_ij > 0,arr.ind = T)
  
  y_ij = matrix(0, N,N)
  # Select upper triangular, non-zero entries in n_ij
  n_ij_upper <- n_ij[uppertri_nonzero_entries]
  p_ij <- p_n_true[uppertri_nonzero_entries]
  y_ij[uppertri_nonzero_entries] <- rbinom(length(uppertri_nonzero_entries)/2, n_ij_upper, p_ij)
  y_ij[lower.tri(y_ij)] <- n_ij[lower.tri(n_ij)] - t(y_ij)[lower.tri(y_ij)]
  
  
  
  
  
  print(paste("total number of matches:", sum(n_ij)/2 - sum(diag(n_ij))/2))
  print(paste("total number of victories:", sum(y_ij)))
  
  return(list(n_ij_true = n_ij, y_ij_true = y_ij, p_ij_true = p_n_true, P_matrix = p, z_true= z))
}



#----------------

simulating_tournament_new_norm<- function(N,alpha, beta_max, K,M, n_ij_max, gamma_vec, model, diag0.5){
  #simulating K from a truncated Poisson(1)
  
  K_true =  K
  labels_available = 1:K
  #expected number of players per block
  while(TRUE){
    expected_players = round(gamma_vec/sum(gamma_vec)*N,0)
    if(sum(expected_players)==N){
      break
    }else if(sum(expected_players) < N){
      expected_players[K] = expected_players[K] +1
      break
    }else{
      expected_players[K] = expected_players[K] -1
      break}
  }
  
  # #simulating z|K from dirichlet multinomial with gamma=1
  # while(TRUE){
  #   KtimesNmatrix = rdirichlet_multinomial(N,K,gamma_vec)
  #   if(sum((rowSums(KtimesNmatrix)>0)) == K){
  #     break
  #   }else{print(paste('Number of sampled blocks =', sum((rowSums(KtimesNmatrix)>0))))}}
  # 
  
  z = vector()
  # for(i in 1:ncol(KtimesNmatrix)){
  #   z[i]= which(KtimesNmatrix[,i] >0)}
  for(i in 1:length(labels_available)){
    z= append(z, rep(i, expected_players[i]))
  }
  
  
  
  if(model == 'POMM'){
    # simulating the matrix p according to the SST property
    truncations= improper_prior5(K,beta_max,alpha = alpha,diag0.5)
    p = simulating_POMM_powerlaw_norm(K = K,alpha = alpha,truncations = truncations,beta_max = beta_max,diag0.5)
  }else if(model == 'Simple'){
    p = matrix(rbeta(K*K,1,1),K,K)
    diag(p) <- rep(0.5,K)
    p=semi_symmetric(p)
  }
  
  #creating the first dataframe
  z_players <- data.frame(id = 1:N, z = z)
  z_mat_true=vec2mat(z)
  matrix_z_p_true = p%*%t(z_mat_true)
  p_n_true = z_mat_true%*%matrix_z_p_true
  
  similarity_plot(p_n_true,z,z)
  
  #stratfified sampling
  block_size = colSums(z_mat_true)
  print(paste("block size",c(1:K),"=",block_size))
  
  n_ij = matrix(0, N,N)
  i=1
  while(i<=M){
    #1) sample the blocks
    block1 = sample(x = labels_available,size = 1,prob = 1/block_size/1/sum(block_size))
    block2 = sample(x = labels_available,size = 1,prob = 1/block_size/1/sum(block_size))
    #2) sample a random player within that block
    # Get player IDs within each block
    players_block1 <- z_players[z_players$z == block1, "id"]
    players_block2 <- z_players[z_players$z == block2, "id"]
    
    # Avoid selecting the same player from both blocks
    repeat {
      pl_1_i <- sample(players_block1, size = 1)
      pl_2_i <- sample(players_block2, size = 1)
      
      if (pl_1_i != pl_2_i)  # Check if the players are different
        break
    }
    
    # Update n_ij matrix
    n_ij[cbind(pl_1_i, pl_2_i)] <- n_ij[cbind(pl_1_i, pl_2_i)] + 1
    n_ij[cbind(pl_2_i, pl_1_i)] <- n_ij[cbind(pl_2_i, pl_1_i)] + 1
    
    if(i %%(M/10) == 0){
      print(paste((i/M*100),'-th percent process complete'))
    }
    i =i+1
  }
  
  uppertri_nonzero_entries = which(upper.tri(n_ij) & n_ij > 0,arr.ind = T)
  
  y_ij = matrix(0, N,N)
  # Select upper triangular, non-zero entries in n_ij
  n_ij_upper <- n_ij[uppertri_nonzero_entries]
  p_ij <- p_n_true[uppertri_nonzero_entries]
  y_ij[uppertri_nonzero_entries] <- rbinom(length(uppertri_nonzero_entries)/2, n_ij_upper, p_ij)
  y_ij[lower.tri(y_ij)] <- n_ij[lower.tri(n_ij)] - t(y_ij)[lower.tri(y_ij)]
  
  
  
  
  
  print(paste("total number of matches:", sum(n_ij)/2 - sum(diag(n_ij))/2))
  print(paste("total number of victories:", sum(y_ij)))
  
  return(list(n_ij_true = n_ij, y_ij_true = y_ij, p_ij_true = p_n_true, P_matrix = p, z_true= z))
}
