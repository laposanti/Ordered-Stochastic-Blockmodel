
adjl_creation <- function(M){
  # create an indicator for all diagonals in the matrix
  d <- (row(M) - col(M))
  d<- upper.tri(d,diag = T)*d + matrix(-1,nrow(d),ncol(d))*upper.tri(d,diag = T)
  
  # use split to group on these values
  
  adjl= data.frame(p_ij =NA, adjl_1=NA,adjl_2=NA )
  L_k = split(M, d)
  for(i in 1:(length(L_k)-2)){
    n = length(as.vector(L_k[[i]]))
    for(j in 1:n){
      adjl = rbind(adjl, data.frame(p_ij =as.vector(L_k[[i]])[j], 
                                    adjl_1 = as.vector(L_k[[i+1]])[j],
                                    adjl_2 = as.vector(L_k[[i+1]])[j+1]))
    }
  }
  adjl = adjl[-1,]
  return(adjl)
}


lambda_0=0.04
#Computing the likelihood
my_df_1 <- adjl_creation(means1)
my_df_1 <- cbind(my_df_1, data.frame(M= apply(my_df_1[,c(2:3)],1,max)) )
log_lik_dtruncbeta(my_df_1$p_ij, 
                   a = my_df_1$M, 
                   b = .99, 
                   alpha = estBetaParams(my_df_1$M,lambda_0)$alpha, 
                   beta = rep(estBetaParams(my_df_1$M,lambda_0)$beta, nrow(my_df_1)))
                   