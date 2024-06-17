#experiment


# Define a function to check the matrix
check_SST <- function(mat) {
  n <- nrow(mat)
  if (n != ncol(mat)) {
    stop("The matrix is not square.")
  }
  
  # Extract upper triangular part (excluding diagonal)
  # Extract upper triangular part (excluding diagonal)
  upper_tri <- mat
  upper_tri[lower.tri(upper_tri, diag = TRUE)] <- NA
  
  
  
  # Check if entries in each column are decreasing
  col_increasing <- sapply(2:n, function(j) all(na.omit(upper_tri[, j]) == sort(na.omit(upper_tri[, j]), decreasing = T)))
  
  # Check if entries in each row are icreasing
  row_decreasing <- sapply(1:(n-1), function(i) all(na.omit(upper_tri[i, ]) == sort(na.omit(upper_tri[i, ]), decreasing = F)))
  
  # Combine the results
  boolean=  all(col_increasing) && all(row_decreasing) && all(mat[upper.tri(mat,diag =F)] >= 0.5)
  return(boolean)
}

check_WST <- function(mat) {
  n <- nrow(mat)
  if (n != ncol(mat)) {
    stop("The matrix is not square.")
  }
  
  
  
  
  # Check if entries in each column are increasing
  greater_than_0.5 = all(mat[upper.tri(mat,diag = F)] >= 0.5)
  
  # Combine the results
  return(greater_than_0.5)
}



experiment_df = data.frame(K = 2:8, check_across_samplesSST = NA, check_across_samplesWST = NA)
for(K in 2:8){
  
  check_across_samplesWST = 0
  check_across_samplesSST = 0
  
  for(i in 1:10000){
    check_across_permutation_WST = 0
    check_across_permutation_SST = 0
    P_sample = matrix(0,K,K)
    P_sample[upper.tri(P_sample,diag =T)]<-runif((K*(K-1))/2+K)
    P_sample[lower.tri(P_sample)] = 1 - P_sample[upper.tri(P_sample)]
    for(i in 1:K){
      
      new_order = c(i:K, 1:(i-1))
      if(i==1){
        new_order = 1:K
      }
      
      mat_permuted <- P_sample[new_order, new_order]
      # if(check_WST(mat_permuted)>0){
      #   print("WST matrix")
      #   print(mat_permuted)
      # }
      # 
      # if(check_SST(mat_permuted)>0){
      #   print("SST matrix")
      #   print(mat_permuted)
      # }
      check_across_permutation_WST = check_across_permutation_WST + check_WST(mat_permuted)
      check_across_permutation_SST=  check_across_permutation_SST+ check_SST(mat_permuted)
      
    }
    if(check_across_permutation_WST > 0 ){
      check_across_samplesWST =check_across_samplesWST+1
    }
    if(check_across_permutation_SST>0){
      check_across_samplesSST=check_across_samplesSST+1
    }
  }
  experiment_df[which(experiment_df$K == K),2] = check_across_samplesSST
  experiment_df[which(experiment_df$K == K),3] = check_across_samplesWST

}
experiment_df
