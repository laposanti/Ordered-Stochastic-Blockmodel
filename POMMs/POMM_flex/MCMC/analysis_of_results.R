
library(label.switching)
library(collpcm)
library(loo)

source('/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/MCMC/Inference_functions.R')

#where the data are stored
data_wd<- "/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/MCMC/results/old_studies"

setwd(data_wd)

getwd()

model<- "POMM"
filename <- list.files(pattern = paste0("True_ModelSimpleEst_model_", model))
print(filename)

#uploading the data
results_gen <- data.frame(
  MAP = 0,
  MINVI = 0,
  MISCLASSERROR = 0,
  WAIC_est = 0,
  WAIC_se = 0
)

#where to save the results
setwd("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/MCMC/results/old_studies/plots")
for(i in 1:length(filename)){
  a<- obtaining_resultsI(filename[i],model = model)
  results_gen<- rbind(results_gen,a)
}

results_gen = results_gen[-1,]
results_gen = round(results_gen,2)


setwd("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/MCMC/results/old_studies/plots")
resultsII_gen = data.frame(alphaest = 0, alpha0.05 = 0 , alpha0.95 =0 , overlapest =0, overlap0.05 = 0 , overlap0.95 =0)
resultsII_gen<- resultsII_gen[-1,]
for(i in 1:length(filename)){
  a<- obtaining_resultsII(filename[i])
  resultsII_gen<- rbind(resultsII_gen,a)
}
resultsII_gen<- round(resultsII_gen,2)
resultsII_gen = results_genII[-1,]

setwd("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/MCMC/results/old_studies/plots")
resultsIII_gen = data.frame(p_est = 0, p_est0.05 = 0 , p_est0.95 =0 , mean_MSE = 0)
for(i in 1:length(filename)){
  a<- obtaining_resultsII(filename[i])
  resultsII_gen<- rbind(resultsII_gen,a)
}
resultsII_gen<- resultsII_gen[-1,]


resultsIII_gen = data.frame(counter5perc = 0, counter1perc =0 , mean_MSE = 0)

for(i in 1:length(filename)){
  a<- obtaining_resultsIII(filename[i],"Simple")
  resultsIII_gen<- rbind(resultsIII_gen,a)
}
resultsIII_gen<- resultsIII_gen[-1,]
resultsIII_gen <- round(resultsIII_gen,2)


obj_POMM<- readRDS(paste0(data_wd,"/",filename[1]))



  
  
  
  
  results = data.frame(counter5perc = 0, counter1perc =0 , mean_MSE = 0)
  obj_POMM<- readRDS(paste0(data_wd,"/",filename[1]))
  print(filename)
  
  K = nrow(obj_POMM$p_true)
  M = sum(obj_POMM$Yij_matrix)
  
  entries<- which(upper.tri(obj_POMM$p_true),arr.ind = T)
  
  
  
  
  