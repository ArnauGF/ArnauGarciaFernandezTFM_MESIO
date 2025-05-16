library(JMbayes2)
library(MASS)
library(lattice)
library(dplyr)
library(tidyr)

pow2 <- function(x){
  x^2
}

# we load the lists of data sets
load("G:/La meva unitat/TFM/MFPCA_simuls/R/datasets_MFPCAvsSL_29apr2025.RData")


repl <- 100
list_DF_randCens <- list()
list_DF_test_randCens <- list()

for(count in 1:repl){
  n <- 175 # number of subjects
  n_test <- 175
  K <- 10 # number of measurements per subject
  
  # we set DF ans DF_test to the corresponding element of the list
  DF <- list_DF[[count]]
  DF_test <- list_DF_test[[count]]
  
  ##saving true times
  trueTimes <- DF[!duplicated(DF$id),]$Time
  
  ##saving true times
  trueTimes_test <- DF_test[!duplicated(DF_test$id),]$Time
  
  #simulating censoring; Random censoring
  Ctimes <- rexp(n, 1/18.5)
  Time <- pmin(trueTimes, Ctimes)
  event <- as.numeric(trueTimes <= Ctimes) # event indicator
  # we keep the longitudinal measurements before the event times
  DF$Time <- Time[DF$id]
  DF$event <- event[DF$id]
  DF <- DF[DF$time <= DF$Time, ]
  
  
  ## Testing data:
  Time <- pmin(trueTimes_test, Ctimes)
  event <- as.numeric(trueTimes_test <= Ctimes) # event indicator
  # we keep the longitudinal measurements before the event times
  DF_test$Time <- Time[DF_test$id]
  DF_test$event <- event[DF_test$id]
  DF_test <- DF_test[DF_test$time <= DF_test$Time, ]
  

  
  ## We save the datasets
  list_DF_randCens <- append(list_DF_randCens, list(DF))
  list_DF_test_randCens <- append(list_DF_test_randCens, list(DF_test))
  
  if(count == repl){
    save(list_DF_randCens, list_DF_test_randCens,
         file="G:/La meva unitat/TFM/MFPCA_simuls/R/datasets_randCens_MFPCA_29apr2025.RData")
  }
  
  print(count)
}