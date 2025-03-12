#######3
# preliminars simuls
###############
library(JMbayes2)
library(MASS)
library(lattice)
library(dplyr)
library(tidyr)

n <- 250 # number of subjects
n_test <- 250
K <- 10 # number of measurements per subject
t_max <- 20 # m

set.seed(1234)
sigmaa <- matrix(c(runif(1,0.25,0.75), runif(7, -0.0005, 0.0005), 
                   runif(1,0,0.15), runif(7, -0.0005, 0.0005),
                   runif(1,0.25,0.75), runif(7, -0.0005, 0.0005),
                   runif(1,0,0.15), runif(7, -0.0005, 0.0005),
                   runif(1,0.25,0.75), runif(7, -0.0005, 0.0005),
                   runif(1,0,0.025), runif(7, -0.0005, 0.0005),
                   runif(1,0,0.05))
                 ,nrow = 7)
treat <- rep(as.factor(sample(c("A", "B"), size = n, replace = TRUE))
             ,each = K)


## now we run the different scripts
source("informativeCens_simul_SLvsmJM.R")

source("randCens_simulSLvsmJM.R")


source("G:/TFM/SLinJointModels/adminCens_simul_SLvsmJM.R")


