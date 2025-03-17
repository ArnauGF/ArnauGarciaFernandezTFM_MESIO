#########
# preliminars simuls
###############
library(JMbayes2)
library(MASS)
library(lattice)
library(dplyr)
library(tidyr)

# we load the lists of data sets
load("G:/TFM/SLinJointModels/simul_SL/datasets_SL_14mar2025.RData")

## now we run the different scripts

# We run the R script with Administrative censoring, 30% of censoring:
source("G:/TFM/SLinJointModels/simul_SL/adminCens_simul_SLvsmJM.R")

# Random censoring, 30% of censoring
source("G:/TFM/SLinJointModels/simul_SL/randCens_simulSLvsmJM.R")

# Informative censoring, 30% of censoring
source("G:/TFM/SLinJointModels/simul_SL/informativeCens_simul_SLvsmJM.R")






