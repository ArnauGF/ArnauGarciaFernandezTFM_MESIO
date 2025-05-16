
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


## now we run the different scripts

# We run the R script with Administrative censoring, 30% of censoring:
#source("G:/TFM/SLinJointModels/simul_SL/adminCens_simul_SLvsmJM.R")

# Random censoring, 30% of censoring
#source("G:/TFM/SLinJointModels/simul_SL/randCens_simulSLvsmJM.R")

# Informative censoring, 30% of censoring
#source("G:/TFM/SLinJointModels/simul_SL/informativeCens_simul_SLvsmJM.R")


# We run the R script with Administrative censoring, 60% of censoring:
source("G:/TFM/SLinJointModels/simul_SL/adminCens_simul_SLvsmJM_60cens.R")

# Random censoring, 60% of censoring
source("G:/TFM/SLinJointModels/simul_SL/randCens_simulSLvsmJM_60cens.R")

# Informative censoring, 60% of censoring
source("G:/TFM/SLinJointModels/simul_SL/informativeCens_simul_SLvsmJM_60cens.R")