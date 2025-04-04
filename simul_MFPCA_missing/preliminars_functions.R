#Simulation study multivariate MM vs MFPCA

###########################################################################
###########################################################################
# General settings
# - L=3 (stan_mvmer() accepts max 3 outcomes)
# - Observations times in a fixed grid every 0.5. Visits every 5 or 3 or 1 year
#   with a noise to have realistic data (each possibility is a different scenario)
#.  by the moment the noise is doing summing U(-0.5,0.5)
# - We focus on dropout as a missingness mechanism, we use 30% and 60% of dropout
# - stan_mvmer() function from rstanarm package to fit multivariate MM under
#   the Bayesian framework
# - MFPCA R package to conduct MFPCA anlaysis
###########################################################################
###########################################################################
library(dplyr)
library(tidyr)
library(MFPCA)
library(rstanarm)
library(MASS)
library(bayesplot)
num_datasets <- 100
n <- 200 # number of subjects
n_test <- 200
K <- 10
#####D matrix for the random effects:
set.seed(12345)
sigmaa <- matrix(c(runif(1,0,1.5), runif(6, -0.005, 0.005), 
                   runif(1,0,1.5), runif(6, -0.005, 0.005),
                   runif(1,0,1.5), runif(6, -0.005, 0.005),
                   runif(1,0,1.5), runif(6, -0.005, 0.005),
                   runif(1,0,1.5), runif(6, -0.005, 0.005),
                   runif(1,0,1.5))
                 ,nrow = 6)
treat <- rep(as.factor(sample(c("A", "B"), size = n, replace = TRUE))
             ,each = K)

####### Function to prepare the longitudinal data in the fixed grid format
grid_longitudinal_data <- function(DF, n, K){
  DF_mfpca <- DF
  DF_mfpca <- DF_mfpca %>%
    mutate(roundtime = round(2 * time) / 2)
  max_len <- as.integer(max(DF_mfpca$roundtime) * 2 + 1)
  obs_time <- seq(0, max_len / 2, by = 0.5)
  
  #we create an empty data frame, full of NA in the longitudinal outcomes
  #we create an empty data frame
  DF_mfpca2 <- data.frame(id = rep(seq_len(n), each = K*2+1),
                          time = c(replicate(n, obs_time)),
                          y1_2 = c(replicate(n, rep(NA, K*2+1))),
                          y3_2 = c(replicate(n, rep(NA, K*2+1))),
                          y5_2 = c(replicate(n, rep(NA, K*2+1))))
  #we put the values
  # we have to group by id and fill the NAs with values whenever
  #there is info in that time of the grid
  df_filled <- DF_mfpca2 %>%
    left_join(DF_mfpca %>% dplyr::select(id, roundtime, y1, y3, y5), by = c("id", "time" = "roundtime")) %>%
    mutate(
      y1_2 = ifelse(is.na(y1_2), y1, y1_2),
      y3_2 = ifelse(is.na(y3_2), y3, y3_2),
      y5_2 = ifelse(is.na(y5_2), y5, y5_2)
    ) %>%
    dplyr::select(-y1, -y3, -y5)
  
  # now we save in 3 different data sets with long format in order to use
  # funData()
  df_y1_wide <- df_filled %>%
    dplyr::select(id, time, y1_2) %>%
    distinct(id, time, .keep_all = TRUE) %>%   
    pivot_wider(names_from = time, 
                values_from = y1_2, 
                names_prefix = "y1_time_")
  
  df_y3_wide <- df_filled %>%
    dplyr::select(id, time, y3_2) %>%
    distinct(id, time, .keep_all = TRUE) %>%   
    pivot_wider(names_from = time, 
                values_from = y3_2, 
                names_prefix = "y3_time_")
  
  df_y5_wide <- df_filled %>%
    dplyr::select(id, time, y5_2) %>%
    distinct(id, time, .keep_all = TRUE) %>%   
    pivot_wider(names_from = time, 
                values_from = y5_2, 
                names_prefix = "y5_time_")
  return(list(df_y1_wide, df_y3_wide, df_y5_wide, obs_time))
}

simulate_dropout <- function(data) {
  # Find the first time when dropout == 1
  dropout_time <- data %>% filter(dropout == 1) %>% slice(1) %>% pull(time)
  
  # If a dropout is found, set all subsequent longitudinal variables to NA
  if (length(dropout_time) > 0) {
    data <- data %>%
      mutate(across(c(y1, y3, y5), ~ ifelse(time >= dropout_time, NA, .)))
  }
  
  return(data)
}

counting_missing_data <- function(data, n, K){
  return(sum(is.na(data$y1))/(n*K))
}
