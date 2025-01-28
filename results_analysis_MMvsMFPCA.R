#######################################################################
#######################################################################
## Download and analyze results of MMvsMFPCA
#######################################################################
#######################################################################
load("D:/La meva unitat/TFM/ResultsMMvsMFPCA/results_50_MCAR_28jan2025.RData")

##let us do some plots to compare predictions

library(ggplot2)
library(tidyr)
library(dplyr)

## Y1 with BAYESIAN (missing data)
# Convert the matrix to a data frame
df_y1_bayes <- as.data.frame(list_preds_tmod_y1_miss2[[23]])

# Add an "id" column to represent each individual
df_y1_bayes$id <- 1:nrow(df_y1_bayes)

# Reshape the data frame to long format (tidy format)
df_long <- df_y1_bayes %>%
  pivot_longer(
    cols = -id,          # All columns except 'id' will be pivoted
    names_to = "time",    # New column for time points
    values_to = "value"   # New column for longitudinal measurements
  )

# Convert "time" to numeric (since it's treated as a factor after pivoting)
df_long$time <- as.numeric(gsub("V", "", df_long$time))

# Plot the longitudinal profiles (spaghetti plot)
ggplot(df_long, aes(x = time, y = value, group = id)) +
  geom_line(alpha = 0.5) +  # Individual profiles
  labs(
    title = "Spaghetti Plot of Longitudinal Profiles",
    x = "Time",
    y = "Longitudinal Measurement"
  ) +
  theme_minimal()

## Y1 with MFPCA (missing data)
# Convert the matrix to a data frame
df_y1_mfpca <- as.data.frame(list_pred_mfpca1_y1[[23]])

# Add an "id" column to represent each individual
df_y1_mfpca$id <- 1:nrow(df_y1_mfpca)

# Reshape the data frame to long format (tidy format)
df_long_mfpca <- df_y1_mfpca %>%
  pivot_longer(
    cols = -id,          # All columns except 'id' will be pivoted
    names_to = "time",    # New column for time points
    values_to = "value"   # New column for longitudinal measurements
  )

# Convert "time" to numeric (since it's treated as a factor after pivoting)
df_long_mfpca$time <- as.numeric(gsub("V", "", df_long_mfpca$time))

# Plot the longitudinal profiles (spaghetti plot)
ggplot(df_long_mfpca, aes(x = time, y = value, group = id)) +
  geom_line(alpha = 0.5) +  # Individual profiles
  labs(
    title = "Spaghetti Plot of Longitudinal Profiles",
    x = "Time",
    y = "Longitudinal Measurement"
  ) +
  theme_minimal()


## Y2 with BAYESIAN (missing data)
# Convert the matrix to a data frame
df_y2_bayes <- as.data.frame(list_preds_tmod_y2_miss2[[13]])

# Add an "id" column to represent each individual
df_y2_bayes$id <- 1:nrow(df_y2_bayes)

# Reshape the data frame to long format (tidy format)
df_long <- df_y2_bayes %>%
  pivot_longer(
    cols = -id,          # All columns except 'id' will be pivoted
    names_to = "time",    # New column for time points
    values_to = "value"   # New column for longitudinal measurements
  )

# Convert "time" to numeric (since it's treated as a factor after pivoting)
df_long$time <- as.numeric(gsub("V", "", df_long$time))

# Plot the longitudinal profiles (spaghetti plot)
ggplot(df_long, aes(x = time, y = value, group = id)) +
  geom_line(alpha = 0.5) +  # Individual profiles
  labs(
    title = "Spaghetti Plot of Longitudinal Profiles",
    x = "Time",
    y = "Longitudinal Measurement"
  ) +
  theme_minimal()

## Y2 with MFPCA (missing data)
# Convert the matrix to a data frame
df_y2_mfpca <- as.data.frame(list_pred_mfpca1_y2[[13]])

# Add an "id" column to represent each individual
df_y2_mfpca$id <- 1:nrow(df_y2_mfpca)

# Reshape the data frame to long format (tidy format)
df_long_mfpca <- df_y2_mfpca %>%
  pivot_longer(
    cols = -id,          # All columns except 'id' will be pivoted
    names_to = "time",    # New column for time points
    values_to = "value"   # New column for longitudinal measurements
  )

# Convert "time" to numeric (since it's treated as a factor after pivoting)
df_long_mfpca$time <- as.numeric(gsub("V", "", df_long_mfpca$time))

# Plot the longitudinal profiles (spaghetti plot)
ggplot(df_long_mfpca, aes(x = time, y = value, group = id)) +
  geom_line(alpha = 0.5) +  # Individual profiles
  labs(
    title = "Spaghetti Plot of Longitudinal Profiles",
    x = "Time",
    y = "Longitudinal Measurement"
  ) +
  theme_minimal()


###########################################################################
#######Computing RMSE
###########################################################################

# we write a matrix with the real values:

## TRAINING

matrix_y1_train_COMPLETE <- matrix(list_DF[[1]]$y1, nrow = 200, ncol = 10, 
                                                     byrow = TRUE)
matrix_y2_train_COMPLETE <- matrix(list_DF[[1]]$y3, nrow = 200, ncol = 10, 
                                   byrow = TRUE)
matrix_y3_train_COMPLETE <- matrix(list_DF[[1]]$y5, nrow = 200, ncol = 10, 
                                   byrow = TRUE)

matrix_y1_train_MISS <- matrix(list_DF_miss[[1]]$y1, nrow = 200, ncol = 10, 
                                   byrow = TRUE)
matrix_y2_train_MISS <- matrix(list_DF_miss[[1]]$y3, nrow = 200, ncol = 10, 
                                   byrow = TRUE)
matrix_y3_train_MISS <- matrix(list_DF_miss[[1]]$y5, nrow = 200, ncol = 10, 
                                   byrow = TRUE)

## TEST

matrix_y1_test_COMPLETE <- matrix(list_DF_test[[1]]$y1, nrow = 200, ncol = 10, 
                                   byrow = TRUE)
matrix_y2_test_COMPLETE <- matrix(list_DF_test[[1]]$y3, nrow = 200, ncol = 10, 
                                   byrow = TRUE)
matrix_y3_test_COMPLETE <- matrix(list_DF_test[[1]]$y5, nrow = 200, ncol = 10, 
                                   byrow = TRUE)

matrix_y1_test_MISS <- matrix(list_DF_test_miss[[1]]$y1, nrow = 200, ncol = 10, 
                               byrow = TRUE)
matrix_y2_test_MISS <- matrix(list_DF_test_miss[[1]]$y3, nrow = 200, ncol = 10, 
                               byrow = TRUE)
matrix_y3_test_MISS <- matrix(list_DF_test_miss[[1]]$y5, nrow = 200, ncol = 10, 
                               byrow = TRUE)

## we create a single data frame will all the info we need to compute RMSEs
list_n <- 50
for(i in 1:list_n){
  #we take the complete data frame TRAING
  data_train <- list_DF[[1]]
  data_train <- data_train %>%
    dplyr::select(id, time, sex, treatment, y1, y3, y5) %>%
    dplyr::rename(
      y1_true_c = y1,
      y2_true_c = y3,
      y3_true_c = y5
    ) %>%
    dplyr::mutate(
      y1_true_m = list_DF_miss[[1]]$y1,
      y2_true_m = list_DF_miss[[1]]$y3,
      y3_true_m = list_DF_miss[[1]]$y5,
      y1_pred_MM = as.vector(t(list_preds_tmod_y1_miss2[[1]])),
      y2_pred_MM = as.vector(t(list_preds_tmod_y2_miss2[[1]])),
      y3_pred_MM = as.vector(t(list_preds_tmod_y3_miss2[[1]]))
    )
  
  #take the vector of rounded times:
  round_times <- round(data_train$time*2)/2
  
  #take the order in the fixed grid
  fix_grid <- seq(from=0, to=10, by=0.5)
  orders <- match(yk, fix_grid)
  #to a matrix
  orders <- matrix(orders, nrow = 200, ncol = 10, 
                   byrow = TRUE)
  #now we transform the matrix of 200x21 (predictions in the fixed grid)
  # to a matrix 200x10, for the real times matching the approx grid
  y1_pred_MFPCA <- numeric()
  y2_pred_MFPCA <- numeric()
  y3_pred_MFPCA <- numeric()
  for(j in 1:200){
    y1_pred_MFPCA <- c(y1_pred_MFPCA, list_pred_mfpca1_y1[[1]][j,][orders[j,]])
    y2_pred_MFPCA <- c(y2_pred_MFPCA, list_pred_mfpca1_y2[[1]][j,][orders[j,]])
    y3_pred_MFPCA <- c(y3_pred_MFPCA, list_pred_mfpca1_y3[[1]][j,][orders[j,]])
  }
  
  data_train <- data_train %>%
    dplyr::mutate(
      y1_pred_MFPCA = y1_pred_MFPCA,
      y2_pred_MFPCA = y2_pred_MFPCA,
      y3_pred_MFPCA = y3_pred_MFPCA
    )


