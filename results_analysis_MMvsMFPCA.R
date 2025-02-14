#######################################################################
#######################################################################
## Download and analyze results of MMvsMFPCA
#######################################################################
#######################################################################
load("D:/La meva unitat/TFM/ResultsMMvsMFPCA/dataframes_MCAR_30.RData")


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


#######################
#TRAIN
#######################

## we create a single data frame will all the info we need to compute RMSEs
list_n <- 100
rmse_y1_MM <- numeric(list_n)
rmse_y2_MM <- numeric(list_n)
rmse_y3_MM <- numeric(list_n)
rmse_y1_MFPCA <- numeric(list_n)
rmse_y2_MFPCA <- numeric(list_n)
rmse_y3_MFPCA <- numeric(list_n)

for(i in 1:list_n){
  #we take the complete data frame TRAING
  data_train <- list_DF[[i]]
  data_train <- data_train %>%
    dplyr::select(id, time, sex, treatment, y1, y3, y5) %>%
    dplyr::rename(
      y1_true_c = y1,
      y2_true_c = y3,
      y3_true_c = y5
    ) %>%
    dplyr::mutate(
      y1_true_m = list_DF_miss[[i]]$y1,
      y2_true_m = list_DF_miss[[i]]$y3,
      y3_true_m = list_DF_miss[[i]]$y5,
      y1_pred_MM = as.vector(t(list_preds_tmod_y1_miss2[[i]])),
      y2_pred_MM = as.vector(t(list_preds_tmod_y2_miss2[[i]])),
      y3_pred_MM = as.vector(t(list_preds_tmod_y3_miss2[[i]]))
    )
  
  #take the vector of rounded times:
  round_times <- round(data_train$time*2)/2
  
  #take the order in the fixed grid
  fix_grid <- seq(from=0, to=10, by=0.5)
  orders <- match(round_times, fix_grid)
  #to a matrix
  orders <- matrix(orders, nrow = 200, ncol = 10, 
                   byrow = TRUE)
  #now we transform the matrix of 200x21 (predictions in the fixed grid)
  # to a matrix 200x10, for the real times matching the approx grid
  y1_pred_MFPCA <- numeric()
  y2_pred_MFPCA <- numeric()
  y3_pred_MFPCA <- numeric()
  for(j in 1:200){
    y1_pred_MFPCA <- c(y1_pred_MFPCA, list_pred_mfpca1_y1[[i]][j,][orders[j,]])
    y2_pred_MFPCA <- c(y2_pred_MFPCA, list_pred_mfpca1_y2[[i]][j,][orders[j,]])
    y3_pred_MFPCA <- c(y3_pred_MFPCA, list_pred_mfpca1_y3[[i]][j,][orders[j,]])
  }
  
  data_train <- data_train %>%
    dplyr::mutate(
      y1_pred_MFPCA = y1_pred_MFPCA,
      y2_pred_MFPCA = y2_pred_MFPCA,
      y3_pred_MFPCA = y3_pred_MFPCA,
      se_y1_MM = (y1_true_c - y1_pred_MM)^2,
      se_y2_MM = (y2_true_c - y2_pred_MM)^2,
      se_y3_MM = (y3_true_c - y3_pred_MM)^2,
      se_y1_MFPCA = (y1_true_c - y1_pred_MFPCA)^2,
      se_y2_MFPCA = (y2_true_c - y2_pred_MFPCA)^2,
      se_y3_MFPCA = (y3_true_c - y3_pred_MFPCA)^2,
    )
  
  ## Now we have all the info in a data frame and we can compute RMSE:
  rmse_y1_MM[i] <- sqrt(mean(data_train$se_y1_MM[is.na(data_train$y1_true_m)]))
  rmse_y2_MM[i] <- sqrt(mean(data_train$se_y2_MM[is.na(data_train$y2_true_m)]))
  rmse_y3_MM[i] <- sqrt(mean(data_train$se_y3_MM[is.na(data_train$y3_true_m)]))
  rmse_y1_MFPCA[i] <- sqrt(mean(data_train$se_y1_MFPCA[is.na(data_train$y1_true_m)]))
  rmse_y2_MFPCA[i] <- sqrt(mean(data_train$se_y2_MFPCA[is.na(data_train$y2_true_m)]))
  rmse_y3_MFPCA[i] <- sqrt(mean(data_train$se_y3_MFPCA[is.na(data_train$y3_true_m)]))

}

################################
##### plotting results: 
#############################

df_results_train <- data.frame("RMSE Y1 MM" = rmse_y1_MM,
                         "RMSE Y2 MM" = rmse_y2_MM,
                         "RMSE Y3 MM" = rmse_y3_MM,
                         "RMSE Y1 MFPCA" = rmse_y1_MFPCA,
                         "RMSE Y2 MFPCA" = rmse_y2_MFPCA,
                         "RMSE Y3 MFPCA" = rmse_y3_MFPCA)

df_res_long_train <- df_results_train %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value")

ploty <- ggplot(df_res_long_train, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot() +
  scale_fill_manual(values = c("RMSE.Y1.MM" = "lightblue",
                               "RMSE.Y1.MFPCA" = "brown3",
                               "RMSE.Y2.MM" = "lightblue",
                               "RMSE.Y2.MFPCA" = "brown3",
                               "RMSE.Y3.MM" = "lightblue",
                               "RMSE.Y3.MFPCA" = "brown3"))  +
  theme_minimal() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 10),
        plot.title = element_text(size = 13),
        legend.position = "none") +
  xlab("") + ylab("RMSE") +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "grey", linewidth = 1)+
  geom_vline(xintercept = 4.5, linetype = "dashed", color = "grey", linewidth = 1)+
  scale_x_discrete(labels = c("RMSE.Y1.MM" = "RMSE Y1\n MM",
                              "RMSE.Y1.MFPCA" = "RMSE Y1\n MFPCA",
                              "RMSE.Y2.MM" = "RMSE Y2\n MM",
                              "RMSE.Y2.MFPCA" = "RMSE Y2\n MFPCA",
                              "RMSE.Y3.MM" = "RMSE Y3\n MM",
                              "RMSE.Y3.MFPCA" = "RMSE Y3\n MFPCA",
                              "eSL.Test.Data" = "eSL\n Test Data"))+
  ggtitle("RMSE MCAR 15% dropout; Training data")

ploty

#######################
#TEST
#######################

## we create a single data frame will all the info we need to compute RMSEs
list_n <- 100
rmse_y1_MM <- numeric(list_n)
rmse_y2_MM <- numeric(list_n)
rmse_y3_MM <- numeric(list_n)
rmse_y1_MFPCA <- numeric(list_n)
rmse_y2_MFPCA <- numeric(list_n)
rmse_y3_MFPCA <- numeric(list_n)

for(i in 1:list_n){
  #we take the complete data frame TESTING DATA
  data_test <- list_DF_test[[i]]
  data_test <- data_test %>%
    dplyr::select(id, time, sex, treatment, y1, y3, y5) %>%
    dplyr::rename(
      y1_true_c = y1,
      y2_true_c = y3,
      y3_true_c = y5
    ) %>%
    dplyr::mutate(
      y1_true_m = list_DF_test_miss[[i]]$y1,
      y2_true_m = list_DF_test_miss[[i]]$y3,
      y3_true_m = list_DF_test_miss[[i]]$y5,
      y1_pred_MM = as.vector(t(list_preds_tmod_y1_miss_test[[i]])),
      y2_pred_MM = as.vector(t(list_preds_tmod_y2_miss_test[[i]])),
      y3_pred_MM = as.vector(t(list_preds_tmod_y3_miss_test[[i]]))
    )
  
  #take the vector of rounded times:
  round_times <- round(data_test$time*2)/2
  
  #take the order in the fixed grid
  fix_grid <- seq(from=0, to=10, by=0.5)
  orders <- match(round_times, fix_grid)
  #to a matrix
  orders <- matrix(orders, nrow = 200, ncol = 10, 
                   byrow = TRUE)
  #now we transform the matrix of 200x21 (predictions in the fixed grid)
  # to a matrix 200x10, for the real times matching the approx grid
  y1_pred_MFPCA <- numeric()
  y2_pred_MFPCA <- numeric()
  y3_pred_MFPCA <- numeric()
  for(j in 1:200){
    y1_pred_MFPCA <- c(y1_pred_MFPCA, list_pred_mfpca1_y1_test[[i]][j,][orders[j,]])
    y2_pred_MFPCA <- c(y2_pred_MFPCA, list_pred_mfpca1_y2_test[[i]][j,][orders[j,]])
    y3_pred_MFPCA <- c(y3_pred_MFPCA, list_pred_mfpca1_y3_test[[i]][j,][orders[j,]])
  }
  
  data_test <- data_test %>%
    dplyr::mutate(
      y1_pred_MFPCA = y1_pred_MFPCA,
      y2_pred_MFPCA = y2_pred_MFPCA,
      y3_pred_MFPCA = y3_pred_MFPCA,
      se_y1_MM = (y1_true_c - y1_pred_MM)^2,
      se_y2_MM = (y2_true_c - y2_pred_MM)^2,
      se_y3_MM = (y3_true_c - y3_pred_MM)^2,
      se_y1_MFPCA = (y1_true_c - y1_pred_MFPCA)^2,
      se_y2_MFPCA = (y2_true_c - y2_pred_MFPCA)^2,
      se_y3_MFPCA = (y3_true_c - y3_pred_MFPCA)^2,
    )
  
  ## Now we have all the info in a data frame and we can compute RMSE:
  rmse_y1_MM[i] <- sqrt(mean(data_test$se_y1_MM[is.na(data_test$y1_true_m)]))
  rmse_y2_MM[i] <- sqrt(mean(data_test$se_y2_MM[is.na(data_test$y2_true_m)]))
  rmse_y3_MM[i] <- sqrt(mean(data_test$se_y3_MM[is.na(data_test$y3_true_m)]))
  rmse_y1_MFPCA[i] <- sqrt(mean(data_test$se_y1_MFPCA[is.na(data_test$y1_true_m)]))
  rmse_y2_MFPCA[i] <- sqrt(mean(data_test$se_y2_MFPCA[is.na(data_test$y2_true_m)]))
  rmse_y3_MFPCA[i] <- sqrt(mean(data_test$se_y3_MFPCA[is.na(data_test$y3_true_m)]))
  
}

################################
##### plotting results: 
#############################

df_results_test <- data.frame("RMSE Y1 MM" = rmse_y1_MM,
                               "RMSE Y2 MM" = rmse_y2_MM,
                               "RMSE Y3 MM" = rmse_y3_MM,
                               "RMSE Y1 MFPCA" = rmse_y1_MFPCA,
                               "RMSE Y2 MFPCA" = rmse_y2_MFPCA,
                               "RMSE Y3 MFPCA" = rmse_y3_MFPCA)

df_res_long_test <- df_results_test %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value")

ploty2 <- ggplot(df_res_long_test, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot() +
  scale_fill_manual(values = c("RMSE.Y1.MM" = "lightblue",
                               "RMSE.Y1.MFPCA" = "brown3",
                               "RMSE.Y2.MM" = "lightblue",
                               "RMSE.Y2.MFPCA" = "brown3",
                               "RMSE.Y3.MM" = "lightblue",
                               "RMSE.Y3.MFPCA" = "brown3"))  +
  theme_minimal() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 10),
        plot.title = element_text(size = 13),
        legend.position = "none") +
  xlab("") + ylab("RMSE") +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "grey", linewidth = 1)+
  geom_vline(xintercept = 4.5, linetype = "dashed", color = "grey", linewidth = 1)+
  scale_x_discrete(labels = c("RMSE.Y1.MM" = "RMSE Y1\n MM",
                              "RMSE.Y1.MFPCA" = "RMSE Y1\n MFPCA",
                              "RMSE.Y2.MM" = "RMSE Y2\n MM",
                              "RMSE.Y2.MFPCA" = "RMSE Y2\n MFPCA",
                              "RMSE.Y3.MM" = "RMSE Y3\n MM",
                              "RMSE.Y3.MFPCA" = "RMSE Y3\n MFPCA",
                              "eSL.Test.Data" = "eSL\n Test Data"))+
  ggtitle("RMSE MCAR 15% dropout; Test data")

ploty2


##################################################################################
#################################################################################
## checks only with MFPCA
##################################################################################
#################################################################################

load("D:/La meva unitat/TFM/ResultsMMvsMFPCA/results_100_MFPCA_MCAR30_4feb2025.RData")

list_n <- 100
rmse_y1_MFPCA <- numeric(list_n)
rmse_y2_MFPCA <- numeric(list_n)
rmse_y3_MFPCA <- numeric(list_n)

for(i in 1:list_n){
  #we take the complete data frame TESTING DATA
  data_test <- list_DF_test[[i]]
  data_test <- data_test %>%
    dplyr::select(id, time, sex, treatment, y1, y3, y5) %>%
    dplyr::rename(
      y1_true_c = y1,
      y2_true_c = y3,
      y3_true_c = y5
    ) %>%
    dplyr::mutate(
      y1_true_m = list_DF_test_miss[[i]]$y1,
      y2_true_m = list_DF_test_miss[[i]]$y3,
      y3_true_m = list_DF_test_miss[[i]]$y5
    )
  
  #take the vector of rounded times:
  round_times <- round(data_test$time*2)/2
  
  #take the order in the fixed grid
  fix_grid <- seq(from=0, to=10, by=0.5)
  orders <- match(round_times, fix_grid)
  #to a matrix
  orders <- matrix(orders, nrow = 200, ncol = 10, 
                   byrow = TRUE)
  #now we transform the matrix of 200x21 (predictions in the fixed grid)
  # to a matrix 200x10, for the real times matching the approx grid
  y1_pred_MFPCA <- numeric()
  y2_pred_MFPCA <- numeric()
  y3_pred_MFPCA <- numeric()
  for(j in 1:200){
    y1_pred_MFPCA <- c(y1_pred_MFPCA, list_pred_mfpca1_y1_test[[i]][j,][orders[j,]])
    y2_pred_MFPCA <- c(y2_pred_MFPCA, list_pred_mfpca1_y2_test[[i]][j,][orders[j,]])
    y3_pred_MFPCA <- c(y3_pred_MFPCA, list_pred_mfpca1_y3_test[[i]][j,][orders[j,]])
  }
  
  data_test <- data_test %>%
    dplyr::mutate(
      y1_pred_MFPCA = y1_pred_MFPCA,
      y2_pred_MFPCA = y2_pred_MFPCA,
      y3_pred_MFPCA = y3_pred_MFPCA,
      se_y1_MFPCA = (y1_true_c - y1_pred_MFPCA)^2,
      se_y2_MFPCA = (y2_true_c - y2_pred_MFPCA)^2,
      se_y3_MFPCA = (y3_true_c - y3_pred_MFPCA)^2,
    )
  
  ## Now we have all the info in a data frame and we can compute RMSE:
  rmse_y1_MFPCA[i] <- sqrt(mean(data_test$se_y1_MFPCA[is.na(data_test$y1_true_m)]))
  rmse_y2_MFPCA[i] <- sqrt(mean(data_test$se_y2_MFPCA[is.na(data_test$y2_true_m)]))
  rmse_y3_MFPCA[i] <- sqrt(mean(data_test$se_y3_MFPCA[is.na(data_test$y3_true_m)]))
  
}

################################
##### plotting results: 
#############################

df_results_test <- data.frame("RMSE Y1 MFPCA" = rmse_y1_MFPCA,
                              "RMSE Y2 MFPCA" = rmse_y2_MFPCA,
                              "RMSE Y3 MFPCA" = rmse_y3_MFPCA)

df_res_long_test <- df_results_test %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value")

ploty2 <- ggplot(df_res_long_test, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot() +
  scale_fill_manual(values = c(
                               "RMSE.Y1.MFPCA" = "brown3",
                               "RMSE.Y2.MFPCA" = "brown3",
                               "RMSE.Y3.MFPCA" = "brown3"))  +
  theme_minimal() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 10),
        plot.title = element_text(size = 13),
        legend.position = "none") +
  xlab("") + ylab("RMSE") +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "grey", linewidth = 1)+
  geom_vline(xintercept = 4.5, linetype = "dashed", color = "grey", linewidth = 1)+
  scale_x_discrete(labels = c(
                              "RMSE.Y1.MFPCA" = "RMSE Y1\n MFPCA",
                              "RMSE.Y2.MFPCA" = "RMSE Y2\n MFPCA",
                              "RMSE.Y3.MFPCA" = "RMSE Y3\n MFPCA",
                              "eSL.Test.Data" = "eSL\n Test Data"))+
  ggtitle("RMSE MFPCA MCAR 30% dropout; Test data")

ploty2
