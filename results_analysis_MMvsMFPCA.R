#######################################################################
#######################################################################
## Download and analyze results of MMvsMFPCA
#######################################################################
#######################################################################
load("D:/La meva unitat/TFM/ResultsMMvsMFPCA/results_10_MCAR_27jan2025.RData")

##let us do some plots to compare predictions

library(ggplot2)
library(tidyr)
library(dplyr)

## Y1 with BAYESIAN (missing data)
# Convert the matrix to a data frame
df_y1_bayes <- as.data.frame(list_preds_tmod_y1_miss2[[1]])

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
df_y1_mfpca <- as.data.frame(list_pred_mfpca1_y1[[1]])

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
df_y2_bayes <- as.data.frame(list_preds_tmod_y2_miss2[[1]])

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
df_y2_mfpca <- as.data.frame(list_pred_mfpca1_y2[[1]])

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


#########################
##Computing R