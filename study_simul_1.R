###############################################################################
###############################################################################
#SMIULATION STUDY
#we follow (https://drizopoulos.github.io/JMbayes2/articles/Non_Gaussian_Mixed_Models.html)

## This is the simulation study to compare dynamic predictions of the multivariate
## JM (which is the true data-generating model), with the SL-based predictions
###############################################################################
###############################################################################
library(JMbayes2)
library(MASS)

###########################################################################
###########################################################################
## We first generate 200 data sets with L=5 longitudinal outcomes. We generate 
## the event times. Then, since for each scenario we have different censoring
## mechanisms, the latter will be simulated in each scenario.
## We use a full vector of random effects to account for the correlation
## within and between longitudinal outcomes
###########################################################################
###########################################################################
num_dat <- 200

#####D matrix for the random effects:
set.seed(12345)
sigmaa <- matrix(c(runif(1,0,1.5), runif(10, -0.005, 0.005), 
                   runif(1,0,1.5), runif(10, -0.005, 0.005),
                   runif(1,0,1.5), runif(10, -0.005, 0.005),
                   runif(1,0,1.5), runif(10, -0.005, 0.005),
                   runif(1,0,1.5), runif(10, -0.005, 0.005),
                   runif(1,0,1.5), runif(10, -0.005, 0.005),
                   runif(1,0,1.5), runif(10, -0.005, 0.005),
                   runif(1,0,1.5), runif(10, -0.005, 0.005),
                   runif(1,0,1.5), runif(10, -0.005, 0.005),
                   runif(1,0,1.5))
                 ,nrow = 10)


list_DF <- list()
list_DF_test <- list()

checkTimes <- checkTimes2 <- numeric(250)

for(count in 1:num_dat){
  ############################
  #We generate the dataset
  ###########################
  n <- 250 # number of subjects
  n_test <- 250
  K <- 10 # number of measurements per subject
  t_max <- 25 # maximum follow-up time
  
  # we construct a data frame with the design:
  # everyone has a baseline measurement, and then measurements at random 
  # follow-up times up to t_max
  DF <- data.frame(id = rep(seq_len(n), each = K),
                   time = c(replicate(n, c(0, sort(runif(K - 1, 0, t_max))))),
                   sex = rep(gl(2, n/2, labels = c("male", "female")), each = K))
  DF_test <- data.frame(id = rep(seq_len(n_test), each = K),
                        time = c(replicate(n_test, c(0, sort(runif(K - 1, 0, t_max))))),
                        sex = rep(gl(2, n_test/2, labels = c("male", "female")), each = K))
  
  
  # design matrices for the fixed and random effects for longitudinal submodels
  X1 <- model.matrix(~ sex * time, data = DF)
  Z1 <- model.matrix(~ time, data = DF)
  
  X2 <- model.matrix(~ sex + time, data = DF )
  Z2 <- model.matrix(~ time, data = DF)
  
  X3 <- model.matrix(~ time, data = DF )
  Z3 <- model.matrix(~ time, data = DF)
  
  X4 <- model.matrix(~ sex + time, data = DF)
  Z4 <- model.matrix(~ time, data = DF)
  
  X5 <- model.matrix(~ time, data = DF)
  Z5 <- model.matrix(~ time, data = DF)
  
  ####for test data set
  X1_test <- model.matrix(~ sex * time, data = DF_test)
  Z1_test <- model.matrix(~ time, data = DF_test)
  
  X2_test <- model.matrix(~ sex + time, data = DF_test )
  Z2_test <- model.matrix(~ time, data = DF_test)
  
  X3_test <- model.matrix(~ time, data = DF_test )
  Z3_test <- model.matrix(~ time, data = DF_test)
  
  X4_test <- model.matrix(~ sex + time, data = DF_test)
  Z4_test <- model.matrix(~ time, data = DF_test)
  
  X5_test <- model.matrix(~ time, data = DF_test)
  Z5_test <- model.matrix(~ time, data = DF_test)
  
  #Simulate random effects
  #################################################
  b <-  mvrnorm(n = n, mu = rep(0,10), Sigma = sigmaa%*%t(sigmaa))
  b_test <-  mvrnorm(n = n, mu = rep(0,10), Sigma = sigmaa%*%t(sigmaa))
  
  
  #Simulate first longitudinal outcome
  ###############################################
  betas1 <- c(-2.2, -0.25, 1.24, -0.05) # fixed effects coefficients
  sigma1 <- 0.125 # errors sd
  
  # random effects
  b1 <-  b[, c(1,2)]
  # linear predictor
  eta_y1 <- as.vector(X1 %*% betas1 + rowSums(Z1 * b1[DF$id, ]))
  # we simulate normal longitudinal data
  DF$y1 <- rnorm(n * K, mean = eta_y1, sd = sigma1)
  
  #Test data:
  # we simulate random effects
  b1_test <- b_test[, c(1,2)] 
  # linear predictor
  eta_y1_test <- as.vector(X1_test %*% betas1 + rowSums(Z1_test * b1_test[DF_test$id, ]))
  # we simulate normal longitudinal data
  DF_test$y1 <- rnorm(n_test * K, mean = eta_y1_test, sd = sigma1)
  
  #Simulate second longitudinal outcome
  ###############################################
  betas2 <- c(-1.8, -0.06, 0.5) # fixed effects coefficients
  betas2 <- c(50, -6, -3.54789) # fixed effects coefficients
  sigma2 <- 0.25 # errors sd
  
  # we simulate random effects
  b2 <- b[, c(3,4)]
  # linear predictor
  eta_y2 <- as.vector(X2 %*% betas2 + rowSums(Z2 * b2[DF$id, ]))
  # we simulate normal longitudinal data
  DF$y2 <- rnorm(n * K, mean = eta_y2, sd = sigma2)
  
  ####Test data
  # we simulate random effects
  b2_test <- b_test[, c(3,4)]
  # linear predictor
  eta_y2_test <- as.vector(X2_test %*% betas2 + rowSums(Z2_test * b2[DF_test$id, ]))
  # we simulate normal longitudinal data
  DF_test$y2 <- rnorm(n_test * K, mean = eta_y2_test, sd = sigma2)
  
  #Simulate third longitudinal outcome
  ###############################################
  betas3 <- c(-2.5, 0.0333) # fixed effects coefficients
  betas3 <- c(-2.5, 2.33) # fixed effects coefficients
  sigma3 <- 0.25 # errors sd
  
  
  # we simulate random effects
  b3 <- b[, c(5,6)]
  # linear predictor
  eta_y3 <- as.vector(X3 %*% betas3 + rowSums(Z3 * b3[DF$id, ]))
  # we simulate normal longitudinal data
  DF$y3 <- rnorm(n * K, mean = eta_y3, sd = sigma3)
  # we assume that values below 0 are not observed, and set equal to 0
  
  
  ########Test data
  # we simulate random effects
  b3_test <- b_test[, c(5,6)]
  # linear predictor
  eta_y3_test <- as.vector(X3_test %*% betas3 + rowSums(Z3_test * b3_test[DF_test$id, ]))
  # we simulate normal longitudinal data
  DF_test$y3 <- rnorm(n_test * K, mean = eta_y3_test, sd = sigma3)
  
  #Simulate forth longitudinal outcome
  ###############################################
  betas4 <- c(0.01, 0.5, -0.31416) # fixed effects coefficients
  
  # we simulate random effects
  b4 <- b[, c(7,8)]
  # linear predictor
  eta_y4 <- as.vector(X4 %*% betas4 + rowSums(Z4 * b4[DF$id, ]))
  # mean of the binomial distribution
  mu_y4 <- plogis(eta_y4)
  # we simulate binomial longitudinal data
  DF$y4 <- rbinom(n * K, size = 1, prob = mu_y4)
  
  #####Test data
  # we simulate random effects
  b4_test <- b_test[, c(7,8)]
  # linear predictor
  eta_y4_test <- as.vector(X4_test %*% betas4 + rowSums(Z4_test * b4_test[DF_test$id, ]))
  # mean of the binomial distribution
  mu_y4_test <- plogis(eta_y4_test)
  # we simulate binomial longitudinal data
  DF_test$y4 <- rbinom(n_test * K, size = 1, prob = mu_y4_test)
  
  
  #Simulate fifth longitudinal outcome
  ###############################################
  betas5 <- c(1, 0.155) # fixed effects coefficients
  
  # we simulate random effects
  b5 <- b[, c(9,10)]
  # linear predictor
  eta_y5 <- as.vector(X5 %*% betas5 + rowSums(Z5 * b5[DF$id, ]))
  # mean of the binomial distribution
  mu_y5 <- plogis(eta_y5)
  # we simulate binomial longitudinal data
  DF$y5 <- rbinom(n * K, size = 1, prob = mu_y5)
  
  ####Test data
  # we simulate random effects
  b5_test <- b_test[, c(9,10)]
  # linear predictor
  eta_y5_test <- as.vector(X5_test %*% betas5 + rowSums(Z5_test * b5_test[DF_test$id, ]))
  # mean of the binomial distribution
  mu_y5_test <- plogis(eta_y5_test)
  # we simulate binomial longitudinal data
  DF_test$y5 <- rbinom(n_test * K, size = 1, prob = mu_y5_test)
  
  #OBS: other glmm distr (beta, gamma, negative binom, etc) can be used
  
  
  #Simulate event times
  ########################################
  alpha <- c(0.8, 0.61, 0.38, 0.222, 0.74) # association coefficients
  alpha <- (-0.0015)*alpha # association coefficients
  #gammas <- c("(Intercept)" = -15, "sex" = -0.5)
  gammas <- c("(Intercept)" = -20, "sex" = -0.05)
  W <- model.matrix(~ sex, data = DF[!duplicated(DF$id), ])
  # linear predictor for the survival model
  eta_t <- as.vector(W %*% gammas)
  # to simulate event times we use inverse transform sampling
  # (https://en.wikipedia.org/wiki/Inverse_transform_sampling). Namely, we want 
  # to find t, such that S(t) = u, where S(.) is the survival function, and u a 
  # number from the Unif(0, 1) distribution. The function below calculates 
  # log(u) - log(S(t)), and for a given u, we want to find t for which it equals
  # zero. We do that below using the uniroot() function
  invS <- function (t, i) {
    # i denotes the subject
    sex_i <- W[i, 2L]
    # h() is the hazard function and we assume a Weibull baseline hazard
    h <- function (s) {
      X1_at_s <- cbind(1, sex_i, s, sex_i * s)
      Z1_at_s <- cbind(1, s)
      X2_at_s <- cbind(1, sex_i, s)
      Z2_at_s <- cbind(1, s)
      X3_at_s <- cbind(1, s)
      Z3_at_s <- cbind(1, s)
      X4_at_s <- cbind(1, sex_i, s)
      Z4_at_s <- cbind(1, s)
      X5_at_s <- cbind(1, s)
      Z5_at_s <- cbind(1, s)
      # the linear predictor from the mixed model evaluated at time s
      f1 <- as.vector(X1_at_s %*% betas1 +
                        rowSums(Z1_at_s * b1[rep(i, nrow(Z1_at_s)), ]))
      f2 <- as.vector(X2_at_s %*% betas2 +
                        rowSums(Z2_at_s * b2[rep(i, nrow(Z2_at_s)), ]))
      f3 <- as.vector(X3_at_s %*% betas3 +
                        rowSums(Z3_at_s * b3[rep(i, nrow(Z3_at_s)), ]))
      f4 <- as.vector(X4_at_s %*% betas4 +
                        rowSums(Z4_at_s * b4[rep(i, nrow(Z4_at_s)), ]))
      f5 <- as.vector(X5_at_s %*% betas5 +
                        rowSums(Z5_at_s * b5[rep(i, nrow(Z5_at_s)), ]))
      exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t[i] 
          + f1*alpha[1] + f2*alpha[2] + f3*alpha[3] + f4*alpha[4] + f5*alpha[5])
    }
    # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
    integrate(h, lower = 0, upper = t)$value + log(u[i])
  }
  
  # we simulate the event times
  u <- runif(n)
  trueTimes <- numeric(n)
  for (i in seq_len(n)) {
    Up <- 100
    Up <- 50
    Root <- try(uniroot(invS, interval = c(1e-05, Up), i = i)$root, TRUE)
    trueTimes[i] <- if (!inherits(Root, "try-error")) Root else 150
  }

  ##saving true times
  checkTimes[count] <- sum(trueTimes==150)/n
  
  #save Times in data frame
  Time <- trueTimes
  DF$Time <- Time[DF$id]
  
  
  #######Test data##################
  W_test <- model.matrix(~ sex, data = DF_test[!duplicated(DF_test$id), ])
  # linear predictor for the survival model
  eta_t_test <- as.vector(W_test %*% gammas)
  # to simulate event times we use inverse transform sampling
  # (https://en.wikipedia.org/wiki/Inverse_transform_sampling). Namely, we want 
  # to find t, such that S(t) = u, where S(.) is the survival function, and u a 
  # number from the Unif(0, 1) distribution. The function below calculates 
  # log(u) - log(S(t)), and for a given u, we want to find t for which it equals
  # zero. We do that below using the uniroot() function
  invS <- function (t, i) {
    # i denotes the subject
    sex_i <- W_test[i, 2L]
    # h() is the hazard function and we assume a Weibull baseline hazard
    h <- function (s) {
      X1_at_s <- cbind(1, sex_i, s, sex_i * s)
      Z1_at_s <- cbind(1, s)
      X2_at_s <- cbind(1, sex_i, s)
      Z2_at_s <- cbind(1, s)
      X3_at_s <- cbind(1, s)
      Z3_at_s <- cbind(1, s)
      X4_at_s <- cbind(1, sex_i, s)
      Z4_at_s <- cbind(1, s)
      X5_at_s <- cbind(1, s)
      Z5_at_s <- cbind(1, s)
      # the linear predictor from the mixed model evaluated at time s
      f1 <- as.vector(X1_at_s %*% betas1 +
                        rowSums(Z1_at_s * b1_test[rep(i, nrow(Z1_at_s)), ]))
      f2 <- as.vector(X2_at_s %*% betas2 +
                        rowSums(Z2_at_s * b2_test[rep(i, nrow(Z2_at_s)), ]))
      f3 <- as.vector(X3_at_s %*% betas3 +
                        rowSums(Z3_at_s * b3_test[rep(i, nrow(Z3_at_s)), ]))
      f4 <- as.vector(X4_at_s %*% betas4 +
                        rowSums(Z4_at_s * b4_test[rep(i, nrow(Z4_at_s)), ]))
      f5 <- as.vector(X5_at_s %*% betas5 +
                        rowSums(Z5_at_s * b5_test[rep(i, nrow(Z5_at_s)), ]))
      exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t_test[i] 
          + f1*alpha[1] + f2*alpha[2] + f3*alpha[3] + f4*alpha[4] + f5*alpha[5])
    }
    # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
    integrate(h, lower = 0, upper = t)$value + log(u[i])
  }
  
  # we simulate the event times
  u <- runif(n)
  trueTimes2 <- numeric(n)
  for (i in seq_len(n)) {
    Up <- 100
    Up <- 50
    Root <- try(uniroot(invS, interval = c(1e-05, Up), i = i)$root, TRUE)
    trueTimes2[i] <- if (!inherits(Root, "try-error")) Root else 150
  }
  
  ##saving truetimes
  checkTimes2[count] <- sum(trueTimes2==150)/n
  
  #save Times in data frame
  Time <- trueTimes2
  DF_test$Time <- Time[DF_test$id]
  
  #we save data frames in a list
  try(list_DF <- append(list_DF, list(DF)))
  try(list_DF_test <- append(list_DF_test, list(DF_test)))

  
  if(count==num_dat){
    setwd("D:/La meva unitat/TFM/Results_simuls_1")
    save(checkTimes, checkTimes2, list_DF, list_DF_test,
         file="data_simul1.RData")
  }
  print(count) 
}


sum(trueTimes==150)/250
sum(trueTimes2==150)/250
mean(trueTimes[trueTimes!=150])
mean(trueTimes2[trueTimes2!=150])


library(lattice)
xyplot(y1 ~ time, groups = id,
       data = DF[1:825,],
       type = "l" ,xlab="Time",ylab="Y1")

xyplot(y2 ~ time, groups = id,
       data = DF[1:825,],
       type = "l" ,xlab="Time",ylab="Y2")

xyplot(y3 ~ time, groups = id,
       data = DF[1:825,],
       type = "l" ,xlab="Time",ylab="Y3")

xyplot(y4 ~ time, groups = id,
       data = DF[1:825,],
       type = "l" ,xlab="Time",ylab="Y4")

xyplot(y5 ~ time, groups = id,
       data = DF[1:825,],
       type = "l" ,xlab="Time",ylab="Y5")

alpha <- c(0.8, 0.61, 0.38, 0.222, 0.74) # association coefficients
alpha <- alpha/2 # association coefficients
gammas <- c("(Intercept)" = -15, "sex" = -0.5)
gammas <- c("(Intercept)" = -20, "sex" = -0.05)
W <- model.matrix(~ sex, data = DF[!duplicated(DF$id), ])
# linear predictor for the survival model
eta_t <- as.vector(W %*% gammas)

i
sex_i <- W[i, 2L]
h <- function (s) {
  X1_at_s <- cbind(1, sex_i, s, sex_i * s)
  Z1_at_s <- cbind(1, s)
  X2_at_s <- cbind(1, sex_i, s)
  Z2_at_s <- cbind(1, s)
  X3_at_s <- cbind(1, s)
  Z3_at_s <- cbind(1, s)
  X4_at_s <- cbind(1, sex_i, s)
  Z4_at_s <- cbind(1, s)
  X5_at_s <- cbind(1, s)
  Z5_at_s <- cbind(1, s)
  # the linear predictor from the mixed model evaluated at time s
  f1 <- as.vector(X1_at_s %*% betas1 +
                    rowSums(Z1_at_s * b1[rep(i, nrow(Z1_at_s)), ]))
  f2 <- as.vector(X2_at_s %*% betas2 +
                    rowSums(Z2_at_s * b2[rep(i, nrow(Z2_at_s)), ]))
  f3 <- as.vector(X3_at_s %*% betas3 +
                    rowSums(Z3_at_s * b3[rep(i, nrow(Z3_at_s)), ]))
  f4 <- as.vector(X4_at_s %*% betas4 +
                    rowSums(Z4_at_s * b4[rep(i, nrow(Z4_at_s)), ]))
  f5 <- as.vector(X5_at_s %*% betas5 +
                    rowSums(Z5_at_s * b5[rep(i, nrow(Z5_at_s)), ]))
  exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t[i] 
      + f1*alpha[1] + f2*alpha[2] + f3*alpha[3] + f4*alpha[4] + f5*alpha[5])
}

h(5)
integrate(h, lower = 0, upper = 8)$value + log(runif(1))

invS <- function (t, i) {
  # i denotes the subject
  sex_i <- W[i, 2L]
  # h() is the hazard function and we assume a Weibull baseline hazard
  
  # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
  integrate(h, lower = 0, upper = t)$value + log(u[i])
}

# we simulate the event times
u <- runif(n)
trueTimes <- numeric(n)
for (i in seq_len(n)) {
  Up <- 100
  Root <- try(uniroot(invS, interval = c(1e-05, Up), i = i)$root, TRUE)
  trueTimes[i] <- if (!inherits(Root, "try-error")) Root else 150
}



