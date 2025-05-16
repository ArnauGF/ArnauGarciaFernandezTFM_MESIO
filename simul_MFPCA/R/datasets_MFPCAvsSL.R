###############################################################################
###############################################################################
## We generate 250 data sets for the simul study MFPCA vs SL
## - mJM is the data-generarting model
## - L=4, 4 LMM
## - Fixed obs points with a noise 
###############################################################################
###############################################################################
library(JMbayes2)
library(MASS)
library(lattice)
library(dplyr)
library(tidyr)

n <- 175 # number of subjects
n_test <- 175
K <- 10 # number of measurements per subject

set.seed(1234)
sigmaa <- matrix(c(runif(1,0.25,0.75), runif(8, -0.0005, 0.0005), 
                   runif(1,0,0.15), runif(8, -0.0005, 0.0005),
                   runif(1,0.25,0.75), runif(8, -0.0005, 0.0005),
                   runif(1,0,0.15), runif(8, -0.0005, 0.0005),
                   runif(1,0.25,0.75), runif(8, -0.0005, 0.0005),
                   runif(1,0,0.025), runif(8, -0.0005, 0.0005),
                   runif(1,0.25,0.75), runif(8, -0.005, 0.005),
                   runif(1,-0.05,0.15))
                 ,nrow = 8)
treat <- rep(as.factor(sample(c("A", "B"), size = n, replace = TRUE))
             ,each = K)

list_DF <- list()
list_DF_test <- list()

repl <- 250
for(count in 1:repl){
  n <- 175 # number of subjects
  n_test <- 175
  K <- 10 # number of measurements per subject
  
  DF <- data.frame(id = rep(seq_len(n), each = K),
                   sex = rep(gl(2, n/2, labels = c("male", "female")), each = K),
                   treat = treat,
                   time = c(replicate(n, c(0, 1:(K-1) + runif(K-1,-0.25,0.25)))))
  DF_test <- data.frame(id = rep(seq_len(n_test), each = K),
                        sex = rep(gl(2, n_test/2, labels = c("male", "female")), each = K),
                        treat = treat,
                        time = c(replicate(n, c(0, 1:(K-1) + runif(K-1,-0.25,0.25)))))
  
  #correlated random effects
  b <-  mvrnorm(n = n, mu = rep(0,8), Sigma = sigmaa%*%t(sigmaa))
  b_test <-  mvrnorm(n = n, mu = rep(0,8), Sigma = sigmaa%*%t(sigmaa))
  
  ## Now we define the model
  
  ##Train:
  X1 <- model.matrix(~ sex + time, data = DF)
  Z1 <- model.matrix(~ time, data = DF)
  
  ##Test:
  X1_test <- model.matrix(~ sex + time, data = DF_test)
  Z1_test <- model.matrix(~ time, data = DF_test)
  
  ##Paramaters specification:
  #b1 <- mvrnorm(n = n, mu = c(0,0), Sigma = matrix(c(0.067, 0.0025, 0.0000178, 0.0008988), 
  #                                                 nrow = 2, ncol = 2))
  b1 <- b[,c(1,2)]
  betas1 <- c(0.123, 0.098765, 1.5)
  sigma1 <- 0.125 #error term sd
  
  ## Longitudinal profile simulation:
  eta_y1 <- as.vector(X1 %*% betas1 + rowSums(Z1 * b1[DF$id, ]))
  DF$y1 <- rnorm(n * K, mean = eta_y1, sd = sigma1)
  
  
  #b1_test <- mvrnorm(n = n_test, mu = c(0,0), Sigma = matrix(c(0.067, 0.0025, 0.0000178, 0.0008988), 
  #                                                           nrow = 2, ncol = 2))
  b1_test <- b_test[,c(1,2)]
  eta_y1_test <- as.vector(X1_test %*% betas1 + rowSums(Z1_test * b1_test[DF_test$id, ]))
  DF_test$y1 <- rnorm(n * K, mean = eta_y1_test, sd = sigma1)
  
  X2 <- model.matrix(~ time, data = DF)
  Z2 <- model.matrix(~ time, data = DF)
  X2_test <- model.matrix(~ time, data = DF_test)
  Z2_test <- model.matrix(~ time, data = DF_test)
  
  ##Paramaters specification:
  #b2 <- mvrnorm(n = n, mu = c(0,0), Sigma = matrix(c(0.067, 0.00025, 0.00178, 0.0008988), 
  #                                                 nrow = 2, ncol = 2))
  b2 <- b[, c(3,4)]
  betas2 <- c(15.641, -1.25918)
  sigma2 <- 0.125 #error term sd
  
  ##Longitudinal profile simulation
  eta_y2 <- as.vector(X2 %*% betas2 + rowSums(Z2 * b2[DF$id, ]))
  DF$y2 <- rnorm(n * K, mean = eta_y2, sd = sigma2)
  
  #b2_test <- mvrnorm(n = n_test, mu = c(0,0), Sigma = matrix(c(0.067, 0.00025, 0.00178, 0.0008988), 
  #                                                           nrow = 2, ncol = 2))
  b2_test <- b_test[, c(3,4)]
  eta_y2_test <- as.vector(X2_test %*% betas2 + rowSums(Z2_test * b2_test[DF_test$id, ]))
  DF_test$y2 <- rnorm(n * K, mean = eta_y2_test, sd = sigma2)
  
  #We define a function to compute the square of a given number
  pow2 <- function(x){
    x^2
  }
  
  X3 <- model.matrix(~ treat + pow2(time), data = DF)
  Z3 <- model.matrix(~ pow2(time), data = DF)
  X3_test <- model.matrix(~ treat + pow2(time), data = DF_test)
  Z3_test <- model.matrix(~ pow2(time), data = DF_test)
  
  ##Paramaters specification:
  #b3 <- mvrnorm(n = n, mu = c(0,0), Sigma = matrix(c(0.009123, 0.000025, 0.0005716, 0.00008988), 
  #                                                 nrow = 2, ncol = 2))
  b3 <- b[, c(5,6)]
  betas3 <- c(1.314159, 0.125,0.08181)
  sigma3 <- 0.125 #error term sd
  
  ##Longitudinal profile simulation
  eta_y3 <- as.vector(X3 %*% betas3 + rowSums(Z3 * b3[DF$id, ]))
  DF$y3 <- rnorm(n * K, mean = eta_y3, sd = sigma2)
  
  
  #b3_test <- mvrnorm(n = n_test, mu = c(0,0), Sigma = matrix(c(0.009123, 0.000025, 0.0005716, 0.00008988), 
  #                                                           nrow = 2, ncol = 2))
  b3_test <- b_test[, c(5,6)]
  eta_y3_test <- as.vector(X3_test %*% betas3 + rowSums(Z3_test * b3_test[DF_test$id, ]))
  DF_test$y3 <- rnorm(n * K, mean = eta_y3_test, sd = sigma3)
  
  ## y4, glmm to model binary long outcome
  X4 <- model.matrix(~ sex + treat + time, data = DF)
  Z4 <- model.matrix(~ time, data = DF)
  X4_test <- model.matrix(~ sex + treat + time, data = DF_test)
  Z4_test <- model.matrix(~ time, data = DF_test)
  
  ##Paramaters specification:
  #b4 <- mvrnorm(n = n, mu = c(0,0), Sigma = matrix(c(0.009123, 0.000025, 0.0005716, 0.00008988), 
  #                                                 nrow = 2, ncol = 2))
  b4 <- b[, c(7,8)]
  betas4 <- c(1.2323, -1.7152, 0.42313, 0.5)
  sigma4 <- 0.125 #error term sd
  
  ##Longitudinal profile simulation
  eta_y4 <- as.vector(X4 %*% betas4 + rowSums(Z4 * b4[DF$id, ]))
  DF$y4 <- rnorm(n * K, mean = eta_y4, sd = sigma4)
  
  
  #b4_test <- mvrnorm(n = n_test, mu = c(0,0), Sigma = matrix(c(0.009123, 0.000025, 0.0005716, 0.00008988), 
  #                                                           nrow = 2, ncol = 2))
  b4_test <- b_test[, c(7,8)]
  eta_y4_test <- as.vector(X4_test %*% betas4 + rowSums(Z4_test * b4_test[DF_test$id, ]))
  DF_test$y4 <- rnorm(n_test * K, mean = eta_y4_test, sd = sigma4)
  
  ##################################
  ### Time-to-event simulation
  ##################################
  
  # We assume baseline hazard to be ctt (exponential distr). We use the inverse 
  # transform sampling method to simulate survival times.
  
  #parameters specification:
  shape_wb <- 1.2 # we assume exp distr and ctt baseline hazard
  #alpha <- c(0.1, -0.2, 0.312) # association coefficients
  #alpha <- c(0.0879, -0.43, 0.2, -0.001)
  alpha <- c(0.0879, -0.43, 0.2, -0.01)
  ##obs: when alpha=0 we know we have convergence
  #alpha <- 2*alpha
  gammas <- c("(Intercept)" = 0.25, "sex" = -0.5)
  #gammas <- c("(Intercept)" = 0.15, "sex" = -0.25)
  #model matrix and linear predictor:
  W <- model.matrix(~ sex, data = DF[!duplicated(DF$id), ])
  treat_ind <- DF[!duplicated(DF$id), ]$treat
  eta_t <- as.vector(W %*% gammas)
  
  invS <- function (t, i) {
    # i denotes the subject
    sex_i <- W[i, 2L]
    treat_i <- treat_ind[i]
    # h() is the hazard function and we assume a Weibull baseline hazard
    h <- function (s) {
      X1_at_s <- cbind(1, sex_i, s)
      Z1_at_s <- cbind(1, s)
      X2_at_s <- cbind(1, s)
      Z2_at_s <- cbind(1, s)
      X3_at_s <- cbind(1, treat_i, pow2(s))
      Z3_at_s <- cbind(1,pow2(s))
      X4_at_s <- cbind(1, sex_i, treat_i, s)
      Z4_at_s <- cbind(1)
      # the linear predictor from the mixed model evaluated at time s
      f1 <- as.vector(X1_at_s %*% betas1 +
                        rowSums(Z1_at_s * b1[rep(i, nrow(Z1_at_s)), ]))
      f2 <- as.vector(X2_at_s %*% betas2 +
                        rowSums(Z2_at_s * b2[rep(i, nrow(Z2_at_s)), ]))
      f3 <- as.vector(X3_at_s %*% betas3 +
                        rowSums(Z3_at_s * b3[rep(i, nrow(Z3_at_s)), ]))
      f4 <- as.vector(X4_at_s %*% betas4 +
                        rowSums(Z4_at_s * b4[rep(i, nrow(Z4_at_s)) ]))
      exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t[i] + f1*alpha[1] + 
            f2*alpha[2] + f3*alpha[3] + f4*alpha[4])
    }
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
  
  
  ## Testing data:
  #model matrix and linear predictor:
  W_test <- model.matrix(~ sex, data = DF_test[!duplicated(DF_test$id), ])
  treat_ind_test <- DF_test[!duplicated(DF_test$id), ]$treat
  eta_t_test <- as.vector(W_test %*% gammas)
  
  invS <- function (t, i) {
    # i denotes the subject
    sex_i <- W_test[i, 2L]
    treat_i <- treat_ind_test[i]
    # h() is the hazard function and we assume a Weibull baseline hazard
    h <- function (s) {
      X1_at_s <- cbind(1, sex_i, s)
      Z1_at_s <- cbind(1, s)
      X2_at_s <- cbind(1, s)
      Z2_at_s <- cbind(1, s)
      X3_at_s <- cbind(1, treat_i, pow2(s))
      Z3_at_s <- cbind(1,pow2(s))
      X4_at_s <- cbind(1, sex_i, treat_i, s)
      Z4_at_s <- cbind(1)
      # the linear predictor from the mixed model evaluated at time s
      f1 <- as.vector(X1_at_s %*% betas1 +
                        rowSums(Z1_at_s * b1_test[rep(i, nrow(Z1_at_s)), ]))
      f2 <- as.vector(X2_at_s %*% betas2 +
                        rowSums(Z2_at_s * b2_test[rep(i, nrow(Z2_at_s)), ]))
      f3 <- as.vector(X3_at_s %*% betas3 +
                        rowSums(Z3_at_s * b3_test[rep(i, nrow(Z3_at_s)), ]))
      f4 <- as.vector(X4_at_s %*% betas4 +
                        rowSums(Z4_at_s * b4_test[rep(i, nrow(Z4_at_s))]))
      exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t_test[i] + f1*alpha[1] + 
            f2*alpha[2] + f3*alpha[3] + f4*alpha[4])
    }
    # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
    integrate(h, lower = 0, upper = t)$value + log(u[i])
  }
  
  # we simulate the event times
  u <- runif(n)
  trueTimes_test <- numeric(n)
  for (i in seq_len(n)) {
    Up <- 100
    Root <- try(uniroot(invS, interval = c(1e-05, Up), i = i)$root, TRUE)
    trueTimes_test[i] <- if (!inherits(Root, "try-error")) Root else 150
  }
  
  ##saving true times
  #checkTimes[count] <- sum(trueTimes==150)/n
  
  #save Times in data frame
  #Time <- trueTimes
  DF$Time <- trueTimes[DF$id]
  
  ##saving true times
  #checkTimes_test[count] <- sum(trueTimes_test==150)/n_test
  
  #save Times in data frame
  #Time <- trueTimes_test
  DF_test$Time <- trueTimes_test[DF_test$id]
  
  # saving datasets in lists
  list_DF <- append(list_DF, list(DF))
  list_DF_test <- append(list_DF_test, list(DF_test))
  
  if(count == repl){
    save(list_DF, list_DF_test,
         file="G:/La meva unitat/TFM/MFPCA_simuls/R/Rdatasets_MFPCAvsSL_29apr2025.RData")
  }
  
  print(count)
  ### end loop  
}


###############################################################################

xyplot(y1 ~ time, groups = id,
       data = DF,
       type = "l" ,xlab="Time",ylab="Y1")

xyplot(y1 ~ time, groups = id,
       data = DF_test,
       type = "l" ,xlab="Time",ylab="Y1")

xyplot(y2 ~ time, groups = id,
       data = DF,
       type = "l" ,xlab="Time",ylab="Y2")

xyplot(y2 ~ time, groups = id,
       data = DF_test,
       type = "l" ,xlab="Time",ylab="Y2")

xyplot(y3 ~ time, groups = id,
       data = DF,
       type = "l" ,xlab="Time",ylab="Y3")

xyplot(y3 ~ time, groups = id,
       data = DF_test,
       type = "l" ,xlab="Time",ylab="Y3")

xyplot(y4 ~ time, groups = id,
       data = DF,
       type = "l" ,xlab="Time",ylab="Y4")

xyplot(y4 ~ time, groups = id,
       data = DF_test,
       type = "l" ,xlab="Time",ylab="Y4")


