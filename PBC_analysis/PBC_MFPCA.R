###########################################################################
###########################################################################
## PBC data set
## - We perform the different case study using PBC data set for the MFPCA
##.  section
###########################################################################
###########################################################################
library(JMbayes2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(MFPCA)

# load the data set pbc2 from JMbayes2 package:
data("pbc2")
data("pbc2.id")

## Survival model
CoxFit <- coxph(Surv(years, status2) ~ drug, data = pbc2.id)

pbc2_noNA <- pbc2[!is.na(pbc2$alkaline),]
pbc2_noNA <- pbc2_noNA[!is.na(pbc2_noNA$platelets),]
## Longitudinal models
LM1 <- lme(log(serBilir) ~ drug + year, data = pbc2_noNA, random = ~ year | id)
LM2 <- lme(albumin ~ sex + year, data = pbc2_noNA, random = ~ year | id)
LM3 <- lme(alkaline ~ year, data = pbc2_noNA, random = ~ year | id)
LM4 <- lme(platelets ~ year, data = pbc2_noNA, random = ~ year | id)

##########################################################################
#Fitting the multivariate JM
multiJM <- jm(CoxFit, list(LM1, LM2, LM3, LM4), time_var = "year",
              n_iter = 12000L, n_burnin = 2000L, n_thin = 5L)

multiJM$statistics$Rhat
summary(multiJM)
# we have a good convergence

t0 <- 5
dt <- 3
## computing IBS (IPCW) and EPCE in (t0, t0+ Dt]=(5, 5+3]
brier_score_multi <- tvBrier(multiJM, newdata = pbc2_noNA, Tstart = t0, Dt = dt, 
                             integrated = TRUE, type_weights = "IPCW")
#EPCE_score_multi <- tvEPCE(multiJM, newdata = pbc2_noNA, Tstart = t0, Dt = dt)

## IBS (IPCW) and EPCE in (t0, t0+ Dt]=(7, 7+3]
brier_score_multi_2 <- tvBrier(multiJM, newdata = pbc2_noNA, Tstart = 7, Dt = 3, 
                               integrated = TRUE, type_weights = "IPCW")
#EPCE_score_multi_2 <- tvEPCE(multiJM, newdata = pbc2_noNA, Tstart = 7, Dt = 3)

##########################################################################
#SuperLearning with the library of models built with the univariate JM

CVdats <- create_folds(pbc2_noNA, V = 5, id_var = "id")

fit_models <- function (data) {
  library("JMbayes2")
  data$status2 <- as.numeric(data$status != "alive")
  data_id <- data[!duplicated(data$id), ]
  CoxFit <- coxph(Surv(years, status2) ~ drug, data = data_id)
  LM1 <- lme(log(serBilir) ~ drug + year, data = data, random = ~ year | id)
  LM2 <- lme(albumin ~ sex + year, data = data, random = ~ year | id)
  LM3 <- lme(alkaline ~ year, data = data, random = ~ year | id)
  LM4 <- lme(platelets ~ year, data = data, random = ~ year | id)
  JM1 <- jm(CoxFit, LM1, time_var = "year")
  JM2 <- jm(CoxFit, LM2, time_var = "year")
  JM3 <- jm(CoxFit, LM3, time_var = "year")
  JM4 <- jm(CoxFit, LM4, time_var = "year")
  out <- list(M1 = JM1, M2 = JM2, M3 = JM3, M4 = JM4)
  class(out) <- "jmList"
  out
}

cl <- parallel::makeCluster(4L)
Models_folds <- parallel::parLapply(cl, CVdats$training, fit_models)
parallel::stopCluster(cl)

# we compute the metrics for the mJM in the CV data in order to do a fair 
# comparison with SL
Brier_multi_CV <- tvBrier(multiJM, newdata = CVdats$testing, 
                          integrated = TRUE, Tstart = t0, Dt = dt,
                          type_weights = "IPCW")

Brier_multi_CV_2 <- tvBrier(multiJM, newdata = CVdats$testing, 
                            integrated = TRUE, Tstart = 7, Dt = 3,
                            type_weights = "IPCW")

EPCE_multi_CV <- tvEPCE(multiJM, newdata = CVdats$testing, 
                        Tstart = t0, Dt = dt)

EPCE_multi_CV_2 <- tvEPCE(multiJM, newdata = CVdats$testing, 
                          Tstart = 7, Dt = 3)

#computing Brier weights
Brier_weights <- tvBrier(Models_folds, newdata = CVdats$testing, 
                         integrated = TRUE, Tstart = t0, Dt = dt,
                         type_weights = "IPCW")
#computing Brier weights in another time window
Brier_weights_2 <- tvBrier(Models_folds, newdata = CVdats$testing, 
                           integrated = TRUE, Tstart = 7, Dt = 3,
                           type_weights = "IPCW")
#computing EPCE weights
#EPCE_weights <- tvEPCE(Models_folds, newdata = CVdats$testing, 
#                       Tstart = t0, Dt = dt)
#computing EPCE weights in another time window
#EPCE_weights_2 <- tvEPCE(Models_folds, newdata = CVdats$testing, 
#                         Tstart = 7, Dt = 3)

#Now we use the weights with the whole data set:
Models <- fit_models(pbc2_noNA)

#we use Brier weights to compute IBS with full data:
bw <- Brier_weights$weights
brier_score_SL_full <- tvBrier(Models, newdata = pbc2_noNA, model_weights = bw, 
                               Tstart = t0, Dt = dt, integrated = TRUE,
                               type_weights = "IPCW")
#we use Brier weights to compute IBS with full data (TIME WINDOW 2):
bw2 <- Brier_weights_2$weights
brier_score_SL_full_2 <- tvBrier(Models, newdata = pbc2_noNA, model_weights = bw2, 
                                 Tstart = 7, Dt = 3, integrated = TRUE,
                                 type_weights = "IPCW")

#we use EPCE weights to compute EPCE with full data:
ew <- EPCE_weights$weights
EPCE_score_SL_full <- tvEPCE(Models, newdata = pbc2_noNA, model_weights = ew,
                             Tstart = t0, Dt = dt)
#we use EPCE weights to compute EPCE with full data (NEW TIME WINDOW):
ew2 <- EPCE_weights_2$weights
EPCE_score_SL_full_2 <- tvEPCE(Models, newdata = pbc2_noNA, model_weights = ew2,
                               Tstart = 7, Dt = 3)


#################################
##Saving results
#################################

## Table for the IBS (IPCW):
data_ibs_results <- data.frame(
  ibs_int1 = c(
               Brier_multi_CV$Brier,
               Brier_weights$Brier,
               Brier_weights$Brier_per_model[1],
               Brier_weights$Brier_per_model[2],
               Brier_weights$Brier_per_model[3],
               Brier_weights$Brier_per_model[4]
  ),
  weights_int1 = c(
                   NA,
                   NA,
                   Brier_weights$weights[1],
                   Brier_weights$weights[2],
                   Brier_weights$weights[3],
                   Brier_weights$weights[4]
  ),
  ibs_int2 = c(
               Brier_multi_CV_2$Brier,
               Brier_weights_2$Brier,
               Brier_weights_2$Brier_per_model[1],
               Brier_weights_2$Brier_per_model[2],
               Brier_weights_2$Brier_per_model[3],
               Brier_weights_2$Brier_per_model[4]
  ),
  weights_int2 = c(
                   NA,
                   NA,
                   Brier_weights_2$weights[1],
                   Brier_weights_2$weights[2],
                   Brier_weights_2$weights[3],
                   Brier_weights_2$weights[4]
  )
)

rownames(data_ibs_results) <- c(
                                "mJM",
                                "SL",
                                "M1",
                                "M2",
                                "M3",
                                "M4")

## Table for the EPCE:
data_epce_results <- data.frame(
  epce_int1 = c(
                EPCE_multi_CV$EPCE,
                EPCE_weights$EPCE,
                EPCE_weights$EPCE_per_model[1],
                EPCE_weights$EPCE_per_model[2],
                EPCE_weights$EPCE_per_model[3],
                EPCE_weights$EPCE_per_model[4]
  ),
  weights_int1 = c(
                   NA,
                   NA,
                   EPCE_weights$weights[1],
                   EPCE_weights$weights[2],
                   EPCE_weights$weights[3],
                   EPCE_weights$weights[4]
  ),
  epce_int2 = c(
                EPCE_multi_CV_2$EPCE,
                EPCE_weights_2$EPCE,
                EPCE_weights_2$EPCE_per_model[1],
                EPCE_weights_2$EPCE_per_model[2],
                EPCE_weights_2$EPCE_per_model[3],
                EPCE_weights_2$EPCE_per_model[4]
  ),
  weights_int2 = c(
                   NA,
                   NA,
                   EPCE_weights_2$weights[1],
                   EPCE_weights_2$weights[2],
                   EPCE_weights_2$weights[3],
                   EPCE_weights_2$weights[4]
  )
)

rownames(data_epce_results) <- c(
                                 "mJM",
                                 "SL",
                                 "M1",
                                 "M2",
                                 "M3",
                                 "M4")

save(data_ibs_results, data_epce_results,
     file = "/Users/arnaugarcia/Desktop/TFM_super/PBC_analysis/results/results_MFPCA.RData")


##############################################################################
# Predictions using MFPCA
#############################################################################
grid_longitudinal_data <- function(DF, n){
  DF_mfpca <- DF
  DF_mfpca$id <- as.integer(DF_mfpca$id)
  DF_mfpca <- DF_mfpca %>%
    mutate(roundtime = round(2 * year) / 2)
  max_len <- as.integer(max(DF_mfpca$roundtime) * 2 + 1)
  obs_time <- seq(0, max_len / 2, by = 0.5)
  
  #we create an empty data frame, full of NA in the longitudinal outcomes
  #we create an empty data frame
  DF_mfpca2 <- data.frame(id = rep(seq_len(n), each = length(obs_time)),
                          year = c(replicate(n, obs_time)),
                          y1_2 = c(replicate(n, rep(NA, length(obs_time)))),
                          y2_2 = c(replicate(n, rep(NA, length(obs_time)))),
                          y3_2 = c(replicate(n, rep(NA, length(obs_time)))),
                          y4_2 = c(replicate(n, rep(NA, length(obs_time)))))
  #we put the values
  # we have to group by id and fill the NAs with values whenever
  #there is info in that time of the grid
  df_filled <- DF_mfpca2 %>%
    left_join(DF_mfpca %>% dplyr::select(id, roundtime, serBilir, albumin, alkaline,
                                         platelets), by = c("id", "year" = "roundtime")) %>%
    mutate(
      y1_2 = ifelse(is.na(y1_2), log(serBilir), y1_2),
      y2_2 = ifelse(is.na(y2_2), albumin, y2_2),
      y3_2 = ifelse(is.na(y3_2), alkaline, y3_2),
      y4_2 = ifelse(is.na(y4_2), platelets, y4_2),
    ) %>%
    dplyr::select(-serBilir, -albumin, -alkaline, -platelets)
  
  # now we save in 3 different data sets with long format in order to use
  # funData()
  df_y1_wide <- df_filled %>%
    dplyr::select(id, year, y1_2) %>%
    distinct(id, year, .keep_all = TRUE) %>%   
    pivot_wider(names_from = year, 
                values_from = y1_2, 
                names_prefix = "y1_year_")
  
  df_y2_wide <- df_filled %>%
    dplyr::select(id, year, y2_2) %>%
    distinct(id, year, .keep_all = TRUE) %>%   
    pivot_wider(names_from = year, 
                values_from = y2_2, 
                names_prefix = "y2_year_")
  
  df_y3_wide <- df_filled %>%
    dplyr::select(id, year, y3_2) %>%
    distinct(id, year, .keep_all = TRUE) %>%   
    pivot_wider(names_from = year, 
                values_from = y3_2, 
                names_prefix = "y3_year_")
  
  df_y4_wide <- df_filled %>%
    dplyr::select(id, year, y4_2) %>%
    distinct(id, year, .keep_all = TRUE) %>%   
    pivot_wider(names_from = year, 
                values_from = y4_2, 
                names_prefix = "y4_year_")
  
  return(list(df_y1_wide, df_y2_wide, df_y3_wide, df_y4_wide, obs_time))
}


#We compute grid longitudinal data:
grid_long <- grid_longitudinal_data(pbc2_noNA, 312)
#univariate functional data
obs_time <- grid_long[[5]]
f1 <- funData(obs_time, as.matrix(grid_long[[1]][,-1]))
f2 <- funData(obs_time, as.matrix(grid_long[[2]][,-1]))
f3 <- funData(obs_time, as.matrix(grid_long[[3]][,-1]))
f4 <- funData(obs_time, as.matrix(grid_long[[4]][,-1]))
#multivariate functional data class
m1 <- multiFunData(list(f1,f2,f3,f4))

m_3_out <- multiFunData(list(f1,f2,f3))


# Grid longitudinal data in TEST:
#grid_long_test <- grid_longitudinal_data(list_DF_test_miss[[count]], n, K)
#univariate functional data TEST case
#obs_time_test <- grid_long_test[[4]]
#f1_test <- funData(obs_time_test, as.matrix(grid_long_test[[1]][,-1]))
#f2_test <- funData(obs_time_test, as.matrix(grid_long_test[[2]][,-1]))
#f3_test <- funData(obs_time_test, as.matrix(grid_long_test[[3]][,-1]))
#multivariate functional data class TEST case
#m1_test <- multiFunData(list(f1_test,f2_test,f3_test))

#now we use MFPCA function from mfpca R package

#We fit with NO weights and pve=0.8
mfpca1 <- MFPCA(m1, M = 3, 
                    uniExpansions = list(list(type="uFPCA", pve = 0.9),
                                         list(type ="uFPCA", pve=0.9),
                                         list(type="uFPCA", pve=0.9),
                                         list(type="uFPCA", pve=0.9)),
                    fit = TRUE)

mfpca_3_out <- MFPCA(m_3_out, M = 3, 
                     uniExpansions = list(list(type="uFPCA", pve = 0.9),
                                          list(type ="uFPCA", pve=0.9),
                                          list(type="uFPCA", pve=0.9)),
                     fit = TRUE)


### predictions:
# predict in training data:
pred_mfpca1 <- predict(mfpca1)
pred_mfpca1_y1 <- pred_mfpca1[[1]]@X
pred_mfpca1_y2 <- pred_mfpca1[[2]]@X
pred_mfpca1_y3 <- pred_mfpca1[[3]]@X
pred_mfpca1_y4 <- pred_mfpca1[[4]]@X


pred_mfpca_3_out <- predict(mfpca_3_out)
pred_mfpca3_y1 <- pred_mfpca_3_out[[1]]@X
pred_mfpca3_y2 <- pred_mfpca_3_out[[2]]@X
pred_mfpca3_y3 <- pred_mfpca_3_out[[3]]@X

## plotting all outcomes together:
plot(pred_mfpca_3_out)

## separately using ggplot2:

set.seed(123456)
individiuals <- sample(1:nrow(pbc2.id), size = 100)

## outcome 1

df_long <- as.data.frame(pred_mfpca3_y1) %>%
  mutate(id = row_number()) %>%     
  pivot_longer(
    cols = -id,                    
    names_to = "time_point",      
    values_to = "measurement"       
  ) %>%
  mutate(time_point = as.numeric(gsub("V", "", time_point)))  

df_long$time_point <- c(replicate(312, obs_time))

df_subs <- df_long %>% 
  filter(id %in% individiuals)

logBilir_mfpca <- ggplot(df_subs, aes(x = time_point, y = measurement, group = id)) +
  geom_line(alpha = 0.5, color = "brown4") +  
  labs(
    x = expression(time),
    y = expression(log(serBilir)(t))
  )  +
  theme_minimal()

## outcome 2

df_long <- as.data.frame(pred_mfpca3_y2) %>%
  mutate(id = row_number()) %>%     
  pivot_longer(
    cols = -id,                    
    names_to = "time_point",      
    values_to = "measurement"       
  ) %>%
  mutate(time_point = as.numeric(gsub("V", "", time_point)))  

df_long$time_point <- c(replicate(312, obs_time))

df_subs <- df_long %>% 
  filter(id %in% individiuals)

albumin_mfpca <- ggplot(df_subs, aes(x = time_point, y = measurement, group = id)) +
  geom_line(alpha = 0.5, color = "brown4") +  
  labs(
    x = expression(time),
    y = expression(albumin(t))
  )  +
  theme_minimal()

#outcome 3 

df_long <- as.data.frame(pred_mfpca3_y3) %>%
  mutate(id = row_number()) %>%     
  pivot_longer(
    cols = -id,                    
    names_to = "time_point",      
    values_to = "measurement"       
  ) %>%
  mutate(time_point = as.numeric(gsub("V", "", time_point)))  

df_long$time_point <- c(replicate(312, obs_time))

df_subs <- df_long %>% 
  filter(id %in% individiuals)

alkaline_mfpca <- ggplot(df_subs, aes(x = time_point, y = measurement, group = id)) +
  geom_line(alpha = 0.5, color = "brown4") +  
  labs(
    x = expression(time),
    y = expression(alkaline(t))
  )  +
  theme_minimal()

### arrenge plot

library(ggpubr)

separator1 <- ggplot(data = NULL) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1) + 
  theme_void() 

arr1 <- ggarrange(logBilir_mfpca, albumin_mfpca, alkaline_mfpca,
                  nrow = 1, ncol = 3)

arr2 <- ggarrange(y1_spagPlot, y2_spagPlot, y3_spagPlot,
                  nrow = 1, ncol = 3)


arreng_plot <- ggarrange(arr2,
                         separator1,
                         arr1,
                         ncol = 1,
                         heights = c(2, 0.1, 2))


pdf("/Users/arnaugarcia/Desktop/TFM_super/PBC_analysis/slides/MFPCA_outcomes.pdf",width=12)
arreng_plot
dev.off()
