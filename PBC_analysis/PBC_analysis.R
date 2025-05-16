###########################################################################
###########################################################################
## PBC data set
## - We do all the graphics from the PBC data set used in the work
## - We perform the different case studies using PBC data set
###########################################################################
###########################################################################
library(JMbayes2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

# load the data set pbc2 from JMbayes2 package:
data("pbc2")
data("pbc2.id")


#######################################################
## Create a table with summary statistics of data
#######################################################
# 1) We have to do the summary for baseline covariates (drug, sex, age)
summary(pbc2.id)
with(pbc2.id, sum(status2==0))

# 2) Summary for time-to-event data (Time, status2)

#follow-up time using reverse KM estimator
sv_fwup <- survfit(Surv(years, status2==0) ~ 1, data = pbc2.id)
summary(sv_fwup)$table["median"]
summary(sv_fwup)$table[c("0.95LCL", "0.95UCL")]

# overall survival time
sv <- survfit(Surv(years, status2) ~ 1, data = pbc2.id)
summary(sv)$table["median"]
summary(sv)$table[c("0.95LCL", "0.95UCL")]

# 3 years survival rate
pbc2.id$year3 <- ifelse(pbc2.id$years <= 3 & pbc2.id$status2 == 1, 1, 0)
table(pbc2.id$year3)
round(table(pbc2.id$year3)/length(pbc2.id$year3)*100, 1)

# 3) Summary for longitudinal data (serBilir, albumin, alkaline, ascites)

summary(pbc2)

pbc2_long <- pbc2 %>%
  select(id, serBilir, albumin, alkaline, ascites, platelets)

obs_counts <- pbc2_long %>%
  group_by(id) %>%
  summarise(
    serBilir_n = sum(!is.na(serBilir)),
    albumin_n = sum(!is.na(albumin)),
    alkaline_n = sum(!is.na(alkaline)),
    ascites_n = sum(!is.na(ascites)),
    platelets_n = sum(!is.na(platelets))
  )


summary(obs_counts)



#######################################################
## Kaplan Meier of the survival
#######################################################
survKM <- survfit(Surv(years, status2) ~ 1, data = pbc2.id)

km_data <- data.frame(
  time = survKM$time,
  surv = survKM$surv,
  lower = survKM$lower,
  upper = survKM$upper
)

kmcurve_p1 <- ggplot(data = km_data, aes(x = time, y = surv)) +
  geom_step(color = "black", linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3,
              fill = "cornflowerblue") +
  labs(
    #title = "",
    x = expression(time),
    y = expression(bolditalic(hat(S)(t)))
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6)
  ) + 
  theme_minimal()

## we plot cumulative hazard via Nelson-Aalen estimator:
nelson_aalen_fit <- survfit(Surv(years, status2) ~ 1, data = pbc2, type = "fh")

nelson_aalen_data <- data.frame(
  time = nelson_aalen_fit$time,         
  cumhaz = nelson_aalen_fit$cumhaz,     
  lower = nelson_aalen_fit$cumhaz - 1.96*nelson_aalen_fit$std.err,
  upper = nelson_aalen_fit$cumhaz + 1.96*nelson_aalen_fit$std.err
)

NA_plot <- ggplot(nelson_aalen_data, aes(x = time, y = cumhaz)) +
  geom_step(color = "black", linewidth = 1) +  # Plot the cumulative hazard as a step function
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "cornflowerblue") +  # Confidence intervals
  labs(#title = "",
       x = expression(time), 
       y = expression(bolditalic(hat(H)(t)))
       ) +
  theme_minimal()

## KM curve for two groups, treated and non-treated
survKM2 <- survfit(Surv(years, status2) ~ drug, data = pbc2.id)

km_data2 <- data.frame(time = survKM2$time, 
                      surv = survKM2$surv, 
                      lower = survKM2$lower, 
                      upper = survKM2$upper, 
                      strata = rep(names(survKM2$strata), survKM2$strata))

# Plot Kaplan-Meier curves using ggplot2
km_curve_drug <- ggplot(km_data2, aes(x = time, y = surv, color = strata)) +
  geom_step() +  # Plot the survival curve as a step function
  labs(x = "Time", y = "Survival Probability", color = "Drug Group") +
  theme_minimal() +
  scale_color_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("blue", "red"))





####################################################
### Spaguetti plots for longitudinal profiles
### Variables of interest:
### - log(serBilir) (serum bilirubin in mg/dl)
### - albumin (albumin in g/dl)
### - alkaline (alkaline phosphatase in U/liter)
### - ascites (a factor with levels No and Yes.)
####################################################

# we select randomly 100 individuals
set.seed(123456)
individiuals <- sample(1:nrow(pbc2.id), size = 100)

#y1 = log(serBilir)
df_subs <- pbc2 %>% 
  filter(id %in% individiuals)

y1_spagPlot <- ggplot(df_subs, aes(x = year, y = log(serBilir), group = id)) +
  geom_line(alpha = 0.5, color = "black") +  
  labs(
    x = expression(time),
    y = expression(log(serBilir)(t))
  ) +
  theme_minimal()


#y2 = albumin
y2_spagPlot <- ggplot(df_subs, aes(x = year, y = albumin, group = id)) +
  geom_line(alpha = 0.5, color = "black") +  
  labs(
    x = expression(time),
    y = expression(albumin(t))
  ) +
  theme_minimal()

#y3 = alkaline
y3_spagPlot <- ggplot(df_subs, aes(x = year, y = alkaline, group = id)) +
  geom_line(alpha = 0.5, color = "black") +
  labs(
    x = expression(time),
    y = expression(alkaline(t))
  ) +
  ylim(c(NA,5000))+
  theme_minimal()

#y4 = ascites (binary outcome)
y4_spagPlot <- ggplot(na.omit(df_subs), aes(x = year, y = ascites, group = id)) +
  geom_line(alpha = 0.5, color = "black") +
  labs(
    x = expression(time),
    y = expression(ascites(t))
  ) +
  theme_minimal()

# another type of visualization for binary longitudinal outcome
y4_linesPlot <- ggplot(na.omit(df_subs), aes(x = year, y = factor(id), color = ascites)) +
  geom_point(size = 2) +    # Plot points based on the ascites status
  geom_line(aes(group = id), color = "black", alpha = 0.8) +  # Add gray lines for the time points per id
  scale_color_manual(values = c("Yes" = "green4", "No" = "red4")) +  # Red for Yes, Green for No
  labs(x = expression(time), y = "Individuals", color = "Ascites Status") +  # Add labels
  theme_minimal() +
  theme(axis.text.y = element_blank(),  # Hide y-axis text for individuals
        panel.grid.major.y = element_blank(),  # Hide major grid lines for y-axis
        panel.grid.minor.y = element_blank())

#y5 = platelets
y5_spagPlot <- ggplot(df_subs, aes(x = year, y = platelets, group = id)) +
  geom_line(alpha = 0.5, color = "black") +
  labs(
    x = expression(time),
    y = expression(platelets(t))
  ) +
  theme_minimal()


#########
# Let us arrenge plots by means of ggarrenge()
#########

arreng_plot <- ggarrange(kmcurve_p1, NA_plot, y1_spagPlot,
                         y2_spagPlot, y3_spagPlot, y4_linesPlot,
                         nrow = 2, ncol = 3)


pdf("/Users/arnaugarcia/Desktop/TFM_super/PBC_analysis/PBC_several_plots.pdf",width=12)
arreng_plot
dev.off()

arreng_plot2 <- ggarrange(kmcurve_p1, NA_plot, y1_spagPlot,
                         y2_spagPlot, y3_spagPlot, y5_spagPlot,
                         nrow = 2, ncol = 3)


pdf("/Users/arnaugarcia/Desktop/TFM_super/PBC_analysis/PBC_several_plots_continous.pdf",width=12)
arreng_plot2
dev.off()
#######################################################################
#######################################################################
## Case study 1: SL applied to PBC data set
## - we use L=4, with variables of interest above, three LMM, one GLMM
## - We fit the multivariate JM and use SL, and we compare
#######################################################################
#######################################################################

## Survival model
CoxFit <- coxph(Surv(years, status2) ~ drug, data = pbc2.id)

pbc2_noNA <- pbc2[!is.na(pbc2$alkaline),]
## Longitudinal models
LM1 <- lme(log(serBilir) ~ drug + year, data = pbc2_noNA, random = ~ year | id)
LM2 <- lme(albumin ~ sex + year, data = pbc2_noNA, random = ~ year | id)
LM3 <- lme(alkaline ~ year, data = pbc2_noNA, random = ~ year | id)
GLM4 <- mixed_model(ascites ~ year, data = pbc2_noNA,
                   random = ~ 1 | id, family = binomial())

##########################################################################
#Fitting the multivariate JM
multiJM <- jm(CoxFit, list(LM1, LM2, LM3, GLM4), time_var = "year",
              n_iter = 12000L, n_burnin = 2000L, n_thin = 5L)

multiJM$statistics$Rhat
summary(multiJM)
# we have a good convergence



######################## 
#### expl on how to compute predictions and plot it
######################## 
t0 <- 5
#We compute a prediction for Patients 25 and 93. 
ND <- pbc2_noNA[pbc2_noNA$id %in% c(25, 93, 111), ]
ND <- ND[ND$year < t0, ]
ND$event <- 0
ND$years <- t0

#In the following example, we calculate predictions from time t0 to time 12
predSurv <- predict(multiJM, newdata = ND, process = "event",
                    times = seq(t0, 12, length.out = 51),
                    return_newdata = TRUE)

png("/Users/arnaugarcia/Desktop/TFM_super/PBC_analysis/slides/indiv_i_dynam_pred.png",
    width = 2000, height = 1600, res = 300)
plot(predSurv, fun_event = function(x) 1-x,
     ylab_event = "Survival probnabilities",
     xlab = "time",
     subject = 111)
dev.off()
######################## 
######################## 

t0 <- 5
dt <- 3
## computing IBS (IPCW) and EPCE in (t0, t0+ Dt]=(5, 5+3]
brier_score_multi <- tvBrier(multiJM, newdata = pbc2_noNA, Tstart = t0, Dt = dt, 
                             integrated = TRUE, type_weights = "IPCW")
EPCE_score_multi <- tvEPCE(multiJM, newdata = pbc2_noNA, Tstart = t0, Dt = dt)

## IBS (IPCW) and EPCE in (t0, t0+ Dt]=(7, 7+3]
brier_score_multi_2 <- tvBrier(multiJM, newdata = pbc2_noNA, Tstart = 7, Dt = 3, 
                               integrated = TRUE, type_weights = "IPCW")
EPCE_score_multi_2 <- tvEPCE(multiJM, newdata = pbc2_noNA, Tstart = 7, Dt = 3)

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
  GLM4 <- mixed_model(ascites ~ year, data = data,
                     random = ~ year | id, family = binomial())
  JM1 <- jm(CoxFit, LM1, time_var = "year")
  JM2 <- jm(CoxFit, LM2, time_var = "year")
  JM3 <- jm(CoxFit, LM3, time_var = "year")
  JM4 <- jm(CoxFit, GLM4, time_var = "year")
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
EPCE_weights <- tvEPCE(Models_folds, newdata = CVdats$testing, 
                       Tstart = t0, Dt = dt)
#computing EPCE weights in another time window
EPCE_weights_2 <- tvEPCE(Models_folds, newdata = CVdats$testing, 
                         Tstart = 7, Dt = 3)

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
  ibs_int1 = c(brier_score_multi$Brier,
               brier_score_SL_full$Brier,
               Brier_multi_CV$Brier,
               Brier_weights$Brier,
               Brier_weights$Brier_per_model[1],
               Brier_weights$Brier_per_model[2],
               Brier_weights$Brier_per_model[3],
               Brier_weights$Brier_per_model[4]
               ),
  weights_int1 = c(NA,
                   NA,
                   NA,
                   NA,
                   Brier_weights$weights[1],
                   Brier_weights$weights[2],
                   Brier_weights$weights[3],
                   Brier_weights$weights[4]
                   ),
  ibs_int2 = c(brier_score_multi_2$Brier,
               brier_score_SL_full_2$Brier,
               Brier_multi_CV_2$Brier,
               Brier_weights_2$Brier,
               Brier_weights_2$Brier_per_model[1],
               Brier_weights_2$Brier_per_model[2],
               Brier_weights_2$Brier_per_model[3],
               Brier_weights_2$Brier_per_model[4]
               ),
  weights_int2 = c(NA,
                   NA,
                   NA,
                   NA,
                   Brier_weights_2$weights[1],
                   Brier_weights_2$weights[2],
                   Brier_weights_2$weights[3],
                   Brier_weights_2$weights[4]
  )
)

rownames(data_ibs_results) <- c("mJM_full",
                                "SL_full",
                                "mJM",
                                "SL",
                                "M1",
                                "M2",
                                "M3",
                                "M4")

## Table for the EPCE:
data_epce_results <- data.frame(
  epce_int1 = c(EPCE_score_multi$EPCE,
               EPCE_score_SL_full$EPCE,
               EPCE_multi_CV$EPCE,
               EPCE_weights$EPCE,
               EPCE_weights$EPCE_per_model[1],
               EPCE_weights$EPCE_per_model[2],
               EPCE_weights$EPCE_per_model[3],
               EPCE_weights$EPCE_per_model[4]
  ),
  weights_int1 = c(NA,
                   NA,
                   NA,
                   NA,
                   EPCE_weights$weights[1],
                   EPCE_weights$weights[2],
                   EPCE_weights$weights[3],
                   EPCE_weights$weights[4]
  ),
  epce_int2 = c(EPCE_score_multi_2$EPCE,
                EPCE_score_SL_full_2$EPCE,
                EPCE_multi_CV_2$EPCE,
                EPCE_weights_2$EPCE,
                EPCE_weights_2$EPCE_per_model[1],
                EPCE_weights_2$EPCE_per_model[2],
                EPCE_weights_2$EPCE_per_model[3],
                EPCE_weights_2$EPCE_per_model[4]
  ),
  weights_int2 = c(NA,
                   NA,
                   NA,
                   NA,
                   EPCE_weights_2$weights[1],
                   EPCE_weights_2$weights[2],
                   EPCE_weights_2$weights[3],
                   EPCE_weights_2$weights[4]
  )
)

rownames(data_epce_results) <- c("mJM_full",
                                "SL_full",
                                "mJM",
                                "SL",
                                "M1",
                                "M2",
                                "M3",
                                "M4")

save(data_ibs_results, data_epce_results,
     file = "/Users/arnaugarcia/Desktop/TFM_super/PBC_analysis/results/results1.RData")



