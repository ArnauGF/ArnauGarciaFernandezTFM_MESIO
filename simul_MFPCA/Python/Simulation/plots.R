#############################################################################
##############################################################################
## Plotting results for MFPCA simulation study
##############################################################################
##############################################################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

## - We define a function that generates boxplots for the results obtained
##.  for both metrics, IBS and EPCE
plot_boxs <- function(string, string_cox, string_deep, sce_cens, n = 100,
                      int1, int2,
                      lims_ibs_1 = c(NA, NA),
                      lims_ibs_2 = c(NA, NA)){
  #Load the RData file with results
  load(string)
  dat_cox <- read.csv(string_cox)
  dat_deep <- read.csv(string_deep)
  
  #########################################
  #### Plots for the interval1
  #########################################
  df_comparing_ibs <- data.frame(IBS_mJM_train = ibs_train[1:n],
                                 IBS_mJM_test = ibs_test[1:n],
                                 IBS_eSL_train = eSL_cv_IBS[1:n],
                                 IBS_eSL_test = eSL_test_IBS[1:n],
                                 IBS_dSL_train = dSL_cv_IBS[1:n],
                                 IBS_dSL_test = dSL_test_IBS[1:n],
                                 IBS_mfpca_cox_train = dat_cox$IBS_train6_7[1:n],
                                 IBS_mfpca_cox_test = dat_cox$IBS_test6_7[1:n],
                                 IBS_mfpca_deepsurv_train = dat_deep$IBS_train6_7[1:n],
                                 IBS_mfpca_deepsurv_test = dat_deep$IBS_test6_7[1:n],
                                 perc_cens_train = perc_cens_train[1:n],
                                 perc_cens_test = perc_cens_test[1:n],
                                 n_event_train = n_event_train[1:n],
                                 n_event_test = n_event_test[1:n],
                                 n_risk_train = n_risk_train[1:n],
                                 n_risk_test = n_risk_test[1:n])
  
  
  #boxplot IBS
  df <- with(df_comparing_ibs, data.frame("MFPCA-Cox Train Data" = IBS_mfpca_cox_train,
                                          "MFPCA-DeepSurv Train Data" = IBS_mfpca_deepsurv_train,
                                          "dSL CV Data" = IBS_dSL_train,
                                          "eSL CV Data" = IBS_eSL_train,
                                          "Multivariate model CV Data" = IBS_mJM_train,
                                          "MFPCA-Cox Test Data" = IBS_mfpca_cox_test,
                                          "MFPCA-DeepSurv Test Data" = IBS_mfpca_deepsurv_test,
                                          "dSL Test Data" = IBS_dSL_test,
                                          "eSL Test Data" = IBS_eSL_test,
                                          "Multivariate model Test Data" = IBS_mJM_test))
  df_long <- df %>%
    pivot_longer(cols = everything(), names_to = "variable", values_to = "value")
  
  titl_ibs <- paste0("IBS; ", int1, ";\n", sce_cens)
  
  df_long$variable <- factor(df_long$variable,
                             levels = c("MFPCA.Cox.Train.Data",
                                        "MFPCA.DeepSurv.Train.Data",
                                        "dSL.CV.Data",
                                        "eSL.CV.Data",
                                        "Multivariate.model.CV.Data",
                                        "MFPCA.Cox.Test.Data",
                                        "MFPCA.DeepSurv.Test.Data",
                                        "dSL.Test.Data",
                                        "eSL.Test.Data",
                                        "Multivariate.model.Test.Data"))
  
  plot_ibs <- ggplot(df_long, aes(x = variable, y = value, fill = variable)) +
    geom_boxplot() +
    scale_fill_manual(values = c(
      "MFPCA.Cox.Train.Data" = "lightgreen",
      "MFPCA.DeepSurv.Train.Data" = "coral3",
      "dSL.CV.Data" = "lightblue",
      "eSL.CV.Data" = "lightblue",
      "Multivariate.model.CV.Data" = "lightblue",
      "MFPCA.Cox.Test.Data" = "lightgreen",
      "MFPCA.DeepSurv.Test.Data" = "coral3",
      "dSL.Test.Data" = "lightblue",
      "eSL.Test.Data" = "lightblue",
      "Multivariate.model.Test.Data" = "lightblue"
    )) +
    theme_minimal() +
    theme(axis.title = element_text(size = 6),
          axis.text = element_text(size = 5),
          plot.title = element_text(size = 8),
          legend.position = "none") +
    xlab("") + ylab("IBS") +
    ylim(lims_ibs_1) + 
    scale_x_discrete(labels = c("MFPCA.Cox.Train.Data" = "MFPCA\n Cox \n Train Data",
                                "MFPCA.DeepSurv.Train.Data" = "MFPCA\n DeepSurv \n Train Data",
                                "dSL.CV.Data" = "dSL\n CV Data",
                                "eSL.CV.Data" = "eSL\n CV Data",
                                "Multivariate.model.CV.Data" = "mJM\n CV Data",
                                "MFPCA.Cox.Test.Data" = "MFPCA\n Cox \n Test Data",
                                "MFPCA.DeepSurv.Test.Data" = "MFPCA\n DeepSurv \n Test Data",
                                "dSL.Test.Data" = "dSL\n Test Data",
                                "eSL.Test.Data" = "eSL\n Test Data",
                                "Multivariate.model.Test.Data" = "mJM\n Test Data"))+
    geom_vline(xintercept = 5.5, linetype = "dashed", color = "grey", linewidth = 0.5)+
    ggtitle(titl_ibs)
  
  
  #########################################
  #### Plots for the interval (4,5.5]
  #########################################
  df_comparing_ibs_2 <- data.frame(IBS_mJM_train = ibs_train_2[1:n],
                                   IBS_mJM_test = ibs_test_2[1:n],
                                   IBS_eSL_train = eSL_cv_IBS_2[1:n],
                                   IBS_eSL_test = eSL_test_IBS_2[1:n],
                                   IBS_dSL_train = dSL_cv_IBS_2[1:n],
                                   IBS_dSL_test = dSL_test_IBS_2[1:n],
                                   IBS_mfpca_cox_train = dat_cox$IBS_train[1:n],
                                   IBS_mfpca_cox_test = dat_cox$IBS_test[1:n],
                                   IBS_mfpca_deepsurv_train = dat_deep$IBS_train[1:n],
                                   IBS_mfpca_deepsurv_test = dat_deep$IBS_test[1:n],
                                   perc_cens_train = perc_cens_train[1:n],
                                   perc_cens_test = perc_cens_test[1:n],
                                   n_event_train = n_event_train_2[1:n],
                                   n_event_test = n_event_test_2[1:n],
                                   n_risk_train = n_risk_train_2[1:n],
                                   n_risk_test = n_risk_test_2[1:n])
  
  
  
  #boxplot IBS
  df_2 <- with(df_comparing_ibs_2, data.frame("MFPCA-Cox Train Data" = IBS_mfpca_cox_train,
                                              "MFPCA-DeepSurv Train Data" = IBS_mfpca_deepsurv_train,
                                              "dSL CV Data" = IBS_dSL_train,
                                              "eSL CV Data" = IBS_eSL_train,
                                              "Multivariate model CV Data" = IBS_mJM_train,
                                              "MFPCA-Cox Test Data" = IBS_mfpca_cox_test,
                                              "MFPCA-DeepSurv Test Data" = IBS_mfpca_deepsurv_test,
                                              "dSL Test Data" = IBS_dSL_test,
                                              "eSL Test Data" = IBS_eSL_test,
                                              "Multivariate model Test Data" = IBS_mJM_test))
  df_long <- df_2 %>%
    pivot_longer(cols = everything(), names_to = "variable", values_to = "value")
  
  titl_ibs <- paste0("IBS; ", int2, ";\n", sce_cens)
  
  df_long$variable <- factor(df_long$variable,
                             levels = c("MFPCA.Cox.Train.Data",
                                        "MFPCA.DeepSurv.Train.Data",
                                        "dSL.CV.Data",
                                        "eSL.CV.Data",
                                        "Multivariate.model.CV.Data",
                                        "MFPCA.Cox.Test.Data",
                                        "MFPCA.DeepSurv.Test.Data",
                                        "dSL.Test.Data",
                                        "eSL.Test.Data",
                                        "Multivariate.model.Test.Data"))
  
  plot_ibs_2 <- ggplot(df_long, aes(x = variable, y = value, fill = variable)) +
    geom_boxplot() +
    scale_fill_manual(values = c(
      "MFPCA.Cox.Train.Data" = "lightgreen",
      "MFPCA.DeepSurv.Train.Data" = "coral3",
      "dSL.CV.Data" = "lightblue",
      "eSL.CV.Data" = "lightblue",
      "Multivariate.model.CV.Data" = "lightblue",
      "MFPCA.Cox.Test.Data" = "lightgreen",
      "MFPCA.DeepSurv.Test.Data" = "coral3",
      "dSL.Test.Data" = "lightblue",
      "eSL.Test.Data" = "lightblue",
      "Multivariate.model.Test.Data" = "lightblue"
    )) +
    theme_minimal() +
    theme(axis.title = element_text(size = 6),
          axis.text = element_text(size = 5),
          plot.title = element_text(size = 8),
          legend.position = "none") +
    xlab("") + ylab("IBS") +
    ylim(lims_ibs_2) + 
    scale_x_discrete(labels = c("MFPCA.Cox.Train.Data" = "MFPCA\n Cox \n Train Data",
                                "MFPCA.DeepSurv.Train.Data" = "MFPCA\n DeepSurv \n Train Data",
                                "dSL.CV.Data" = "dSL\n CV Data",
                                "eSL.CV.Data" = "eSL\n CV Data",
                                "Multivariate.model.CV.Data" = "mJM\n CV Data",
                                "MFPCA.Cox.Test.Data" = "MFPCA\n Cox \n Test Data",
                                "MFPCA.DeepSurv.Test.Data" = "MFPCA\n DeepSurv \n Test Data",
                                "dSL.Test.Data" = "dSL\n Test Data",
                                "eSL.Test.Data" = "eSL\n Test Data",
                                "Multivariate.model.Test.Data" = "mJM\n Test Data"))+
    geom_vline(xintercept = 5.5, linetype = "dashed", color = "grey", linewidth = 0.5)+
    ggtitle(titl_ibs)
  
  
  return(list(plot_ibs,
              plot_ibs_2))
}


## an expl on how to use the function above with preliminar results:
box_randCens <- plot_boxs("/Users/arnaugarcia/Desktop/TFM_super/results_SL/ibsipcw_epce_simul_4outcomes_correl_randCens_SL_08mar2025.RData",
                          "(6, 7.5]; Random Censoring")
box_infoCens <- plot_boxs("/Users/arnaugarcia/Desktop/TFM_super/results_SL/ibsipcw_epce_simul_4outcomes_correl_infoCens_SL_08mar2025.RData",
                          "(6, 7.5]; Informative Censoring")

box_adminCens <- plot_boxs("/Users/arnaugarcia/Desktop/TFM_super/results_SL/ibsipcw_epce_simul_4outcomes_correl_typ1_SL_06mar2025.RData",
                           "(4,5.5]; Admin Censoring")

## we can use ggarrenge() function to put boxplots together:
library(ggpubr)

# Create a blank plot with a red dashed horizontal line (simulating a separator between rows)
separator <- ggplot() + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1) + 
  theme_void() 

boxs_ibs <- ggarrange(box_adminCens[[1]], box_randCens[[1]], box_infoCens[[1]],
                      nrow = 1, ncol = 3)

boxs_epce <- ggarrange(box_adminCens[[2]], box_randCens[[2]], box_infoCens[[2]],
                       nrow = 1, ncol = 3)


boxplots_prev <- ggarrange(boxs_ibs, 
                           separator,
                           boxs_epce, 
                           ncol = 1,
                           heights = c(2, 0.1, 2))

setwd("/Users/arnaugarcia/Desktop/TFM_super")
pdf("boxplots_SL.pdf",width=12)
boxplots_prev
dev.off()



#################################################################
## Plots for def results
#################################################################
box_adminCens <- plot_boxs("/Users/arnaugarcia/Desktop/TFM_super/ML/Simulation/results/repl100_adminCens_30_2apr2025.RData",
                           string_cox = "/Users/arnaugarcia/Desktop/TFM_super/ML/Simulation/results/Cox_adminCens.csv",
                           string_deep = "/Users/arnaugarcia/Desktop/TFM_super/ML/Simulation/results/DeepSurv_adminCens_9may2025.csv",
                           sce_cens = "Administrative Censoring",
                           int1 = "(6,7.5]",
                           int2 = "(4,5.5]",
                           n = 100,
                           lims_ibs_1 = c(0.1, 0.5),
                           lims_ibs_2 = c(0.015, 0.11))

box_adminCens[[1]]
box_adminCens[[2]]

box_randCens <- plot_boxs("/Users/arnaugarcia/Desktop/TFM_super/ML/Simulation/results/repl100_randCens_30_2apr2025.RData",
                          string_cox = "/Users/arnaugarcia/Desktop/TFM_super/ML/Simulation/results/Cox_randCens.csv",
                          string_deep = "/Users/arnaugarcia/Desktop/TFM_super/ML/Simulation/results/DeepSurv_randCens_09may2025.csv",
                          sce_cens = "Random Censoring",
                          int1 = "(6,7.5]",
                          int2 = "(4,5.5]",
                          n = 100,
                          lims_ibs_1 = c(0.1, 0.5),
                          lims_ibs_2 = c(0.015, 0.11))

box_randCens[[1]]
box_randCens[[2]]


box_informCens <- plot_boxs("/Users/arnaugarcia/Desktop/TFM_super/ML/Simulation/results/repl100_informativeCens_30_2apr2025.RData",
                            string_cox = "/Users/arnaugarcia/Desktop/TFM_super/ML/Simulation/results/Cox_infCens.csv",
                            string_deep = "/Users/arnaugarcia/Desktop/TFM_super/ML/Simulation/results/DeepSurv_infCens_9may2025.csv",
                            sce_cens = "Informative Censoring",
                            int1 = "(6,7.5]",
                            int2 = "(4,5.5]",
                            n = 100,
                            lims_ibs_1 = c(0.1, 0.5),
                            lims_ibs_2 = c(0.015, 0.11))

box_informCens[[1]]
box_informCens[[2]]


## we can use ggarrenge() function to put boxplots together:
library(ggpubr)

# Create a blank plot with a red dashed horizontal line (simulating a separator between rows)
separator <- ggplot() + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1) + 
  theme_void() 

separator2 <- ggplot() + 
  geom_hline(yintercept = 0, color = "black", linewidth = 1) + 
  theme_void() 

r1 <- ggarrange(box_adminCens[[1]], box_randCens[[1]], box_informCens[[1]],
                      nrow = 1, ncol = 3)

r2 <- ggarrange(box_adminCens[[2]], box_randCens[[2]], box_informCens[[2]],
                nrow = 1, ncol = 3)

## first plot with separator
big_r1 <- ggarrange(r2,
                    separator2,
                    r1,
                    ncol = 1,
                    heights = c(2, 0.1, 2))

# w = 12, h = 15
pdf("/Users/arnaugarcia/Desktop/TFM_super/ML/Simulation/plot_MFPCA_simul_12_def.pdf",
    width = 12)
big_r1
dev.off()

###########################################
### Plots for 60% of censoring
############################################
#0.25
box_adminCens_60 <- plot_boxs("/Users/arnaugarcia/Desktop/TFM_super/ML/Simulation/results/repl100_adminCens_60_08may2025.RData",
                           string_cox = "/Users/arnaugarcia/Desktop/TFM_super/ML/Simulation/results/Cox_adminCens_60_13May2025.csv",
                           string_deep = "/Users/arnaugarcia/Desktop/TFM_super/ML/Simulation/results/DeepSurv_adminCens_60_15may2025.csv",
                           sce_cens = "Administrative Censoring; 60%",
                           int1 = "(6,6.49]",
                           int2 = "(4,5.5]",
                           n = 100,
                           lims_ibs_1 = c(0.01, 0.51),
                           lims_ibs_2 = c(0.015, 0.11))

box_adminCens_60[[1]]
box_adminCens_60[[2]]

box_randCens_60 <- plot_boxs("/Users/arnaugarcia/Desktop/TFM_super/ML/Simulation/results/repl100_randCens_60_08may2025.RData",
                          string_cox = "/Users/arnaugarcia/Desktop/TFM_super/ML/Simulation/results/Cox_randCens_60_13may.csv",
                          string_deep = "/Users/arnaugarcia/Desktop/TFM_super/ML/Simulation/results/DeepSurv_randCen_60_15may2025.csv",
                          sce_cens = "Random Censoring; 60%", 
                          int1 = "(6,7.5]",
                          int2 = "(4,5.5]",
                          n = 100,
                          lims_ibs_1 = c(0.1, 0.51),
                          lims_ibs_2 = c(0.015, 0.11))

box_randCens_60[[1]]
box_randCens_60[[2]]


box_informCens_60 <- plot_boxs("/Users/arnaugarcia/Desktop/TFM_super/ML/Simulation/results/repl100_informativeCens_60_08may2025.RData",
                            string_cox = "/Users/arnaugarcia/Desktop/TFM_super/ML/Simulation/results/Cox_infCens_60_13may2025.csv",
                            string_deep = "/Users/arnaugarcia/Desktop/TFM_super/ML/Simulation/results/DeepSurv_infCens_60_15may2025.csv",
                            sce_cens = "Informative Censoring; 60%",
                            int1 = "(6,7.5]",
                            int2 = "(4,5.5]",
                            n = 100,
                            lims_ibs_1 = c(0.1, 0.51),
                            lims_ibs_2 = c(0.015, 0.11))

box_informCens_60[[1]]
box_informCens_60[[2]]

## using ggarrenge to put plots together

# Create a blank plot with a red dashed horizontal line (simulating a separator between rows)
separator <- ggplot() + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1) + 
  theme_void() 

separator2 <- ggplot() + 
  geom_hline(yintercept = 0, color = "black", linewidth = 1) + 
  theme_void() 

r1 <- ggarrange(box_adminCens_60[[1]], box_randCens_60[[1]], box_informCens_60[[1]],
                nrow = 1, ncol = 3)

r2 <- ggarrange(box_adminCens_60[[2]], box_randCens_60[[2]], box_informCens_60[[2]],
                nrow = 1, ncol = 3)

## first plot with separator
big_r1 <- ggarrange(r2,
                    separator2,
                    r1,
                    ncol = 1,
                    heights = c(2, 0.1, 2))

# w = 12, h = 15
pdf("/Users/arnaugarcia/Desktop/TFM_super/ML/Simulation/plot_MFPCA_simul_60cens_12_def.pdf",
    width = 12)
big_r1
dev.off()


##############################################
## Computing summary statistics on the weights
## obtained for SL in different scenarios
##############################################

############################
##### Admin censoring, 30%
############################
load("/Users/arnaugarcia/Desktop/TFM_super/ML/Simulation/results/repl100_adminCens_30_2apr2025.RData")

## We build a data frame with the weights per model:

## IBS interval (6,7.5]
w_ibs_int1 <- as.data.frame(IBS_w)
names(w_ibs_int1) <- c("M1", "M2", "M3", "M4")
summary(w_ibs_int1)

## IBS interval (4,5.5]
w_ibs_int2 <- as.data.frame(IBS_w_2)
names(w_ibs_int2) <- c("M1", "M2", "M3", "M4")
summary(w_ibs_int2)


############################
##### Random censoring, 30%
############################
load("/Users/arnaugarcia/Desktop/TFM_super/ML/Simulation/results/repl100_randCens_30_2apr2025.RData")

## We build a data frame with the weights per model:

## IBS interval (6,7.5]
w_ibs_int1 <- as.data.frame(IBS_w)
names(w_ibs_int1) <- c("M1", "M2", "M3", "M4")
summary(w_ibs_int1)

## IBS interval (4,5.5]
w_ibs_int2 <- as.data.frame(IBS_w_2)
names(w_ibs_int2) <- c("M1", "M2", "M3", "M4")
summary(w_ibs_int2)

############################
##### Informative censoring, 30%
############################
load("/Users/arnaugarcia/Desktop/TFM_super/ML/Simulation/results/repl100_informativeCens_30_2apr2025.RData")

## We build a data frame with the weights per model:

## IBS interval (6,7.5]
w_ibs_int1 <- as.data.frame(IBS_w)
names(w_ibs_int1) <- c("M1", "M2", "M3", "M4")
summary(w_ibs_int1)

## IBS interval (4,5.5]
w_ibs_int2 <- as.data.frame(IBS_w_2)
names(w_ibs_int2) <- c("M1", "M2", "M3", "M4")
summary(w_ibs_int2)

######################################################
## Computing other metrics:
## - dSL selected per scenario
######################################################

## Admin cens 30%
load("/Users/arnaugarcia/Desktop/TFM_super/ML/Simulation/results/repl100_adminCens_30_2apr2025.RData")
table(disSL_ibs_2)
table(disSL_ibs)

## Random cens 30%
load("/Users/arnaugarcia/Desktop/TFM_super/ML/Simulation/results/repl100_randCens_30_2apr2025.RData")
table(disSL_ibs_2)
table(disSL_ibs)


## Info cens 30%
load("/Users/arnaugarcia/Desktop/TFM_super/ML/Simulation/results/repl100_informativeCens_30_2apr2025.RData")
table(disSL_ibs_2)
table(disSL_ibs)

## Admin cens 60%
load("/Users/arnaugarcia/Desktop/TFM_super/ML/Simulation/results/repl100_adminCens_60_08may2025.RData")
table(disSL_ibs_2)
table(disSL_ibs)

## Random cens 60%
load("/Users/arnaugarcia/Desktop/TFM_super/ML/Simulation/results/repl100_randCens_60_08may2025.RData")
table(disSL_ibs_2)
table(disSL_ibs)

## Info cens 60%
load("/Users/arnaugarcia/Desktop/TFM_super/ML/Simulation/results/repl100_informativeCens_60_08may2025.RData")
table(disSL_ibs_2)
table(disSL_ibs)
