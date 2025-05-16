##############################################################################
##############################################################################
## Plotting results for SL simulation study
##############################################################################
##############################################################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

## - We define a function that generates boxplots for the results obtained
##.  for both metrics, IBS and EPCE
plot_boxs <- function(string, sce_cens, n = 100,
                      int1, int2,
                      lims_ibs_1 = c(NA, NA),
                      lims_ibs_2 = c(NA, NA), 
                      lims_epce_1 = c(NA, NA),
                      lims_epce_2 = c(NA, NA)){
  #Load the RData file with results
  load(string)

  #########################################
  #### Plots for the interval (4,5.5]
  #########################################
  df_comparing_ibs <- data.frame(IBS_mJM_train = ibs_train[1:n],
                                 IBS_mJM_test = ibs_test[1:n],
                                 IBS_eSL_train = eSL_cv_IBS[1:n],
                                 IBS_eSL_test = eSL_test_IBS[1:n],
                                 IBS_dSL_train = dSL_cv_IBS[1:n],
                                 IBS_dSL_test = dSL_test_IBS[1:n],
                                 perc_cens_train = perc_cens_train[1:n],
                                 perc_cens_test = perc_cens_test[1:n],
                                 n_event_train = n_event_train[1:n],
                                 n_event_test = n_event_test[1:n],
                                 n_risk_train = n_risk_train[1:n],
                                 n_risk_test = n_risk_test[1:n])
  
  
  
  df_comparing_epce <- data.frame(EPCE_mJM_train = epce_train[1:n],
                                  EPCE_mJM_test = epce_test[1:n],
                                  EPCE_eSL_train = eSL_cv_EPCE[1:n],
                                  EPCE_eSL_test = eSL_test_EPCE[1:n],
                                  EPCE_dSL_train = dSL_cv_EPCE[1:n],
                                  EPCE_dSL_test = dSL_test_EPCE[1:n],
                                  perc_cens_train = perc_cens_train[1:n],
                                  perc_cens_test = perc_cens_test[1:n],
                                  n_event_train = n_event_train[1:n],
                                  n_event_test = n_event_test[1:n],
                                  n_risk_train = n_risk_train[1:n],
                                  n_risk_test = n_risk_test[1:n])
  
  
  #boxplot IBS
  df <- with(df_comparing_ibs, data.frame("dSL CV Data" = IBS_dSL_train,
                                          "eSL CV Data" = IBS_eSL_train,
                                          "Multivariate model CV Data" = IBS_mJM_train,
                                          "dSL Test Data" = IBS_dSL_test,
                                          "eSL Test Data" = IBS_eSL_test,
                                          "Multivariate model Test Data" = IBS_mJM_test))
  df_long <- df %>%
    pivot_longer(cols = everything(), names_to = "variable", values_to = "value")
  
  titl_ibs <- paste0("IBS; ", int1, ";\n", sce_cens)
  
  df_long$variable <- factor(df_long$variable,
                             levels = c("dSL.CV.Data",
                                        "eSL.CV.Data",
                                        "Multivariate.model.CV.Data",
                                        "dSL.Test.Data",
                                        "eSL.Test.Data",
                                        "Multivariate.model.Test.Data"))
  
  plot_ibs <- ggplot(df_long, aes(x = variable, y = value)) +
    geom_boxplot(fill = "lightblue") +
    theme_minimal() +
    theme(axis.title = element_text(size = 6),
          axis.text = element_text(size = 7.5),
          plot.title = element_text(size = 8)) +
    xlab("") + ylab("IBS") +
    ylim(lims_ibs_1) + 
    scale_x_discrete(labels = c("dSL.CV.Data" = "dSL\n CV Data",
                                "eSL.CV.Data" = "eSL\n CV Data",
                                "Multivariate.model.CV.Data" = "mJM\n CV Data",
                                "dSL.Test.Data" = "dSL\n Test Data",
                                "eSL.Test.Data" = "eSL\n Test Data",
                                "Multivariate.model.Test.Data" = "mJM\n Test Data"))+
    geom_vline(xintercept = (3+4)/2, linetype = "dashed", color = "grey", linewidth = 0.5)+
    ggtitle(titl_ibs)
  
  
  ## EPCE plot
  df2 <- with(df_comparing_epce, data.frame("dSL CV Data" = EPCE_dSL_train,
                                            "eSL CV Data" = EPCE_eSL_train,
                                            "Multivariate model CV Data" = EPCE_mJM_train,
                                            "dSL Test Data" = EPCE_dSL_test,
                                            "eSL Test Data" = EPCE_eSL_test,
                                            "Multivariate model Test Data" = EPCE_mJM_test))
  df_long2 <- df2 %>%
    pivot_longer(cols = everything(), names_to = "variable", values_to = "value")
  
  df_long2$variable <- factor(df_long2$variable,
                              levels = c("dSL.CV.Data",
                                         "eSL.CV.Data",
                                         "Multivariate.model.CV.Data",
                                         "dSL.Test.Data",
                                         "eSL.Test.Data",
                                         "Multivariate.model.Test.Data"))
  
  titl_epce <- paste0("EPCE; ", int1, ";\n", sce_cens)
  
  plot_epce <- ggplot(df_long2, aes(x = variable, y = value)) +
    geom_boxplot(fill = "lightblue") +
    theme_minimal() +
    theme(axis.title = element_text(size = 6),
          axis.text = element_text(size = 7.5),
          plot.title = element_text(size = 8)) +
    xlab("") + ylab("EPCE") +
    ylim(lims_epce_1) + 
    scale_x_discrete(labels = c("dSL.CV.Data" = "dSL\n CV Data",
                                "eSL.CV.Data" = "eSL\n CV Data",
                                "Multivariate.model.CV.Data" = "mJM\n CV Data",
                                "dSL.Test.Data" = "dSL\n Test Data",
                                "eSL.Test.Data" = "eSL\n Test Data",
                                "Multivariate.model.Test.Data" = "mJM\n Test Data"))+
    geom_vline(xintercept = (3+4)/2, linetype = "dashed", color = "grey", linewidth = 0.5)+
    ggtitle(titl_epce)
  
  #########################################
  #### Plots for the interval (4,5.5]
  #########################################
  df_comparing_ibs_2 <- data.frame(IBS_mJM_train = ibs_train_2[1:n],
                                 IBS_mJM_test = ibs_test_2[1:n],
                                 IBS_eSL_train = eSL_cv_IBS_2[1:n],
                                 IBS_eSL_test = eSL_test_IBS_2[1:n],
                                 IBS_dSL_train = dSL_cv_IBS_2[1:n],
                                 IBS_dSL_test = dSL_test_IBS_2[1:n],
                                 perc_cens_train = perc_cens_train[1:n],
                                 perc_cens_test = perc_cens_test[1:n],
                                 n_event_train = n_event_train_2[1:n],
                                 n_event_test = n_event_test_2[1:n],
                                 n_risk_train = n_risk_train_2[1:n],
                                 n_risk_test = n_risk_test_2[1:n])
  
  
  
  df_comparing_epce_2 <- data.frame(EPCE_mJM_train = epce_train_2[1:n],
                                  EPCE_mJM_test = epce_test_2[1:n],
                                  EPCE_eSL_train = eSL_cv_EPCE_2[1:n],
                                  EPCE_eSL_test = eSL_test_EPCE_2[1:n],
                                  EPCE_dSL_train = dSL_cv_EPCE_2[1:n],
                                  EPCE_dSL_test = dSL_test_EPCE_2[1:n],
                                  perc_cens_train = perc_cens_train[1:n],
                                  perc_cens_test = perc_cens_test[1:n],
                                  n_event_train = n_event_train_2[1:n],
                                  n_event_test = n_event_test_2[1:n],
                                  n_risk_train = n_risk_train_2[1:n],
                                  n_risk_test = n_risk_test_2[1:n])
  
  
  #boxplot IBS
  df_2 <- with(df_comparing_ibs_2, data.frame("dSL CV Data" = IBS_dSL_train,
                                          "eSL CV Data" = IBS_eSL_train,
                                          "Multivariate model CV Data" = IBS_mJM_train,
                                          "dSL Test Data" = IBS_dSL_test,
                                          "eSL Test Data" = IBS_eSL_test,
                                          "Multivariate model Test Data" = IBS_mJM_test))
  df_long <- df_2 %>%
    pivot_longer(cols = everything(), names_to = "variable", values_to = "value")
  
  titl_ibs <- paste0("IBS; ", int2, ";\n", sce_cens)
  
  df_long$variable <- factor(df_long$variable,
                             levels = c("dSL.CV.Data",
                                        "eSL.CV.Data",
                                        "Multivariate.model.CV.Data",
                                        "dSL.Test.Data",
                                        "eSL.Test.Data",
                                        "Multivariate.model.Test.Data"))
  
  plot_ibs_2 <- ggplot(df_long, aes(x = variable, y = value)) +
    geom_boxplot(fill = "lightblue") +
    theme_minimal() +
    theme(axis.title = element_text(size = 6),
          axis.text = element_text(size = 7.5),
          plot.title = element_text(size = 8)) +
    xlab("") + ylab("IBS") +
    ylim(lims_ibs_2) + 
    scale_x_discrete(labels = c("dSL.CV.Data" = "dSL\n CV Data",
                                "eSL.CV.Data" = "eSL\n CV Data",
                                "Multivariate.model.CV.Data" = "mJM\n CV Data",
                                "dSL.Test.Data" = "dSL\n Test Data",
                                "eSL.Test.Data" = "eSL\n Test Data",
                                "Multivariate.model.Test.Data" = "mJM\n Test Data"))+
    geom_vline(xintercept = (3+4)/2, linetype = "dashed", color = "grey", linewidth = 0.5)+
    ggtitle(titl_ibs)
  
  
  ## EPCE plot
  df2_2 <- with(df_comparing_epce_2, data.frame("dSL CV Data" = EPCE_dSL_train,
                                            "eSL CV Data" = EPCE_eSL_train,
                                            "Multivariate model CV Data" = EPCE_mJM_train,
                                            "dSL Test Data" = EPCE_dSL_test,
                                            "eSL Test Data" = EPCE_eSL_test,
                                            "Multivariate model Test Data" = EPCE_mJM_test))
  df_long2 <- df2_2 %>%
    pivot_longer(cols = everything(), names_to = "variable", values_to = "value")
  
  df_long2$variable <- factor(df_long2$variable,
                              levels = c("dSL.CV.Data",
                                         "eSL.CV.Data",
                                         "Multivariate.model.CV.Data",
                                         "dSL.Test.Data",
                                         "eSL.Test.Data",
                                         "Multivariate.model.Test.Data"))
  
  titl_epce <- paste0("EPCE; ", int2, ";\n", sce_cens)
  
  plot_epce_2 <- ggplot(df_long2, aes(x = variable, y = value)) +
    geom_boxplot(fill = "lightblue") +
    theme_minimal() +
    theme(axis.title = element_text(size = 6),
          axis.text = element_text(size = 7.5),
          plot.title = element_text(size = 8)) +
    xlab("") + ylab("EPCE") +
    ylim(lims_epce_2) + 
    scale_x_discrete(labels = c("dSL.CV.Data" = "dSL\n CV Data",
                                "eSL.CV.Data" = "eSL\n CV Data",
                                "Multivariate.model.CV.Data" = "mJM\n CV Data",
                                "dSL.Test.Data" = "dSL\n Test Data",
                                "eSL.Test.Data" = "eSL\n Test Data",
                                "Multivariate.model.Test.Data" = "mJM\n Test Data"))+
    geom_vline(xintercept = (3+4)/2, linetype = "dashed", color = "grey", linewidth = 0.5)+
    ggtitle(titl_epce)
  
  
  return(list(plot_ibs, plot_epce,
              plot_ibs_2, plot_epce_2))
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
box_adminCens <- plot_boxs("/Users/arnaugarcia/Desktop/TFM_super/simul_SL/results/repl100_adminCens_30_19mar2025.RData",
                           sce_cens = "Administrative Censoring",
                           int1 = "(6,7.5]",
                           int2 = "(4,5.5]",
                           n = 100,
                           lims_ibs_1 = c(0.12, 0.22),
                           lims_ibs_2 = c(0.02, 0.10),
                           lims_epce_1 = c(0.6, 1.1),
                           lims_epce_2 = c(0.3, 0.6))

box_adminCens[[1]]
box_adminCens[[2]]
box_adminCens[[3]]
box_adminCens[[4]]

box_randCens <- plot_boxs("/Users/arnaugarcia/Desktop/TFM_super/simul_SL/results/repl100_randCens_30_25mar2025.RData",
                          sce_cens = "Random Censoring",
                          int1 = "(6,7.5]",
                          int2 = "(4,5.5]",
                          n = 100,
                          lims_ibs_1 = c(0.12, 0.22),
                          lims_ibs_2 = c(0.02, 0.10),
                          lims_epce_1 = c(0.6, 1.1),
                          lims_epce_2 = c(0.3, 0.6))

box_randCens[[1]]
box_randCens[[2]]
box_randCens[[3]]
box_randCens[[4]]

box_informCens <- plot_boxs("/Users/arnaugarcia/Desktop/TFM_super/simul_SL/results/repl100_informativeCens_30_31mar2025.RData",
                          sce_cens = "Informative Censoring",
                          int1 = "(6,7.5]",
                          int2 = "(4,5.5]",
                          n = 100,
                          lims_ibs_1 = c(0.12, 0.22),
                          lims_ibs_2 = c(0.02, 0.10),
                          lims_epce_1 = c(0.6, 1.1),
                          lims_epce_2 = c(0.3, 0.6))

box_informCens[[1]]
box_informCens[[2]]
box_informCens[[3]]
box_informCens[[4]]

### Results 60% censoring
box_adminCens_60 <- plot_boxs("/Users/arnaugarcia/Desktop/TFM_super/simul_SL/results/repl50_adminCens_60_07apr2025.RData",
                             sce_cens = "Administrative Censoring; 60%",
                             int1 = "(6,7.5]",
                             int2 = "(4,5.5]",
                             n = 50,
                             lims_ibs_1 = c(0, 0.1),
                             lims_ibs_2 = c(0.02, 0.10),
                             lims_epce_1 = c(0.2, 0.5),
                             lims_epce_2 = c(0.2, 0.65))

box_adminCens_60[[1]]
box_adminCens_60[[2]]
box_adminCens_60[[3]]
box_adminCens_60[[4]]

box_randCens_60 <- plot_boxs("/Users/arnaugarcia/Desktop/TFM_super/simul_SL/results/repl100_randCens_60_07apr2025.RData",
                            sce_cens = "Random Censoring; 60%",
                            int1 = "(6,7.5]",
                            int2 = "(4,5.5]",
                            n = 100,
                            lims_ibs_1 = c(0.12, 0.22),
                            lims_ibs_2 = c(0.02, 0.10),
                            lims_epce_1 = c(0.5, 1.05),
                            lims_epce_2 = c(0.2, 0.65))

box_randCens_60[[1]]
box_randCens_60[[2]]
box_randCens_60[[3]]
box_randCens_60[[4]]

box_informCens_60 <- plot_boxs("/Users/arnaugarcia/Desktop/TFM_super/simul_SL/results/repl100_informativeCens_60_07apr2025.RData",
                            sce_cens = "Informative Censoring; 60%",
                            int1 = "(6,7.5]",
                            int2 = "(4,5.5]",
                            n = 100,
                            lims_ibs_1 = c(0.12, 0.22),
                            lims_ibs_2 = c(0.02, 0.10),
                            lims_epce_1 = c(0.5, 1.05),
                            lims_epce_2 = c(0.2, 0.65))

box_informCens_60[[1]]
box_informCens_60[[2]]
box_informCens_60[[3]]
box_informCens_60[[4]]


## setting lines to separate the boxplots

separator1 <- ggplot() + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1) + 
  theme_void() 

separator2 <- ggplot() + 
  geom_hline(yintercept = 0, color = "black", linewidth = 1) + 
  theme_void() 

### Plot tfm:

r1 <- ggarrange(box_adminCens[[1]], box_randCens[[1]], box_informCens[[1]],
                nrow = 1, ncol = 3)

r2 <- ggarrange(box_adminCens[[2]], box_randCens[[2]], box_informCens[[2]],
                nrow = 1, ncol = 3)

r3 <- ggarrange(box_adminCens[[3]], box_randCens[[3]], box_informCens[[3]],
                nrow = 1, ncol = 3)

r4 <- ggarrange(box_adminCens[[4]], box_randCens[[4]], box_informCens[[4]],
                nrow = 1, ncol = 3)

## first plot with separator
big_r1 <- ggarrange(r1,
                    separator1,
                    r2,
                    ncol = 1,
                    heights = c(2, 0.1, 2))

big_r2 <- ggarrange(r3,
                    separator1,
                    r4,
                    ncol = 1,
                    heights = c(2, 0.1, 2))

final_plot <- ggarrange(big_r2,
                        separator2,
                        big_r1,
                        ncol = 1,
                        heights = c(2, 0.1, 2))


# w = 12, h = 15
pdf("/Users/arnaugarcia/Desktop/TFM_super/simul_SL/plots/plot_SL_simul_12x15_def.pdf",
    width = 12, height = 15)
final_plot
dev.off()

#######################
## Plot tfm 60%:

r1 <- ggarrange(box_adminCens_60[[1]], box_randCens_60[[1]], box_informCens_60[[1]],
                nrow = 1, ncol = 3)

r2 <- ggarrange(box_adminCens_60[[2]], box_randCens_60[[2]], box_informCens_60[[2]],
                nrow = 1, ncol = 3)

r3 <- ggarrange(box_adminCens_60[[3]], box_randCens_60[[3]], box_informCens_60[[3]],
                nrow = 1, ncol = 3)

r4 <- ggarrange(box_adminCens_60[[4]], box_randCens_60[[4]], box_informCens_60[[4]],
                nrow = 1, ncol = 3)

## first plot with separator
big_r1 <- ggarrange(r1,
                    separator1,
                    r2,
                    ncol = 1,
                    heights = c(2, 0.1, 2))

big_r2 <- ggarrange(r3,
                    separator1,
                    r4,
                    ncol = 1,
                    heights = c(2, 0.1, 2))

final_plot <- ggarrange(big_r2,
                        separator2,
                        big_r1,
                        ncol = 1,
                        heights = c(2, 0.1, 2))

# w = 12, h = 15
pdf("/Users/arnaugarcia/Desktop/TFM_super/simul_SL/plots/plot_SL_simul_60Cens_12x15_def.pdf",
    width = 12, height = 15)
final_plot
dev.off()


##############################################
### plot extended abstract CEB:
##############################################

ploti1 <- ggarrange(box_adminCens[[1]], box_adminCens[[3]],
                    box_randCens[[1]], box_randCens[[3]],
                    nrow = 1, ncol = 4)

ploti2 <- ggarrange(box_adminCens[[2]], box_adminCens[[4]],
                    box_randCens[[2]], box_randCens[[4]],
                    nrow = 1, ncol = 4)

ploti_ceb <- ggarrange(ploti1,
                       separator1,
                       ploti2,
                       ncol = 1,
                       heights = c(2, 0.1, 2))

ggarrange(box_adminCens[[1]], box_adminCens[[3]],
          box_randCens[[1]], box_randCens[[3]],
          box_adminCens[[2]], box_adminCens[[4]],
          box_randCens[[2]], box_randCens[[4]],
          nrow = 1, ncol = 8)


pdf("/Users/arnaugarcia/Desktop/TFM_super/simul_SL/plot_ceb4.pdf",width = 10, height = 5)
ploti_ceb
dev.off()


### OBS: we have to compute indiv at risk before the point on avg, etc...
###.     so we should compute all the relevant metrics on avg

load("/Users/arnaugarcia/Desktop/TFM_super/simul_SL/results/repl50_adminCens_60_07apr2025.RData")

n<-50
df_comparing_ibs <- data.frame(IBS_mJM_train = ibs_train[1:n],
                               IBS_mJM_test = ibs_test[1:n],
                               IBS_eSL_train = eSL_cv_IBS[1:n],
                               IBS_eSL_test = eSL_test_IBS[1:n],
                               IBS_dSL_train = dSL_cv_IBS[1:n],
                               IBS_dSL_test = dSL_test_IBS[1:n],
                               perc_cens_train = perc_cens_train[1:n],
                               perc_cens_test = perc_cens_test[1:n],
                               n_event_train = n_event_train[1:n],
                               n_event_test = n_event_test[1:n],
                               n_risk_train = n_risk_train[1:n],
                               n_risk_test = n_risk_test[1:n])

df_comparing_ibs2 <- data.frame(IBS_mJM_train = ibs_train_2[1:n],
                               IBS_mJM_test = ibs_test_2[1:n],
                               IBS_eSL_train = eSL_cv_IBS_2[1:n],
                               IBS_eSL_test = eSL_test_IBS_2[1:n],
                               IBS_dSL_train = dSL_cv_IBS_2[1:n],
                               IBS_dSL_test = dSL_test_IBS_2[1:n],
                               perc_cens_train = perc_cens_train[1:n],
                               perc_cens_test = perc_cens_test[1:n],
                               n_event_train = n_event_train_2[1:n],
                               n_event_test = n_event_test_2[1:n],
                               n_risk_train = n_risk_train_2[1:n],
                               n_risk_test = n_risk_test_2[1:n])

df_comparing_epce <- data.frame(EPCE_mJM_train = epce_train[1:n],
                                EPCE_mJM_test = epce_test[1:n],
                                EPCE_eSL_train = eSL_cv_EPCE[1:n],
                                EPCE_eSL_test = eSL_test_EPCE[1:n],
                                EPCE_dSL_train = dSL_cv_EPCE[1:n],
                                EPCE_dSL_test = dSL_test_EPCE[1:n],
                                perc_cens_train = perc_cens_train[1:n],
                                perc_cens_test = perc_cens_test[1:n],
                                n_event_train = n_event_train[1:n],
                                n_event_test = n_event_test[1:n],
                                n_risk_train = n_risk_train[1:n],
                                n_risk_test = n_risk_test[1:n])

df_comparing_epce2 <- data.frame(EPCE_mJM_train = epce_train_2[1:n],
                                EPCE_mJM_test = epce_test_2[1:n],
                                EPCE_eSL_train = eSL_cv_EPCE_2[1:n],
                                EPCE_eSL_test = eSL_test_EPCE_2[1:n],
                                EPCE_dSL_train = dSL_cv_EPCE_2[1:n],
                                EPCE_dSL_test = dSL_test_EPCE_2[1:n],
                                perc_cens_train = perc_cens_train[1:n],
                                perc_cens_test = perc_cens_test[1:n],
                                n_event_train = n_event_train_2[1:n],
                                n_event_test = n_event_test_2[1:n],
                                n_risk_train = n_risk_train_2[1:n],
                                n_risk_test = n_risk_test_2[1:n])


##############################################
## Computing summary statistics on the weights
## obtained for SL in different scenarios
##############################################

############################
##### Admin censoring, 30%
############################
load("/Users/arnaugarcia/Desktop/TFM_super/simul_SL/results/repl100_adminCens_30_19mar2025.RData")

## We build a data frame with the weights per model:

## IBS interval (6,7.5]
w_ibs_int1 <- as.data.frame(IBS_w)
names(w_ibs_int1) <- c("M1", "M2", "M3", "M4")
summary(w_ibs_int1)

## IBS interval (4,5.5]
w_ibs_int2 <- as.data.frame(IBS_w_2)
names(w_ibs_int2) <- c("M1", "M2", "M3", "M4")
summary(w_ibs_int2)

## EPCE interval (6,7.5]
w_epce_int1 <- as.data.frame(EPCE_w)
names(w_epce_int1) <- c("M1", "M2", "M3", "M4")
summary(w_epce_int1)

## EPCE interval (4,5.5]
w_epce_int2 <- as.data.frame(EPCE_w_2)
names(w_epce_int2) <- c("M1", "M2", "M3", "M4")
summary(w_epce_int2)

############################
##### Random censoring, 30%
############################
load("/Users/arnaugarcia/Desktop/TFM_super/simul_SL/results/repl100_randCens_30_25mar2025.RData")

## We build a data frame with the weights per model:

## IBS interval (6,7.5]
w_ibs_int1 <- as.data.frame(IBS_w)
names(w_ibs_int1) <- c("M1", "M2", "M3", "M4")
summary(w_ibs_int1)

## IBS interval (4,5.5]
w_ibs_int2 <- as.data.frame(IBS_w_2)
names(w_ibs_int2) <- c("M1", "M2", "M3", "M4")
summary(w_ibs_int2)

## EPCE interval (6,7.5]
w_epce_int1 <- as.data.frame(EPCE_w)
names(w_epce_int1) <- c("M1", "M2", "M3", "M4")
summary(w_epce_int1)

## EPCE interval (4,5.5]
w_epce_int2 <- as.data.frame(EPCE_w_2)
names(w_epce_int2) <- c("M1", "M2", "M3", "M4")
summary(w_epce_int2)

############################
##### Informative censoring, 30%
############################
load("/Users/arnaugarcia/Desktop/TFM_super/simul_SL/results/repl100_informativeCens_30_31mar2025.RData")

## We build a data frame with the weights per model:

## IBS interval (6,7.5]
w_ibs_int1 <- as.data.frame(IBS_w)
names(w_ibs_int1) <- c("M1", "M2", "M3", "M4")
summary(w_ibs_int1)

## IBS interval (4,5.5]
w_ibs_int2 <- as.data.frame(IBS_w_2)
names(w_ibs_int2) <- c("M1", "M2", "M3", "M4")
summary(w_ibs_int2)

## EPCE interval (6,7.5]
w_epce_int1 <- as.data.frame(EPCE_w)
names(w_epce_int1) <- c("M1", "M2", "M3", "M4")
summary(w_epce_int1)

## EPCE interval (4,5.5]
w_epce_int2 <- as.data.frame(EPCE_w_2)
names(w_epce_int2) <- c("M1", "M2", "M3", "M4")
summary(w_epce_int2)

######################################################
## Computing other metrics:
## - % of censoring (on avg) per scenario
## - dSL selected per scenario
######################################################

## Admin cens 30%
load("/Users/arnaugarcia/Desktop/TFM_super/simul_SL/results/repl100_adminCens_30_19mar2025.RData")
mean(perc_cens_train)
mean(perc_cens_test)
table(disSL_ibs_2)
table(disSL_ibs)
table(disSL_epce_2)
table(disSL_epce)

## Random cens 30%
load("/Users/arnaugarcia/Desktop/TFM_super/simul_SL/results/repl100_randCens_30_25mar2025.RData")
mean(perc_cens_train)
mean(perc_cens_test)
table(disSL_ibs_2)
table(disSL_ibs)
table(disSL_epce_2)
table(disSL_epce)

## Info cens 30%
load("/Users/arnaugarcia/Desktop/TFM_super/simul_SL/results/repl100_informativeCens_30_31mar2025.RData")
mean(perc_cens_train)
mean(perc_cens_test)
table(disSL_ibs_2)
table(disSL_ibs)
table(disSL_epce_2)
table(disSL_epce)

## Admin cens 60%
load("/Users/arnaugarcia/Desktop/TFM_super/simul_SL/results/repl75_adminCens_60_07apr2025.RData")
mean(perc_cens_train[1:75])
mean(perc_cens_test[1:75])
table(disSL_ibs_2)
table(disSL_ibs)
table(disSL_epce_2)
table(disSL_epce)

load("/Users/arnaugarcia/Desktop/TFM_super/simul_SL/results/repl75to100_adminCens_60_29apr2025.RData")
mean(perc_cens_train)
mean(perc_cens_test)
table(disSL_ibs_2)
table(disSL_ibs)
table(disSL_epce_2)
table(disSL_epce)





## Random cens 60%
load("/Users/arnaugarcia/Desktop/TFM_super/simul_SL/results/repl100_randCens_60_07apr2025.RData")
mean(perc_cens_train)
mean(perc_cens_test)
table(disSL_ibs_2)
table(disSL_ibs)
table(disSL_epce_2)
table(disSL_epce)


## Info cens 60%
load("/Users/arnaugarcia/Desktop/TFM_super/simul_SL/results/repl100_informativeCens_60_07apr2025.RData")
mean(perc_cens_train)
mean(perc_cens_test)
table(disSL_ibs_2)
table(disSL_ibs)
table(disSL_epce_2)
table(disSL_epce)




