#### Comparing FR Predictions Alternative Models Fits ####
# Using Jeffrey-Raftery's (1995) guidelines evidence in favor of the model with the smaller BIC is considered ‘weak’ for a ΔBIC of 0-2; ‘positive’ for a ΔBIC of 2-6; and ‘strong’ for a ΔBIC of 6 -10
# We compare the predictions for all models within 6 BIC units for the 'best' model

## Colipidium 
# 1. Type 2 - linear a, ΔBIC = 0
# 2. Generalized - linear a, ΔBIC = 2.5 
# 3. Type 2 - linear a & linear h, ΔBIC = 5.5 


## Dexiostoma 
# 1. Generalized - linear a, ΔBIC = 0
# 2. Generalized - linear a & linear h, ΔBIC = 2.2 
# 3. Generalized - van veen (non-linear a), ΔBIC = 3.7 


#### Visualize Fits
library(tidyverse)
library(here)
library(lme4)
library(bbmle)
library(ggpubr)

library(MASS)

library(readr)
library(readxl)
library(broom)
library(janitor)

#### Functions ####
source(here("02 - Functions", "Function-T2-a1-h1-LV.R"))
source(here("02 - Functions", "Function-T3-a1-h1-LV.R"))
source(here("02 - Functions", "Function-T3-h1-vv-LV.R"))

## Plotting Functions ###
t2a12_mod <- function(N0, M0, P, a, a12, h){ (((a + a12 * M0) * N0) / (1 + (a + a12 * M0) * h * N0)) * P}
t2a12h12_mod <- function(N0, M0, P, a, a12, h, h12){ (((a + a12 * M0) * N0) / (1 + (a + a12 * M0) * (h + h12 * M0) * N0)) * P}
t3a12_mod <- function(N0, M0, P, a, a12, kr, h){ ((((a + a12 * M0) * N0) / (kr + N0))  * N0) / (1 + (((a + a12 * M0) * N0) / (kr + N0)) * h * N0) * P}
t3a12h12_mod <- function(N0, M0, P, a, a12, kr, h, h12){ ((((a + a12 * M0) * N0) / (kr + N0))  * N0) / (1 + (((a + a12 * M0) * N0) / (kr + N0)) * (h + h12 * M0) * N0) * P}
t3vv_mod <- function(N0, M0, P, a, kr, h, w){ (((a * N0) / (kr + N0))  * N0) / (1 + ((a * N0) / (kr + N0)) * h * N0 + w * M0) * P}

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

quant = function(col) quantile(col, c(0.025,0.975)) # 95% percentiles

#### Data ####
mort <- 0.13333
conversion <- exp(-15)
nresamp=1000

TIM_data <- loadRData(here("01 - Data", "TIM_data.RData")) # Data

##### Colpidium #####
colp_para_wrk <- TIM_data %>% filter(Exp == "CP")
CP_best_mod <- loadRData(here("05 - Output", "Model Fits", "Colpidium", "CP_t2.a12_fixLV_mod.Rdata"))
CP_mod2 <- loadRData(here("05 - Output", "Model Fits", "Colpidium", "CP_t3.a12_fixLV_mod.Rdata"))
CP_mod3 <- loadRData(here("05 - Output", "Model Fits", "Colpidium", "CP_t2.a12.h12_fixLV_mod.Rdata"))

##### Dexiostoma ##### 
dexio_para_wrk <- TIM_data %>% filter(Exp == "DP")
DP_best_mod <- loadRData(here("05 - Output", "Model Fits", "Dexiostoma", "DP_t3.a12_fixLV_mod.Rdata"))
DP_mod2 <- loadRData(here("05 - Output", "Model Fits", "Dexiostoma", "DP_t3.a12.h12_fixLV_mod.Rdata"))
DP_mod3 <- loadRData(here("05 - Output", "Model Fits", "Dexiostoma", "DP_t3.vv_fixLV_mod.Rdata"))

#### Preditions - Colpidium ####
##### Best Mod ####
Colp_fit1 <- tidy(CP_best_mod)
Colp_fit1_list = list()

for (i in unique(colp_para_wrk$Para_T0.dens.ml)){
  CaP_FR_temp <- sapply(seq(1,1000, length.out=100), function(x) t2a12_mod(N0 = x,
                                                                           M0 = i,
                                                                           
                                                                           a = exp(Colp_fit1[Colp_fit1$term == "a_log",]$estimate),
                                                                           a12 = Colp_fit1[Colp_fit1$term == "a12",]$estimate,
                                                                           h = exp(Colp_fit1[Colp_fit1$term == "h_log",]$estimate), 
                                                                           
                                                                           P = 5)
                        )
  CaP_FR <- append(CaP_FR_temp, i)
  Colp_fit1_list[[i]] <- CaP_FR
}

(Colp_fit1_predict <- as.data.frame(t(do.call(rbind, Colp_fit1_list))))

colnames(Colp_fit1_predict)<-unique(colp_para_wrk$Para_T0.dens.ml); tail(Colp_fit1_predict)
Colp_fit1_predict <- Colp_fit1_predict[1:100,]
Colp_fit1_predict$N0s<-seq(1,1000, length.out=100)
Colp_fit1_predict_long<-Colp_fit1_predict %>% pivot_longer(cols = !N0s, names_to = c("Para_group"), values_to = "Prey_eaten")
Colp_fit1_predict_long$model <-as.factor("M2")

## certainity interval ##
set.seed(1001)

pars.picked = mvrnorm(1000, mu = Colp_fit1$estimate, Sigma = vcov(CP_best_mod)) # pick new parameter values by sampling from multivariate normal distribution based on fit
Colp_fit1_pe_confint <- apply(pars.picked,2,quant)

Colp_fit1_lci_list = list() ; Colp_fit1_uci_list = list() 
for (m in unique(colp_para_wrk$Para_T0.dens.ml)){
  FR_yvals = matrix(0, nrow = nresamp, ncol = nrow(Colp_fit1_predict))
  
  for (i in 1:nresamp) {
    FR_yvals[i,] = sapply(seq(1,1000, length.out=100), function(x) t2a12_mod(N0 = x, 
                                                                             M0 = m,
                                                                             
                                                                             a =  exp(pars.picked[i, 1]),
                                                                             a12 = pars.picked[i, 2],
                                                                             h = exp(pars.picked[i, 3]),
                                                                             
                                                                             P = 5
    ))
  }
  
  FR_conflims = apply(FR_yvals,2,quant) # 95% confidence intervals
  
  Colp_fit1_lci_list[[m]] <- FR_conflims[1,]
  Colp_fit1_uci_list[[m]] <- FR_conflims[2,]
}
#lower confidence interval
(Colp_fit1_LCI <- as.data.frame(t(do.call(rbind, Colp_fit1_lci_list))))
colnames(Colp_fit1_LCI)<-unique(colp_para_wrk$Para_T0.dens.ml); Colp_fit1_LCI$N0s<-seq(1,1000, length.out=100)
Colp_fit1_LCI_long<-Colp_fit1_LCI %>%
  pivot_longer(cols = !N0s, names_to = c("Para_group"), values_to = "LCI")
Colp_fit1_LCI_long$LCI[ Colp_fit1_LCI_long$LCI < 0 ] <- 0


#upper confidence interval
(Colp_fit1_UCI <- as.data.frame(t(do.call(rbind, Colp_fit1_uci_list))))
colnames(Colp_fit1_UCI)<-unique(colp_para_wrk$Para_T0.dens.ml); Colp_fit1_UCI$N0s<-seq(1,1000, length.out=100)
Colp_fit1_UCI_long<-Colp_fit1_UCI %>%
  pivot_longer(cols = !N0s, names_to = c("Para_group"), values_to = "UCI")

#combine
Colp_fit1_CIs <- full_join(Colp_fit1_LCI_long,Colp_fit1_UCI_long)
Colp_fit1_CIs$Para_group <- factor(Colp_fit1_CIs$Para_group, levels=c("1", "25", "40", "60", "90","140", "210", "325", "500"))


##### Fit # 2 ####
Colp_fit2 <- tidy(CP_mod2)
Colp_fit2_list = list()

for (i in unique(colp_para_wrk$Para_T0.dens.ml)){
  CaP_FR_temp <- sapply(seq(1,1000, length.out=100), function(x) t3a12_mod(N0 = x,
                                                                           M0 = i,
                                                                           
                                                                           a = exp(Colp_fit2[Colp_fit2$term == "a_log",]$estimate),
                                                                           a12 = Colp_fit2[Colp_fit2$term == "a12",]$estimate,
                                                                           kr = exp(Colp_fit2[Colp_fit2$term == "kr_log",]$estimate),
                                                                           h = exp(Colp_fit2[Colp_fit2$term == "h_log",]$estimate), 
                                                                           
                                                                           P = 5)
  )
  CaP_FR <- append(CaP_FR_temp, i)
  Colp_fit2_list[[i]] <- CaP_FR
}

(Colp_fit2_predict <- as.data.frame(t(do.call(rbind, Colp_fit2_list))))

colnames(Colp_fit2_predict)<-unique(colp_para_wrk$Para_T0.dens.ml); tail(Colp_fit2_predict)
Colp_fit2_predict <- Colp_fit2_predict[1:100,]
Colp_fit2_predict$N0s<-seq(1,1000, length.out=100)
Colp_fit2_predict_long<-Colp_fit2_predict %>% pivot_longer(cols = !N0s, names_to = c("Para_group"), values_to = "Prey_eaten")
Colp_fit2_predict_long$model <-as.factor("M9")

##### Fit # 3 ####
Colp_fit3 <- tidy(CP_mod3)
Colp_fit3_list = list()

for (i in unique(colp_para_wrk$Para_T0.dens.ml)){
  CaP_FR_temp <- sapply(seq(1,1000, length.out=100), function(x) t2a12h12_mod(N0 = x,
                                                                              M0 = i,
                                                                              
                                                                              a = exp(Colp_fit3[Colp_fit3$term == "a_log",]$estimate),
                                                                              a12 = Colp_fit3[Colp_fit3$term == "a12",]$estimate,
                                                                              h = exp(Colp_fit3[Colp_fit3$term == "h_log",]$estimate),
                                                                              h12 = Colp_fit3[Colp_fit3$term == "h12",]$estimate,
                                                                              
                                                                              P = 5)
  )
  CaP_FR <- append(CaP_FR_temp, i)
  Colp_fit3_list[[i]] <- CaP_FR
}

(Colp_fit3_predict <- as.data.frame(t(do.call(rbind, Colp_fit3_list))))

colnames(Colp_fit3_predict)<-unique(colp_para_wrk$Para_T0.dens.ml); tail(Colp_fit3_predict)
Colp_fit3_predict <- Colp_fit3_predict[1:100,]
Colp_fit3_predict$N0s<-seq(1,1000, length.out=100)
Colp_fit3_predict_long<-Colp_fit3_predict %>% pivot_longer(cols = !N0s, names_to = c("Para_group"), values_to = "Prey_eaten")
Colp_fit3_predict_long$model <-as.factor("M4")


##### Plotting #####
colp_prediction <- rbind(Colp_fit1_predict_long, Colp_fit2_predict_long, Colp_fit3_predict_long)
colp_prediction$Para_group <- factor(colp_prediction$Para_group, levels=c("1", "25", "40", "60", "90","140", "210", "325", "500"))

density_labels <- c("1" = "0 ind/mL", "25" = "25 ind/mL", "40" = "40 ind/mL", "90" = "90 ind/mL", "140" = "140 ind/mL", "210" = "210 ind/mL",  "500" = "500 ind/mL")

## Plot
Colp_alt_fig <- colp_prediction %>% filter(Para_group != "60" & Para_group != "90" & Para_group != "325") %>% 
  ggplot() + 
  xlab("Colpidium Starting Density (ind/mL)") + ylab(" ") +
  ylab(bquote('Colpidium Consumed (ind ml'^-1*' day '^-1*')')) +
  facet_wrap(Para_group ~ ., labeller=labeller(Para_group = density_labels), nrow = 1) +
  geom_line(aes(x=N0s, y=Prey_eaten, group=model, color=model), lty=1, linewidth = 1) +
  geom_ribbon(data = filter(Colp_fit1_CIs, Para_group != "60" & Para_group != "90" & Para_group != "325"), aes(x = N0s, ymin=LCI, ymax=UCI), fill = "#33A02C", linetype = 1, alpha=0.2) +
  scale_color_manual(values=c("#33A02C", "#a1dab4", "#2c7fb8")) +
  coord_cartesian(ylim=c(0, NA)) +
  theme_bw() +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=18), strip.text = element_text(size=16), plot.title = element_text(size=15),
          axis.text.x = element_text(angle = 40, vjust = 1, hjust=1),
          panel.spacing = unit(5, "mm"), panel.border = element_rect(colour = "black", fill=NA, linewidth=2), strip.background = element_rect(colour="white", fill="white"),
          plot.margin = margin(3, 3, 5, 10, "mm"), 
        legend.position="bottom", legend.title=element_text(size=15, face = "bold"), legend.text=element_text(size=15))

tiff(filename="TIM_S3_1.tiff", height=4500/4, width=3500, units="px", res=200, compression="lzw")
Colp_alt_fig
dev.off()

#### * ####
#### * ####
#### Predictions - Dexiostoma ####
##### Best Mod ####
Dexio_fit1 <- tidy(DP_best_mod)
Dexio_fit1_list = list()

for (i in unique(dexio_para_wrk$Para_T0.dens.ml)){
  DaP_FR_temp <- sapply(seq(1,3250, length.out=100), function(x) t3a12_mod(N0 = x,
                                                                           M0 = i,
                                                                           
                                                                           a = exp(Dexio_fit1[Dexio_fit1$term == "a_log",]$estimate),
                                                                           a12 = Dexio_fit1[Dexio_fit1$term == "a12",]$estimate,
                                                                           kr = exp(Dexio_fit1[Dexio_fit1$term == "kr_log",]$estimate),
                                                                           h = exp(Dexio_fit1[Dexio_fit1$term == "h_log",]$estimate), 
                                                                           
                                                                           P = 5)
  )
  DaP_FR <- append(DaP_FR_temp, i)
  Dexio_fit1_list[[i]] <- DaP_FR
}

(Dexio_fit1_predict <- as.data.frame(t(do.call(rbind, Dexio_fit1_list))))

colnames(Dexio_fit1_predict)<-Dexio_fit1_predict[101,]; tail(Dexio_fit1_predict)
Dexio_fit1_predict <- Dexio_fit1_predict[1:100,]
Dexio_fit1_predict$N0s<-seq(1,3250, length.out=100)
Dexio_fit1_predict_long<-Dexio_fit1_predict %>% pivot_longer(cols = !N0s, names_to = c("Para_group"), values_to = "Prey_eaten")
Dexio_fit1_predict_long$model <-as.factor("M9")

## certainity interval ##
set.seed(1001)

pars.picked = mvrnorm(1000, mu = Dexio_fit1$estimate, Sigma = vcov(DP_best_mod)) # pick new parameter values by sampling from multivariate normal distribution based on fit
Dexio_fit1_pe_confint <- apply(pars.picked,2,quant)

Dexio_fit1_lci_list = list() ; Dexio_fit1_uci_list = list() 
for (m in unique(dexio_para_wrk$Para_T0.dens.ml)){
  FR_yvals = matrix(0, nrow = nresamp, ncol = nrow(Dexio_fit1_predict))
  
  for (i in 1:nresamp) {
    FR_yvals[i,] = sapply(seq(1,3250, length.out=100), function(x) t2a12_mod(N0 = x, 
                                                                             M0 = m,
                                                                             
                                                                             a =  exp(pars.picked[i, 1]),
                                                                             a12 = pars.picked[i, 2],
                                                                             h = exp(pars.picked[i, 3]),
                                                                             
                                                                             P = 5
    ))
  }
  
  FR_conflims = apply(FR_yvals,2,quant) # 95% confidence intervals
  
  Dexio_fit1_lci_list[[m]] <- FR_conflims[1,]
  Dexio_fit1_uci_list[[m]] <- FR_conflims[2,]
}
#lower confidence interval
(Dexio_fit1_LCI <- as.data.frame(t(do.call(rbind, Dexio_fit1_lci_list))))
colnames(Dexio_fit1_LCI)<-colnames(Dexio_fit1_predict)[1:9]; Dexio_fit1_LCI$N0s<-seq(1,3250, length.out=100)
Dexio_fit1_LCI_long<-Dexio_fit1_LCI %>%
  pivot_longer(cols = !N0s, names_to = c("Para_group"), values_to = "LCI")
Dexio_fit1_LCI_long$LCI[ Dexio_fit1_LCI_long$LCI < 0 ] <- 0

#upper confidence interval
(Dexio_fit1_UCI <- as.data.frame(t(do.call(rbind, Dexio_fit1_uci_list))))
colnames(Dexio_fit1_UCI)<-colnames(Dexio_fit1_predict)[1:9]; Dexio_fit1_UCI$N0s<-seq(1,3250, length.out=100)
Dexio_fit1_UCI_long<-Dexio_fit1_UCI %>%
  pivot_longer(cols = !N0s, names_to = c("Para_group"), values_to = "UCI")

#combine
Dexio_fit1_CIs <- full_join(Dexio_fit1_LCI_long,Dexio_fit1_UCI_long)
Dexio_fit1_CIs$Para_group <- factor(Dexio_fit1_CIs$Para_group, levels=c("1", "25", "40", "60", "90","140", "210", "325", "500"))

##### Fit # 2 ####
Dexio_fit2 <- tidy(DP_mod2)
Dexio_fit2_list = list()

for (i in unique(dexio_para_wrk$Para_T0.dens.ml)){
  DaP_FR_temp <- sapply(seq(1,3250, length.out=100), function(x) t3a12h12_mod(N0 = x,
                                                                              M0 = i,
                                                                              
                                                                              a = exp(Dexio_fit2[Dexio_fit2$term == "a_log",]$estimate),
                                                                              a12 = Dexio_fit2[Dexio_fit2$term == "a12",]$estimate,
                                                                              kr = exp(Dexio_fit2[Dexio_fit2$term == "kr_log",]$estimate),
                                                                              h = exp(Dexio_fit2[Dexio_fit2$term == "h_log",]$estimate), 
                                                                              h12 = Dexio_fit2[Dexio_fit2$term == "h12",]$estimate,
                                                                              
                                                                              P = 5)
  )
  DaP_FR <- append(DaP_FR_temp, i)
  Dexio_fit2_list[[i]] <- DaP_FR
}

(Dexio_fit2_predict <- as.data.frame(t(do.call(rbind, Dexio_fit2_list))))

colnames(Dexio_fit2_predict)<-Dexio_fit2_predict[101,]; tail(Dexio_fit2_predict)
Dexio_fit2_predict <- Dexio_fit2_predict[1:100,]
Dexio_fit2_predict$N0s<-seq(1,3250, length.out=100)
Dexio_fit2_predict_long<-Dexio_fit2_predict %>% pivot_longer(cols = !N0s, names_to = c("Para_group"), values_to = "Prey_eaten")
Dexio_fit2_predict_long$model <-as.factor("M11")

##### Fit # 3 ####
Dexio_fit3 <- tidy(DP_mod3)
Dexio_fit3_list = list()

for (i in unique(dexio_para_wrk$Para_T0.dens.ml)){
  DaP_FR_temp <- sapply(seq(1,3250, length.out=100), function(x) t3vv_mod(N0 = x,
                                                                          M0 = i,
                                                                          
                                                                          a = exp(Dexio_fit3[Dexio_fit3$term == "a_log",]$estimate),
                                                                          kr = exp(Dexio_fit3[Dexio_fit3$term == "kr_log",]$estimate),
                                                                          h = exp(Dexio_fit3[Dexio_fit3$term == "h_log",]$estimate), 
                                                                          w = exp(Dexio_fit3[Dexio_fit3$term == "log_w",]$estimate), 
                                                                          
                                                                          P = 5)
  )
  DaP_FR <- append(DaP_FR_temp, i)
  Dexio_fit3_list[[i]] <- DaP_FR
}

(Dexio_fit3_predict <- as.data.frame(t(do.call(rbind, Dexio_fit3_list))))

colnames(Dexio_fit3_predict)<-Dexio_fit3_predict[101,]; tail(Dexio_fit3_predict)
Dexio_fit3_predict <- Dexio_fit3_predict[1:100,]
Dexio_fit3_predict$N0s<-seq(1,3250, length.out=100)
Dexio_fit3_predict_long<-Dexio_fit3_predict %>% pivot_longer(cols = !N0s, names_to = c("Para_group"), values_to = "Prey_eaten")
Dexio_fit3_predict_long$model <-as.factor("M12")

##### Fit # 4 ####
Dexio_fit4 <- tidy(DP_mod4)
Dexio_fit4_list = list()

for (i in unique(dexio_para_wrk$Para_T0.dens.ml)){
  DaP_FR_temp <- sapply(seq(1,3250, length.out=100), function(x) t3cm_mod(N0 = x,
                                                                          M0 = i,
                                                                          
                                                                          a = exp(Dexio_fit4[Dexio_fit4$term == "a_log",]$estimate),
                                                                          kr = exp(Dexio_fit4[Dexio_fit4$term == "kr_log",]$estimate),
                                                                          h = exp(Dexio_fit4[Dexio_fit4$term == "h_log",]$estimate), 
                                                                          w = exp(Dexio_fit4[Dexio_fit4$term == "log_w",]$estimate), 
                                                                          
                                                                          P = 5)
  )
  DaP_FR <- append(DaP_FR_temp, i)
  Dexio_fit4_list[[i]] <- DaP_FR
}

(Dexio_fit4_predict <- as.data.frame(t(do.call(rbind, Dexio_fit4_list))))

colnames(Dexio_fit4_predict)<-Dexio_fit4_predict[101,]; tail(Dexio_fit4_predict)
Dexio_fit4_predict <- Dexio_fit4_predict[1:100,]
Dexio_fit4_predict$N0s<-seq(1,3250, length.out=100)
Dexio_fit4_predict_long<-Dexio_fit4_predict %>% pivot_longer(cols = !N0s, names_to = c("Para_group"), values_to = "Prey_eaten")
Dexio_fit4_predict_long$model <-as.factor("Generalized non-linear-a & non-linear-h")

##### Plotting #####
dexio_prediction <- rbind(Dexio_fit1_predict_long, Dexio_fit2_predict_long, Dexio_fit3_predict_long)
dexio_prediction$Para_group <- factor(dexio_prediction$Para_group, levels=c("1", "25", "40", "60", "90","140", "210", "325", "500"))

density_labels <- c("1" = "0 ind/mL", "25" = "25 ind/mL", "40" = "40 ind/mL", "90" = "90 ind/mL", "140" = "140 ind/mL", "210" = "210 ind/mL",  "500" = "500 ind/mL")

## Plot
Dexio_alt_fig <- dexio_prediction %>% filter(Para_group != "60" & Para_group != "90" & Para_group != "325") %>% 
  ggplot() + 
  xlab("Dexiostoma Starting Density (ind/mL)") + ylab(" ") +
  ylab(bquote('Dexiostoma Consumed (ind ml'^-1*' day '^-1*')')) +
  facet_wrap(Para_group ~ ., labeller=labeller(Para_group = density_labels), nrow = 1) +
  geom_line(aes(x=N0s, y=Prey_eaten, group=model, color=model), lty=1, linewidth = 1) +
  geom_ribbon(data = filter(Dexio_fit1_CIs, Para_group != "60" & Para_group != "90" & Para_group != "325"), aes(x = N0s, ymin=LCI, ymax=UCI), fill = "#1f78b4", linetype = 1, alpha=0.2) +
  scale_color_manual(values=c("#1F78B4", "#a6cee3", "#d7b5d8")) +
  coord_cartesian(ylim=c(0, 400)) +
  theme_bw() +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=18), strip.text = element_text(size=16), plot.title = element_text(size=15),
        axis.text.x = element_text(angle = 40, vjust = 1, hjust=1),
        panel.spacing = unit(5, "mm"), panel.border = element_rect(colour = "black", fill=NA, linewidth=2), strip.background = element_rect(colour="white", fill="white"),
        plot.margin = margin(3, 3, 5, 10, "mm"), 
        legend.position="bottom", legend.title=element_text(size=15, face = "bold"), legend.text=element_text(size=15))

tiff(filename="TIM_S3_2.tiff", height=4500/4, width=3500, units="px", res=200, compression="lzw")
Dexio_alt_fig
dev.off()
