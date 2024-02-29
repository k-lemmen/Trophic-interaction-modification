## Testing Dynamic Models comparing with and without TIM ### 

library(bbmle)
library(broom)
library(emdbook)
library(ggpubr)
library(here)
library(janitor)
library(MASS)
library(tidyverse)
library(readr)
library(readxl)
library(scales)

#### Functions ####
source(here("02 - Functions", "Function-T2-a1-h1-LV.R"))
source(here("02 - Functions", "Function-T3-a1-h1-LV.R"))


loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

find_peaks <- function (x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] < x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}


find_valley <- function (x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  valley <- sapply(which(shape > 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] >= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  valley <- unlist(valley)
  valley
}

quant = function(col) quantile(col, c(0.025,0.975)) # 95% percentiles

##### Data ####
mort <- 0.13333
conversion <- exp(-15)
TIM_data <- loadRData(here("01 - Data", "TIM_data.RData")) # Data
CPcomp_pe_mod <- loadRData(here("05 - Output", "Model Fits", "Colpidium", "colp_para_comp_fit.RData"))
DPcomp_pe_mod <- loadRData(here("05 - Output", "Model Fits", "Dexiostoma","dexio_para_comp_fit.RData"))

dexio_para_wrk <- TIM_data %>% filter(Exp == "DP")
colp_para_wrk <- TIM_data %>% filter(Exp == "CP")

colp_TIM_mod <- loadRData(here("05 - Output", "Model Fits", "Colpidium", "CP_t2.a12_fixLV_mod.Rdata"))
colp_T2_mod <- loadRData(here("05 - Output", "Model Fits", "Colpidium", "CP_t2_fixLV_mod.Rdata"))

dexio_TIM_mod  <- loadRData(here("05 - Output", "Model Fits", "Dexiostoma","DP_t3.a12_fixLV_mod.Rdata"))
dexio_T3_mod  <- loadRData(here("05 - Output", "Model Fits", "Dexiostoma","DP_t3_fixLV_mod.Rdata"))

#### Colpidium ####
##### Dynamic Mod Start Point ####
# Start each species at 1/2 K
pred <- 5
colp_mod <- 550/2
colp_prey <- 750/2
Days <- 1100
pred_c <- exp(-5)

###### TIM Model - Type II + linear a TIM ####
(colp_TIM_pe <- tidy(colp_TIM_mod))
colp_TIM_vcov <-vcov(colp_TIM_mod)

colp_TIM_dyn <- eq.odeint.t2.a1.h1.lv(start.v = c(pred, colp_prey, colp_mod),
                                      a = exp(colp_TIM_pe[colp_TIM_pe$term == "a_log",]$estimate),
                                      a12 = colp_TIM_pe[colp_TIM_pe$term == "a12",]$estimate,
                                      h = exp(colp_TIM_pe[colp_TIM_pe$term == "h_log",]$estimate),
                                      h12 = 0, 
                                      
                                      r1 = exp(coef(CPcomp_pe_mod)[1]), 
                                      alpha11 = coef(CPcomp_pe_mod)[2],
                                      alpha12 = coef(CPcomp_pe_mod)[3],
                                      
                                      r2 = exp(coef(CPcomp_pe_mod)[4]), 
                                      alpha22 = coef(CPcomp_pe_mod)[5],
                                      alpha21 = coef(CPcomp_pe_mod)[6],
                                      
                                      c = pred_c, 
                                      m = mort,
                                      Tt = Days,
                                      timesteplength = 0.1
)

colnames(colp_TIM_dyn) <- c('Day', "Predator", "Prey", "Modifier")

###### T2 Model - NO TIM ####
(colp_T2_pe <- tidy(colp_T2_mod))
colp_T2_vcov <-vcov(colp_T2_mod)

### Dynamic Mod Start Point ###
colp_T2_dyn <- eq.odeint.t2.a1.h1.lv(start.v = c(pred, colp_prey, colp_mod),
                                     a = exp(colp_T2_pe[colp_T2_pe$term == "a_log",]$estimate),
                                     a12 = 0,
                                     h = exp(colp_T2_pe[colp_T2_pe$term == "h_log",]$estimate),
                                     h12 = 0, 
                                     
                                     r1 = exp(coef(CPcomp_pe_mod)[1]), 
                                     alpha11 = coef(CPcomp_pe_mod)[2],
                                     alpha12 = coef(CPcomp_pe_mod)[3],
                                     
                                     r2 = exp(coef(CPcomp_pe_mod)[4]), 
                                     alpha22 = coef(CPcomp_pe_mod)[5],
                                     alpha21 = coef(CPcomp_pe_mod)[6],
                                     
                                     c = pred_c, 
                                     m = mort,
                                     Tt = Days,
                                     timesteplength = 0.1)
colnames(colp_T2_dyn) <- c('Day', "Predator", "Prey", "Modifier")

##### Generating 95% CIs #####
set.seed(1001); nresamp=100

###### Colp TIM 95% CI ####
colp_TIM_pars.picked = mvrnorm(1000, mu = colp_TIM_pe$estimate, Sigma = vcov(colp_TIM_mod)) # pick new parameter values by sampling from 

colp_TIM_prey_yvals = matrix(0, nrow = nresamp, ncol = (Days/0.1))
colp_TIM_mod_yvals = matrix(0, nrow = nresamp, ncol = (Days/0.1))
colp_TIM_pred_yvals = matrix(0, nrow = nresamp, ncol = (Days/0.1))

for (i in 1:nresamp) {
  mod_output <- eq.odeint.t2.a1.h1.lv(start.v = c(pred, colp_prey, colp_mod),
                                      a = exp(colp_TIM_pars.picked[i, 1]),
                                      a12 = colp_TIM_pars.picked[i, 2],
                                      h = exp(colp_TIM_pars.picked[i, 3]),
                                      h12 = 0, 
                                      
                                      r1 = exp(coef(CPcomp_pe_mod)[1]), alpha11 = coef(CPcomp_pe_mod)[2], alpha12 = coef(CPcomp_pe_mod)[3],
                                      r2 = exp(coef(CPcomp_pe_mod)[4]), alpha22 = coef(CPcomp_pe_mod)[5], alpha21 = coef(CPcomp_pe_mod)[6],
                                      c = pred_c, 
                                      m = mort,
                                      Tt = Days,
                                      timesteplength = 0.1
  )
  colnames(mod_output) <- c('Day', "Predator", "Prey", "Modifier")
  if(nrow(mod_output) > Days/0.1){
    mod_output <- head(mod_output, - 1) }
  colp_TIM_pred_yvals[i,] <- mod_output[[2]]
  colp_TIM_prey_yvals[i,] <- mod_output[[3]]
  colp_TIM_mod_yvals[i,] <- mod_output[[4]]
  
}

colp_TIM_prey_conflims = as.data.frame(t(apply(colp_TIM_prey_yvals,2,quant))); colp_TIM_prey_conflims$Day <- seq(0,Days-0.1, by=0.1); colnames(colp_TIM_prey_conflims) <- c("LCI","UCI","Day")
colp_TIM_mod_conflims =  as.data.frame(t(apply(colp_TIM_mod_yvals,2,quant))); colp_TIM_mod_conflims$Day <- seq(0,Days-0.1, by=0.1); colnames(colp_TIM_mod_conflims) <- c("LCI","UCI","Day")
colp_TIM_pred_conflims =  as.data.frame(t(apply(colp_TIM_pred_yvals,2,quant))); colp_TIM_pred_conflims$Day <- seq(0,Days-0.1, by=0.1); colnames(colp_TIM_pred_conflims) <- c("LCI","UCI","Day")


###### Colp T2 95% CI ####
colp_T2_pars.picked = mvrnorm(1000, mu = colp_T2_pe$estimate, Sigma = vcov(colp_T2_mod)) # pick new parameter values by sampling from 

colp_T2_prey_yvals = matrix(0, nrow = nresamp, ncol = (Days/0.1))
colp_T2_mod_yvals = matrix(0, nrow = nresamp, ncol = (Days/0.1 ))
colp_T2_pred_yvals = matrix(0, nrow = nresamp, ncol = (Days/0.1))

for (i in 1:nresamp) {
  mod_output <- eq.odeint.t2.a1.h1.lv(start.v = c(pred, colp_prey, colp_mod),
                                      a = exp(colp_T2_pars.picked[i, 1]),
                                      a12 = 0,
                                      h = exp(colp_T2_pars.picked[i, 2]),
                                      h12 = 0, 
                                      
                                      r1 = exp(coef(CPcomp_pe_mod)[1]), alpha11 = coef(CPcomp_pe_mod)[2], alpha12 = coef(CPcomp_pe_mod)[3],
                                      r2 = exp(coef(CPcomp_pe_mod)[4]), alpha22 = coef(CPcomp_pe_mod)[5], alpha21 = coef(CPcomp_pe_mod)[6],
                                      c = pred_c, 
                                      m = mort,
                                      Tt = Days,
                                      timesteplength = 0.1
  )
  colnames(mod_output) <- c('Day', "Predator", "Prey", "Modifier")
  if(nrow(mod_output) > Days/0.1){
    mod_output <- head(mod_output, - 1) }
  colp_T2_pred_yvals[i,] <- mod_output[[2]]
  colp_T2_prey_yvals[i,] <- mod_output[[3]]
  colp_T2_mod_yvals[i,] <- mod_output[[4]]
  
}

####
colp_T2_prey_conflims = as.data.frame(t(apply(colp_T2_prey_yvals,2,quant))); colp_T2_prey_conflims$Day <- seq(0,Days-0.1, by=0.1); colnames(colp_T2_prey_conflims) <- c("LCI","UCI","Day")
colp_T2_mod_conflims =  as.data.frame(t(apply(colp_T2_mod_yvals,2,quant))); colp_T2_mod_conflims$Day <- seq(0,Days-0.1, by=0.1); colnames(colp_T2_mod_conflims) <- c("LCI","UCI","Day")
colp_T2_pred_conflims =  as.data.frame(t(apply(colp_T2_pred_yvals,2,quant))); colp_T2_pred_conflims$Day <- seq(0,Days-0.1, by=0.1); colnames(colp_T2_pred_conflims) <- c("LCI","UCI","Day")

#####-- Time to Extinction --#####
mean_prey_TIM_extinction <- ceiling(head(subset(colp_TIM_dyn, 0.1 > Prey), 1)$Day); 
LCI_prey_TIM_ex <- ceiling(head(subset(colp_TIM_prey_conflims, 0.1 > LCI), 1)$Day); UCI_prey_TIM_ex <- ceiling(head(subset(colp_TIM_prey_conflims, 0.1 > UCI), 1)$Day)

mean_prey_T2_extinction <- ceiling(head(subset(colp_T2_dyn, 0.1 > Prey), 1)$Day)
LCI_prey_T2_ex <- ceiling(head(subset(colp_T2_prey_conflims, 0.1 > LCI), 1)$Day); UCI_prey_T2_ex <- ceiling(head(subset(colp_T2_prey_conflims, 0.1 > UCI), 1)$Day)

mean_pred_TIM_extinction <- ceiling(head(subset(colp_TIM_dyn, 0.1 > Predator), 1)$Day); 
LCI_pred_TIM_ex <- ceiling(head(subset(colp_TIM_pred_conflims, 0.1 > LCI), 1)$Day); UCI_pred_TIM_ex <- ceiling(head(subset(colp_TIM_pred_conflims, 0.1 > UCI), 1)$Day)

mean_pred_T2_extinction <- ceiling(head(subset(colp_T2_dyn, 0.1 > Predator), 1)$Day)
LCI_pred_T2_ex <- ceiling(head(subset(colp_T2_pred_conflims, 0.1 > LCI), 1)$Day); UCI_pred_T2_ex <- ceiling(head(subset(colp_T2_pred_conflims, 0.1 > UCI), 1)$Day)

Spp <- as.factor(c("Predator", "Predator", "Prey", "Prey"))
FR <-  as.factor(c("Pairwise", "TIM", "Pairwise", "TIM"))
TTE <- as.numeric(c(mean_pred_T2_extinction, mean_pred_TIM_extinction, mean_prey_T2_extinction, mean_prey_TIM_extinction))
TTE_LCI <- as.numeric(c(LCI_pred_T2_ex, LCI_pred_TIM_ex, LCI_prey_T2_ex, LCI_prey_TIM_ex))
TTE_UCI <- as.numeric(c(UCI_pred_T2_ex, UCI_pred_TIM_ex, UCI_prey_T2_ex, UCI_prey_TIM_ex))

TTE_pred_prey <- as.data.frame(cbind(Spp, FR, TTE, TTE_LCI, TTE_UCI))
TTE_pred_prey$Spp <- as.factor(TTE_pred_prey$Spp); TTE_pred_prey$FR <- as.factor(TTE_pred_prey$FR)
TTE_pred_prey$FR <- recode(TTE_pred_prey$FR, "1" = "Pairwise", "2" = "TIM")

TTE_pred_prey

#### Dexiostoma ####

##### Dynamic Mod Start Point ####
# If NO predators two species co-exist --> dexio @ 1170 and paramecium at 151 individuals
#

pred <- 5
dexio_mod <- 151/2
dexio_prey <- 1180/2
Days <- 1100
pred_c <- exp(-5) # conversion rate of predator impacts visibility of TIM

###### TIM Model - Type III + linear a TIM ####
(dexio_TIM_pe <- tidy(dexio_TIM_mod))
dexio_TIM_vcov <-vcov(dexio_TIM_mod)

dexio_TIM_dyn <- eq.odeint.t3.a1.h1.lv(start.v = c(pred, dexio_prey, dexio_mod),
                                       a = exp(dexio_TIM_pe[dexio_TIM_pe$term == "a_log",]$estimate),
                                       a12 = dexio_TIM_pe[dexio_TIM_pe$term == "a12",]$estimate,
                                       h = exp(dexio_TIM_pe[dexio_TIM_pe$term == "h_log",]$estimate),
                                       kr = exp(dexio_TIM_pe[dexio_TIM_pe$term == "kr_log",]$estimate),
                                       h12 = 0, 
                                       
                                       r1 = exp(coef(DPcomp_pe_mod)[1]), alpha11 = coef(DPcomp_pe_mod)[2], alpha12 = coef(DPcomp_pe_mod)[3],
                                       r2 = exp(coef(DPcomp_pe_mod)[4]), alpha22 = coef(DPcomp_pe_mod)[5], alpha21 = coef(DPcomp_pe_mod)[6],
                                       
                                       c = pred_c, 
                                       m = mort,
                                       Tt = Days,
                                       timesteplength = 0.1
)

colnames(dexio_TIM_dyn) <- c('Day', "Predator", "Prey", "Modifier")
tail(dexio_TIM_dyn)

###### T3 Model - NO TIM ####
(dexio_T3_pe <- tidy(dexio_T3_mod))
dexio_T3_vcov <-vcov(dexio_T3_mod)

### Dynamic Mod Start Point ###
dexio_T3_dyn <- eq.odeint.t3.a1.h1.lv(start.v = c(pred, dexio_prey, dexio_mod),
                                      a = exp(dexio_T3_pe[dexio_T3_pe$term == "a_log",]$estimate),
                                      a12 = 0,
                                      h = exp(dexio_T3_pe[dexio_T3_pe$term == "h_log",]$estimate),
                                      kr = exp(dexio_T3_pe[dexio_T3_pe$term == "kr_log",]$estimate),
                                      h12 = 0, 
                                      
                                      r1 = exp(coef(DPcomp_pe_mod)[1]), alpha11 = coef(DPcomp_pe_mod)[2], alpha12 = coef(DPcomp_pe_mod)[3],
                                      r2 = exp(coef(DPcomp_pe_mod)[4]), alpha22 = coef(DPcomp_pe_mod)[5], alpha21 = coef(DPcomp_pe_mod)[6],
                                      
                                      c = pred_c, 
                                      m = mort,
                                      Tt = Days,
                                      timesteplength = 0.1)
colnames(dexio_T3_dyn) <- c('Day', "Predator", "Prey", "Modifier")

##### Generating 95% CIs #####
set.seed(1001)
nresamp=100

###### dexio TIM 95% CI ####
dexio_TIM_pars.picked = mvrnorm(1000, mu = dexio_TIM_pe$estimate, Sigma = vcov(dexio_TIM_mod)) # pick new parameter values by sampling from 

dexio_TIM_prey_yvals = matrix(0, nrow = nresamp, ncol = (Days/0.1)+1)
dexio_TIM_mod_yvals = matrix(0, nrow = nresamp, ncol = (Days/0.1)+1)
dexio_TIM_pred_yvals = matrix(0, nrow = nresamp, ncol = (Days/0.1)+1)

for (i in 1:nresamp) {
  mod_output <- eq.odeint.t3.a1.h1.lv(start.v = c(pred, dexio_prey, dexio_mod),
                                      a = exp(dexio_TIM_pars.picked[i, 1]),
                                      a12 = dexio_TIM_pars.picked[i, 2],
                                      h = exp(dexio_TIM_pars.picked[i, 3]),
                                      kr = exp(dexio_TIM_pars.picked[i, 4]),
                                      h12 = 0, 
                                      
                                      r1 = exp(coef(DPcomp_pe_mod)[1]), alpha11 = coef(DPcomp_pe_mod)[2], alpha12 = coef(DPcomp_pe_mod)[3],
                                      r2 = exp(coef(DPcomp_pe_mod)[4]), alpha22 = coef(DPcomp_pe_mod)[5], alpha21 = coef(DPcomp_pe_mod)[6],
                                      c = pred_c, 
                                      m = mort,
                                      Tt = Days,
                                      timesteplength = 0.1
  )
  colnames(mod_output) <- c('Day', "Predator", "Prey", "Modifier")
  dexio_TIM_pred_yvals[i,] <- mod_output[[2]]
  dexio_TIM_prey_yvals[i,] <- mod_output[[3]]
  dexio_TIM_mod_yvals[i,] <- mod_output[[4]]
  
}
dexio_TIM_prey_conflims = as.data.frame(t(apply(dexio_TIM_prey_yvals,2,quant))); dexio_TIM_prey_conflims$Day <- seq(0,Days, by=0.1); colnames(dexio_TIM_prey_conflims) <- c("LCI","UCI","Day")
dexio_TIM_mod_conflims =  as.data.frame(t(apply(dexio_TIM_mod_yvals,2,quant))); dexio_TIM_mod_conflims$Day <- seq(0,Days, by=0.1); colnames(dexio_TIM_mod_conflims) <- c("LCI","UCI","Day")
dexio_TIM_pred_conflims =  as.data.frame(t(apply(dexio_TIM_pred_yvals,2,quant))); dexio_TIM_pred_conflims$Day <- seq(0,Days, by=0.1); colnames(dexio_TIM_pred_conflims) <- c("LCI","UCI","Day")


###### dexio T3 95% CI ####
dexio_T3_pars.picked = mvrnorm(1000, mu = dexio_T3_pe$estimate, Sigma = vcov(dexio_T3_mod)) # pick new parameter values by sampling from 

dexio_T3_prey_yvals = matrix(0, nrow = nresamp, ncol = (Days/0.1)+1)
dexio_T3_mod_yvals = matrix(0, nrow = nresamp, ncol = (Days/0.1)+1)
dexio_T3_pred_yvals = matrix(0, nrow = nresamp, ncol = (Days/0.1)+1)

for (i in 1:nresamp) {
  mod_output <- eq.odeint.t3.a1.h1.lv(start.v = c(pred, dexio_prey, dexio_mod),
                                      a = exp(dexio_T3_pars.picked[i, 1]),
                                      a12 = 0,
                                      h = exp(dexio_T3_pars.picked[i, 2]),
                                      kr = exp(dexio_T3_pars.picked[i, 3]),
                                      h12 = 0, 
                                      
                                      r1 = exp(coef(DPcomp_pe_mod)[1]), alpha11 = coef(DPcomp_pe_mod)[2], alpha12 = coef(DPcomp_pe_mod)[3],
                                      r2 = exp(coef(DPcomp_pe_mod)[4]), alpha22 = coef(DPcomp_pe_mod)[5], alpha21 = coef(DPcomp_pe_mod)[6],
                                      c = pred_c, 
                                      m = mort,
                                      Tt = Days,
                                      timesteplength = 0.1
  )
  colnames(mod_output) <- c('Day', "Predator", "Prey", "Modifier")
  dexio_T3_pred_yvals[i,] <- mod_output[[2]]
  dexio_T3_prey_yvals[i,] <- mod_output[[3]]
  dexio_T3_mod_yvals[i,] <- mod_output[[4]]
  
}

####
dexio_T3_prey_conflims = as.data.frame(t(apply(dexio_T3_prey_yvals,2,quant))); dexio_T3_prey_conflims$Day <- seq(0,Days, by=0.1); colnames(dexio_T3_prey_conflims) <- c("LCI","UCI","Day")
dexio_T3_mod_conflims =  as.data.frame(t(apply(dexio_T3_mod_yvals,2,quant))); dexio_T3_mod_conflims$Day <- seq(0,Days, by=0.1); colnames(dexio_T3_mod_conflims) <- c("LCI","UCI","Day")
dexio_T3_pred_conflims =  as.data.frame(t(apply(dexio_T3_pred_yvals,2,quant))); dexio_T3_pred_conflims$Day <- seq(0,Days, by=0.1); colnames(dexio_T3_pred_conflims) <- c("LCI","UCI","Day")

##### Cycle info #####

###### TIM Model - Type III + linear a TIM ####
#Estimate
(TIM_peak <- dexio_TIM_dyn[find_peaks(dexio_TIM_dyn$Prey, m = 10),]); (TIM_Tpeak <- TIM_peak[TIM_peak$Prey > (max(TIM_peak$Prey)-100) ,])
(TIM_valley <- dexio_TIM_dyn[find_valley(dexio_TIM_dyn$Prey, m = 10),]); (TIM_Tvalley <- TIM_valley[TIM_valley$Prey < (min(TIM_valley$Prey)+0.1) ,])

TIM_PvV <- as.numeric(nrow(TIM_Tpeak) - nrow(TIM_Tvalley))

TIM_period <- mean(diff(round(TIM_Tpeak$Day)))

if(TIM_PvV == 0) {
  TIM_amp <- mean(TIM_Tpeak$Prey - TIM_Tvalley$Prey)/2
} else if(TIM_PvV < 0) {
  TIM_amp <- mean(TIM_Tpeak$Prey - TIM_Tvalley$Prey[-1])/2
} else {
  TIM_amp <- mean(TIM_Tpeak$Prey[-1] - TIM_Tvalley$Prey)/2
}

#Confidence interval
set.seed(3716)
dexio_TIM_pars.picked = mvrnorm(1000, mu = dexio_TIM_pe$estimate, Sigma = vcov(dexio_TIM_mod)) # pick new parameter values by sampling from 

dexio_TIM_prey_yvals = matrix(0, nrow = nresamp, ncol = (Days/1 + 1))

dexio_TIM_period_mean = list();
dexio_TIM_amp_mean = list();

for (i in 1:nresamp) {
  mod_output <- eq.odeint.t3.a1.h1.lv(start.v = c(pred, dexio_prey, dexio_mod),
                                      a = exp(dexio_TIM_pars.picked[i, 1]),
                                      a12 = dexio_TIM_pars.picked[i, 2],
                                      h = exp(dexio_TIM_pars.picked[i, 3]),
                                      kr = exp(dexio_TIM_pars.picked[i, 4]),
                                      h12 = 0, 
                                      
                                      r1 = exp(coef(DPcomp_pe_mod)[1]), alpha11 = coef(DPcomp_pe_mod)[2], alpha12 = coef(DPcomp_pe_mod)[3],
                                      r2 = exp(coef(DPcomp_pe_mod)[4]), alpha22 = coef(DPcomp_pe_mod)[5], alpha21 = coef(DPcomp_pe_mod)[6],
                                      c = pred_c, 
                                      m = mort,
                                      Tt = Days,
                                      timesteplength = 1
  )
  colnames(mod_output) <- c('Day', "Predator", "Prey", "Modifier")
  dexio_TIM_prey_yvals[i,] <- mod_output[[3]]
  
  (TIM_peak <- mod_output[find_peaks(mod_output$Prey, m = 10),]); (TIM_Tpeak <- TIM_peak[TIM_peak$Prey > (max(TIM_peak$Prey)-100) ,])
  (TIM_valley <- mod_output[find_valley(mod_output$Prey, m = 10),]); (TIM_Tvalley <- TIM_valley[TIM_valley$Prey < (min(TIM_valley$Prey)+0.4) ,])
  
  (TIM_PvV <- as.numeric(nrow(TIM_Tpeak) - nrow(TIM_Tvalley)))
  
  TIM_period <- mean(diff(round(TIM_Tpeak$Day)))
  
  if(TIM_PvV == 0) {
    TIM_amp <- mean(TIM_Tpeak$Prey - TIM_Tvalley$Prey)/2
  } else if(TIM_PvV < 0) {
    TIM_amp <- mean(TIM_Tpeak$Prey - TIM_Tvalley$Prey[-1])/2
  } else {
    TIM_amp <- mean(head(TIM_Tpeak$Prey,-1) - TIM_Tvalley$Prey)/2
  }
  
  dexio_TIM_period_mean[[i]] <- TIM_period
  dexio_TIM_amp_mean[[i]] <- TIM_amp
  
}
TIM_period_999<-as.numeric(unlist(dexio_TIM_period_mean))
TIM_period_qurt<-quant(TIM_period_999)

TIM_amp_999<-as.numeric(unlist(dexio_TIM_amp_mean))
TIM_amp_qurt<-quant(TIM_amp_999)

##### Pairwise #####
set.seed(7619)
dexio_T3_pars.picked = mvrnorm(1000, mu = dexio_T3_pe$estimate, Sigma = vcov(dexio_T3_mod))

dexio_T3_prey_yvals = matrix(0, nrow = nresamp, ncol = (Days/1 +1))

dexio_PW_period_mean = list();
dexio_PW_amp_mean = list();

for (i in 1:nresamp) {
  mod_output <- eq.odeint.t3.a1.h1.lv(start.v = c(pred, dexio_prey, dexio_mod),
                                      a = exp(dexio_T3_pars.picked[i, 1]),
                                      a12 = 0,
                                      h = exp(dexio_T3_pars.picked[i, 2]),
                                      kr = exp(dexio_T3_pars.picked[i, 3]),
                                      h12 = 0, 
                                      
                                      r1 = exp(coef(DPcomp_pe_mod)[1]), alpha11 = coef(DPcomp_pe_mod)[2], alpha12 = coef(DPcomp_pe_mod)[3],
                                      r2 = exp(coef(DPcomp_pe_mod)[4]), alpha22 = coef(DPcomp_pe_mod)[5], alpha21 = coef(DPcomp_pe_mod)[6],
                                      c = pred_c, 
                                      m = mort,
                                      Tt = Days,
                                      timesteplength = 1
  )
  colnames(mod_output) <- c('Day', "Predator", "Prey", "Modifier")
  dexio_T3_prey_yvals[i,] <- mod_output[[3]]
  
  (PW_peak <- mod_output[find_peaks(mod_output$Prey, m = 10),]); (PW_Tpeak <- PW_peak[PW_peak$Prey > (max(PW_peak$Prey)-100) ,])
  (PW_valley <- mod_output[find_valley(mod_output$Prey, m = 10),]); (PW_Tvalley <- PW_valley[PW_valley$Prey < (min(PW_valley$Prey)+0.4) ,])
  
  (PW_PvV <- as.numeric(nrow(PW_Tpeak) - nrow(PW_Tvalley)))
  
  PW_period <- mean(diff(round(PW_Tpeak$Day)))
  
  if(PW_PvV == 0) {
    PW_amp <- mean(PW_Tpeak$Prey - PW_Tvalley$Prey)/2
  } else if(PW_PvV < 0) {
    PW_amp <- mean(PW_Tpeak$Prey - PW_Tvalley$Prey[-1])/2
  } else {
    PW_amp <- mean(head(PW_Tpeak$Prey,-1) - PW_Tvalley$Prey)/2
  }
  
  dexio_PW_period_mean[[i]] <- PW_period
  dexio_PW_amp_mean[[i]] <- PW_amp
  
}
PW_period_999<-as.numeric(unlist(dexio_PW_period_mean))
(PW_period_qurt<-quant(PW_period_999))

PW_amp_999<-as.numeric(unlist(dexio_PW_amp_mean))
(PW_amp_qurt<-quant(PW_amp_999))


Period_pt <- c(TIM_period, PW_period); Period_LCI <- c(TIM_period_qurt[[1]], PW_period_qurt[[1]]); Period_UCI <- c(TIM_period_qurt[[2]], PW_period_qurt[[2]])
period_df <- as.data.frame(cbind(Period_pt, Period_LCI, Period_UCI))
period_df$Mod_ID <- as.factor(c("TIM","Pairwise"))
period_df

Amp_pt <- c(TIM_amp, PW_amp); Amp_LCI <- c(TIM_amp_qurt[[1]], PW_amp_qurt[[1]]); Amp_UCI <- c(TIM_amp_qurt[[2]], PW_amp_qurt[[2]])
amp_df <- as.data.frame(cbind(Amp_pt, Amp_LCI, Amp_UCI))
amp_df$Mod_ID <- as.factor(c("TIM","Pairwise"))
amp_df
