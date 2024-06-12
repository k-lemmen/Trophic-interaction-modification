# Model comparisons #
library(tidyverse)
library(here)
library(bbmle)
library(ggpubr)

library(doParallel)
library(foreach)
library(bbmle)

library(readr)
library(readxl)
library(broom)
library(janitor)

#### Functions ####
source(here("02 - Functions", "Function-LV.R"))
source(here("02 - Functions", "Function-T2-a1-h1-LV.R"))
source(here("02 - Functions", "Function-T2-h1-vv-LV.R"))
source(here("02 - Functions", "Function-T2-cm-LV.R"))
source(here("02 - Functions", "Function-T3-a1-h1-LV.R"))
source(here("02 - Functions", "Function-T3-h1-vv-LV.R"))
source(here("02 - Functions", "Function-T3-cm-LV.R"))


loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

##### Data ####
mort <- 0.13333
conversion <- exp(-15)
TIM_data <- loadRData(here("01 - Data", "TIM_data.RData")) # Data
#DPcomp_pe_mod <- loadRData(here("05 - Output", "Model Fits", "Dexiostoma","dexio_para_comp_fit.RData")) # if not running competition fits in this script load here

dexio_para_wrk <- TIM_data %>% filter(Exp == "DP")


#### Model Fits ####

#### 1. Type II #####
##### Search Grid #####
t2.pre_sg <- expand.grid(a_log = seq(-3, -1, length.out = 10), 
                         h_log = seq(-5, -1, length.out = 10), 
                         r1_log = c(-2, -1, -0.5), 
                         r2_log = c(-3, -1, -0.5), 
                         alpha11 = c(0.0005, 0.001, 0.005), 
                         alpha22 = c(0.0005, 0.001, 0.005),  
                         alpha12 = c(0.0005, 0.001, 0.005), 
                         alpha21 = c(0.0005, 0.001, 0.005),
                         c_log = c(-15, -10, -5),
                         sigma = c(0.2),
                         sigma2 = c(0.2)
)

t2.pre_sg$strtpt <- seq(1:nrow(t2.pre_sg))

set.seed(39514)
t2.sim_grid <- t2.pre_sg[sample(t2.pre_sg$strtpt, 19, replace = F), ]

t2.sim_grid$error <- NA

t2.sim_grid$a_log.est <- NA
t2.sim_grid$h_log.est <- NA

t2.sim_grid$r1_log.est <- NA
t2.sim_grid$r2_log.est <- NA

t2.sim_grid$alpha11.est <- NA
t2.sim_grid$alpha12.est <- NA
t2.sim_grid$alpha22.est <- NA
t2.sim_grid$alpha21.est <- NA

t2.sim_grid$c_log.est <- NA
t2.sim_grid$sigma.est <- NA
t2.sim_grid$sigma2.est <- NA

t2.sim_grid$converge <- NA
t2.sim_grid$BIC <- NA
t2.sim_grid$logLik <- NA


# Set the computer to run in Parallel #
bigstatsr::nb_cores()

my.cluster <- parallel::makeCluster(
  bigstatsr::nb_cores(),
  type = "FORK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

t_start <- Sys.time()
# Run Grid #
DP_t2_flexLV_pe <- NA # name of file where parallel runs will be save
mcoptions <- list(preschedule=FALSE)

DP_t2_flexLV_pe <- foreach(i = 1:nrow(t2.sim_grid),
                           .combine = 'rbind',
                           .options.multicore = mcoptions,
                           .packages = c('bbmle', 'odeintr', 'tidyverse', 'R.utils'),
                           .errorhandling = "pass") %dopar%
  {
    compile_sys("FRS_t2_a1_h1_lv", FRS.t2.a1.h1.lv, pars = c("a", "a12", "h", "h12", "r1", "r2","alpha11", "alpha22", "alpha12", "alpha21", "c", "m"), method = "rk54")
    
    r<-tryCatch({ 
      setTimeLimit(35*60)
      fit = mle2(minuslogl = nll.odeint.t2.a1.h1.lv,
                 control = list(maxit = 5000), 
                 method = "Nelder-Mead",
                 
                 start = list(a_log = t2.sim_grid$a_log[[i]], 
                              h_log = t2.sim_grid$h_log[[i]],
                              
                              r1_log = t2.sim_grid$r1_log[[i]], 
                              alpha11 = t2.sim_grid$alpha11[[i]],
                              alpha12 = t2.sim_grid$alpha12[[i]],
                              r2_log = t2.sim_grid$r2_log[[i]], 
                              alpha22 = t2.sim_grid$alpha22[[i]],
                              alpha21 = t2.sim_grid$alpha21[[i]],
                              
                              c_log = t2.sim_grid$c_log[[i]], 
                              sigma = t2.sim_grid$sigma[[i]],
                              sigma2 = t2.sim_grid$sigma2[[i]]
                 ),
                 
                 data = with(list(N0 = Dexio_T0.dens.ml,
                                  Ndead = Dexio_T0.dens.ml - Dexio_T24.dens.ml,
                                  P = Pred_T0.dens.ml, 
                                  P.end = Pred_T24.dens.ml,
                                  M0 = Para_T0.dens.ml,
                                  M0.end = Para_T24.dens.ml,
                                  Tt = Time.day),
                             data = filter(dexio_para_wrk)),
                 
                 fixed = list(m = mort,
                              a12 = 0,
                              h12 = 0) 
      )
      
      if(fit@details$convergence==0) {
        t2.sim_grid$a_log.est[i] <- fit@coef[[1]]
        t2.sim_grid$h_log.est[i] <- fit@coef[[2]]
        
        t2.sim_grid$r1_log.est[i] <- fit@coef[[3]]
        t2.sim_grid$alpha11.est[i] <- fit@coef[[4]]
        t2.sim_grid$alpha12.est[i] <- fit@coef[[5]]
        t2.sim_grid$r2_log.est[i] <- fit@coef[[6]]
        t2.sim_grid$alpha22.est[i] <- fit@coef[[7]]
        t2.sim_grid$alpha21.est[i] <- fit@coef[[8]]
        
        t2.sim_grid$c_log.est[i] <- fit@coef[[9]]
        t2.sim_grid$sigma.est [i] <- fit@coef[[10]]
        t2.sim_grid$sigma2.est [i] <- fit@coef[[11]]
        
        t2.sim_grid$converge[i] <- fit@details$convergence
        t2.sim_grid$BIC[i] <- BIC(fit)
        t2.sim_grid$logLik[i] <- logLik(fit)[[1]]
      } else {
        t2.sim_grid$a_log.est[i] <- NaN
        t2.sim_grid$h_log.est[i] <- NaN
        
        t2.sim_grid$r1_log.est[i] <- NaN
        t2.sim_grid$alpha11.est[i] <- NaN
        t2.sim_grid$alpha12.est[i] <- NaN
        t2.sim_grid$r2_log.est[i] <- NaN
        t2.sim_grid$alpha22.est[i] <- NaN
        t2.sim_grid$alpha21.est[i] <- NaN
        
        t2.sim_grid$c_log.est[i] <- NaN
        t2.sim_grid$sigma.est [i] <- NaN
        t2.sim_grid$sigma2.est [i] <- NaN
        
        t2.sim_grid$converge[i] <- NaN
        t2.sim_grid$BIC[i] <- NaN
        t2.sim_grid$logLik[i] <- NaN
      }
      
      #Remove fit for next iteration
      rm(fit)
    }, error=function(e) list(message = e$message)) # to print error in console use {cat("ERROR :", paste(conditionMessage(e), "Combo", (i)), "\n")}
    
    if(length(r) > 0){
      t2.sim_grid$error[i] <- r 
    } else {
      t2.sim_grid$error[i] <- NaN
    }  
    
    return(t2.sim_grid[i,])
    
  }

parallel::stopCluster(cl = my.cluster)
stopImplicitCluster()

DP_t2_flexLV_pe %>% arrange(BIC)

save(DP_t2_flexLV_pe, file = "DP_t2_flexLV_pe.RData")

###### Best Model #####
(DP_T2_PECs <- DP_t2_flexLV_pe %>% 
   filter(converge == 0) %>% 
   filter(BIC <= min(DP_t2_flexLV_pe$BIC, na.rm = T)+2)) #check for similarity in parameter estimates within 2 AIC of min

min_T2 <- DP_T2_PECs %>% filter(BIC == min(DP_T2_PECs$BIC, na.rm = T)) %>% select(c(a_log.est:sigma2.est, BIC))

if(nrow(min_T2)>1){
  min_T2 <- DP_T2_PECs %>% filter(strtpt %in% sample(filter(DP_T2_PECs, BIC <= min(DP_T2_PECs$BIC, na.rm = T)+0.1)$strtpt,1))
}

DP_min_T2_flexLV_mod =  mle2(minuslogl = nll.odeint.t2.a1.h1.lv,
                             method= "Nelder-Mead",
                             control = list(maxit = 5000),
                             
                             start = list(a_log = min_T2$a_log.est, 
                                          h_log = min_T2$h_log.est,
                                          
                                          r1_log = min_T2$r1_log.est, 
                                          alpha11 = min_T2$alpha11.est,
                                          alpha12 = min_T2$alpha12.est,
                                          r2_log = min_T2$r2_log.est, 
                                          alpha22 = min_T2$alpha22.est,
                                          alpha21 = min_T2$alpha21.est,
                                          
                                          c_log = min_T2$c_log.est, 
                                          sigma = min_T2$sigma.est,
                                          sigma2 = min_T2$sigma2.est
                             ),
                             
                             data = with(list(N0 = Dexio_T0.dens.ml,
                                              Ndead = Dexio_T0.dens.ml - Dexio_T24.dens.ml,
                                              P = Pred_T0.dens.ml, 
                                              P.end = Pred_T24.dens.ml,
                                              M0 = Para_T0.dens.ml,
                                              M0.end = Para_T24.dens.ml,
                                              Tt = Time.day),
                                         data = filter(dexio_para_wrk)),
                             
                             fixed = list(m = mort,
                                          
                                          a12 = 0,
                                          h12 = 0)
)
save(DP_min_T2_flexLV_mod, file = "DP_T2_flexLV_mod.Rdata")
tidy(DP_min_T2_flexLV_mod)

(nrow(DP_T2_PECs)/nrow(DP_t2_flexLV_pe %>% filter(converge == 0))) * 100 # percent of starting combinations within 2 AIC that 
BIC(DP_min_T2_flexLV_mod)

t_end <- Sys.time()
mod1_time <- t_end-t_start

#### 2. Type II + linear a TIM #####
##### Search Grid #####
t2.a12.pre_sg <- expand.grid(a_log = seq(-3, -1, length.out = 5), 
                             a12 = seq(-0.0005, 0.0001, length.out = 5),
                             h_log = seq(-5, -1, length.out = 5),
                             r1_log = c(-2, -1, -0.5), 
                             r2_log = c(-3, -1, -0.5), 
                             alpha11 = c(0.0005, 0.001, 0.005), 
                             alpha22 = c(0.0005, 0.001, 0.005),  
                             alpha12 = c(0.0005, 0.001, 0.005), 
                             alpha21 = c(0.0005, 0.001, 0.005),
                             c_log = c(-15, -10, -5),
                             sigma = c(0.2), 
                             sigma2 = c(0.2)
)
t2.a12.pre_sg$strtpt <- seq(1:nrow(t2.a12.pre_sg))

set.seed(54791)
t2.a12.sim_grid <- t2.a12.pre_sg[sample(t2.a12.pre_sg$strtpt, 98, replace = F), ]

t2.a12.sim_grid$error <- NA

t2.a12.sim_grid$a_log.est <- NA
t2.a12.sim_grid$a12.est <- NA
t2.a12.sim_grid$h_log.est <- NA

t2.a12.sim_grid$r1_log.est <- NA
t2.a12.sim_grid$r2_log.est <- NA

t2.a12.sim_grid$alpha11.est <- NA
t2.a12.sim_grid$alpha12.est <- NA
t2.a12.sim_grid$alpha22.est <- NA
t2.a12.sim_grid$alpha21.est <- NA

t2.a12.sim_grid$c_log.est <- NA
t2.a12.sim_grid$sigma.est <- NA
t2.a12.sim_grid$sigma2.est <- NA

t2.a12.sim_grid$converge <- NA
t2.a12.sim_grid$BIC <- NA
t2.a12.sim_grid$logLik <- NA


# Set the computer to run in Parallel #
bigstatsr::nb_cores()

my.cluster <- parallel::makeCluster(
  bigstatsr::nb_cores(),
  type = "FORK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

t_start <- Sys.time()
# Run Grid #
DP_t2.a12_flexLV_pe <- NA # name of file where parallel runs will be save
mcoptions <- list(preschedule=FALSE)

DP_t2.a12_flexLV_pe <- foreach(i = 1:nrow(t2.a12.sim_grid),
                               .combine = 'rbind',
                               .options.multicore = mcoptions,
                               .packages = c('bbmle', 'odeintr', 'tidyverse', 'R.utils'),
                               .errorhandling = "pass") %dopar%
  {
    compile_sys("FRS_t2_a1_h1_lv", FRS.t2.a1.h1.lv, pars = c("a", "a12", "h", "h12", "r1", "r2","alpha11", "alpha22", "alpha12", "alpha21", "c", "m"), method = "rk54")
    
    r<-tryCatch({ 
      setTimeLimit(35*60)
      fit = mle2(minuslogl = nll.odeint.t2.a1.h1.lv,
                 control = list(maxit = 5000), 
                 method = "Nelder-Mead",
                 
                 start = list(a_log = t2.a12.sim_grid$a_log[[i]], 
                              a12 = t2.a12.sim_grid$a12[[i]], 
                              h_log = t2.a12.sim_grid$h_log[[i]],
                              
                              r1_log = t2.a12.sim_grid$r1_log[[i]], 
                              alpha11 = t2.a12.sim_grid$alpha11[[i]],
                              alpha12 = t2.a12.sim_grid$alpha12[[i]],
                              r2_log = t2.a12.sim_grid$r2_log[[i]], 
                              alpha22 = t2.a12.sim_grid$alpha22[[i]],
                              alpha21 = t2.a12.sim_grid$alpha21[[i]],
                              
                              c_log = t2.a12.sim_grid$c_log[[i]], 
                              sigma = t2.a12.sim_grid$sigma[[i]],
                              sigma2 = t2.a12.sim_grid$sigma2[[i]]
                 ),
                 
                 data = with(list(N0 = Dexio_T0.dens.ml,
                                  Ndead = Dexio_T0.dens.ml - Dexio_T24.dens.ml,
                                  P = Pred_T0.dens.ml, 
                                  P.end = Pred_T24.dens.ml,
                                  M0 = Para_T0.dens.ml,
                                  M0.end = Para_T24.dens.ml,
                                  Tt = Time.day),
                             data = filter(dexio_para_wrk)),
                 
                 fixed = list(m = mort,
                              h12 = 0) 
      )
      
      if(fit@details$convergence==0) {
        t2.a12.sim_grid$a_log.est[i] <- fit@coef[[1]]
        t2.a12.sim_grid$a12.est[i] <- fit@coef[[2]]
        t2.a12.sim_grid$h_log.est[i] <- fit@coef[[3]]
        
        t2.a12.sim_grid$r1_log.est[i] <- fit@coef[[4]]
        t2.a12.sim_grid$alpha11.est[i] <- fit@coef[[5]]
        t2.a12.sim_grid$alpha12.est[i] <- fit@coef[[6]]
        t2.a12.sim_grid$r2_log.est[i] <- fit@coef[[7]]
        t2.a12.sim_grid$alpha22.est[i] <- fit@coef[[8]]
        t2.a12.sim_grid$alpha21.est[i] <- fit@coef[[9]]
        
        t2.a12.sim_grid$c_log.est[i] <- fit@coef[[10]]
        t2.a12.sim_grid$sigma.est [i] <- fit@coef[[11]]
        t2.a12.sim_grid$sigma2.est [i] <- fit@coef[[12]]
        
        t2.a12.sim_grid$converge[i] <- fit@details$convergence
        t2.a12.sim_grid$BIC[i] <- BIC(fit)
        t2.a12.sim_grid$logLik[i] <- logLik(fit)[[1]]
      } else {
        t2.a12.sim_grid$a_log.est[i] <- NaN
        t2.a12.sim_grid$a12.est[i] <- NaN
        t2.a12.sim_grid$h_log.est[i] <- NaN
        
        t2.a12.sim_grid$r1_log.est[i] <- NaN
        t2.a12.sim_grid$alpha11.est[i] <- NaN
        t2.a12.sim_grid$alpha12.est[i] <- NaN
        t2.a12.sim_grid$r2_log.est[i] <- NaN
        t2.a12.sim_grid$alpha22.est[i] <- NaN
        t2.a12.sim_grid$alpha21.est[i] <- NaN
        
        t2.a12.sim_grid$c_log.est[i] <- NaN
        t2.a12.sim_grid$sigma.est [i] <- NaN
        t2.a12.sim_grid$sigma2.est [i] <- NaN
        
        t2.a12.sim_grid$converge[i] <- NaN
        t2.a12.sim_grid$BIC[i] <- NaN
        t2.a12.sim_grid$logLik[i] <- NaN
      }
      
      #Remove fit for next iteration
      rm(fit)
    }, error=function(e) list(message = e$message)) # to print error in console use {cat("ERROR :", paste(conditionMessage(e), "Combo", (i)), "\n")}
    
    if(length(r) > 0){
      t2.a12.sim_grid$error[i] <- r 
    } else {
      t2.a12.sim_grid$error[i] <- NaN
    }  
    
    return(t2.a12.sim_grid[i,])
    
  }

parallel::stopCluster(cl = my.cluster)
stopImplicitCluster()

DP_t2.a12_flexLV_pe %>% arrange(BIC)
save(DP_t2.a12_flexLV_pe, file = "DP_t2.a12_flexLV_pe.RData")

###### Best Model #####
(DP_t2.a12_PECs <- DP_t2.a12_flexLV_pe %>% 
   filter(converge == 0) %>% 
   filter(BIC <= min(DP_t2.a12_flexLV_pe$BIC, na.rm = T)+2)) #check for similarity in parameter estimates within 2 AIC of min

min_t2.a12 <- DP_t2.a12_PECs %>% filter(BIC == min(DP_t2.a12_PECs$BIC, na.rm = T)) %>% select(c(a_log.est:sigma2.est, BIC))

if(nrow(min_t2.a12)>1){
  min_t2.a12 <- DP_t2.a12_PECs %>% filter(strtpt %in% sample(filter(DP_t2.a12_PECs, BIC <= min(DP_t2.a12_PECs$BIC, na.rm = T)+0.1)$strtpt,1))
}

DP_min_t2.a12_flexLV_mod =  mle2(minuslogl = nll.odeint.t2.a1.h1.lv,
                                 method= "Nelder-Mead",
                                 control = list(maxit = 5000),
                                 
                                 start = list(a_log = min_t2.a12$a_log.est, 
                                              h_log = min_t2.a12$h_log.est,
                                              a12 = min_t2.a12$a12.est,
                                              
                                              r1_log = min_t2.a12$r1_log.est, 
                                              alpha11 = min_t2.a12$alpha11.est,
                                              alpha12 = min_t2.a12$alpha12.est,
                                              r2_log = min_t2.a12$r2_log.est, 
                                              alpha22 = min_t2.a12$alpha22.est,
                                              alpha21 = min_t2.a12$alpha21.est,
                                              
                                              c_log = min_t2.a12$c_log.est, 
                                              sigma = min_t2.a12$sigma.est,
                                              sigma2 = min_t2.a12$sigma2.est
                                 ),
                                 
                                 data = with(list(N0 = Dexio_T0.dens.ml,
                                                  Ndead = Dexio_T0.dens.ml - Dexio_T24.dens.ml,
                                                  P = Pred_T0.dens.ml, 
                                                  P.end = Pred_T24.dens.ml,
                                                  M0 = Para_T0.dens.ml,
                                                  M0.end = Para_T24.dens.ml,
                                                  Tt = Time.day),
                                             data = filter(dexio_para_wrk)),
                                 
                                 fixed = list(m = mort,
                                              
                                              h12 = 0)
)
save(DP_min_t2.a12_flexLV_mod, file = "DP_t2.a12_flexLV_mod.Rdata")
tidy(DP_min_t2.a12_flexLV_mod)

(nrow(DP_t2.a12_PECs)/nrow(DP_t2.a12_flexLV_pe %>% filter(converge == 0))) * 100 # percent of starting combinations within 2 AIC that 
BIC(DP_min_t2.a12_flexLV_mod)

t_end <- Sys.time()
mod2_time <- t_end-t_start

#### 3. Type II + linear h TIM #####
##### Search Grid #####
t2.h12.pre_sg <- expand.grid(a_log = seq(-3, -1, length.out = 5), 
                             h12 = seq(-0.0005, 0.0001, length.out = 5),
                             h_log = seq(-5, -1, length.out = 5),
                             r1_log = c(-2, -1, -0.5), 
                             r2_log = c(-3, -1, -0.5), 
                             alpha11 = c(0.0005, 0.001, 0.005), 
                             alpha22 = c(0.0005, 0.001, 0.005),  
                             alpha12 = c(0.0005, 0.001, 0.005), 
                             alpha21 = c(0.0005, 0.001, 0.005),
                             c_log = c(-15, -10, -5),
                             sigma = c(0.2), 
                             sigma2 = c(0.2)
)
t2.h12.pre_sg$strtpt <- seq(1:nrow(t2.h12.pre_sg))

set.seed(64715)
t2.h12.sim_grid <- t2.h12.pre_sg[sample(t2.h12.pre_sg$strtpt, 98, replace = F), ]

t2.h12.sim_grid$error <- NA

t2.h12.sim_grid$a_log.est <- NA
t2.h12.sim_grid$h_log.est <- NA
t2.h12.sim_grid$h12.est <- NA

t2.h12.sim_grid$r1_log.est <- NA
t2.h12.sim_grid$r2_log.est <- NA

t2.h12.sim_grid$alpha11.est <- NA
t2.h12.sim_grid$alpha12.est <- NA
t2.h12.sim_grid$alpha22.est <- NA
t2.h12.sim_grid$alpha21.est <- NA

t2.h12.sim_grid$c_log.est <- NA
t2.h12.sim_grid$sigma.est <- NA
t2.h12.sim_grid$sigma2.est <- NA

t2.h12.sim_grid$converge <- NA
t2.h12.sim_grid$BIC <- NA
t2.h12.sim_grid$logLik <- NA


# Set the computer to run in Parallel #
bigstatsr::nb_cores()

my.cluster <- parallel::makeCluster(
  bigstatsr::nb_cores(),
  type = "FORK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

t_start <- Sys.time()
# Run Grid #
DP_t2.h12_flexLV_pe <- NA # name of file where parallel runs will be save
mcoptions <- list(preschedule=FALSE)

DP_t2.h12_flexLV_pe <- foreach(i = 1:nrow(t2.h12.sim_grid),
                               .combine = 'rbind',
                               .options.multicore = mcoptions,
                               .packages = c('bbmle', 'odeintr', 'tidyverse', 'R.utils'),
                               .errorhandling = "pass") %dopar%
  {
    compile_sys("FRS_t2_a1_h1_lv", FRS.t2.a1.h1.lv, pars = c("a", "a12", "h", "h12", "r1", "r2","alpha11", "alpha22", "alpha12", "alpha21", "c", "m"), method = "rk54")
    
    r<-tryCatch({ 
      setTimeLimit(35*60)
      fit = mle2(minuslogl = nll.odeint.t2.a1.h1.lv,
                 control = list(maxit = 5000), 
                 method = "Nelder-Mead",
                 
                 start = list(a_log = t2.h12.sim_grid$a_log[[i]], 
                              h12 = t2.h12.sim_grid$h12[[i]], 
                              h_log = t2.h12.sim_grid$h_log[[i]],
                              
                              r1_log = t2.h12.sim_grid$r1_log[[i]], 
                              alpha11 = t2.h12.sim_grid$alpha11[[i]],
                              alpha12 = t2.h12.sim_grid$alpha12[[i]],
                              r2_log = t2.h12.sim_grid$r2_log[[i]], 
                              alpha22 = t2.h12.sim_grid$alpha22[[i]],
                              alpha21 = t2.h12.sim_grid$alpha21[[i]],
                              
                              c_log = t2.h12.sim_grid$c_log[[i]], 
                              sigma = t2.h12.sim_grid$sigma[[i]],
                              sigma2 = t2.h12.sim_grid$sigma2[[i]]
                 ),
                 
                 data = with(list(N0 = Dexio_T0.dens.ml,
                                  Ndead = Dexio_T0.dens.ml - Dexio_T24.dens.ml,
                                  P = Pred_T0.dens.ml, 
                                  P.end = Pred_T24.dens.ml,
                                  M0 = Para_T0.dens.ml,
                                  M0.end = Para_T24.dens.ml,
                                  Tt = Time.day),
                             data = filter(dexio_para_wrk)),
                 
                 fixed = list(m = mort,
                              a12 = 0) 
      )
      
      if(fit@details$convergence==0) {
        t2.h12.sim_grid$a_log.est[i] <- fit@coef[[1]]
        t2.h12.sim_grid$h_log.est[i] <- fit@coef[[2]]
        t2.h12.sim_grid$h12.est[i] <- fit@coef[[3]]
        
        t2.h12.sim_grid$r1_log.est[i] <- fit@coef[[4]]
        t2.h12.sim_grid$alpha11.est[i] <- fit@coef[[5]]
        t2.h12.sim_grid$alpha12.est[i] <- fit@coef[[6]]
        t2.h12.sim_grid$r2_log.est[i] <- fit@coef[[7]]
        t2.h12.sim_grid$alpha22.est[i] <- fit@coef[[8]]
        t2.h12.sim_grid$alpha21.est[i] <- fit@coef[[9]]
        
        t2.h12.sim_grid$c_log.est[i] <- fit@coef[[10]]
        t2.h12.sim_grid$sigma.est [i] <- fit@coef[[11]]
        t2.h12.sim_grid$sigma2.est [i] <- fit@coef[[12]]
        
        t2.h12.sim_grid$converge[i] <- fit@details$convergence
        t2.h12.sim_grid$BIC[i] <- BIC(fit)
        t2.h12.sim_grid$logLik[i] <- logLik(fit)[[1]]
      } else {
        t2.h12.sim_grid$a_log.est[i] <- NaN
        t2.h12.sim_grid$h_log.est[i] <- NaN
        t2.h12.sim_grid$h12.est[i] <- NaN
        
        t2.h12.sim_grid$r1_log.est[i] <- NaN
        t2.h12.sim_grid$alpha11.est[i] <- NaN
        t2.h12.sim_grid$alpha12.est[i] <- NaN
        t2.h12.sim_grid$r2_log.est[i] <- NaN
        t2.h12.sim_grid$alpha22.est[i] <- NaN
        t2.h12.sim_grid$alpha21.est[i] <- NaN
        
        t2.h12.sim_grid$c_log.est[i] <- NaN
        t2.h12.sim_grid$sigma.est [i] <- NaN
        t2.h12.sim_grid$sigma2.est [i] <- NaN
        
        t2.h12.sim_grid$converge[i] <- NaN
        t2.h12.sim_grid$BIC[i] <- NaN
        t2.h12.sim_grid$logLik[i] <- NaN
      }
      
      #Remove fit for next iteration
      rm(fit)
    }, error=function(e) list(message = e$message)) # to print error in console use {cat("ERROR :", paste(conditionMessage(e), "Combo", (i)), "\n")}
    
    if(length(r) > 0){
      t2.h12.sim_grid$error[i] <- r 
    } else {
      t2.h12.sim_grid$error[i] <- NaN
    }  
    
    return(t2.h12.sim_grid[i,])
    
  }

parallel::stopCluster(cl = my.cluster)
stopImplicitCluster()

DP_t2.h12_flexLV_pe %>% arrange(BIC)

save(DP_t2.h12_flexLV_pe, file = "DP_t2.h12_flexLV_pe.RData")

###### Best Model #####
(DP_t2.h12_PECs <- DP_t2.h12_flexLV_pe %>% 
   filter(converge == 0) %>% 
   filter(BIC <= min(DP_t2.h12_flexLV_pe$BIC, na.rm = T)+2)) #check for similarity in parameter estimates within 2 AIC of min

min_t2.h12 <- DP_t2.h12_PECs %>% filter(BIC == min(DP_t2.h12_PECs$BIC, na.rm = T)) %>% select(c(a_log.est:sigma2.est, BIC))

if(nrow(min_t2.h12)>1){
  min_t2.h12 <- DP_t2.h12_PECs %>% filter(strtpt %in% sample(filter(DP_t2.h12_PECs, BIC <= min(DP_t2.h12_PECs$BIC, na.rm = T)+0.1)$strtpt,1))
}

DP_min_t2.h12_flexLV_mod =  mle2(minuslogl = nll.odeint.t2.a1.h1.lv,
                                 method= "Nelder-Mead",
                                 control = list(maxit = 5000),
                                 
                                 start = list(a_log = min_t2.h12$a_log.est, 
                                              h_log = min_t2.h12$h_log.est,
                                              h12 = min_t2.h12$h12.est,
                                              
                                              r1_log = min_t2.h12$r1_log.est, 
                                              alpha11 = min_t2.h12$alpha11.est,
                                              alpha12 = min_t2.h12$alpha12.est,
                                              r2_log = min_t2.h12$r2_log.est, 
                                              alpha22 = min_t2.h12$alpha22.est,
                                              alpha21 = min_t2.h12$alpha21.est,
                                              
                                              c_log = min_t2.h12$c_log.est, 
                                              sigma = min_t2.h12$sigma.est,
                                              sigma2 = min_t2.h12$sigma2.est
                                 ),
                                 
                                 data = with(list(N0 = Dexio_T0.dens.ml,
                                                  Ndead = Dexio_T0.dens.ml - Dexio_T24.dens.ml,
                                                  P = Pred_T0.dens.ml, 
                                                  P.end = Pred_T24.dens.ml,
                                                  M0 = Para_T0.dens.ml,
                                                  M0.end = Para_T24.dens.ml,
                                                  Tt = Time.day),
                                             data = filter(dexio_para_wrk)),
                                 
                                 fixed = list(m = mort,
                                              a12=0)
)
save(DP_min_t2.h12_flexLV_mod, file = "DP_t2.h12_flexLV_mod.Rdata")
tidy(DP_min_t2.h12_flexLV_mod)

(nrow(DP_t2.h12_PECs)/nrow(DP_t2.h12_flexLV_pe %>% filter(converge == 0))) * 100 # percent of starting combinations within 2 AIC that 
BIC(DP_min_t2.h12_flexLV_mod)

t_end <- Sys.time()
mod3_time <- t_end-t_start

#### 4. Type II + linear a TIM  + linear h TIM #####
##### Search Grid #####
t2.a12.h12.pre_sg <- expand.grid(a_log = seq(-3, -1, length.out = 5), 
                                 a12 = seq(-0.0005, 0.0001, length.out = 5),
                                 h12 = seq(-0.0005, 0.0001, length.out = 5),
                                 h_log = seq(-5, -1, length.out = 5),
                                 r1_log = c(-2, -1, -0.5), 
                                 r2_log = c(-3, -1, -0.5), 
                                 alpha11 = c(0.0005, 0.001, 0.005), 
                                 alpha22 = c(0.0005, 0.001, 0.005),  
                                 alpha12 = c(0.0005, 0.001, 0.005), 
                                 alpha21 = c(0.0005, 0.001, 0.005),
                                 c_log = c(-15, -10, -5),
                                 sigma = c(0.2), 
                                 sigma2 = c(0.2)
)
t2.a12.h12.pre_sg$strtpt <- seq(1:nrow(t2.a12.h12.pre_sg))

set.seed(47815)
t2.a12.h12.sim_grid <- t2.a12.h12.pre_sg[sample(t2.a12.h12.pre_sg$strtpt, 98, replace = F), ]

t2.a12.h12.sim_grid$error <- NA

t2.a12.h12.sim_grid$a_log.est <- NA
t2.a12.h12.sim_grid$a12.est <- NA
t2.a12.h12.sim_grid$h_log.est <- NA
t2.a12.h12.sim_grid$h12.est <- NA

t2.a12.h12.sim_grid$r1_log.est <- NA
t2.a12.h12.sim_grid$r2_log.est <- NA

t2.a12.h12.sim_grid$alpha11.est <- NA
t2.a12.h12.sim_grid$alpha12.est <- NA
t2.a12.h12.sim_grid$alpha22.est <- NA
t2.a12.h12.sim_grid$alpha21.est <- NA

t2.a12.h12.sim_grid$c_log.est <- NA
t2.a12.h12.sim_grid$sigma.est <- NA
t2.a12.h12.sim_grid$sigma2.est <- NA

t2.a12.h12.sim_grid$converge <- NA
t2.a12.h12.sim_grid$BIC <- NA
t2.a12.h12.sim_grid$logLik <- NA


# Set the computer to run in Parallel #
bigstatsr::nb_cores()

my.cluster <- parallel::makeCluster(
  bigstatsr::nb_cores(),
  type = "FORK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

t_start <- Sys.time()
# Run Grid #
DP_t2.a12.h12_flexLV_pe <- NA # name of file where parallel runs will be save
mcoptions <- list(preschedule=FALSE)

DP_t2.a12.h12_flexLV_pe <- foreach(i = 1:nrow(t2.a12.h12.sim_grid),
                                   .combine = 'rbind',
                                   .options.multicore = mcoptions,
                                   .packages = c('bbmle', 'odeintr', 'tidyverse', 'R.utils'),
                                   .errorhandling = "pass") %dopar%
  {
    compile_sys("FRS_t2_a1_h1_lv", FRS.t2.a1.h1.lv, pars = c("a", "a12", "h", "h12", "r1", "r2","alpha11", "alpha22", "alpha12", "alpha21", "c", "m"), method = "rk54")
    
    r<-tryCatch({ 
      setTimeLimit(35*60)
      fit = mle2(minuslogl = nll.odeint.t2.a1.h1.lv,
                 control = list(maxit = 5000), 
                 method = "Nelder-Mead",
                 
                 start = list(a_log = t2.a12.h12.sim_grid$a_log[[i]], 
                              a12 = t2.a12.h12.sim_grid$a12[[i]], 
                              h_log = t2.a12.h12.sim_grid$h_log[[i]],
                              h12 = t2.a12.h12.sim_grid$h12[[i]], 
                              
                              r1_log = t2.a12.h12.sim_grid$r1_log[[i]], 
                              alpha11 = t2.a12.h12.sim_grid$alpha11[[i]],
                              alpha12 = t2.a12.h12.sim_grid$alpha12[[i]],
                              r2_log = t2.a12.h12.sim_grid$r2_log[[i]], 
                              alpha22 = t2.a12.h12.sim_grid$alpha22[[i]],
                              alpha21 = t2.a12.h12.sim_grid$alpha21[[i]],
                              
                              c_log = t2.a12.h12.sim_grid$c_log[[i]], 
                              sigma = t2.a12.h12.sim_grid$sigma[[i]],
                              sigma2 = t2.a12.h12.sim_grid$sigma2[[i]]
                 ),
                 
                 data = with(list(N0 = Dexio_T0.dens.ml,
                                  Ndead = Dexio_T0.dens.ml - Dexio_T24.dens.ml,
                                  P = Pred_T0.dens.ml, 
                                  P.end = Pred_T24.dens.ml,
                                  M0 = Para_T0.dens.ml,
                                  M0.end = Para_T24.dens.ml,
                                  Tt = Time.day),
                             data = filter(dexio_para_wrk)),
                 
                 fixed = list(m = mort) 
      )
      
      if(fit@details$convergence==0) {
        t2.a12.h12.sim_grid$a_log.est[i] <- fit@coef[[1]]
        t2.a12.h12.sim_grid$a12.est[i] <- fit@coef[[2]]
        t2.a12.h12.sim_grid$h_log.est[i] <- fit@coef[[3]]
        t2.a12.h12.sim_grid$h12.est[i] <- fit@coef[[4]]
        
        t2.a12.h12.sim_grid$r1_log.est[i] <- fit@coef[[5]]
        t2.a12.h12.sim_grid$alpha11.est[i] <- fit@coef[[6]]
        t2.a12.h12.sim_grid$alpha12.est[i] <- fit@coef[[7]]
        t2.a12.h12.sim_grid$r2_log.est[i] <- fit@coef[[8]]
        t2.a12.h12.sim_grid$alpha22.est[i] <- fit@coef[[9]]
        t2.a12.h12.sim_grid$alpha21.est[i] <- fit@coef[[10]]
        
        t2.a12.h12.sim_grid$c_log.est[i] <- fit@coef[[11]]
        t2.a12.h12.sim_grid$sigma.est [i] <- fit@coef[[12]]
        t2.a12.h12.sim_grid$sigma2.est [i] <- fit@coef[[13]]
        
        t2.a12.h12.sim_grid$converge[i] <- fit@details$convergence
        t2.a12.h12.sim_grid$BIC[i] <- BIC(fit)
        t2.a12.h12.sim_grid$logLik[i] <- logLik(fit)[[1]]
      } else {
        t2.a12.h12.sim_grid$a_log.est[i] <- NaN
        t2.a12.h12.sim_grid$a12.est[i] <- NaN
        t2.a12.h12.sim_grid$h_log.est[i] <- NaN
        t2.a12.h12.sim_grid$h12.est[i] <- NaN
        
        t2.a12.h12.sim_grid$r1_log.est[i] <- NaN
        t2.a12.h12.sim_grid$alpha11.est[i] <- NaN
        t2.a12.h12.sim_grid$alpha12.est[i] <- NaN
        t2.a12.h12.sim_grid$r2_log.est[i] <- NaN
        t2.a12.h12.sim_grid$alpha22.est[i] <- NaN
        t2.a12.h12.sim_grid$alpha21.est[i] <- NaN
        
        t2.a12.h12.sim_grid$c_log.est[i] <- NaN
        t2.a12.h12.sim_grid$sigma.est [i] <- NaN
        t2.a12.h12.sim_grid$sigma2.est [i] <- NaN
        
        t2.a12.h12.sim_grid$converge[i] <- NaN
        t2.a12.h12.sim_grid$BIC[i] <- NaN
        t2.a12.h12.sim_grid$logLik[i] <- NaN
      }
      
      #Remove fit for next iteration
      rm(fit)
    }, error=function(e) list(message = e$message)) # to print error in console use {cat("ERROR :", paste(conditionMessage(e), "Combo", (i)), "\n")}
    
    if(length(r) > 0){
      t2.a12.h12.sim_grid$error[i] <- r 
    } else {
      t2.a12.h12.sim_grid$error[i] <- NaN
    }  
    
    return(t2.a12.h12.sim_grid[i,])
    
  }

parallel::stopCluster(cl = my.cluster)
stopImplicitCluster()

DP_t2.a12.h12_flexLV_pe %>% arrange(BIC)

save(DP_t2.a12.h12_flexLV_pe, file = "DP_t2.a12.h12_flexLV_pe.RData")

###### Best Model #####
(DP_t2.a12.h12_PECs <- DP_t2.a12.h12_flexLV_pe %>% 
   filter(converge == 0) %>% 
   filter(BIC <= min(DP_t2.a12.h12_flexLV_pe$BIC, na.rm = T)+2)) #check for similarity in parameter estimates within 2 AIC of min

min_t2.a12.h12 <- DP_t2.a12.h12_PECs %>% filter(BIC == min(DP_t2.a12.h12_PECs$BIC, na.rm = T)) %>% select(c(a_log.est:sigma2.est, BIC))

if(nrow(min_t2.a12.h12)>1){
  min_t2.a12.h12 <- DP_t2.a12.h12_PECs %>% filter(strtpt %in% sample(filter(DP_t2.a12.h12_PECs, BIC <= min(DP_t2.a12.h12_PECs$BIC, na.rm = T)+0.1)$strtpt,1))
}

DP_min_t2.a12.h12_flexLV_mod =  mle2(minuslogl = nll.odeint.t2.a1.h1.lv,
                                     method= "Nelder-Mead",
                                     control = list(maxit = 5000),
                                     
                                     start = list(a_log = min_t2.a12.h12$a_log.est, 
                                                  a12 = min_t2.a12.h12$a12.est,
                                                  h_log = min_t2.a12.h12$h_log.est,
                                                  h12 = min_t2.a12.h12$h12.est,
                                                  
                                                  r1_log = min_t2.a12.h12$r1_log.est, 
                                                  alpha11 = min_t2.a12.h12$alpha11.est,
                                                  alpha12 = min_t2.a12.h12$alpha12.est,
                                                  r2_log = min_t2.a12.h12$r2_log.est, 
                                                  alpha22 = min_t2.a12.h12$alpha22.est,
                                                  alpha21 = min_t2.a12.h12$alpha21.est,
                                                  
                                                  c_log = min_t2.a12.h12$c_log.est, 
                                                  sigma = min_t2.a12.h12$sigma.est,
                                                  sigma2 = min_t2.a12.h12$sigma2.est
                                     ),
                                     
                                     data = with(list(N0 = Dexio_T0.dens.ml,
                                                      Ndead = Dexio_T0.dens.ml - Dexio_T24.dens.ml,
                                                      P = Pred_T0.dens.ml, 
                                                      P.end = Pred_T24.dens.ml,
                                                      M0 = Para_T0.dens.ml,
                                                      M0.end = Para_T24.dens.ml,
                                                      Tt = Time.day),
                                                 data = filter(dexio_para_wrk)),
                                     
                                     fixed = list(m = mort)
)
save(DP_min_t2.a12.h12_flexLV_mod, file = "DP_t2.a12.h12_flexLV_mod.Rdata")
tidy(DP_min_t2.a12.h12_flexLV_mod)

(nrow(DP_t2.a12.h12_PECs)/nrow(DP_t2.a12.h12_flexLV_pe %>% filter(converge == 0))) * 100 # percent of starting combinations within 2 AIC that 
BIC(DP_min_t2.a12.h12_flexLV_mod)

t_end <- Sys.time()
mod4_time <- t_end-t_start

#### 5. Type II + van Veen #####
##### Search Grid #####
t2.vv.pre_sg <- expand.grid(a_log = seq(-3, -1, length.out = 5), 
                            h_log = seq(-5, -1, length.out = 5), 
                            log_w = seq(-10, -5, length.out = 5),
                            r1_log = c(-2, -1, -0.5), 
                            r2_log = c(-3, -1, -0.5), 
                            alpha11 = c(0.0005, 0.001, 0.005), 
                            alpha22 = c(0.0005, 0.001, 0.005),  
                            alpha12 = c(0.0005, 0.001, 0.005), 
                            alpha21 = c(0.0005, 0.001, 0.005),
                            c_log = c(-15, -10, -5),
                            sigma = c(0.2), 
                            sigma2 = c(0.2)
)
t2.vv.pre_sg$strtpt <- seq(1:nrow(t2.vv.pre_sg))

set.seed(13875)
t2.vv.sim_grid <- t2.vv.pre_sg[sample(t2.vv.pre_sg$strtpt, 98, replace = F), ]

t2.vv.sim_grid$error <- NA

t2.vv.sim_grid$a_log.est <- NA
t2.vv.sim_grid$h_log.est <- NA
t2.vv.sim_grid$log_w.est <- NA

t2.vv.sim_grid$r1_log.est <- NA
t2.vv.sim_grid$r2_log.est <- NA

t2.vv.sim_grid$alpha11.est <- NA
t2.vv.sim_grid$alpha12.est <- NA
t2.vv.sim_grid$alpha22.est <- NA
t2.vv.sim_grid$alpha21.est <- NA

t2.vv.sim_grid$c_log.est <- NA
t2.vv.sim_grid$sigma.est <- NA
t2.vv.sim_grid$sigma2.est <- NA

t2.vv.sim_grid$converge <- NA
t2.vv.sim_grid$BIC <- NA
t2.vv.sim_grid$logLik <- NA

# Set the computer to run in Parallel #
bigstatsr::nb_cores()

my.cluster <- parallel::makeCluster(
  bigstatsr::nb_cores(),
  type = "FORK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

t_start <- Sys.time()
# Run Grid #
DP_t2.vv_flexLV_pe <- NA # name of file where parallel runs will be save
mcoptions <- list(preschedule=FALSE)

DP_t2.vv_flexLV_pe <- foreach(i = 1:nrow(t2.vv.sim_grid),
                              .combine = 'rbind',
                              .options.multicore = mcoptions,
                              .packages = c('bbmle', 'odeintr', 'tidyverse', 'R.utils'),
                              .errorhandling = "pass") %dopar%
  {
    compile_sys("FRS_t2_h1_VV_lv", FRS.t2.h1.VV.lv, pars = c("a", "h", "h12", "w", "r1", "r2","alpha11", "alpha22", "alpha12", "alpha21", "c", "m"), method = "rk54")
    
    r<-tryCatch({ 
      setTimeLimit(35*60)
      fit = mle2(minuslogl = nll.odeint.t2.h1.VV.lv,
                 control = list(maxit = 5000), 
                 method = "Nelder-Mead",
                 
                 start = list(a_log = t2.vv.sim_grid$a_log[[i]], 
                              h_log = t2.vv.sim_grid$h_log[[i]],
                              log_w = t2.vv.pre_sg$log_w[[i]],
                              
                              r1_log = t2.vv.sim_grid$r1_log[[i]], 
                              alpha11 = t2.vv.sim_grid$alpha11[[i]],
                              alpha12 = t2.vv.sim_grid$alpha12[[i]],
                              r2_log = t2.vv.sim_grid$r2_log[[i]], 
                              alpha22 = t2.vv.sim_grid$alpha22[[i]],
                              alpha21 = t2.vv.sim_grid$alpha21[[i]],
                              
                              c_log = t2.vv.sim_grid$c_log[[i]],
                              sigma = t2.vv.sim_grid$sigma[[i]],
                              sigma2 = t2.vv.sim_grid$sigma2[[i]]
                 ),
                 
                 data = with(list(N0 = Dexio_T0.dens.ml,
                                  Ndead = Dexio_T0.dens.ml - Dexio_T24.dens.ml,
                                  P = Pred_T0.dens.ml, 
                                  P.end = Pred_T24.dens.ml,
                                  M0 = Para_T0.dens.ml,
                                  M0.end = Para_T24.dens.ml,
                                  Tt = Time.day),
                             data = filter(dexio_para_wrk)),
                 
                 fixed = list(m = mort,
                              h12 = 0) 
      )
      
      if(fit@details$convergence==0) {
        t2.vv.sim_grid$a_log.est[i] <- fit@coef[[1]]
        t2.vv.sim_grid$h_log.est[i] <- fit@coef[[2]]
        t2.vv.sim_grid$log_w.est[i] <- fit@coef[[3]]
        
        t2.vv.sim_grid$r1_log.est[i] <- fit@coef[[4]]
        t2.vv.sim_grid$alpha11.est[i] <- fit@coef[[5]]
        t2.vv.sim_grid$alpha12.est[i] <- fit@coef[[6]]
        t2.vv.sim_grid$r2_log.est[i] <- fit@coef[[7]]
        t2.vv.sim_grid$alpha22.est[i] <- fit@coef[[8]]
        t2.vv.sim_grid$alpha21.est[i] <- fit@coef[[9]]
        
        t2.vv.sim_grid$c_log.est[i] <- fit@coef[[10]]
        t2.vv.sim_grid$sigma.est [i] <- fit@coef[[11]]
        t2.vv.sim_grid$sigma2.est [i] <- fit@coef[[12]]
        
        t2.vv.sim_grid$converge[i] <- fit@details$convergence
        t2.vv.sim_grid$BIC[i] <- BIC(fit)
        t2.vv.sim_grid$logLik[i] <- logLik(fit)[[1]]
      } else {
        t2.vv.sim_grid$a_log.est[i] <- NaN
        t2.vv.sim_grid$h_log.est[i] <- NaN
        t2.vv.sim_grid$log_w.est[i] <- NaN
        
        t2.vv.sim_grid$r1_log.est[i] <- NaN
        t2.vv.sim_grid$alpha11.est[i] <- NaN
        t2.vv.sim_grid$alpha12.est[i] <- NaN
        t2.vv.sim_grid$r2_log.est[i] <- NaN
        t2.vv.sim_grid$alpha22.est[i] <- NaN
        t2.vv.sim_grid$alpha21.est[i] <- NaN
        
        t2.vv.sim_grid$c_log.est[i] <- NaN
        t2.vv.sim_grid$sigma.est [i] <- NaN
        t2.vv.sim_grid$sigma2.est [i] <- NaN
        
        t2.vv.sim_grid$converge[i] <- NaN
        t2.vv.sim_grid$BIC[i] <- NaN
        t2.vv.sim_grid$logLik[i] <- NaN
      }
      
      #Remove fit for next iteration
      rm(fit)
    }, error=function(e) list(message = e$message)) # to print error in console use {cat("ERROR :", paste(conditionMessage(e), "Combo", (i)), "\n")}
    
    if(length(r) > 0){
      t2.vv.sim_grid$error[i] <- r 
    } else {
      t2.vv.sim_grid$error[i] <- NaN
    }  
    
    return(t2.vv.sim_grid[i,])
    
  }

parallel::stopCluster(cl = my.cluster)
stopImplicitCluster()

DP_t2.vv_flexLV_pe %>% arrange(BIC)

save(DP_t2.vv_flexLV_pe, file = "DP_t2.vv_flexLV_pe.RData")

###### Best Model #####
(DP_t2.vv_PECs <- DP_t2.vv_flexLV_pe %>% 
   filter(converge == 0) %>% 
   filter(BIC <= min(DP_t2.vv_flexLV_pe$BIC, na.rm = T)+2)) #check for similarity in parameter estimates within 2 AIC of min

set.seed(273) 
min_t2.vv <- DP_t2.vv_PECs %>% filter(BIC == min(DP_t2.vv_PECs$BIC, na.rm = T)) %>% select(c(a_log.est:sigma2.est, BIC))

if(nrow(min_t2.vv)>1){
  min_t2.vv <- DP_t2.vv_PECs %>% filter(strtpt %in% sample(filter(DP_t2.vv_PECs, BIC <= min(DP_t2.vv_PECs$BIC, na.rm = T)+0.1)$strtpt,1))
}

DP_min_t2.vv_flexLV_mod =  mle2(minuslogl = nll.odeint.t2.h1.VV.lv,
                                method= "Nelder-Mead",
                                control = list(maxit = 5000),
                                
                                start = list(a_log = min_t2.vv$a_log.est, 
                                             h_log = min_t2.vv$h_log.est,
                                             log_w = min_t2.vv$log_w.est,
                                             
                                             r1_log = min_t2.vv$r1_log.est, 
                                             alpha11 = min_t2.vv$alpha11.est,
                                             alpha12 = min_t2.vv$alpha12.est,
                                             r2_log = min_t2.vv$r2_log.est, 
                                             alpha22 = min_t2.vv$alpha22.est,
                                             alpha21 = min_t2.vv$alpha21.est,
                                             
                                             c_log = min_t2.vv$c_log.est, 
                                             sigma = min_t2.vv$sigma.est,
                                             sigma2 = min_t2.vv$sigma2.est
                                ),
                                
                                data = with(list(N0 = Dexio_T0.dens.ml,
                                                 Ndead = Dexio_T0.dens.ml - Dexio_T24.dens.ml,
                                                 P = Pred_T0.dens.ml, 
                                                 P.end = Pred_T24.dens.ml,
                                                 M0 = Para_T0.dens.ml,
                                                 M0.end = Para_T24.dens.ml,
                                                 Tt = Time.day),
                                            data = filter(dexio_para_wrk)),
                                
                                fixed = list(m = mort,
                                             h12=0)
)
save(DP_min_t2.vv_flexLV_mod, file = "DP_t2.vv_flexLV_mod.Rdata")
tidy(DP_min_t2.vv_flexLV_mod)

(nrow(DP_t2.vv_PECs)/nrow(DP_t2.vv_flexLV_pe %>% filter(converge == 0))) * 100 # percent of starting combinations within 2 AIC that 
BIC(DP_min_t2.vv_flexLV_mod)

t_end <- Sys.time()
mod5_time <- NA #t_end-t_start

#### 6. Type II + van Veen + linear h TIM #####
##### Search Grid #####
t2.vv.h12.pre_sg <- expand.grid(a_log = seq(-3, -1, length.out = 5), 
                                h_log = seq(-5, -1, length.out = 5), 
                                h12 = seq(-0.0005, 0.0001, length.out = 5),
                                log_w = seq(-10, -5, length.out = 5),
                                r1_log = c(-2, -1, -0.5), 
                                r2_log = c(-3, -1, -0.5), 
                                alpha11 = c(0.0005, 0.001, 0.005), 
                                alpha22 = c(0.0005, 0.001, 0.005),  
                                alpha12 = c(0.0005, 0.001, 0.005), 
                                alpha21 = c(0.0005, 0.001, 0.005),
                                c_log = c(-15, -10, -5),
                                sigma = c(0.2), 
                                sigma2 = c(0.2)
)
t2.vv.h12.pre_sg$strtpt <- seq(1:nrow(t2.vv.h12.pre_sg))

set.seed(97464)
t2.vv.h12.sim_grid <- t2.vv.h12.pre_sg[sample(t2.vv.h12.pre_sg$strtpt, 98, replace = F), ]

t2.vv.h12.sim_grid$error <- NA

t2.vv.h12.sim_grid$a_log.est <- NA
t2.vv.h12.sim_grid$h_log.est <- NA
t2.vv.h12.sim_grid$h12.est <- NA
t2.vv.h12.sim_grid$log_w.est <- NA

t2.vv.h12.sim_grid$r1_log.est <- NA
t2.vv.h12.sim_grid$r2_log.est <- NA

t2.vv.h12.sim_grid$alpha11.est <- NA
t2.vv.h12.sim_grid$alpha12.est <- NA
t2.vv.h12.sim_grid$alpha22.est <- NA
t2.vv.h12.sim_grid$alpha21.est <- NA

t2.vv.h12.sim_grid$c_log.est <- NA
t2.vv.h12.sim_grid$sigma.est <- NA
t2.vv.h12.sim_grid$sigma2.est <- NA

t2.vv.h12.sim_grid$converge <- NA
t2.vv.h12.sim_grid$BIC <- NA
t2.vv.h12.sim_grid$logLik <- NA


# Set the computer to run in Parallel #
bigstatsr::nb_cores()

my.cluster <- parallel::makeCluster(
  bigstatsr::nb_cores(),
  type = "FORK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

t_start <- Sys.time()
# Run Grid #
DP_t2.vv.h12_flexLV_pe <- NA # name of file where parallel runs will be save
mcoptions <- list(preschedule=FALSE)

DP_t2.vv.h12_flexLV_pe <- foreach(i = 1:nrow(t2.vv.h12.sim_grid),
                                  .combine = 'rbind',
                                  .options.multicore = mcoptions,
                                  .packages = c('bbmle', 'odeintr', 'tidyverse', 'R.utils'),
                                  .errorhandling = "pass") %dopar%
  {
    compile_sys("FRS_t2_h1_vv.h12_lv", FRS.t2.h1.VV.lv, pars = c("a", "h", "h12", "w", "r1", "r2","alpha11", "alpha22", "alpha12", "alpha21", "c", "m"), method = "rk54")
    
    r<-tryCatch({ 
      setTimeLimit(35*60)
      fit = mle2(minuslogl = nll.odeint.t2.h1.VV.lv,
                 control = list(maxit = 5000), 
                 method = "Nelder-Mead",
                 
                 start = list(a_log = t2.vv.h12.sim_grid$a_log[[i]], 
                              h_log = t2.vv.h12.sim_grid$h_log[[i]],
                              h12 = t2.vv.h12.sim_grid$h12[[i]],
                              log_w = t2.vv.h12.pre_sg$log_w[[i]],
                              
                              r1_log = t2.vv.h12.sim_grid$r1_log[[i]], 
                              alpha11 = t2.vv.h12.sim_grid$alpha11[[i]],
                              alpha12 = t2.vv.h12.sim_grid$alpha12[[i]],
                              r2_log = t2.vv.h12.sim_grid$r2_log[[i]], 
                              alpha22 = t2.vv.h12.sim_grid$alpha22[[i]],
                              alpha21 = t2.vv.h12.sim_grid$alpha21[[i]],
                              
                              c_log = t2.vv.h12.sim_grid$c_log[[i]], 
                              sigma = t2.vv.h12.sim_grid$sigma[[i]],
                              sigma2 = t2.vv.h12.sim_grid$sigma2[[i]]
                 ),
                 
                 data = with(list(N0 = Dexio_T0.dens.ml,
                                  Ndead = Dexio_T0.dens.ml - Dexio_T24.dens.ml,
                                  P = Pred_T0.dens.ml, 
                                  P.end = Pred_T24.dens.ml,
                                  M0 = Para_T0.dens.ml,
                                  M0.end = Para_T24.dens.ml,
                                  Tt = Time.day),
                             data = filter(dexio_para_wrk)),
                 
                 fixed = list(m = mort) 
      )
      
      if(fit@details$convergence==0) {
        t2.vv.h12.sim_grid$a_log.est[i] <- fit@coef[[1]]
        t2.vv.h12.sim_grid$h_log.est[i] <- fit@coef[[2]]
        t2.vv.h12.sim_grid$h12.est[i] <- fit@coef[[3]]
        t2.vv.h12.sim_grid$log_w.est[i] <- fit@coef[[4]]
        
        t2.vv.h12.sim_grid$r1_log.est[i] <- fit@coef[[5]]
        t2.vv.h12.sim_grid$alpha11.est[i] <- fit@coef[[6]]
        t2.vv.h12.sim_grid$alpha12.est[i] <- fit@coef[[7]]
        t2.vv.h12.sim_grid$r2_log.est[i] <- fit@coef[[8]]
        t2.vv.h12.sim_grid$alpha22.est[i] <- fit@coef[[9]]
        t2.vv.h12.sim_grid$alpha21.est[i] <- fit@coef[[10]]
        
        t2.vv.h12.sim_grid$c_log.est[i] <- fit@coef[[11]]
        t2.vv.h12.sim_grid$sigma.est [i] <- fit@coef[[12]]
        t2.vv.h12.sim_grid$sigma2.est [i] <- fit@coef[[13]]
        
        t2.vv.h12.sim_grid$converge[i] <- fit@details$convergence
        t2.vv.h12.sim_grid$BIC[i] <- BIC(fit)
        t2.vv.h12.sim_grid$logLik[i] <- logLik(fit)[[1]]
      } else {
        t2.vv.h12.sim_grid$a_log.est[i] <- NaN
        t2.vv.h12.sim_grid$h_log.est[i] <- NaN
        t2.vv.h12.sim_grid$h12.est[i] <- NaN
        t2.vv.h12.sim_grid$log_w.est[i] <- NaN
        
        t2.vv.h12.sim_grid$r1_log.est[i] <- NaN
        t2.vv.h12.sim_grid$alpha11.est[i] <- NaN
        t2.vv.h12.sim_grid$alpha12.est[i] <- NaN
        t2.vv.h12.sim_grid$r2_log.est[i] <- NaN
        t2.vv.h12.sim_grid$alpha22.est[i] <- NaN
        t2.vv.h12.sim_grid$alpha21.est[i] <- NaN
        
        t2.vv.h12.sim_grid$c_log.est[i] <- NaN
        t2.vv.h12.sim_grid$sigma.est [i] <- NaN
        t2.vv.h12.sim_grid$sigma2.est [i] <- NaN
        
        t2.vv.h12.sim_grid$converge[i] <- NaN
        t2.vv.h12.sim_grid$BIC[i] <- NaN
        t2.vv.h12.sim_grid$logLik[i] <- NaN
      }
      
      #Remove fit for next iteration
      rm(fit)
    }, error=function(e) list(message = e$message)) # to print error in console use {cat("ERROR :", paste(conditionMessage(e), "Combo", (i)), "\n")}
    
    if(length(r) > 0){
      t2.vv.h12.sim_grid$error[i] <- r 
    } else {
      t2.vv.h12.sim_grid$error[i] <- NaN
    }  
    
    return(t2.vv.h12.sim_grid[i,])
    
  }

parallel::stopCluster(cl = my.cluster)
stopImplicitCluster()

DP_t2.vv.h12_flexLV_pe %>% arrange(BIC)

save(DP_t2.vv.h12_flexLV_pe, file = "DP_t2.vv.h12_flexLV_pe.RData")

###### Best Model #####
(DP_t2.vv.h12_PECs <- DP_t2.vv.h12_flexLV_pe %>% 
   filter(converge == 0) %>% 
   filter(BIC <= min(DP_t2.vv.h12_flexLV_pe$BIC, na.rm = T)+2)) #check for similarity in parameter estimates within 2 AIC of min

min_t2.vv.h12 <- DP_t2.vv.h12_PECs %>% filter(BIC == min(DP_t2.vv.h12_PECs$BIC, na.rm = T)) %>% select(c(a_log.est:sigma2.est, BIC))

if(nrow(min_t2.vv.h12)>1){
  min_t2.vv.h12 <- DP_t2.vv.h12_PECs %>% filter(strtpt %in% sample(filter(DP_t2.vv.h12_PECs, BIC <= min(DP_t2.vv.h12_PECs$BIC, na.rm = T)+0.1)$strtpt,1))
}

DP_min_t2.vv.h12_flexLV_mod =  mle2(minuslogl = nll.odeint.t2.h1.VV.lv,
                                    method= "Nelder-Mead",
                                    control = list(maxit = 5000),
                                    
                                    start = list(a_log = min_t2.vv.h12$a_log.est, 
                                                 h_log = min_t2.vv.h12$h_log.est,
                                                 h12 = min_t2.vv.h12$h12.est,
                                                 log_w = min_t2.vv.h12$log_w.est,
                                                 
                                                 r1_log = min_t2.vv.h12$r1_log.est, 
                                                 alpha11 = min_t2.vv.h12$alpha11.est,
                                                 alpha12 = min_t2.vv.h12$alpha12.est,
                                                 r2_log = min_t2.vv.h12$r2_log.est, 
                                                 alpha22 = min_t2.vv.h12$alpha22.est,
                                                 alpha21 = min_t2.vv.h12$alpha21.est,
                                                 
                                                 c_log = min_t2.vv.h12$c_log.est, 
                                                 sigma = min_t2.vv.h12$sigma.est,
                                                 sigma2 = min_t2.vv.h12$sigma2.est
                                    ),
                                    
                                    data = with(list(N0 = Dexio_T0.dens.ml,
                                                     Ndead = Dexio_T0.dens.ml - Dexio_T24.dens.ml,
                                                     P = Pred_T0.dens.ml, 
                                                     P.end = Pred_T24.dens.ml,
                                                     M0 = Para_T0.dens.ml,
                                                     M0.end = Para_T24.dens.ml,
                                                     Tt = Time.day),
                                                data = filter(dexio_para_wrk)),
                                    
                                    fixed = list(m = mort)
)
save(DP_min_t2.vv.h12_flexLV_mod, file = "DP_t2.vv.h12_flexLV_mod.Rdata")
tidy(DP_min_t2.vv.h12_flexLV_mod)

(nrow(DP_t2.vv.h12_PECs)/nrow(DP_t2.vv.h12_flexLV_pe %>% filter(converge == 0))) * 100 # percent of starting combinations within 2 AIC that 
BIC(DP_min_t2.vv.h12_flexLV_mod)

t_end <- Sys.time()
mod6_time <- t_end-t_start

#### 7. Type II + Crowley-Martin #####
##### Search Grid #####
t2.cm.pre_sg <- expand.grid(a_log = seq(-3, -1, length.out = 5),
                            h_log = seq(-5, -1, length.out = 5),
                            log_w = seq(-10, -5, length.out = 5),
                            r1_log = c(-2, -1, -0.5), 
                            r2_log = c(-3, -1, -0.5), 
                            alpha11 = c(0.0005, 0.001, 0.005), 
                            alpha22 = c(0.0005, 0.001, 0.005),  
                            alpha12 = c(0.0005, 0.001, 0.005), 
                            alpha21 = c(0.0005, 0.001, 0.005),
                            c_log = c(-15, -10, -5),
                            sigma = c(0.2), 
                            sigma2 = c(0.2)
)
t2.cm.pre_sg$strtpt <- seq(1:nrow(t2.cm.pre_sg))

set.seed(10640)
t2.cm.sim_grid <- t2.cm.pre_sg[sample(t2.cm.pre_sg$strtpt, 98, replace = F), ]

t2.cm.sim_grid$error <- NA

t2.cm.sim_grid$a_log.est <- NA
t2.cm.sim_grid$h_log.est <- NA
t2.cm.sim_grid$log_w.est <- NA

t2.cm.sim_grid$r1_log.est <- NA
t2.cm.sim_grid$r2_log.est <- NA

t2.cm.sim_grid$alpha11.est <- NA
t2.cm.sim_grid$alpha12.est <- NA
t2.cm.sim_grid$alpha22.est <- NA
t2.cm.sim_grid$alpha21.est <- NA

t2.cm.sim_grid$c_log.est <- NA
t2.cm.sim_grid$sigma.est <- NA
t2.cm.sim_grid$sigma2.est <- NA

t2.cm.sim_grid$converge <- NA
t2.cm.sim_grid$BIC <- NA
t2.cm.sim_grid$logLik <- NA


# Set the computer to run in Parallel #
bigstatsr::nb_cores()

my.cluster <- parallel::makeCluster(
  bigstatsr::nb_cores(),
  type = "FORK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

t_start <- Sys.time()
# Run Grid #
DP_t2.cm_flexLV_pe <- NA # name of file where parallel runs will be save
mcoptions <- list(preschedule=FALSE)

DP_t2.cm_flexLV_pe <- foreach(i = 1:nrow(t2.cm.sim_grid),
                              .combine = 'rbind',
                              .options.multicore = mcoptions,
                              .packages = c('bbmle', 'odeintr', 'tidyverse', 'R.utils'),
                              .errorhandling = "pass") %dopar%
  {
    compile_sys("FRS_t2_CM_lv", FRS.t2.CM.lv, pars = c("a", "h", "w", "r1", "r2","alpha11", "alpha22", "alpha12", "alpha21", "c", "m"), method = "rk54")
    
    r<-tryCatch({ 
      setTimeLimit(35*60)
      fit = mle2(minuslogl = nll.odeint.t2.CM.lv,
                 control = list(maxit = 5000), 
                 method = "Nelder-Mead",
                 
                 start = list(a_log = t2.cm.sim_grid$a_log[[i]], 
                              h_log = t2.cm.sim_grid$h_log[[i]],
                              log_w = t2.cm.pre_sg$log_w[[i]],
                              
                              r1_log = t2.cm.sim_grid$r1_log[[i]], 
                              alpha11 = t2.cm.sim_grid$alpha11[[i]],
                              alpha12 = t2.cm.sim_grid$alpha12[[i]],
                              r2_log = t2.cm.sim_grid$r2_log[[i]], 
                              alpha22 = t2.cm.sim_grid$alpha22[[i]],
                              alpha21 = t2.cm.sim_grid$alpha21[[i]],
                              
                              c_log = t2.cm.sim_grid$c_log[[i]], 
                              sigma = t2.cm.sim_grid$sigma[[i]],
                              sigma2 = t2.cm.sim_grid$sigma2[[i]]
                 ),
                 
                 data = with(list(N0 = Dexio_T0.dens.ml,
                                  Ndead = Dexio_T0.dens.ml - Dexio_T24.dens.ml,
                                  P = Pred_T0.dens.ml, 
                                  P.end = Pred_T24.dens.ml,
                                  M0 = Para_T0.dens.ml,
                                  M0.end = Para_T24.dens.ml,
                                  Tt = Time.day),
                             data = filter(dexio_para_wrk)),
                 
                 fixed = list(m = mort) 
      )
      
      if(fit@details$convergence==0) {
        t2.cm.sim_grid$a_log.est[i] <- fit@coef[[1]]
        t2.cm.sim_grid$h_log.est[i] <- fit@coef[[2]]
        t2.cm.sim_grid$log_w.est[i] <- fit@coef[[3]]
        
        t2.cm.sim_grid$r1_log.est[i] <- fit@coef[[4]]
        t2.cm.sim_grid$alpha11.est[i] <- fit@coef[[5]]
        t2.cm.sim_grid$alpha12.est[i] <- fit@coef[[6]]
        t2.cm.sim_grid$r2_log.est[i] <- fit@coef[[7]]
        t2.cm.sim_grid$alpha22.est[i] <- fit@coef[[8]]
        t2.cm.sim_grid$alpha21.est[i] <- fit@coef[[9]]
        
        t2.cm.sim_grid$c_log.est[i] <- fit@coef[[10]]
        t2.cm.sim_grid$sigma.est [i] <- fit@coef[[11]]
        t2.cm.sim_grid$sigma2.est [i] <- fit@coef[[12]]
        
        t2.cm.sim_grid$converge[i] <- fit@details$convergence
        t2.cm.sim_grid$BIC[i] <- BIC(fit)
        t2.cm.sim_grid$logLik[i] <- logLik(fit)[[1]]
      } else {
        t2.cm.sim_grid$a_log.est[i] <- NaN
        t2.cm.sim_grid$h_log.est[i] <- NaN
        t2.cm.sim_grid$log_w.est[i] <- NaN
        
        t2.cm.sim_grid$r1_log.est[i] <- NaN
        t2.cm.sim_grid$alpha11.est[i] <- NaN
        t2.cm.sim_grid$alpha12.est[i] <- NaN
        t2.cm.sim_grid$r2_log.est[i] <- NaN
        t2.cm.sim_grid$alpha22.est[i] <- NaN
        t2.cm.sim_grid$alpha21.est[i] <- NaN
        
        t2.cm.sim_grid$c_log.est[i] <- NaN
        t2.cm.sim_grid$sigma.est [i] <- NaN
        t2.cm.sim_grid$sigma2.est [i] <- NaN
        
        t2.cm.sim_grid$converge[i] <- NaN
        t2.cm.sim_grid$BIC[i] <- NaN
        t2.cm.sim_grid$logLik[i] <- NaN
      }
      
      #Remove fit for next iteration
      rm(fit)
    }, error=function(e) list(message = e$message)) # to print error in console use {cat("ERROR :", paste(conditionMessage(e), "Combo", (i)), "\n")}
    
    if(length(r) > 0){
      t2.cm.sim_grid$error[i] <- r 
    } else {
      t2.cm.sim_grid$error[i] <- NaN
    }  
    
    return(t2.cm.sim_grid[i,])
    
  }

parallel::stopCluster(cl = my.cluster)
stopImplicitCluster()

DP_t2.cm_flexLV_pe %>% arrange(BIC)

save(DP_t2.cm_flexLV_pe, file = "DP_t2.cm_flexLV_pe.RData")

###### Best Model #####
(DP_t2.cm_PECs <- DP_t2.cm_flexLV_pe %>% 
   filter(converge == 0) %>% 
   filter(BIC <= min(DP_t2.cm_flexLV_pe$BIC, na.rm = T)+2)) #check for similarity in parameter estimates within 2 AIC of min

min_t2.cm <- DP_t2.cm_PECs %>% filter(BIC == min(DP_t2.cm_PECs$BIC, na.rm = T)) %>% select(c(a_log.est:sigma2.est, BIC))

if(nrow(min_t2.cm)>1){
  min_t2.cm <- DP_t2.cm_PECs %>% filter(strtpt %in% sample(filter(DP_t2.cm_PECs, BIC <= min(DP_t2.cm_PECs$BIC, na.rm = T)+0.1)$strtpt,1))
}

DP_min_t2.cm_flexLV_mod =  mle2(minuslogl = nll.odeint.t2.CM.lv,
                                method= "Nelder-Mead",
                                control = list(maxit = 5000),
                                
                                start = list(a_log = min_t2.cm$a_log.est, 
                                             h_log = min_t2.cm$h_log.est,
                                             log_w = min_t2.cm$log_w.est,
                                             
                                             r1_log = min_t2.cm$r1_log.est, 
                                             alpha11 = min_t2.cm$alpha11.est,
                                             alpha12 = min_t2.cm$alpha12.est,
                                             r2_log = min_t2.cm$r2_log.est, 
                                             alpha22 = min_t2.cm$alpha22.est,
                                             alpha21 = min_t2.cm$alpha21.est,
                                             
                                             c_log = min_t2.cm$c_log.est, 
                                             sigma = min_t2.cm$sigma.est,
                                             sigma2 = min_t2.cm$sigma2.est
                                ),
                                
                                data = with(list(N0 = Dexio_T0.dens.ml,
                                                 Ndead = Dexio_T0.dens.ml - Dexio_T24.dens.ml,
                                                 P = Pred_T0.dens.ml, 
                                                 P.end = Pred_T24.dens.ml,
                                                 M0 = Para_T0.dens.ml,
                                                 M0.end = Para_T24.dens.ml,
                                                 Tt = Time.day),
                                            data = filter(dexio_para_wrk)),
                                
                                fixed = list(m = mort)
)
save(DP_min_t2.cm_flexLV_mod, file = "DP_t2.cm_flexLV_mod.Rdata")
tidy(DP_min_t2.cm_flexLV_mod)

(nrow(DP_t2.cm_PECs)/nrow(DP_t2.cm_flexLV_pe %>% filter(converge == 0))) * 100 # percent of starting combinations within 2 AIC that 
BIC(DP_min_t2.cm_flexLV_mod)

t_end <- Sys.time()
mod7_time <- t_end-t_start

#### 8. Type III #####
##### Search Grid #####
t3.pre_sg <- expand.grid(a_log = seq(-3, -1, length.out = 10), 
                         h_log = seq(-5, -1, length.out = 10),
                         kr_log = seq(-3, -0.5, length.out = 5),
                         r1_log = c(-2, -1, -0.5), 
                         r2_log = c(-3, -1, -0.5), 
                         alpha11 = c(0.0005, 0.001, 0.005), 
                         alpha22 = c(0.0005, 0.001, 0.005),  
                         alpha12 = c(0.0005, 0.001, 0.005), 
                         alpha21 = c(0.0005, 0.001, 0.005),
                         c_log = c(-15, -10, -5),
                         sigma = c(0.2), 
                         sigma2 = c(0.2)
)
t3.pre_sg$strtpt <- seq(1:nrow(t3.pre_sg))

set.seed(52882)
t3.sim_grid <- t3.pre_sg[sample(t3.pre_sg$strtpt, 98, replace = F), ]

t3.sim_grid$error <- NA

t3.sim_grid$a_log.est <- NA
t3.sim_grid$h_log.est <- NA
t3.sim_grid$kr_log.est <- NA

t3.sim_grid$r1_log.est <- NA
t3.sim_grid$r2_log.est <- NA

t3.sim_grid$alpha11.est <- NA
t3.sim_grid$alpha12.est <- NA
t3.sim_grid$alpha22.est <- NA
t3.sim_grid$alpha21.est <- NA

t3.sim_grid$c_log.est <- NA
t3.sim_grid$sigma.est <- NA
t3.sim_grid$sigma2.est <- NA

t3.sim_grid$converge <- NA
t3.sim_grid$BIC <- NA
t3.sim_grid$logLik <- NA



# Set the computer to run in Parallel #
bigstatsr::nb_cores()

my.cluster <- parallel::makeCluster(
  bigstatsr::nb_cores(),
  type = "FORK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

t_start <- Sys.time()
# Run Grid #
DP_t3_flexLV_pe <- NA # name of file where parallel runs will be save
mcoptions <- list(preschedule=FALSE)

DP_t3_flexLV_pe <- foreach(i = 1:nrow(t3.sim_grid),
                           .combine = 'rbind',
                           .options.multicore = mcoptions,
                           .packages = c('bbmle', 'odeintr', 'tidyverse', 'R.utils'),
                           .errorhandling = "pass") %dopar%
  {
    compile_sys("FRS_t3_a1_h1_lv", FRS.t3.a1.h1.lv, pars = c("a", "a12", "kr", "h", "h12", "r1", "r2","alpha11", "alpha22", "alpha12", "alpha21", "c", "m"), method = "rk54")
    
    r<-tryCatch({ 
      setTimeLimit(35*60)
      fit = mle2(minuslogl = nll.odeint.t3.a1.h1.lv,
                 control = list(maxit = 5000), 
                 method = "Nelder-Mead",
                 
                 start = list(a_log = t3.sim_grid$a_log[[i]], 
                              h_log = t3.sim_grid$h_log[[i]],
                              kr_log = t3.sim_grid$kr_log[[i]],
                              
                              r1_log = t3.sim_grid$r1_log[[i]], 
                              alpha11 = t3.sim_grid$alpha11[[i]],
                              alpha12 = t3.sim_grid$alpha12[[i]],
                              r2_log = t3.sim_grid$r2_log[[i]], 
                              alpha22 = t3.sim_grid$alpha22[[i]],
                              alpha21 = t3.sim_grid$alpha21[[i]],
                              
                              c_log = t3.sim_grid$c_log[[i]], 
                              sigma = t3.sim_grid$sigma[[i]],
                              sigma2 = t3.sim_grid$sigma2[[i]]
                 ),
                 
                 data = with(list(N0 = Dexio_T0.dens.ml,
                                  Ndead = Dexio_T0.dens.ml - Dexio_T24.dens.ml,
                                  P = Pred_T0.dens.ml, 
                                  P.end = Pred_T24.dens.ml,
                                  M0 = Para_T0.dens.ml,
                                  M0.end = Para_T24.dens.ml,
                                  Tt = Time.day),
                             data = filter(dexio_para_wrk)),
                 
                 fixed = list(m = mort,
                              a12 = 0,
                              h12 = 0) 
      )
      
      if(fit@details$convergence==0) {
        t3.sim_grid$a_log.est[i] <- fit@coef[[1]]
        t3.sim_grid$h_log.est[i] <- fit@coef[[2]]
        t3.sim_grid$kr_log.est[i] <- fit@coef[[3]]
        
        t3.sim_grid$r1_log.est[i] <- fit@coef[[4]]
        t3.sim_grid$alpha11.est[i] <- fit@coef[[5]]
        t3.sim_grid$alpha12.est[i] <- fit@coef[[6]]
        t3.sim_grid$r2_log.est[i] <- fit@coef[[7]]
        t3.sim_grid$alpha22.est[i] <- fit@coef[[8]]
        t3.sim_grid$alpha21.est[i] <- fit@coef[[9]]
        
        t3.sim_grid$c_log.est[i] <- fit@coef[[10]]
        t3.sim_grid$sigma.est [i] <- fit@coef[[11]]
        t3.sim_grid$sigma2.est [i] <- fit@coef[[12]]
        
        t3.sim_grid$converge[i] <- fit@details$convergence
        t3.sim_grid$BIC[i] <- BIC(fit)
        t3.sim_grid$logLik[i] <- logLik(fit)[[1]]
      } else {
        t3.sim_grid$a_log.est[i] <- NaN
        t3.sim_grid$h_log.est[i] <- NaN
        t3.sim_grid$kr_log.est[i] <- NaN
        
        t3.sim_grid$r1_log.est[i] <- NaN
        t3.sim_grid$alpha11.est[i] <- NaN
        t3.sim_grid$alpha12.est[i] <- NaN
        t3.sim_grid$r2_log.est[i] <- NaN
        t3.sim_grid$alpha22.est[i] <- NaN
        t3.sim_grid$alpha21.est[i] <- NaN
        
        t3.sim_grid$c_log.est[i] <- NaN
        t3.sim_grid$sigma.est [i] <- NaN
        t3.sim_grid$sigma2.est [i] <- NaN
        
        t3.sim_grid$converge[i] <- NaN
        t3.sim_grid$BIC[i] <- NaN
        t3.sim_grid$logLik[i] <- NaN
      }
      
      #Remove fit for next iteration
      rm(fit)
    }, error=function(e) list(message = e$message)) # to print error in console use {cat("ERROR :", paste(conditionMessage(e), "Combo", (i)), "\n")}
    
    if(length(r) > 0){
      t3.sim_grid$error[i] <- r 
    } else {
      t3.sim_grid$error[i] <- NaN
    }  
    
    return(t3.sim_grid[i,])
    
  }

parallel::stopCluster(cl = my.cluster)
stopImplicitCluster()

DP_t3_flexLV_pe %>% arrange(BIC)

save(DP_t3_flexLV_pe, file = "DP_t3_flexLV_pe.RData")

###### Best Model #####
(DP_t3_PECs <- DP_t3_flexLV_pe %>% 
   filter(converge == 0) %>% 
   filter(BIC <= min(DP_t3_flexLV_pe$BIC, na.rm = T)+2)) #check for similarity in parameter estimates within 2 AIC of min

min_T3 <- DP_t3_PECs %>% filter(BIC == min(DP_t3_PECs$BIC, na.rm = T)) %>% select(c(a_log.est:sigma2.est, BIC))

if(nrow(min_T3)>1){
  min_T3 <- DP_t3_PECs %>% filter(strtpt %in% sample(filter(DP_t3_PECs, BIC <= min(DP_t3_PECs$BIC, na.rm = T)+0.1)$strtpt,1))
}

DP_min_t3_flexLV_mod =  mle2(minuslogl = nll.odeint.t3.a1.h1.lv,
                             method= "Nelder-Mead",
                             control = list(maxit = 5000),
                             
                             start = list(a_log = min_T3$a_log.est, 
                                          h_log = min_T3$h_log.est,
                                          kr_log = min_T3$kr_log.est,
                                          
                                          r1_log = min_T2$r1_log.est, 
                                          alpha11 = min_T2$alpha11.est,
                                          alpha12 = min_T2$alpha12.est,
                                          r2_log = min_T2$r2_log.est, 
                                          alpha22 = min_T2$alpha22.est,
                                          alpha21 = min_T2$alpha21.est,
                                          
                                          c_log = min_T2$c_log.est, 
                                          sigma = min_T3$sigma.est,
                                          sigma2 = min_T3$sigma2.est
                             ),
                             
                             data = with(list(N0 = Dexio_T0.dens.ml,
                                              Ndead = Dexio_T0.dens.ml - Dexio_T24.dens.ml,
                                              P = Pred_T0.dens.ml, 
                                              P.end = Pred_T24.dens.ml,
                                              M0 = Para_T0.dens.ml,
                                              M0.end = Para_T24.dens.ml,
                                              Tt = Time.day),
                                         data = filter(dexio_para_wrk)),
                             
                             fixed = list(m = mort,
                                          a12 = 0,
                                          h12 = 0)
)
save(DP_min_t3_flexLV_mod, file = "DP_t3_flexLV_mod.Rdata")
tidy(DP_min_t3_flexLV_mod)

(nrow(DP_t3_PECs)/nrow(DP_t3_flexLV_pe %>% filter(converge == 0))) * 100 # percent of starting combinations within 2 AIC that 
BIC(DP_min_t3_flexLV_mod)

t_end <- Sys.time()
mod8_time <- t_end-t_start

#### 9. Type III + linear a TIM #####
##### Search Grid #####
t3.a12.pre_sg <- expand.grid(a_log = seq(-3, -1, length.out = 5), 
                             a12 = seq(-0.0005, 0.0001, length.out = 5),
                             h_log = seq(-5, -1, length.out = 5),
                             kr_log = seq(-3, -0.5, length.out = 5),
                             r1_log = c(-2, -1, -0.5), 
                             r2_log = c(-3, -1, -0.5), 
                             alpha11 = c(0.0005, 0.001, 0.005), 
                             alpha22 = c(0.0005, 0.001, 0.005),  
                             alpha12 = c(0.0005, 0.001, 0.005), 
                             alpha21 = c(0.0005, 0.001, 0.005),
                             c_log = c(-15, -10, -5),
                             sigma = c(0.2), 
                             sigma2 = c(0.2)
)
t3.a12.pre_sg$strtpt <- seq(1:nrow(t3.a12.pre_sg))

set.seed(22761)
t3.a12.sim_grid <- t3.a12.pre_sg[sample(t3.a12.pre_sg$strtpt, 98, replace = F), ]

t3.a12.sim_grid$error <- NA

t3.a12.sim_grid$a_log.est <- NA
t3.a12.sim_grid$a12.est <- NA
t3.a12.sim_grid$h_log.est <- NA
t3.a12.sim_grid$kr_log.est <- NA

t3.a12.sim_grid$r1_log.est <- NA
t3.a12.sim_grid$r2_log.est <- NA

t3.a12.sim_grid$alpha11.est <- NA
t3.a12.sim_grid$alpha12.est <- NA
t3.a12.sim_grid$alpha22.est <- NA
t3.a12.sim_grid$alpha21.est <- NA

t3.a12.sim_grid$c_log.est <- NA
t3.a12.sim_grid$sigma.est <- NA
t3.a12.sim_grid$sigma2.est <- NA

t3.a12.sim_grid$converge <- NA
t3.a12.sim_grid$BIC <- NA
t3.a12.sim_grid$logLik <- NA


# Set the computer to run in Parallel #
bigstatsr::nb_cores()

my.cluster <- parallel::makeCluster(
  bigstatsr::nb_cores(),
  type = "FORK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

t_start <- Sys.time()
# Run Grid #
DP_t3.a12_flexLV_pe <- NA # name of file where parallel runs will be save
mcoptions <- list(preschedule=FALSE)

DP_t3.a12_flexLV_pe <- foreach(i = 1:nrow(t3.a12.sim_grid),
                               .combine = 'rbind',
                               .options.multicore = mcoptions,
                               .packages = c('bbmle', 'odeintr', 'tidyverse', 'R.utils'),
                               .errorhandling = "pass") %dopar%
  {
    compile_sys("FRS_t3_a1_h1_lv", FRS.t3.a1.h1.lv, pars = c("a", "a12", "kr", "h", "h12", "r1", "r2","alpha11", "alpha22", "alpha12", "alpha21", "c", "m"), method = "rk54")
    
    
    r<-tryCatch({ 
      setTimeLimit(35*60)
      fit = mle2(minuslogl = nll.odeint.t3.a1.h1.lv,
                 control = list(maxit = 5000), 
                 method = "Nelder-Mead",
                 
                 start = list(a_log = t3.a12.sim_grid$a_log[[i]], 
                              a12 = t3.a12.sim_grid$a12[[i]], 
                              h_log = t3.a12.sim_grid$h_log[[i]],
                              kr_log = t3.a12.sim_grid$kr_log[[i]],
                              
                              r1_log = t3.a12.sim_grid$r1_log[[i]], 
                              alpha11 = t3.a12.sim_grid$alpha11[[i]],
                              alpha12 = t3.a12.sim_grid$alpha12[[i]],
                              r2_log = t3.a12.sim_grid$r2_log[[i]], 
                              alpha22 = t3.a12.sim_grid$alpha22[[i]],
                              alpha21 = t3.a12.sim_grid$alpha21[[i]],
                              
                              c_log = t3.a12.sim_grid$c_log[[i]], 
                              sigma = t3.a12.sim_grid$sigma[[i]],
                              sigma2 = t3.a12.sim_grid$sigma2[[i]]
                 ),
                 
                 data = with(list(N0 = Dexio_T0.dens.ml,
                                  Ndead = Dexio_T0.dens.ml - Dexio_T24.dens.ml,
                                  P = Pred_T0.dens.ml, 
                                  P.end = Pred_T24.dens.ml,
                                  M0 = Para_T0.dens.ml,
                                  M0.end = Para_T24.dens.ml,
                                  Tt = Time.day),
                             data = filter(dexio_para_wrk)),
                 
                 fixed = list(m = mort,
                              h12 = 0) 
      )
      
      if(fit@details$convergence==0) {
        t3.a12.sim_grid$a_log.est[i] <- fit@coef[[1]]
        t3.a12.sim_grid$a12.est[i] <- fit@coef[[2]]
        t3.a12.sim_grid$h_log.est[i] <- fit@coef[[3]]
        t3.a12.sim_grid$kr_log.est[i] <- fit@coef[[4]]
        
        t3.a12.sim_grid$r1_log.est[i] <- fit@coef[[5]]
        t3.a12.sim_grid$alpha11.est[i] <- fit@coef[[6]]
        t3.a12.sim_grid$alpha12.est[i] <- fit@coef[[7]]
        t3.a12.sim_grid$r2_log.est[i] <- fit@coef[[8]]
        t3.a12.sim_grid$alpha22.est[i] <- fit@coef[[9]]
        t3.a12.sim_grid$alpha21.est[i] <- fit@coef[[10]]
        
        t3.a12.sim_grid$c_log.est[i] <- fit@coef[[11]]
        t3.a12.sim_grid$sigma.est [i] <- fit@coef[[12]]
        t3.a12.sim_grid$sigma2.est [i] <- fit@coef[[13]]
        
        t3.a12.sim_grid$converge[i] <- fit@details$convergence
        t3.a12.sim_grid$BIC[i] <- BIC(fit)
        t3.a12.sim_grid$logLik[i] <- logLik(fit)[[1]]
      } else {
        t3.a12.sim_grid$a_log.est[i] <- NaN
        t3.a12.sim_grid$a12.est[i] <- NaN
        t3.a12.sim_grid$h_log.est[i] <- NaN
        t3.a12.sim_grid$kr_log.est[i] <- NaN
        
        
        t3.a12.sim_grid$r1_log.est[i] <- NaN
        t3.a12.sim_grid$alpha11.est[i] <- NaN
        t3.a12.sim_grid$alpha12.est[i] <- NaN
        t3.a12.sim_grid$r2_log.est[i] <- NaN
        t3.a12.sim_grid$alpha22.est[i] <- NaN
        t3.a12.sim_grid$alpha21.est[i] <- NaN
        
        t3.a12.sim_grid$c_log.est[i] <- NaN
        t3.a12.sim_grid$sigma.est [i] <- NaN
        t3.a12.sim_grid$sigma2.est [i] <- NaN
        
        t3.a12.sim_grid$converge[i] <- NaN
        t3.a12.sim_grid$BIC[i] <- NaN
        t3.a12.sim_grid$logLik[i] <- NaN
      }
      
      #Remove fit for next iteration
      rm(fit)
    }, error=function(e) list(message = e$message)) # to print error in console use {cat("ERROR :", paste(conditionMessage(e), "Combo", (i)), "\n")}
    
    if(length(r) > 0){
      t3.a12.sim_grid$error[i] <- r 
    } else {
      t3.a12.sim_grid$error[i] <- NaN
    }  
    
    return(t3.a12.sim_grid[i,])
    
  }

parallel::stopCluster(cl = my.cluster)
stopImplicitCluster()

DP_t3.a12_flexLV_pe %>% arrange(BIC)
save(DP_t3.a12_flexLV_pe, file = "DP_t3.a12_flexLV_pe.RData")

###### Best Model #####
(DP_t3.a12_PECs <- DP_t3.a12_flexLV_pe %>% 
   filter(converge == 0) %>% 
   filter(BIC <= min(DP_t3.a12_flexLV_pe$BIC, na.rm = T)+2)) #check for similarity in parameter estimates within 2 AIC of min

min_t3.a12 <- DP_t3.a12_PECs %>% filter(BIC == min(DP_t3.a12_PECs$BIC, na.rm = T)) %>% select(c(a_log.est:sigma2.est, BIC))

if(nrow(min_t3.a12)>1){
  min_t3.a12 <- DP_t3.a12_PECs %>% filter(strtpt %in% sample(filter(DP_t3.a12_PECs, BIC <= min(DP_t3.a12_PECs$BIC, na.rm = T)+0.1)$strtpt,1))
}

DP_min_t3.a12_flexLV_mod =  mle2(minuslogl = nll.odeint.t3.a1.h1.lv,
                                 method= "Nelder-Mead",
                                 control = list(maxit = 5000),
                                 
                                 start = list(a_log = min_t3.a12$a_log.est, 
                                              h_log = min_t3.a12$h_log.est,
                                              a12 = min_t3.a12$a12.est,
                                              kr_log = min_t3.a12$kr_log.est,
                                              
                                              r1_log = min_t3.a12$r1_log.est, 
                                              alpha11 = min_t3.a12$alpha11.est,
                                              alpha12 = min_t3.a12$alpha12.est,
                                              r2_log = min_t3.a12$r2_log.est, 
                                              alpha22 = min_t3.a12$alpha22.est,
                                              alpha21 = min_t3.a12$alpha21.est,
                                              
                                              c_log = min_t3.a12$c_log.est, 
                                              sigma = min_t3.a12$sigma.est,
                                              sigma2 = min_t3.a12$sigma2.est
                                 ),
                                 
                                 data = with(list(N0 = Dexio_T0.dens.ml,
                                                  Ndead = Dexio_T0.dens.ml - Dexio_T24.dens.ml,
                                                  P = Pred_T0.dens.ml, 
                                                  P.end = Pred_T24.dens.ml,
                                                  M0 = Para_T0.dens.ml,
                                                  M0.end = Para_T24.dens.ml,
                                                  Tt = Time.day),
                                             data = filter(dexio_para_wrk)),
                                 
                                 fixed = list(m = mort,
                                              h12 = 0)
)
save(DP_min_t3.a12_flexLV_mod, file = "DP_t3.a12_flexLV_mod.Rdata")
tidy(DP_min_t3.a12_flexLV_mod)

(nrow(DP_t3.a12_PECs)/nrow(DP_t3.a12_flexLV_pe %>% filter(converge == 0))) * 100 # percent of starting combinations within 2 AIC that 
BIC(DP_min_t3.a12_flexLV_mod)

t_end <- Sys.time()
mod9_time <- t_end-t_start

#### 10. Type III + linear h TIM #####
##### Search Grid #####
t3.h12.pre_sg <- expand.grid(a_log = seq(-3, -1, length.out = 5), 
                             h12 = seq(-0.0005, 0.0001, length.out = 5),
                             h_log = seq(-5, -1, length.out = 5),
                             kr_log = seq(-3, -0.5, length.out = 5),
                             r1_log = c(-2, -1, -0.5), 
                             r2_log = c(-3, -1, -0.5), 
                             alpha11 = c(0.0005, 0.001, 0.005), 
                             alpha22 = c(0.0005, 0.001, 0.005),  
                             alpha12 = c(0.0005, 0.001, 0.005), 
                             alpha21 = c(0.0005, 0.001, 0.005),
                             c_log = c(-15, -10, -5),
                             sigma = c(0.2), 
                             sigma2 = c(0.2)
)
t3.h12.pre_sg$strtpt <- seq(1:nrow(t3.h12.pre_sg))

set.seed(34432)
t3.h12.sim_grid <- t3.h12.pre_sg[sample(t3.h12.pre_sg$strtpt, 98, replace = F), ]

t3.h12.sim_grid$error <- NA

t3.h12.sim_grid$a_log.est <- NA
t3.h12.sim_grid$h_log.est <- NA
t3.h12.sim_grid$h12.est <- NA
t3.h12.sim_grid$kr_log.est <- NA

t3.h12.sim_grid$r1_log.est <- NA
t3.h12.sim_grid$r2_log.est <- NA

t3.h12.sim_grid$alpha11.est <- NA
t3.h12.sim_grid$alpha12.est <- NA
t3.h12.sim_grid$alpha22.est <- NA
t3.h12.sim_grid$alpha21.est <- NA

t3.h12.sim_grid$c_log.est <- NA
t3.h12.sim_grid$sigma.est <- NA
t3.h12.sim_grid$sigma2.est <- NA

t3.h12.sim_grid$converge <- NA
t3.h12.sim_grid$BIC <- NA
t3.h12.sim_grid$logLik <- NA


# Set the computer to run in Parallel #
bigstatsr::nb_cores()

my.cluster <- parallel::makeCluster(
  bigstatsr::nb_cores(),
  type = "FORK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

t_start <- Sys.time()
# Run Grid #
DP_t3.h12_flexLV_pe <- NA # name of file where parallel runs will be save
mcoptions <- list(preschedule=FALSE)

DP_t3.h12_flexLV_pe <- foreach(i = 1:nrow(t3.h12.sim_grid),
                               .combine = 'rbind',
                               .options.multicore = mcoptions,
                               .packages = c('bbmle', 'odeintr', 'tidyverse', 'R.utils'),
                               .errorhandling = "pass") %dopar%
  {
    compile_sys("FRS_t3_a1_h1_lv", FRS.t3.a1.h1.lv, pars = c("a", "a12", "kr", "h", "h12", "r1", "r2","alpha11", "alpha22", "alpha12", "alpha21", "c", "m"), method = "rk54")
    
    r<-tryCatch({ 
      setTimeLimit(35*60)
      fit = mle2(minuslogl = nll.odeint.t3.a1.h1.lv,
                 control = list(maxit = 5000), 
                 method = "Nelder-Mead",
                 
                 start = list(a_log = t3.h12.sim_grid$a_log[[i]], 
                              h12 = t3.h12.sim_grid$h12[[i]], 
                              h_log = t3.h12.sim_grid$h_log[[i]],
                              kr_log = t3.h12.sim_grid$kr_log[[i]],
                              
                              r1_log = t3.h12.sim_grid$r1_log[[i]], 
                              alpha11 = t3.h12.sim_grid$alpha11[[i]],
                              alpha12 = t3.h12.sim_grid$alpha12[[i]],
                              r2_log = t3.h12.sim_grid$r2_log[[i]], 
                              alpha22 = t3.h12.sim_grid$alpha22[[i]],
                              alpha21 = t3.h12.sim_grid$alpha21[[i]],
                              
                              c_log = t3.h12.sim_grid$c_log[[i]],
                              sigma = t3.h12.sim_grid$sigma[[i]],
                              sigma2 = t3.h12.sim_grid$sigma2[[i]]
                 ),
                 
                 data = with(list(N0 = Dexio_T0.dens.ml,
                                  Ndead = Dexio_T0.dens.ml - Dexio_T24.dens.ml,
                                  P = Pred_T0.dens.ml, 
                                  P.end = Pred_T24.dens.ml,
                                  M0 = Para_T0.dens.ml,
                                  M0.end = Para_T24.dens.ml,
                                  Tt = Time.day),
                             data = filter(dexio_para_wrk)),
                 
                 fixed = list(m = mort,
                              a12 = 0) 
      )
      
      if(fit@details$convergence==0) {
        t3.h12.sim_grid$a_log.est[i] <- fit@coef[[1]]
        t3.h12.sim_grid$h_log.est[i] <- fit@coef[[2]]
        t3.h12.sim_grid$h12.est[i] <- fit@coef[[3]]
        t3.h12.sim_grid$kr_log.est[i] <- fit@coef[[4]]
        
        t3.h12.sim_grid$r1_log.est[i] <- fit@coef[[5]]
        t3.h12.sim_grid$alpha11.est[i] <- fit@coef[[6]]
        t3.h12.sim_grid$alpha12.est[i] <- fit@coef[[7]]
        t3.h12.sim_grid$r2_log.est[i] <- fit@coef[[8]]
        t3.h12.sim_grid$alpha22.est[i] <- fit@coef[[9]]
        t3.h12.sim_grid$alpha21.est[i] <- fit@coef[[10]]
        
        t3.h12.sim_grid$c_log.est[i] <- fit@coef[[11]]
        t3.h12.sim_grid$sigma.est [i] <- fit@coef[[12]]
        t3.h12.sim_grid$sigma2.est [i] <- fit@coef[[13]]
        
        t3.h12.sim_grid$converge[i] <- fit@details$convergence
        t3.h12.sim_grid$BIC[i] <- BIC(fit)
        t3.h12.sim_grid$logLik[i] <- logLik(fit)[[1]]
      } else {
        t3.h12.sim_grid$a_log.est[i] <- NaN
        t3.h12.sim_grid$h_log.est[i] <- NaN
        t3.h12.sim_grid$h12.est[i] <- NaN
        
        t3.h12.sim_grid$r1_log.est[i] <- NaN
        t3.h12.sim_grid$alpha11.est[i] <- NaN
        t3.h12.sim_grid$alpha12.est[i] <- NaN
        t3.h12.sim_grid$r2_log.est[i] <- NaN
        t3.h12.sim_grid$alpha22.est[i] <- NaN
        t3.h12.sim_grid$alpha21.est[i] <- NaN
        
        t3.h12.sim_grid$c_log.est[i] <- NaN
        t3.h12.sim_grid$sigma.est [i] <- NaN
        t3.h12.sim_grid$sigma2.est [i] <- NaN
        
        t3.h12.sim_grid$converge[i] <- NaN
        t3.h12.sim_grid$BIC[i] <- NaN
        t3.h12.sim_grid$logLik[i] <- NaN
      }
      
      #Remove fit for next iteration
      rm(fit)
    }, error=function(e) list(message = e$message)) # to print error in console use {cat("ERROR :", paste(conditionMessage(e), "Combo", (i)), "\n")}
    
    if(length(r) > 0){
      t3.h12.sim_grid$error[i] <- r 
    } else {
      t3.h12.sim_grid$error[i] <- NaN
    }  
    
    return(t3.h12.sim_grid[i,])
    
  }

parallel::stopCluster(cl = my.cluster)
stopImplicitCluster()

DP_t3.h12_flexLV_pe %>% arrange(BIC)

save(DP_t3.h12_flexLV_pe, file = "DP_t3.h12_flexLV_pe.RData")

###### Best Model #####
(DP_t3.h12_PECs <- DP_t3.h12_flexLV_pe %>% 
   filter(converge == 0) %>% 
   filter(BIC <= min(DP_t3.h12_flexLV_pe$BIC, na.rm = T)+2)) #check for similarity in parameter estimates within 2 AIC of min

min_t3.h12 <- DP_t3.h12_PECs %>% filter(BIC == min(DP_t3.h12_PECs$BIC, na.rm = T)) %>% select(c(a_log.est:sigma2.est, BIC))

if(nrow(min_t3.h12)>1){
  min_t3.h12 <- DP_t3.h12_PECs %>% filter(strtpt %in% sample(filter(DP_t3.h12_PECs, BIC <= min(DP_t3.h12_PECs$BIC, na.rm = T)+0.1)$strtpt,1))
}

DP_min_t3.h12_flexLV_mod =  mle2(minuslogl = nll.odeint.t3.a1.h1.lv,
                                 method= "Nelder-Mead",
                                 control = list(maxit = 5000),
                                 
                                 start = list(a_log = min_t3.h12$a_log.est, 
                                              h_log = min_t3.h12$h_log.est,
                                              h12 = min_t3.h12$h12.est,
                                              kr_log = min_t3.h12$kr_log.est,
                                              
                                              r1_log = min_t3.h12$r1_log.est, 
                                              alpha11 = min_t3.h12$alpha11.est,
                                              alpha12 = min_t3.h12$alpha12.est,
                                              r2_log = min_t3.h12$r2_log.est, 
                                              alpha22 = min_t3.h12$alpha22.est,
                                              alpha21 = min_t3.h12$alpha21.est,
                                              
                                              c_log = min_t3.h12$c_log.est, 
                                              sigma = min_t3.h12$sigma.est,
                                              sigma2 = min_t3.h12$sigma2.est
                                 ),
                                 
                                 data = with(list(N0 = Dexio_T0.dens.ml,
                                                  Ndead = Dexio_T0.dens.ml - Dexio_T24.dens.ml,
                                                  P = Pred_T0.dens.ml, 
                                                  P.end = Pred_T24.dens.ml,
                                                  M0 = Para_T0.dens.ml,
                                                  M0.end = Para_T24.dens.ml,
                                                  Tt = Time.day),
                                             data = filter(dexio_para_wrk)),
                                 
                                 fixed = list(m = mort,
                                              a12=0)
)
save(DP_min_t3.h12_flexLV_mod, file = "DP_t3.h12_flexLV_mod.Rdata")
tidy(DP_min_t3.h12_flexLV_mod)

(nrow(DP_t3.h12_PECs)/nrow(DP_t3.h12_flexLV_pe %>% filter(converge == 0))) * 100 # percent of starting combinations within 2 AIC that 
BIC(DP_min_t3.h12_flexLV_mod)

t_end <- Sys.time()
mod10_time <- t_end-t_start

#### 11. Type III + linear a TIM  + linear h TIM #####
##### Search Grid #####
t3.a12.h12.pre_sg <- expand.grid(a_log = seq(-3, -1, length.out = 5), 
                                 a12 = seq(-0.0005, 0.0001, length.out = 5),
                                 h12 = seq(-0.0005, 0.0001, length.out = 5),
                                 h_log = seq(-5, -1, length.out = 5),
                                 kr_log = seq(-3, -0.5, length.out = 5),
                                 r1_log = c(-2, -1, -0.5), 
                                 r2_log = c(-3, -1, -0.5), 
                                 alpha11 = c(0.0005, 0.001, 0.005), 
                                 alpha22 = c(0.0005, 0.001, 0.005),  
                                 alpha12 = c(0.0005, 0.001, 0.005), 
                                 alpha21 = c(0.0005, 0.001, 0.005),
                                 c_log = c(-15, -10, -5),
                                 sigma = c(0.2), 
                                 sigma2 = c(0.2)
)
t3.a12.h12.pre_sg$strtpt <- seq(1:nrow(t3.a12.h12.pre_sg))

set.seed(31394)
t3.a12.h12.sim_grid <- t3.a12.h12.pre_sg[sample(t3.a12.h12.pre_sg$strtpt, 98, replace = F), ]

t3.a12.h12.sim_grid$error <- NA

t3.a12.h12.sim_grid$a_log.est <- NA
t3.a12.h12.sim_grid$a12.est <- NA
t3.a12.h12.sim_grid$h_log.est <- NA
t3.a12.h12.sim_grid$h12.est <- NA
t3.a12.h12.sim_grid$kr_log.est <- NA

t3.a12.h12.sim_grid$r1_log.est <- NA
t3.a12.h12.sim_grid$r2_log.est <- NA

t3.a12.h12.sim_grid$alpha11.est <- NA
t3.a12.h12.sim_grid$alpha12.est <- NA
t3.a12.h12.sim_grid$alpha22.est <- NA
t3.a12.h12.sim_grid$alpha21.est <- NA

t3.a12.h12.sim_grid$c_log.est <- NA
t3.a12.h12.sim_grid$sigma.est <- NA
t3.a12.h12.sim_grid$sigma2.est <- NA

t3.a12.h12.sim_grid$converge <- NA
t3.a12.h12.sim_grid$BIC <- NA
t3.a12.h12.sim_grid$logLik <- NA



# Set the computer to run in Parallel #
bigstatsr::nb_cores()

my.cluster <- parallel::makeCluster(
  bigstatsr::nb_cores(),
  type = "FORK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

t_start <- Sys.time()
# Run Grid #
DP_t3.a12.h12_flexLV_pe <- NA # name of file where parallel runs will be save
mcoptions <- list(preschedule=FALSE)

DP_t3.a12.h12_flexLV_pe <- foreach(i = 1:nrow(t3.a12.h12.sim_grid),
                                   .combine = 'rbind',
                                   .options.multicore = mcoptions,
                                   .packages = c('bbmle', 'odeintr', 'tidyverse', 'R.utils'),
                                   .errorhandling = "pass") %dopar%
  {
    compile_sys("FRS_t3_a1_h1_lv", FRS.t3.a1.h1.lv, pars = c("a", "a12", "kr", "h", "h12", "r1", "r2","alpha11", "alpha22", "alpha12", "alpha21", "c", "m"), method = "rk54")
    
    r<-tryCatch({ 
      setTimeLimit(35*60)
      fit = mle2(minuslogl = nll.odeint.t3.a1.h1.lv,
                 control = list(maxit = 5000), 
                 method = "Nelder-Mead",
                 
                 start = list(a_log = t3.a12.h12.sim_grid$a_log[[i]], 
                              a12 = t3.a12.h12.sim_grid$a12[[i]], 
                              h_log = t3.a12.h12.sim_grid$h_log[[i]],
                              h12 = t3.a12.h12.sim_grid$h12[[i]], 
                              kr_log = t3.a12.h12.sim_grid$kr_log[[i]],
                              
                              r1_log = t3.a12.h12.sim_grid$r1_log[[i]], 
                              alpha11 = t3.a12.h12.sim_grid$alpha11[[i]],
                              alpha12 = t3.a12.h12.sim_grid$alpha12[[i]],
                              r2_log = t3.a12.h12.sim_grid$r2_log[[i]], 
                              alpha22 = t3.a12.h12.sim_grid$alpha22[[i]],
                              alpha21 = t3.a12.h12.sim_grid$alpha21[[i]],
                              
                              c_log = t3.a12.h12.sim_grid$c_log[[i]], 
                              sigma = t3.a12.h12.sim_grid$sigma[[i]],
                              sigma2 = t3.a12.h12.sim_grid$sigma2[[i]]
                 ),
                 
                 data = with(list(N0 = Dexio_T0.dens.ml,
                                  Ndead = Dexio_T0.dens.ml - Dexio_T24.dens.ml,
                                  P = Pred_T0.dens.ml, 
                                  P.end = Pred_T24.dens.ml,
                                  M0 = Para_T0.dens.ml,
                                  M0.end = Para_T24.dens.ml,
                                  Tt = Time.day),
                             data = filter(dexio_para_wrk)),
                 
                 fixed = list(m = mort) 
      )
      
      if(fit@details$convergence==0) {
        t3.a12.h12.sim_grid$a_log.est[i] <- fit@coef[[1]]
        t3.a12.h12.sim_grid$a12.est[i] <- fit@coef[[2]]
        t3.a12.h12.sim_grid$h_log.est[i] <- fit@coef[[3]]
        t3.a12.h12.sim_grid$h12.est[i] <- fit@coef[[4]]
        t3.a12.h12.sim_grid$kr_log.est[i] <- fit@coef[[5]]
        
        t3.a12.h12.sim_grid$r1_log.est[i] <- fit@coef[[6]]
        t3.a12.h12.sim_grid$alpha11.est[i] <- fit@coef[[7]]
        t3.a12.h12.sim_grid$alpha12.est[i] <- fit@coef[[8]]
        t3.a12.h12.sim_grid$r2_log.est[i] <- fit@coef[[9]]
        t3.a12.h12.sim_grid$alpha22.est[i] <- fit@coef[[10]]
        t3.a12.h12.sim_grid$alpha21.est[i] <- fit@coef[[11]]
        
        t3.a12.h12.sim_grid$c_log.est[i] <- fit@coef[[12]]
        t3.a12.h12.sim_grid$sigma.est [i] <- fit@coef[[13]]
        t3.a12.h12.sim_grid$sigma2.est [i] <- fit@coef[[14]]
        
        t3.a12.h12.sim_grid$converge[i] <- fit@details$convergence
        t3.a12.h12.sim_grid$BIC[i] <- BIC(fit)
        t3.a12.h12.sim_grid$logLik[i] <- logLik(fit)[[1]]
      } else {
        t3.a12.h12.sim_grid$a_log.est[i] <- NaN
        t3.a12.h12.sim_grid$a12.est[i] <- NaN
        t3.a12.h12.sim_grid$h_log.est[i] <- NaN
        t3.a12.h12.sim_grid$a12.est[i] <- NaN
        t3.a12.h12.sim_grid$kr_log.est[i] <- NaN
        
        t3.a12.h12.sim_grid$r1_log.est[i] <- NaN
        t3.a12.h12.sim_grid$alpha11.est[i] <- NaN
        t3.a12.h12.sim_grid$alpha12.est[i] <- NaN
        t3.a12.h12.sim_grid$r2_log.est[i] <- NaN
        t3.a12.h12.sim_grid$alpha22.est[i] <- NaN
        t3.a12.h12.sim_grid$alpha21.est[i] <- NaN
        
        t3.a12.h12.sim_grid$c_log.est[i] <- NaN
        t3.a12.h12.sim_grid$sigma.est [i] <- NaN
        t3.a12.h12.sim_grid$sigma2.est [i] <- NaN
        
        t3.a12.h12.sim_grid$converge[i] <- NaN
        t3.a12.h12.sim_grid$BIC[i] <- NaN
        t3.a12.h12.sim_grid$logLik[i] <- NaN
      }
      
      #Remove fit for next iteration
      rm(fit)
    }, error=function(e) list(message = e$message)) # to print error in console use {cat("ERROR :", paste(conditionMessage(e), "Combo", (i)), "\n")}
    
    if(length(r) > 0){
      t3.a12.h12.sim_grid$error[i] <- r 
    } else {
      t3.a12.h12.sim_grid$error[i] <- NaN
    }  
    
    return(t3.a12.h12.sim_grid[i,])
    
  }

parallel::stopCluster(cl = my.cluster)
stopImplicitCluster()

DP_t3.a12.h12_flexLV_pe %>% arrange(BIC)

save(DP_t3.a12.h12_flexLV_pe, file = "DP_t3.a12.h12_flexLV_pe.RData")

###### Best Model #####
(DP_t3.a12.h12_PECs <- DP_t3.a12.h12_flexLV_pe %>% 
   filter(converge == 0) %>% 
   filter(BIC <= min(DP_t3.a12.h12_flexLV_pe$BIC, na.rm = T)+2)) #check for similarity in parameter estimates within 2 AIC of min

min_t3.a12.h12 <- DP_t3.a12.h12_PECs %>% filter(BIC == min(DP_t3.a12.h12_PECs$BIC, na.rm = T)) %>% select(c(a_log.est:sigma2.est, BIC))

if(nrow(min_t3.a12.h12)>1){
  min_t3.a12.h12 <- DP_t3.a12.h12_PECs %>% filter(strtpt %in% sample(filter(DP_t3.a12.h12_PECs, BIC <= min(DP_t3.a12.h12_PECs$BIC, na.rm = T)+0.1)$strtpt,1))
}

DP_min_t3.a12.h12_flexLV_mod =  mle2(minuslogl = nll.odeint.t3.a1.h1.lv,
                                     method= "Nelder-Mead",
                                     control = list(maxit = 5000),
                                     
                                     start = list(a_log = min_t3.a12.h12$a_log.est, 
                                                  a12 = min_t3.a12.h12$a12.est,
                                                  h_log = min_t3.a12.h12$h_log.est,
                                                  h12 = min_t3.a12.h12$h12.est,
                                                  kr_log = min_t3.a12.h12$kr_log.est,
                                                  
                                                  r1_log = min_t3.a12.h12$r1_log.est, 
                                                  alpha11 = min_t3.a12.h12$alpha11.est,
                                                  alpha12 = min_t3.a12.h12$alpha12.est,
                                                  r2_log = min_t3.a12.h12$r2_log.est, 
                                                  alpha22 = min_t3.a12.h12$alpha22.est,
                                                  alpha21 = min_t3.a12.h12$alpha21.est,
                                                  
                                                  c_log = min_t3.a12.h12$c_log.est, 
                                                  sigma = min_t3.a12.h12$sigma.est,
                                                  sigma2 = min_t3.a12.h12$sigma2.est
                                     ),
                                     
                                     data = with(list(N0 = Dexio_T0.dens.ml,
                                                      Ndead = Dexio_T0.dens.ml - Dexio_T24.dens.ml,
                                                      P = Pred_T0.dens.ml, 
                                                      P.end = Pred_T24.dens.ml,
                                                      M0 = Para_T0.dens.ml,
                                                      M0.end = Para_T24.dens.ml,
                                                      Tt = Time.day),
                                                 data = filter(dexio_para_wrk)),
                                     
                                     fixed = list(m = mort)
)
save(DP_min_t3.a12.h12_flexLV_mod, file = "DP_t3.a12.h12_flexLV_mod.Rdata")
tidy(DP_min_t3.a12.h12_flexLV_mod)

(nrow(DP_t3.a12.h12_PECs)/nrow(DP_t3.a12.h12_flexLV_pe %>% filter(converge == 0))) * 100 # percent of starting combinations within 2 AIC that 
BIC(DP_min_t3.a12.h12_flexLV_mod)

t_end <- Sys.time()
mod11_time <- t_end-t_start

#### 12. Type III + van Veen #####
##### Search Grid #####
t3.vv.pre_sg <- expand.grid(a_log = seq(-3, -1, length.out = 5), 
                            h_log = seq(-5, -1, length.out = 5), 
                            log_w = seq(-10, -5, length.out = 5), 
                            kr_log = seq(-3, -0.5, length.out = 5),
                            r1_log = c(-2, -1, -0.5), 
                            r2_log = c(-3, -1, -0.5), 
                            alpha11 = c(0.0005, 0.001, 0.005), 
                            alpha22 = c(0.0005, 0.001, 0.005),  
                            alpha12 = c(0.0005, 0.001, 0.005), 
                            alpha21 = c(0.0005, 0.001, 0.005),
                            c_log = c(-15, -10, -5),
                            sigma = c(0.2), 
                            sigma2 = c(0.2)
)
t3.vv.pre_sg$strtpt <- seq(1:nrow(t3.vv.pre_sg))

set.seed(96542)
t3.vv.sim_grid <- t3.vv.pre_sg[sample(t3.vv.pre_sg$strtpt, 98, replace = F), ]

t3.vv.sim_grid$error <- NA

t3.vv.sim_grid$a_log.est <- NA
t3.vv.sim_grid$h_log.est <- NA
t3.vv.sim_grid$log_w.est <- NA
t3.vv.sim_grid$kr_log.est <- NA

t3.vv.sim_grid$r1_log.est <- NA
t3.vv.sim_grid$r2_log.est <- NA

t3.vv.sim_grid$alpha11.est <- NA
t3.vv.sim_grid$alpha12.est <- NA
t3.vv.sim_grid$alpha22.est <- NA
t3.vv.sim_grid$alpha21.est <- NA

t3.vv.sim_grid$c_log.est <- NA
t3.vv.sim_grid$sigma.est <- NA
t3.vv.sim_grid$sigma2.est <- NA

t3.vv.sim_grid$converge <- NA
t3.vv.sim_grid$BIC <- NA
t3.vv.sim_grid$logLik <- NA

# Set the computer to run in Parallel #
bigstatsr::nb_cores()

my.cluster <- parallel::makeCluster(
  bigstatsr::nb_cores(),
  type = "FORK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

t_start <- Sys.time()
# Run Grid #
DP_t3.vv_flexLV_pe <- NA # name of file where parallel runs will be save
mcoptions <- list(preschedule=FALSE)

DP_t3.vv_flexLV_pe <- foreach(i = 1:nrow(t3.vv.sim_grid),
                              .combine = 'rbind',
                              .options.multicore = mcoptions,
                              .packages = c('bbmle', 'odeintr', 'tidyverse', 'R.utils'),
                              .errorhandling = "pass") %dopar%
  {
    compile_sys("FRS_T3_h1_VV_lv", FRS.T3.h1.VV.lv, pars = c("a", "kr", "h", "h12", "w", "r1", "r2","alpha11", "alpha22", "alpha12", "alpha21", "c", "m"), method = "rk54")
    
    r<-tryCatch({ 
      setTimeLimit(35*60)
      fit = mle2(minuslogl = nll.odeint.T3.h1.VV.lv,
                 control = list(maxit = 5000), 
                 method = "Nelder-Mead",
                 
                 start = list(a_log = t3.vv.sim_grid$a_log[[i]], 
                              h_log = t3.vv.sim_grid$h_log[[i]],
                              log_w = t3.vv.pre_sg$log_w[[i]],
                              kr_log = t3.vv.sim_grid$kr_log[[i]],
                              
                              r1_log = t3.vv.sim_grid$r1_log[[i]], 
                              alpha11 = t3.vv.sim_grid$alpha11[[i]],
                              alpha12 = t3.vv.sim_grid$alpha12[[i]],
                              r2_log = t3.vv.sim_grid$r2_log[[i]], 
                              alpha22 = t3.vv.sim_grid$alpha22[[i]],
                              alpha21 = t3.vv.sim_grid$alpha21[[i]],
                              
                              c_log = t3.vv.sim_grid$c_log[[i]], 
                              sigma = t3.vv.sim_grid$sigma[[i]],
                              sigma2 = t3.vv.sim_grid$sigma2[[i]]
                 ),
                 
                 data = with(list(N0 = Dexio_T0.dens.ml,
                                  Ndead = Dexio_T0.dens.ml - Dexio_T24.dens.ml,
                                  P = Pred_T0.dens.ml, 
                                  P.end = Pred_T24.dens.ml,
                                  M0 = Para_T0.dens.ml,
                                  M0.end = Para_T24.dens.ml,
                                  Tt = Time.day),
                             data = filter(dexio_para_wrk)),
                 
                 fixed = list(m = mort,
                              h12 = 0) 
      )
      
      if(fit@details$convergence==0) {
        t3.vv.sim_grid$a_log.est[i] <- fit@coef[[1]]
        t3.vv.sim_grid$h_log.est[i] <- fit@coef[[2]]
        t3.vv.sim_grid$kr_log.est[i] <- fit@coef[[3]]
        t3.vv.sim_grid$log_w.est[i] <- fit@coef[[4]]
        
        
        t3.vv.sim_grid$r1_log.est[i] <- fit@coef[[5]]
        t3.vv.sim_grid$alpha11.est[i] <- fit@coef[[6]]
        t3.vv.sim_grid$alpha12.est[i] <- fit@coef[[7]]
        t3.vv.sim_grid$r2_log.est[i] <- fit@coef[[8]]
        t3.vv.sim_grid$alpha22.est[i] <- fit@coef[[9]]
        t3.vv.sim_grid$alpha21.est[i] <- fit@coef[[10]]
        
        t3.vv.sim_grid$c_log.est[i] <- fit@coef[[11]]
        t3.vv.sim_grid$sigma.est [i] <- fit@coef[[12]]
        t3.vv.sim_grid$sigma2.est [i] <- fit@coef[[13]]
        
        t3.vv.sim_grid$converge[i] <- fit@details$convergence
        t3.vv.sim_grid$BIC[i] <- BIC(fit)
        t3.vv.sim_grid$logLik[i] <- logLik(fit)[[1]]
      } else {
        t3.vv.sim_grid$a_log.est[i] <- NaN
        t3.vv.sim_grid$h_log.est[i] <- NaN
        t3.vv.sim_grid$log_w.est[i] <- NaN
        t3.vv.sim_grid$kr_log.est[i] <- NaN
        
        t3.vv.sim_grid$r1_log.est[i] <- NaN
        t3.vv.sim_grid$alpha11.est[i] <- NaN
        t3.vv.sim_grid$alpha12.est[i] <- NaN
        t3.vv.sim_grid$r2_log.est[i] <- NaN
        t3.vv.sim_grid$alpha22.est[i] <- NaN
        t3.vv.sim_grid$alpha21.est[i] <- NaN
        
        t3.vv.sim_grid$c_log.est[i] <- NaN
        t3.vv.sim_grid$sigma.est [i] <- NaN
        t3.vv.sim_grid$sigma2.est [i] <- NaN
        
        t3.vv.sim_grid$converge[i] <- NaN
        t3.vv.sim_grid$BIC[i] <- NaN
        t3.vv.sim_grid$logLik[i] <- NaN
      }
      
      #Remove fit for next iteration
      rm(fit)
    }, error=function(e) list(message = e$message)) # to print error in console use {cat("ERROR :", paste(conditionMessage(e), "Combo", (i)), "\n")}
    
    if(length(r) > 0){
      t3.vv.sim_grid$error[i] <- r 
    } else {
      t3.vv.sim_grid$error[i] <- NaN
    }  
    
    return(t3.vv.sim_grid[i,])
    
  }

parallel::stopCluster(cl = my.cluster)
stopImplicitCluster()

DP_t3.vv_flexLV_pe %>% arrange(BIC)

save(DP_t3.vv_flexLV_pe, file = "DP_t3.vv_flexLV_pe.RData")

###### Best Model #####
(DP_t3.vv_PECs <- DP_t3.vv_flexLV_pe %>% 
   filter(converge == 0) %>% 
   filter(BIC <= min(DP_t3.vv_flexLV_pe$BIC, na.rm = T)+2)) #check for similarity in parameter estimates within 2 AIC of min

min_t3.vv <- DP_t3.vv_PECs %>% filter(BIC == min(DP_t3.vv_PECs$BIC, na.rm = T)) %>% select(c(a_log.est:sigma2.est, BIC))

if(nrow(min_t3.vv)>1){
  min_t3.vv <- DP_t3.vv_PECs %>% filter(strtpt %in% sample(filter(DP_t3.vv_PECs, BIC <= min(DP_t3.vv_PECs$BIC, na.rm = T)+0.1)$strtpt,1))
}

DP_min_t3.vv_flexLV_mod =  mle2(minuslogl = nll.odeint.T3.h1.VV.lv,
                                method= "Nelder-Mead",
                                control = list(maxit = 5000),
                                
                                start = list(a_log = min_t3.vv$a_log.est, 
                                             h_log = min_t3.vv$h_log.est,
                                             log_w = min_t3.vv$log_w.est,
                                             kr_log = min_t3.vv$kr_log.est,
                                             
                                             r1_log = min_t3.vv$r1_log.est, 
                                             alpha11 = min_t3.vv$alpha11.est,
                                             alpha12 = min_t3.vv$alpha12.est,
                                             r2_log = min_t3.vv$r2_log.est, 
                                             alpha22 = min_t3.vv$alpha22.est,
                                             alpha21 = min_t3.vv$alpha21.est,
                                             
                                             c_log = min_t3.vv$c_log.est, sigma = min_t3.vv$sigma.est,
                                             sigma2 = min_t3.vv$sigma2.est
                                ),
                                
                                data = with(list(N0 = Dexio_T0.dens.ml,
                                                 Ndead = Dexio_T0.dens.ml - Dexio_T24.dens.ml,
                                                 P = Pred_T0.dens.ml, 
                                                 P.end = Pred_T24.dens.ml,
                                                 M0 = Para_T0.dens.ml,
                                                 M0.end = Para_T24.dens.ml,
                                                 Tt = Time.day),
                                            data = filter(dexio_para_wrk)),
                                
                                fixed = list(m = mort,
                                             h12=0)
)
save(DP_min_t3.vv_flexLV_mod, file = "DP_t3.vv_flexLV_mod.Rdata")
tidy(DP_min_t3.vv_flexLV_mod)

(nrow(DP_t3.vv_PECs)/nrow(DP_t3.vv_flexLV_pe %>% filter(converge == 0))) * 100 # percent of starting combinations within 2 AIC that 
BIC(DP_min_t3.vv_flexLV_mod)

t_end <- Sys.time()
mod12_time <- t_end-t_start

#### 13. Type III + Crowley-Martin #####
##### Search Grid #####
t3.cm.pre_sg <- expand.grid(a_log = seq(-3, -1, length.out = 5),
                            h_log = seq(-5, -1, length.out = 5),
                            log_w = seq(-10, -5, length.out = 5),
                            kr_log = seq(-3, -0.5, length.out = 5),
                            r1_log = c(-2, -1, -0.5), 
                            r2_log = c(-3, -1, -0.5), 
                            alpha11 = c(0.0005, 0.001, 0.005), 
                            alpha22 = c(0.0005, 0.001, 0.005),  
                            alpha12 = c(0.0005, 0.001, 0.005), 
                            alpha21 = c(0.0005, 0.001, 0.005),
                            c_log = c(-15, -10, -5),
                            sigma = c(0.2), 
                            sigma2 = c(0.2)
)
t3.cm.pre_sg$strtpt <- seq(1:nrow(t3.cm.pre_sg))

set.seed(23396)
t3.cm.sim_grid <- t3.cm.pre_sg[sample(t3.cm.pre_sg$strtpt, 98, replace = F), ]

t3.cm.sim_grid$error <- NA

t3.cm.sim_grid$a_log.est <- NA
t3.cm.sim_grid$h_log.est <- NA
t3.cm.sim_grid$log_w.est <- NA
t3.cm.sim_grid$kr_log.est <- NA

t3.cm.sim_grid$r1_log.est <- NA
t3.cm.sim_grid$r2_log.est <- NA

t3.cm.sim_grid$alpha11.est <- NA
t3.cm.sim_grid$alpha12.est <- NA
t3.cm.sim_grid$alpha22.est <- NA
t3.cm.sim_grid$alpha21.est <- NA

t3.cm.sim_grid$c_log.est <- NA
t3.cm.sim_grid$sigma.est <- NA
t3.cm.sim_grid$sigma2.est <- NA

t3.cm.sim_grid$converge <- NA
t3.cm.sim_grid$BIC <- NA
t3.cm.sim_grid$logLik <- NA


# Set the computer to run in Parallel #
bigstatsr::nb_cores()

my.cluster <- parallel::makeCluster(
  bigstatsr::nb_cores(),
  type = "FORK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

t_start <- Sys.time()
# Run Grid #
DP_t3.cm_flexLV_pe <- NA # name of file where parallel runs will be save
mcoptions <- list(preschedule=FALSE)

DP_t3.cm_flexLV_pe <- foreach(i = 1:nrow(t3.cm.sim_grid),
                              .combine = 'rbind',
                              .options.multicore = mcoptions,
                              .packages = c('bbmle', 'odeintr', 'tidyverse', 'R.utils'),
                              .errorhandling = "pass") %dopar%
  {
    compile_sys("FRS_T3_CM_lv", FRS.T3.CM.lv, pars = c("a", "kr", "h", "w", "r1", "r2","alpha11", "alpha22", "alpha12", "alpha21", "c", "m"), method = "rk54")
    
    r<-tryCatch({ 
      setTimeLimit(35*60)
      fit = mle2(minuslogl = nll.odeint.T3.CM.lv,
                 control = list(maxit = 5000), 
                 method = "Nelder-Mead",
                 
                 start = list(a_log = t3.cm.sim_grid$a_log[[i]], 
                              h_log = t3.cm.sim_grid$h_log[[i]],
                              log_w = t3.cm.pre_sg$log_w[[i]],
                              kr_log = t3.cm.sim_grid$kr_log[[i]],
                              
                              
                              r1_log = t3.cm.sim_grid$r1_log[[i]], 
                              alpha11 = t3.cm.sim_grid$alpha11[[i]],
                              alpha12 = t3.cm.sim_grid$alpha12[[i]],
                              r2_log = t3.cm.sim_grid$r2_log[[i]], 
                              alpha22 = t3.cm.sim_grid$alpha22[[i]],
                              alpha21 = t3.cm.sim_grid$alpha21[[i]],
                              
                              c_log = t3.cm.sim_grid$c_log[[i]],
                              sigma = t3.cm.sim_grid$sigma[[i]],
                              sigma2 = t3.cm.sim_grid$sigma2[[i]]
                 ),
                 
                 data = with(list(N0 = Dexio_T0.dens.ml,
                                  Ndead = Dexio_T0.dens.ml - Dexio_T24.dens.ml,
                                  P = Pred_T0.dens.ml, 
                                  P.end = Pred_T24.dens.ml,
                                  M0 = Para_T0.dens.ml,
                                  M0.end = Para_T24.dens.ml,
                                  Tt = Time.day),
                             data = filter(dexio_para_wrk)),
                 
                 fixed = list(m = mort) 
      )
      
      if(fit@details$convergence==0) {
        t3.cm.sim_grid$a_log.est[i] <- fit@coef[[1]]
        t3.cm.sim_grid$h_log.est[i] <- fit@coef[[2]]
        t3.cm.sim_grid$kr_log.est[i] <- fit@coef[[3]]
        t3.cm.sim_grid$log_w.est[i] <- fit@coef[[4]]
        
        t3.cm.sim_grid$r1_log.est[i] <- fit@coef[[5]]
        t3.cm.sim_grid$alpha11.est[i] <- fit@coef[[6]]
        t3.cm.sim_grid$alpha12.est[i] <- fit@coef[[7]]
        t3.cm.sim_grid$r2_log.est[i] <- fit@coef[[8]]
        t3.cm.sim_grid$alpha22.est[i] <- fit@coef[[9]]
        t3.cm.sim_grid$alpha21.est[i] <- fit@coef[[10]]
        
        t3.cm.sim_grid$c_log.est[i] <- fit@coef[[11]]
        t3.cm.sim_grid$sigma.est [i] <- fit@coef[[12]]
        t3.cm.sim_grid$sigma2.est [i] <- fit@coef[[13]]
        
        t3.cm.sim_grid$converge[i] <- fit@details$convergence
        t3.cm.sim_grid$BIC[i] <- BIC(fit)
        t3.cm.sim_grid$logLik[i] <- logLik(fit)[[1]]
      } else {
        t3.cm.sim_grid$a_log.est[i] <- NaN
        t3.cm.sim_grid$h_log.est[i] <- NaN
        t3.cm.sim_grid$log_w.est[i] <- NaN
        t3.cm.sim_grid$kr_log.est[i] <- NaN
        
        t3.cm.sim_grid$r1_log.est[i] <- NaN
        t3.cm.sim_grid$alpha11.est[i] <- NaN
        t3.cm.sim_grid$alpha12.est[i] <- NaN
        t3.cm.sim_grid$r2_log.est[i] <- NaN
        t3.cm.sim_grid$alpha22.est[i] <- NaN
        t3.cm.sim_grid$alpha21.est[i] <- NaN
        
        t3.cm.sim_grid$c_log.est[i] <- NaN
        t3.cm.sim_grid$sigma.est [i] <- NaN
        t3.cm.sim_grid$sigma2.est [i] <- NaN
        
        t3.cm.sim_grid$converge[i] <- NaN
        t3.cm.sim_grid$BIC[i] <- NaN
        t3.cm.sim_grid$logLik[i] <- NaN
      }
      
      #Remove fit for next iteration
      rm(fit)
    }, error=function(e) list(message = e$message)) # to print error in console use {cat("ERROR :", paste(conditionMessage(e), "Combo", (i)), "\n")}
    
    if(length(r) > 0){
      t3.cm.sim_grid$error[i] <- r 
    } else {
      t3.cm.sim_grid$error[i] <- NaN
    }  
    
    return(t3.cm.sim_grid[i,])
    
  }

parallel::stopCluster(cl = my.cluster)
stopImplicitCluster()

DP_t3.cm_flexLV_pe %>% arrange(BIC)

save(DP_t3.cm_flexLV_pe, file = "DP_t3.cm_flexLV_pe.RData")

###### Best Model #####
(DP_t3.cm_PECs <- DP_t3.cm_flexLV_pe %>% 
   filter(converge == 0) %>% 
   filter(BIC <= min(DP_t3.cm_flexLV_pe$BIC, na.rm = T)+2)) #check for similarity in parameter estimates within 2 AIC of min

min_t3.cm <- DP_t3.cm_PECs %>% filter(BIC == min(DP_t3.cm_PECs$BIC, na.rm = T)) %>% select(c(a_log.est:sigma2.est, BIC))

if(nrow(min_t3.cm)>1){
  min_t3.cm <- DP_t3.cm_PECs %>% filter(strtpt %in% sample(filter(DP_t3.cm_PECs, BIC <= min(DP_t3.cm_PECs$BIC, na.rm = T)+0.1)$strtpt,1))
}

DP_min_t3.cm_flexLV_mod =  mle2(minuslogl = nll.odeint.T3.CM.lv,
                                method= "Nelder-Mead",
                                control = list(maxit = 5000),
                                
                                start = list(a_log = min_t3.cm$a_log.est, 
                                             h_log = min_t3.cm$h_log.est,
                                             log_w = min_t3.cm$log_w.est,
                                             kr_log = min_t3.cm$kr_log.est,
                                             
                                             r1_log = min_t3.cm$r1_log.est, 
                                             alpha11 = min_t3.cm$alpha11.est,
                                             alpha12 = min_t3.cm$alpha12.est,
                                             r2_log = min_t3.cm$r2_log.est, 
                                             alpha22 = min_t3.cm$alpha22.est,
                                             alpha21 = min_t3.cm$alpha21.est,
                                             
                                             c_log = min_t3.cm$c_log.est, 
                                             sigma = min_t3.cm$sigma.est,
                                             sigma2 = min_t3.cm$sigma2.est
                                ),
                                
                                data = with(list(N0 = Dexio_T0.dens.ml,
                                                 Ndead = Dexio_T0.dens.ml - Dexio_T24.dens.ml,
                                                 P = Pred_T0.dens.ml, 
                                                 P.end = Pred_T24.dens.ml,
                                                 M0 = Para_T0.dens.ml,
                                                 M0.end = Para_T24.dens.ml,
                                                 Tt = Time.day),
                                            data = filter(dexio_para_wrk)),
                                
                                fixed = list(m = mort)
)
save(DP_min_t3.cm_flexLV_mod, file = "DP_t3.cm_flexLV_mod.Rdata")
tidy(DP_min_t3.cm_flexLV_mod)

(nrow(DP_t3.cm_PECs)/nrow(DP_t3.cm_flexLV_pe %>% filter(converge == 0))) * 100 # percent of starting combinations within 2 AIC that 
BIC(DP_min_t3.cm_flexLV_mod)

t_end <- Sys.time()
mod13_time <- t_end-t_start

#### BIC Table ####
DP_BIC<-as.data.frame(BIC(DP_min_T2_flexLV_mod,
                          DP_min_t2.a12_flexLV_mod,
                          DP_min_t2.h12_flexLV_mod,
                          DP_min_t2.a12.h12_flexLV_mod,
                          DP_min_t2.vv_flexLV_mod,
                          DP_min_t2.vv.h12_flexLV_mod,
                          DP_min_t2.cm_flexLV_mod,
                          DP_min_t3_flexLV_mod,
                          DP_min_t3.a12_flexLV_mod,
                          DP_min_t3.h12_flexLV_mod,
                          DP_min_t3.a12.h12_flexLV_mod,
                          DP_min_t3.vv_flexLV_mod,
                          DP_min_t3.cm_flexLV_mod))
DP_BIC$model<-seq(1, nrow(DP_BIC), by =1); DP_BIC$dBIC<-round(DP_BIC$BIC-min(DP_BIC$BIC, na.rm = T),1)

DP_BIC %>% arrange(BIC)

(total_time <-sum(mod1_time, mod2_time, mod3_time, mod4_time, mod5_time, mod6_time, mod7_time, mod8_time, mod9_time, mod10_time, mod11_time, mod12_time, mod13_time))

save.image(file = "dexioTIM_fit_flex.RData")

#### Parameter Table
pe_list <- list(coef(DP_min_T2_flexLV_mod),
                coef(DP_min_t2.a12_flexLV_mod),
                coef(DP_min_t2.h12_flexLV_mod),
                coef(DP_min_t2.a12.h12_flexLV_mod),
                coef(DP_min_t2.vv_flexLV_mod),
                coef(DP_min_t2.vv.h12_flexLV_mod),
                coef(DP_min_t2.cm_flexLV_mod),
                coef(DP_min_t3_flexLV_mod),
                coef(DP_min_t3.a12_flexLV_mod),
                coef(DP_min_t3.h12_flexLV_mod),
                coef(DP_min_t3.a12.h12_flexLV_mod),
                coef(DP_min_t3.vv_flexLV_mod),
                coef(DP_min_t3.cm_flexLV_mod)
)

pe_table <- as.data.frame(do.call(rbind, lapply(lapply(pe_list, unlist), "[",
                                                unique(unlist(c(sapply(pe_list,names)))))))

pe_table$model<-seq(1, nrow(pe_table), by =1)
pe_table

