# Model comparisons #
library(tidyverse)
library(here)
library(bbmle)
library(ggpubr)

library(MASS)

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

quant = function(col) quantile(col, c(0.025,0.975)) # 95% percentiles

##### Data ####
mort <- 0.13333
conversion <- exp(-15)
nresamp=1000

TIM_data <- loadRData(here("01 - Data", "TIM_data.RData")) # Data
#CPcomp_pe_mod <- loadRData(here("05 - Output", "Model Fits", "Colpidium", "colp_para_comp_fit.RData")) # if not running competition fits in this script load here

colp_para_wrk <- TIM_data %>% filter(Exp == "CP")

#### Model Fits ####


#### 0. Competition #####
##### Search Grid #####
comp.pre_sg <- expand.grid(r1_log = c(-2, -1, -0.5, 0.5), 
                           r2_log = c(-5, -3, -1, -0.5), 
                           alpha11 = c(0.0005, 0.001, 0.005, 0.001), 
                           alpha22 = c(0.0005, 0.001, 0.005, 0.001),  
                           alpha12 = c(0.0005, 0.001, 0.005, 0.001), 
                           alpha21 = c(0.0005, 0.001, 0.005, 0.001),  
                           sigma = c(0.2), 
                           sigma2 = c(0.2))


comp.pre_sg$error <- NA

comp.pre_sg$r1_log.est <- NA
comp.pre_sg$r2_log.est <- NA

comp.pre_sg$alpha11.est <- NA
comp.pre_sg$alpha12.est <- NA
comp.pre_sg$alpha22.est <- NA
comp.pre_sg$alpha21.est <- NA

comp.pre_sg$sigma.est <- NA
comp.pre_sg$sigma2.est <- NA

comp.pre_sg$converge <- NA
comp.pre_sg$BIC <- NA
comp.pre_sg$logLik <- NA
comp.pre_sg$strtpt <- seq(1:nrow(comp.pre_sg))

set.seed(67808)
comp.sim_grid <- comp.pre_sg[sample(comp.pre_sg$strtpt, 5*98, replace = F), ]


# Set the computer to run in Parallel #
bigstatsr::nb_cores()

my.cluster <- parallel::makeCluster(
  bigstatsr::nb_cores(),
  type = "FORK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()


# Run Grid #
colp_para_comp_flex_pe <- NA # name of file where parallel runs will be save

(t0<-Sys.time()) # timer
mcoptions <- list(preschedule=FALSE)

colp_para_comp_flex_pe <- foreach(i = 1:nrow(comp.sim_grid),
                                   .combine = 'rbind',
                                   .options.multicore = mcoptions,
                                   .packages = c('bbmle', 'odeintr', 'tidyverse', 'R.utils'),
                                   .errorhandling = "pass") %dopar%
  {
    compile_sys("FRS_lv", FRS.lv, pars = c("r1", "r2","alpha11", "alpha22", "alpha12", "alpha21"), method = "rk54")
    
    r<-tryCatch({ 
      setTimeLimit(35*60)
      fit = mle2(minuslogl = nll.odeint.lv,
                 browse_obj=F,
                 skip.hessian=F, 
                 trace = T,
                 control = list(trace=1, REPORT=1, maxit = 5000),
                 start = list(r1_log = comp.sim_grid$r1_log[[i]],
                              alpha11 = comp.sim_grid$alpha11[[i]],
                              alpha12 = comp.sim_grid$alpha12[[i]],
                              
                              r2_log = comp.sim_grid$r2_log[[i]],
                              alpha22 = comp.sim_grid$alpha22[[i]],
                              alpha21 = comp.sim_grid$alpha21[[i]],
                              
                              sigma = comp.sim_grid$sigma[[i]],
                              sigma2 = comp.sim_grid$sigma2[[i]]),
                 
                 data = with(list(N0 = Colp_T0.dens.ml,
                                  Ndead = Colp_T0.dens.ml - Colp_T24.dens.ml,
                                  M0 = Para_T0.dens.ml,
                                  M0.end = Para_T24.dens.ml,
                                  Tt = Time.day),
                             data = filter(colp_para_wrk, Type =="Comp_Control" | Type =="Colp_Con" | Type =="Para_Con"))
      ) 
      
      if(fit@details$convergence==0) {
        comp.sim_grid$r1_log.est[i] <- fit@coef[[1]]
        comp.sim_grid$alpha11.est[i] <- fit@coef[[2]]
        comp.sim_grid$alpha12.est[i] <- fit@coef[[3]]
        
        comp.sim_grid$r2_log.est[i] <- fit@coef[[4]]
        comp.sim_grid$alpha22.est[i] <- fit@coef[[5]]
        comp.sim_grid$alpha21.est[i] <- fit@coef[[6]]
        
        comp.sim_grid$sigma.est [i] <- fit@coef[[7]]
        comp.sim_grid$sigma2.est [i] <- fit@coef[[8]]
        
        comp.sim_grid$converge[i] <- fit@details[[4]]
        comp.sim_grid$BIC[i] <- BIC(fit)
        comp.sim_grid$logLik[i] <- logLik(fit)[[1]]
      } else {
        comp.sim_grid$r1_log.est[i] <- NaN
        comp.sim_grid$alpha11.est[i] <- NaN
        comp.sim_grid$alpha12.est[i] <- NaN
        
        comp.sim_grid$r2_log.est[i] <- NaN
        comp.sim_grid$alpha22.est[i] <- NaN
        comp.sim_grid$alpha21.est[i] <- NaN
        
        comp.sim_grid$sigma.est [i] <- NaN
        comp.sim_grid$sigma2.est [i] <- NaN
        
        comp.sim_grid$converge[i] <- NaN
        comp.sim_grid$iteration[i] <- NaN
        comp.sim_grid$BIC[i] <- NaN
        comp.sim_grid$logLik[i] <- NaN
      }
      
      #Remove fit for next iteration
      rm(fit)
    }, error=function(e) list(message = e$message)) # to print error in console use {cat("ERROR :", paste(conditionMessage(e), "Combo", (i)), "\n")}
    
    if(length(r) > 0){
      comp.sim_grid$error[i] <- r 
    } else {
      comp.sim_grid$error[i] <- NaN
    }  
    
    return(comp.sim_grid[i,])
    
  }

parallel::stopCluster(cl = my.cluster)
stopImplicitCluster()

colp_para_comp_flex_pe %>% arrange(BIC)
#save(colp_para_comp_flex_pe, file = "CP_comp_pe.RData")

###### Best Model #####
CP_LV_pe <- colp_para_comp_flex_pe
(CP_comp_PECs <- CP_LV_pe %>% 
   filter(converge == 0) %>% 
   filter(BIC <= min(CP_LV_pe$BIC, na.rm = T)+2)) #check for similarity in parameter estimates within 2 AIC of min

min_comp <- CP_comp_PECs %>% filter(BIC <= min(CP_comp_PECs$BIC, na.rm = T)+0.0001)

if(nrow(min_comp)>1){
  min_comp <- CP_comp_PECs %>% filter(strtpt %in% sample(filter(CP_comp_PECs, BIC <= min(CP_comp_PECs$BIC, na.rm = T)+0.0001)$strtpt,1))
}

CPcomp_pe_mod =  mle2(minuslogl = nll.odeint.lv,
                      method= "Nelder-Mead",
                      control = list(maxit = 5000),
                      
                      start = list(r1_log = min_comp$r1_log.est,
                                   alpha11 = min_comp$alpha11.est,
                                   alpha12 = min_comp$alpha12.est,
                                   
                                   r2_log = min_comp$r2_log.est,
                                   alpha22 = min_comp$alpha22.est,
                                   alpha21 = min_comp$alpha21.est,
                                   
                                   sigma = min_comp$sigma.est,
                                   sigma2 = min_comp$sigma2.est
                      ),
                      
                      data = with(list(N0 = Colp_T0.dens.ml,
                                       Ndead = Colp_T0.dens.ml - Colp_T24.dens.ml,
                                       M0 = Para_T0.dens.ml,
                                       M0.end = Para_T24.dens.ml,
                                       Tt = Time.day),
                                  data = filter(colp_para_wrk, Type =="Comp_Control" | Type =="Colp_Con" | Type =="Para_Con"))
                      
)

#save(CPcomp_pe_mod, file = "colp_para_comp_fit.RData")
comp_tidy<-tidy(CPcomp_pe_mod)

(nrow(CP_comp_PECs)/nrow(CP_LV_pe %>% filter(converge == 0))) * 100 # percent of starting combinations within 2 AIC that 


#### 1. Type II #####
##### Search Grid #####
t2.pre_sg <- expand.grid(a_log = seq(-3, -1, length.out = 10), 
                         h_log = seq(-5, -1, length.out = 10), 
                         sigma = c(0.2), 
                         sigma2 = c(0.2)
)

t2.pre_sg$error <- NA

t2.pre_sg$a_log.est <- NA
t2.pre_sg$h_log.est <- NA

t2.pre_sg$sigma.est <- NA
t2.pre_sg$sigma2.est <- NA

t2.pre_sg$converge <- NA
t2.pre_sg$BIC <- NA
t2.pre_sg$logLik <- NA
t2.pre_sg$strtpt <- seq(1:nrow(t2.pre_sg))

set.seed(39514)
t2.sim_grid <- t2.pre_sg[sample(t2.pre_sg$strtpt, 98, replace = F), ]


# Set the computer to run in Parallel #
bigstatsr::nb_cores()

my.cluster <- parallel::makeCluster(
  bigstatsr::nb_cores(),
  type = "FORK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()


# Run Grid #
CP_t2_fixLV_pe <- NA # name of file where parallel runs will be save
mcoptions <- list(preschedule=FALSE)

CP_t2_fixLV_pe <- foreach(i = 1:nrow(t2.sim_grid),
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
                              
                              sigma = t2.sim_grid$sigma[[i]],
                              sigma2 = t2.sim_grid$sigma2[[i]]
                 ),
                 
                 data = with(list(N0 = Colp_T0.dens.ml,
                                  Ndead = Colp_T0.dens.ml - Colp_T24.dens.ml,
                                  P = Pred_T0.dens.ml, 
                                  P.end = Pred_T24.dens.ml,
                                  M0 = Para_T0.dens.ml,
                                  M0.end = Para_T24.dens.ml,
                                  Tt = Time.day),
                             data = filter(colp_para_wrk)),
                 
                 fixed = list(m = mort,
                              c_log = -15,
                              a12 = 0,
                              h12 = 0,
                              
                              r1_log = coef(CPcomp_pe_mod)[1], 
                              alpha11 = coef(CPcomp_pe_mod)[2],
                              alpha12 = coef(CPcomp_pe_mod)[3],
                              r2_log = coef(CPcomp_pe_mod)[4], 
                              alpha22 = coef(CPcomp_pe_mod)[5],
                              alpha21 = coef(CPcomp_pe_mod)[6]) 
      )
      
      if(fit@details$convergence==0) {
        t2.sim_grid$a_log.est[i] <- fit@coef[[1]]
        t2.sim_grid$h_log.est[i] <- fit@coef[[2]]
        
        t2.sim_grid$sigma.est [i] <- fit@coef[[3]]
        t2.sim_grid$sigma2.est [i] <- fit@coef[[4]]
        
        t2.sim_grid$converge[i] <- fit@details$convergence
        t2.sim_grid$BIC[i] <- BIC(fit)
        t2.sim_grid$logLik[i] <- logLik(fit)[[1]]
      } else {
        t2.sim_grid$a_log.est[i] <- NaN
        
        t2.sim_grid$h_log.est[i] <- NaN
        
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

CP_t2_fixLV_pe %>% arrange(BIC)

#save(CP_t2_fixLV_pe, file = "CP_t2_fixLV_pe.RData")

###### Best Model #####
(CP_T2_PECs <- CP_t2_fixLV_pe %>% 
   filter(converge == 0) %>% 
   filter(BIC <= min(CP_t2_fixLV_pe$BIC, na.rm = T)+2)) #check for similarity in parameter estimates within 2 AIC of min

min_T2 <- CP_T2_PECs %>% filter(BIC <= min(CP_T2_PECs$BIC, na.rm = T)+0.0001) %>% dplyr::select(c(a_log.est:sigma2.est, BIC))

if(nrow(min_T2)>1){
  min_T2 <- CP_T2_PECs %>% filter(strtpt %in% sample(filter(CP_T2_PECs, BIC <= min(CP_T2_PECs$BIC, na.rm = T)+0.0001)$strtpt,1))
}

CP_min_T2_fixLV_mod =  mle2(minuslogl = nll.odeint.t2.a1.h1.lv,
                            method= "Nelder-Mead",
                            control = list(maxit = 5000),
                            
                            start = list(a_log = min_T2$a_log.est, 
                                         h_log = min_T2$h_log.est,
                                         
                                         sigma = min_T2$sigma.est,
                                         sigma2 = min_T2$sigma2.est
                            ),
                            
                            data = with(list(N0 = Colp_T0.dens.ml,
                                             Ndead = Colp_T0.dens.ml - Colp_T24.dens.ml,
                                             P = Pred_T0.dens.ml, 
                                             P.end = Pred_T24.dens.ml,
                                             M0 = Para_T0.dens.ml,
                                             M0.end = Para_T24.dens.ml,
                                             Tt = Time.day),
                                        data = filter(colp_para_wrk)),
                            
                            fixed = list(m = mort,
                                         c_log = -15,
                                         a12 = 0,
                                         h12 = 0,
                                         r1_log = coef(CPcomp_pe_mod)[1], 
                                         alpha11 = coef(CPcomp_pe_mod)[2],
                                         alpha12 = coef(CPcomp_pe_mod)[3],
                                         r2_log = coef(CPcomp_pe_mod)[4], 
                                         alpha22 = coef(CPcomp_pe_mod)[5],
                                         alpha21 = coef(CPcomp_pe_mod)[6])
)
#save(CP_min_T2_fixLV_mod, file = "CP_T2_fixLV_mod.Rdata")
(Colp_fit1 <- tidy(CP_min_T2_fixLV_mod))

(nrow(CP_T2_PECs)/nrow(CP_t2_fixLV_pe %>% filter(converge == 0))) * 100 # percent of starting combinations within 2 AIC that 
BIC(CP_min_T2_fixLV_mod)

set.seed(1001)
pars.picked = mvrnorm(1000, mu = Colp_fit1$estimate, Sigma = vcov(CP_min_T2_fixLV_mod)) # pick new parameter values by sampling from multivariate normal distribution based on fit
Colp_fit1_pe_confint <- apply(pars.picked,2,quant)

#### 2. Type II + linear a TIM #####
##### Search Grid #####
t2.a12.pre_sg <- expand.grid(a_log = seq(-3, -1, length.out = 5), 
                             a12 = seq(-0.0005, 0.0001, length.out = 5),
                             h_log = seq(-5, -1, length.out = 5), 
                             sigma = c(0.2), 
                             sigma2 = c(0.2)
)

t2.a12.pre_sg$error <- NA

t2.a12.pre_sg$a_log.est <- NA
t2.a12.pre_sg$a12.est <- NA
t2.a12.pre_sg$h_log.est <- NA

t2.a12.pre_sg$sigma.est <- NA
t2.a12.pre_sg$sigma2.est <- NA

t2.a12.pre_sg$converge <- NA
t2.a12.pre_sg$BIC <- NA
t2.a12.pre_sg$logLik <- NA
t2.a12.pre_sg$strtpt <- seq(1:nrow(t2.a12.pre_sg))

set.seed(54791)
t2.a12.sim_grid <- t2.a12.pre_sg[sample(t2.a12.pre_sg$strtpt, 98, replace = F), ]


# Set the computer to run in Parallel #
bigstatsr::nb_cores()

my.cluster <- parallel::makeCluster(
  bigstatsr::nb_cores(),
  type = "FORK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()


# Run Grid #
CP_t2.a12_fixLV_pe <- NA # name of file where parallel runs will be save
mcoptions <- list(preschedule=FALSE)

CP_t2.a12_fixLV_pe <- foreach(i = 1:nrow(t2.a12.sim_grid),
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
                              
                              sigma = t2.a12.sim_grid$sigma[[i]],
                              sigma2 = t2.a12.sim_grid$sigma2[[i]]
                 ),
                 
                 data = with(list(N0 = Colp_T0.dens.ml,
                                  Ndead = Colp_T0.dens.ml - Colp_T24.dens.ml,
                                  P = Pred_T0.dens.ml, 
                                  P.end = Pred_T24.dens.ml,
                                  M0 = Para_T0.dens.ml,
                                  M0.end = Para_T24.dens.ml,
                                  Tt = Time.day),
                             data = filter(colp_para_wrk)),
                 
                 fixed = list(m = mort,
                              c_log = -15,
                              h12 = 0,
                              
                              r1_log = coef(CPcomp_pe_mod)[1], 
                              alpha11 = coef(CPcomp_pe_mod)[2],
                              alpha12 = coef(CPcomp_pe_mod)[3],
                              r2_log = coef(CPcomp_pe_mod)[4], 
                              alpha22 = coef(CPcomp_pe_mod)[5],
                              alpha21 = coef(CPcomp_pe_mod)[6]) 
      )
      
      if(fit@details$convergence==0) {
        t2.a12.sim_grid$a_log.est[i] <- fit@coef[[1]]
        t2.a12.sim_grid$a12.est[i] <- fit@coef[[2]]
        t2.a12.sim_grid$h_log.est[i] <- fit@coef[[3]]
        
        t2.a12.sim_grid$sigma.est [i] <- fit@coef[[4]]
        t2.a12.sim_grid$sigma2.est [i] <- fit@coef[[5]]
        
        t2.a12.sim_grid$converge[i] <- fit@details$convergence
        t2.a12.sim_grid$BIC[i] <- BIC(fit)
        t2.a12.sim_grid$logLik[i] <- logLik(fit)[[1]]
      } else {
        t2.a12.sim_grid$a_log.est[i] <- NaN
        t2.a12.sim_grid$a12.est[i] <- NaN
        t2.a12.sim_grid$h_log.est[i] <- NaN
        
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

CP_t2.a12_fixLV_pe %>% arrange(BIC)
#save(CP_t2.a12_fixLV_pe, file = "CP_t2.a12_fixLV_pe.RData")

###### Best Model #####
(CP_t2.a12_PECs <- CP_t2.a12_fixLV_pe %>% 
   filter(converge == 0) %>% 
   filter(BIC <= min(CP_t2.a12_fixLV_pe$BIC, na.rm = T)+2)) #check for similarity in parameter estimates within 2 AIC of min

min_t2.a12 <- CP_t2.a12_PECs %>% filter(BIC <= min(CP_t2.a12_PECs$BIC, na.rm = T)+0.0001) %>% dplyr::select(c(a_log.est:sigma2.est, BIC))

if(nrow(min_t2.a12)>1){
  min_t2.a12 <- CP_t2.a12_PECs %>% filter(strtpt %in% sample(filter(CP_t2.a12_PECs, BIC <= min(CP_t2.a12_PECs$BIC, na.rm = T)+0.0001)$strtpt,1))
}

CP_min_t2.a12_fixLV_mod =  mle2(minuslogl = nll.odeint.t2.a1.h1.lv,
                                method= "Nelder-Mead",
                                control = list(maxit = 5000),
                                
                                start = list(a_log = min_t2.a12$a_log.est, 
                                             h_log = min_t2.a12$h_log.est,
                                             a12 = min_t2.a12$a12.est,
                                             
                                             sigma = min_t2.a12$sigma.est,
                                             sigma2 = min_t2.a12$sigma2.est
                                ),
                                
                                data = with(list(N0 = Colp_T0.dens.ml,
                                                 Ndead = Colp_T0.dens.ml - Colp_T24.dens.ml,
                                                 P = Pred_T0.dens.ml, 
                                                 P.end = Pred_T24.dens.ml,
                                                 M0 = Para_T0.dens.ml,
                                                 M0.end = Para_T24.dens.ml,
                                                 Tt = Time.day),
                                            data = filter(colp_para_wrk)),
                                
                                fixed = list(m = mort,
                                             c_log = -15,
                                             h12 = 0, 
                                             r1_log = coef(CPcomp_pe_mod)[1], 
                                             alpha11 = coef(CPcomp_pe_mod)[2],
                                             alpha12 = coef(CPcomp_pe_mod)[3],
                                             r2_log = coef(CPcomp_pe_mod)[4], 
                                             alpha22 = coef(CPcomp_pe_mod)[5],
                                             alpha21 = coef(CPcomp_pe_mod)[6])
)
#save(CP_min_t2.a12_fixLV_mod, file = "CP_t2.a12_fixLV_mod.Rdata")
Colp_fit2 <- tidy(CP_min_t2.a12_fixLV_mod)

(nrow(CP_t2.a12_PECs)/nrow(CP_t2.a12_fixLV_pe %>% filter(converge == 0))) * 100 # percent of starting combinations within 2 AIC that 
BIC(CP_min_t2.a12_fixLV_mod)

set.seed(1001)
pars.picked = mvrnorm(1000, mu = Colp_fit2$estimate, Sigma = vcov(CP_min_t2.a12_fixLV_mod)) # pick new parameter values by sampling from multivariate normal distribution based on fit
Colp_fit2_pe_confint <- apply(pars.picked,2,quant)

#### 3. Type II + linear h TIM #####
##### Search Grid #####
t2.h12.pre_sg <- expand.grid(a_log = seq(-3, -1, length.out = 5), 
                             h12 = seq(-0.0005, 0.0001, length.out = 5),
                             h_log = seq(-5, -1, length.out = 5), 
                             sigma = c(0.2), 
                             sigma2 = c(0.2)
)

t2.h12.pre_sg$error <- NA

t2.h12.pre_sg$a_log.est <- NA
t2.h12.pre_sg$h_log.est <- NA
t2.h12.pre_sg$h12.est <- NA

t2.h12.pre_sg$sigma.est <- NA
t2.h12.pre_sg$sigma2.est <- NA

t2.h12.pre_sg$converge <- NA
t2.h12.pre_sg$BIC <- NA
t2.h12.pre_sg$logLik <- NA
t2.h12.pre_sg$strtpt <- seq(1:nrow(t2.h12.pre_sg))

set.seed(64715)
t2.h12.sim_grid <- t2.h12.pre_sg[sample(t2.h12.pre_sg$strtpt, 98, replace = F), ]


# Set the computer to run in Parallel #
bigstatsr::nb_cores()

my.cluster <- parallel::makeCluster(
  bigstatsr::nb_cores(),
  type = "FORK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()


# Run Grid #
CP_t2.h12_fixLV_pe <- NA # name of file where parallel runs will be save
mcoptions <- list(preschedule=FALSE)

CP_t2.h12_fixLV_pe <- foreach(i = 1:nrow(t2.h12.sim_grid),
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
                              
                              sigma = t2.h12.sim_grid$sigma[[i]],
                              sigma2 = t2.h12.sim_grid$sigma2[[i]]
                 ),
                 
                 data = with(list(N0 = Colp_T0.dens.ml,
                                  Ndead = Colp_T0.dens.ml - Colp_T24.dens.ml,
                                  P = Pred_T0.dens.ml, 
                                  P.end = Pred_T24.dens.ml,
                                  M0 = Para_T0.dens.ml,
                                  M0.end = Para_T24.dens.ml,
                                  Tt = Time.day),
                             data = filter(colp_para_wrk)),
                 
                 fixed = list(m = mort,
                              c_log = -15,
                              a12 = 0,
                              
                              r1_log = coef(CPcomp_pe_mod)[1], 
                              alpha11 = coef(CPcomp_pe_mod)[2],
                              alpha12 = coef(CPcomp_pe_mod)[3],
                              r2_log = coef(CPcomp_pe_mod)[4], 
                              alpha22 = coef(CPcomp_pe_mod)[5],
                              alpha21 = coef(CPcomp_pe_mod)[6]) 
      )
      
      if(fit@details$convergence==0) {
        t2.h12.sim_grid$a_log.est[i] <- fit@coef[[1]]
        t2.h12.sim_grid$h_log.est[i] <- fit@coef[[2]]
        t2.h12.sim_grid$h12.est[i] <- fit@coef[[3]]
        
        t2.h12.sim_grid$sigma.est [i] <- fit@coef[[4]]
        t2.h12.sim_grid$sigma2.est [i] <- fit@coef[[5]]
        
        t2.h12.sim_grid$converge[i] <- fit@details$convergence
        t2.h12.sim_grid$BIC[i] <- BIC(fit)
        t2.h12.sim_grid$logLik[i] <- logLik(fit)[[1]]
      } else {
        t2.h12.sim_grid$a_log.est[i] <- NaN
        t2.h12.sim_grid$h_log.est[i] <- NaN
        t2.h12.sim_grid$h12.est[i] <- NaN
        
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

CP_t2.h12_fixLV_pe %>% arrange(BIC)

#save(CP_t2.h12_fixLV_pe, file = "CP_t2.h12_fixLV_pe.RData")

###### Best Model #####
(CP_t2.h12_PECs <- CP_t2.h12_fixLV_pe %>% 
   filter(converge == 0) %>% 
   filter(BIC <= min(CP_t2.h12_fixLV_pe$BIC, na.rm = T)+2)) #check for similarity in parameter estimates within 2 AIC of min

min_t2.h12 <- CP_t2.h12_PECs %>% filter(BIC <= min(CP_t2.h12_PECs$BIC, na.rm = T)+0.0001) %>% dplyr::select(c(a_log.est:sigma2.est, BIC))

if(nrow(min_t2.h12)>1){
  min_t2.h12 <- CP_t2.h12_PECs %>% filter(strtpt %in% sample(filter(CP_t2.h12_PECs, BIC <= min(CP_t2.h12_PECs$BIC, na.rm = T)+0.0001)$strtpt,1))
}

CP_min_t2.h12_fixLV_mod =  mle2(minuslogl = nll.odeint.t2.a1.h1.lv,
                                method= "Nelder-Mead",
                                control = list(maxit = 5000),
                                
                                start = list(a_log = min_t2.h12$a_log.est, 
                                             h_log = min_t2.h12$h_log.est,
                                             h12 = min_t2.h12$h12.est,
                                             sigma = min_t2.h12$sigma.est,
                                             sigma2 = min_t2.h12$sigma2.est
                                ),
                                
                                data = with(list(N0 = Colp_T0.dens.ml,
                                                 Ndead = Colp_T0.dens.ml - Colp_T24.dens.ml,
                                                 P = Pred_T0.dens.ml, 
                                                 P.end = Pred_T24.dens.ml,
                                                 M0 = Para_T0.dens.ml,
                                                 M0.end = Para_T24.dens.ml,
                                                 Tt = Time.day),
                                            data = filter(colp_para_wrk)),
                                
                                fixed = list(m = mort,
                                             c_log = -15,
                                             a12=0,
                                             r1_log = coef(CPcomp_pe_mod)[1], 
                                             alpha11 = coef(CPcomp_pe_mod)[2],
                                             alpha12 = coef(CPcomp_pe_mod)[3],
                                             r2_log = coef(CPcomp_pe_mod)[4], 
                                             alpha22 = coef(CPcomp_pe_mod)[5],
                                             alpha21 = coef(CPcomp_pe_mod)[6])
)
#save(CP_min_t2.h12_fixLV_mod, file = "CP_t2.h12_fixLV_mod.Rdata")
Colp_fit3 <- tidy(CP_min_t2.h12_fixLV_mod)

(nrow(CP_t2.h12_PECs)/nrow(CP_t2.h12_fixLV_pe %>% filter(converge == 0))) * 100 # percent of starting combinations within 2 AIC that 
BIC(CP_min_t2.h12_fixLV_mod)

set.seed(1001)
pars.picked = mvrnorm(1000, mu = Colp_fit3$estimate, Sigma = vcov(CP_min_t2.h12_fixLV_mod)) # pick new parameter values by sampling from multivariate normal distribution based on fit
Colp_fit3_pe_confint <- apply(pars.picked,2,quant)

#### 4. Type II + linear a TIM  + linear h TIM #####
##### Search Grid #####
t2.a12.h12.pre_sg <- expand.grid(a_log = seq(-3, -1, length.out = 5), 
                                 a12 = seq(-0.0005, 0.0001, length.out = 5),
                                 h12 = seq(-0.0005, 0.0001, length.out = 5),
                                 h_log = seq(-5, -1, length.out = 5), 
                                 sigma = c(0.2), 
                                 sigma2 = c(0.2)
)

t2.a12.h12.pre_sg$error <- NA

t2.a12.h12.pre_sg$a_log.est <- NA
t2.a12.h12.pre_sg$a12.est <- NA
t2.a12.h12.pre_sg$h_log.est <- NA
t2.a12.h12.pre_sg$h12.est <- NA

t2.a12.h12.pre_sg$sigma.est <- NA
t2.a12.h12.pre_sg$sigma2.est <- NA

t2.a12.h12.pre_sg$converge <- NA
t2.a12.h12.pre_sg$BIC <- NA
t2.a12.h12.pre_sg$logLik <- NA
t2.a12.h12.pre_sg$strtpt <- seq(1:nrow(t2.a12.h12.pre_sg))

set.seed(47815)
t2.a12.h12.sim_grid <- t2.a12.h12.pre_sg[sample(t2.a12.h12.pre_sg$strtpt, 98, replace = F), ]


# Set the computer to run in Parallel #
bigstatsr::nb_cores()

my.cluster <- parallel::makeCluster(
  bigstatsr::nb_cores(),
  type = "FORK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()


# Run Grid #
CP_t2.a12.h12_fixLV_pe <- NA # name of file where parallel runs will be save
mcoptions <- list(preschedule=FALSE)

CP_t2.a12.h12_fixLV_pe <- foreach(i = 1:nrow(t2.a12.h12.sim_grid),
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
                              
                              sigma = t2.a12.h12.sim_grid$sigma[[i]],
                              sigma2 = t2.a12.h12.sim_grid$sigma2[[i]]
                 ),
                 
                 data = with(list(N0 = Colp_T0.dens.ml,
                                  Ndead = Colp_T0.dens.ml - Colp_T24.dens.ml,
                                  P = Pred_T0.dens.ml, 
                                  P.end = Pred_T24.dens.ml,
                                  M0 = Para_T0.dens.ml,
                                  M0.end = Para_T24.dens.ml,
                                  Tt = Time.day),
                             data = filter(colp_para_wrk)),
                 
                 fixed = list(m = mort,
                              c_log = -15,
                              
                              r1_log = coef(CPcomp_pe_mod)[1], 
                              alpha11 = coef(CPcomp_pe_mod)[2],
                              alpha12 = coef(CPcomp_pe_mod)[3],
                              r2_log = coef(CPcomp_pe_mod)[4], 
                              alpha22 = coef(CPcomp_pe_mod)[5],
                              alpha21 = coef(CPcomp_pe_mod)[6]) 
      )
      
      if(fit@details$convergence==0) {
        t2.a12.h12.sim_grid$a_log.est[i] <- fit@coef[[1]]
        t2.a12.h12.sim_grid$a12.est[i] <- fit@coef[[2]]
        t2.a12.h12.sim_grid$h_log.est[i] <- fit@coef[[3]]
        t2.a12.h12.sim_grid$h12.est[i] <- fit@coef[[4]]
        
        t2.a12.h12.sim_grid$sigma.est [i] <- fit@coef[[5]]
        t2.a12.h12.sim_grid$sigma2.est [i] <- fit@coef[[6]]
        
        t2.a12.h12.sim_grid$converge[i] <- fit@details$convergence
        t2.a12.h12.sim_grid$BIC[i] <- BIC(fit)
        t2.a12.h12.sim_grid$logLik[i] <- logLik(fit)[[1]]
      } else {
        t2.a12.h12.sim_grid$a_log.est[i] <- NaN
        t2.a12.h12.sim_grid$a12.est[i] <- NaN
        t2.a12.h12.sim_grid$h_log.est[i] <- NaN
        t2.a12.h12.sim_grid$h12.est[i] <- NaN
        
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

CP_t2.a12.h12_fixLV_pe %>% arrange(BIC)

#save(CP_t2.a12.h12_fixLV_pe, file = "CP_t2.a12.h12_fixLV_pe.RData")

###### Best Model #####
(CP_t2.a12.h12_PECs <- CP_t2.a12.h12_fixLV_pe %>% 
   filter(converge == 0) %>% 
   filter(BIC <= min(CP_t2.a12.h12_fixLV_pe$BIC, na.rm = T)+2)) #check for similarity in parameter estimates within 2 AIC of min

min_t2.a12.h12 <- CP_t2.a12.h12_PECs %>% filter(BIC <= min(CP_t2.a12.h12_PECs$BIC, na.rm = T)+0.0001) %>% dplyr::select(c(a_log.est:sigma2.est, BIC))

if(nrow(min_t2.a12.h12)>1){
  min_t2.a12.h12 <- CP_t2.a12.h12_PECs %>% filter(strtpt %in% sample(filter(CP_t2.a12.h12_PECs, BIC <= min(CP_t2.a12.h12_PECs$BIC, na.rm = T)+0.0001)$strtpt,1))
}

CP_min_t2.a12.h12_fixLV_mod =  mle2(minuslogl = nll.odeint.t2.a1.h1.lv,
                                    method= "Nelder-Mead",
                                    control = list(maxit = 5000),
                                    
                                    start = list(a_log = min_t2.a12.h12$a_log.est, 
                                                 a12 = min_t2.a12.h12$a12.est,
                                                 h_log = min_t2.a12.h12$h_log.est,
                                                 h12 = min_t2.a12.h12$h12.est,
                                                 sigma = min_t2.a12.h12$sigma.est,
                                                 sigma2 = min_t2.a12.h12$sigma2.est
                                    ),
                                    
                                    data = with(list(N0 = Colp_T0.dens.ml,
                                                     Ndead = Colp_T0.dens.ml - Colp_T24.dens.ml,
                                                     P = Pred_T0.dens.ml, 
                                                     P.end = Pred_T24.dens.ml,
                                                     M0 = Para_T0.dens.ml,
                                                     M0.end = Para_T24.dens.ml,
                                                     Tt = Time.day),
                                                data = filter(colp_para_wrk)),
                                    
                                    fixed = list(m = mort,
                                                 c_log = -15,
                                                 r1_log = coef(CPcomp_pe_mod)[1], 
                                                 alpha11 = coef(CPcomp_pe_mod)[2],
                                                 alpha12 = coef(CPcomp_pe_mod)[3],
                                                 r2_log = coef(CPcomp_pe_mod)[4], 
                                                 alpha22 = coef(CPcomp_pe_mod)[5],
                                                 alpha21 = coef(CPcomp_pe_mod)[6])
)
#save(CP_min_t2.a12.h12_fixLV_mod, file = "CP_t2.a12.h12_fixLV_mod.Rdata")
Colp_fit4 <- tidy(CP_min_t2.a12.h12_fixLV_mod)

(nrow(CP_t2.a12.h12_PECs)/nrow(CP_t2.a12.h12_fixLV_pe %>% filter(converge == 0))) * 100 # percent of starting combinations within 2 AIC that 
BIC(CP_min_t2.a12.h12_fixLV_mod)

set.seed(1001)
pars.picked = mvrnorm(1000, mu = Colp_fit4$estimate, Sigma = vcov(CP_min_t2.a12.h12_fixLV_mod)) # pick new parameter values by sampling from multivariate normal distribution based on fit
Colp_fit4_pe_confint <- apply(pars.picked,2,quant)


#### 5. Type II + van Veen #####
##### Search Grid #####
t2.vv.pre_sg <- expand.grid(a_log = seq(-3, -1, length.out = 5), 
                            h_log = seq(-5, -1, length.out = 5), 
                            log_w = seq(-10, -5, length.out = 5), 
                            sigma = c(0.2), 
                            sigma2 = c(0.2)
)

t2.vv.pre_sg$error <- NA

t2.vv.pre_sg$a_log.est <- NA
t2.vv.pre_sg$h_log.est <- NA
t2.vv.pre_sg$log_w.est <- NA

t2.vv.pre_sg$sigma.est <- NA
t2.vv.pre_sg$sigma2.est <- NA

t2.vv.pre_sg$converge <- NA
t2.vv.pre_sg$BIC <- NA
t2.vv.pre_sg$logLik <- NA
t2.vv.pre_sg$strtpt <- seq(1:nrow(t2.vv.pre_sg))

set.seed(13875)
t2.vv.sim_grid <- t2.vv.pre_sg[sample(t2.vv.pre_sg$strtpt, 98, replace = F), ]


# Set the computer to run in Parallel #
bigstatsr::nb_cores()

my.cluster <- parallel::makeCluster(
  bigstatsr::nb_cores(),
  type = "FORK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()


# Run Grid #
CP_t2.vv_fixLV_pe <- NA # name of file where parallel runs will be save
mcoptions <- list(preschedule=FALSE)

CP_t2.vv_fixLV_pe <- foreach(i = 1:nrow(t2.vv.sim_grid),
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
                              
                              sigma = t2.vv.sim_grid$sigma[[i]],
                              sigma2 = t2.vv.sim_grid$sigma2[[i]]
                 ),
                 
                 data = with(list(N0 = Colp_T0.dens.ml,
                                  Ndead = Colp_T0.dens.ml - Colp_T24.dens.ml,
                                  P = Pred_T0.dens.ml, 
                                  P.end = Pred_T24.dens.ml,
                                  M0 = Para_T0.dens.ml,
                                  M0.end = Para_T24.dens.ml,
                                  Tt = Time.day),
                             data = filter(colp_para_wrk)),
                 
                 fixed = list(m = mort,
                              c_log = -15,
                              h12 = 0,
                              
                              r1_log = coef(CPcomp_pe_mod)[1], 
                              alpha11 = coef(CPcomp_pe_mod)[2],
                              alpha12 = coef(CPcomp_pe_mod)[3],
                              r2_log = coef(CPcomp_pe_mod)[4], 
                              alpha22 = coef(CPcomp_pe_mod)[5],
                              alpha21 = coef(CPcomp_pe_mod)[6]) 
      )
      
      if(fit@details$convergence==0) {
        t2.vv.sim_grid$a_log.est[i] <- fit@coef[[1]]
        t2.vv.sim_grid$h_log.est[i] <- fit@coef[[2]]
        t2.vv.sim_grid$log_w.est[i] <- fit@coef[[3]]
        
        t2.vv.sim_grid$sigma.est [i] <- fit@coef[[4]]
        t2.vv.sim_grid$sigma2.est [i] <- fit@coef[[5]]
        
        t2.vv.sim_grid$converge[i] <- fit@details$convergence
        t2.vv.sim_grid$BIC[i] <- BIC(fit)
        t2.vv.sim_grid$logLik[i] <- logLik(fit)[[1]]
      } else {
        t2.vv.sim_grid$a_log.est[i] <- NaN
        t2.vv.sim_grid$h_log.est[i] <- NaN
        t2.vv.sim_grid$log_w.est[i] <- NaN
        
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

CP_t2.vv_fixLV_pe %>% arrange(BIC)

#save(CP_t2.vv_fixLV_pe, file = "CP_t2.vv_fixLV_pe.RData")

###### Best Model #####
(CP_t2.vv_PECs <- CP_t2.vv_fixLV_pe %>% 
   filter(converge == 0) %>% 
   filter(BIC <= min(CP_t2.vv_fixLV_pe$BIC, na.rm = T)+2)) #check for similarity in parameter estimates within 2 AIC of min

min_t2.vv <- CP_t2.vv_PECs %>% filter(strtpt %in% sample(filter(CP_t2.vv_PECs, BIC <= min(CP_t2.vv_PECs$BIC, na.rm = T)+0.0001)$strtpt,1))

if(nrow(min_t2.vv)>1){
  min_t2.vv <- CP_t2.vv_PECs %>% filter(strtpt %in% sample(filter(CP_t2.vv_PECs, BIC <= min(CP_t2.vv_PECs$BIC, na.rm = T)+0.0001)$strtpt,1))
}

CP_min_t2.vv_fixLV_mod =  mle2(minuslogl = nll.odeint.t2.h1.VV.lv,
                               method= "Nelder-Mead",
                               control = list(maxit = 5000),
                               
                               start = list(a_log = min_t2.vv$a_log.est, 
                                            h_log = min_t2.vv$h_log.est,
                                            log_w = min_t2.vv$log_w.est,
                                            sigma = min_t2.vv$sigma.est,
                                            sigma2 = min_t2.vv$sigma2.est
                               ),
                               
                               data = with(list(N0 = Colp_T0.dens.ml,
                                                Ndead = Colp_T0.dens.ml - Colp_T24.dens.ml,
                                                P = Pred_T0.dens.ml, 
                                                P.end = Pred_T24.dens.ml,
                                                M0 = Para_T0.dens.ml,
                                                M0.end = Para_T24.dens.ml,
                                                Tt = Time.day),
                                           data = filter(colp_para_wrk)),
                               
                               fixed = list(m = mort,
                                            c_log = -15,
                                            h12=0,
                                            
                                            r1_log = coef(CPcomp_pe_mod)[1], 
                                            alpha11 = coef(CPcomp_pe_mod)[2],
                                            alpha12 = coef(CPcomp_pe_mod)[3],
                                            r2_log = coef(CPcomp_pe_mod)[4], 
                                            alpha22 = coef(CPcomp_pe_mod)[5],
                                            alpha21 = coef(CPcomp_pe_mod)[6])
)
#save(CP_min_t2.vv_fixLV_mod, file = "CP_t2.vv_fixLV_mod.Rdata")
Colp_fit5 <- tidy(CP_min_t2.vv_fixLV_mod)

(nrow(CP_t2.vv_PECs)/nrow(CP_t2.vv_fixLV_pe %>% filter(converge == 0))) * 100 # percent of starting combinations within 2 AIC that 
BIC(CP_min_t2.vv_fixLV_mod)

set.seed(1001)
pars.picked = mvrnorm(1000, mu = Colp_fit5$estimate, Sigma = vcov(CP_min_t2.vv_fixLV_mod)) # pick new parameter values by sampling from multivariate normal distribution based on fit
Colp_fit5_pe_confint <- apply(pars.picked,2,quant)

#### 6. Type II + van Veen + linear h TIM #####
##### Search Grid #####
t2.vv.h12.pre_sg <- expand.grid(a_log = seq(-3, -1, length.out = 5), 
                                h_log = seq(-5, -1, length.out = 5), 
                                h12 = seq(-0.0005, 0.0001, length.out = 5),
                                log_w = seq(-10, -5, length.out = 5),
                                
                                sigma = c(0.2), 
                                sigma2 = c(0.2)
)

t2.vv.h12.pre_sg$error <- NA

t2.vv.h12.pre_sg$a_log.est <- NA
t2.vv.h12.pre_sg$h_log.est <- NA
t2.vv.h12.pre_sg$h12.est <- NA
t2.vv.h12.pre_sg$log_w.est <- NA

t2.vv.h12.pre_sg$sigma.est <- NA
t2.vv.h12.pre_sg$sigma2.est <- NA

t2.vv.h12.pre_sg$converge <- NA
t2.vv.h12.pre_sg$BIC <- NA
t2.vv.h12.pre_sg$logLik <- NA
t2.vv.h12.pre_sg$strtpt <- seq(1:nrow(t2.vv.h12.pre_sg))

set.seed(97464)
t2.vv.h12.sim_grid <- t2.vv.h12.pre_sg[sample(t2.vv.h12.pre_sg$strtpt, 98, replace = F), ]


# Set the computer to run in Parallel #
bigstatsr::nb_cores()

my.cluster <- parallel::makeCluster(
  bigstatsr::nb_cores(),
  type = "FORK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()


# Run Grid #
CP_t2.vv.h12_fixLV_pe <- NA # name of file where parallel runs will be save
mcoptions <- list(preschedule=FALSE)

CP_t2.vv.h12_fixLV_pe <- foreach(i = 1:nrow(t2.vv.h12.sim_grid),
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
                              
                              sigma = t2.vv.h12.sim_grid$sigma[[i]],
                              sigma2 = t2.vv.h12.sim_grid$sigma2[[i]]
                 ),
                 
                 data = with(list(N0 = Colp_T0.dens.ml,
                                  Ndead = Colp_T0.dens.ml - Colp_T24.dens.ml,
                                  P = Pred_T0.dens.ml, 
                                  P.end = Pred_T24.dens.ml,
                                  M0 = Para_T0.dens.ml,
                                  M0.end = Para_T24.dens.ml,
                                  Tt = Time.day),
                             data = filter(colp_para_wrk)),
                 
                 fixed = list(m = mort,
                              c_log = -15,
                              
                              r1_log = coef(CPcomp_pe_mod)[1], 
                              alpha11 = coef(CPcomp_pe_mod)[2],
                              alpha12 = coef(CPcomp_pe_mod)[3],
                              r2_log = coef(CPcomp_pe_mod)[4], 
                              alpha22 = coef(CPcomp_pe_mod)[5],
                              alpha21 = coef(CPcomp_pe_mod)[6]) 
      )
      
      if(fit@details$convergence==0) {
        t2.vv.h12.sim_grid$a_log.est[i] <- fit@coef[[1]]
        t2.vv.h12.sim_grid$h_log.est[i] <- fit@coef[[2]]
        t2.vv.h12.sim_grid$h12.est[i] <- fit@coef[[3]]
        t2.vv.h12.sim_grid$log_w.est[i] <- fit@coef[[4]]
        
        t2.vv.h12.sim_grid$sigma.est [i] <- fit@coef[[5]]
        t2.vv.h12.sim_grid$sigma2.est [i] <- fit@coef[[6]]
        
        t2.vv.h12.sim_grid$converge[i] <- fit@details$convergence
        t2.vv.h12.sim_grid$BIC[i] <- BIC(fit)
        t2.vv.h12.sim_grid$logLik[i] <- logLik(fit)[[1]]
      } else {
        t2.vv.h12.sim_grid$a_log.est[i] <- NaN
        t2.vv.h12.sim_grid$h_log.est[i] <- NaN
        t2.vv.h12.sim_grid$h12.est[i] <- NaN
        t2.vv.h12.sim_grid$log_w.est[i] <- NaN
        
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

CP_t2.vv.h12_fixLV_pe %>% arrange(BIC)

#save(CP_t2.vv.h12_fixLV_pe, file = "CP_t2.vv.h12_fixLV_pe.RData")

###### Best Model #####
(CP_t2.vv.h12_PECs <- CP_t2.vv.h12_fixLV_pe %>% 
   filter(converge == 0) %>% 
   filter(BIC <= min(CP_t2.vv.h12_fixLV_pe$BIC, na.rm = T)+2)) #check for similarity in parameter estimates within 2 AIC of min

min_t2.vv.h12 <- CP_t2.vv.h12_PECs %>% filter(BIC <= min(CP_t2.vv.h12_PECs$BIC, na.rm = T)+0.0001) %>% dplyr::select(c(a_log.est:sigma2.est, BIC))

if(nrow(min_t2.vv.h12)>1){
  min_t2.vv.h12 <- CP_t2.vv.h12_PECs %>% filter(strtpt %in% sample(filter(CP_t2.vv.h12_PECs, BIC <= min(CP_t2.vv.h12_PECs$BIC, na.rm = T)+0.0001)$strtpt,1))
}

CP_min_t2.vv.h12_fixLV_mod =  mle2(minuslogl = nll.odeint.t2.h1.VV.lv,
                                   method= "Nelder-Mead",
                                   control = list(maxit = 5000),
                                   
                                   start = list(a_log = min_t2.vv.h12$a_log.est, 
                                                h_log = min_t2.vv.h12$h_log.est,
                                                h12 = min_t2.vv.h12$h12.est,
                                                log_w = min_t2.vv.h12$log_w.est,
                                                sigma = min_t2.vv.h12$sigma.est,
                                                sigma2 = min_t2.vv.h12$sigma2.est
                                   ),
                                   
                                   data = with(list(N0 = Colp_T0.dens.ml,
                                                    Ndead = Colp_T0.dens.ml - Colp_T24.dens.ml,
                                                    P = Pred_T0.dens.ml, 
                                                    P.end = Pred_T24.dens.ml,
                                                    M0 = Para_T0.dens.ml,
                                                    M0.end = Para_T24.dens.ml,
                                                    Tt = Time.day),
                                               data = filter(colp_para_wrk)),
                                   
                                   fixed = list(m = mort,
                                                c_log = -15,
                                                
                                                r1_log = coef(CPcomp_pe_mod)[1], 
                                                alpha11 = coef(CPcomp_pe_mod)[2],
                                                alpha12 = coef(CPcomp_pe_mod)[3],
                                                r2_log = coef(CPcomp_pe_mod)[4], 
                                                alpha22 = coef(CPcomp_pe_mod)[5],
                                                alpha21 = coef(CPcomp_pe_mod)[6])
)
#save(CP_min_t2.vv.h12_fixLV_mod, file = "CP_t2.vv.h12_fixLV_mod.Rdata")
Colp_fit6 <- tidy(CP_min_t2.vv.h12_fixLV_mod)

(nrow(CP_t2.vv.h12_PECs)/nrow(CP_t2.vv.h12_fixLV_pe %>% filter(converge == 0))) * 100 # percent of starting combinations within 2 AIC that 
BIC(CP_min_t2.vv.h12_fixLV_mod)

set.seed(1001)
pars.picked = mvrnorm(1000, mu = Colp_fit6$estimate, Sigma = vcov(CP_min_t2.vv.h12_fixLV_mod)) # pick new parameter values by sampling from multivariate normal distribution based on fit
Colp_fit6_pe_confint <- apply(pars.picked,2,quant)


#### 7. Type II + Crowley-Martin #####
##### Search Grid #####
t2.cm.pre_sg <- expand.grid(a_log = seq(-3, -1, length.out = 5),
                            h_log = seq(-5, -1, length.out = 5),
                            log_w = seq(-10, -5, length.out = 5), 
                            sigma = c(0.2), 
                            sigma2 = c(0.2)
)

t2.cm.pre_sg$error <- NA

t2.cm.pre_sg$a_log.est <- NA
t2.cm.pre_sg$h_log.est <- NA
t2.cm.pre_sg$log_w.est <- NA

t2.cm.pre_sg$sigma.est <- NA
t2.cm.pre_sg$sigma2.est <- NA

t2.cm.pre_sg$converge <- NA
t2.cm.pre_sg$BIC <- NA
t2.cm.pre_sg$logLik <- NA
t2.cm.pre_sg$strtpt <- seq(1:nrow(t2.cm.pre_sg))

set.seed(10640)
t2.cm.sim_grid <- t2.cm.pre_sg[sample(t2.cm.pre_sg$strtpt, 98, replace = F), ]


# Set the computer to run in Parallel #
bigstatsr::nb_cores()

my.cluster <- parallel::makeCluster(
  bigstatsr::nb_cores(),
  type = "FORK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()


# Run Grid #
CP_t2.cm_fixLV_pe <- NA # name of file where parallel runs will be save
mcoptions <- list(preschedule=FALSE)

CP_t2.cm_fixLV_pe <- foreach(i = 1:nrow(t2.cm.sim_grid),
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
                              
                              sigma = t2.cm.sim_grid$sigma[[i]],
                              sigma2 = t2.cm.sim_grid$sigma2[[i]]
                 ),
                 
                 data = with(list(N0 = Colp_T0.dens.ml,
                                  Ndead = Colp_T0.dens.ml - Colp_T24.dens.ml,
                                  P = Pred_T0.dens.ml, 
                                  P.end = Pred_T24.dens.ml,
                                  M0 = Para_T0.dens.ml,
                                  M0.end = Para_T24.dens.ml,
                                  Tt = Time.day),
                             data = filter(colp_para_wrk)),
                 
                 fixed = list(m = mort,
                              c_log = -15,
                              
                              r1_log = coef(CPcomp_pe_mod)[1], 
                              alpha11 = coef(CPcomp_pe_mod)[2],
                              alpha12 = coef(CPcomp_pe_mod)[3],
                              r2_log = coef(CPcomp_pe_mod)[4], 
                              alpha22 = coef(CPcomp_pe_mod)[5],
                              alpha21 = coef(CPcomp_pe_mod)[6]) 
      )
      
      if(fit@details$convergence==0) {
        t2.cm.sim_grid$a_log.est[i] <- fit@coef[[1]]
        t2.cm.sim_grid$h_log.est[i] <- fit@coef[[2]]
        t2.cm.sim_grid$log_w.est[i] <- fit@coef[[3]]
        
        t2.cm.sim_grid$sigma.est [i] <- fit@coef[[4]]
        t2.cm.sim_grid$sigma2.est [i] <- fit@coef[[5]]
        
        t2.cm.sim_grid$converge[i] <- fit@details$convergence
        t2.cm.sim_grid$BIC[i] <- BIC(fit)
        t2.cm.sim_grid$logLik[i] <- logLik(fit)[[1]]
      } else {
        t2.cm.sim_grid$a_log.est[i] <- NaN
        t2.cm.sim_grid$h_log.est[i] <- NaN
        t2.cm.sim_grid$log_w.est[i] <- NaN
        
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

CP_t2.cm_fixLV_pe %>% arrange(BIC)

#save(CP_t2.cm_fixLV_pe, file = "CP_t2.cm_fixLV_pe.RData")

###### Best Model #####
(CP_t2.cm_PECs <- CP_t2.cm_fixLV_pe %>% 
   filter(converge == 0) %>% 
   filter(BIC <= min(CP_t2.cm_fixLV_pe$BIC, na.rm = T)+2)) #check for similarity in parameter estimates within 2 AIC of min

min_t2.cm <- CP_t2.cm_PECs %>% filter(BIC <= min(CP_t2.cm_PECs$BIC, na.rm = T)+0.0001) %>% dplyr::select(c(a_log.est:sigma2.est, BIC))

if(nrow(min_t2.cm)>1){
  min_t2.cm <- CP_t2.cm_PECs %>% filter(strtpt %in% sample(filter(CP_t2.cm_PECs, BIC <= min(CP_t2.cm_PECs$BIC, na.rm = T)+0.0001)$strtpt,1))
}

CP_min_t2.cm_fixLV_mod =  mle2(minuslogl = nll.odeint.t2.CM.lv,
                               method= "Nelder-Mead",
                               control = list(maxit = 5000),
                               
                               start = list(a_log = min_t2.cm$a_log.est, 
                                            h_log = min_t2.cm$h_log.est,
                                            log_w = min_t2.cm$log_w.est,
                                            sigma = min_t2.cm$sigma.est,
                                            sigma2 = min_t2.cm$sigma2.est
                               ),
                               
                               data = with(list(N0 = Colp_T0.dens.ml,
                                                Ndead = Colp_T0.dens.ml - Colp_T24.dens.ml,
                                                P = Pred_T0.dens.ml, 
                                                P.end = Pred_T24.dens.ml,
                                                M0 = Para_T0.dens.ml,
                                                M0.end = Para_T24.dens.ml,
                                                Tt = Time.day),
                                           data = filter(colp_para_wrk)),
                               
                               fixed = list(m = mort,
                                            c_log = -15,
                                            
                                            r1_log = coef(CPcomp_pe_mod)[1], 
                                            alpha11 = coef(CPcomp_pe_mod)[2],
                                            alpha12 = coef(CPcomp_pe_mod)[3],
                                            r2_log = coef(CPcomp_pe_mod)[4], 
                                            alpha22 = coef(CPcomp_pe_mod)[5],
                                            alpha21 = coef(CPcomp_pe_mod)[6])
)
#save(CP_min_t2.cm_fixLV_mod, file = "CP_t2.cm_fixLV_mod.Rdata")
Colp_fit7 <- tidy(CP_min_t2.cm_fixLV_mod)

(nrow(CP_t2.cm_PECs)/nrow(CP_t2.cm_fixLV_pe %>% filter(converge == 0))) * 100 # percent of starting combinations within 2 AIC that 
BIC(CP_min_t2.cm_fixLV_mod)

set.seed(1001)
pars.picked = mvrnorm(1000, mu = Colp_fit7$estimate, Sigma = vcov(CP_min_t2.cm_fixLV_mod)) # pick new parameter values by sampling from multivariate normal distribution based on fit
Colp_fit7_pe_confint <- apply(pars.picked,2,quant)


#### 8. Type III #####
##### Search Grid #####
t3.pre_sg <- expand.grid(a_log = seq(-3, -1, length.out = 10), 
                         h_log = seq(-5, -1, length.out = 10),
                         kr_log = seq(-3, -0.5, length.out = 5),
                         sigma = c(0.2), 
                         sigma2 = c(0.2)
)

t3.pre_sg$error <- NA

t3.pre_sg$a_log.est <- NA
t3.pre_sg$h_log.est <- NA
t3.pre_sg$kr_log.est <- NA

t3.pre_sg$sigma.est <- NA
t3.pre_sg$sigma2.est <- NA

t3.pre_sg$converge <- NA
t3.pre_sg$BIC <- NA
t3.pre_sg$logLik <- NA
t3.pre_sg$strtpt <- seq(1:nrow(t3.pre_sg))

set.seed(52882)
t3.sim_grid <- t3.pre_sg[sample(t3.pre_sg$strtpt, 98, replace = F), ]


# Set the computer to run in Parallel #
bigstatsr::nb_cores()

my.cluster <- parallel::makeCluster(
  bigstatsr::nb_cores(),
  type = "FORK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()


# Run Grid #
CP_t3_fixLV_pe <- NA # name of file where parallel runs will be save
mcoptions <- list(preschedule=FALSE)

CP_t3_fixLV_pe <- foreach(i = 1:nrow(t3.sim_grid),
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
                              
                              sigma = t3.sim_grid$sigma[[i]],
                              sigma2 = t3.sim_grid$sigma2[[i]]
                 ),
                 
                 data = with(list(N0 = Colp_T0.dens.ml,
                                  Ndead = Colp_T0.dens.ml - Colp_T24.dens.ml,
                                  P = Pred_T0.dens.ml, 
                                  P.end = Pred_T24.dens.ml,
                                  M0 = Para_T0.dens.ml,
                                  M0.end = Para_T24.dens.ml,
                                  Tt = Time.day),
                             data = filter(colp_para_wrk)),
                 
                 fixed = list(m = mort,
                              c_log = -15,
                              a12 = 0,
                              h12 = 0,
                              
                              r1_log = coef(CPcomp_pe_mod)[1], 
                              alpha11 = coef(CPcomp_pe_mod)[2],
                              alpha12 = coef(CPcomp_pe_mod)[3],
                              r2_log = coef(CPcomp_pe_mod)[4], 
                              alpha22 = coef(CPcomp_pe_mod)[5],
                              alpha21 = coef(CPcomp_pe_mod)[6]) 
      )
      
      if(fit@details$convergence==0) {
        t3.sim_grid$a_log.est[i] <- fit@coef[[1]]
        t3.sim_grid$h_log.est[i] <- fit@coef[[2]]
        t3.sim_grid$kr_log.est[i] <- fit@coef[[3]]
        
        t3.sim_grid$sigma.est [i] <- fit@coef[[4]]
        t3.sim_grid$sigma2.est [i] <- fit@coef[[5]]
        
        t3.sim_grid$converge[i] <- fit@details$convergence
        t3.sim_grid$BIC[i] <- BIC(fit)
        t3.sim_grid$logLik[i] <- logLik(fit)[[1]]
      } else {
        t3.sim_grid$a_log.est[i] <- NaN
        t3.sim_grid$h_log.est[i] <- NaN
        t3.sim_grid$kr_log.est[i] <- NaN
        
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

CP_t3_fixLV_pe %>% arrange(BIC)

#save(CP_t3_fixLV_pe, file = "CP_t3_fixLV_pe.RData")

###### Best Model #####
(CP_t3_PECs <- CP_t3_fixLV_pe %>% 
   filter(converge == 0) %>% 
   filter(BIC <= min(CP_t3_fixLV_pe$BIC, na.rm = T)+2)) #check for similarity in parameter estimates within 2 AIC of min

min_T3 <- CP_t3_PECs %>% filter(BIC <= min(CP_t3_PECs$BIC, na.rm = T)+0.0001) %>% dplyr::select(c(a_log.est:sigma2.est, BIC))

if(nrow(min_T3)>1){
  min_T3 <- CP_t3_PECs %>% filter(strtpt %in% sample(filter(CP_t3_PECs, BIC <= min(CP_t3_PECs$BIC, na.rm = T)+0.0001)$strtpt,1))
}

CP_min_t3_fixLV_mod =  mle2(minuslogl = nll.odeint.t3.a1.h1.lv,
                            method= "Nelder-Mead",
                            control = list(maxit = 5000),
                            
                            start = list(a_log = min_T3$a_log.est, 
                                         h_log = min_T3$h_log.est,
                                         kr_log = min_T3$kr_log.est,
                                         sigma = min_T3$sigma.est,
                                         sigma2 = min_T3$sigma2.est
                            ),
                            
                            data = with(list(N0 = Colp_T0.dens.ml,
                                             Ndead = Colp_T0.dens.ml - Colp_T24.dens.ml,
                                             P = Pred_T0.dens.ml, 
                                             P.end = Pred_T24.dens.ml,
                                             M0 = Para_T0.dens.ml,
                                             M0.end = Para_T24.dens.ml,
                                             Tt = Time.day),
                                        data = filter(colp_para_wrk)),
                            
                            fixed = list(m = mort,
                                         c_log = -15,
                                         a12 = 0,
                                         h12 = 0,
                                         r1_log = coef(CPcomp_pe_mod)[1], 
                                         alpha11 = coef(CPcomp_pe_mod)[2],
                                         alpha12 = coef(CPcomp_pe_mod)[3],
                                         r2_log = coef(CPcomp_pe_mod)[4], 
                                         alpha22 = coef(CPcomp_pe_mod)[5],
                                         alpha21 = coef(CPcomp_pe_mod)[6])
)
#save(CP_min_t3_fixLV_mod, file = "CP_t3_fixLV_mod.Rdata")
Colp_fit8 <- tidy(CP_min_t3_fixLV_mod)

(nrow(CP_t3_PECs)/nrow(CP_t3_fixLV_pe %>% filter(converge == 0))) * 100 # percent of starting combinations within 2 AIC that 
BIC(CP_min_t3_fixLV_mod)

set.seed(1001)
pars.picked = mvrnorm(1000, mu = Colp_fit8$estimate, Sigma = vcov(CP_min_t3_fixLV_mod)) # pick new parameter values by sampling from multivariate normal distribution based on fit
Colp_fit8_pe_confint <- apply(pars.picked,2,quant)

#### 9. Type III + linear a TIM #####
##### Search Grid #####
t3.a12.pre_sg <- expand.grid(a_log = seq(-3, -1, length.out = 5), 
                             a12 = seq(-0.0005, 0.0001, length.out = 5),
                             h_log = seq(-5, -1, length.out = 5),
                             kr_log = seq(-3, -0.5, length.out = 5),
                             sigma = c(0.2), 
                             sigma2 = c(0.2)
)

t3.a12.pre_sg$error <- NA

t3.a12.pre_sg$a_log.est <- NA
t3.a12.pre_sg$a12.est <- NA
t3.a12.pre_sg$h_log.est <- NA
t3.a12.pre_sg$kr_log.est <- NA

t3.a12.pre_sg$sigma.est <- NA
t3.a12.pre_sg$sigma2.est <- NA

t3.a12.pre_sg$converge <- NA
t3.a12.pre_sg$BIC <- NA
t3.a12.pre_sg$logLik <- NA
t3.a12.pre_sg$strtpt <- seq(1:nrow(t3.a12.pre_sg))

set.seed(22761)
t3.a12.sim_grid <- t3.a12.pre_sg[sample(t3.a12.pre_sg$strtpt, 98, replace = F), ]


# Set the computer to run in Parallel #
bigstatsr::nb_cores()

my.cluster <- parallel::makeCluster(
  bigstatsr::nb_cores(),
  type = "FORK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()


# Run Grid #
CP_t3.a12_fixLV_pe <- NA # name of file where parallel runs will be save
mcoptions <- list(preschedule=FALSE)

CP_t3.a12_fixLV_pe <- foreach(i = 1:nrow(t3.a12.sim_grid),
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
                              
                              sigma = t3.a12.sim_grid$sigma[[i]],
                              sigma2 = t3.a12.sim_grid$sigma2[[i]]
                 ),
                 
                 data = with(list(N0 = Colp_T0.dens.ml,
                                  Ndead = Colp_T0.dens.ml - Colp_T24.dens.ml,
                                  P = Pred_T0.dens.ml, 
                                  P.end = Pred_T24.dens.ml,
                                  M0 = Para_T0.dens.ml,
                                  M0.end = Para_T24.dens.ml,
                                  Tt = Time.day),
                             data = filter(colp_para_wrk)),
                 
                 fixed = list(m = mort,
                              c_log = -15,
                              h12 = 0,
                              
                              r1_log = coef(CPcomp_pe_mod)[1], 
                              alpha11 = coef(CPcomp_pe_mod)[2],
                              alpha12 = coef(CPcomp_pe_mod)[3],
                              r2_log = coef(CPcomp_pe_mod)[4], 
                              alpha22 = coef(CPcomp_pe_mod)[5],
                              alpha21 = coef(CPcomp_pe_mod)[6]) 
      )
      
      if(fit@details$convergence==0) {
        t3.a12.sim_grid$a_log.est[i] <- fit@coef[[1]]
        t3.a12.sim_grid$a12.est[i] <- fit@coef[[2]]
        t3.a12.sim_grid$h_log.est[i] <- fit@coef[[3]]
        t3.a12.sim_grid$kr_log.est[i] <- fit@coef[[4]]
        
        t3.a12.sim_grid$sigma.est [i] <- fit@coef[[5]]
        t3.a12.sim_grid$sigma2.est [i] <- fit@coef[[6]]
        
        t3.a12.sim_grid$converge[i] <- fit@details$convergence
        t3.a12.sim_grid$BIC[i] <- BIC(fit)
        t3.a12.sim_grid$logLik[i] <- logLik(fit)[[1]]
      } else {
        t3.a12.sim_grid$a_log.est[i] <- NaN
        t3.a12.sim_grid$a12.est[i] <- NaN
        t3.a12.sim_grid$h_log.est[i] <- NaN
        t3.a12.sim_grid$kr_log.est[i] <- NaN
        
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

CP_t3.a12_fixLV_pe %>% arrange(BIC)
#save(CP_t3.a12_fixLV_pe, file = "CP_t3.a12_fixLV_pe.RData")

###### Best Model #####
(CP_t3.a12_PECs <- CP_t3.a12_fixLV_pe %>% 
   filter(converge == 0) %>% 
   filter(BIC <= min(CP_t3.a12_fixLV_pe$BIC, na.rm = T)+2)) #check for similarity in parameter estimates within 2 AIC of min

min_t3.a12 <- CP_t3.a12_PECs %>% filter(BIC <= min(CP_t3.a12_PECs$BIC, na.rm = T)+0.0001) %>% dplyr::select(c(a_log.est:sigma2.est, BIC))

if(nrow(min_t3.a12)>1){
  min_t3.a12 <- CP_t3.a12_PECs %>% filter(strtpt %in% sample(filter(CP_t3.a12_PECs, BIC <= min(CP_t3.a12_PECs$BIC, na.rm = T)+0.0001)$strtpt,1))
}

CP_min_t3.a12_fixLV_mod =  mle2(minuslogl = nll.odeint.t3.a1.h1.lv,
                                method= "Nelder-Mead",
                                control = list(maxit = 5000),
                                
                                start = list(a_log = min_t3.a12$a_log.est, 
                                             h_log = min_t3.a12$h_log.est,
                                             a12 = min_t3.a12$a12.est,
                                             kr_log = min_t3.a12$kr_log.est,
                                             sigma = min_t3.a12$sigma.est,
                                             sigma2 = min_t3.a12$sigma2.est
                                ),
                                
                                data = with(list(N0 = Colp_T0.dens.ml,
                                                 Ndead = Colp_T0.dens.ml - Colp_T24.dens.ml,
                                                 P = Pred_T0.dens.ml, 
                                                 P.end = Pred_T24.dens.ml,
                                                 M0 = Para_T0.dens.ml,
                                                 M0.end = Para_T24.dens.ml,
                                                 Tt = Time.day),
                                            data = filter(colp_para_wrk)),
                                
                                fixed = list(m = mort,
                                             c_log = -15,
                                             h12 = 0, 
                                             r1_log = coef(CPcomp_pe_mod)[1], 
                                             alpha11 = coef(CPcomp_pe_mod)[2],
                                             alpha12 = coef(CPcomp_pe_mod)[3],
                                             r2_log = coef(CPcomp_pe_mod)[4], 
                                             alpha22 = coef(CPcomp_pe_mod)[5],
                                             alpha21 = coef(CPcomp_pe_mod)[6])
)
#save(CP_min_t3.a12_fixLV_mod, file = "CP_t3.a12_fixLV_mod.Rdata")
Colp_fit9 <- tidy(CP_min_t3.a12_fixLV_mod)

(nrow(CP_t3.a12_PECs)/nrow(CP_t3.a12_fixLV_pe %>% filter(converge == 0))) * 100 # percent of starting combinations within 2 AIC that 
BIC(CP_min_t3.a12_fixLV_mod)

set.seed(1001)
pars.picked = mvrnorm(1000, mu = Colp_fit9$estimate, Sigma = vcov(CP_min_t3.a12_fixLV_mod)) # pick new parameter values by sampling from multivariate normal distribution based on fit
Colp_fit9_pe_confint <- apply(pars.picked,2,quant)


#### 10. Type III + linear h TIM #####
##### Search Grid #####
t3.h12.pre_sg <- expand.grid(a_log = seq(-3, -1, length.out = 5), 
                             h12 = seq(-0.0005, 0.0001, length.out = 5),
                             h_log = seq(-5, -1, length.out = 5),
                             kr_log = seq(-3, -0.5, length.out = 5),
                             sigma = c(0.2), 
                             sigma2 = c(0.2)
)

t3.h12.pre_sg$error <- NA

t3.h12.pre_sg$a_log.est <- NA
t3.h12.pre_sg$h_log.est <- NA
t3.h12.pre_sg$h12.est <- NA
t3.h12.pre_sg$kr_log.est <- NA

t3.h12.pre_sg$sigma.est <- NA
t3.h12.pre_sg$sigma2.est <- NA

t3.h12.pre_sg$converge <- NA
t3.h12.pre_sg$BIC <- NA
t3.h12.pre_sg$logLik <- NA
t3.h12.pre_sg$strtpt <- seq(1:nrow(t3.h12.pre_sg))

set.seed(34432)
t3.h12.sim_grid <- t3.h12.pre_sg[sample(t3.h12.pre_sg$strtpt, 98, replace = F), ]


# Set the computer to run in Parallel #
bigstatsr::nb_cores()

my.cluster <- parallel::makeCluster(
  bigstatsr::nb_cores(),
  type = "FORK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()


# Run Grid #
CP_t3.h12_fixLV_pe <- NA # name of file where parallel runs will be save
mcoptions <- list(preschedule=FALSE)

CP_t3.h12_fixLV_pe <- foreach(i = 1:nrow(t3.h12.sim_grid),
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
                              
                              sigma = t3.h12.sim_grid$sigma[[i]],
                              sigma2 = t3.h12.sim_grid$sigma2[[i]]
                 ),
                 
                 data = with(list(N0 = Colp_T0.dens.ml,
                                  Ndead = Colp_T0.dens.ml - Colp_T24.dens.ml,
                                  P = Pred_T0.dens.ml, 
                                  P.end = Pred_T24.dens.ml,
                                  M0 = Para_T0.dens.ml,
                                  M0.end = Para_T24.dens.ml,
                                  Tt = Time.day),
                             data = filter(colp_para_wrk)),
                 
                 fixed = list(m = mort,
                              c_log = -15,
                              a12 = 0,
                              
                              r1_log = coef(CPcomp_pe_mod)[1], 
                              alpha11 = coef(CPcomp_pe_mod)[2],
                              alpha12 = coef(CPcomp_pe_mod)[3],
                              r2_log = coef(CPcomp_pe_mod)[4], 
                              alpha22 = coef(CPcomp_pe_mod)[5],
                              alpha21 = coef(CPcomp_pe_mod)[6]) 
      )
      
      if(fit@details$convergence==0) {
        t3.h12.sim_grid$a_log.est[i] <- fit@coef[[1]]
        t3.h12.sim_grid$h_log.est[i] <- fit@coef[[2]]
        t3.h12.sim_grid$h12.est[i] <- fit@coef[[3]]
        t3.h12.sim_grid$kr_log.est[i] <- fit@coef[[4]]
        
        t3.h12.sim_grid$sigma.est [i] <- fit@coef[[5]]
        t3.h12.sim_grid$sigma2.est [i] <- fit@coef[[6]]
        
        t3.h12.sim_grid$converge[i] <- fit@details$convergence
        t3.h12.sim_grid$BIC[i] <- BIC(fit)
        t3.h12.sim_grid$logLik[i] <- logLik(fit)[[1]]
      } else {
        t3.h12.sim_grid$a_log.est[i] <- NaN
        t3.h12.sim_grid$h_log.est[i] <- NaN
        t3.h12.sim_grid$h12.est[i] <- NaN
        
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

CP_t3.h12_fixLV_pe %>% arrange(BIC)

#save(CP_t3.h12_fixLV_pe, file = "CP_t3.h12_fixLV_pe.RData")

###### Best Model #####
(CP_t3.h12_PECs <- CP_t3.h12_fixLV_pe %>% 
   filter(converge == 0) %>% 
   filter(BIC <= min(CP_t3.h12_fixLV_pe$BIC, na.rm = T)+2)) #check for similarity in parameter estimates within 2 AIC of min

min_t3.h12 <- CP_t3.h12_PECs %>% filter(BIC <= min(CP_t3.h12_PECs$BIC, na.rm = T)+0.0001) %>% dplyr::select(c(a_log.est:sigma2.est, BIC))

if(nrow(min_t3.h12)>1){
  min_t3.h12 <- CP_t3.h12_PECs %>% filter(strtpt %in% sample(filter(CP_t3.h12_PECs, BIC <= min(CP_t3.h12_PECs$BIC, na.rm = T)+0.0001)$strtpt,1))
}

CP_min_t3.h12_fixLV_mod =  mle2(minuslogl = nll.odeint.t3.a1.h1.lv,
                                method= "Nelder-Mead",
                                control = list(maxit = 5000),
                                
                                start = list(a_log = min_t3.h12$a_log.est, 
                                             h_log = min_t3.h12$h_log.est,
                                             h12 = min_t3.h12$h12.est,
                                             kr_log = min_t3.h12$kr_log.est,
                                             sigma = min_t3.h12$sigma.est,
                                             sigma2 = min_t3.h12$sigma2.est
                                ),
                                
                                data = with(list(N0 = Colp_T0.dens.ml,
                                                 Ndead = Colp_T0.dens.ml - Colp_T24.dens.ml,
                                                 P = Pred_T0.dens.ml, 
                                                 P.end = Pred_T24.dens.ml,
                                                 M0 = Para_T0.dens.ml,
                                                 M0.end = Para_T24.dens.ml,
                                                 Tt = Time.day),
                                            data = filter(colp_para_wrk)),
                                
                                fixed = list(m = mort,
                                             c_log = -15,
                                             a12=0,
                                             r1_log = coef(CPcomp_pe_mod)[1], 
                                             alpha11 = coef(CPcomp_pe_mod)[2],
                                             alpha12 = coef(CPcomp_pe_mod)[3],
                                             r2_log = coef(CPcomp_pe_mod)[4], 
                                             alpha22 = coef(CPcomp_pe_mod)[5],
                                             alpha21 = coef(CPcomp_pe_mod)[6])
)
#save(CP_min_t3.h12_fixLV_mod, file = "CP_t3.h12_fixLV_mod.Rdata")
Colp_fit10 <- tidy(CP_min_t3.h12_fixLV_mod)

(nrow(CP_t3.h12_PECs)/nrow(CP_t3.h12_fixLV_pe %>% filter(converge == 0))) * 100 # percent of starting combinations within 2 AIC that 
BIC(CP_min_t3.h12_fixLV_mod)

set.seed(1001)
pars.picked = mvrnorm(1000, mu = Colp_fit10$estimate, Sigma = vcov(CP_min_t3.h12_fixLV_mod)) # pick new parameter values by sampling from multivariate normal distribution based on fit
Colp_fit10_pe_confint <- apply(pars.picked,2,quant)

#### 11. Type III + linear a TIM  + linear h TIM #####
##### Search Grid #####
t3.a12.h12.pre_sg <- expand.grid(a_log = seq(-3, -1, length.out = 5), 
                                 a12 = seq(-0.0005, 0.0001, length.out = 5),
                                 h12 = seq(-0.0005, 0.0001, length.out = 5),
                                 h_log = seq(-5, -1, length.out = 5),
                                 kr_log = seq(-3, -0.5, length.out = 5),
                                 sigma = c(0.2), 
                                 sigma2 = c(0.2)
)

t3.a12.h12.pre_sg$error <- NA

t3.a12.h12.pre_sg$a_log.est <- NA
t3.a12.h12.pre_sg$a12.est <- NA
t3.a12.h12.pre_sg$h_log.est <- NA
t3.a12.h12.pre_sg$h12.est <- NA
t3.a12.h12.pre_sg$kr_log.est <- NA

t3.a12.h12.pre_sg$sigma.est <- NA
t3.a12.h12.pre_sg$sigma2.est <- NA

t3.a12.h12.pre_sg$converge <- NA
t3.a12.h12.pre_sg$BIC <- NA
t3.a12.h12.pre_sg$logLik <- NA
t3.a12.h12.pre_sg$strtpt <- seq(1:nrow(t3.a12.h12.pre_sg))

set.seed(31394)
t3.a12.h12.sim_grid <- t3.a12.h12.pre_sg[sample(t3.a12.h12.pre_sg$strtpt, 98, replace = F), ]


# Set the computer to run in Parallel #
bigstatsr::nb_cores()

my.cluster <- parallel::makeCluster(
  bigstatsr::nb_cores(),
  type = "FORK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()


# Run Grid #
CP_t3.a12.h12_fixLV_pe <- NA # name of file where parallel runs will be save
mcoptions <- list(preschedule=FALSE)

CP_t3.a12.h12_fixLV_pe <- foreach(i = 1:nrow(t3.a12.h12.sim_grid),
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
                              
                              sigma = t3.a12.h12.sim_grid$sigma[[i]],
                              sigma2 = t3.a12.h12.sim_grid$sigma2[[i]]
                 ),
                 
                 data = with(list(N0 = Colp_T0.dens.ml,
                                  Ndead = Colp_T0.dens.ml - Colp_T24.dens.ml,
                                  P = Pred_T0.dens.ml, 
                                  P.end = Pred_T24.dens.ml,
                                  M0 = Para_T0.dens.ml,
                                  M0.end = Para_T24.dens.ml,
                                  Tt = Time.day),
                             data = filter(colp_para_wrk)),
                 
                 fixed = list(m = mort,
                              c_log = -15,
                              
                              r1_log = coef(CPcomp_pe_mod)[1], 
                              alpha11 = coef(CPcomp_pe_mod)[2],
                              alpha12 = coef(CPcomp_pe_mod)[3],
                              r2_log = coef(CPcomp_pe_mod)[4], 
                              alpha22 = coef(CPcomp_pe_mod)[5],
                              alpha21 = coef(CPcomp_pe_mod)[6]) 
      )
      
      if(fit@details$convergence==0) {
        t3.a12.h12.sim_grid$a_log.est[i] <- fit@coef[[1]]
        t3.a12.h12.sim_grid$a12.est[i] <- fit@coef[[2]]
        t3.a12.h12.sim_grid$h_log.est[i] <- fit@coef[[3]]
        t3.a12.h12.sim_grid$h12.est[i] <- fit@coef[[4]]
        t3.a12.h12.sim_grid$kr_log.est[i] <- fit@coef[[5]]
        
        t3.a12.h12.sim_grid$sigma.est [i] <- fit@coef[[6]]
        t3.a12.h12.sim_grid$sigma2.est [i] <- fit@coef[[7]]
        
        t3.a12.h12.sim_grid$converge[i] <- fit@details$convergence
        t3.a12.h12.sim_grid$BIC[i] <- BIC(fit)
        t3.a12.h12.sim_grid$logLik[i] <- logLik(fit)[[1]]
      } else {
        t3.a12.h12.sim_grid$a_log.est[i] <- NaN
        t3.a12.h12.sim_grid$a12.est[i] <- NaN
        t3.a12.h12.sim_grid$h_log.est[i] <- NaN
        t3.a12.h12.sim_grid$a12.est[i] <- NaN
        t3.a12.h12.sim_grid$kr_log.est[i] <- NaN
        
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

CP_t3.a12.h12_fixLV_pe %>% arrange(BIC)

#save(CP_t3.a12.h12_fixLV_pe, file = "CP_t3.a12.h12_fixLV_pe.RData")

###### Best Model #####
(CP_t3.a12.h12_PECs <- CP_t3.a12.h12_fixLV_pe %>% 
   filter(converge == 0) %>% 
   filter(BIC <= min(CP_t3.a12.h12_fixLV_pe$BIC, na.rm = T)+2)) #check for similarity in parameter estimates within 2 AIC of min

min_t3.a12.h12 <- CP_t3.a12.h12_PECs %>% filter(BIC <= min(CP_t3.a12.h12_PECs$BIC, na.rm = T)+0.0001) %>% dplyr::select(c(a_log.est:sigma2.est, BIC))

if(nrow(min_t3.a12.h12)>1){
  min_t3.a12.h12 <- CP_t3.a12.h12_PECs %>% filter(strtpt %in% sample(filter(CP_t3.a12.h12_PECs, BIC <= min(CP_t3.a12.h12_PECs$BIC, na.rm = T)+0.0001)$strtpt,1))
}

CP_min_t3.a12.h12_fixLV_mod =  mle2(minuslogl = nll.odeint.t3.a1.h1.lv,
                                    method= "Nelder-Mead",
                                    control = list(maxit = 5000),
                                    
                                    start = list(a_log = min_t3.a12.h12$a_log.est, 
                                                 a12 = min_t3.a12.h12$a12.est,
                                                 h_log = min_t3.a12.h12$h_log.est,
                                                 h12 = min_t3.a12.h12$h12.est,
                                                 kr_log = min_t3.a12.h12$kr_log.est,
                                                 sigma = min_t3.a12.h12$sigma.est,
                                                 sigma2 = min_t3.a12.h12$sigma2.est
                                    ),
                                    
                                    data = with(list(N0 = Colp_T0.dens.ml,
                                                     Ndead = Colp_T0.dens.ml - Colp_T24.dens.ml,
                                                     P = Pred_T0.dens.ml, 
                                                     P.end = Pred_T24.dens.ml,
                                                     M0 = Para_T0.dens.ml,
                                                     M0.end = Para_T24.dens.ml,
                                                     Tt = Time.day),
                                                data = filter(colp_para_wrk)),
                                    
                                    fixed = list(m = mort,
                                                 c_log = -15,
                                                 r1_log = coef(CPcomp_pe_mod)[1], 
                                                 alpha11 = coef(CPcomp_pe_mod)[2],
                                                 alpha12 = coef(CPcomp_pe_mod)[3],
                                                 r2_log = coef(CPcomp_pe_mod)[4], 
                                                 alpha22 = coef(CPcomp_pe_mod)[5],
                                                 alpha21 = coef(CPcomp_pe_mod)[6])
)
#save(CP_min_t3.a12.h12_fixLV_mod, file = "CP_t3.a12.h12_fixLV_mod.Rdata")
Colp_fit11 <- tidy(CP_min_t3.a12.h12_fixLV_mod)

(nrow(CP_t3.a12.h12_PECs)/nrow(CP_t3.a12.h12_fixLV_pe %>% filter(converge == 0))) * 100 # percent of starting combinations within 2 AIC that 
BIC(CP_min_t3.a12.h12_fixLV_mod)

set.seed(1001)
pars.picked = mvrnorm(1000, mu = Colp_fit11$estimate, Sigma = vcov(CP_min_t3.a12.h12_fixLV_mod)) # pick new parameter values by sampling from multivariate normal distribution based on fit
Colp_fit11_pe_confint <- apply(pars.picked,2,quant)


#### 12. Type III + van Veen #####
##### Search Grid #####
t3.vv.pre_sg <- expand.grid(a_log = seq(-3, -1, length.out = 5), 
                            h_log = seq(-5, -1, length.out = 5), 
                            log_w = seq(-10, -5, length.out = 5), 
                            kr_log = seq(-3, -0.5, length.out = 5),
                            sigma = c(0.2), 
                            sigma2 = c(0.2)
)

t3.vv.pre_sg$error <- NA

t3.vv.pre_sg$a_log.est <- NA
t3.vv.pre_sg$h_log.est <- NA
t3.vv.pre_sg$log_w.est <- NA
t3.vv.pre_sg$kr_log.est <- NA

t3.vv.pre_sg$sigma.est <- NA
t3.vv.pre_sg$sigma2.est <- NA

t3.vv.pre_sg$converge <- NA
t3.vv.pre_sg$BIC <- NA
t3.vv.pre_sg$logLik <- NA
t3.vv.pre_sg$strtpt <- seq(1:nrow(t3.vv.pre_sg))

set.seed(96542)
t3.vv.sim_grid <- t3.vv.pre_sg[sample(t3.vv.pre_sg$strtpt, 98, replace = F), ]


# Set the computer to run in Parallel #
bigstatsr::nb_cores()

my.cluster <- parallel::makeCluster(
  bigstatsr::nb_cores(),
  type = "FORK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()


# Run Grid #
CP_t3.vv_fixLV_pe <- NA # name of file where parallel runs will be save
mcoptions <- list(preschedule=FALSE)

CP_t3.vv_fixLV_pe <- foreach(i = 1:nrow(t3.vv.sim_grid),
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
                              
                              sigma = t3.vv.sim_grid$sigma[[i]],
                              sigma2 = t3.vv.sim_grid$sigma2[[i]]
                 ),
                 
                 data = with(list(N0 = Colp_T0.dens.ml,
                                  Ndead = Colp_T0.dens.ml - Colp_T24.dens.ml,
                                  P = Pred_T0.dens.ml, 
                                  P.end = Pred_T24.dens.ml,
                                  M0 = Para_T0.dens.ml,
                                  M0.end = Para_T24.dens.ml,
                                  Tt = Time.day),
                             data = filter(colp_para_wrk)),
                 
                 fixed = list(m = mort,
                              c_log = -15,
                              h12 = 0,
                              
                              r1_log = coef(CPcomp_pe_mod)[1], 
                              alpha11 = coef(CPcomp_pe_mod)[2],
                              alpha12 = coef(CPcomp_pe_mod)[3],
                              r2_log = coef(CPcomp_pe_mod)[4], 
                              alpha22 = coef(CPcomp_pe_mod)[5],
                              alpha21 = coef(CPcomp_pe_mod)[6]) 
      )
      
      if(fit@details$convergence==0) {
        t3.vv.sim_grid$a_log.est[i] <- fit@coef[[1]]
        t3.vv.sim_grid$h_log.est[i] <- fit@coef[[2]]
        t3.vv.sim_grid$log_w.est[i] <- fit@coef[[3]]
        t3.vv.sim_grid$kr_log.est[i] <- fit@coef[[4]]
        
        t3.vv.sim_grid$sigma.est [i] <- fit@coef[[5]]
        t3.vv.sim_grid$sigma2.est [i] <- fit@coef[[6]]
        
        t3.vv.sim_grid$converge[i] <- fit@details$convergence
        t3.vv.sim_grid$BIC[i] <- BIC(fit)
        t3.vv.sim_grid$logLik[i] <- logLik(fit)[[1]]
      } else {
        t3.vv.sim_grid$a_log.est[i] <- NaN
        t3.vv.sim_grid$h_log.est[i] <- NaN
        t3.vv.sim_grid$log_w.est[i] <- NaN
        t3.vv.sim_grid$kr_log.est[i] <- NaN
        
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

CP_t3.vv_fixLV_pe %>% arrange(BIC)

#save(CP_t3.vv_fixLV_pe, file = "CP_t3.vv_fixLV_pe.RData")

###### Best Model #####
(CP_t3.vv_PECs <- CP_t3.vv_fixLV_pe %>% 
   filter(converge == 0) %>% 
   filter(BIC <= min(CP_t3.vv_fixLV_pe$BIC, na.rm = T)+2)) #check for similarity in parameter estimates within 2 AIC of min

min_t3.vv <- CP_t3.vv_PECs %>% filter(BIC <= min(CP_t3.vv_PECs$BIC, na.rm = T)+0.0001) %>% dplyr::select(c(a_log.est:sigma2.est, BIC))

if(nrow(min_t3.vv)>1){
  min_t3.vv <- CP_t3.vv_PECs %>% filter(strtpt %in% sample(filter(CP_t3.vv_PECs, BIC <= min(CP_t3.vv_PECs$BIC, na.rm = T)+0.0001)$strtpt,1))
}

CP_min_t3.vv_fixLV_mod =  mle2(minuslogl = nll.odeint.T3.h1.VV.lv,
                               method= "Nelder-Mead",
                               control = list(maxit = 5000),
                               
                               start = list(a_log = min_t3.vv$a_log.est, 
                                            h_log = min_t3.vv$h_log.est,
                                            log_w = min_t3.vv$log_w.est,
                                            kr_log = min_t3.vv$kr_log.est,
                                            sigma = min_t3.vv$sigma.est,
                                            sigma2 = min_t3.vv$sigma2.est
                               ),
                               
                               data = with(list(N0 = Colp_T0.dens.ml,
                                                Ndead = Colp_T0.dens.ml - Colp_T24.dens.ml,
                                                P = Pred_T0.dens.ml, 
                                                P.end = Pred_T24.dens.ml,
                                                M0 = Para_T0.dens.ml,
                                                M0.end = Para_T24.dens.ml,
                                                Tt = Time.day),
                                           data = filter(colp_para_wrk)),
                               
                               fixed = list(m = mort,
                                            c_log = -15,
                                            h12=0,
                                            
                                            r1_log = coef(CPcomp_pe_mod)[1], 
                                            alpha11 = coef(CPcomp_pe_mod)[2],
                                            alpha12 = coef(CPcomp_pe_mod)[3],
                                            r2_log = coef(CPcomp_pe_mod)[4], 
                                            alpha22 = coef(CPcomp_pe_mod)[5],
                                            alpha21 = coef(CPcomp_pe_mod)[6])
)
#save(CP_min_t3.vv_fixLV_mod, file = "CP_t3.vv_fixLV_mod.Rdata")
Colp_fit12 <- tidy(CP_min_t3.vv_fixLV_mod)

(nrow(CP_t3.vv_PECs)/nrow(CP_t3.vv_fixLV_pe %>% filter(converge == 0))) * 100 # percent of starting combinations within 2 AIC that 
BIC(CP_min_t3.vv_fixLV_mod)

set.seed(1001)
pars.picked = mvrnorm(1000, mu = Colp_fit12$estimate, Sigma = vcov(CP_min_t3.vv_fixLV_mod)) # pick new parameter values by sampling from multivariate normal distribution based on fit
Colp_fit12_pe_confint <- apply(pars.picked,2,quant)


#### 13. Type III + Crowley-Martin #####
##### Search Grid #####
t3.cm.pre_sg <- expand.grid(a_log = seq(-3, -1, length.out = 5),
                            h_log = seq(-5, -1, length.out = 5),
                            log_w = seq(-10, -5, length.out = 5),
                            kr_log = seq(-3, -0.5, length.out = 5),
                            sigma = c(0.2), 
                            sigma2 = c(0.2)
)

t3.cm.pre_sg$error <- NA

t3.cm.pre_sg$a_log.est <- NA
t3.cm.pre_sg$h_log.est <- NA
t3.cm.pre_sg$log_w.est <- NA
t3.cm.pre_sg$kr_log.est <- NA

t3.cm.pre_sg$sigma.est <- NA
t3.cm.pre_sg$sigma2.est <- NA

t3.cm.pre_sg$converge <- NA
t3.cm.pre_sg$BIC <- NA
t3.cm.pre_sg$logLik <- NA
t3.cm.pre_sg$strtpt <- seq(1:nrow(t3.cm.pre_sg))

set.seed(23396)
t3.cm.sim_grid <- t3.cm.pre_sg[sample(t3.cm.pre_sg$strtpt, 98, replace = F), ]


# Set the computer to run in Parallel #
bigstatsr::nb_cores()

my.cluster <- parallel::makeCluster(
  bigstatsr::nb_cores(),
  type = "FORK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()


# Run Grid #
CP_t3.cm_fixLV_pe <- NA # name of file where parallel runs will be save
mcoptions <- list(preschedule=FALSE)

CP_t3.cm_fixLV_pe <- foreach(i = 1:nrow(t3.cm.sim_grid),
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
                              
                              sigma = t3.cm.sim_grid$sigma[[i]],
                              sigma2 = t3.cm.sim_grid$sigma2[[i]]
                 ),
                 
                 data = with(list(N0 = Colp_T0.dens.ml,
                                  Ndead = Colp_T0.dens.ml - Colp_T24.dens.ml,
                                  P = Pred_T0.dens.ml, 
                                  P.end = Pred_T24.dens.ml,
                                  M0 = Para_T0.dens.ml,
                                  M0.end = Para_T24.dens.ml,
                                  Tt = Time.day),
                             data = filter(colp_para_wrk)),
                 
                 fixed = list(m = mort,
                              c_log = -15,
                              
                              r1_log = coef(CPcomp_pe_mod)[1], 
                              alpha11 = coef(CPcomp_pe_mod)[2],
                              alpha12 = coef(CPcomp_pe_mod)[3],
                              r2_log = coef(CPcomp_pe_mod)[4], 
                              alpha22 = coef(CPcomp_pe_mod)[5],
                              alpha21 = coef(CPcomp_pe_mod)[6]) 
      )
      
      if(fit@details$convergence==0) {
        t3.cm.sim_grid$a_log.est[i] <- fit@coef[[1]]
        t3.cm.sim_grid$h_log.est[i] <- fit@coef[[2]]
        t3.cm.sim_grid$log_w.est[i] <- fit@coef[[3]]
        t3.cm.sim_grid$kr_log.est[i] <- fit@coef[[4]]
        
        t3.cm.sim_grid$sigma.est [i] <- fit@coef[[5]]
        t3.cm.sim_grid$sigma2.est [i] <- fit@coef[[6]]
        
        t3.cm.sim_grid$converge[i] <- fit@details$convergence
        t3.cm.sim_grid$BIC[i] <- BIC(fit)
        t3.cm.sim_grid$logLik[i] <- logLik(fit)[[1]]
      } else {
        t3.cm.sim_grid$a_log.est[i] <- NaN
        t3.cm.sim_grid$h_log.est[i] <- NaN
        t3.cm.sim_grid$log_w.est[i] <- NaN
        t3.cm.sim_grid$kr_log.est[i] <- NaN
        
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

CP_t3.cm_fixLV_pe %>% arrange(BIC)

#save(CP_t3.cm_fixLV_pe, file = "CP_t3.cm_fixLV_pe.RData")

###### Best Model #####
(CP_t3.cm_PECs <- CP_t3.cm_fixLV_pe %>% 
   filter(converge == 0) %>% 
   filter(BIC <= min(CP_t3.cm_fixLV_pe$BIC, na.rm = T)+2)) #check for similarity in parameter estimates within 2 AIC of min

min_t3.cm <- CP_t3.cm_PECs %>% filter(BIC <= min(CP_t3.cm_PECs$BIC, na.rm = T)+0.0001) %>% dplyr::select(c(a_log.est:sigma2.est, BIC))

if(nrow(min_t3.cm)>1){
  min_t3.cm <- CP_t3.cm_PECs %>% filter(strtpt %in% sample(filter(CP_t3.cm_PECs, BIC <= min(CP_t3.cm_PECs$BIC, na.rm = T)+0.0001)$strtpt,1))
}

CP_min_t3.cm_fixLV_mod =  mle2(minuslogl = nll.odeint.T3.CM.lv,
                               method= "Nelder-Mead",
                               control = list(maxit = 5000),
                               
                               start = list(a_log = min_t3.cm$a_log.est, 
                                            h_log = min_t3.cm$h_log.est,
                                            log_w = min_t3.cm$log_w.est,
                                            kr_log = min_t3.cm$kr_log.est,
                                            sigma = min_t3.cm$sigma.est,
                                            sigma2 = min_t3.cm$sigma2.est
                               ),
                               
                               data = with(list(N0 = Colp_T0.dens.ml,
                                                Ndead = Colp_T0.dens.ml - Colp_T24.dens.ml,
                                                P = Pred_T0.dens.ml, 
                                                P.end = Pred_T24.dens.ml,
                                                M0 = Para_T0.dens.ml,
                                                M0.end = Para_T24.dens.ml,
                                                Tt = Time.day),
                                           data = filter(colp_para_wrk)),
                               
                               fixed = list(m = mort,
                                            c_log = -15,
                                            
                                            r1_log = coef(CPcomp_pe_mod)[1], 
                                            alpha11 = coef(CPcomp_pe_mod)[2],
                                            alpha12 = coef(CPcomp_pe_mod)[3],
                                            r2_log = coef(CPcomp_pe_mod)[4], 
                                            alpha22 = coef(CPcomp_pe_mod)[5],
                                            alpha21 = coef(CPcomp_pe_mod)[6])
)
#save(CP_min_t3.cm_fixLV_mod, file = "CP_t3.cm_fixLV_mod.Rdata")
Colp_fit13 <- tidy(CP_min_t3.cm_fixLV_mod)

(nrow(CP_t3.cm_PECs)/nrow(CP_t3.cm_fixLV_pe %>% filter(converge == 0))) * 100 # percent of starting combinations within 2 AIC that 
BIC(CP_min_t3.cm_fixLV_mod)

set.seed(1001)
pars.picked = mvrnorm(1000, mu = Colp_fit13$estimate, Sigma = vcov(CP_min_t3.cm_fixLV_mod)) # pick new parameter values by sampling from multivariate normal distribution based on fit
Colp_fit13_pe_confint <- apply(pars.picked,2,quant)

#### BIC Table ####
CP_BIC<-as.data.frame(BIC(CP_min_T2_fixLV_mod,
                          CP_min_t2.a12_fixLV_mod,
                          CP_min_t2.h12_fixLV_mod,
                          CP_min_t2.a12.h12_fixLV_mod,
                          CP_min_t2.vv_fixLV_mod,
                          CP_min_t2.vv.h12_fixLV_mod,
                          CP_min_t2.cm_fixLV_mod,
                          CP_min_t3_fixLV_mod,
                          CP_min_t3.a12_fixLV_mod,
                          CP_min_t3.h12_fixLV_mod,
                          CP_min_t3.a12.h12_fixLV_mod,
                          CP_min_t3.vv_fixLV_mod,
                          CP_min_t3.cm_fixLV_mod))
CP_BIC$model<-seq(1, nrow(CP_BIC), by =1); CP_BIC$dBIC<-round(CP_BIC$BIC-min(CP_BIC$BIC, na.rm = T),1)

CP_BIC %>% arrange(BIC)


#### Parameter Table
pe_list <- list(coef(CP_min_T2_fixLV_mod),
                coef(CP_min_t2.a12_fixLV_mod),
                coef(CP_min_t2.h12_fixLV_mod),
                coef(CP_min_t2.a12.h12_fixLV_mod),
                coef(CP_min_t2.vv_fixLV_mod),
                coef(CP_min_t2.vv.h12_fixLV_mod),
                coef(CP_min_t2.cm_fixLV_mod),
                coef(CP_min_t3_fixLV_mod),
                coef(CP_min_t3.a12_fixLV_mod),
                coef(CP_min_t3.h12_fixLV_mod),
                coef(CP_min_t3.a12.h12_fixLV_mod),
                coef(CP_min_t3.vv_fixLV_mod),
                coef(CP_min_t3.cm_fixLV_mod)
)

pe_table <- as.data.frame(do.call(rbind, lapply(lapply(pe_list, unlist), "[",
                                                unique(unlist(c(sapply(pe_list,names)))))))
pe_table$model<-seq(1, nrow(pe_table), by =1)
pe_table

#### Variance table 
LCI_list <- list(Colp_fit1_pe_confint[1,],
                 Colp_fit2_pe_confint[1,],
                 Colp_fit3_pe_confint[1,],
                 Colp_fit4_pe_confint[1,],
                 Colp_fit5_pe_confint[1,],
                 Colp_fit6_pe_confint[1,],
                 Colp_fit7_pe_confint[1,],
                 Colp_fit8_pe_confint[1,],
                 Colp_fit9_pe_confint[1,],
                 Colp_fit10_pe_confint[1,],
                 Colp_fit11_pe_confint[1,],
                 Colp_fit12_pe_confint[1,])
LCI_table <- as.data.frame(do.call(rbind, lapply(lapply(LCI_list, unlist), "[",
                                                unique(unlist(c(sapply(LCI_list,names)))))))
colnames(LCI_table) <- c("a_log","h_log","sigma","sigma2","a12","h12","log_w","kr_log")
LCI_table$model<-seq(1, nrow(LCI_table), by =1)
LCI_table

UCI_list <- list(Colp_fit1_pe_confint[2,],
                 Colp_fit2_pe_confint[2,],
                 Colp_fit3_pe_confint[2,],
                 Colp_fit4_pe_confint[2,],
                 Colp_fit5_pe_confint[2,],
                 Colp_fit6_pe_confint[2,],
                 Colp_fit7_pe_confint[2,],
                 Colp_fit8_pe_confint[2,],
                 Colp_fit9_pe_confint[2,],
                 Colp_fit10_pe_confint[2,],
                 Colp_fit11_pe_confint[2,],
                 Colp_fit12_pe_confint[2,])
UCI_table <- as.data.frame(do.call(rbind, lapply(lapply(UCI_list, unlist), "[",
                                                 unique(unlist(c(sapply(UCI_list,names)))))))
colnames(UCI_table) <- c("a_log","h_log","sigma","sigma2","a12","h12","log_w","kr_log")
UCI_table$model<-seq(1, nrow(UCI_table), by =1)
UCI_table
