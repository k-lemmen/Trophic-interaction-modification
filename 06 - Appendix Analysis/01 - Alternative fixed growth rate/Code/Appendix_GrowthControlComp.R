library(tidyverse)
library(here)
library(lme4)
library(ggpubr)
library(odeintr)

library(doParallel)
library(foreach)
library(bbmle)

library(readr)
library(readxl)
library(broom)
library(janitor)

#### Functions ####
source(here("02 - Functions", "Function_LVgrowth.R"))

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

##### Data ####
mort <- 0.13333
conversion <- exp(-15)
TIM_data <- loadRData(here("01 - Data", "TIM_data.RData")) # Data

colp_para_wrk <- TIM_data %>% filter(Exp == "CP")
dexio_para_wrk <- TIM_data %>% filter(Exp == "DP")

#### Colpidium ####
G1_colp_growth <- colp_para_wrk %>% 
  filter(Type == "Colp_Con") %>% 
  filter(Date =="23.11.2021" | Date =="28.11.2021" | Date =="2.12.2021" | Date =="4.12.2021" )

G2_colp_growth <- colp_para_wrk %>% 
  filter(Type == "Colp_Con") %>% 
  filter(Date !="23.11.2021" & Date !="28.11.2021" & Date !="2.12.2021" & Date !="4.12.2021" )

##### Grid Search #####
pre_sg <- expand.grid(r1_log = c(-4, -2, -1.5, -1, -0.5, -0.1, 0), 
                      alpha11 = c(0.0001, 0.001, 0.003, 0.005, 0.01, 0.1), 
                      sigma = 0.2)

pre_sg$error <- NA

pre_sg$r1_log.est <- NA
pre_sg$alpha11.est <- NA
pre_sg$sigma.est <- NA

pre_sg$converge <- NA
pre_sg$iteration <- NA
pre_sg$AICc <- NA
pre_sg$logLik <- NA
pre_sg$strtpt <- seq(1:nrow(pre_sg))

sim_grid <- pre_sg #another approach is to create very large search grid and then subset down to 100 starting PECs

###### Group 1 ######
bigstatsr::nb_cores()

my.cluster <- parallel::makeCluster(
  bigstatsr::nb_cores(), 
  type = "FORK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

###### Run Grid ######
G1_colp_PE<-NA # name of file where parallel runs will be save

(t0<-Sys.time()) # timer
mcoptions <- list(preschedule=FALSE)

G1_colp_PE <- foreach(i = 1:nrow(sim_grid),
                          .combine = 'rbind',
                          .options.multicore = mcoptions,
                          .packages = c('bbmle', 'odeintr', 'tidyverse', 'R.utils'),
                          .errorhandling = "pass") %dopar%
  {
    compile_sys("FRS_lvgrowth", FRS.lvgrowth, pars = c("r1", "alpha11"), method = "rk54")
    
    r<-tryCatch({ 
      setTimeLimit(35*60)
      fit = mle2(minuslogl = nll.odeint.lvgrowth,
                 method="BFGS",
                 control = list(maxit = 5000),
                 start = list(r1_log = sim_grid$r1_log[[i]],
                              alpha11 = sim_grid$alpha11[[i]],
                              sigma = sim_grid$sigma[[i]]),
                 
                 data = with(list(N0 = Colp_T0.dens.ml,
                                  Ndead = Colp_T0.dens.ml - Colp_T24.dens.ml,
                                  Tt = Time.day),
                             data = G1_colp_growth)
      ) 
      
      if(fit@details$convergence==0) {
        sim_grid$r1_log.est[i] <- fit@coef[[1]]
        sim_grid$alpha11.est[i] <- fit@coef[[2]]
        
        sim_grid$sigma.est [i] <- fit@coef[[3]]
        
        sim_grid$iteration[i] <- fit@details[[3]][[2]]
        sim_grid$converge[i] <- fit@details[[4]]
        sim_grid$AICc[i] <- AICc(fit)
        sim_grid$logLik[i] <- logLik(fit)[[1]]
      } else {
        sim_grid$r1_log.est[i] <- NaN
        sim_grid$alpha11.est[i] <- NaN
        
        sim_grid$converge[i] <- NaN
        sim_grid$iteration[i] <- NaN
        sim_grid$AICc[i] <- NaN
        sim_grid$logLik[i] <- NaN
      }
      
      #Remove fit for next iteration
      rm(fit)
    }, error=function(e) list(message = e$message)) # to print error in console use {cat("ERROR :", paste(conditionMessage(e), "Combo", (i)), "\n")}
    
    if(length(r) > 0){
      sim_grid$error[i] <- r 
    } else {
      sim_grid$error[i] <- NaN
    }  
    
    return(sim_grid[i,])
    
  }

parallel::stopCluster(cl = my.cluster)
stopImplicitCluster()

G1_colp_PE %>%arrange(AICc)

(G1_colp_best <- G1_colp_PE %>% filter(AICc == min(G1_colp_PE$AICc, na.rm = T)))

G1_colp_r_best.mod = mle2(minuslogl = nll.odeint.lvgrowth,
                       method="BFGS",
                       control = list(maxit = 5000),
                       start = list(r1_log = G1_colp_best$r1_log.est,
                                    alpha11 = G1_colp_best$alpha11.est,
                                    sigma = G1_colp_best$sigma.est),
                       
                       data = with(list(N0 = Colp_T0.dens.ml,
                                        Ndead = Colp_T0.dens.ml - Colp_T24.dens.ml,
                                        Tt = Time.day),
                                   data = G1_colp_growth)
)
G1_colp_mod_tidy <- tidy(G1_colp_r_best.mod)

### 95% Confidence intervals
library(MASS); set.seed(1001)
nresamp=1000
quant = function(col) quantile(col, c(0.025,0.975)) # 95% percentiles 
pars.picked = mvrnorm(1000, mu = G1_colp_mod_tidy$estimate, Sigma = vcov(G1_colp_r_best.mod)) # pick new parameter values by sampling from multivariate normal distribution based on fit
(G1_colp_pe_confint <- apply(pars.picked,2,quant))


###### Group 2 ######
bigstatsr::nb_cores()

my.cluster <- parallel::makeCluster(
  bigstatsr::nb_cores(), 
  type = "FORK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

###### Run Grid ######
G2_colp_PE<-NA # name of file where parallel runs will be save

(t0<-Sys.time()) # timer
mcoptions <- list(preschedule=FALSE)

G2_colp_PE <- foreach(i = 1:nrow(sim_grid),
                      .combine = 'rbind',
                      .options.multicore = mcoptions,
                      .packages = c('bbmle', 'odeintr', 'tidyverse', 'R.utils'),
                      .errorhandling = "pass") %dopar%
  {
    compile_sys("FRS_lvgrowth", FRS.lvgrowth, pars = c("r1", "alpha11"), method = "rk54")
    
    r<-tryCatch({ 
      setTimeLimit(35*60)
      fit = mle2(minuslogl = nll.odeint.lvgrowth,
                 method="BFGS",
                 control = list(maxit = 5000),
                 start = list(r1_log = sim_grid$r1_log[[i]],
                              alpha11 = sim_grid$alpha11[[i]],
                              sigma = sim_grid$sigma[[i]]),
                 
                 data = with(list(N0 = Colp_T0.dens.ml,
                                  Ndead = Colp_T0.dens.ml - Colp_T24.dens.ml,
                                  Tt = Time.day),
                             data = G2_colp_growth)
      ) 
      
      if(fit@details$convergence==0) {
        sim_grid$r1_log.est[i] <- fit@coef[[1]]
        sim_grid$alpha11.est[i] <- fit@coef[[2]]
        
        sim_grid$sigma.est [i] <- fit@coef[[3]]
        
        sim_grid$iteration[i] <- fit@details[[3]][[2]]
        sim_grid$converge[i] <- fit@details[[4]]
        sim_grid$AICc[i] <- AICc(fit)
        sim_grid$logLik[i] <- logLik(fit)[[1]]
      } else {
        sim_grid$r1_log.est[i] <- NaN
        sim_grid$alpha11.est[i] <- NaN
        
        sim_grid$converge[i] <- NaN
        sim_grid$iteration[i] <- NaN
        sim_grid$AICc[i] <- NaN
        sim_grid$logLik[i] <- NaN
      }
      
      #Remove fit for next iteration
      rm(fit)
    }, error=function(e) list(message = e$message)) # to print error in console use {cat("ERROR :", paste(conditionMessage(e), "Combo", (i)), "\n")}
    
    if(length(r) > 0){
      sim_grid$error[i] <- r 
    } else {
      sim_grid$error[i] <- NaN
    }  
    
    return(sim_grid[i,])
    
  }

parallel::stopCluster(cl = my.cluster)
stopImplicitCluster()

G2_colp_PE %>%arrange(AICc)

(G2_colp_best <- G2_colp_PE %>% filter(AICc == min(G2_colp_PE$AICc, na.rm = T)))

G2_colp_r_best.mod = mle2(minuslogl = nll.odeint.lvgrowth,
                          method="BFGS",
                          control = list(maxit = 5000),
                          start = list(r1_log = G2_colp_best$r1_log.est,
                                       alpha11 = G2_colp_best$alpha11.est,
                                       sigma = G2_colp_best$sigma.est),
                          
                          data = with(list(N0 = Colp_T0.dens.ml,
                                           Ndead = Colp_T0.dens.ml - Colp_T24.dens.ml,
                                           Tt = Time.day),
                                      data = G2_colp_growth)
)
G2_colp_mod_tidy <- tidy(G2_colp_r_best.mod)

### 95% Confidence intervals
library(MASS); set.seed(1001)
nresamp=1000
quant = function(col) quantile(col, c(0.025,0.975)) # 95% percentiles 
pars.picked = mvrnorm(1000, mu = G2_colp_mod_tidy$estimate, Sigma = vcov(G2_colp_r_best.mod)) # pick new parameter values by sampling from multivariate normal distribution based on fit
(G2_colp_pe_confint <- apply(pars.picked,2,quant))

##### Comparison #####
Group <- c("G1", "G1", "G1", "G2", "G2", "G2")
Term <- c(G1_colp_mod_tidy$term, G2_colp_mod_tidy$term)
Estimate <- c(G1_colp_mod_tidy$estimate, G2_colp_mod_tidy$estimate)
LCI <-c(G1_colp_pe_confint[1,], G2_colp_pe_confint[1,])
UCI <-c(G1_colp_pe_confint[2,], G2_colp_pe_confint[2,])

colp_comp <- as.data.frame(cbind(Group, Term, Estimate, LCI, UCI))
colp_comparison <- type_convert(colp_comp)
colp_comparison <- type_convert(colp_comp)

Colp_growth_comp <- colp_comparison %>% filter(Term != "sigma") %>% 
  ggplot() +
  geom_point(aes(Group, Estimate)) +
  geom_errorbar(aes(x=Group, ymin=LCI, ymax=UCI), width=.1) +
  facet_wrap(Term ~ ., scales = "free") +
  labs(title = "Colpidium Block 1 vs Block 2", y="Parameter Estimate") +
  theme_bw () +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), strip.text = element_text(size=14),plot.title = element_text(size=16, hjust = 0.5),
        panel.spacing = unit(5, "mm"), panel.border = element_rect(colour = "black", fill=NA, size=2), strip.background = element_rect(colour="white", fill="white"),
        legend.position="none")

#row.names(colp_comparison) <- c("colp_logr_B1","colp_alpha11_B1","sigma_B1","colp_logr_B2","colp_alpha11_B2","sigma_B2")
#save(colp_comparison, file = "CP_growth_pe.RData")

#### Dexiostoma ####
G1_dexio_growth <- dexio_para_wrk %>% 
  filter(Type == "Dexio_Con") %>% 
  filter(Date =="7.12.2021" | Date =="8.12.2021" | Date =="10.12.2021" | Date =="11.12.2021" | Date =="13.12.2021" | Date =="15.12.2021")

G2_dexio_growth <- dexio_para_wrk %>% 
  filter(Type == "Dexio_Con") %>% 
  filter(Date !="7.12.2021" & Date !="8.12.2021" & Date !="10.12.2021" & Date !="11.12.2021" & Date !="13.12.2021" & Date !="15.12.2021")

##### Grid Search #####
pre_sg <- expand.grid(r1_log = seq(-8.5,-7.5, by =0.1), 
                      alpha11 = c(0.00001, 0.0001, 0.0003, 0.0005, 0.00075,0.001, 0.005, 0.01, 0.1), 
                      sigma = 0.2)

pre_sg$error <- NA

pre_sg$r1_log.est <- NA
pre_sg$alpha11.est <- NA
pre_sg$sigma.est <- NA

pre_sg$converge <- NA
pre_sg$iteration <- NA
pre_sg$AICc <- NA
pre_sg$logLik <- NA
pre_sg$strtpt <- seq(1:nrow(pre_sg))

sim_grid <- pre_sg #another approach is to create very large search grid and then subset down to 100 starting PECs

###### Group 1 ######
bigstatsr::nb_cores()

my.cluster <- parallel::makeCluster(
  bigstatsr::nb_cores(), 
  type = "FORK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

###### Run Grid ######
G1_dexio_PE<-NA # name of file where parallel runs will be save

(t0<-Sys.time()) # timer
mcoptions <- list(preschedule=FALSE)

G1_dexio_PE <- foreach(i = 1:nrow(sim_grid),
                       .combine = 'rbind',
                       .options.multicore = mcoptions,
                       .packages = c('bbmle', 'odeintr', 'tidyverse', 'R.utils'),
                       .errorhandling = "pass") %dopar%
  {
    compile_sys("FRS_lvgrowth", FRS.lvgrowth, pars = c("r1", "alpha11"), method = "rk54")
    
    r<-tryCatch({ 
      setTimeLimit(35*60)
      fit = mle2(minuslogl = nll.odeint.lvgrowth,
                 method="BFGS",
                 control = list(maxit = 5000),
                 start = list(r1_log = sim_grid$r1_log[[i]],
                              alpha11 = sim_grid$alpha11[[i]],
                              sigma = sim_grid$sigma[[i]]),
                 
                 data = with(list(N0 = Dexio_T0.dens.ml,
                                  Ndead = Dexio_T0.dens.ml - Dexio_T24.dens.ml,
                                  Tt = Time.day),
                             data = G1_dexio_growth)
      ) 
      
      if(fit@details$convergence==0) {
        sim_grid$r1_log.est[i] <- fit@coef[[1]]
        sim_grid$alpha11.est[i] <- fit@coef[[2]]
        
        sim_grid$sigma.est [i] <- fit@coef[[3]]
        
        sim_grid$iteration[i] <- fit@details[[3]][[2]]
        sim_grid$converge[i] <- fit@details[[4]]
        sim_grid$AICc[i] <- AICc(fit)
        sim_grid$logLik[i] <- logLik(fit)[[1]]
      } else {
        sim_grid$r1_log.est[i] <- NaN
        sim_grid$alpha11.est[i] <- NaN
        
        sim_grid$sigma.est [i] <- NaN
        
        sim_grid$converge[i] <- NaN
        sim_grid$iteration[i] <- NaN
        sim_grid$AICc[i] <- NaN
        sim_grid$logLik[i] <- NaN
      }
      
      #Remove fit for next iteration
      rm(fit)
    }, error=function(e) list(message = e$message)) # to print error in console use {cat("ERROR :", paste(conditionMessage(e), "Combo", (i)), "\n")}
    
    if(length(r) > 0){
      sim_grid$error[i] <- r 
    } else {
      sim_grid$error[i] <- NaN
    }  
    
    return(sim_grid[i,])
    
  }

parallel::stopCluster(cl = my.cluster)
stopImplicitCluster()

G1_dexio_PE %>%arrange(AICc)

(G1_dexio_best <- G1_dexio_PE %>% filter(AICc == min(G1_dexio_PE$AICc, na.rm = T)))

G1_dexio_r_best.mod = mle2(minuslogl = nll.odeint.lvgrowth,
                           method="BFGS",
                           control = list(maxit = 5000),
                           start = list(r1_log = G1_dexio_best$r1_log.est,
                                        alpha11 = G1_dexio_best$alpha11.est,
                                        sigma = G1_dexio_best$sigma.est),
                           
                           data = with(list(N0 = Dexio_T0.dens.ml,
                                            Ndead = Dexio_T0.dens.ml - Dexio_T24.dens.ml,
                                            Tt = Time.day),
                                       data = G1_dexio_growth)
)
G1_dexio_mod_tidy <- tidy(G1_dexio_r_best.mod)

### 95% Confidence intervals
library(MASS); set.seed(1001)
nresamp=1000
quant = function(col) quantile(col, c(0.025,0.975)) # 95% percentiles 
pars.picked = mvrnorm(1000, mu = G1_dexio_mod_tidy$estimate, Sigma = vcov(G1_dexio_r_best.mod)) # pick new parameter values by sampling from multivariate normal distribution based on fit
(G1_dexio_pe_confint <- apply(pars.picked,2,quant))


###### Group 2 ######
bigstatsr::nb_cores()

my.cluster <- parallel::makeCluster(
  bigstatsr::nb_cores(), 
  type = "FORK"
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

###### Run Grid ######
G2_dexio_PE<-NA # name of file where parallel runs will be save

(t0<-Sys.time()) # timer
mcoptions <- list(preschedule=FALSE)

G2_dexio_PE <- foreach(i = 1:nrow(sim_grid),
                       .combine = 'rbind',
                       .options.multicore = mcoptions,
                       .packages = c('bbmle', 'odeintr', 'tidyverse', 'R.utils'),
                       .errorhandling = "pass") %dopar%
  {
    compile_sys("FRS_lvgrowth", FRS.lvgrowth, pars = c("r1", "alpha11"), method = "rk54")
    
    r<-tryCatch({ 
      setTimeLimit(35*60)
      fit = mle2(minuslogl = nll.odeint.lvgrowth,
                 method="BFGS",
                 control = list(maxit = 5000),
                 start = list(r1_log = sim_grid$r1_log[[i]],
                              alpha11 = sim_grid$alpha11[[i]],
                              sigma = sim_grid$sigma[[i]]),
                 
                 data = with(list(N0 = Dexio_T0.dens.ml,
                                  Ndead = Dexio_T0.dens.ml - Dexio_T24.dens.ml,
                                  Tt = Time.day),
                             data = G2_dexio_growth)
      ) 
      
      if(fit@details$convergence==0) {
        sim_grid$r1_log.est[i] <- fit@coef[[1]]
        sim_grid$alpha11.est[i] <- fit@coef[[2]]
        
        sim_grid$sigma.est [i] <- fit@coef[[3]]
        
        sim_grid$iteration[i] <- fit@details[[3]][[2]]
        sim_grid$converge[i] <- fit@details[[4]]
        sim_grid$AICc[i] <- AICc(fit)
        sim_grid$logLik[i] <- logLik(fit)[[1]]
      } else {
        sim_grid$r1_log.est[i] <- NaN
        sim_grid$alpha11.est[i] <- NaN
        
        sim_grid$converge[i] <- NaN
        sim_grid$iteration[i] <- NaN
        sim_grid$AICc[i] <- NaN
        sim_grid$logLik[i] <- NaN
      }
      
      #Remove fit for next iteration
      rm(fit)
    }, error=function(e) list(message = e$message)) # to print error in console use {cat("ERROR :", paste(conditionMessage(e), "Combo", (i)), "\n")}
    
    if(length(r) > 0){
      sim_grid$error[i] <- r 
    } else {
      sim_grid$error[i] <- NaN
    }  
    
    return(sim_grid[i,])
    
  }

parallel::stopCluster(cl = my.cluster)
stopImplicitCluster()

G2_dexio_PE %>%arrange(AICc)

(G2_dexio_best <- G2_dexio_PE %>% filter(AICc == min(G2_dexio_PE$AICc, na.rm = T)))

G2_dexio_r_best.mod = mle2(minuslogl = nll.odeint.lvgrowth,
                           method="BFGS",
                           control = list(maxit = 5000),
                           start = list(r1_log = G2_dexio_best$r1_log.est,
                                        alpha11 = G2_dexio_best$alpha11.est,
                                        sigma = G2_dexio_best$sigma.est),
                           
                           data = with(list(N0 = Dexio_T0.dens.ml,
                                            Ndead = Dexio_T0.dens.ml - Dexio_T24.dens.ml,
                                            Tt = Time.day),
                                       data = G2_dexio_growth)
)
G2_dexio_mod_tidy <- tidy(G2_dexio_r_best.mod)

### 95% Confidence intervals
library(MASS); set.seed(1001)
nresamp=1000
quant = function(col) quantile(col, c(0.025,0.975)) # 95% percentiles 
pars.picked = mvrnorm(1000, mu = G2_dexio_mod_tidy$estimate, Sigma = vcov(G2_dexio_r_best.mod)) # pick new parameter values by sampling from multivariate normal distribution based on fit
(G2_dexio_pe_confint <- apply(pars.picked,2,quant))

##### Comparison #####
Group <- c("G1", "G1", "G1", "G2", "G2", "G2")
Term <- c(G1_dexio_mod_tidy$term, G2_dexio_mod_tidy$term)
D_Estimate <- c(G1_dexio_mod_tidy$estimate, G2_dexio_mod_tidy$estimate)
D_LCI <-c(G1_dexio_pe_confint[1,], G2_dexio_pe_confint[1,])
D_UCI <-c(G1_dexio_pe_confint[2,], G2_dexio_pe_confint[2,])

dexio_comp <- as.data.frame(cbind(Group, Term, D_Estimate, D_LCI, D_UCI))
dexio_comparison <- type_convert(dexio_comp)

Dexio_growth_comp <- dexio_comparison %>% filter(Term != "sigma") %>% 
  ggplot() +
  geom_point(aes(Group, D_Estimate)) +
  geom_errorbar(aes(x=Group, ymin=D_LCI, ymax=D_UCI), width=.1) +
  facet_wrap(Term ~ ., scales = "free") +
  labs(title = "Dexiostoma Block 1 vs Block 2", y="Parameter Estimate") +
  theme_bw () +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), strip.text = element_text(size=14),plot.title = element_text(size=16, hjust = 0.5),
        #axis.text.x = element_text(angle = 40, vjust = 1, hjust=1),
        panel.spacing = unit(5, "mm"), panel.border = element_rect(colour = "black", fill=NA, size=2), strip.background = element_rect(colour="white", fill="white"),
        legend.position="none")

#row.names(dexio_comparison) <- c("dexio_logr_B1","dexio_alpha11_B1","sigma_B1","dexio_logr_B2","dexio_alpha11_B2","sigma_B2")
#save(dexio_comparison, file = "DP_growth_pe.RData")

#### Figure ####
(Growth_Comp <- ggarrange(Dexio_growth_comp, Colp_growth_comp,
                                  ncol = 2, nrow = 1,
                                  label.x = 0, labels = c("a","b"), font.label = list(size = 28, face = "bold"),
                                  common.legend=F) 
)  

png(filename="Prey_Spp_Growth_Control_Comp.png", height=500, width=1500, units="px", res=75)
Growth_Comp
dev.off()
