#### Effect of modifier on Space clearance constant ####

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
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

quant = function(col) quantile(col, c(0.025,0.975)) # 95% percentiles

#### Data ####
nresamp=1000
mod_density <- seq(1,500, length.out = 1000)

TIM_data <- loadRData(here("01 - Data", "TIM_data.RData")) # Data

##### Colpidium: Type 2 - linear a #####
colp_para_wrk <- TIM_data %>% filter(Exp == "CP")
CP_best_mod <- loadRData(here("05 - Output", "Model Fits", "Colpidium", "CP_t2.a12_fixLV_mod.Rdata"))

##### Dexiostoma: Generalized - linear a ##### 
dexio_para_wrk <- TIM_data %>% filter(Exp == "DP")
DP_best_mod <- loadRData(here("05 - Output", "Model Fits", "Dexiostoma", "DP_t3.a12_fixLV_mod.Rdata"))

#### Colpidium - effect of modifier on space clearance constant ####
Colp_PE <- tidy(CP_best_mod)

set.seed(1001)
pars.picked = mvrnorm(1000, mu = Colp_PE$estimate, Sigma = vcov(CP_best_mod)) # pick new parameter values by sampling from multivariate normal distribution based on fit
Colp_pe_confint <- apply(pars.picked,2,quant)

realized_SCR <- exp(Colp_PE[Colp_PE$term == "a_log",]$estimate) + Colp_PE[Colp_PE$term == "a12",]$estimate * mod_density
LCI_rSCR <- exp(Colp_pe_confint[1,1]) + Colp_pe_confint[1,2] * mod_density
LCI_rSCR[LCI_rSCR< 0] <- 0

UCI_rSCR <- exp(Colp_pe_confint[2,1]) + Colp_pe_confint[2,2] * mod_density
  
colp_real_SCR <- as.data.frame(cbind(mod_density,realized_SCR, LCI_rSCR, UCI_rSCR))

colp_rscr_viz <- colp_real_SCR %>% ggplot() +
  ggtitle("Colpidium") +
  xlab("Modifier Density (ind/mL)") + ylab(" ") +
  ylab(bquote('Realized Space Clearance Constant (a)')) +
  geom_ribbon(aes(x = mod_density, ymin=LCI_rSCR, ymax=UCI_rSCR), fill = "red", linetype = 1, alpha=0.2) +
  geom_line(aes(x=mod_density, y=realized_SCR), lty=1, linewidth = 1) +
  theme_bw() +
  coord_cartesian(ylim=c(0, 0.31)) +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=18), strip.text = element_text(size=16), plot.title = element_text(size=25),
        plot.margin = margin(3, 3, 5, 10, "mm"), 
        legend.position="none")

#### Dexiostoma - effect of modifier on space clearance constant ####
Dexio_PE <- tidy(DP_best_mod)

set.seed(1001)
pars.picked = mvrnorm(1000, mu = Dexio_PE$estimate, Sigma = vcov(DP_best_mod)) # pick new parameter values by sampling from multivariate normal distribution based on fit
Dexio_pe_confint <- apply(pars.picked,2,quant)

realized_SCR <- exp(Dexio_PE[Dexio_PE$term == "a_log",]$estimate) + Dexio_PE[Dexio_PE$term == "a12",]$estimate * mod_density
LCI_rSCR <- exp(Dexio_pe_confint[1,1]) + Dexio_pe_confint[1,2] * mod_density
LCI_rSCR[LCI_rSCR< 0] <- 0

UCI_rSCR <- exp(Dexio_pe_confint[2,1]) + Dexio_pe_confint[2,2] * mod_density

dexio_real_SCR <- as.data.frame(cbind(mod_density,realized_SCR, LCI_rSCR, UCI_rSCR))

dexio_rscr_viz <- dexio_real_SCR %>% ggplot() +
  ggtitle("Dexiostoma") +
  xlab("Modifier Density (ind/mL)") + ylab(" ") +
  ylab(bquote('Realized Space Clearance Constant (a)')) +
  geom_ribbon(aes(x = mod_density, ymin=LCI_rSCR, ymax=UCI_rSCR), fill = "red", linetype = 1, alpha=0.2) +
  geom_line(aes(x=mod_density, y=realized_SCR), lty=1, linewidth = 1) +
  coord_cartesian(ylim=c(0, 0.31)) +
  theme_bw() +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=18), strip.text = element_text(size=16), plot.title = element_text(size=25),
        plot.margin = margin(3, 3, 5, 10, "mm"), 
        legend.position="none")  

ggarrange(colp_rscr_viz, dexio_rscr_viz, ncol = 2, nrow = 1)
