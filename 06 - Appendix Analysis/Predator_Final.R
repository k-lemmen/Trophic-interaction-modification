# Model comparisons #
library(plyr)
library(tidyverse)
library(here)
library(ggpubr)


#### Functions ####
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#### data ####
TIM_data <- loadRData(here("01 - Data", "TIM_data.RData")) # Data

#### plot ####
plot_data <- TIM_data %>% filter(Pred_YN == "Y") %>% mutate_if(is.character, as.factor)
plot_data$Trt <- fct_recode(plot_data$Type, "Prey Only" = "Colp_FR", "Prey Only" = "Dexio_FR", "Modifier Only" = "Para_FR", "Prey and Modifier" = "MSFR")

exp_labels <- c("DP" = "Dexiostoma-Paramecium", "CP" = "Colpidium-Paramecium")

(Colp_pred <- plot_data %>% filter(Exp == "CP") %>% 
    ggplot(aes(x= as.factor(Trt), y = Pred_T24.dens.ml)) +
    stat_summary(fun.data = "mean_cl_boot",  width = .1, size = 1) +
    geom_jitter(size = 1, alpha = 0.25, width = 0.2, height = 0) +
  facet_wrap(Exp ~ ., labeller=labeller(Exp = exp_labels), nrow = 2) +
  geom_abline(slope=0, intercept=5, lty=3) +
  xlab("\nTreatment") +
  ylab(bquote('Final Predator Density (ind ml'^-1*')')) +
  theme_bw() +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=18), strip.text = element_text(size=16), plot.title = element_text(size=15),
        panel.spacing = unit(10, "mm"), panel.border = element_rect(colour = "black", fill=NA, size=2), strip.background = element_rect(colour="white", fill="white"),
        plot.margin = margin(3, 3, 5, 10, "mm"), legend.position="none")
)

(Dexio_pred <- plot_data %>% filter(Exp == "DP") %>% 
    ggplot(aes(x= as.factor(Trt), y = Pred_T24.dens.ml)) +
    stat_summary(fun.data = "mean_cl_boot",  width = .1, size = 1) +
    geom_jitter(size = 1, alpha = 0.25, width = 0.2, height = 0) +
    facet_wrap(Exp ~ ., labeller=labeller(Exp = exp_labels), nrow = 2) +
    geom_abline(slope=0, intercept=5, lty=3) +
    xlab("\nTreatment") +
    ylab(bquote('Final Predator Density (ind ml'^-1*')')) +
    theme_bw() +
    theme(axis.text=element_text(size=16), axis.title=element_text(size=18), strip.text = element_text(size=16), plot.title = element_text(size=15),
          panel.spacing = unit(10, "mm"), panel.border = element_rect(colour = "black", fill=NA, size=2), strip.background = element_rect(colour="white", fill="white"),
          plot.margin = margin(3, 3, 5, 10, "mm"), legend.position="none")
)

FigS5_1 <- ggarrange(Colp_pred, Dexio_pred,
                       ncol = 1, nrow = 2, 
                       label.x = 0, labels = c("a","b"), font.label = list(size = 28, face = "bold"),
                       common.legend=F) 


tiff(filename="TIM_FigureS5-1.tiff", height=2000, width=2000, units="px", res=200, compression="lzw")
FigS5_1
dev.off()
