library(tidyverse)
library(latex2exp)

source("ldsc_reg/plotfunctions/helperfuncs.R")


log_files <- Sys.glob("ldsc_reg/multigen_sim/05012020/[0-9]*.log")
simdf <- readlogfiles(log_files)

true_sd <- sd(simdf$r)
true_se <- true_sd/sqrt(nrow(simdf))

simdf %>% 
    ggplot() +
    geom_point(aes(traitcode, r_invhse, color = "Inverse Hessian")) +
    geom_point(aes(traitcode, r_jkse, color = "Block Jack Knife")) +
    geom_hline(yintercept = true_sd) +
    geom_text(aes(x = 50, y = true_sd + 0.0006, label = "True SD")) +
    geom_hline(yintercept = true_se) +
    geom_text(aes(x = 50, y = true_se + 0.0006, label = "True SE")) +
    scale_color_brewer(palette = "Dark2") +
    theme_minimal() +
    labs(x = "Run Number", y = "Standard Error Estimate",
         color = "")
ggsave("ldsc_reg/multigen_sim/05012020/comparing_se.png",
       height = 6, width = 9)   


simdf %>% 
    ggplot() +
    geom_point(aes(traitcode, r, color = "Estimates")) +
    guides(color = F) +
    theme_minimal() +
    scale_color_brewer(palette = "Dark2") +
    labs(x = "Run", y = TeX("$\\textit{r_g}$ Direct + Population"))

ggsave("ldsc_reg/multigen_sim/05012020/rg_estimate.png",
       height = 6, width = 9)   