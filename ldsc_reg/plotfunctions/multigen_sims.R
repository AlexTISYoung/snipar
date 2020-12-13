library(tidyverse)
library(readxl)
theme_set(theme_minimal())
ldsc_path <- "C:/Users/Hariharan/Documents/git_repos/SNIPar/ldsc_reg/"
source(paste0(ldsc_path, "plotfunctions/helperfuncs.R"))


# reading in mutligen sim data
log_files <- Sys.glob(paste0(ldsc_path, "multigen_sim/multisims/[0-9]*.log"))
multigen <- readlogfiles(log_files)

r_sd <- sd(multigen$r)
N <- length(multigen$r)
r_se <- r_sd#/sqrt(N)

# plotting SEs
multigen %>%
    ggplot() +
    geom_hline(aes(color = "True Standard Error", 
                   yintercept = r_se)) +
    geom_point(aes(x = 1:nrow(multigen), y = r_invhse,
                   color = "Inverse Hessian")) +
    geom_smooth(aes(x = 1:nrow(multigen), y = r_invhse,
                    color = "Inverse Hessian"), 
                method = lm, se= FALSE) +
    geom_point(aes(x = 1:nrow(multigen), y = r_jkse,
                   color = "Block Jack Knife")) +
    geom_smooth(aes(x = 1:nrow(multigen), y = r_jkse,
                    color = "Block Jack Knife"), 
                method = lm, se= FALSE) +
    scale_color_brewer(palette = "Set1") +
    labs(color = "",
         x = "Run", y = "",
         title = "Jack Knife vs Inverse Hessian")


ggsave(paste0(ldsc_path, "multigen_sim/multisims/comparing_se.png"),
       height = 6,
       width = 9) 
