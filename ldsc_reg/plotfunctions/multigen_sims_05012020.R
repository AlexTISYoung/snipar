library(tidyverse)
library(latex2exp)

source("ldsc_reg/plotfunctions/helperfuncs.R")


# Reading in Data ---------------------------------------------------------


log_files_avgparental <- Sys.glob("ldsc_reg/multigen_sim/05012020/*_avgparental.log")
avg_parental <- readlogfiles(log_files_avgparental)

log_files_pop <- Sys.glob("ldsc_reg/multigen_sim/05012020/*_pop.log")
pop <- readlogfiles(log_files_pop)

ldsc_estimates <- read_csv("ldsc_reg/multigen_sim/05012020/ldsc_estimates.csv")

true_sd_avgpar_snipar <- sd(avg_parental$r)
true_sd_pop_snipar <- sd(pop$r)

true_sd_avgpar_ldsc <- ldsc_estimates %>% 
    filter(effect_estimated == "Average Parental") %>% 
    pull(estimate) %>% 
    sd
true_sd_pop_ldsc <- ldsc_estimates %>% 
    filter(effect_estimated == "Population") %>% 
    pull(estimate) %>% 
    sd


# Plotting ----------------------------------------------------------------

# avg parental
ggplot() +
    geom_point(data = avg_parental, 
               aes(traitcode, r_invhse, color = "Inverse Hessian")) +
    geom_point(data = avg_parental, 
               aes(traitcode, r_jkse, color = "Block Jack Knife")) +
    geom_point(data = ldsc_estimates %>% filter(effect_estimated == "Average Parental"), 
               aes(run, se, color = "LDSC Block Jack Knife")) +
    geom_hline(yintercept = true_sd_avgpar_snipar) +
    geom_text(aes(x = 50, y = true_sd_avgpar_snipar - 0.005, label = "True SD - SNIPar")) +
    geom_hline(yintercept = true_sd_avgpar_ldsc) +
    geom_text(aes(x = 50, y = true_sd_avgpar_ldsc + 0.005, label = "True SD - LDSC")) +
    scale_color_brewer(palette = "Dark2") +
    theme_minimal() +
    labs(x = "Run Number", y = "Standard Error Estimate",
         color = "")

ggsave("ldsc_reg/multigen_sim/05012020/se_avgparental.png",
       height = 6, width = 9)   

ggplot() +
    geom_point(data = pop, 
               aes(traitcode, r_invhse, color = "Inverse Hessian")) +
    geom_point(data = pop, 
               aes(traitcode, r_jkse, color = "Block Jack Knife")) +
    geom_point(data = ldsc_estimates %>% filter(effect_estimated == "Population"), 
               aes(run, se, color = "LDSC Block Jack Knife")) +
    geom_hline(yintercept = true_sd_pop_snipar) +
    geom_text(aes(x = 50, y = true_sd_pop_snipar + 0.002, label = "True SD - SNIPar")) +
    geom_hline(yintercept = true_sd_pop_ldsc) +
    geom_text(aes(x = 50, y = true_sd_pop_ldsc - 0.002, label = "True SD - LDSC")) +
    scale_color_brewer(palette = "Dark2") +
    theme_minimal() +
    labs(x = "Run Number", y = "Standard Error Estimate",
         color = "")

ggsave("ldsc_reg/multigen_sim/05012020/se_pop.png",
       height = 6, width = 9)   

# Comparing estimates
ggplot() +
    geom_point(data = pop, 
               aes(traitcode, r, color = "SNIPar")) +
    geom_point(data = ldsc_estimates %>% filter(effect_estimated == "Population"), 
               aes(run, estimate, color = "LDSC")) +
    geom_hline(yintercept = 1.0) +
    scale_color_brewer(palette = "Dark2") +
    theme_minimal() +
    labs(x = "Run Number", y = "Estimate",
         color = "")

ggsave("ldsc_reg/multigen_sim/05012020/estimates_pop.png",
       height = 6, width = 9)  

ggplot() +
    geom_point(data = avg_parental, 
               aes(traitcode, r, color = "SNIPar")) +
    geom_point(data = ldsc_estimates %>% filter(effect_estimated == "Average Parental"), 
               aes(run, estimate, color = "LDSC")) +
    geom_hline(yintercept = 0.5) +
    scale_color_brewer(palette = "Dark2") +
    theme_minimal() +
    labs(x = "Run Number", y = "Estimate",
         color = "")

ggsave("ldsc_reg/multigen_sim/05012020/estimates_avgparental.png",
       height = 6, width = 9)   
