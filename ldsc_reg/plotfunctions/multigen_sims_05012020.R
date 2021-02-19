library(tidyverse)
library(latex2exp)

source("ldsc_reg/plotfunctions/helperfuncs.R")


# Reading in Data ---------------------------------------------------------


log_files_avgparental <- Sys.glob("ldsc_reg/multigen_sim/05012020/*_avgparental.log")
avg_parental <- readlogfiles(log_files_avgparental)

log_files_pop <- Sys.glob("ldsc_reg/multigen_sim/05012020/*_pop.log")
pop <- readlogfiles(log_files_pop)

ldsc_estimates <- read_csv("ldsc_reg/multigen_sim/05012020/ldsc_estimates.csv")

true_sd_r_avgpar_snipar <- sd(avg_parental$r)
true_sd_r_pop_snipar <- sd(pop$r)
true_sd_h21_avgpar_snipar <- sd(avg_parental$v1)
true_sd_h21_pop_snipar <- sd(pop$v1)
true_sd_h22_avgpar_snipar <- sd(avg_parental$v2)
true_sd_h22_pop_snipar <- sd(pop$v2)

true_sd_r_avgpar_ldsc <- ldsc_estimates %>% 
    filter(effect_estimated == "Average Parental") %>% 
    pull(r) %>% 
    sd
true_sd_h21_avgpar_ldsc <- ldsc_estimates %>% 
    filter(effect_estimated == "Average Parental") %>% 
    pull(h2_dir) %>% 
    sd
true_sd_h22_avgpar_ldsc <- ldsc_estimates %>% 
    filter(effect_estimated == "Average Parental") %>% 
    pull(h2_second) %>% 
    sd
true_sd_r_pop_ldsc <- ldsc_estimates %>% 
    filter(effect_estimated == "Population") %>% 
    pull(r) %>% 
    sd
true_sd_h21_pop_ldsc <- ldsc_estimates %>% 
    filter(effect_estimated == "Population") %>% 
    pull(h2_dir) %>% 
    sd
true_sd_h22_pop_ldsc <- ldsc_estimates %>% 
    filter(effect_estimated == "Population") %>% 
    pull(h2_second) %>% 
    sd


# Plotting ----------------------------------------------------------------

# avg parental
ggplot() +
    geom_density(data = avg_parental, 
               aes(r_invhse, color = "Inverse Hessian")) +
    geom_density(data = avg_parental, 
               aes(r_jkse, color = "Block Jack Knife")) +
    geom_density(data = ldsc_estimates %>% filter(effect_estimated == "Average Parental"), 
               aes(r_se, color = "LDSC Block Jack Knife")) +
    geom_vline(aes(xintercept = true_sd_r_avgpar_snipar, linetype = "True SD SNIPar")) +
    geom_vline(aes(xintercept = true_sd_r_avgpar_ldsc, linetype = "True SD LDSC")) +
    scale_color_brewer(palette = "Dark2") +
    theme_minimal() +
    labs(color = "",
         x = "SE",
         title = "Genetic Correlation - Average Parental")

ggsave("ldsc_reg/multigen_sim/05012020/se_avgparental.png",
       height = 6, width = 9)   

ggplot() +
    geom_density(data = pop, 
               aes(r_invhse, color = "Inverse Hessian")) +
    geom_density(data = pop, 
               aes(r_jkse, color = "Block Jack Knife")) +
    geom_density(data = ldsc_estimates %>% filter(effect_estimated == "Population"), 
               aes(r_se, color = "LDSC Block Jack Knife")) +
    geom_vline(aes(xintercept = true_sd_r_pop_snipar, linetype = "True SD SNIPar")) +
    geom_vline(aes(xintercept = true_sd_r_pop_ldsc, linetype = "True SD LDSC")) +
    scale_color_brewer(palette = "Dark2") +
    theme_minimal() +
    labs(color = "",
         x = "SE",
         title = "Genetic Correlation - Population")

ggsave("ldsc_reg/multigen_sim/05012020/se_pop.png",
       height = 6, width = 9)   

# Comparing h2
ggplot() +
    geom_density(data = pop,
                   aes(v1_invhse, color = "Inverse Hessian")) +
    geom_density(data = pop,
                 aes(v1_jkse, color = "Block Jack Knife")) +
    geom_density(data = ldsc_estimates %>% filter(effect_estimated == "Population"), 
               aes(h2_dir_se, color = "LDSC Block Jack Knife")) +
    geom_vline(aes(xintercept = true_sd_h21_pop_snipar, linetype = "SNIPar - True SD")) +
    geom_vline(aes(xintercept = true_sd_h21_pop_ldsc, linetype = "LDSC - True SD")) +
    scale_color_brewer(palette = "Dark2") +
    theme_minimal() +
    labs(color = "",
         linetype = "",
         x = "SE",
         title = "Heritability of Direct Effect - Population")



ggplot() +
    geom_density(data = avg_parental,
                 aes(v1_invhse, color = "Inverse Hessian")) +
    geom_density(data = avg_parental,
                 aes(v1_jkse, color = "Block Jack Knife")) +
    geom_density(data = ldsc_estimates %>% filter(effect_estimated == "Average Parental"), 
                 aes(h2_dir_se, color = "LDSC Block Jack Knife")) +
    geom_vline(aes(xintercept = true_sd_h21_pop_snipar, linetype = "SNIPAR - True SD")) +
    geom_vline(aes(xintercept = true_sd_h21_pop_ldsc, linetype = "LDSC - True SD")) +
    scale_color_brewer(palette = "Dark2") +
    theme_minimal() +
    labs(color = "",
         linetype = "",
         x = "SE",
         title = "Heritability of Direct Effect - Average Parental")


# Comparing estimates
ggplot() +
    geom_density(data = pop, 
               aes(r, color = "SNIPar")) +
    geom_density(data = ldsc_estimates %>% filter(effect_estimated == "Population"), 
               aes(r, color = "LDSC")) +
    geom_vline(xintercept = 1.0) +
    scale_color_brewer(palette = "Dark2") +
    theme_minimal() +
    labs(x = "Run Number", y = "Estimate",
         color = "")

ggsave("ldsc_reg/multigen_sim/05012020/estimates_pop.png",
       height = 6, width = 9)  

ggplot() +
    geom_density(data = avg_parental, 
               aes(r, color = "SNIPar")) +
    geom_density(data = ldsc_estimates %>% filter(effect_estimated == "Average Parental"), 
               aes(r, color = "LDSC")) +
    geom_vline(xintercept = 0.5) +
    scale_color_brewer(palette = "Dark2") +
    theme_minimal() +
    labs(color = "")

ggsave("ldsc_reg/multigen_sim/05012020/estimates_avgparental.png",
       height = 6, width = 9)   


ggplot() +
    geom_density(data = avg_parental, 
               aes(v1, color = "SNIPar")) +
    geom_density(data = ldsc_estimates %>% filter(effect_estimated == "Average Parental"), 
               aes(h2_dir, color = "LDSC")) +
    scale_color_brewer(palette = "Dark2") +
    theme_minimal() +
    labs(color = "")

ggplot() +
    geom_density(data = pop, 
               aes(v2, color = "SNIPar")) +
    geom_density(data = ldsc_estimates %>% filter(effect_estimated == "Population"), 
               aes(h2_second, color = "LDSC")) +
    scale_color_brewer(palette = "Dark2") +
    theme_minimal() +
    labs(color = "")

