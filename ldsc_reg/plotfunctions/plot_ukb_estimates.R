library(tidyverse)
theme_set(theme_minimal())

ldsc_path <- "C:/Users/Hariharan/Documents/git_repos/SNIPar/ldsc_reg/"
source(paste0(ldsc_path, "plotfunctions/helperfuncs.R"))

# reading in UKB data
log_files <- Sys.glob(paste0(ldsc_path, "ukb/[0-9]*.log"))
log_files_rbound <- log_files[!str_detect(log_files, "norbound")]
log_files_norbound <- log_files[str_detect(log_files, "norbound")]


# getting trait names
traitnames <- read_delim(paste0(ldsc_path, "ukb/traits.txt"), 
                         delim = " ",
                         col_names = c("traitcode", "traitname")) %>%
    mutate(traitname = str_to_title(traitname),
           traitname = str_replace(traitname, "\\.", " "),
           traitname = case_when(traitname == "Hdl" ~ "HDL",
                                 traitname == "Aafb" ~ "AAFB",
                                 traitname == "Fev1" ~ "FEV1",
                                 traitname == "Ea" ~ "EA",
                                 traitname == "Bmi" ~ "BMI",
                                 TRUE ~ traitname))


# reading in data
estimates_df_rbound <- readlogfiles(log_files_rbound)
estimates_df_norbound <- readlogfiles(log_files_norbound)

estimates_df_rbound <- estimates_df_rbound %>%
    left_join(traitnames, by = "traitcode")

estimates_df_norbound <- estimates_df_norbound %>%
    left_join(traitnames, by = "traitcode")


# making confidence intervals
estimates_df_rbound <- estimates_df_rbound %>%
    mutate(r_ci_lo_invh = r - 1.96 * r_invhse,
           r_ci_hi_invh = r + 1.96 * r_invhse)

estimates_df_rbound <- estimates_df_rbound %>%
    mutate(r_ci_lo_jk = r - 1.96 * r_jkse,
           r_ci_hi_jk = r + 1.96 * r_jkse)

estimates_df_norbound <- estimates_df_norbound %>%
    mutate(r_ci_lo_invh = r - 1.96 * r_invhse,
           r_ci_hi_invh = r + 1.96 * r_invhse)

estimates_df_norbound <- estimates_df_norbound %>%
    mutate(r_ci_lo_jk = r - 1.96 * r_jkse,
           r_ci_hi_jk = r + 1.96 * r_jkse)

# plotting
estimates_df_rbound %>%
    ggplot() +
    geom_hline(yintercept = 1.0, color = 'darkgray') +
    geom_point(aes(reorder(traitname, r), 
                   r, color = traitname),
               size = 3) +
    geom_errorbar(aes(x = traitname, 
                      ymin = r_ci_lo_invh, 
                      ymax = r_ci_hi_invh,
                      color = traitname,
                      linetype = "Inverse Hessian")) +
    geom_errorbar(aes(x = traitname, 
                      ymin = r_ci_lo_jk, 
                      ymax = r_ci_hi_jk,
                      color = traitname,
                      linetype = "Block Jack Knife")) +
    labs(x = "Trait", y = "Genetic Correlation", linetype = "") +
    guides(color = FALSE) +
    scale_linetype_manual(values = c("solid",
                                     "dashed")) +
    theme(legend.position = c(0.8, 0.2))

ggsave(paste0(ldsc_path, "ukb/ukb_rbound.png"),
       height = 6,
       width = 9)    


estimates_df_rbound %>%
    ggplot() +
    geom_point(aes(r_invhse, r_jkse)) +
    lims(x = c(0, 0.2), y = c(0, 0.2)) +
    geom_abline(aes(intercept = 0, slope = 1, color = "45 Degree Line")) +
    annotate(geom="text", x=0.12, y=0.13, label="45 Degree Line", color="red", angle = 45) +
    labs(x = "Inverse Hessian SE", y = "Block Jack Knife SE") +
    guides(color = FALSE)

ggsave(paste0(ldsc_path, "ukb/ukb_ses_rbound.png"),
       height = 9,
       width = 9)    



estimates_df_norbound %>%
    filter(traitcode != 15) %>%
    ggplot() +
    geom_hline(yintercept = 1.0, color = 'darkgray') +
    geom_point(aes(reorder(traitname, r), 
                   r, color = traitname),
               size = 3) +
    geom_errorbar(aes(x = traitname, 
                      ymin = r_ci_lo_invh, 
                      ymax = r_ci_hi_invh,
                      color = traitname,
                      linetype = "Inverse Hessian")) +
    geom_errorbar(aes(x = traitname, 
                      ymin = r_ci_lo_jk, 
                      ymax = r_ci_hi_jk,
                      color = traitname,
                      linetype = "Block Jack Knife")) +
    labs(x = "Trait", y = "Genetic Correlation", linetype = "") +
    guides(color = FALSE) +
    scale_linetype_manual(values = c("solid",
                                     "dashed")) +
    theme(legend.position = c(0.8, 0.2))

ggsave(paste0(ldsc_path, "ukb/ukb_norbound.png"),
       height = 6,
       width = 9)    


estimates_df_norbound %>%
    filter(traitcode != 15) %>%
    ggplot() +
    geom_point(aes(r_invhse, r_jkse)) +
    lims(x = c(0, 1), y = c(0, 1)) +
    geom_abline(aes(intercept = 0, slope = 1, color = "45 Degree Line")) +
    annotate(geom="text", x=0.12, y=0.13, label="45 Degree Line", color="red", angle = 45) +
    labs(x = "Inverse Hessian SE", y = "Block Jack Knife SE") +
    guides(color = FALSE)

ggsave(paste0(ldsc_path, "ukb/ukb_ses_norbound.png"),
       height = 9,
       width = 9)    
