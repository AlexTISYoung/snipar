library(tidyverse)
library(readxl)
theme_set(theme_minimal())
ldsc_path <- "C:/Users/Hariharan/Documents/git_repos/SNIPar/ldsc_reg/"
source(paste0(ldsc_path, "plotfunctions/helperfuncs.R"))


# reading in UKB data
log_files <- Sys.glob(paste0(ldsc_path, "ukb/[0-9]*.log"))
log_files_rbound <- log_files[!str_detect(log_files, "norbound")]
log_files_avgparental <- log_files[str_detect(log_files, "norbound_avgparental")]
log_files_norbound <- log_files[str_detect(log_files, "norbound") & !str_detect(log_files, "avgparental")]



# getting trait names
traitnames <- read_delim(paste0(ldsc_path, "ukb/traits.txt"), 
                         delim = " ",
                         col_names = c("traitcode", "traitname"))  %>%
    mutate(traitname = str_to_title(traitname),
           traitname = str_replace(traitname, "\\.", " "),
           traitname = case_when(traitname == "Hdl" ~ "HDL",
                                 traitname == "Aafb" ~ "AAFB",
                                 traitname == "Fev1" ~ "FEV1",
                                 traitname == "Ea" ~ "EA",
                                 traitname == "Bmi" ~ "BMI",
                                 TRUE ~ traitname))


# reading in data
estimates_df_rbound <- readlogfiles(log_files_rbound) %>%
    left_join(traitnames, by = "traitcode")
estimates_df_norbound <- readlogfiles(log_files_norbound) %>%
    left_join(traitnames, by = "traitcode")
estimates_df_avgparental <- readlogfiles(log_files_avgparental) %>%
    left_join(traitnames, by = "traitcode")


# making confidence intervals
estimates_df_rbound <- estimates_df_rbound %>%
    mutate(r_ci_lo_invh = r - 1.96 * r_invhse,
           r_ci_hi_invh = r + 1.96 * r_invhse,
           r_ci_lo_jk = r - 1.96 * r_jkse,
           r_ci_hi_jk = r + 1.96 * r_jkse)


estimates_df_norbound <- estimates_df_norbound %>%
    mutate(r_ci_lo_invh = r - 1.96 * r_invhse,
           r_ci_hi_invh = r + 1.96 * r_invhse,
           r_ci_lo_jk = r - 1.96 * r_jkse,
           r_ci_hi_jk = r + 1.96 * r_jkse)

estimates_df_avgparental <- estimates_df_avgparental %>%
    mutate(r_ci_lo_invh = r - 1.96 * r_invhse,
           r_ci_hi_invh = r + 1.96 * r_invhse,
           r_ci_lo_jk = r - 1.96 * r_jkse,
           r_ci_hi_jk = r + 1.96 * r_jkse)


# reading in ldsc estimates for horserace
ldsc <- read_excel(paste0(ldsc_path, "ukb/ldsc_estimates.xlsx")) %>%
    rename(traitcode = trait)

ldsc <- ldsc %>%
    mutate(r_ci_lo = r - 1.96 * rse,
           r_ci_hi = r + 1.96 * rse)

ldsc <- ldsc %>%
    left_join(traitnames, by = "traitcode")

ldsc_restricted <- ldsc %>%
    filter(restricted == "yes")

# plots
ggplot() +
    geom_hline(yintercept = 1.0, color = 'darkgray') +
    geom_point(data = estimates_df_rbound,
               aes(reorder(traitname, r), 
                   r, color = traitname,
                   shape = "SNIPAR"),
               size = 3, position = position_nudge(x = -0.1)) +
    geom_errorbar(data = estimates_df_rbound,
                  aes(x = traitname, 
                      ymin = r_ci_lo_invh, 
                      ymax = r_ci_hi_invh,
                      color = traitname,
                      linetype = "SNIPAR"),
                  width = 0.2,
                  position = position_nudge(x = -0.1)) +
    geom_point(data = ldsc_restricted %>% filter(traitcode %in% estimates_df_rbound$traitcode &
                                                 effect_estimated == "population"),
               aes(reorder(traitname, estimates_df_rbound$r), 
                   r, color = traitname,
                   shape = "LDSC"),
               size = 3, position = position_nudge(x = 0.1)) +
    geom_errorbar(data = ldsc_restricted %>% filter(traitcode %in% estimates_df_rbound$traitcode &
                                                        effect_estimated == "population"),
                  aes(x = traitname, 
                      ymin = r_ci_lo, 
                      ymax = r_ci_hi,
                      color = traitname,
                      linetype = "LDSC"),
                  width = 0.2,
                  position = position_nudge(x = 0.1)) +
    guides(color = FALSE) +
    labs(linetype = "", shape = "", 
         y = "",
         x = "",
         title = "Horserace - Genetic Correlation Estimates",
         caption = "Dataset: UKB") +
    scale_linetype_manual(values = c("longdash", "solid")) +
    scale_shape_manual(values = c(17, 16)) + 
    theme(legend.position = c(0.8, 0.2),
          plot.caption = element_text(face = "italic", color = "gray", hjust = 0))

ggsave(paste0(ldsc_path, "ukb/ukb_horserace_rbound.png"),
       height = 6,
       width = 9)  



ggplot() +
    geom_hline(yintercept = 1.0, color = 'darkgray') +
    geom_point(data = estimates_df_norbound,
               aes(reorder(traitname, r), 
                   r, color = traitname,
                   shape = "SNIPAR"),
               size = 3, position = position_nudge(x = -0.1)) +
    geom_errorbar(data = estimates_df_norbound,
                  aes(x = traitname, 
                      ymin = r_ci_lo_invh, 
                      ymax = r_ci_hi_invh,
                      color = traitname,
                      linetype = "SNIPAR"),
                  width = 0.2,
                  position = position_nudge(x = -0.1)) +
    geom_point(data = ldsc_restricted %>% filter(traitcode %in% estimates_df_norbound$traitcode &
                                                     effect_estimated == "population"),
               aes(reorder(traitname, r), 
                   r, color = traitname,
                   shape = "LDSC"),
               size = 3, position = position_nudge(x = 0.1)) +
    geom_errorbar(data = ldsc_restricted %>% filter(traitcode %in% estimates_df_norbound$traitcode &
                                                        effect_estimated == "population"),
                  aes(x = traitname, 
                      ymin = r_ci_lo, 
                      ymax = r_ci_hi,
                      color = traitname,
                      linetype = "LDSC"),
                  width = 0.2,
                  position = position_nudge(x = 0.1)) +
    guides(color = FALSE) +
    labs(linetype = "", shape = "", 
         y = "",
         x = "",
         title = "Horserace - Genetic Correlation Estimates",
         caption = "Dataset: UKB") +
    scale_linetype_manual(values = c("longdash", "solid")) +
    scale_shape_manual(values = c(17, 16)) + 
    theme(legend.position = c(0.8, 0.2),
          plot.caption = element_text(face = "italic", color = "gray", hjust = 0))

ggsave(paste0(ldsc_path, "ukb/ukb_horserace_norbound.png"),
       height = 6,
       width = 9)  


ggplot() +
    geom_hline(yintercept = 1.0, color = 'darkgray') +
    geom_hline(yintercept = -1.0, color = 'darkgray') +
    geom_point(data = estimates_df_avgparental,
               aes(reorder(traitname, r), 
                   r, color = traitname,
                   shape = "SNIPAR"),
               size = 3, position = position_nudge(x = -0.1)) +
    geom_errorbar(data = estimates_df_avgparental,
                  aes(x = traitname, 
                      ymin = r_ci_lo_invh, 
                      ymax = r_ci_hi_invh,
                      color = traitname,
                      linetype = "SNIPAR"),
                  width = 0.2,
                  position = position_nudge(x = -0.1)) +
    geom_point(data = ldsc_restricted  %>% filter(traitcode %in% estimates_df_avgparental$traitcode &
                                                      effect_estimated == "averageparental"),
               aes(reorder(traitname, r), 
                   r, color = traitname,
                   shape = "LDSC"),
               size = 3, position = position_nudge(x = 0.1)) +
    geom_errorbar(data = ldsc_restricted  %>% filter(traitcode %in% estimates_df_avgparental$traitcode &
                                                         effect_estimated == "averageparental"),
                  aes(x = traitname, 
                      ymin = r_ci_lo, 
                      ymax = r_ci_hi,
                      color = traitname,
                      linetype = "LDSC"),
                  width = 0.2,
                  position = position_nudge(x = 0.1)) +
    guides(color = FALSE) +
    labs(linetype = "", shape = "", 
         y = "",
         x = "",
         title = "Horserace - Genetic Correlation Estimates",
         caption = "Dataset: UKB") +
    scale_linetype_manual(values = c("longdash", "solid")) +
    scale_shape_manual(values = c(17, 16)) + 
    theme(legend.position = c(0.2, 0.8),
          plot.caption = element_text(face = "italic", color = "gray", hjust = 0))

ggsave(paste0(ldsc_path, "ukb/ukb_horserace_norbound_avgparental.png"),
       height = 6,
       width = 9)  

# with JKSE CI

ggplot() +
    geom_hline(yintercept = 1.0, color = 'darkgray') +
    geom_point(data = estimates_df_rbound,
               aes(reorder(traitname, r), 
                   r, color = traitname,
                   shape = "SNIPAR"),
               size = 3, position = position_nudge(x = -0.1)) +
    geom_errorbar(data = estimates_df_rbound,
                  aes(x = traitname, 
                      ymin = r_ci_lo_jk, 
                      ymax = r_ci_hi_jk,
                      color = traitname,
                      linetype = "SNIPAR"),
                  width = 0.2,
                  position = position_nudge(x = -0.1)) +
    geom_point(data = ldsc_restricted %>% filter(traitcode %in% estimates_df_rbound$traitcode &
                                                     effect_estimated == "population"),
               aes(reorder(traitname, estimates_df_rbound$r), 
                   r, color = traitname,
                   shape = "LDSC"),
               size = 3, position = position_nudge(x = 0.1)) +
    geom_errorbar(data = ldsc_restricted %>% filter(traitcode %in% estimates_df_rbound$traitcode &
                                                        effect_estimated == "population"),
                  aes(x = traitname, 
                      ymin = r_ci_lo, 
                      ymax = r_ci_hi,
                      color = traitname,
                      linetype = "LDSC"),
                  width = 0.2,
                  position = position_nudge(x = 0.1)) +
    guides(color = FALSE) +
    labs(linetype = "", shape = "", 
         y = "",
         x = "",
         title = "Horserace - Genetic Correlation Estimates",
         caption = "Dataset: UKB") +
    scale_linetype_manual(values = c("longdash", "solid")) +
    scale_shape_manual(values = c(17, 16)) + 
    theme(legend.position = c(0.8, 0.2),
          plot.caption = element_text(face = "italic", color = "gray", hjust = 0))


ggsave(paste0(ldsc_path, "ukb/ukb_horserace_rbound_jkse.png"),
       height = 6,
       width = 9)  

ggplot() +
    geom_hline(yintercept = 1.0, color = 'darkgray') +
    geom_point(data = estimates_df_norbound,
               aes(reorder(traitname, r), 
                   r, color = traitname,
                   shape = "SNIPAR"),
               size = 3, position = position_nudge(x = -0.1)) +
    geom_errorbar(data = estimates_df_norbound,
                  aes(x = traitname, 
                      ymin = r_ci_lo_jk, 
                      ymax = r_ci_hi_jk,
                      color = traitname,
                      linetype = "SNIPAR"),
                  width = 0.2,
                  position = position_nudge(x = -0.1)) +
    geom_point(data = ldsc_restricted %>% filter(traitcode %in% estimates_df_rbound$traitcode &
                                                     effect_estimated == "population"),
               aes(reorder(traitname, estimates_df_rbound$r), 
                   r, color = traitname,
                   shape = "LDSC"),
               size = 3, position = position_nudge(x = 0.1)) +
    geom_errorbar(data = ldsc_restricted %>% filter(traitcode %in% estimates_df_rbound$traitcode &
                                                        effect_estimated == "population"),
                  aes(x = traitname, 
                      ymin = r_ci_lo, 
                      ymax = r_ci_hi,
                      color = traitname,
                      linetype = "LDSC"),
                  width = 0.2,
                  position = position_nudge(x = 0.1)) +
    guides(color = FALSE) +
    labs(linetype = "", shape = "", 
         y = "",
         x = "",
         title = "Horserace - Genetic Correlation Estimates",
         caption = "Dataset: UKB") +
    scale_linetype_manual(values = c("longdash", "solid")) +
    scale_shape_manual(values = c(17, 16)) + 
    theme(legend.position = c(0.8, 0.2),
          plot.caption = element_text(face = "italic", color = "gray", hjust = 0))

ggsave(paste0(ldsc_path, "ukb/ukb_horserace_norbound_jkse.png"),
       height = 6,
       width = 9)  

ggplot() +
    geom_hline(yintercept = 1.0, color = 'darkgray') +
    geom_point(data = estimates_df_norbound %>% filter(!(traitcode %in% c(12, 4))),
               aes(reorder(traitname, r), 
                   r, color = traitname,
                   shape = "SNIPAR"),
               size = 3, position = position_nudge(x = -0.1)) +
    geom_errorbar(data = estimates_df_norbound %>% filter(!(traitcode %in% c(12, 4))),
                  aes(x = traitname, 
                      ymin = r_ci_lo_jk, 
                      ymax = r_ci_hi_jk,
                      color = traitname,
                      linetype = "SNIPAR"),
                  width = 0.2,
                  position = position_nudge(x = -0.1)) +
    geom_point(data = ldsc_restricted %>% filter(traitcode %in% estimates_df_rbound$traitcode &
                                                     effect_estimated == "population"),
               aes(reorder(traitname, estimates_df_rbound$r), 
                   r, color = traitname,
                   shape = "LDSC"),
               size = 3, position = position_nudge(x = 0.1)) +
    geom_errorbar(data = ldsc_restricted %>% filter(traitcode %in% estimates_df_rbound$traitcode &
                                                        effect_estimated == "population"),
                  aes(x = traitname, 
                      ymin = r_ci_lo, 
                      ymax = r_ci_hi,
                      color = traitname,
                      linetype = "LDSC"),
                  width = 0.2,
                  position = position_nudge(x = 0.1)) +
    guides(color = FALSE) +
    labs(linetype = "", shape = "", 
         y = "",
         x = "",
         title = "Horserace - Genetic Correlation Estimates",
         caption = "Dataset: UKB") +
    scale_linetype_manual(values = c("longdash", "solid")) +
    scale_shape_manual(values = c(17, 16)) + 
    theme(legend.position = c(0.8, 0.2),
          plot.caption = element_text(face = "italic", color = "gray", hjust = 0))

ggsave(paste0(ldsc_path, "ukb/ukb_horserace_norbound_jkse_constrained.png"),
       height = 6,
       width = 9)  


# SE Comparison

estimates_df_norbound %>% left_join(ldsc_restricted %>% filter(traitcode %in% estimates_df_norbound$traitcode &
                                                                   effect_estimated == "population"), 
                                    by = "traitcode", suffix = c(".snipar", ".ldsc")) %>%
    ggplot() +
    geom_hline(yintercept = 1.0, color = 'darkgray') +
    geom_point(aes(reorder(traitname.snipar, r.snipar), rse^2/r_invhse^2)) +
    labs(x = "", y = "(LDSC SE)^2/(SNIPAR SE)^2",
         title = "Comparing Standard Errors",
         subtitle = "R not Bounded")

ggsave(paste0(ldsc_path, "ukb/ukb_comparingse_norbound.png"),
       height = 6,
       width = 9)  

estimates_df_avgparental %>% left_join(ldsc_restricted %>% filter(traitcode %in% estimates_df_rbound$traitcode &
                                                                   effect_estimated == "averageparental"), 
                                    by = "traitcode", suffix = c(".snipar", ".ldsc")) %>%
    ggplot() +
    geom_hline(yintercept = 1.0, color = 'darkgray') +
    geom_point(aes(reorder(traitname.snipar, r.snipar), rse^2/r_invhse^2)) +
    labs(x = "", y = "(LDSC SE)^2/(SNIPAR SE)^2",
         title = "Comparing Standard Errors",
         subtitle = "R Bounded")

ggsave(paste0(ldsc_path, "ukb/ukb_comparingse_rbound.png"),
       height = 6,
       width = 9)  

