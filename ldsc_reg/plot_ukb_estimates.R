library(tidyverse)
library(jsonlite)
theme_set(theme_minimal())


# reading in UKB data
ldsc_path <- "C:/Users/Hariharan/Documents/git_repos/SNIPar/ldsc_reg/"
log_files <- Sys.glob(paste0(ldsc_path, "ukb/*.log"))


estimates_df <- tibble(
    
    traitcode = numeric(),
    v1 = numeric(),
    v2 = numeric(),
    r = numeric(),
    v1_invhse = numeric(),
    v2_invhse = numeric(),
    r_invhse = numeric(),
    v1_jkse = numeric(),
    v2_jkse = numeric(),
    r_jkse = numeric()
    
)
for (file in log_files) {
    
    logdata <- read_delim(file, delim = "\n")
    
    traitno <- str_sub(file, -6) %>% 
        str_replace(".log", "") %>%
        as.numeric()
    
    estimates <- logdata[11, ] %>% pull
    estimates <- substr(estimates, 12, nchar(estimates))
    estimates <- str_replace_all(estimates, "\'", "")
    estimates <- str_replace_all(estimates, "\\{", "")
    estimates <- str_replace_all(estimates, "\\}", "")
    estimates <- str_replace_all(estimates, "  ", " ")
    estimates <- str_split(estimates, ",")[[1]]
    v1 <- str_split(estimates[1], ":")[[1]][2] %>% as.numeric()
    v2 <- str_split(estimates[2], ":")[[1]][2] %>% as.numeric()
    r <- str_split(estimates[3], ":")[[1]][2] %>% as.numeric()
    
    invh_se <- logdata[12, ] %>% pull
    invh_se <- str_split(invh_se, ":")[[1]][2]
    invh_se <- str_sub(invh_se, 2)
    invh_se <- str_replace_all(invh_se, "\\[", "")
    invh_se <- str_replace_all(invh_se, "\\]", "")
    invh_se <- str_replace_all(invh_se, "  ", " ")
    invh_se <- str_split(invh_se, " ")[[1]]
    v1_invhse <- invh_se[1] %>% as.numeric()
    v2_invhse <- invh_se[2] %>% as.numeric()
    r_invhse <- invh_se[3] %>% as.numeric()
    
    jkse <- logdata[25, ] %>% pull
    jkse <- str_split(jkse, ":")[[1]][2]
    jkse <- str_sub(jkse, 2)
    jkse <- str_replace_all(jkse, "\\[", "")
    jkse <- str_replace_all(jkse, "\\]", "")
    jkse <- str_replace_all(jkse, "  ", " ")
    jkse <- str_split(jkse, " ")[[1]]
    v1_jkse <- jkse[1] %>% as.numeric()
    v2_jkse <- jkse[2] %>% as.numeric()
    r_jkse <- jkse[3] %>% as.numeric()
    
    
    estimates_toappend <- tibble(
        
        traitcode = traitno,
        v1 = v1,
        v2 = v2,
        r = r,
        v1_invhse = v1_invhse,
        v2_invhse = v2_invhse,
        r_invhse = r_invhse,
        v1_jkse = v1_jkse,
        v2_jkse = v2_jkse,
        r_jkse = r_jkse
        
    )
    
    estimates_df <- estimates_df %>%
        bind_rows(estimates_toappend)
    
    
}

# getting trait names
traitnames <- read_delim(paste0(ldsc_path, "ukb/traits.txt"), 
                         delim = " ",
                         col_names = c("traitcode", "traitname"))
estimates_df <- estimates_df %>%
    left_join(traitnames, by = "traitcode")

# cleaning trait names
estimates_df <- estimates_df %>%
    mutate(traitname = str_to_title(traitname),
           traitname = str_replace(traitname, "\\.", " "),
           traitname = case_when(traitname == "Hdl" ~ "HDL",
                                 traitname == "Aafb" ~ "AAFB",
                                 traitname == "Fev1" ~ "FEV1",
                                 traitname == "Ea" ~ "EA",
                                 traitname == "Bmi" ~ "BMI",
                             TRUE ~ traitname))


# making confidence intervals
estimates_df <- estimates_df %>%
    mutate(r_ci_lo_invh = r - 1.96 * r_invhse,
           r_ci_hi_invh = r + 1.96 * r_invhse)

estimates_df <- estimates_df %>%
    mutate(r_ci_lo_jk = r - 1.96 * r_jkse,
           r_ci_hi_jk = r + 1.96 * r_jkse)

# plotting
estimates_df %>%
    ggplot() +
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

ggsave(paste0(ldsc_path, "ukb/ukb_estimates.png"),
       height = 6,
       width = 9)    


estimates_df %>%
    ggplot() +
    geom_point(aes(r_invhse, r_jkse)) +
    lims(x = c(0, 0.2), y = c(0, 0.2)) +
    geom_abline(aes(intercept = 0, slope = 1, color = "45 Degree Line")) +
    annotate(geom="text", x=0.12, y=0.13, label="45 Degree Line", color="red", angle = 45) +
    labs(x = "Inverse Hessian SE", y = "Block Jack Knife SE") +
    guides(color = FALSE)

ggsave(paste0(ldsc_path, "ukb/ukb_ses.png"),
       height = 9,
       width = 9)    
