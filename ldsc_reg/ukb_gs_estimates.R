library(tidyverse)
theme_set(theme_minimal())


# reading in data
estimates <- read_csv("ukb_gs_estimates.csv")

# making confidence intervals
estimates <- estimates %>%
    mutate(r_ci_lo = r - 1.96 * se_r,
           r_ci_hi = r + 1.96 * se_r)

# cleaning trait names
estimates <- estimates %>%
    mutate(Trait = str_to_title(Trait),
           Trait = str_replace(Trait, "\\.", " "),
           Trait = case_when(Trait == "Hdl" ~ "HDL",
                             Trait == "Aafb" ~ "AAFB",
                             Trait == "Fev1" ~ "FEV1",
                             Trait == "Ea" ~ "EA",
                             Trait == "Bmi" ~ "BMI",
                             TRUE ~ Trait))

estimates <- estimates %>%
    mutate(Dataset = factor(Dataset,
                            levels = c("UKB", "Generation Scotland")))

# plotting
estimates %>%
    ggplot() +
    geom_point(aes(Trait, r, 
                   color = Trait, shape = Dataset),
               size = 3,
               position = position_dodge2(padding = 0, preserve = "single", width = 0.6)) +
    geom_errorbar(aes(x = Trait, ymin = r_ci_lo, ymax = r_ci_hi,
                      color = Trait),
                  width = 0.6, position = position_dodge2(width=0, 
                                                          preserve = "single",
                                                          padding = 0.5)) +
    labs(x = "Trait", y = "Genetic Correlation") +
    guides(color = FALSE) +
    theme(legend.position = c(0.3, 0.9))

ggsave("ukb_gs_estimates.png",
       height = 6,
       width = 13)
