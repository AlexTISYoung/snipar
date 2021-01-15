library(dplyr)
library(readr)
library(stringr)

paths <- "/disk/genetics/ukb/alextisyoung/haplotypes/simulated_pops_large/from_chr1_to_chr23_start0_endNone_run0_p0-0_ab_corr0-5_vb0-25_length2/phenotype_dir_par_corr_0.5/*"

read_estimates <- function(folders, estimate_filename){

    estimate_filename <- paste0("/", estimate_filename)

    dfout <- tibble(run = integer(), estimate = numeric(), se = numeric())
    for (folder in folders){

        file <- paste0(folder, estimate_filename)
        est <- read_delim(file, delim = "\n")

        cor_line <- est[50,] %>% pull
        cor_est <- str_sub(cor_line, 20, str_length(cor_line))
        estimate <- str_extract(cor_est, "\\d+.\\d+") %>% 
            as.numeric()
        se <- str_extract(cor_est, "\\(\\d+.\\d+\\)") %>% 
            str_remove("\\(") %>% 
            str_remove("\\)") %>% 
            as.numeric()

        df_toappend <- tibble(run = str_sub(folder, 162, str_length(folder)) %>% as.numeric(),
                                estimate = estimate,
                                se = se)

        dfout <- bind_rows(dfout, df_toappend)
        
    }

    return(dfout)

}

replicates <- Sys.glob(paths)

avgparental <- read_estimates(replicates, "rg_delta_eta_par_constrain.log") %>%
    mutate(effect_estimated = "Average Parental")
population <- read_estimates(replicates, "rg_delta_beta_constrain.log") %>%
    mutate(effect_estimated = "Population")

ldsc_estimates <- bind_rows(avgparental, population)

ldsc_estimates %>% 
    write_csv("ldsc_reg/multigen_sim/05012020/ldsc_estimates.csv")


