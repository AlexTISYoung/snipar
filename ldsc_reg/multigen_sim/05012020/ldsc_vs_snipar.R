library(dplyr)
library(readr)
library(stringr)

paths <- "/disk/genetics/ukb/alextisyoung/haplotypes/simulated_pops_large/from_chr1_to_chr23_start0_endNone_run0_p0-0_ab_corr0-5_vb0-25_length2/phenotype_dir_par_corr_0.5/*"

read_estimates <- function(folders, estimate_filename){

    estimate_filename <- paste0("/", estimate_filename)

    dfout <- tibble(run = integer())
    for (folder in folders){

        file <- paste0(folder, estimate_filename)
        est <- read_delim(file, delim = "\n")


        h21_line <- est[33,] %>% pull
        h21_est <- str_sub(h21_line, 24, str_length(h21_line))
        h21 <- str_extract(h21_est, "\\d+.\\d+") %>% 
            as.numeric()
        h21_se <- str_extract(h21_est, "\\(\\d+.\\d+\\)") %>% 
            str_remove("\\(") %>% 
            str_remove("\\)") %>% 
            as.numeric()

        h22_line <- est[39,] %>% pull
        h22_est <- str_sub(h22_line, 24, str_length(h22_line))
        h22 <- str_extract(h22_est, "\\d+.\\d+") %>% 
            as.numeric()
        h22_se <- str_extract(h22_est, "\\(\\d+.\\d+\\)") %>% 
            str_remove("\\(") %>% 
            str_remove("\\)") %>% 
            as.numeric()

        cor_line <- est[50,] %>% pull
        cor_est <- str_sub(cor_line, 20, str_length(cor_line))
        cor <- str_extract(cor_est, "\\d+.\\d+") %>% 
            as.numeric()
        cor_se <- str_extract(cor_est, "\\(\\d+.\\d+\\)") %>% 
            str_remove("\\(") %>% 
            str_remove("\\)") %>% 
            as.numeric()


        df_toappend <- tibble(run = str_sub(folder, 162, str_length(folder)) %>% as.numeric(),
                                h2_dir = h21,
                                h2_dir_se = h21_se,
                                h2_second = h22,
                                h2_second_se = h22_se,
                                r = cor,
                                r_se = cor_se)

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


