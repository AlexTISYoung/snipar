readlogfiles <- function(log_files){
    
    # Reads in logfiles and converts them into
    # a tibble
    

    
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
        
        print(file)
        logdata <- read_delim(file, delim = "\n")
        logdata <- logdata[[1]]
        
        traitno <- str_split(file, "/")[[1]]
        traitno <- traitno[length(traitno)]
        traitno <- str_extract(traitno, "\\d+")
        traitno <- as.numeric(traitno)
        
        estimates <- logdata[str_which(logdata, "Estimates")]
        if (length(estimates) != 0){
            estimates <- substr(estimates, 12, nchar(estimates))
            estimates <- str_replace_all(estimates, "\'", "")
            estimates <- str_replace_all(estimates, "\\{", "")
            estimates <- str_replace_all(estimates, "\\}", "")
            estimates <- str_replace_all(estimates, "  ", " ")
            estimates <- str_split(estimates, ",")[[1]]
            v1 <- str_split(estimates[1], ":")[[1]][2] %>% as.numeric()
            v2 <- str_split(estimates[2], ":")[[1]][2] %>% as.numeric()
            r <- str_split(estimates[3], ":")[[1]][2] %>% as.numeric()
        
            invh_se <- logdata[str_which(logdata, "Standard Errors:")]
            invh_se <- str_extract_all(invh_se, "\\d+.\\d+|nan")[[1]]
            v1_invhse <- invh_se[1] %>% as.numeric()
            v2_invhse <- invh_se[2] %>% as.numeric()
            r_invhse <- invh_se[3] %>% as.numeric()
            
            jkse <- logdata[str_which(logdata, "Block Jack Knife Standard Errors")] 
            jkse <- str_extract_all(jkse, "\\d+.\\d+")[[1]]
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
        
        
    }
    
    return(estimates_df)
}