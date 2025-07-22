
library(tidyverse)
library(janitor)


# Connect to SQL server ---------------------------------------------------

source(here::here("scripts/connect_to_sql_server.R"))

# Function ----------------------------------------------------------------

widen_by_test <- function(df = eval_hrd, 
                          test_string,
                          prefix_string) {
  
  output <- eval_hrd |> 
    filter(test_name == test_string) |> 
    select(-nhsn_umber) |> 
    pivot_wider(id_cols = c(referral_number, 
                            test_order_date, 
                            test_identifier),
                names_from = field,
                values_from = data_value,
                names_prefix = prefix_string) |> 
    clean_names()
  
  if(nrow(output) == 0){
    stop("Output dataframe is empty")
  }
  
  return(output)
  
}

# Extract tumour and germline variant results -----------------------------

glvar_igene <- widen_by_test(df = eval_hrd,
                              test_string = "PANEL: R207.1 - Inherited ovarian cancer (without breast cancer) v4.0 (ICP)", 
                              prefix_string = "glvar")

tvar_igene <- widen_by_test(df = eval_hrd,
                             test_string = "PANEL: M2_tBRCA_PS", 
                             prefix_string = "tvar") 

# Export results ----------------------------------------------------------

write_csv(glvar_igene,
          paste0(config::get("data_folderpath"),
                 "01_initial/",
                 "glvar_igene_extracted.csv"))

write_csv(tvar_igene,
          paste0(config::get("data_folderpath"),
                 "01_initial/",
                 "tvar_igene_extracted.csv"))
