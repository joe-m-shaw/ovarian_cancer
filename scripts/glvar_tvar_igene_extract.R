library(tidyverse, quietly  = TRUE)
library(janitor, quietly  = TRUE)

# Connect to SQL server ---------------------------------------------------

source(here::here("scripts/connect_to_sql_server.R"))

test_names <- c("PANEL: M2.5 - SeqOne HRD Status",
                "PANEL: R207.1 - Inherited ovarian cancer (without breast cancer) v4.0 (ICP)",
                "PANEL: M2_tBRCA_PS",
                "PANEL: M2_tBRCA_PS v2.a",
                "PANEL: M2_tBRCA_PS v2.b")

stopifnot(nrow(eval_hrd) != 0)
stopifnot(length(setdiff(unique(eval_hrd$test_name), test_names)) == 0)
stopifnot(ncol(eval_hrd) == 8)

# Function ----------------------------------------------------------------

message("Extracting glvar and tvar results from iGene")

widen_by_test <- function(df = eval_hrd, 
                          test_vector,
                          prefix_string) {
  
  output <- eval_hrd |> 
    filter(test_name %in% test_vector) |> 
    pivot_wider(id_cols = c(referral_number, 
                            nhsn_umber,
                            test_order_date, 
                            test_identifier),
                names_from = box,
                values_from = data_value,
                names_prefix = prefix_string) |> 
    clean_names() |> 
    rename(nhsno = nhsn_umber)
  
  if(nrow(output) == 0){
    stop("Output dataframe is empty")
  }
  
  return(output)
  
}

# Extract tumour and germline variant results -----------------------------

glvar_igene <- widen_by_test(df = eval_hrd,
                              test_vector = c("PANEL: R207.1 - Inherited ovarian cancer (without breast cancer) v4.0 (ICP)"), 
                              prefix_string = "glvar")

stopifnot(ncol(glvar_igene) == 23)

message(paste0(nrow(glvar_igene), " germline tests identified for these panels: ",
               "PANEL: R207.1 - Inherited ovarian cancer (without breast cancer) v4.0 (ICP)"))

tvar_igene <- widen_by_test(df = eval_hrd,
                             test_vector = c("PANEL: M2_tBRCA_PS",
                                             "PANEL: M2_tBRCA_PS v2.a",
                                             "PANEL: M2_tBRCA_PS v2.b"), 
                             prefix_string = "tvar")

message(paste0(nrow(tvar_igene), " tumour tests identified for these panels: ", 
               paste(c("PANEL: M2_tBRCA_PS",
                       "PANEL: M2_tBRCA_PS v2.a",
                       "PANEL: M2_tBRCA_PS v2.b"),
                     collapse = ", ")))

gi_igene <- widen_by_test(df = eval_hrd,
                          test_vector = c("PANEL: M2.5 - SeqOne HRD Status"), 
                          prefix_string = "gi")

message(paste0(nrow(gi_igene), " SeqOne GI tests identified for these panels: ",
               "PANEL: M2.5 - SeqOne HRD Status"))

# Perform checks ----------------------------------------------------------

message("Checking extracted glvar and tvar iGene data")

stopifnot(nrow(glvar_igene) != 0)
stopifnot(nrow(tvar_igene) != 0)
stopifnot(nrow(gi_igene) != 0)

tests_in_eval_hrd <- unique(eval_hrd$test_identifier)

tests_in_wide_tables <- c(unique(glvar_igene$test_identifier),
                          unique(tvar_igene$test_identifier),
                          unique(gi_igene$test_identifier))

# Check no tests have been lost from the dataset
stopifnot(length(setdiff(tests_in_eval_hrd, tests_in_wide_tables)) == 0)

message(paste0("There are ", 
               length(tests_in_eval_hrd), 
               " tests in eval_hrd and ",
               length(tests_in_wide_tables), 
               " tests in the reformatted data."))

# Remove entries without headline results ---------------------------------

glvar_igene_na <- glvar_igene |> 
  filter(is.na(glvar_headline_result))

glvar_igene_no_na <- glvar_igene |> 
  filter(!is.na(glvar_headline_result))

message(paste0(nrow(glvar_igene_na),
               " germline entries with no headline result were removed."))

tvar_igene_na <- tvar_igene |> 
  filter(is.na(tvar_headline_result))

tvar_igene_no_na <- tvar_igene |> 
  filter(!is.na(tvar_headline_result))

message(paste0(nrow(tvar_igene_na),
               " tumour entries with no headline result were removed."))

# Check dates of entries --------------------------------------------------

message(paste0("iGene germline variant data ranges from ",
               format.Date(x = min(glvar_igene_no_na$test_order_date), 
                           format = "%d %B %Y"),
               " to ",
               format.Date(x = max(glvar_igene_no_na$test_order_date), 
                           format = "%d %B %Y")))

message(paste0("iGene tumour variant data ranges from ",
               format.Date(x = min(tvar_igene_no_na$test_order_date), 
                           format = "%d %B %Y"),
               " to ",
               format.Date(x = max(tvar_igene_no_na$test_order_date), 
                           format = "%d %B %Y")))

# Export results ----------------------------------------------------------

message("Exporting iGene results")

write_csv(glvar_igene_no_na,
          paste0(config::get("data_folderpath"),
                 "01_initial/",
                 "glvar_igene_extracted.csv"))

write_csv(tvar_igene_no_na,
          paste0(config::get("data_folderpath"),
                 "01_initial/",
                 "tvar_igene_extracted.csv"))
