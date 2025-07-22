library(tidyverse)

# Load initial data -------------------------------------------------------

glvar_dnadb_extracted <- read_csv(
  paste0(config::get("data_folderpath"),
         "01_initial/",
         "glvar_dnadb_extracted.csv"))

tvar_dnadb_extracted <- read_csv(
  paste0(config::get("data_folderpath"),
         "01_initial/",
         "tvar_dnadb_extracted.csv"))

# Add headline to germline results ----------------------------------------

glvar_dnadb_cleaned <- glvar_dnadb_extracted |> 
  mutate(glvar_headline_result = case_when(
    genotype %in% c("No pathogenic variant identified; No relevant CNV identified",
                    "No pathogenic variant identified; CNV Analysis Failed") ~"No reportable variant(s) detected",
    genotype %in% c("Fail") ~"Analysis failed (see quality score)",
    TRUE ~"Reportable variant(s) detected"))

# Add headline to tumour variant results ----------------------------------

no_path_var_strings <- unique(grep(pattern = "no\\spathogenic",
            x = tvar_dnadb_extracted$genotype,
            ignore.case = TRUE,
            value = TRUE))

fail_strings <- c("Analysis failed", "Fail")

tvar_dnadb_cleaned <- tvar_dnadb_extracted |> 
  mutate(tvar_headline_result = case_when(
    genotype %in% no_path_var_strings ~"No reportable variant(s) detected",
    genotype %in% fail_strings ~"Analysis failed (see quality score)",
    TRUE ~"Reportable variant(s) detected"
  ))

# Export results ----------------------------------------------------------

write_csv(glvar_dnadb_cleaned,
          paste0(config::get("data_folderpath"),
                 "02_cleaned/",
                 "glvar_dnadb_cleaned.csv"))

write_csv(tvar_dnadb_cleaned,
          paste0(config::get("data_folderpath"),
                 "02_cleaned/",
                 "tvar_dnadb_cleaned.csv"))
