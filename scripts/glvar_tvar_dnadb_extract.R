library(tidyverse)

source(here::here("scripts/connect_to_sql_server.R"))

gi_csv_cleaned <- read_csv(file = paste0(config::get(
  "data_folderpath"),
  "02_cleaned/",
  "gi_csv_cleaned.csv"),
  col_types = list(
    nhsno = col_character()
  ))

gi_nhsnos <- gi_csv_cleaned$nhsno

# Identify all lab numbers from patients ----------------------------------

labno_df <- sample_tbl |> 
  filter(nhsno %in% gi_nhsnos) |> 
  select(labno, nhsno) |> 
  collect()

labno_query <- unique(labno_df$labno)

# Extract all results for lab numbers -------------------------------------

dnadb_results <- results_tbl |> 
  filter(labno %in% labno_query) |> 
  select(labno, pcrid, test, genotype, genotype2,
         genocomm) |> 
  collect()


# Extract germline variant results ----------------------------------------

glvar_dnadb_results <- dnadb_results |> 
  filter(test %in% unique(grep(pattern = "hs2(\\s|)icp", 
                               x = dnadb_results$test, 
                               ignore.case = TRUE, 
                               value = TRUE))) |> 
  left_join(labno_df, by = "labno") |> 
  relocate(nhsno)

# Extract tumour variant results ------------------------------------------

tvar_dnadb_results <- dnadb_results |> 
  filter(test %in% unique(grep(pattern = "seq\\span", 
                               x = dnadb_results$test, 
                               ignore.case = TRUE, 
                               value = TRUE))) |> 
  left_join(labno_df, by = "labno") |> 
  relocate(nhsno)

# Export data -------------------------------------------------------------

write_csv(glvar_dnadb_results,
          paste0(config::get(
            "data_folderpath"),
            "01_initial/",
            "glvar_dnadb_extracted.csv"))

write_csv(tvar_dnadb_results,
          paste0(config::get(
            "data_folderpath"),
            "01_initial/",
            "tvar_dnadb_extracted.csv"))
