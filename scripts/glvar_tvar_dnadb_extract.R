library(tidyverse)

source(here::here("scripts/connect_to_sql_server.R"))

gi_csv_cleaned <- read_csv(file = paste0(config::get(
  "data_folderpath"),
  "02_cleaned/",
  "gi_csv_cleaned_orpp.csv"),
  col_types = list(
    nhsno = col_character()
  ))

stopifnot(nrow(gi_csv_cleaned) != 0)
stopifnot(anyNA(gi_csv_cleaned$nhsno) == FALSE)

gi_nhsnos <- gi_csv_cleaned$nhsno

# Identify all lab numbers from patients ----------------------------------

message("Finding all DNA Database lab numbers for patients with GI results")

labno_df <- sample_tbl |> 
  filter(nhsno %in% gi_nhsnos) |> 
  select(labno, nhsno) |> 
  collect()

stopifnot(length(setdiff(unique(labno_df$nhsno), gi_nhsnos)) == 0)

labno_query <- unique(labno_df$labno)

# Extract all results for lab numbers -------------------------------------

message("Extracting DNA Database results for all lab numbers")

dnadb_results <- results_tbl |> 
  filter(labno %in% labno_query) |> 
  select(labno, genodate, pcrid, test, genotype, genotype2,
         genocomm) |> 
  collect()

stopifnot(nrow(dnadb_results) != 0)

# Extract germline variant results ----------------------------------------

message("Finding germline variant DNA Database results")

icp_test_strings <- unique(grep(pattern = "hs2(\\s|)icp", 
            x = dnadb_results$test, 
            ignore.case = TRUE, 
            value = TRUE))

glvar_dnadb_results <- dnadb_results |> 
  filter(test %in% c(icp_test_strings,
                     "NGS SSXT ICP",
                     "ICP PANEL",
                     "SSXT ICP NGS",
                     "ICP SSXT NGS",
                     "Panel re-analysis of 24043064 from WS144546",
                     "NGS SSXT ICPv4",
                     "SSXTHS2 ICPv4",
                     "ICPv4 NGS SSXT HS2")) |> 
  left_join(labno_df, by = "labno") |> 
  relocate(nhsno)

# Extract tumour variant results ------------------------------------------

message("Finding tumour variant DNA Database results")

pansolid_test_strings <- unique(grep(pattern = "seq\\span", 
                                     x = dnadb_results$test, 
                                     ignore.case = TRUE, 
                                     value = TRUE))

tvar_dnadb_results <- dnadb_results |> 
  filter(test %in% c(pansolid_test_strings,
                     "NGS Pansolid", "NGS PanSolid QIAseq")) |> 
  left_join(labno_df, by = "labno") |> 
  relocate(nhsno)

# Check date range of data ------------------------------------------------

message(paste0("dnadb germline variant data ranges from ",
               format.Date(x = min(glvar_dnadb_results$genodate), 
                           format = "%d %B %Y"),
               " to ",
               format.Date(x = max(glvar_dnadb_results$genodate), 
                           format = "%d %B %Y")))
message(paste0("dnadb tumour variant data ranges from ",
               format.Date(x = min(tvar_dnadb_results$genodate), 
                           format = "%d %B %Y"),
               " to ",
               format.Date(x = max(tvar_dnadb_results$genodate), 
                           format = "%d %B %Y")))

# Export data -------------------------------------------------------------

message("Exporting DNA database tvar and glvar data")

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
