# Add identifiers to SeqOne genomic instability files

library(tidyverse)

# Connect to SQL server ---------------------------------------------------

source(here::here("scripts/connect_to_sql_server.R"))

# Load collated csv data --------------------------------------------------

gi_csv_collated <- read_csv(file = paste0(
  config::get("data_folderpath"),
  "01_initial/",
  "gi_csv_collated.csv"),
  col_types = list(
    labno = col_character(),
    worksheet = col_character(),
    sample = col_character(),
    analysis_date = col_character(),
    date = col_datetime(),
    somahrd_version = col_character(),
    lga = col_integer(),
    lpc = col_integer(),
    score = col_number(),
    status = col_character(),
    brca_status = col_logical(),
    brca_mutation = col_logical(),
    ccne1_cn = col_number(),
    rad51b_cn = col_number(),
    coverage = col_number(),
    pct_mapped_reads = col_number(),
    pct_tum_cell = col_number(),
    gi_confidence = col_number(),
    low_tumor_fraction = col_number(),
    filepath = col_character()
  ))

# Get sample identifiers --------------------------------------------------

gi_labnos <- unique(gi_csv_collated$labno)

gi_sample_info <- sample_tbl |> 
  filter(labno %in% gi_labnos) |> 
  select(labno, firstname, surname, nhsno) |> 
  collect()

stopifnot(anyNA.data.frame(gi_sample_info |> 
                   select(-nhsno)) == FALSE)

# Add sample identifiers --------------------------------------------------

gi_csv_collated_patient_info <- gi_csv_collated |> 
  left_join(gi_sample_info, by = "labno") |> 
  relocate(labno, surname, firstname, nhsno)

stopifnot(anyNA(gi_csv_collated_patient_info$firstname) == FALSE)

# Annotate validation data ------------------------------------------------

DOC6192_validation_worksheets <- c("WS133557", "WS134687", "WS134928", 
                                   "WS135001", "WS135498")

DOC6255_validation_worksheets <- c("WS136827", "WS138201", "WS138439", 
                                   "WS138627")

DOC6588_validation_worksheets <- c("WS147582", "WS149085", "WS149086")

validation_worksheets <- c(DOC6192_validation_worksheets,
                           DOC6255_validation_worksheets,
                           DOC6588_validation_worksheets)

gi_csv_collated_validation_info <- gi_csv_collated_patient_info |> 
  mutate(service_validation = case_when(
    worksheet %in% validation_worksheets ~"validation",
    # 2 validation samples were on WS138061. The other samples on this
    # worksheet were live clinical samples
    (worksheet == "WS138061" & labno %in% c("23047082", "23053359")) ~"validation",
    TRUE ~"service"
  ),
  patient_non_patient = case_when(
    surname %in% c("Seraseq", "GenQA") ~"non-patient",
    TRUE ~"patient"
  ))

# Remove duplicates -------------------------------------------------------

gi_csv_cleaned_orpp <- gi_csv_collated_validation_info |> 
  filter(!is.na(nhsno) &
           service_validation == "service" &
           patient_non_patient == "patient") |> 
  # Some patient have multiple samples. To select samples with conclusive 
  # results (positive or negative statuses), arrange by NHS number and 
  # status and then
  # remove NHS number duplicates
  arrange(nhsno, status) |> 
  filter(!duplicated(nhsno))

stopifnot(unique(gi_csv_cleaned_orpp$patient_non_patient) == "patient")
stopifnot(unique(gi_csv_cleaned_orpp$service_validation) == "service")
stopifnot(anyDuplicated(gi_csv_cleaned_orpp$labno) == 0)
stopifnot(anyDuplicated(gi_csv_cleaned_orpp$nhsno) == 0)

gi_csv_cleaned_orps <- gi_csv_collated_validation_info |> 
  filter(service_validation == "service" &
           patient_non_patient == "patient") |> 
  arrange(labno, status) |> 
  filter(!duplicated(labno))

# Export cleaned data -----------------------------------------------------

write_csv(gi_csv_cleaned_orpp,
          paste0(config::get("data_folderpath"),
                 "02_cleaned/",
                 "gi_csv_cleaned_orpp.csv"))

write_csv(gi_csv_cleaned_orps,
          paste0(config::get("data_folderpath"),
                 "02_cleaned/",
                 "gi_csv_cleaned_orps.csv"))
