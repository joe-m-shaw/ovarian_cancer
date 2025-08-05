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
    status = col_factor(levels = c("Positive",
                                   "Negative",
                                   "Non-conclusive"))
    ))

stopifnot(nrow(gi_csv_collated) != 0)
stopifnot(anyNA(gi_csv_collated$status) == 0)

# Get sample identifiers --------------------------------------------------

gi_labnos <- unique(gi_csv_collated$labno)

message("Adding identifiers for ",
        length(gi_labnos),
        " samples to GI csv data")

gi_sample_info <- sample_tbl |> 
  filter(labno %in% gi_labnos) |> 
  select(labno, i_gene_r_no, i_gene_s_no, 
         firstname, surname, nhsno, pathno) |> 
  collect()

stopifnot(anyNA.data.frame(gi_sample_info |> 
                   select(-nhsno)) == FALSE)

stopifnot(length(setdiff(gi_sample_info$labno, gi_labnos)) == 0)

stopifnot(anyDuplicated(gi_sample_info$labno) == 0)

# Add sample identifiers --------------------------------------------------

gi_csv_collated_patient_info <- gi_csv_collated |> 
  left_join(gi_sample_info, by = "labno",
            # One lab number can have multiple GI results, but should only
            # have one NHS number
            relationship = "many-to-one") |> 
  relocate(labno, i_gene_r_no, i_gene_s_no,
           surname, firstname, nhsno, pathno) |> 
  rename(rno = i_gene_r_no,
         sno = i_gene_s_no)

stopifnot(anyNA(gi_csv_collated_patient_info$firstname) == FALSE)

# Annotate validation data ------------------------------------------------

message("Annotating GI csv data with validation information")

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
  ),
  full_name = paste0(firstname, " ", surname))

stopifnot(anyNA(gi_csv_collated_validation_info$service_validation) == FALSE)
stopifnot(anyNA(gi_csv_collated_validation_info$patient_non_patient) == FALSE)

# Check name fields do not contain numbers
stopifnot(length(grep(pattern = "[[:digit:]]",
     gi_csv_collated_validation_info$full_name)) == 0)

# Filter to one result per patient ----------------------------------------

message("Filtering GI csv data to one result per patient")

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

# Check the lab number with a conclusive result has been selected for 
# patients with multiple samples including non-conclusive results
stopifnot(nrow(gi_csv_cleaned_orpp[gi_csv_cleaned_orpp$labno == "24044667" &
                           gi_csv_cleaned_orpp$status == "Negative",]) == 1)

stopifnot(nrow(gi_csv_cleaned_orpp[gi_csv_cleaned_orpp$labno == "24032495" &
                           gi_csv_cleaned_orpp$status == "Negative",]) == 1)

stopifnot(nrow(gi_csv_cleaned_orpp[gi_csv_cleaned_orpp$labno == "24045338" &
                          gi_csv_cleaned_orpp$status == "Negative",]) == 1)

stopifnot(nrow(gi_csv_cleaned_orpp[gi_csv_cleaned_orpp$labno == "24064756" &
                          gi_csv_cleaned_orpp$status == "Positive",]) == 1)

stopifnot(nrow(gi_csv_cleaned_orpp[gi_csv_cleaned_orpp$labno == "25030664" &
                          gi_csv_cleaned_orpp$status == "Positive",]) == 1)

if(length(anyDuplicated(gi_csv_cleaned_orpp$full_name)) != 0){
  message("Warning: there are patients with identical names and different NHS numbers")
} else{
  message("All patients have unique names")
}

# Filter to one result per sample -----------------------------------------

message("Filtering GI csv data to one result per sample")

gi_csv_cleaned_orps <- gi_csv_collated_validation_info |> 
  filter(service_validation == "service" &
           patient_non_patient == "patient") |> 
  arrange(labno, status) |> 
  filter(!duplicated(labno))

stopifnot(anyDuplicated(gi_csv_cleaned_orps$labno) == 0)

stopifnot(nrow(gi_csv_cleaned_orps[gi_csv_cleaned_orps$labno == "24045060" &
                      gi_csv_cleaned_orps$status == "Positive",]) == 1)

# Export cleaned data -----------------------------------------------------

message("Exporting cleaned GI csv data")

write_csv(gi_csv_collated_validation_info,
          paste0(config::get("data_folderpath"),
                 "02_cleaned/",
                 "gi_csv_cleaned.csv"))

write_csv(gi_csv_cleaned_orpp,
          paste0(config::get("data_folderpath"),
                 "02_cleaned/",
                 "gi_csv_cleaned_orpp.csv"))

write_csv(gi_csv_cleaned_orps,
          paste0(config::get("data_folderpath"),
                 "02_cleaned/",
                 "gi_csv_cleaned_orps.csv"))

