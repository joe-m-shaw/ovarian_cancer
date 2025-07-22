# Collate SeqOne Genomic Instability csv files

# DNA Database connection -------------------------------------------------

message("Connecting to DNA database")

source(here::here("scripts/connect_to_sql_server.R"))
source(here::here("functions/functions.R"))

# Find SeqOne GI worksheets -----------------------------------------------

all_worksheets <- dna_db_worksheets |> 
  select(pcrid, date, description) |> 
  collect() |> 
  mutate(ws = paste0("WS", pcrid))

stopifnot(nrow(all_worksheets) > 0)

gi_ws_info <- all_worksheets |> 
  filter(grepl(pattern = "seqone|seq\\sone|seq_one",
               x = description,
               ignore.case = TRUE))

stopifnot(nrow(gi_ws_info) > 0)

# Find csvs in worksheet folders ------------------------------------------

message("Locating SeqOne csv files")

gi_ws_list <- list(gi_ws_info$ws)

gi_ws_filepaths <- gi_ws_list |> 
  map(\(gi_ws_list) find_ws_filepaths(worksheet = gi_ws_list,
                                      pattern = "hrd-results.*csv")) |> 
  flatten()

gi_filepath_df <- tibble(
  filepath = unlist(gi_ws_filepaths)) |> 
  mutate(filename = str_extract(string = filepath,
                                pattern = "WS\\d{6}/(.*hrd-results.*\\.csv)",
                                group = 1))

stopifnot(anyNA(gi_filepath_df$filename) == FALSE)

# Identify new files ------------------------------------------------------

initial_folder <-  paste0(config::get("data_folderpath"), 
                          "01_initial/")

gi_csv_archive_folder <- paste0(initial_folder,
                                "gi_csv_archive")

gi_csv_new_folder <- paste0(initial_folder,
                            "gi_csv_new")

if(length(list.files(gi_csv_new_folder)) != 0){
  stop("New data folder is not empty")
} else {
  message("New data folder is empty")
}

gi_archive_file_df <- tibble(
  filename = list.files(gi_csv_archive_folder,
                                full.names = FALSE))

gi_new_filepath_df <- gi_filepath_df |> 
  filter(!filename %in% gi_archive_file_df$filename)

if(length(gi_new_filepath_df) > 0) {
  message(paste0(length(gi_new_filepath_df$filepath),
                 " new files identified"))
} else {
  stop("No new files identified")
}

gi_filepaths_to_copy <- gi_new_filepath_df$filepath

# Copy to data folder -----------------------------------------------------

message("Copying csv files to directory")

file.copy(from = gi_filepaths_to_copy,
          to = gi_csv_new_folder)

# Collate csv files -------------------------------------------------------

message("Collating csv files")

new_filepaths_to_collate <- list.files(
  gi_csv_new_folder,
  full.names = TRUE,
  recursive = FALSE,
  pattern = "hrd-results.*csv")

gi_collated_new <- new_filepaths_to_collate |> 
  map(\(new_filepaths_to_collate) read_seqone_gi_csv(new_filepaths_to_collate)) |> 
  list_rbind() |> 
  mutate(date = parse_date_time(x = analysis_date, 
                                orders = c("dmy", "ymd")),
         worksheet = str_extract(string = sample,
                                 pattern = "(WS[0-9]{6})_(\\d{8})",
                                 group = 1),
         labno = str_extract(string = sample,
                             pattern = "(WS[0-9]{6})_(\\d{8})",
                             group = 2),
         status = factor(status,
                         levels = c("Positive", "Negative",
                                    "Non-conclusive"))) |> 
  relocate(date, .after = analysis_date) |> 
  relocate(labno, .before = sample) |> 
  relocate(worksheet, .after = labno)

stopifnot(ncol(gi_collated_new) == 20)

# Load existing collated data ---------------------------------------------

gi_collated_old <- read_csv(paste0(config::get("data_folderpath"), 
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

# Add new data ------------------------------------------------------------

gi_collated_updated <- rbind(gi_collated_old,
                             gi_collated_new)

# Check data --------------------------------------------------------------

message("Checking data")

if(anyNA.data.frame(gi_collated_updated |> 
                    select(-c(brca_status, brca_mutation)))){
  warning("NAs in gi_collated_updated")
}

if((min(gi_collated_updated$lga, na.rm = TRUE) < 0) |
   (max(gi_collated_updated$lga, na.rm = TRUE) > 100)){
  stop("LGA should be in range 0-100")
} else {
  message("LGA check passed")
}

if((min(gi_collated_updated$lpc, na.rm = TRUE) < 0)|
   (max(gi_collated_updated$lpc, na.rm = TRUE) > 100)){
  stop("LPC should be in range 0-100")
} else {
  message("LPC check passed")
}

if((min(gi_collated_updated$score, na.rm = TRUE) < 0)|
   (max(gi_collated_updated$score, na.rm = TRUE) > 1)){
  stop("SeqOne score should be between -10")
} else{
  message("SeqOne score check passed")
}

stopifnot(levels(gi_collated_updated$status) == 
            c("Positive", "Negative", "Non-conclusive"))

if(min(gi_collated_updated$ccne1_cn, na.rm = TRUE) < 0 |
   max(gi_collated_updated$ccne1_cn, na.rm = TRUE) > 100){
  stop("Check CCNE1 copy number")
} else {
  message("CCNE1 copy number check passed")
}

if(min(gi_collated_updated$rad51b_cn, na.rm = TRUE) < 0 |
   max(gi_collated_updated$rad51b_cn, na.rm = TRUE) > 100){
  stop("Check RAD51B copy number")
} else {
  message("RAD51B copy number check passed")
}

if(min(gi_collated_updated$coverage) < 0 |
   max(gi_collated_updated$coverage) > 7){
  stop("Check coverage")
} else {
  message("Coverage check passed")
}

if(min(gi_collated_updated$pct_tum_cell) < 0 |
   max(gi_collated_updated$pct_tum_cell) > 1){
  stop("Check percentage tumour cells")
} else {
  message("Percentage tumour cells check passed")
}

if(min(gi_collated_updated$gi_confidence, na.rm = TRUE) < 0 |
   max(gi_collated_updated$gi_confidence, na.rm = TRUE) > 1){
  stop("Check GI confidence")
} else {
  message("GI confidence check passed")
}

message("Data check complete")

# Export data -------------------------------------------------------------

write_csv(gi_collated_updated,
          paste0(initial_folder, 
                 "gi_csv_collated.csv"))

# Copy new data to archive ------------------------------------------------

file.copy(from = new_filepaths_to_collate,
          to = gi_csv_archive_folder)

# Delete files from new data folder ---------------------------------------

file.remove(new_filepaths_to_collate)

if(length(list.files(gi_csv_new_folder)) == 0){
  message("New file folder is empty")
}
