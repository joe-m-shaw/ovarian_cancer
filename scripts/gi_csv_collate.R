# Collate new SeqOne Genomic Instability csv files

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
