# Identify new SeqOne Genomic Instability csv files

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
  filter(grepl(pattern = "seqone|seq\\sone|seq_one|SSXT\\ssWGS\\sHRD",
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

# Issue 1: the same file may be saved in multiple folders, so filepath cannot 
# be used as a unique identifier
# Issue 2: the folder structure is not consistent. Some files are saved in a
# folder called the worksheet name, others are in sub-folders. This makes
# parsing filenames from filepaths with regex more complicated.

filename_regex <- regex(
  r"[
  /           # The last forward slash in the filepath
  ([^/]*      # Variable name before hrd-results but must not include forward slash
  hrd-results # All filenames contain hrd-results somewhere in the string
  .*          # Variable name after hrd-results
  \.csv)      # File type
  ]",
  comments = TRUE
)

gi_filepath_df <- tibble(
  filepath = unlist(gi_ws_filepaths)) |> 
  mutate(filename = str_extract(string = filepath,
                                pattern = filename_regex,
                                group = 1),
         filename_length = str_length(filename))

stopifnot(anyNA(gi_filepath_df$filename) == FALSE)
stopifnot(max(gi_filepath_df$filename_length) < 50)

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

if(nrow(gi_new_filepath_df) == 0) {
  message("No new files identified.")
}
if(nrow(gi_new_filepath_df) > 0) {
  
  message(paste0(length(gi_new_filepath_df$filepath),
                 " new files identified"))
  
  gi_filepaths_to_copy <- gi_new_filepath_df$filepath
  
  message("Copying csv files to directory")
  
  file.copy(from = gi_filepaths_to_copy,
            to = gi_csv_new_folder)
  
  source(here::here("scripts/gi_csv_collate.R"))
} 