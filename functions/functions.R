
find_ws_filepaths <- function(worksheet,
                              path = config::get("ws_folderpath"),
                              pattern,
                              recursive = TRUE,
                              full_names = TRUE) {
  
  output <- list.files(str_c(path, worksheet, "/"),
                       full.names = full_names,
                       recursive = recursive,
                       pattern = pattern)
  
  return(output)
  
}

read_seqone_gi_csv <- function(file) {
  
  output <- readr::read_csv(file, 
                     n_max = 1,
                     col_types = list(
                       sample = col_character(),
                       analysis_date = col_character(),
                       somahrd_version = col_character(),
                       LGA = col_integer(),
                       LPC = col_integer(),
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
                       low_tumor_fraction = col_number()
                     )) |> 
    janitor::clean_names() |> 
    dplyr::mutate(filepath = file)
  
  return(output)
  
}
