# Clean results from DNA Database

library(tidyverse)

# Load initial data -------------------------------------------------------

glvar_dnadb_extracted <- read_csv(
  paste0(config::get("data_folderpath"),
         "01_initial/",
         "glvar_dnadb_extracted.csv"),
  show_col_types = FALSE)

stopifnot(ncol(glvar_dnadb_extracted) == 8)
stopifnot(nrow(glvar_dnadb_extracted) != 0)

tvar_dnadb_extracted <- read_csv(
  paste0(config::get("data_folderpath"),
         "01_initial/",
         "tvar_dnadb_extracted.csv"),
  show_col_types = FALSE)

stopifnot(ncol(tvar_dnadb_extracted) == 8)
stopifnot(nrow(tvar_dnadb_extracted) != 0)

# Annotate variants in germline results -----------------------------------

glvar_classifications <- read_csv(paste0(config::get("data_folderpath"),
                                         "01_initial/",
                                         "glvar_classifications.csv"),
                                  show_col_types = FALSE)

stopifnot(anyDuplicated(glvar_classifications$glvar_seqv1_description) == 0)

message("Adding ",
        nrow(glvar_classifications),
        " variant classifications to germline variant DNA database results")

glvar_dnadb_classifications <- glvar_dnadb_extracted |> 
  mutate(
    # Split the "genotype" field into SNV and CNV results
    glvar_snv_result = str_extract(string = genotype,
                                        pattern = "(.*);(.*)",
                                        group = 1),
    glvar_cnv_result = str_extract(string = genotype,
                                        pattern = "(.*);(.*)",
                                        group = 2),
    # Correct input for sample 18028742 - checked against report
    glvar_snv_result = case_when(
      glvar_snv_result == "BRCA1 Exon 3 3 Copies" ~"No pathogenic variant identified",
      TRUE ~glvar_snv_result
    ),
    # Separate HGVS nomenclature into glvar_seqv1_description column
    glvar_seqv1_description = str_extract(string = glvar_snv_result,
                                         pattern = "(.*)\\s\\d{1,3}%",
                                         group = 1),
    # Remove non-pathogenic BRCA1 variant from sample 24030433 which also
    # has a pathogenic BRCA2 variant
    glvar_seqv1_description = case_when(
      genotype == "BRCA2 c.4276dup p.(Thr1426AsnfsTer12) 48%; BRCA1 c.5074+13C>A 43%; CNV Analysis Failed" ~"BRCA2 c.4276dup p.(Thr1426AsnfsTer12)",
      TRUE ~glvar_seqv1_description
    ),
    # glvar_icnv1_description is used in the iGene results to describe CNVs
    glvar_icnv1_description = str_extract(string = glvar_cnv_result,
                                  pattern = ".*Ex.*")
    ) |> 
  left_join(glvar_classifications,
             by = "glvar_seqv1_description") |> 
  mutate(glvar_icnv1_classification = case_when(
    glvar_icnv1_description == " BRCA1 Ex02 deletion" ~"Pathogenic",
    TRUE ~NA
  ))

samples_with_gl_sn_variants <- glvar_dnadb_classifications |> 
  filter(!is.na(glvar_seqv1_description))

samples_with_gl_cn_variants <- glvar_dnadb_classifications |> 
  filter(!is.na(glvar_icnv1_description))

# Check all variants have a classification
stopifnot(anyNA(samples_with_gl_sn_variants$glvar_seqv1_classification) == FALSE)

stopifnot(anyNA(samples_with_gl_cn_variants$glvar_icnv1_classification) == FALSE)

# Add headline to germline results ----------------------------------------

message("Adding headline results to germline variant DNA Database results")

glvar_dnadb_cleaned <- glvar_dnadb_classifications |> 
  mutate(glvar_headline_result = case_when(
    genotype %in% c("No pathogenic variant identified; No relevant CNV identified",
                    "No pathogenic variant identified; CNV Analysis Failed") ~"No reportable variant(s) detected",
    genotype %in% c("Fail") ~"Analysis failed (see quality score)",
    glvar_seqv1_classification %in% c("Pathogenic", "Likely pathogenic",
                                "Uncertain significance") ~"Reportable variant(s) detected",
    glvar_seqv1_classification == "Not reported" ~"No reportable variant(s) detected",
    glvar_icnv1_classification == "Pathogenic" ~"Reportable variant(s) detected",
    # Specify headline for 18028742
    (glvar_snv_result == "No pathogenic variant identified" &
      glvar_cnv_result == " No pathogenic variant identified") ~"No reportable variant(s) detected"))

stopifnot(anyNA(glvar_dnadb_cleaned$glvar_headline_result) == FALSE)

# Add classification to tumour variant results ----------------------------

tvar_classifications <- read_csv(paste0(config::get("data_folderpath"),
                                        "01_initial/",
                                        "tvar_classifications.csv"),
                                 show_col_types = FALSE)

stopifnot(anyDuplicated(tvar_classifications$tvar_seqv1_description) == 0)

message("Adding ",
        nrow(tvar_classifications),
        " variant classifications to tumour variant DNA Database results")

tvar_dnadb_classifications <- tvar_dnadb_extracted |> 
  mutate(tvar_seqv1_description = str_extract(string = genotype,
                                             pattern = "(.*)\\s\\d{1,3}%",
                                             group = 1),
         # Add variant percentage
         tvar_seqv1_state = str_extract(string = genotype,
                                        pattern = ".*\\s(\\d{1,3})%",
                                        group = 1)) |> 
  left_join(tvar_classifications, by = "tvar_seqv1_description")

samples_with_tumour_variants <- tvar_dnadb_classifications |> 
  filter(!is.na(tvar_seqv1_description))

stopifnot(anyNA(samples_with_tumour_variants$tvar_seqv1_description) == FALSE)

# Add headline to tumour variant results ----------------------------------

message("Adding headline results to tumour variant DNA Database results")

no_path_var_strings <- unique(grep(pattern = "no\\spathogenic",
            x = tvar_dnadb_extracted$genotype,
            ignore.case = TRUE,
            value = TRUE))

fail_strings <- c("Analysis failed", "Fail")

tvar_dnadb_cleaned <- tvar_dnadb_classifications |> 
  mutate(tvar_headline_result = case_when(
    genotype %in% no_path_var_strings ~"No reportable variant(s) detected",
    genotype %in% fail_strings ~"Analysis failed (see quality score)",
    tvar_seqv1_classification %in% c("Pathogenic", "Likely pathogenic",
                                "Uncertain significance") ~"Reportable variant(s) detected",
    tvar_seqv1_classification == "Not reported" ~"No reportable variant(s) detected"))

stopifnot(anyNA(tvar_dnadb_cleaned$tvar_headline_result) == FALSE)

# Export results ----------------------------------------------------------

message("Exporting DNA database cleaned results")

write_csv(glvar_dnadb_cleaned,
          paste0(config::get("data_folderpath"),
                 "02_cleaned/",
                 "glvar_dnadb_cleaned.csv"))

write_csv(tvar_dnadb_cleaned,
          paste0(config::get("data_folderpath"),
                 "02_cleaned/",
                 "tvar_dnadb_cleaned.csv"))
