library(tidyverse)

# Load initial data -------------------------------------------------------

glvar_dnadb_extracted <- read_csv(
  paste0(config::get("data_folderpath"),
         "01_initial/",
         "glvar_dnadb_extracted.csv"),
  show_col_types = FALSE)

tvar_dnadb_extracted <- read_csv(
  paste0(config::get("data_folderpath"),
         "01_initial/",
         "tvar_dnadb_extracted.csv"),
  show_col_types = FALSE)

# Annotate variants in germline results -----------------------------------

message("Adding variant classifications to germline variant DNA database results")

glvar_classifications <- read_csv(paste0(config::get("data_folderpath"),
                                         "01_initial/",
                                         "glvar_classifications.csv"),
                                  show_col_types = FALSE)

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
    # Separate HGVS nomenclature into glvar_hgvs_description column
    glvar_hgvs_description = str_extract(string = glvar_snv_result,
                                         pattern = "(.*)\\s\\d{1,3}%",
                                         group = 1),
    # Remove non-pathogenic BRCA1 variant from sample 24030433 which also
    # has a pathogenic BRCA2 variant
    glvar_hgvs_description = case_when(
      genotype == "BRCA2 c.4276dup p.(Thr1426AsnfsTer12) 48%; BRCA1 c.5074+13C>A 43%; CNV Analysis Failed" ~"BRCA2 c.4276dup p.(Thr1426AsnfsTer12)",
      TRUE ~glvar_hgvs_description
    ),
    # glvar_description is used in the iGene results to describe CNVs
    glvar_description = str_extract(string = glvar_cnv_result,
                                  pattern = ".*Ex.*")
    ) |> 
  left_join(glvar_classifications,
             by = "glvar_hgvs_description") |> 
  mutate(glvar_classification = case_when(
    glvar_description == " BRCA1 Ex02 deletion" ~"Pathogenic",
    TRUE ~glvar_classification
  ))

samples_with_gl_variants <- glvar_dnadb_classifications |> 
  filter(!is.na(glvar_hgvs_description) |
           !is.na(glvar_description))

# Check all variants have a classification
stopifnot(anyNA(samples_with_gl_variants$glvar_classification) == FALSE)

# Add headline to germline results ----------------------------------------

glvar_dnadb_cleaned <- glvar_dnadb_classifications |> 
  mutate(glvar_headline_result = case_when(
    genotype %in% c("No pathogenic variant identified; No relevant CNV identified",
                    "No pathogenic variant identified; CNV Analysis Failed") ~"No reportable variant(s) detected",
    genotype %in% c("Fail") ~"Analysis failed (see quality score)",
    glvar_classification %in% c("Pathogenic", "Likely pathogenic",
                                "Uncertain significance") ~"Reportable variant(s) detected",
    glvar_classification == "Not reported" ~"No reportable variant(s) detected",
    !is.na(glvar_description) ~"Reportable variant(s) detected",
    # Specify headline for 18028742
    (glvar_snv_result == "No pathogenic variant identified" &
      glvar_cnv_result == " No pathogenic variant identified") ~"No reportable variant(s) detected"))

stopifnot(anyNA(glvar_dnadb_cleaned$glvar_headline_result) == FALSE)

# Add classification to tumour variant results ----------------------------

message("Adding classifications to tumour variant DNA Database results")

tvar_classifications <- read_csv(paste0(config::get("data_folderpath"),
                                        "01_initial/",
                                        "tvar_classifications.csv"),
                                 show_col_types = FALSE)


tvar_dnadb_classifications <- tvar_dnadb_extracted |> 
  mutate(tvar_hgvs_description = str_extract(string = genotype,
                                             pattern = "(.*)\\s\\d{1,3}%",
                                             group = 1)) |> 
  left_join(tvar_classifications, by = "tvar_hgvs_description")

samples_with_tumour_variants <- tvar_dnadb_classifications |> 
  filter(!is.na(tvar_hgvs_description))

stopifnot(anyNA(samples_with_tumour_variants$tvar_classification) == FALSE)

# Add headline to tumour variant results ----------------------------------

no_path_var_strings <- unique(grep(pattern = "no\\spathogenic",
            x = tvar_dnadb_extracted$genotype,
            ignore.case = TRUE,
            value = TRUE))

fail_strings <- c("Analysis failed", "Fail")

tvar_dnadb_cleaned <- tvar_dnadb_classifications |> 
  mutate(tvar_headline_result = case_when(
    genotype %in% no_path_var_strings ~"No reportable variant(s) detected",
    genotype %in% fail_strings ~"Analysis failed (see quality score)",
    tvar_classification %in% c("Pathogenic", "Likely pathogenic",
                                "Uncertain significance") ~"Reportable variant(s) detected",
    tvar_classification == "Not reported" ~"No reportable variant(s) detected"))

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
