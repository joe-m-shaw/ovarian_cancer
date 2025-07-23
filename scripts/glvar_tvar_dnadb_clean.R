library(tidyverse)

# Load initial data -------------------------------------------------------

glvar_dnadb_extracted <- read_csv(
  paste0(config::get("data_folderpath"),
         "01_initial/",
         "glvar_dnadb_extracted.csv"))

tvar_dnadb_extracted <- read_csv(
  paste0(config::get("data_folderpath"),
         "01_initial/",
         "tvar_dnadb_extracted.csv"))


# Annotate variants in germline results -----------------------------------

glvar_classifications <- read_csv(paste0(config::get("data_folderpath"),
                                         "01_initial/",
                                         "glvar_classifications.csv"))

glvar_dnadb_classifications <- glvar_dnadb_extracted |> 
  mutate(
    # Split the "genotype" field into SNV and CNV results
    glvar_snv_result = str_extract(string = genotype,
                                        pattern = "(.*);(.*)",
                                        group = 1),
    glvar_cnv_result = str_extract(string = genotype,
                                        pattern = "(.*);(.*)",
                                        group = 2),
    # Separate HGVS nomenclature into glvar_hgvs_description column
    glvar_hgvs_description = str_extract(string = glvar_snv_result,
                                         pattern = "(.*)\\s\\d{1,3}%",
                                         group = 1),
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

# Add headline to germline results ----------------------------------------

glvar_dnadb_cleaned <- glvar_dnadb_classifications |> 
  mutate(glvar_headline_result = case_when(
    genotype %in% c("No pathogenic variant identified; No relevant CNV identified",
                    "No pathogenic variant identified; CNV Analysis Failed") ~"No reportable variant(s) detected",
    genotype %in% c("Fail") ~"Analysis failed (see quality score)",
    glvar_classification %in% c("Pathogenic", "Likely pathogenic",
                                "Uncertain significance") ~"Reportable variant(s) detected",
    glvar_classification == "Not reported" ~"No reportable variant(s) detected",
    !is.na(glvar_description) ~"Reportable variant(s) detected"))

stopifnot(anyNA(glvar_dnadb_cleaned$glvar_headline_result) == FALSE)

# Add headline to tumour variant results ----------------------------------

no_path_var_strings <- unique(grep(pattern = "no\\spathogenic",
            x = tvar_dnadb_extracted$genotype,
            ignore.case = TRUE,
            value = TRUE))

fail_strings <- c("Analysis failed", "Fail")

tvar_dnadb_cleaned <- tvar_dnadb_extracted |> 
  mutate(tvar_headline_result = case_when(
    genotype %in% no_path_var_strings ~"No reportable variant(s) detected",
    genotype %in% fail_strings ~"Analysis failed (see quality score)",
    TRUE ~"Reportable variant(s) detected"
  ))

# Export results ----------------------------------------------------------

write_csv(glvar_dnadb_cleaned,
          paste0(config::get("data_folderpath"),
                 "02_cleaned/",
                 "glvar_dnadb_cleaned.csv"))

write_csv(tvar_dnadb_cleaned,
          paste0(config::get("data_folderpath"),
                 "02_cleaned/",
                 "tvar_dnadb_cleaned.csv"))
