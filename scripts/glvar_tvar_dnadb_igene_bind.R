library(tidyverse)

source(here::here("scripts/connect_to_sql_server.R"))

# Load results ------------------------------------------------------------

glvar_igene_extracted <- read_csv(file = 
                                    paste0(config::get("data_folderpath"),
                                           "01_initial/",
                                           "glvar_igene_extracted.csv"))

tvar_igene_extracted <- read_csv(file = 
                                    paste0(config::get("data_folderpath"),
                                           "01_initial/",
                                           "tvar_igene_extracted.csv"))

glvar_dnadb_cleaned <- read_csv(file = 
                                  paste0(config::get("data_folderpath"),
                                         "02_cleaned/",
                                         "glvar_dnadb_cleaned.csv"),
                                col_types = list(
                                  labno = col_character()
                                ))

tvar_dnadb_cleaned <- read_csv(file = 
                                  paste0(config::get("data_folderpath"),
                                         "02_cleaned/",
                                         "tvar_dnadb_cleaned.csv"),
                                col_types = list(
                                  labno = col_character()
                                ))

# Add patient identifiers to iGene data -----------------------------------

message("Adding patient identifiers before bind")

patient_id_df <- sample_tbl |> 
  select(nhsno, i_gene_r_no, labno) |> 
  collect() |> 
  rename(rno = i_gene_r_no)

patient_rno_nhsno_df <- patient_id_df |> 
  select(-labno) |> 
  filter(!duplicated(rno)) |> 
  filter((!is.na(rno) &
            !is.na(nhsno)))

glvar_igene_for_bind <- glvar_igene_extracted |>
  rename(rno = referral_number) |> 
  mutate(data_source = "igene",
         genotype = "",
         labno = "") |> 
  left_join(patient_rno_nhsno_df, by = "rno",
            relationship = "many-to-one") |> 
  select(nhsno, labno, rno, data_source, 
         glvar_headline_result, genotype, 
         glvar_panel_coverage, glvar_reflex_test,
         glvar_quality_score, glvar_zygosity,
         glvar_hgvs_description, glvar_classification,
         glvar_inheritance, glvar_genomic_coordinates,
         glvar_incidental_finding, glvar_description,
         glvar_copy_number_state, glvar_checker_comments)

stopifnot(nrow(glvar_igene_for_bind) == nrow(glvar_igene_extracted))

tvar_igene_for_bind <- tvar_igene_extracted |> 
  rename(rno = referral_number) |> 
  mutate(data_source = "igene",
         genotype = "",
         labno = "") |> 
  left_join(patient_rno_nhsno_df, by = "rno",
            relationship = "many-to-one") |> 
  select(nhsno, labno, rno, data_source, 
         tvar_headline_result, genotype, 
         tvar_reflex_test, tvar_hgvs_description, 
         tvar_failed_hotspots, tvar_vaf_percent,
         tvar_classification, tvar_quality_score,
         tvar_genomic_coordinates, tvar_panel_coverage,
         tvar_analyst_comments, tvar_checker_comments,
         tvar_incidental_finding, tvar_evidence)

stopifnot(nrow(tvar_igene_for_bind) == nrow(tvar_igene_extracted))

# Add patient identifiers to DNA database data ----------------------------

glvar_dnadb_for_bind <- glvar_dnadb_cleaned |> 
  select(-nhsno) |> 
  mutate(data_source = "dnadb") |> 
  left_join(patient_id_df, by = "labno",
            relationship = "many-to-one") |> 
  mutate(
    glvar_panel_coverage = "", 
    glvar_reflex_test = "",
    glvar_quality_score = "",
    glvar_zygosity = "",
    glvar_inheritance = "", 
    glvar_genomic_coordinates = "",
    glvar_incidental_finding = "", 
    glvar_copy_number_state = "", 
    glvar_checker_comments = "") |> 
  select(nhsno, labno, rno, data_source, 
         glvar_headline_result, genotype, 
         glvar_panel_coverage, glvar_reflex_test,
         glvar_quality_score, glvar_zygosity,
         glvar_hgvs_description, glvar_classification,
         glvar_inheritance, glvar_genomic_coordinates,
         glvar_incidental_finding, glvar_description,
         glvar_copy_number_state, glvar_checker_comments)

stopifnot(nrow(glvar_dnadb_for_bind) == nrow(glvar_dnadb_cleaned))

tvar_dnadb_for_bind <- tvar_dnadb_cleaned |> 
  select(-nhsno) |> 
  mutate(data_source = "dnadb") |> 
  left_join(patient_id_df, by = "labno",
            relationship = "many-to-one") |> 
  mutate(
    tvar_reflex_test = "",
    genotype = "", 
    tvar_reflex_test = "", 
    tvar_failed_hotspots = "",
    tvar_vaf_percent = "",
    tvar_classification = "", 
    tvar_quality_score = "",
    tvar_genomic_coordinates = "", 
    tvar_panel_coverage = "",
    tvar_analyst_comments = "", 
    tvar_checker_comments = "",
    tvar_incidental_finding = "", 
    tvar_evidence = "") |> 
  select(nhsno, labno, rno, data_source, 
         tvar_headline_result, genotype, 
         tvar_reflex_test, tvar_hgvs_description, 
         tvar_failed_hotspots, tvar_vaf_percent,
         tvar_classification, tvar_quality_score,
         tvar_genomic_coordinates, tvar_panel_coverage,
         tvar_analyst_comments, tvar_checker_comments,
         tvar_incidental_finding, tvar_evidence)

stopifnot(nrow(tvar_dnadb_for_bind) == nrow(tvar_dnadb_cleaned))

# Bind germline variant data ----------------------------------------------

message("Binding iGene and DNA Database data")

glvar_dnadb_igene_bound <- rbind(glvar_dnadb_for_bind,
                                    glvar_igene_for_bind)

glvar_dnadb_igene_bound_orpp <- glvar_dnadb_igene_bound |> 
  filter(!is.na(nhsno)) |> 
  mutate(glvar_headline_result = factor(glvar_headline_result,
                                        levels = c("Reportable variant(s) detected",
                                                   "No reportable variant(s) detected",
                                                   "Analysis failed (see quality score)"))) |> 
  arrange(nhsno, glvar_headline_result) |> 
  filter(!duplicated(nhsno))

anyNA(glvar_dnadb_igene_bound_orpp$nhsno)

# Some samples have multiple results with inconclusive results. 
# Check that the conclusive results have been selected for 3 samples.

stopifnot(nrow(glvar_dnadb_igene_bound_orpp |> 
                 filter(labno == "24024388" &
                          glvar_headline_result == "No reportable variant(s) detected")) == 1)

stopifnot(nrow(glvar_dnadb_igene_bound_orpp |> 
  filter(rno == "R24-1J8H" &
           glvar_headline_result == "No reportable variant(s) detected")) == 1)

stopifnot(nrow(glvar_dnadb_igene_bound_orpp |> 
                 filter(labno == "24030686" &
                          glvar_headline_result == "Reportable variant(s) detected")) == 1)

# Bind tumour variant data ------------------------------------------------

tvar_dnadb_igene_bound <- rbind(tvar_dnadb_for_bind,
                                tvar_igene_for_bind)

tvar_dnadb_igene_bound_orpp <- tvar_dnadb_igene_bound |> 
  filter(!is.na(nhsno)) |> 
  mutate(tvar_headline_result = factor(tvar_headline_result,
                                        levels = c("Reportable variant(s) detected",
                                                   "No reportable variant(s) detected",
                                                   "Analysis failed (see quality score)"))) |> 
  arrange(nhsno, tvar_headline_result) |> 
  filter(!duplicated(nhsno)) 

# Checks for samples with multiple results including inconclusive results

stopifnot(nrow(tvar_dnadb_igene_bound_orpp |> 
       filter(labno == "24009901" &
                tvar_headline_result == "No reportable variant(s) detected")) == 1)

stopifnot(nrow(tvar_dnadb_igene_bound_orpp |> 
                 filter(rno == "R24-1E83" &
                          tvar_headline_result == "Reportable variant(s) detected")) == 1)

stopifnot(nrow(tvar_dnadb_igene_bound_orpp |> 
                 filter(labno == "23060393" &
                          tvar_headline_result == "Reportable variant(s) detected")) == 1)

# Export data -------------------------------------------------------------

message("Exporting bound iGene and DNA Database results")

write_csv(glvar_dnadb_igene_bound_orpp,
          paste0(config::get("data_folderpath"),
                 "02_cleaned/",
                 "glvar_dnadb_igene_bound_orpp.csv"))

write_csv(tvar_dnadb_igene_bound_orpp,
          paste0(config::get("data_folderpath"),
                 "02_cleaned/",
                 "tvar_dnadb_igene_bound_orpp.csv"))
