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

patient_id_df <- sample_tbl |> 
  select(nhsno, i_gene_r_no, labno) |> 
  collect() |> 
  rename(rno = i_gene_r_no)

glvar_igene_for_collation <- glvar_igene_extracted |>
  rename(rno = referral_number) |> 
  mutate(data_source = "igene",
         genotype = "") |> 
  left_join(patient_id_df, by = "rno") |> 
  select(nhsno, labno, rno, data_source, 
         glvar_headline_result, genotype, 
         glvar_panel_coverage, glvar_reflex_test,
         glvar_quality_score, glvar_zygosity,
         glvar_hgvs_description, glvar_classification,
         glvar_inheritance, glvar_genomic_coordinates,
         glvar_incidental_finding, glvar_description,
         glvar_copy_number_state, glvar_checker_comments)

tvar_igene_for_bind <- tvar_igene_extracted |> 
  rename(rno = referral_number) |> 
  mutate(data_source = "igene",
         genotype = "") |> 
  left_join(patient_id_df, by = "rno") |> 
  select(nhsno, labno, rno, data_source, 
         tvar_headline_result, genotype, 
         tvar_reflex_test, tvar_hgvs_description, 
         tvar_failed_hotspots, tvar_vaf_percent,
         tvar_classification, tvar_quality_score,
         tvar_genomic_coordinates, tvar_panel_coverage,
         tvar_analyst_comments, tvar_checker_comments,
         tvar_incidental_finding, tvar_evidence)

# Add patient identifiers to DNA database data ----------------------------

glvar_dnadb_for_collation <- glvar_dnadb_cleaned |> 
  select(-nhsno) |> 
  mutate(data_source = "dnadb") |> 
  left_join(patient_id_df, by = "labno") |> 
  mutate(
    glvar_panel_coverage = "", 
    glvar_reflex_test = "",
    glvar_quality_score = "",
    glvar_zygosity = "",
    glvar_hgvs_description = "", 
    glvar_classification = "",
    glvar_inheritance = "", 
    glvar_genomic_coordinates = "",
    glvar_incidental_finding = "", 
    glvar_description = "",
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

tvar_dnadb_for_bind <- tvar_dnadb_cleaned |> 
  select(-nhsno) |> 
  mutate(data_source = "dnadb") |> 
  left_join(patient_id_df, by = "labno") |> 
  mutate(
    tvar_reflex_test = "",
    genotype = "", 
    tvar_reflex_test = "", 
    tvar_hgvs_description = "", 
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

# Join germline variant data ----------------------------------------------

glvar_dnadb_igene_bound <- rbind(glvar_dnadb_for_collation,
                                    glvar_igene_for_collation)

glvar_dnadb_igene_bound_orpp <- glvar_dnadb_igene_bound |> 
  mutate(glvar_headline_result = factor(glvar_headline_result,
                                        levels = c("Reportable variant(s) detected",
                                                   "No reportable variant(s) detected",
                                                   "Analysis failed (see quality score)"))) |> 
  arrange(nhsno, glvar_headline_result) |> 
  filter(!duplicated(nhsno))

# Join tumour variant data ------------------------------------------------

tvar_dnadb_igene_bound <- rbind(tvar_dnadb_for_bind,
                                tvar_igene_for_bind)

tvar_dnadb_igene_bound_orpp <- tvar_dnadb_igene_bound |> 
  mutate(tvar_headline_result = factor(tvar_headline_result,
                                        levels = c("Reportable variant(s) detected",
                                                   "No reportable variant(s) detected",
                                                   "Analysis failed (see quality score)"))) |> 
  arrange(nhsno, tvar_headline_result) |> 
  filter(!duplicated(nhsno))

# Export data -------------------------------------------------------------

write_csv(glvar_dnadb_igene_bound_orpp,
          paste0(config::get("data_folderpath"),
                 "02_cleaned/",
                 "glvar_dnadb_igene_bound_orpp.csv"))

write_csv(tvar_dnadb_igene_bound_orpp,
          paste0(config::get("data_folderpath"),
                 "02_cleaned/",
                 "tvar_dnadb_igene_bound_orpp.csv"))
