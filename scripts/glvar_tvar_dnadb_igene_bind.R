library(tidyverse)

source(here::here("scripts/connect_to_sql_server.R"))

# Load results ------------------------------------------------------------

glvar_igene_extracted <- read_csv(file = 
                                    paste0(config::get("data_folderpath"),
                                           "01_initial/",
                                           "glvar_igene_extracted.csv"),
                                  show_col_types = FALSE)

tvar_igene_extracted <- read_csv(file = 
                                    paste0(config::get("data_folderpath"),
                                           "01_initial/",
                                           "tvar_igene_extracted.csv"),
                                 show_col_types = FALSE)

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

# Prepare iGene data for bind ---------------------------------------------

glvar_igene_for_bind <- glvar_igene_extracted |>
  rename(rno = referral_number) |> 
  mutate(data_source = "igene",
         genotype = "",
         labno = "") |> 
  relocate(nhsno, labno, rno, data_source, .before = test_order_date)

stopifnot(nrow(glvar_igene_for_bind) == nrow(glvar_igene_extracted))

tvar_igene_for_bind <- tvar_igene_extracted |> 
  rename(rno = referral_number) |> 
  mutate(data_source = "igene",
         genotype = "",
         labno = "") |> 
  relocate(nhsno, labno, rno, data_source, .before = test_order_date)

stopifnot(nrow(tvar_igene_for_bind) == nrow(tvar_igene_extracted))

# Prepare DNA Database germline data for bind -----------------------------

patient_id_df <- sample_tbl |> 
  select(i_gene_r_no, labno) |> 
  collect() |> 
  rename(rno = i_gene_r_no)

glvar_dnadb_for_bind <- glvar_dnadb_cleaned |> 
  left_join(patient_id_df, by = "labno",
            relationship = "many-to-one") |> 
  mutate(data_source = "dnadb",
         # Add empty columns for bind
         test_order_date = NA,
         test_identifier = NA,
         glvar_seqv1_state = NA,
         glvar_seqv1_genomic_coordinates = NA,
         glvar_seqv1_evidence = NA,
         glvar_reflex_test = NA,
         glvar_analyst_comments = NA,
         glvar_incidental_finding = NA,
         glvar_panel_coverage = NA,
         glvar_seqv1_inheritance = NA,
         glvar_quality_score = NA,
         glvar_checker_comments = NA,
         glvar_icnv1_evidence = NA,
         glvar_icnv1_state = NA,
         glvar_icnv1_inheritance = NA,
         glvar_icnv1_genomic_coordinates = NA
         ) |> 
  select(
    # Identifiers
    nhsno, labno, rno, 
    # Test information
    data_source, test_order_date, test_identifier, 
    # Result information
    glvar_seqv1_state, glvar_seqv1_genomic_coordinates, glvar_seqv1_evidence, 
    glvar_reflex_test, glvar_analyst_comments, glvar_headline_result, 
    glvar_incidental_finding, glvar_panel_coverage, 
    glvar_seqv1_inheritance, glvar_seqv1_description, 
    glvar_seqv1_classification, glvar_quality_score, 
    glvar_checker_comments, glvar_icnv1_classification, 
    glvar_icnv1_evidence, glvar_icnv1_description, glvar_icnv1_state, 
    glvar_icnv1_inheritance, glvar_icnv1_genomic_coordinates, genotype)

stopifnot(nrow(glvar_dnadb_for_bind) == nrow(glvar_dnadb_cleaned))

# Prepare DNA Database tumour data for bind -------------------------------

tvar_dnadb_for_bind <- tvar_dnadb_cleaned |> 
  left_join(patient_id_df, by = "labno",
            relationship = "many-to-one") |> 
  mutate(data_source = "dnadb",
         test_order_date = NA, 
         test_identifier = NA, 
         tvar_analyst_comments = NA, 
         tvar_failed_hotspots = NA, 
         tvar_panel_coverage = NA, 
         tvar_reflex_test = NA, 
         tvar_quality_score = NA, 
         tvar_seqv1_genomic_coordinates = NA, 
         tvar_seqv1_state = NA, 
         tvar_checker_comments = NA, 
         tvar_icnv1_classification = NA, 
         tvar_icnv1_description = NA, 
         tvar_icnv1_state = NA, 
         tvar_icnv1_genomic_coordinates = NA, 
         tvar_seqv2_description = NA, 
         tvar_seqv2_classification = NA, 
         tvar_seqv2_state = NA, 
         tvar_seqv2_genomic_coordinates = NA, 
         tvar_incidental_finding = NA, 
         tvar_seqv1_evidence = NA,
         tvar_sv1_classification = NA,
         tvar_sv1_evidence = NA,          
         tvar_sv1_state = NA,
         tvar_sv1_description = NA,       
         tvar_sv1_genomic_coordinates = NA) |> 
  select(nhsno, labno, rno, data_source, test_order_date, 
         test_identifier, tvar_panel_coverage, tvar_failed_hotspots, 
         tvar_quality_score, tvar_headline_result, tvar_reflex_test, 
         tvar_analyst_comments, tvar_checker_comments, tvar_seqv2_state, 
         tvar_seqv1_state, tvar_seqv1_classification, tvar_seqv2_description,
         tvar_seqv2_classification, tvar_seqv1_description, 
         tvar_seqv1_genomic_coordinates, tvar_icnv1_genomic_coordinates, 
         tvar_icnv1_description, tvar_icnv1_classification, 
         tvar_icnv1_state, tvar_incidental_finding, 
         tvar_seqv2_genomic_coordinates, tvar_seqv1_evidence, 
         tvar_sv1_classification, tvar_sv1_evidence, tvar_sv1_state, 
         tvar_sv1_description, tvar_sv1_genomic_coordinates, genotype)

stopifnot(nrow(tvar_dnadb_for_bind) == nrow(tvar_dnadb_cleaned))

# Bind germline variant data ----------------------------------------------

message("Binding germline iGene and DNA Database data")

glvar_dnadb_igene_bound <- rbind(glvar_dnadb_for_bind,
                                    glvar_igene_for_bind)

# Annotate germline genes -------------------------------------------------

message("Annotating germline gene information")

gl_genes <- c("BRCA1", "BRCA2", "BRIP1", "PALB2", "RAD51D", 
              "MSH2", "MSH6", "CHEK2")

gl_gene_regex <- paste0(".*(", paste0(gl_genes, collapse = "|"), ").*")

glvar_dnadb_igene_bound_genes <- glvar_dnadb_igene_bound |> 
  mutate(gl_snv_gene = str_extract(string = glvar_seqv1_description,
                                 pattern = gl_gene_regex,
                                 group = 1),
       gl_cnv_gene = str_extract(string = glvar_icnv1_description,
                                 pattern = gl_gene_regex,
                                 group = 1),
       gl_gene = case_when(
         is.na(gl_snv_gene) & !is.na(gl_cnv_gene) ~gl_cnv_gene,
         !is.na(gl_snv_gene) & is.na(gl_cnv_gene) ~gl_snv_gene
       ),
       glvar_seqv1_classification = case_when(
         # Likely pathogenic reduced penetrance
         glvar_seqv1_description == "NM_000059.3(BRCA2):c.9302T>G p.(Leu3101Arg)" ~"Likely pathogenic",
         TRUE ~glvar_seqv1_classification
       ),
       # Add column to summarise SNV and CNV classifications
       glvar_headline_classification = case_when(
         !is.na(glvar_seqv1_classification) ~glvar_seqv1_classification,
         !is.na(glvar_icnv1_classification) ~glvar_icnv1_classification,
         TRUE ~NA))

samples_with_gl_variants <- glvar_dnadb_igene_bound_genes |> 
  filter(glvar_headline_result == "Reportable variant(s) detected")

stopifnot(anyNA(samples_with_gl_variants$gl_gene) == FALSE)
stopifnot(anyNA(samples_with_gl_variants$glvar_headline_classification) == FALSE)

# Filter germline data to one result per patient --------------------------

glvar_dnadb_igene_bound_orpp <- glvar_dnadb_igene_bound_genes |> 
  filter(!is.na(nhsno)) |> 
  mutate(glvar_headline_result = factor(glvar_headline_result,
                                        levels = c("Reportable variant(s) detected",
                                                   "No reportable variant(s) detected",
                                                   "Analysis failed (see quality score)"))) |> 
  arrange(nhsno, glvar_headline_result) |> 
  filter(!duplicated(nhsno))

stopifnot(anyNA(glvar_dnadb_igene_bound_orpp$nhsno) == FALSE)

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

message("Binding tumour iGene and DNA Database data")

tvar_dnadb_igene_bound <- rbind(tvar_dnadb_for_bind,
                                tvar_igene_for_bind)

# Annotate tumour variant genes -------------------------------------------

message("Annotating tumour variant gene information")

tvar_genes <- c("BRCA1", "BRCA2", "KRAS", "BRAF", "ERBB2")

tvar_gene_regex <- paste0(".*(", paste0(tvar_genes, collapse = "|"), ").*")

tvar_dnadb_igene_bound_genes <- tvar_dnadb_igene_bound |> 
  mutate(
    tvar_gene = str_extract(string = tvar_seqv1_description,
                            pattern = tvar_gene_regex,
                            group = 1),
    tvar_headline_classification = case_when(
      !is.na(tvar_seqv1_classification) ~tvar_seqv1_classification,
      !is.na(tvar_icnv1_classification) ~tvar_icnv1_classification,
      TRUE ~NA
    ))

samples_with_t_variants <- tvar_dnadb_igene_bound_genes |> 
  filter(!is.na(tvar_seqv1_description))

stopifnot(anyNA(samples_with_t_variants$tvar_gene) == FALSE)

# Filter tumour variant data to one result per patient --------------------

tvar_dnadb_igene_bound_orpp <- tvar_dnadb_igene_bound_genes |> 
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
