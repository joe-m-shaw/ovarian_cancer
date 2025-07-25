library(tidyverse)

# Load data ---------------------------------------------------------------

gi_csv_cleaned_orpp <- read_csv(paste0(config::get("data_folderpath"),
                                       "02_cleaned/",
                                       "gi_csv_cleaned_orpp.csv"),
                                col_types = list(
                                  nhsno = col_character()
                                ))

glvar_dnadb_igene_bound_orpp <- read_csv(paste0(config::get("data_folderpath"),
                                                "02_cleaned/",
                                                "glvar_dnadb_igene_bound_orpp.csv"),
                                         col_types = list(
                                           nhsno = col_character()
                                         ))

tvar_dnadb_igene_bound_orpp <- read_csv(paste0(config::get("data_folderpath"),
                                               "02_cleaned/",
                                               "tvar_dnadb_igene_bound_orpp.csv"),
                                        col_types = list(
                                          nhsno = col_character()
                                        ))

# Join GI and glvar data --------------------------------------------------

message("Joining GI and germline variant data")

gi_csv_cleaned_orpp_for_join <- gi_csv_cleaned_orpp |> 
  rename(gi_labno = labno) |> 
  select(nhsno, gi_labno, firstname, surname, lga, lpc, score,
         status, gi_confidence)

glvar_dnadb_igene_bound_orpp_for_join <- glvar_dnadb_igene_bound_orpp |> 
  rename(glvar_labno = labno,
         glvar_genotype = genotype) |> 
  select(nhsno, glvar_labno, glvar_headline_result,
         glvar_genotype, glvar_hgvs_description,
         glvar_classification, glvar_description, 
         gl_snv_gene, gl_cnv_gene, gl_gene)

gi_glvar_joined <- gi_csv_cleaned_orpp_for_join |> 
  inner_join(glvar_dnadb_igene_bound_orpp_for_join,
            by = "nhsno",
            relationship = "one-to-one") 

# Join GI and tvar data ---------------------------------------------------

message("Joining GI and tumour variant data")

tvar_dnadb_igene_bound_orpp_for_join <- tvar_dnadb_igene_bound_orpp |> 
  rename(tvar_labno = labno,
         tvar_genotype = genotype) |> 
  select(nhsno, tvar_labno, tvar_headline_result,
         tvar_genotype, tvar_hgvs_description, tvar_classification,
         tvar_gene)

gi_tvar_joined <- gi_csv_cleaned_orpp_for_join |> 
  inner_join(tvar_dnadb_igene_bound_orpp_for_join,
             by = "nhsno",
             relationship = "one-to-one") 

# Export data -------------------------------------------------------------

message("Exporting joined data")

write_csv(gi_glvar_joined,
          paste0(config::get("data_folderpath"),
                 "03_joined/",
                 "gi_glvar_joined.csv"))

write_csv(gi_tvar_joined,
          paste0(config::get("data_folderpath"),
                 "03_joined/",
                 "gi_tvar_joined.csv"))
