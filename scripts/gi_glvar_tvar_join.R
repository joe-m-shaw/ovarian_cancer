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

# Join data ---------------------------------------------------------------

gi_csv_cleaned_orpp_for_join <- gi_csv_cleaned_orpp |> 
  rename(gi_labno = labno) |> 
  select(nhsno, gi_labno, firstname, surname, lga, lpc, score,
         status, gi_confidence)

glvar_dnadb_igene_bound_orpp_for_join <- glvar_dnadb_igene_bound_orpp |> 
  rename(glvar_labno = labno,
         glvar_genotype = genotype) |> 
  select(nhsno, glvar_labno, glvar_headline_result,
         glvar_genotype, glvar_hgvs_description,
         glvar_classification, glvar_description)

tvar_dnadb_igene_bound_orpp_for_join <- tvar_dnadb_igene_bound_orpp |> 
  rename(tvar_labno = labno,
         tvar_genotype = genotype) |> 
  select(nhsno, tvar_labno, tvar_headline_result,
         tvar_genotype, tvar_hgvs_description)

gi_glvar_joined <- gi_csv_cleaned_orpp_for_join |> 
  filter(!is.na(nhsno)) |> 
  inner_join(glvar_dnadb_igene_bound_orpp_for_join,
            by = "nhsno") 

# Export data -------------------------------------------------------------

write_csv(gi_glvar_joined,
          paste0(config::get("data_folderpath"),
                 "03_joined/",
                 "gi_glvar_joined.csv"))
