library(tidyverse)

# Load data ---------------------------------------------------------------

gi_csv_cleaned_orpp <- read_csv(paste0(config::get("data_folderpath"),
                                       "02_cleaned/",
                                       "gi_csv_cleaned_orpp.csv"),
                                col_types = list(
                                  nhsno = col_character()
                                ))

gi_csv_cleaned_orps <- read_csv(paste0(config::get("data_folderpath"),
                                       "02_cleaned/",
                                       "gi_csv_cleaned_orps.csv"),
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

message("Joining GI and germline variant data by NHS number")

gi_csv_cleaned_orpp_for_join <- gi_csv_cleaned_orpp |> 
  rename(gi_labno = labno,
         gi_pathno = pathno) |> 
  select(nhsno, gi_labno, gi_pathno, firstname, surname, lga, lpc, score,
         status, gi_confidence)

glvar_dnadb_igene_bound_orpp_for_join <- glvar_dnadb_igene_bound_orpp |> 
  rename(glvar_labno = labno,
         glvar_genotype = genotype)

gi_glvar_joined <- gi_csv_cleaned_orpp_for_join |> 
  inner_join(glvar_dnadb_igene_bound_orpp_for_join,
            by = "nhsno",
            relationship = "one-to-one") 

# Join GI and tvar data ---------------------------------------------------

message("Joining GI and tumour variant data by labno")

tvar_dnadb_igene_bound_orpp_for_join <- tvar_dnadb_igene_bound_orpp |> 
  rename(tvar_labno = labno,
         tvar_genotype = genotype,
         tvar_pathno = pathno)

gi_tvar_joined_by_labno <- gi_csv_cleaned_orpp_for_join |> 
  inner_join(tvar_dnadb_igene_bound_orpp_for_join,
             join_by("gi_labno" == "tvar_labno"),
             relationship = "one-to-one") 

gi_tvar_joined_by_nhsno <- gi_csv_cleaned_orpp_for_join |> 
  inner_join(tvar_dnadb_igene_bound_orpp_for_join,
             by = "nhsno",
             relationship = "one-to-one") 




# Join tvar and glvar data ------------------------------------------------

glvar_tvar_joined <- glvar_dnadb_igene_bound_orpp |> 
  inner_join(tvar_dnadb_igene_bound_orpp,
             by = "nhsno",
             relationship = "one-to-one")

# Export data -------------------------------------------------------------

message("Exporting joined data")

write_csv(gi_glvar_joined,
          paste0(config::get("data_folderpath"),
                 "03_joined/",
                 "gi_glvar_joined.csv"))

write_csv(gi_tvar_joined_by_labno,
          paste0(config::get("data_folderpath"),
                 "03_joined/",
                 "gi_tvar_joined_by_labno.csv"))

write_csv(gi_tvar_joined_by_nhsno,
          paste0(config::get("data_folderpath"),
                 "03_joined/",
                 "gi_tvar_joined_by_nhsno.csv"))

write_csv(glvar_tvar_joined,
          paste0(config::get("data_folderpath"),
                 "03_joined/",
                 "glvar_tvar_joined.csv"))
