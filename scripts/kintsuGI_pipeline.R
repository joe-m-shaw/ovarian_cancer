# kintsuGI pipeline

message("Running kintsuGI pipeline")

source(here::here("scripts/gi_csv_collate.R"))

source(here::here("scripts/gi_csv_annotate.R"))
      
source(here::here("scripts/glvar_tvar_dnadb_extract.R"))

source(here::here("scripts/glvar_tvar_igene_extract.R"))

source(here::here("scripts/glvar_tvar_dnadb_clean.R"))

source(here::here("scripts/glvar_tvar_dnadb_igene_bind.R"))

source(here::here("scripts/gi_glvar_tvar_join.R"))

message("kintsuGI pipeline finished")
