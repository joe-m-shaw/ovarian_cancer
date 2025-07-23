# kintsuGI pipeline

library(here)

source(here("scripts/gi_csv_collate.R"))

source(here("scripts/gi_csv_annotate.R"))
      
source(here("scripts/glvar_tvar_dnadb_extract.R"))

source(here("scripts/glvar_tvar_igene_extract.R"))

source(here("scripts/glvar_tvar_dnadb_clean.R"))

source(here("scripts/glvar_tvar_dnadb_igene_bind.R"))

source(here("scripts/gi_glvar_tvar_join.R"))
