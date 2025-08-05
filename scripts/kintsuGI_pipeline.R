# kintsuGI pipeline

message("Running kintsuGI pipeline")

source(here::here("scripts/gi_csv_identify.R"))

rm(list=ls())

source(here::here("scripts/gi_csv_annotate.R"))

rm(list=ls())
      
source(here::here("scripts/glvar_tvar_dnadb_extract.R"))

rm(list=ls())

source(here::here("scripts/glvar_tvar_igene_extract.R"))

rm(list=ls())

source(here::here("scripts/glvar_tvar_dnadb_clean.R"))

rm(list=ls())

source(here::here("scripts/glvar_tvar_dnadb_igene_bind.R"))

rm(list=ls())

source(here::here("scripts/gi_glvar_tvar_join.R"))

rm(list=ls())

message("Data processing complete. Rendering kintsuGI_report.qmd")

quarto::quarto_render("documents/kintsuGI_report.qmd")

message("kintsuGI pipeline finished")
