# kintsuGI pipeline

quarto::quarto_render(here::here("scripts/gi_csv_identify.qmd"))

rm(list=ls())

quarto::quarto_render(here::here("scripts/gi_csv_annotate.qmd"))

rm(list=ls())

quarto::quarto_render(here::here("scripts/glvar_tvar_dnadb_extract.qmd"))
    
rm(list=ls())

quarto::quarto_render(here::here("scripts/glvar_tvar_igene_extract.qmd"))

rm(list=ls())

quarto::quarto_render(here::here("scripts/glvar_tvar_dnadb_clean.qmd"))

rm(list=ls())

quarto::quarto_render(here::here("scripts/glvar_tvar_dnadb_igene_bind.qmd"))

rm(list=ls())

quarto::quarto_render(here::here("scripts/gi_glvar_tvar_join.qmd"))

rm(list=ls())

quarto::quarto_render("documents/kintsuGI_report.qmd")
