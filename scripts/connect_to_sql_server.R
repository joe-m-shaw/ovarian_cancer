
library(tidyverse)
library(odbc)
library(DBI)
library(dbplyr)

# Microsoft SQL server connection -----------------------------------------

dbi_con <- DBI::dbConnect(
  drv = odbc::odbc(),
  dsn = "moldb")

# iGene table -------------------------------------------------------------

eval_hrd_lazy <- tbl(dbi_con, dbplyr::in_catalog(catalog = "MolecularDB",
                                                 schema = "GenExport",
                                                 table = "Eval_HRD")) |> 
  janitor::clean_names() 

eval_hrd <- eval_hrd_lazy |> 
  # For some reason you have to select test_order_date first
  select(c(test_order_date, test_name, test_identifier, referral_number, 
           nhsn_umber, test_name, field, data_value, box)) |> 
  collect()

# DNA Database tables -----------------------------------------------------

sample_tbl <- tbl(dbi_con, dbplyr::in_catalog(catalog = "MolecularDB",
                                              schema = "dbo",
                                              table = "Samples")) |> 
  janitor::clean_names()

results_tbl <- tbl(dbi_con, 
                   dbplyr::in_catalog(
                     catalog = "MolecularDB",
                     schema = "dbo",
                     table = "ResultsAccess")) |> 
  janitor::clean_names() |> 
  rename(pcrid = resultsid)

dna_db_worksheets <- tbl(dbi_con, dbplyr::in_catalog(catalog = "MolecularDB",
                                                     schema = "dbo",
                                                     table = "PCR_New"))|> 
  janitor::clean_names()

dna_db_pcr_records <- tbl(dbi_con, dbplyr::in_catalog(catalog = "MolecularDB",
                                                      schema = "dbo",
                                                      table = "PCR_Records"))|> 
  janitor::clean_names()
