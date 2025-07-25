# kintsuGI: joining ovarian cancer results

## Purpose

The purpose of this repo is to extract, collate, clean, analyse and visualise genomic instability (GI) data and other genetic test results for ovarian cancer patients at the North West Genomic Laboratory Hub.
Since December 2023, the data for 3 different tests (GI testing, tumour variant testing and germline variant testing) has been saved in 3 different locations (2 databases: DNA Database and iGene; and individual csv files on a shared drive) in different formats.
The aim of this project is to bring all the data together to audit the ovarian cancer service.

## Documentation

Documentation is saved in the "documents" folder and includes information about the naming conventions I have used, and notes about the data.

## Data

Data is stored in an internal Microsoft SQL server, and results are exported to a data folder on the internal S drive.
**No data should be available on this public Github repo.** The file path for the data folder is specified in a `config.yml` file which I have not committed to Github.

## Analysis

The analysis can be performed by running the `kintsuGi_pipeline.R` script.
A final report is then rendered using Quarto (`kintsuGI_report.qmd`).

## What does kintsugi mean?

Kintsugi is the Japanese art of repairing broken pottery, literally meaning "golden joinery".
The philosophy of kintsugi posits that breaking and rejoining are "part of the history of an object, rather than something to disguise" ([Wikipedia](https://en.wikipedia.org/wiki/Kintsugi)).
This project aims to join data from multiple aspects of ovarian cancer testing, without trying to disguise this.
As the last two letters of kintsugi can also stand for "Genomic Instability" I have capitalised them in the repo name.

## 
