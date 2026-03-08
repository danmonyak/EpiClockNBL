#!/bin/bash

Rscript -e "rmarkdown::render('Data_Prep.Rmd', output_format = 'html_document', output_file = paste0('Data_Prep ', Sys.time(), '.html'))"