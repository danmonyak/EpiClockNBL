#!/bin/bash

Rscript -e "rmarkdown::render('GMM_STAN.Rmd', output_format = 'html_document', output_file = paste0('GMM_STAN ', Sys.time(), '.html'), params = list(dataset='Henrich', beta_vals_filename='Henrich.methyl.TARGET_Clock_sites.tsv'))"