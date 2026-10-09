#!/bin/bash
# run from the parent directory, where GMM_STAN.Rmd is
cd "$(dirname "$0")/.."
Rscript -e "rmarkdown::render('GMM_STAN.Rmd', output_format = 'html_document', output_file = paste0('GMM_STAN_', format(Sys.time(), '%Y-%m-%d_%H-%M-%S'), '.html'), params = list(dataset='TARGET', beta_vals_filename='cohort1.methyl.clock_sites.tsv'))"