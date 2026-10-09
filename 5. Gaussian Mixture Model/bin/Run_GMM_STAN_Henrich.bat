@echo off
rem run from the parent directory, where GMM_STAN.Rmd is
cd /d "%~dp0.."
Rscript -e "rmarkdown::render('GMM_STAN.Rmd', output_format = 'html_document', output_file = paste0('GMM_STAN_', format(Sys.time(), '%%Y-%%m-%%d_%%H-%%M-%%S'), '.html'), params = list(dataset='Henrich', beta_vals_filename='Henrich.methyl.TARGET_Clock_sites.tsv'))"