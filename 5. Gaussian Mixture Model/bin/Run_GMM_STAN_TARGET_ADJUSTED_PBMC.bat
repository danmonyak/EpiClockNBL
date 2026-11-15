@echo off
rem run from the parent directory, where GMM_STAN.Rmd is
cd /d "%~dp0.."
Rscript -e "rmarkdown::render('GMM_STAN.Rmd', output_format = 'html_document', output_file = paste0('ADJUSTED_PBMC_GMM_STAN_', format(Sys.time(), '%%Y-%%m-%%d_%%H-%%M-%%S'), '.html'), params = list(prefix='ADJUSTED_PBMC', dataset='TARGET', beta_vals_filename='cohort1.analysis_tumors.methyl_adjusted_PBMC.clock_sites.tsv'))"