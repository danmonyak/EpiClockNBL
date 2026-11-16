@echo off
rem Runs the "adjusted" GSEA analysis: preprocessing, GSEA, then the figure
rem run from the parent directory, where the GSEA Rmds are
cd /d "%~dp0.."
Rscript -e "rmarkdown::render('GSEA_Data_Preprocessing.Rmd', output_format = 'html_document', output_file = paste0('ADJUSTED_GSEA_Data_Preprocessing_', format(Sys.time(), '%%Y-%%m-%%d_%%H-%%M-%%S'), '.html'), params = list(analysis='adjusted'))" || exit /b 1
Rscript Run_GSEA.R adjusted || exit /b 1
Rscript -e "rmarkdown::render('GSEA_Figure.Rmd', output_format = 'html_document', output_file = paste0('ADJUSTED_GSEA_Figure_', format(Sys.time(), '%%Y-%%m-%%d_%%H-%%M-%%S'), '.html'), params = list(analysis='adjusted'))"