@echo off
rem Runs the "split1" GSEA analysis: preprocessing, GSEA, then the figure
rem run from the parent directory, where the GSEA Rmds are
cd /d "%~dp0.."
Rscript -e "rmarkdown::render('GSEA_Data_Preprocessing.Rmd', output_format = 'html_document', output_file = paste0('SPLIT1_GSEA_Data_Preprocessing_', format(Sys.time(), '%%Y-%%m-%%d_%%H-%%M-%%S'), '.html'), params = list(analysis='split1'))" || exit /b 1
Rscript Run_GSEA.R split1 || exit /b 1
Rscript -e "rmarkdown::render('GSEA_Figure.Rmd', output_format = 'html_document', output_file = paste0('SPLIT1_GSEA_Figure_', format(Sys.time(), '%%Y-%%m-%%d_%%H-%%M-%%S'), '.html'), params = list(analysis='split1', xlab_annotations_ypos=-5.3))"