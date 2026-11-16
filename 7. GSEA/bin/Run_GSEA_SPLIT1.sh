#!/bin/bash
# Runs the "split1" GSEA analysis: preprocessing, GSEA, then the figure
# run from the parent directory, where the GSEA Rmds are
set -e
cd "$(dirname "$0")/.."
Rscript -e "rmarkdown::render('GSEA_Data_Preprocessing.Rmd', output_format = 'html_document', output_file = paste0('SPLIT1_GSEA_Data_Preprocessing_', format(Sys.time(), '%Y-%m-%d_%H-%M-%S'), '.html'), params = list(analysis='split1'))"
Rscript Run_GSEA.R split1
Rscript -e "rmarkdown::render('GSEA_Figure.Rmd', output_format = 'html_document', output_file = paste0('SPLIT1_GSEA_Figure_', format(Sys.time(), '%Y-%m-%d_%H-%M-%S'), '.html'), params = list(analysis='split1'))"
