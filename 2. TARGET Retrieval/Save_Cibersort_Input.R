# Save the TPM gene expression data of the analysis cohort tumors as a tsv file
# This is the mixture file used as input for CIBERSORTx
#
# Run after Data_Processing_Pipeline.ipynb, from inside this directory:
#   Rscript Save_Cibersort_Input.R

library(jsonlite)

source('util.R')

save_dir <- "TARGET"
gene_filename <- "cohort1.rnaseq"
processed_clinical_filename <- "clinical.processed.tsv"

consts <- read_json(file.path(repo_dir, 'config.json'))

if (consts[["Windows"]]) {
  save_path <- file.path("P:", save_dir)
} else {
  save_path <- file.path(consts[['official_indir']], save_dir)
}

tpm <- readRDS(file.path(save_path, paste0(gene_filename, "_tpm.rds")))

clinical <- read.table(file.path(save_path, processed_clinical_filename),
                       sep = "\t", header = TRUE)
clinical$in_analysis_dataset <- as.logical(clinical$in_analysis_dataset)
analysis_tumors <- clinical$sampleID[clinical$in_analysis_dataset]

# keep the gene names and the tumors in the analysis cohort
is_tumor <- colnames(tpm) != "Gene"
keep <- !is_tumor | getTumorIDs(colnames(tpm)) %in% analysis_tumors
cat("Tumors with gene expression data:", sum(is_tumor), "\n")
cat("Tumors in the analysis cohort:", length(analysis_tumors), "\n")
cat("Tumors in the analysis cohort with gene expression data:", sum(keep & is_tumor), "\n")

tpm <- tpm[, keep]

outfile_path <- file.path(save_path, "cohort1.analysis_tumors.rnaseq_tpm.tsv")
write.table(tpm, file = outfile_path,
            sep = "\t", quote = F, row.names = F)
cat("Saved", outfile_path, "\n")
