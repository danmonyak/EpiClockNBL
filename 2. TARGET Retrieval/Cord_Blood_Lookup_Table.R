# Builds the cord blood cell-type methylation lookup table, used by Beta_adjustment.Rmd
# when immune_reference is "CordBlood"
#
# The reference is FlowSorted.CordBloodCombined.450k (Gervin et al., 2019, Clinical Epigenetics 11:125):
# sorted cord blood cell types from four studies, measured on the 450K and EPIC arrays.
# For each cell type, the average beta value (after noob preprocessing) is computed at each CpG site
# of the GSE35069 lookup table, so the two tables have the same rows. Sites that are not in the
# cord blood data take the adult GSE35069 values.
#
# Run from inside "2. TARGET Retrieval":
#   Rscript Cord_Blood_Lookup_Table.R
# Requires the GSE35069 lookup table from WBC_Lookup_Table.ipynb.

suppressMessages({
  library(jsonlite)
  library(minfi)
  library(FlowSorted.Blood.EPIC)
})

source(file.path(repo_dir, 'save_util.R'))

consts <- read_json(file.path(repo_dir, 'config.json'))
other_dir <- file.path(consts[['official_indir']], 'Other')
adult_lookup_filepath <- file.path(other_dir, 'GSE35069_cell_type_mean_beta.tsv')
output_filepath <- file.path(other_dir, 'CordBloodCombined_cell_type_mean_beta.tsv')

# Cord blood cell types, named like the columns of the GSE35069 lookup table
# Nucleated red blood cells (nRBC) are left out: they are only found in cord blood
cell_type_names <- c(CD4T = "CD4+ T cells",
                     CD8T = "CD8+ T cells",
                     Mono = "CD14+ Monocytes",
                     Bcell = "CD19+ B cells",
                     NK = "CD56+ NK cells",
                     Gran = "Granulocytes")

# Raw data of the reference samples (downloaded from ExperimentHub the first time)
cat("Loading FlowSorted.CordBloodCombined.450k...\n")
rgset <- libraryDataGet('FlowSorted.CordBloodCombined.450k')
cat("Samples per cell type in the package:\n")
print(table(rgset$CellType))

rgset <- rgset[, rgset$CellType %in% names(cell_type_names)]

# Beta values after noob background correction and dye-bias normalization
cat("Preprocessing with noob...\n")
beta <- getBeta(preprocessNoob(rgset))

# Rows: the CpG sites of the GSE35069 lookup table, in the same order
adult_lookup <- read.table(adult_lookup_filepath, sep = "\t", header = TRUE, check.names = FALSE)
idx <- match(adult_lookup$ID_REF, rownames(beta))
cat("CpG sites of the GSE35069 table found in the cord blood data:",
    sum(!is.na(idx)), "of", length(idx), "\n")

# Average beta value of each cell type, using the samples measured at each site
lookup <- sapply(names(cell_type_names), function(cell_type) {
  rowMeans(beta[idx, rgset$CellType == cell_type, drop = FALSE], na.rm = TRUE)
})
lookup[is.nan(lookup)] <- NA
colnames(lookup) <- cell_type_names[colnames(lookup)]
lookup <- data.frame(ID_REF = adult_lookup$ID_REF, lookup, check.names = FALSE)

# How different the cord and adult references are at the fCpG sites
clock_cpgs <- read.table(file.path(repo_dir, "3. Select fCpGs", "outputs", "NBL_Clock_CpGs.txt"))$V1
at_clock <- match(clock_cpgs, lookup$ID_REF)
cat("fCpG sites with a cord blood value, for every cell type:",
    sum(complete.cases(lookup[at_clock, -1])), "of", length(clock_cpgs), "\n")
cat("Mean |cord - adult| beta value at the fCpG sites:\n")
for (cell_type in setdiff(colnames(lookup), "ID_REF")) {
  diff <- abs(lookup[at_clock, cell_type] - adult_lookup[at_clock, cell_type])
  cat(sprintf("  %-16s %.3f (sites differing by more than 0.1: %d)\n",
              cell_type, mean(diff, na.rm = TRUE), sum(diff > 0.1, na.rm = TRUE)))
}

# CpG sites that are not in the cord blood data (450K-only probes, since the combined data
# only has the probes shared by the 450K and EPIC arrays) take the adult GSE35069 values,
# so the cord blood reference only changes the sites it has data for
not_in_cord <- is.na(idx)
for (cell_type in setdiff(colnames(lookup), "ID_REF")) {
  lookup[not_in_cord, cell_type] <- adult_lookup[not_in_cord, cell_type]
}
cat("CpG sites given the adult GSE35069 values:", sum(not_in_cord),
    "(fCpG sites:", sum(not_in_cord[at_clock]), ")\n")

# Save, but don't overwrite an existing file that has different values
writeWithoutOverwrite(
  output_filepath,
  lookup,
  writeFunc = function(filepath, data) write.table(data, filepath, sep = "\t", quote = FALSE, row.names = FALSE),
  readFunc = function(filepath) read.table(filepath, sep = "\t", header = TRUE, check.names = FALSE),
  compareFunc = function(data, existing_data) {
    identical(data$ID_REF, existing_data$ID_REF) &&
      identical(colnames(data), colnames(existing_data)) &&
      isTRUE(all.equal(data[, -1], existing_data[, -1], tolerance = 1e-6, check.attributes = FALSE))
  }
)
