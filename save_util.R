# Helper functions for saving output files without overwriting existing data
# Source from an Rmd with: source(file.path(repo_dir, 'save_util.R'))

# Stop right away if any of the output files already exist
# Use this at the start of a long-running Rmd, so that it fails before the long computation
# instead of when it tries to save at the end
#   filepaths: paths of the output files
#   advice: extra text added to the error message (for example, how to rerun anyway)
stopIfOutputsExist <- function(filepaths, advice = NULL) {
  existing <- filepaths[file.exists(filepaths)]
  if (length(existing) > 0) {
    stop(paste0(
      "These output files already exist:\n",
      paste0("  ", existing, collapse = "\n"),
      "\nDelete them to run this again",
      if (is.null(advice)) "." else paste0(", ", advice)
    ), call. = FALSE)
  }
  invisible(NULL)
}

# Save a file, but don't overwrite an existing file at the same path with different data
# (same behavior as writeWithoutOverwrite in EpiClockNBL/util.py)
#   filepath: path of the file to save
#   data: object to save
#   writeFunc: function(filepath, data) that saves data to filepath
#   readFunc: function(filepath) that reads an existing file back in
#   compareFunc: function(data, existing_data) that returns TRUE if they match;
#                if NULL, all.equal() is used
#   verbose: print what happened
# If the file exists and the data differs, an error of class "overwriteError" is raised.
# Delete the old file to save the new data.
writeWithoutOverwrite <- function(filepath, data, writeFunc, readFunc,
                                  compareFunc = NULL, verbose = TRUE) {
  if (file.exists(filepath)) {
    existing_data <- readFunc(filepath)
    if (is.null(compareFunc)) {
      compareFunc <- function(data, existing_data) isTRUE(all.equal(data, existing_data))
    }
    comparison <- tryCatch(isTRUE(compareFunc(data, existing_data)),
                           error = function(e) FALSE)
    if (!comparison) {
      stop(structure(
        class = c("overwriteError", "error", "condition"),
        list(message = paste("Current output would overwrite alternative data at", filepath),
             call = sys.call(-1),
             filepath = filepath)
      ))
    }
    if (verbose) cat("Current output matches existing file.\n")
  } else {
    writeFunc(filepath, data)
    if (verbose) cat("Saved new file.\n")
  }
  if (verbose) cat(filepath, "\n")
  invisible(NULL)
}

# Save a table as one tab of the supplementary data workbook, keeping its other tabs
# (same behavior as writeSupplementaryTab in EpiClockNBL/util.py, which writes the other tabs)
#   sheet_name: name of the tab, for example "Data S2"
#   data: data frame to save (without row names)
#   verbose: print what happened
# If the tab exists with different contents, an error of class "overwriteError" is raised.
writeSupplementaryTab <- function(sheet_name, data, verbose = TRUE) {
  filepath <- file.path(repo_dir, "Supplementary Data", "Supplementary_Data.xlsx")

  # Compare as text, since the workbook is read back as text
  data <- as.data.frame(lapply(data, as.character), check.names = FALSE)
  rownames(data) <- NULL

  sheets <- list()
  if (file.exists(filepath)) {
    for (name in readxl::excel_sheets(filepath)) {
      sheets[[name]] <- as.data.frame(readxl::read_excel(filepath, sheet = name, col_types = "text"),
                                      check.names = FALSE)
    }
  }

  if (sheet_name %in% names(sheets)) {
    if (!isTRUE(all.equal(data, sheets[[sheet_name]], check.attributes = FALSE)) ||
        !identical(colnames(data), colnames(sheets[[sheet_name]]))) {
      stop(structure(
        class = c("overwriteError", "error", "condition"),
        list(message = paste0("Current output would overwrite alternative data in tab \"", sheet_name, "\" of ", filepath),
             call = sys.call(-1),
             filepath = filepath)
      ))
    }
    if (verbose) cat("Current output matches tab \"", sheet_name, "\" of existing file.\n", sep = "")
  } else {
    sheets[[sheet_name]] <- data
    dir.create(dirname(filepath), recursive = TRUE, showWarnings = FALSE)
    writexl::write_xlsx(sheets[sort(names(sheets))], filepath)
    if (verbose) cat("Saved new tab \"", sheet_name, "\".\n", sep = "")
  }
  if (verbose) cat(filepath, "\n")
  invisible(NULL)
}
