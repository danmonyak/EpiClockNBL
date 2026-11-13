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
