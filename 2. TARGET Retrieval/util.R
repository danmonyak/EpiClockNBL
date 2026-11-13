# Shared helper functions for the Rmd files in this directory

keepN_fields <- function(x, n) {
  fields <- strsplit(x, "-")[[1]]
  paste(fields[1:n], collapse = "-")
}
getParticipantIDs <- function(fullTumorIDs) {
  sapply(fullTumorIDs, function (x) keepN_fields(x, 3), USE.NAMES=F)
}
getTumorIDs <- function(fullTumorIDs) {
  sapply(fullTumorIDs, function (x) keepN_fields(x, 4), USE.NAMES=F)
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
