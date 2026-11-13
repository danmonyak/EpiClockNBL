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
