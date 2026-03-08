# This package streamlines the data pull from GDC
if (!require('remotes', quietly=TRUE)) {
  install.packages('remotes')
}

########################################################################

dash.sep <- strrep('-', 200)
dash.prefix <- paste0(strrep('-', 40), '\t')

writeWithPref <- function(x) {
  writeLines(paste0(dash.prefix, x))
}

TCGAbiolinks_version <- read.table('TCGAbiolinks_version.txt')$V1

########################################################################

writeLines(dash.sep)
writeLines(dash.sep)
if (require('TCGAbiolinks', quietly=TRUE)) {
  if (as.character(packageVersion('TCGAbiolinks')) != TCGAbiolinks_version) {
    writeWithPref('Wrong version of TCGAbiolinks installed, need forked version...')
    writeWithPref('Removing current installation of TCGAbiolinks...')
    remove.packages('TCGAbiolinks')
    unloadNamespace('TCGAbiolinks')
    install_forked = TRUE
  } else {
    writeWithPref('Correct version of TCGAbiolinks already installed!')
    install_forked = FALSE
  }
} else {
  writeWithPref('TCGAbiolinks not installed yet...')
  install_forked = TRUE
}
writeLines(dash.sep)
writeLines(dash.sep)

########################################################################

if (install_forked) {
  writeWithPref('Installing forked version of TCGAbiolinks from danmonyak/TCGAbiolinks')
  writeWithPref('Refrain from updating other packages (skip udates) unless necessary due to errors...')
  writeLines(dash.sep)
  writeLines(dash.sep)
  remotes::install_github("danmonyak/TCGAbiolinks", upgrade="never")
  writeLines(dash.sep)
  writeLines(dash.sep)
  writeWithPref('Installed!')
  writeLines(dash.sep)
  writeLines(dash.sep)
}