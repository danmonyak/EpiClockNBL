# Runs GSEA (the GSEA desktop app's command-line version) on the .gct and .cls files
# made by GSEA_Data_Preprocessing.Rmd, with the same settings as the runs made in the app.
#
# Run from inside "7. GSEA":
#   Rscript Run_GSEA.R <analysis>
# where <analysis> is "original", "split1", "split2" or "adjusted".
#
# Requires GSEA 4.4.0 (https://www.gsea-msigdb.org/gsea/downloads.jsp), with its install
# folder set as "gsea_dir" in config.json. GSEA downloads the gene sets and chip file from
# the Broad Institute, so this needs an internet connection.
# The results are saved in a new folder named <analysis>.Gsea.<timestamp> in
# official_indir/TARGET/GSEA, which GSEA_Figure.Rmd then finds by itself.

library(jsonlite)

analysis <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(!is.na(analysis), analysis %in% c("original", "split1", "split2", "adjusted"))

consts <- read_json(file.path(repo_dir, "config.json"))
if (is.null(consts[["gsea_dir"]])) stop("Set gsea_dir in config.json to the GSEA install folder")
gsea_dir <- normalizePath(consts[["gsea_dir"]], mustWork = TRUE)
gsea_data_dir <- normalizePath(file.path(consts[["official_indir"]], "TARGET", "GSEA"), mustWork = TRUE)

# Input files from GSEA_Data_Preprocessing.Rmd (same names it uses)
suffix <- if (analysis == "original") "" else paste0(".", analysis)
phenotype <- paste0("neuroblastoma_gmm_phi", suffix)
gct_name <- paste0("neuroblastoma_deseq_scale_factor", suffix, ".gct")
cls_name <- paste0(phenotype, ".cls")
stopifnot(file.exists(file.path(gsea_data_dir, c(gct_name, cls_name))))

# GSEA reads its file paths from the command line, which on Windows can't hold every
# character (official_indir may contain non-ASCII characters). So GSEA runs in a temporary
# folder, and the results folder is copied to official_indir afterwards.
work_dir <- normalizePath(tempfile("gsea_"), mustWork = FALSE)
dir.create(work_dir)
stopifnot(file.copy(file.path(gsea_data_dir, c(gct_name, cls_name)), work_dir))

windows <- .Platform$OS.type == "windows"
pathArg <- function(path) shQuote(path, type = if (windows) "cmd" else "sh")

# GSEA settings, as recorded in the .rpt files of the runs made in the GSEA app
gsea_args <- c(
  "-res", pathArg(file.path(work_dir, gct_name)),
  "-cls", pathArg(paste0(file.path(work_dir, cls_name), "#", phenotype)),
  "-gmx", "ftp.broadinstitute.org://pub/gsea/msigdb/human/gene_sets/h.all.v2026.1.Hs.symbols.gmt",
  "-chip", "ftp.broadinstitute.org://pub/gsea/msigdb/human/annotations/Human_Ensembl_Gene_ID_MSigDB.v2026.1.Hs.chip",
  "-collapse", "Collapse",
  "-mode", "Max_probe",
  "-norm", "meandiv",
  "-nperm", "1000",
  "-permute", "gene_set",
  "-rnd_seed", "0",
  "-rnd_type", "no_balance",
  "-scoring_scheme", "weighted",
  "-metric", "Spearman",
  "-sort", "real",
  "-order", "descending",
  "-include_only_symbols", "true",
  "-make_sets", "true",
  "-median", "false",
  "-num", "100",
  "-plot_top_x", "20",
  "-save_rnd_lists", "false",
  "-set_max", "500",
  "-set_min", "15",
  "-create_gcts", "false",
  "-create_svgs", "false",
  "-zip_report", "false",
  "-rpt_label", analysis,
  "-out", pathArg(work_dir)
)

cat("Running GSEA for", analysis, "...\n")
old_wd <- setwd(gsea_dir)
if (windows) {
  # gsea-cli.bat starts Java in the background and returns right away, so call Java
  # directly with the same arguments, which waits for GSEA to finish
  status <- system2(
    utils::shortPathName(file.path(gsea_dir, "jdk", "bin", "java.exe")),
    c("--module-path=modules", "-Xmx4g", "-Djava.awt.headless=true",
      "-Djava.util.logging.config.file=logging.properties", "@gsea.args",
      "--module=org.gsea_msigdb.gsea/xapps.gsea.CLI", "GSEA", gsea_args)
  )
} else {
  status <- system2(file.path(gsea_dir, "gsea-cli.sh"), c("GSEA", gsea_args))
}
setwd(old_wd)

result_dir <- list.files(work_dir, pattern = paste0("^", analysis, "\\.Gsea\\.[0-9]+$"))
if (status != 0 || length(result_dir) != 1) {
  stop("GSEA did not finish (exit status ", status, "). Its output is in ", work_dir)
}

# Copy the results folder to official_indir
stopifnot(!dir.exists(file.path(gsea_data_dir, result_dir)))
stopifnot(file.copy(file.path(work_dir, result_dir), gsea_data_dir, recursive = TRUE))
unlink(work_dir, recursive = TRUE)
cat("GSEA results saved in", file.path(gsea_data_dir, result_dir), "\n")
