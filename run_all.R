# Run from the release directory: Rscript run_all.R
# Optional: Rscript run_all.R --stage=01   or   Rscript run_all.R --check-only
run_release <- function(root = getwd(), stages = NULL, check_only = FALSE) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  required <- c("readr", "readxl", "dplyr", "tidyr", "purrr", "stringr",
                "tibble", "ggplot2", "scales", "patchwork", "mixOmics")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
  files <- list.files(file.path(root, "scripts"), pattern = "^[0-9]{2}_.*\\.R$", full.names = TRUE)
  if (!is.null(stages)) files <- files[substr(basename(files), 1, 2) %in% stages]
  if (!length(files)) stop("No matching analysis stages.")
  invisible(lapply(files, parse, encoding = "UTF-8"))
  data_file <- function(name) {
    path <- file.path(root, "data", name)
    if (!file.exists(path)) stop("Missing input file: ", path)
    path
  }
  inputs <- c("ms_formula_assignments.xlsx", "MataData.csv", "total_featuretable.csv",
              "WandB_featuretable.csv", "picrust2.group.pathways.relative.xls",
              "picrust2.sample.pathways.relative.xls")
  invisible(lapply(inputs, data_file))
  if (check_only) return(invisible(TRUE))
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  for (script in files) {
    stage <- tools::file_path_sans_ext(basename(script))
    output_dir <- file.path(root, "outputs", stage)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    env <- new.env(parent = globalenv())
    env$data_file <- data_file
    env$output_dir <- env$out_dir <- output_dir
    env$safe_out <- function(filename) file.path(output_dir, filename)
    env$otu_file <- data_file("total_featuretable.csv")
    env$meta_file <- data_file("MataData.csv")
    env$file <- env$ms_file <- data_file("ms_formula_assignments.xlsx")
    setwd(output_dir)
    set.seed(20260921) # Reproducible jitter and any stochastic plotting steps.
    message("Running ", stage)
    grDevices::cairo_pdf(file.path(output_dir, "diagnostic_plots.pdf"), width = 12, height = 9)
    device <- grDevices::dev.cur()
    tryCatch({
      source(script, local = env, encoding = "UTF-8", keep.source = TRUE)
      if (exists("diablo_pool_res", envir = env, inherits = FALSE)) {
        saveRDS(env$diablo_pool_res, file.path(output_dir, "proxy_diablo_model.rds"))
      }
      if (exists("diablo_res", envir = env, inherits = FALSE)) {
        saveRDS(env$diablo_res, file.path(output_dir, "formula_diablo_model.rds"))
      }
      figure_map <- read.csv(file.path(root, "figure_map.csv"), stringsAsFactors = FALSE)
      figure_map <- figure_map[figure_map$stage == stage, ]
      figure_dir <- file.path(root, "figures")
      dir.create(figure_dir, showWarnings = FALSE)
      for (i in seq_len(nrow(figure_map))) {
        spec <- figure_map[i, ]
        obj <- get(spec$object, envir = env, inherits = FALSE)
        base <- file.path(figure_dir, spec$figure)
        if (spec$figure == "FigS6") {
          draw_scores <- function() mixOmics::plotIndiv(obj, group = env$Y_diablo,
            ind.names = TRUE, legend = TRUE, ellipse = TRUE, title = "DIABLO: MS + 16S")
          grDevices::cairo_pdf(paste0(base, ".pdf"), width = spec$width_in, height = spec$height_in)
          draw_scores()
          grDevices::dev.off()
          grDevices::png(paste0(base, ".png"), width = spec$width_in, height = spec$height_in,
                        units = "in", res = 300)
          draw_scores()
          grDevices::dev.off()
        } else {
          ggplot2::ggsave(paste0(base, ".pdf"), obj, device = grDevices::cairo_pdf,
                         width = spec$width_in, height = spec$height_in, bg = "white")
          ggplot2::ggsave(paste0(base, ".png"), obj, width = spec$width_in,
                         height = spec$height_in, dpi = 300, bg = "white")
        }
      }
      saveRDS(as.list(env, all.names = TRUE)[vapply(as.list(env, all.names = TRUE),
        function(x) is.data.frame(x) || is.matrix(x), logical(1))],
        file.path(output_dir, "analysis_tables.rds"))
    }, finally = {if (device %in% grDevices::dev.list()) grDevices::dev.off(device)})
    writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
  }
  message("Completed. Outputs: ", file.path(root, "outputs"))
}

args <- commandArgs(trailingOnly = TRUE)
stage_arg <- grep("^--stage=", args, value = TRUE)
stage_ids <- if (length(stage_arg)) strsplit(sub("^--stage=", "", stage_arg[1]), ",")[[1]] else NULL
run_release(stages = stage_ids, check_only = "--check-only" %in% args)
