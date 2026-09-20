# Optional upstream preparation from vendor formula-assignment CSVs, not raw spectra.
# Run from the release root: Rscript optional/prepare_ms_assignments.R
# The verified manuscript workbook is NEVER overwritten by this script.
required <- c("readr", "dplyr", "writexl")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))
library(dplyr)
input_dir <- file.path(getwd(), "data", "vendor_assignments")
output_dir <- file.path(getwd(), "outputs", "00_optional_ms_preparation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
files <- file.path(input_dir, paste0("sample ", 1:10, ".csv"))
stopifnot(all(file.exists(files)))

process_dom_file <- function(file_path) {
  sample_name <- tools::file_path_sans_ext(basename(file_path))
  message("Processing: ", sample_name)
  df <- readr::read_csv(file_path, show_col_types = FALSE) %>%
    dplyr::select(where(~ !all(is.na(.)))) %>%
    dplyr::select(-matches("^Unnamed")) %>%
    dplyr::select(-matches("^\\.\\.\\."))
  for (element in c("N", "S", "P")) {
    if (!element %in% names(df)) df[[element]] <- 0
  }
  # Preserve the original H + 1 convention and formula-selection priority.
  filtered <- df %>%
    mutate(
      N = ifelse(is.na(N), 0, N), S = ifelse(is.na(S), 0, S),
      P = ifelse(is.na(P), 0, P), H_plus_1 = H + 1,
      HC_ratio = H_plus_1 / C, OC_ratio = O / C,
      abs_error = abs(`err ppm`), DBE = 1 + (2 * C - H_plus_1 + N) / 2,
      AI_denominator = C - 0.5 * O - S - N - P,
      AI_mod = ifelse(AI_denominator > 0,
        (1 + C - 0.5 * O - S - 0.5 * (H_plus_1 + N + P)) / AI_denominator,
        NA_real_)
    ) %>%
    filter(`Observed m/z` >= 100, `Observed m/z` <= 800,
           HC_ratio >= 0.2, HC_ratio <= 2.3, OC_ratio >= 0, OC_ratio <= 1.2) %>%
    mutate(sum_NS = N + S) %>%
    arrange(`Observed m/z`, sum_NS, abs_error) %>%
    distinct(`Observed m/z`, .keep_all = TRUE) %>%
    mutate(Sample = sample_name, .before = 1)
  writexl::write_xlsx(filtered, file.path(output_dir, paste0(sample_name, "_filtered.xlsx")))
  filtered
}
result_list <- lapply(files, process_dom_file)
names(result_list) <- paste("sample", 1:10)
writexl::write_xlsx(bind_rows(result_list), file.path(output_dir, "all_samples_filtered.xlsx"))
writexl::write_xlsx(result_list, file.path(output_dir, "regenerated_ms_assignments.xlsx"))
message("Preparation complete. Compare regenerated data with the manuscript workbook before substitution.")
