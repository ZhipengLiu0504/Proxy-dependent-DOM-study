# Run through run_all.R; see README.md for inputs and figure mapping.
# Figure S6 uses the formula-level model, not the proxy-level model in Figure 3.
# Preserve the original top-variance prefilters and five selected features per block.
# Input signals are relative-intensity/log1p transformed; genus preprocessing differs
# from the proxy workflow and is deliberately retained for figure reproduction.
# Ellipses and score separation are descriptive, with no held-out validation here.
library(readr)
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(ggplot2)
library(mixOmics)
# 1.2 Helper functions for robust column reading
# ============================================================

num_col <- function(df, candidates, default = NA_real_) {
  candidates <- candidates[candidates %in% names(df)]
  if (length(candidates) == 0) {
    return(rep(default, nrow(df)))
  } else {
    return(suppressWarnings(as.numeric(df[[candidates[1]]])))
  }
}

char_col <- function(df, candidates, default = NA_character_) {
  candidates <- candidates[candidates %in% names(df)]
  if (length(candidates) == 0) {
    return(rep(default, nrow(df)))
  } else {
    return(as.character(df[[candidates[1]]]))
  }
}

first_non_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA_real_)
  } else {
    return(x[1])
  }
}

# ============================================================
# 2. Read and clean metadata
# ============================================================

metadata <- readr::read_csv(meta_file, show_col_types = FALSE)

# Remove spaces from column names and ID values
colnames(metadata) <- trimws(colnames(metadata))
metadata$ID <- trimws(metadata$ID)

matched_meta <- metadata %>%
  dplyr::mutate(
    `sample ID in MS` = na_if(`sample ID in MS`, "NA"),
    MS_sample = `sample ID in MS`
  ) %>%
  dplyr::filter(!is.na(MS_sample)) %>%
  dplyr::mutate(
    MS_sample = trimws(MS_sample),
    MS_order = as.numeric(stringr::str_extract(MS_sample, "\\d+")),
    environment = tolower(trimws(as.character(environment))),
    environment = dplyr::recode(
      environment,
      "water" = "Bulk water",
      "biofilm" = "Biofilm"
    ),
    environment = factor(environment, levels = c("Bulk water", "Biofilm"))
  ) %>%
  dplyr::arrange(MS_order)

print(matched_meta %>% dplyr::select(ID, environment, reservoir, site, MS_sample, MS_order))
print(table(matched_meta$environment))

# ============================================================
# 3. Read and process 16S genus table
# ============================================================

otu_raw <- readr::read_csv(otu_file, show_col_types = FALSE)

# Clean column names and sample IDs
colnames(otu_raw) <- trimws(colnames(otu_raw))
matched_meta$ID <- trimws(matched_meta$ID)

# Remove unwanted taxonomy
otu_clean <- otu_raw %>%
  dplyr::filter(!is.na(Taxonomy)) %>%
  dplyr::mutate(Taxonomy = trimws(Taxonomy)) %>%
  dplyr::filter(Taxonomy != "Others") %>%
  dplyr::filter(
    !stringr::str_detect(
      Taxonomy,
      stringr::regex("chloroplast|mitochondria", ignore_case = TRUE)
    )
  )

# Check whether the 10 matched 16S IDs exist in the 16S feature table
missing_cols <- setdiff(matched_meta$ID, colnames(otu_clean))

if (length(missing_cols) > 0) {
  stop(
    paste0(
      "These 16S sample IDs in metadata are missing from the 16S table: ",
      paste(missing_cols, collapse = ", ")
    )
  )
}

# Keep only the 10 matched 16S samples
otu_10 <- otu_clean %>%
  dplyr::select(
    Taxonomy,
    tidyselect::all_of(matched_meta$ID)
  )

# Convert genus x sample into sample x genus
# Also merge duplicated genus names by summing abundance
otu_mat <- otu_10 %>%
  tidyr::pivot_longer(
    cols = -Taxonomy,
    names_to = "ID",
    values_to = "Abundance"
  ) %>%
  dplyr::mutate(
    Abundance = as.numeric(Abundance)
  ) %>%
  dplyr::left_join(
    matched_meta %>%
      dplyr::select(ID, MS_sample, MS_order, environment),
    by = "ID"
  ) %>%
  dplyr::filter(!is.na(MS_sample)) %>%
  dplyr::group_by(MS_sample, MS_order, Taxonomy) %>%
  dplyr::summarise(
    Abundance = sum(Abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(MS_order) %>%
  dplyr::select(MS_sample, Taxonomy, Abundance) %>%
  tidyr::pivot_wider(
    names_from = Taxonomy,
    values_from = Abundance,
    values_fill = 0
  ) %>%
  tibble::column_to_rownames("MS_sample") %>%
  as.matrix()

# Remove possible NA values
otu_mat[is.na(otu_mat)] <- 0

cat("\n16S matrix dimension:\n")
print(dim(otu_mat))
print(rownames(otu_mat))
print(otu_mat[1:min(5, nrow(otu_mat)), 1:min(5, ncol(otu_mat))])

# ============================================================
# 4. Read and process MS data
#    AImod is calculated here.
#    Condensed aromatics are annotated as AI_mod >= 0.67.
# ============================================================

# Check sheet names
all_sheets <- readxl::excel_sheets(ms_file)
print(all_sheets)

# We only need sample 1 to sample 10
wanted_sheets_lower <- paste("sample", 1:10)
sheet_match <- match(tolower(wanted_sheets_lower), tolower(all_sheets))

if (any(is.na(sheet_match))) {
  missing_sheets <- wanted_sheets_lower[is.na(sheet_match)]
  stop(
    paste0(
      "These MS sheets are missing from Excel file: ",
      paste(missing_sheets, collapse = ", ")
    )
  )
}

ms_sheets <- all_sheets[sheet_match]
names(ms_sheets) <- wanted_sheets_lower

ms_long <- lapply(seq_along(ms_sheets), function(i) {

  sh <- ms_sheets[i]
  ms_sample_name <- names(ms_sheets)[i]

  df <- readxl::read_excel(ms_file, sheet = sh)
  colnames(df) <- trimws(colnames(df))

  df2 <- tibble::tibble(
    MS_sample = ms_sample_name,
    Formula = stringr::str_remove_all(
      char_col(df, c("sum formula", "Neutral_Formula", "Formula")),
      "\\s+"
    ),
    Intensity = num_col(df, c("Observed Intens", "Observed_Intens", "Intensity")),

    C = num_col(df, c("C")),
    H_raw = num_col(df, c("H_neutral", "H")),
    H_plus_1_existing = num_col(df, c("H_plus_1", "H.plus.1"), default = NA_real_),
    O = num_col(df, c("O")),
    N = num_col(df, c("N"), default = 0),
    S = num_col(df, c("S"), default = 0),
    P = num_col(df, c("P"), default = 0),

    H_C_existing = num_col(df, c("H/C", "HC_ratio", "H_C"), default = NA_real_),
    O_C_existing = num_col(df, c("O/C", "OC_ratio", "O_C"), default = NA_real_),
    DBE_existing = num_col(df, c("DBE"), default = NA_real_),
    AI_mod_existing = num_col(df, c("AI_mod", "AImod", "AI.mod"), default = NA_real_)
  ) %>%
    dplyr::mutate(
      N = tidyr::replace_na(N, 0),
      S = tidyr::replace_na(S, 0),
      P = tidyr::replace_na(P, 0),

      # Use H_plus_1 if it exists. Otherwise use H + 1.
      # This keeps the AImod calculation consistent with your filtered MS table.
      H_for_AI = dplyr::case_when(
        is.finite(H_plus_1_existing) ~ H_plus_1_existing,
        is.finite(H_raw) ~ H_raw + 1,
        TRUE ~ NA_real_
      ),

      H_C = dplyr::case_when(
        is.finite(H_C_existing) ~ H_C_existing,
        is.finite(H_for_AI) & is.finite(C) & C > 0 ~ H_for_AI / C,
        TRUE ~ NA_real_
      ),
      O_C = dplyr::case_when(
        is.finite(O_C_existing) ~ O_C_existing,
        is.finite(O) & is.finite(C) & C > 0 ~ O / C,
        TRUE ~ NA_real_
      ),
      DBE = dplyr::case_when(
        is.finite(DBE_existing) ~ DBE_existing,
        is.finite(C) & is.finite(H_for_AI) & is.finite(N) ~
          1 + (2 * C - H_for_AI + N) / 2,
        TRUE ~ NA_real_
      ),

      # Modified aromaticity index
      AI_denominator = C - 0.5 * O - S - N - P,
      AI_mod = dplyr::case_when(
        is.finite(AI_mod_existing) ~ AI_mod_existing,
        is.finite(AI_denominator) & AI_denominator > 0 ~
          (1 + C - 0.5 * O - S - 0.5 * (H_for_AI + N + P)) / AI_denominator,
        TRUE ~ NA_real_
      ),

      # New requirement:
      # condensed aromatics are annotated by AImod >= 0.67
      Formula_class_AImod = dplyr::case_when(
        is.finite(AI_mod) & AI_mod >= 0.67 ~ "Condensed_Aromatics",
        TRUE ~ "Other_or_unclassified"
      )
    ) %>%
    dplyr::filter(
      !is.na(Formula),
      Formula != "",
      !is.na(Intensity),
      is.finite(Intensity),
      Intensity > 0
    ) %>%
    dplyr::select(
      MS_sample,
      Formula,
      Intensity,
      C,
      H_for_AI,
      O,
      N,
      S,
      P,
      H_C,
      O_C,
      DBE,
      AI_denominator,
      AI_mod,
      Formula_class_AImod
    )

  return(df2)

}) %>%
  dplyr::bind_rows()

cat("\nFormula class counts based on AImod:\n")
print(table(ms_long$Formula_class_AImod, useNA = "ifany"))

# Save formula-level annotation so selected MS features can be interpreted later
write_csv(
  ms_long,
  safe_out("MS_formula_AImod_annotation_by_sample.csv")
)

# Merge duplicated formulas within the same sample
ms_long_sum <- ms_long %>%
  dplyr::group_by(MS_sample, Formula) %>%
  dplyr::summarise(
    Intensity = sum(Intensity, na.rm = TRUE),
    C = first_non_na(C),
    H_for_AI = first_non_na(H_for_AI),
    O = first_non_na(O),
    N = first_non_na(N),
    S = first_non_na(S),
    P = first_non_na(P),
    H_C = first_non_na(H_C),
    O_C = first_non_na(O_C),
    DBE = first_non_na(DBE),
    AI_denominator = first_non_na(AI_denominator),
    AI_mod = ifelse(all(is.na(AI_mod)), NA_real_, mean(AI_mod, na.rm = TRUE)),
    Formula_class_AImod = ifelse(
      any(Formula_class_AImod == "Condensed_Aromatics", na.rm = TRUE),
      "Condensed_Aromatics",
      "Other_or_unclassified"
    ),
    .groups = "drop"
  )

write_csv(
  ms_long_sum,
  safe_out("MS_formula_AImod_annotation_summed_by_sample.csv")
)

# Formula-level annotation across all samples
ms_formula_annotation <- ms_long_sum %>%
  dplyr::group_by(Formula) %>%
  dplyr::summarise(
    C = first_non_na(C),
    H_for_AI = first_non_na(H_for_AI),
    O = first_non_na(O),
    N = first_non_na(N),
    S = first_non_na(S),
    P = first_non_na(P),
    H_C = first_non_na(H_C),
    O_C = first_non_na(O_C),
    DBE = first_non_na(DBE),
    AI_mod = ifelse(all(is.na(AI_mod)), NA_real_, mean(AI_mod, na.rm = TRUE)),
    Formula_class_AImod = ifelse(
      any(Formula_class_AImod == "Condensed_Aromatics", na.rm = TRUE),
      "Condensed_Aromatics",
      "Other_or_unclassified"
    ),
    .groups = "drop"
  )

write_csv(
  ms_formula_annotation,
  safe_out("MS_formula_AImod_annotation_unique_formula.csv")
)

# Convert intensity into relative intensity within each sample
ms_long_rel <- ms_long_sum %>%
  dplyr::group_by(MS_sample) %>%
  dplyr::mutate(
    RelIntensity = Intensity / sum(Intensity, na.rm = TRUE)
  ) %>%
  dplyr::ungroup()

# Convert long table into sample x formula matrix
ms_mat <- ms_long_rel %>%
  dplyr::select(MS_sample, Formula, RelIntensity) %>%
  tidyr::pivot_wider(
    names_from = Formula,
    values_from = RelIntensity,
    values_fill = 0
  ) %>%
  dplyr::mutate(
    MS_order = as.numeric(stringr::str_extract(MS_sample, "\\d+"))
  ) %>%
  dplyr::arrange(MS_order) %>%
  dplyr::select(-MS_order) %>%
  tibble::column_to_rownames("MS_sample") %>%
  as.matrix()

# Remove possible NA values
ms_mat[is.na(ms_mat)] <- 0

cat("\nMS matrix dimension:\n")
print(dim(ms_mat))
print(rownames(ms_mat))
print(ms_mat[1:min(5, nrow(ms_mat)), 1:min(5, ncol(ms_mat))])

# ============================================================
# 5. Align MS and 16S matrices
# ============================================================

cat("\nCheck sample names before alignment:\n")
print(rownames(ms_mat))
print(rownames(otu_mat))

# Reorder 16S matrix according to MS matrix
otu_mat <- otu_mat[rownames(ms_mat), , drop = FALSE]

cat("\nAre MS and 16S samples aligned?\n")
print(all(rownames(ms_mat) == rownames(otu_mat)))

if (!all(rownames(ms_mat) == rownames(otu_mat))) {
  stop("MS and 16S rownames are still not matched. Please check metadata.")
}

# ============================================================
# 6. Prepare Y group information
# ============================================================

Y <- matched_meta %>%
  dplyr::arrange(MS_order) %>%
  dplyr::pull(environment)

names(Y) <- matched_meta %>%
  dplyr::arrange(MS_order) %>%
  dplyr::pull(MS_sample)

# Reorder Y according to MS matrix
Y <- Y[rownames(ms_mat)]

cat("\nGroup information:\n")
print(Y)
print(table(Y))

if (!all(rownames(ms_mat) == names(Y))) {
  stop("Y names are not matched with MS matrix rownames.")
}

# ============================================================
# 7. Log transformation
# ============================================================

ms_mat_log  <- log1p(ms_mat)
otu_mat_log <- log1p(otu_mat)

# ============================================================
# 8. Feature filtering
# Because you only have 10 samples, do not use all features.
# Keep high-variance features only.
# ============================================================

# ---------- MS filtering ----------
top_n_ms <- 200

ms_var <- apply(ms_mat_log, 2, var, na.rm = TRUE)

# Remove zero-variance or NA-variance features
ms_var <- ms_var[!is.na(ms_var) & ms_var > 0]

ms_keep <- names(sort(ms_var, decreasing = TRUE))[1:min(top_n_ms, length(ms_var))]

ms_mat_filt <- ms_mat_log[, ms_keep, drop = FALSE]

# ---------- 16S filtering ----------
top_n_otu <- 30

otu_var <- apply(otu_mat_log, 2, var, na.rm = TRUE)

# Remove zero-variance or NA-variance features
otu_var <- otu_var[!is.na(otu_var) & otu_var > 0]

otu_keep <- names(sort(otu_var, decreasing = TRUE))[1:min(top_n_otu, length(otu_var))]

otu_mat_filt <- otu_mat_log[, otu_keep, drop = FALSE]

cat("\nFiltered MS matrix dimension:\n")
print(dim(ms_mat_filt))

cat("\nFiltered 16S matrix dimension:\n")
print(dim(otu_mat_filt))

# ============================================================
# 9. Prepare DIABLO input
# ============================================================

data_diablo <- list(
  MS = ms_mat_filt,
  Microbe = otu_mat_filt
)

Y_diablo <- Y

cat("\nFinal checks before DIABLO:\n")
print(all(rownames(data_diablo$MS) == rownames(data_diablo$Microbe)))
print(all(rownames(data_diablo$MS) == names(Y_diablo)))
print(table(Y_diablo))
print(dim(data_diablo$MS))
print(dim(data_diablo$Microbe))

stopifnot(all(rownames(data_diablo$MS) == rownames(data_diablo$Microbe)))
stopifnot(all(rownames(data_diablo$MS) == names(Y_diablo)))
stopifnot(length(Y_diablo) == nrow(data_diablo$MS))

# ============================================================
# 10. Save processed files
# ============================================================

write_csv(
  as.data.frame(ms_mat_filt) %>%
    tibble::rownames_to_column("MS_sample"),
  safe_out("MS_matrix_for_DIABLO.csv")
)

write_csv(
  as.data.frame(otu_mat_filt) %>%
    tibble::rownames_to_column("MS_sample"),
  safe_out("16S_matrix_for_DIABLO.csv")
)

write_csv(
  matched_meta,
  safe_out("matched_metadata_for_DIABLO.csv")
)

# ============================================================
# 11. Run first DIABLO model
# ============================================================

# Design matrix controls how strongly MS and Microbe blocks are connected.
# 0.1 is conservative and suitable for small sample size.
design <- matrix(
  c(0, 0.1,
    0.1, 0),
  nrow = 2,
  byrow = TRUE
)

rownames(design) <- colnames(design) <- c("MS", "Microbe")
diablo_res <- block.splsda(
  X = data_diablo,
  Y = Y_diablo,
  ncomp = 2,
  keepX = list(
    MS = c(5, 5),
    Microbe = c(5, 5)
  ),
  design = design
)

print(diablo_res)

# ============================================================
# 12. Basic plots
# ============================================================

# Sample separation plot
plotIndiv(
  diablo_res,
  group = Y_diablo,
  ind.names = TRUE,
  legend = TRUE,
  ellipse = TRUE,
  title = "DIABLO: MS + 16S"
)

