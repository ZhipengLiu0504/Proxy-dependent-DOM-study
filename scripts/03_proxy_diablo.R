# Run through run_all.R; see README.md for inputs and figure mapping.
# Figures 3-5 and S7: preserve the proxy-level model and published plotting choices.
# Fig. 4a/b: mean +/- SD, five samples per habitat, unpaired Welch tests.
# Fig. 4c/d: three matched pairs in reservoir 1, paired t tests.
# Stars and Fig. 5 filtering use raw p values; exported BH values are additional context.
# DIABLO loadings are model coefficients, NOT fold changes or significance tests.
# Loading signs can reverse across fits; the original sign-to-habitat colors are retained.
# DIABLO separation here is in-sample, not an independently validated prediction result.
library(readr)
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(ggplot2)
library(scales)
library(mixOmics)
theme_ms_plot <- function(base_size = 16, base_family = "Arial") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(size = base_size + 3, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = base_size - 1, face = "bold", hjust = 0.5),
      axis.title = element_text(size = base_size + 1, face = "bold", colour = "black"),
      axis.text = element_text(size = base_size - 2, face = "bold", colour = "black"),
      strip.text = element_text(size = base_size - 1, face = "bold", colour = "black"),
      legend.title = element_text(size = base_size - 1, face = "bold"),
      legend.text = element_text(size = base_size - 2, face = "bold"),
      plot.caption = element_text(size = base_size - 3, face = "bold")
    )
}

# ============================================================
# 2. Helper functions
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

# ============================================================
# 3. Read and clean metadata
# ============================================================

metadata <- readr::read_csv(meta_file, show_col_types = FALSE)

colnames(metadata) <- trimws(colnames(metadata))
metadata$ID <- trimws(metadata$ID)

matched_meta <- metadata %>%
  dplyr::mutate(
    `sample ID in MS` = na_if(`sample ID in MS`, "NA"),
    MS_sample_original = `sample ID in MS`
  ) %>%
  dplyr::filter(!is.na(MS_sample_original)) %>%
  dplyr::mutate(
    MS_order = as.numeric(stringr::str_extract(MS_sample_original, "\\d+")),
    MS_sample = paste("sample", MS_order),
    environment = tolower(trimws(as.character(environment))),
    environment = dplyr::recode(
      environment,
      "water" = "Bulk water",
      "biofilm" = "Biofilm"
    ),
    environment = factor(environment, levels = c("Bulk water", "Biofilm"))
  ) %>%
  dplyr::arrange(MS_order)

cat("\nMatched metadata:\n")
print(matched_meta %>% dplyr::select(ID, environment, reservoir, site, MS_sample, MS_order))
print(table(matched_meta$environment))

# ============================================================
# 4. Read and process 16S genus table
# ============================================================

otu_raw <- readr::read_csv(otu_file, show_col_types = FALSE)

colnames(otu_raw) <- trimws(colnames(otu_raw))
matched_meta$ID <- trimws(matched_meta$ID)

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

missing_cols <- setdiff(matched_meta$ID, colnames(otu_clean))

if (length(missing_cols) > 0) {
  stop(
    paste0(
      "These 16S sample IDs in metadata are missing from the 16S table: ",
      paste(missing_cols, collapse = ", ")
    )
  )
}

otu_10 <- otu_clean %>%
  dplyr::select(
    Taxonomy,
    tidyselect::all_of(matched_meta$ID)
  )

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

otu_mat[is.na(otu_mat)] <- 0

cat("\n16S matrix dimension:\n")
print(dim(otu_mat))
print(rownames(otu_mat))

# ============================================================
# 5. Read MS data and classify DOM formulas into proxy pools
#    Condensed_Aromatics is now defined by AImod >= 0.67
# ============================================================

all_sheets <- readxl::excel_sheets(ms_file)

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

pool_levels <- c(
  "CRAM",
  "Tannins",
  "Condensed_Aromatics",
  "Lignins",
  "Unsaturated_Hydrocarbons",
  "Carbohydrates",
  "Proteins",
  "Lipids",
  "Other"
)

ms_pool_long <- lapply(seq_along(ms_sheets), function(i) {
  
  sh <- ms_sheets[i]
  ms_sample_name <- names(ms_sheets)[i]
  
  df <- readxl::read_excel(ms_file, sheet = sh)
  colnames(df) <- trimws(colnames(df))
  
  if ("Keep" %in% names(df)) {
    df <- df %>%
      dplyr::filter(
        is.na(Keep) |
          Keep == TRUE |
          Keep == "TRUE" |
          Keep == "True" |
          Keep == 1
      )
  }
  
  df2 <- tibble::tibble(
    MS_sample = ms_sample_name,
    Formula   = stringr::str_remove_all(
      char_col(df, c("sum formula", "Neutral_Formula", "Formula")),
      "\\s+"
    ),
    Intensity = num_col(df, c("Observed Intens", "Observed_Intens", "Intensity")),
    
    H_C = num_col(df, c("H/C", "HC_ratio", "H_C")),
    O_C = num_col(df, c("O/C", "OC_ratio", "O_C")),
    
    C = num_col(df, c("C")),
    H_raw = num_col(df, c("H_neutral", "H")),
    H_plus_1_existing = num_col(df, c("H_plus_1", "H.plus.1"), default = NA_real_),
    O = num_col(df, c("O")),
    N = num_col(df, c("N"), default = 0),
    S = num_col(df, c("S"), default = 0),
    P = num_col(df, c("P"), default = 0),
    DBE = num_col(df, c("DBE")),
    AI_mod_existing = num_col(df, c("AI_mod", "AImod", "AI.mod"), default = NA_real_)
  ) %>%
    dplyr::mutate(
      N = tidyr::replace_na(N, 0),
      S = tidyr::replace_na(S, 0),
      P = tidyr::replace_na(P, 0),
      
      # Use existing H_plus_1 if available.
      # Otherwise use H_raw + 1, consistent with your previous filtering code.
      H_for_AI = dplyr::case_when(
        is.finite(H_plus_1_existing) ~ H_plus_1_existing,
        is.finite(H_raw) ~ H_raw + 1,
        TRUE ~ NA_real_
      ),
      
      # If H/C or O/C were not present, calculate them here.
      H_C = dplyr::case_when(
        is.finite(H_C) ~ H_C,
        is.finite(H_for_AI) & is.finite(C) & C > 0 ~ H_for_AI / C,
        TRUE ~ NA_real_
      ),
      O_C = dplyr::case_when(
        is.finite(O_C) ~ O_C,
        is.finite(O) & is.finite(C) & C > 0 ~ O / C,
        TRUE ~ NA_real_
      ),
      
      # If DBE was not present, calculate it using H_for_AI.
      DBE = dplyr::case_when(
        is.finite(DBE) ~ DBE,
        is.finite(C) & is.finite(H_for_AI) & is.finite(N) ~
          1 + (2 * C - H_for_AI + N) / 2,
        TRUE ~ NA_real_
      ),
      
      # AImod calculation
      AI_denominator = C - 0.5 * O - S - N - P,
      AI_mod = dplyr::case_when(
        is.finite(AI_mod_existing) ~ AI_mod_existing,
        is.finite(AI_denominator) & AI_denominator > 0 ~
          (1 + C - 0.5 * O - S - 0.5 * (H_for_AI + N + P)) / AI_denominator,
        TRUE ~ NA_real_
      )
    ) %>%
    dplyr::filter(
      !is.na(Intensity),
      Intensity > 0,
      is.finite(Intensity),
      is.finite(C),
      is.finite(O),
      is.finite(H_for_AI)
    ) %>%
    dplyr::mutate(
      DBE_C = ifelse(C > 0, DBE / C, NA_real_),
      DBE_H = ifelse(H_for_AI > 0, DBE / H_for_AI, NA_real_),
      DBE_O = ifelse(O > 0, DBE / O, NA_real_),
      N_C   = ifelse(C > 0, N / C, NA_real_),
      P_C   = ifelse(C > 0, P / C, NA_real_)
    ) %>%
    dplyr::mutate(
      pool = dplyr::case_when(
        
        # RDOC proxy: CRAM
        is.finite(DBE_C) & is.finite(DBE_H) & is.finite(DBE_O) &
          DBE_C >= 0.30 & DBE_C <= 0.68 &
          DBE_H >= 0.20 & DBE_H <= 0.95 &
          DBE_O >= 0.77 & DBE_O <= 1.75 ~ "CRAM",
        
        # New standard: AImod-based condensed aromatics
        # AImod >= 0.67 is commonly used for condensed aromatic / black-carbon-like formulas.
        is.finite(AI_mod) & AI_mod >= 0.67 ~ "Condensed_Aromatics",
        
        # Carbohydrate-like
        is.finite(O_C) & is.finite(H_C) &
          O_C >= 0.8 &
          H_C >= 1.65 & H_C < 2.7 &
          N == 0 ~ "Carbohydrates",
        
        # Protein-like
        is.finite(O_C) & is.finite(H_C) & is.finite(N_C) & is.finite(P_C) &
          (
            (
              O_C > 0.12 & O_C <= 0.6 &
                H_C > 0.9 & H_C < 2.5 &
                N_C >= 0.126 & N_C <= 0.7 &
                P_C < 0.17
            ) |
              (
                O_C > 0.6 & O_C <= 1.0 &
                  H_C > 1.2 & H_C < 2.5 &
                  N_C > 0.2 & N_C <= 0.7 &
                  P_C < 0.17
              )
          ) ~ "Proteins",
        
        # Lipid-like
        is.finite(O_C) & is.finite(H_C) & is.finite(N_C) & is.finite(P_C) &
          O_C <= 0.6 &
          H_C >= 1.32 &
          N_C <= 0.125 &
          P_C < 0.35 ~ "Lipids",
        
        # Tannin-like
        is.finite(O_C) & is.finite(H_C) &
          O_C >= 0.67 & O_C <= 0.97 &
          H_C >= 0.53 & H_C <= 1.5 ~ "Tannins",
        
        # Lignin-like
        is.finite(O_C) & is.finite(H_C) &
          O_C >= 0.25 & O_C <= 0.67 &
          H_C >= 0.75 & H_C <= 1.5 ~ "Lignins",
        
        # Unsaturated hydrocarbons
        is.finite(O_C) & is.finite(H_C) &
          O_C >= 0.00 & O_C <= 0.29 &
          H_C >= 1.0 & H_C <= 1.6 ~ "Unsaturated_Hydrocarbons",
        
        TRUE ~ "Other"
      )
    ) %>%
    dplyr::select(
      MS_sample,
      Formula,
      Intensity,
      H_C,
      O_C,
      C,
      H_for_AI,
      O,
      N,
      S,
      P,
      DBE,
      AI_denominator,
      AI_mod,
      DBE_C,
      DBE_H,
      DBE_O,
      N_C,
      P_C,
      pool
    )
  
  return(df2)
  
}) %>%
  dplyr::bind_rows()

cat("\nDOM pool counts after AImod-based classification:\n")
print(table(ms_pool_long$pool))

write_csv(
  ms_pool_long,
  "MS_formula_pool_classification_with_AImod.csv"
)

# ============================================================
# 6. Calculate DOM pool relative abundance per sample
# ============================================================

df_comp <- ms_pool_long %>%
  dplyr::mutate(
    pool = factor(pool, levels = pool_levels)
  ) %>%
  dplyr::group_by(MS_sample, pool) %>%
  dplyr::summarise(
    I = sum(Intensity, na.rm = TRUE),
    Formula_count = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::group_by(MS_sample) %>%
  dplyr::mutate(
    relI = I / sum(I, na.rm = TRUE)
  ) %>%
  dplyr::ungroup()

# 7. Convert DOM pools into sample x proxy matrix
# ============================================================

ms_pool_mat <- df_comp %>%
  dplyr::select(MS_sample, pool, relI) %>%
  tidyr::pivot_wider(
    names_from = pool,
    values_from = relI,
    values_fill = 0
  ) %>%
  dplyr::mutate(
    MS_order = as.numeric(stringr::str_extract(MS_sample, "\\d+"))
  ) %>%
  dplyr::arrange(MS_order) %>%
  dplyr::select(-MS_order) %>%
  tibble::column_to_rownames("MS_sample") %>%
  as.data.frame()

for (p in pool_levels) {
  if (!p %in% colnames(ms_pool_mat)) {
    ms_pool_mat[[p]] <- 0
  }
}

ms_pool_mat <- ms_pool_mat[, pool_levels, drop = FALSE]

ms_pool_mat <- ms_pool_mat %>%
  dplyr::mutate(
    # CRAM itself is the RDOC proxy
    BDOC_PCL = Proteins + Carbohydrates + Lipids,
    RDOC_to_BDOC = ifelse(BDOC_PCL > 0, CRAM / BDOC_PCL, 0)
  ) %>%
  as.matrix()

ms_pool_mat[is.na(ms_pool_mat)] <- 0

cat("\nDOM proxy matrix dimension:\n")
print(dim(ms_pool_mat))
print(rownames(ms_pool_mat))
print(ms_pool_mat)

# ============================================================
# 8. Align DOM proxy matrix and 16S matrix
# ============================================================

otu_mat <- otu_mat[rownames(ms_pool_mat), , drop = FALSE]

cat("\nAre DOM proxy and 16S samples aligned?\n")
print(all(rownames(ms_pool_mat) == rownames(otu_mat)))

if (!all(rownames(ms_pool_mat) == rownames(otu_mat))) {
  stop("DOM proxy matrix and 16S matrix rownames are not matched.")
}

# ============================================================
# 9. Prepare Y group information
# ============================================================

Y_diablo <- matched_meta %>%
  dplyr::arrange(MS_order) %>%
  dplyr::pull(environment)

names(Y_diablo) <- matched_meta %>%
  dplyr::arrange(MS_order) %>%
  dplyr::pull(MS_sample)

Y_diablo <- Y_diablo[rownames(ms_pool_mat)]

cat("\nGroup information:\n")
print(Y_diablo)
print(table(Y_diablo))

if (!all(rownames(ms_pool_mat) == names(Y_diablo))) {
  stop("Y names are not matched with DOM proxy matrix rownames.")
}

# ============================================================
# 10. Preprocess matrices for DIABLO
# ============================================================

# ---------- DOM proxy block ----------
# Remove zero-variance DOM proxies, then scale
dom_var <- apply(ms_pool_mat, 2, var, na.rm = TRUE)
dom_var <- dom_var[!is.na(dom_var) & dom_var > 0]

dom_keep <- names(dom_var)

ms_pool_scaled <- scale(ms_pool_mat[, dom_keep, drop = FALSE])
ms_pool_scaled[is.na(ms_pool_scaled)] <- 0

# ---------- Microbe block ----------
# Normalize genus abundance within each sample, then log-transform
otu_row_sum <- rowSums(otu_mat, na.rm = TRUE)

otu_mat_rel <- sweep(otu_mat, 1, otu_row_sum, FUN = "/")
otu_mat_rel[is.na(otu_mat_rel)] <- 0
otu_mat_rel[is.infinite(otu_mat_rel)] <- 0

otu_mat_log <- log1p(otu_mat_rel)

top_n_otu <- 30

otu_var <- apply(otu_mat_log, 2, var, na.rm = TRUE)
otu_var <- otu_var[!is.na(otu_var) & otu_var > 0]

otu_keep <- names(sort(otu_var, decreasing = TRUE))[1:min(top_n_otu, length(otu_var))]

otu_mat_filt <- otu_mat_log[, otu_keep, drop = FALSE]

# Re-align after filtering
otu_mat_filt <- otu_mat_filt[rownames(ms_pool_scaled), , drop = FALSE]

cat("\nFiltered DOM proxy matrix dimension:\n")
print(dim(ms_pool_scaled))

cat("\nFiltered 16S matrix dimension:\n")
print(dim(otu_mat_filt))

# ============================================================
# 11. Prepare DIABLO input
# ============================================================

data_diablo_pool <- list(
  DOM_proxy = ms_pool_scaled,
  Microbe = otu_mat_filt
)

cat("\nFinal checks before proxy-level DIABLO:\n")
print(all(rownames(data_diablo_pool$DOM_proxy) == rownames(data_diablo_pool$Microbe)))
print(all(rownames(data_diablo_pool$DOM_proxy) == names(Y_diablo)))
print(table(Y_diablo))
print(dim(data_diablo_pool$DOM_proxy))
print(dim(data_diablo_pool$Microbe))

stopifnot(all(rownames(data_diablo_pool$DOM_proxy) == rownames(data_diablo_pool$Microbe)))
stopifnot(all(rownames(data_diablo_pool$DOM_proxy) == names(Y_diablo)))
stopifnot(length(Y_diablo) == nrow(data_diablo_pool$DOM_proxy))

# ============================================================
# 12. Save processed proxy-level matrices
# ============================================================

write_csv(
  as.data.frame(ms_pool_mat) %>%
    tibble::rownames_to_column("MS_sample"),
  "DOM_proxy_matrix_raw_for_DIABLO.csv"
)

write_csv(
  as.data.frame(ms_pool_scaled) %>%
    tibble::rownames_to_column("MS_sample"),
  "DOM_proxy_matrix_scaled_for_DIABLO.csv"
)

write_csv(
  as.data.frame(otu_mat_filt) %>%
    tibble::rownames_to_column("MS_sample"),
  "16S_genus_matrix_for_proxy_DIABLO.csv"
)

write_csv(
  matched_meta,
  "matched_metadata_for_proxy_DIABLO.csv"
)

# ============================================================
# 13. Run proxy-level DIABLO model
# ============================================================

design_pool <- matrix(
  c(0, 0.1,
    0.1, 0),
  nrow = 2,
  byrow = TRUE
)

rownames(design_pool) <- colnames(design_pool) <- c("DOM_proxy", "Microbe")

keep_dom <- ncol(data_diablo_pool$DOM_proxy)
keep_microbe <- 10

diablo_pool_res <- block.splsda(
  X = data_diablo_pool,
  Y = Y_diablo,
  ncomp = 2,
  keepX = list(
    DOM_proxy = c(keep_dom, keep_dom),
    Microbe = c(keep_microbe, keep_microbe)
  ),
  design = design_pool
)

print(diablo_pool_res)

# ============================================================
# 14. Basic plots - softer color palette
# ============================================================
# ============================================================
# Softer color palette for mixOmics plotLoadings
# ============================================================

pal_env <- c(
  "Bulk water" = "#3D5A80",  # muted dark blue
  "Biofilm"    = "#C8553D"   # muted brick red/orange
)

Y_diablo <- factor(Y_diablo, levels = c("Bulk water", "Biofilm"))
library(ggplot2)
library(dplyr)
library(tibble)

plot_loading_gg <- function(res, block_name, title_text, filename) {
  
  loading_df <- res$loadings[[block_name]][, 1, drop = FALSE] %>%
    as.data.frame() %>%
    rownames_to_column("Feature") %>%
    rename(Loading = 2) %>%
    filter(Loading != 0) %>%
    mutate(
      Outcome = ifelse(Loading > 0, "Biofilm", "Bulk water"),
      Feature = factor(Feature, levels = Feature[order(Loading)])
    )
  
  p <- ggplot(loading_df, aes(x = Loading, y = Feature, fill = Outcome)) +
    geom_col(width = 0.65) +
    scale_fill_manual(
      values = c("Bulk water" = "#3D5A80", "Biofilm" = "#C8553D")
    ) +
    labs(
      title = title_text,
      x = "",
      y = "",
      fill = "Outcome"
    ) +
    theme_classic(base_size = 18) +
    theme(
      plot.title = element_text(size = 26, face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 18, face = "bold", color = "black"),
      axis.text.x = element_text(size = 16, face = "bold", color = "black"),
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 15, face = "bold"),
      legend.position = "right"
    )
  
  print(p)
  
  ggsave(
    filename = filename,
    plot = p,
    width = 10,
    height = 7,
    dpi = 300,
    bg = "white"
  )
  return(p)
}

p_loading_dom <- plot_loading_gg(
  res = diablo_pool_res,
  block_name = "DOM_proxy",
  title_text = "Selected DOM proxy pools",
  filename = "DIABLO_selected_DOM_proxy_pools_loading_bold.png"
)

p_loading_microbe <- plot_loading_gg(
  res = diablo_pool_res,
  block_name = "Microbe",
  title_text = "Selected microbial genera",
  filename = "DIABLO_selected_microbial_genera_loading_bold.png"
)
# ============================================================
# 15. Extract selected variables
# ============================================================

dom_selected <- names(which(diablo_pool_res$loadings$DOM_proxy[, 1] != 0))
microbe_selected <- names(which(diablo_pool_res$loadings$Microbe[, 1] != 0))

cat("\nSelected DOM proxies:\n")
print(dom_selected)

cat("\nSelected microbial genera:\n")
print(microbe_selected)

selected_features_pool <- data.frame(
  Block = c(
    rep("DOM_proxy", length(dom_selected)),
    rep("Microbe", length(microbe_selected))
  ),
  Feature = c(dom_selected, microbe_selected)
)

write.csv(
  selected_features_pool,
  "DIABLO_proxy_selected_features_comp1.csv",
  row.names = FALSE
)

# 17. Group mean for selected DOM proxies
# Use raw relative abundance / proxy values for interpretation
# ============================================================

dom_sel_raw <- ms_pool_mat[, dom_selected, drop = FALSE]

dom_group_mean <- as.data.frame(dom_sel_raw) %>%
  tibble::rownames_to_column("MS_sample") %>%
  dplyr::mutate(environment = Y_diablo[MS_sample]) %>%
  tidyr::pivot_longer(
    cols = -c(MS_sample, environment),
    names_to = "DOM_proxy",
    values_to = "Value"
  ) %>%
  dplyr::group_by(DOM_proxy, environment) %>%
  dplyr::summarise(
    Mean = mean(Value, na.rm = TRUE),
    SD = sd(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = environment,
    values_from = c(Mean, SD)
  ) %>%
  dplyr::mutate(
    Higher_in = dplyr::case_when(
      `Mean_Biofilm` > `Mean_Bulk water` ~ "Biofilm",
      `Mean_Biofilm` < `Mean_Bulk water` ~ "Bulk water",
      TRUE ~ "similar"
    )
  )

print(dom_group_mean)

write.csv(
  dom_group_mean,
  "DIABLO_selected_DOM_proxy_group_mean.csv",
  row.names = FALSE
)

# ============================================================
# 18. Group mean for selected microbial genera
# Use relative abundance for interpretation
# ============================================================

microbe_sel_raw <- otu_mat_rel[, microbe_selected, drop = FALSE]

microbe_group_mean <- as.data.frame(microbe_sel_raw) %>%
  tibble::rownames_to_column("MS_sample") %>%
  dplyr::mutate(environment = Y_diablo[MS_sample]) %>%
  tidyr::pivot_longer(
    cols = -c(MS_sample, environment),
    names_to = "Microbe",
    values_to = "Relative_abundance"
  ) %>%
  dplyr::group_by(Microbe, environment) %>%
  dplyr::summarise(
    Mean = mean(Relative_abundance, na.rm = TRUE),
    SD = sd(Relative_abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = environment,
    values_from = c(Mean, SD)
  ) %>%
  dplyr::mutate(
    Higher_in = dplyr::case_when(
      `Mean_Biofilm` > `Mean_Bulk water` ~ "Biofilm",
      `Mean_Biofilm` < `Mean_Bulk water` ~ "Bulk water",
      TRUE ~ "similar"
    )
  )

print(microbe_group_mean)

write.csv(
  microbe_group_mean,
  "DIABLO_selected_Microbe_group_mean_proxy_model.csv",
  row.names = FALSE
)

# ============================================================
# 19. Bar plot of selected DOM proxies by group
# With Welch t-test significance stars
# ============================================================

sig_label_fun <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE ~ "ns"
  )
}

show_only_significant <- FALSE

# -----------------------------
# 19.1 Prepare DOM plot dataframe
# -----------------------------

dom_plot_df <- as.data.frame(dom_sel_raw) %>%
  tibble::rownames_to_column("MS_sample") %>%
  dplyr::mutate(
    environment = Y_diablo[MS_sample],
    environment = factor(environment, levels = c("Bulk water", "Biofilm"))
  ) %>%
  tidyr::pivot_longer(
    cols = -c(MS_sample, environment),
    names_to = "DOM_proxy",
    values_to = "Value"
  ) %>%
  dplyr::mutate(
    DOM_proxy = factor(DOM_proxy, levels = dom_selected)
  )

# -----------------------------
# 19.2 Welch t-test for DOM proxies
# -----------------------------

dom_ttest <- dom_plot_df %>%
  dplyr::group_by(DOM_proxy) %>%
  dplyr::summarise(
    p_value = tryCatch(
      t.test(Value ~ environment)$p.value,
      error = function(e) NA_real_
    ),
    bulk_water_mean = mean(Value[environment == "Bulk water"], na.rm = TRUE),
    biofilm_mean = mean(Value[environment == "Biofilm"], na.rm = TRUE),
    log2FC_Biofilm_vs_Bulk_water = log2((biofilm_mean + 1e-9) / (bulk_water_mean + 1e-9)),
    Higher_in = dplyr::case_when(
      biofilm_mean > bulk_water_mean ~ "Biofilm",
      biofilm_mean < bulk_water_mean ~ "Bulk water",
      TRUE ~ "similar"
    ),
    sig_label = sig_label_fun(p_value),
    .groups = "drop"
  )

print(dom_ttest)

write.csv(
  dom_ttest,
  "DIABLO_selected_DOM_proxy_relative_ttest.csv",
  row.names = FALSE
)

# -----------------------------
# 19.3 DOM bar summary: mean +/- SD
# -----------------------------

dom_bar_df <- dom_plot_df %>%
  dplyr::group_by(DOM_proxy, environment) %>%
  dplyr::summarise(
    Mean = mean(Value, na.rm = TRUE),
    SD = sd(Value, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------
# 19.4 DOM significance annotation position
# -----------------------------

dom_sig_df <- dom_bar_df %>%
  dplyr::mutate(
    upper = Mean + SD
  ) %>%
  dplyr::group_by(DOM_proxy) %>%
  dplyr::summarise(
    y_max_bar = max(upper, na.rm = TRUE),
    y_min_bar = min(Mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    dom_plot_df %>%
      dplyr::group_by(DOM_proxy) %>%
      dplyr::summarise(
        y_max_point = max(Value, na.rm = TRUE),
        y_min_point = min(Value, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "DOM_proxy"
  ) %>%
  dplyr::left_join(dom_ttest, by = "DOM_proxy") %>%
  dplyr::mutate(
    y_max = pmax(y_max_bar, y_max_point, na.rm = TRUE),
    y_min = pmin(y_min_bar, y_min_point, na.rm = TRUE),
    y_range = y_max - y_min,
    y_range = ifelse(is.na(y_range) | y_range == 0, abs(y_max), y_range),
    y_range = ifelse(is.na(y_range) | y_range == 0, 1, y_range),
    y_bracket = y_max + 0.15 * y_range,
    y_text = y_max + 0.25 * y_range,
    x_start = 1,
    x_end = 2,
    x_mid = 1.5
  )

if (show_only_significant) {
  dom_sig_df <- dom_sig_df %>%
    dplyr::filter(!is.na(p_value), p_value < 0.05)
}

print(dom_sig_df %>% dplyr::select(DOM_proxy, p_value, sig_label))

# -----------------------------
# 19.5 DOM bar plot with significance stars
# -----------------------------

p_dom_selected_bar <- ggplot(
  dom_bar_df,
  aes(x = environment, y = Mean, fill = environment)
) +
  geom_col(width = 0.6, alpha = 0.8, color = "black", linewidth = 0.3) +
  geom_errorbar(
    aes(ymin = Mean - SD, ymax = Mean + SD),
    width = 0.18,
    linewidth = 0.5
  ) +
  geom_jitter(
    data = dom_plot_df,
    aes(x = environment, y = Value),
    inherit.aes = FALSE,
    width = 0.08,
    size = 2,
    color = "black"
  ) +
  
  geom_segment(
    data = dom_sig_df,
    aes(x = x_start, xend = x_end, y = y_bracket, yend = y_bracket),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_segment(
    data = dom_sig_df,
    aes(x = x_start, xend = x_start, y = y_bracket, yend = y_bracket - 0.04 * y_range),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_segment(
    data = dom_sig_df,
    aes(x = x_end, xend = x_end, y = y_bracket, yend = y_bracket - 0.04 * y_range),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_text(
    data = dom_sig_df,
    aes(x = x_mid, y = y_text, label = sig_label),
    inherit.aes = FALSE,
    size = 5
  ) +

  facet_wrap(~ DOM_proxy, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("Bulk water" = "#F8766D", "Biofilm" = "#00BFC4")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.30))) +
  theme_ms_plot() +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey85", color = "black"),

    panel.spacing = grid::unit(1.1, "lines"),

    axis.text.x = element_text(size = 13, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    title = "Selected DOM proxy pools from DIABLO",
    x = "",

    y = ""
  )

print(p_dom_selected_bar)

ggsave(
  "DIABLO_selected_DOM_proxy_barplot_with_significance.png",
  p_dom_selected_bar,
  width = 12,
  height = 9,
  dpi = 300
)

# ============================================================
# 20. Bar plot of selected microbial genera by group
# With Welch t-test significance stars
# 20.1 Prepare microbe plot dataframe
# -----------------------------

microbe_plot_df <- as.data.frame(microbe_sel_raw) %>%
  tibble::rownames_to_column("MS_sample") %>%
  dplyr::mutate(
    environment = Y_diablo[MS_sample],
    environment = factor(environment, levels = c("Bulk water", "Biofilm"))
  ) %>%
  tidyr::pivot_longer(
    cols = -c(MS_sample, environment),
    names_to = "Microbe",
    values_to = "Relative_abundance"
  ) %>%
  dplyr::mutate(
    Microbe = factor(Microbe, levels = microbe_selected)
  )

# -----------------------------
# 20.2 Welch t-test for microbial genera
# -----------------------------

microbe_ttest <- microbe_plot_df %>%
  dplyr::group_by(Microbe) %>%
  dplyr::summarise(
    p_value = tryCatch(
      t.test(Relative_abundance ~ environment)$p.value,
      error = function(e) NA_real_
    ),
    bulk_water_mean = mean(Relative_abundance[environment == "Bulk water"], na.rm = TRUE),
    biofilm_mean = mean(Relative_abundance[environment == "Biofilm"], na.rm = TRUE),
    log2FC_Biofilm_vs_Bulk_water = log2((biofilm_mean + 1e-9) / (bulk_water_mean + 1e-9)),
    Higher_in = dplyr::case_when(
      biofilm_mean > bulk_water_mean ~ "Biofilm",
      biofilm_mean < bulk_water_mean ~ "Bulk water",
      TRUE ~ "similar"
    ),
    sig_label = sig_label_fun(p_value),
    .groups = "drop"
  )

print(microbe_ttest)

write.csv(
  microbe_ttest,
  "DIABLO_selected_Microbe_relative_ttest.csv",
  row.names = FALSE
)

# -----------------------------
# 20.3 Microbe bar summary: mean +/- SD
# -----------------------------

microbe_bar_df <- microbe_plot_df %>%
  dplyr::group_by(Microbe, environment) %>%
  dplyr::summarise(
    Mean = mean(Relative_abundance, na.rm = TRUE),
    SD = sd(Relative_abundance, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------
# 20.4 Microbe significance annotation position
# -----------------------------

microbe_sig_df <- microbe_bar_df %>%
  dplyr::mutate(
    upper = Mean + SD
  ) %>%
  dplyr::group_by(Microbe) %>%
  dplyr::summarise(
    y_max_bar = max(upper, na.rm = TRUE),
    y_min_bar = min(Mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    microbe_plot_df %>%
      dplyr::group_by(Microbe) %>%
      dplyr::summarise(
        y_max_point = max(Relative_abundance, na.rm = TRUE),
        y_min_point = min(Relative_abundance, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "Microbe"
  ) %>%
  dplyr::left_join(microbe_ttest, by = "Microbe") %>%
  dplyr::mutate(
    y_max = pmax(y_max_bar, y_max_point, na.rm = TRUE),
    y_min = pmin(y_min_bar, y_min_point, na.rm = TRUE),
    y_range = y_max - y_min,
    y_range = ifelse(is.na(y_range) | y_range == 0, abs(y_max), y_range),
    y_range = ifelse(is.na(y_range) | y_range == 0, 1, y_range),
    y_bracket = y_max + 0.15 * y_range,
    y_text = y_max + 0.25 * y_range,
    x_start = 1,
    x_end = 2,
    x_mid = 1.5
  )

if (show_only_significant) {
  microbe_sig_df <- microbe_sig_df %>%
    dplyr::filter(!is.na(p_value), p_value < 0.05)
}

print(microbe_sig_df %>% dplyr::select(Microbe, p_value, sig_label))

# -----------------------------
# 20.5 Microbe bar plot with significance stars
# -----------------------------

p_microbe_selected_bar <- ggplot(
  microbe_bar_df,
  aes(x = environment, y = Mean, fill = environment)
) +
  geom_col(width = 0.6, alpha = 0.8, color = "black", linewidth = 0.3) +
  geom_errorbar(
    aes(ymin = Mean - SD, ymax = Mean + SD),
    width = 0.18,
    linewidth = 0.5
  ) +
  geom_jitter(
    data = microbe_plot_df,
    aes(x = environment, y = Relative_abundance),
    inherit.aes = FALSE,
    width = 0.08,
    size = 2,
    color = "black"
  ) +
  
  geom_segment(
    data = microbe_sig_df,
    aes(x = x_start, xend = x_end, y = y_bracket, yend = y_bracket),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_segment(
    data = microbe_sig_df,
    aes(x = x_start, xend = x_start, y = y_bracket, yend = y_bracket - 0.04 * y_range),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_segment(
    data = microbe_sig_df,
    aes(x = x_end, xend = x_end, y = y_bracket, yend = y_bracket - 0.04 * y_range),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_text(
    data = microbe_sig_df,
    aes(x = x_mid, y = y_text, label = sig_label),
    inherit.aes = FALSE,
    size = 5
  ) +

  facet_wrap(~ Microbe, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("Bulk water" = "#F8766D", "Biofilm" = "#00BFC4")) +
  scale_y_continuous(

    labels = scales::label_number(accuracy = 0.01),
    expand = expansion(mult = c(0.05, 0.30))
  ) +
  theme_ms_plot() +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey85", color = "black"),

    panel.spacing = grid::unit(1.1, "lines"),

    axis.text.x = element_text(size = 13, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    title = "Selected microbial genera from DIABLO",
    x = "",

    y = ""
  )

print(p_microbe_selected_bar)

ggsave(
  "DIABLO_selected_Microbe_barplot_with_significance.png",
  p_microbe_selected_bar,
  width = 12,
  height = 10,
  dpi = 300
)

# ============================================================
# 22.7 Check total ion intensity for each sample
# With Welch t-test significance label
# ------------------------------------------------------------

tic_df <- ms_pool_long %>%
  dplyr::group_by(MS_sample) %>%
  dplyr::summarise(
    Total_intensity = sum(Intensity, na.rm = TRUE),
    Formula_count = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    environment = Y_diablo[MS_sample],
    environment = factor(environment, levels = c("Bulk water", "Biofilm")),
    MS_order = as.numeric(stringr::str_extract(MS_sample, "\\d+"))
  ) %>%
  dplyr::arrange(MS_order)

tic_p <- tryCatch(
  t.test(Total_intensity ~ environment, data = tic_df)$p.value,
  error = function(e) NA_real_
)

tic_sig_label <- sig_label_fun(tic_p)

tic_y_max <- max(tic_df$Total_intensity, na.rm = TRUE)
tic_y_min <- min(tic_df$Total_intensity, na.rm = TRUE)
tic_y_range <- tic_y_max - tic_y_min

if (tic_y_range == 0 | is.na(tic_y_range)) {
  tic_y_range <- tic_y_max
}

tic_sig_df <- data.frame(
  x_start = 1,
  x_end = 2,
  x_mid = 1.5,
  y_bracket = tic_y_max + 0.12 * tic_y_range,
  y_text = tic_y_max + 0.20 * tic_y_range,
  y_range = tic_y_range,
  sig_label = tic_sig_label,
  p_value = tic_p
)

print(tic_sig_df)

p_tic <- ggplot(
  tic_df,
  aes(x = environment, y = Total_intensity, fill = environment)
) +
  geom_boxplot(width = 0.6, alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.08, size = 2.5, color = "black") +
  
  geom_segment(
    data = tic_sig_df,
    aes(x = x_start, xend = x_end, y = y_bracket, yend = y_bracket),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_segment(
    data = tic_sig_df,
    aes(x = x_start, xend = x_start, y = y_bracket, yend = y_bracket - 0.04 * y_range),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_segment(
    data = tic_sig_df,
    aes(x = x_end, xend = x_end, y = y_bracket, yend = y_bracket - 0.04 * y_range),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_text(
    data = tic_sig_df,
    aes(x = x_mid, y = y_text, label = sig_label),
    inherit.aes = FALSE,
    size = 5
  ) +
  
  scale_fill_manual(values = c("Bulk water" = "#F8766D", "Biofilm" = "#00BFC4")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  theme_ms_plot() +
  theme(legend.position = "none") +
  labs(
    title = "Total FT-ICR-MS\nsignal intensity",
    x = "",
    y = "Total summed peak intensity"
  )

print(p_tic)

ggsave(
  "Total_FTICRMS_signal_intensity_boxplot_with_significance.png",
  p_tic,
  width = 5,
  height = 5,
  dpi = 300
)

write.csv(
  tic_df,
  "Total_FTICRMS_signal_intensity_by_sample.csv",
  row.names = FALSE
)
# Extra analysis: only sample 1/3/5 Bulk water vs sample 2/4/6 biofilm
# Reservoir 1 paired comparison
# ============================================================

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(scales)

# ------------------------------------------------------------
# 1. Define subset samples
# ------------------------------------------------------------

subset_samples <- paste("sample", 1:6)

# Check whether these samples exist
stopifnot(all(subset_samples %in% rownames(ms_pool_mat)))
stopifnot(all(subset_samples %in% rownames(otu_mat_rel)))
stopifnot(all(subset_samples %in% names(Y_diablo)))

# Metadata for sample 1-6
meta_sub <- matched_meta %>%
  dplyr::filter(MS_sample %in% subset_samples) %>%
  dplyr::arrange(MS_order) %>%
  dplyr::mutate(
    environment = factor(environment, levels = c("Bulk water", "Biofilm")),
    
    # pair_id connects Bulk water-Biofilm from the same location:
    # sample 1-2, sample 3-4, sample 5-6
    pair_id = paste0("Reservoir_", reservoir, "_Site_", site)
  )

print(meta_sub %>% dplyr::select(MS_sample, ID, environment, reservoir, site, pair_id))

# Check group numbers
print(table(meta_sub$environment))

# ------------------------------------------------------------
# 2. Helper function for significance labels
# ------------------------------------------------------------

sig_label_fun <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE ~ "ns"
  )
}

# ============================================================
# Part A. DOM proxy analysis using only sample 1-6
# ============================================================

# Use DIABLO-selected DOM proxies if available.
# Since your DOM block is small, this usually includes most/all DOM proxies.
dom_features_use <- intersect(dom_selected, colnames(ms_pool_mat))

# Fallback: use all DOM proxies if dom_selected is empty
if (length(dom_features_use) == 0) {
  dom_features_use <- colnames(ms_pool_mat)
}

dom_sub_mat <- ms_pool_mat[subset_samples, dom_features_use, drop = FALSE]

dom_sub_long <- as.data.frame(dom_sub_mat) %>%
  tibble::rownames_to_column("MS_sample") %>%
  dplyr::left_join(
    meta_sub %>%
      dplyr::select(MS_sample, environment, reservoir, site, pair_id),
    by = "MS_sample"
  ) %>%
  tidyr::pivot_longer(
    cols = all_of(dom_features_use),
    names_to = "DOM_proxy",
    values_to = "Value"
  ) %>%
  dplyr::mutate(
    environment = factor(environment, levels = c("Bulk water", "Biofilm")),
    DOM_proxy = factor(DOM_proxy, levels = dom_features_use)
  )

# ------------------------------------------------------------
# Paired t-test for DOM proxies
# ------------------------------------------------------------

dom_pair_wide <- dom_sub_long %>%
  dplyr::select(pair_id, DOM_proxy, environment, Value) %>%
  tidyr::pivot_wider(
    names_from = environment,
    values_from = Value
  )

dom_paired_stats <- dom_pair_wide %>%
  dplyr::group_by(DOM_proxy) %>%
  dplyr::summarise(
    n_pairs = sum(!is.na(`Bulk water`) & !is.na(Biofilm)),
    bulk_water_mean = mean(`Bulk water`, na.rm = TRUE),
    biofilm_mean = mean(Biofilm, na.rm = TRUE),
    mean_delta_Biofilm_minus_Bulk_water = mean(Biofilm - `Bulk water`, na.rm = TRUE),
    log2FC_Biofilm_vs_Bulk_water = log2((biofilm_mean + 1e-9) / (bulk_water_mean + 1e-9)),
    p_value = tryCatch(
      t.test(Biofilm, `Bulk water`, paired = TRUE)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    sig_label = sig_label_fun(p_value),
    Higher_in = dplyr::case_when(
      biofilm_mean > bulk_water_mean ~ "Biofilm",
      biofilm_mean < bulk_water_mean ~ "Bulk water",
      TRUE ~ "similar"
    )
  ) %>%
  dplyr::arrange(p_value)

print(dom_paired_stats)

write.csv(
  dom_paired_stats,
  "Subset_sample1to6_DOM_proxy_paired_ttest.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# DOM summary for plot
# ------------------------------------------------------------

dom_bar_sub <- dom_sub_long %>%
  dplyr::group_by(DOM_proxy, environment) %>%
  dplyr::summarise(
    Mean = mean(Value, na.rm = TRUE),
    SD = sd(Value, na.rm = TRUE),
    .groups = "drop"
  )

dom_sig_sub <- dom_sub_long %>%
  dplyr::group_by(DOM_proxy) %>%
  dplyr::summarise(
    y_max = max(Value, na.rm = TRUE),
    y_min = min(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(dom_paired_stats, by = "DOM_proxy") %>%
  dplyr::mutate(
    y_range = y_max - y_min,
    y_range = ifelse(is.na(y_range) | y_range == 0, abs(y_max), y_range),
    y_range = ifelse(is.na(y_range) | y_range == 0, 1, y_range),
    y_bracket = y_max + 0.15 * y_range,
    y_text = y_max + 0.25 * y_range,
    x_start = 1,
    x_end = 2,
    x_mid = 1.5
  )

# ------------------------------------------------------------
# DOM paired bar plot
# ------------------------------------------------------------

p_dom_sub <- ggplot() +
  geom_col(
    data = dom_bar_sub,
    aes(x = environment, y = Mean, fill = environment),
    width = 0.6,
    alpha = 0.8,
    color = "black",
    linewidth = 0.3
  ) +
  geom_errorbar(
    data = dom_bar_sub,
    aes(x = environment, ymin = Mean - SD, ymax = Mean + SD),
    width = 0.18,
    linewidth = 0.5
  ) +
  geom_line(
    data = dom_sub_long,
    aes(x = environment, y = Value, group = pair_id),
    color = "grey40",
    linewidth = 0.5,
    alpha = 0.8
  ) +
  geom_point(
    data = dom_sub_long,
    aes(x = environment, y = Value),
    color = "black",
    size = 2
  ) +
  geom_segment(
    data = dom_sig_sub,
    aes(x = x_start, xend = x_end, y = y_bracket, yend = y_bracket),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_segment(
    data = dom_sig_sub,
    aes(x = x_start, xend = x_start, y = y_bracket, yend = y_bracket - 0.04 * y_range),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_segment(
    data = dom_sig_sub,
    aes(x = x_end, xend = x_end, y = y_bracket, yend = y_bracket - 0.04 * y_range),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_text(
    data = dom_sig_sub,
    aes(x = x_mid, y = y_text, label = sig_label),
    inherit.aes = FALSE,
    size = 5
  ) +

  facet_wrap(~ DOM_proxy, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("Bulk water" = "#F8766D", "Biofilm" = "#00BFC4")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.30))) +
  theme_ms_plot() +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey85", color = "black"),

    panel.spacing = grid::unit(1.1, "lines"),
    axis.text.x = element_text(size = 13, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    title = "DOM proxy pools: Reservoir 1 paired comparison",
    subtitle = "sample 1/3/5 = Bulk water; sample 2/4/6 = Biofilm",
    x = "",

    y = ""
  )

print(p_dom_sub)

ggsave(
  "Subset_sample1to6_DOM_proxy_paired_barplot.png",
  p_dom_sub,
  width = 12,
  height = 9,
  dpi = 300
)

# ============================================================
# Part B. Microbial genera analysis using only sample 1-6
# ============================================================

# Use DIABLO-selected genera
microbe_features_use <- intersect(microbe_selected, colnames(otu_mat_rel))

if (length(microbe_features_use) == 0) {
  stop("No selected microbial genera were found in otu_mat_rel.")
}

microbe_sub_mat <- otu_mat_rel[subset_samples, microbe_features_use, drop = FALSE]

microbe_sub_long <- as.data.frame(microbe_sub_mat) %>%
  tibble::rownames_to_column("MS_sample") %>%
  dplyr::left_join(
    meta_sub %>%
      dplyr::select(MS_sample, environment, reservoir, site, pair_id),
    by = "MS_sample"
  ) %>%
  tidyr::pivot_longer(
    cols = all_of(microbe_features_use),
    names_to = "Microbe",
    values_to = "Relative_abundance"
  ) %>%
  dplyr::mutate(
    environment = factor(environment, levels = c("Bulk water", "Biofilm")),
    Microbe = factor(Microbe, levels = microbe_features_use)
  )

# ------------------------------------------------------------
# Paired t-test for microbial genera
# ------------------------------------------------------------

microbe_pair_wide <- microbe_sub_long %>%
  dplyr::select(pair_id, Microbe, environment, Relative_abundance) %>%
  tidyr::pivot_wider(
    names_from = environment,
    values_from = Relative_abundance
  )

microbe_paired_stats <- microbe_pair_wide %>%
  dplyr::group_by(Microbe) %>%
  dplyr::summarise(
    n_pairs = sum(!is.na(`Bulk water`) & !is.na(Biofilm)),
    bulk_water_mean = mean(`Bulk water`, na.rm = TRUE),
    biofilm_mean = mean(Biofilm, na.rm = TRUE),
    mean_delta_Biofilm_minus_Bulk_water = mean(Biofilm - `Bulk water`, na.rm = TRUE),
    log2FC_Biofilm_vs_Bulk_water = log2((biofilm_mean + 1e-9) / (bulk_water_mean + 1e-9)),
    p_value = tryCatch(
      t.test(Biofilm, `Bulk water`, paired = TRUE)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    sig_label = sig_label_fun(p_value),
    Higher_in = dplyr::case_when(
      biofilm_mean > bulk_water_mean ~ "Biofilm",
      biofilm_mean < bulk_water_mean ~ "Bulk water",
      TRUE ~ "similar"
    )
  ) %>%
  dplyr::arrange(p_value)

print(microbe_paired_stats)

write.csv(
  microbe_paired_stats,
  "Subset_sample1to6_Microbe_paired_ttest.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# Microbe summary for plot
# ------------------------------------------------------------

microbe_bar_sub <- microbe_sub_long %>%
  dplyr::group_by(Microbe, environment) %>%
  dplyr::summarise(
    Mean = mean(Relative_abundance, na.rm = TRUE),
    SD = sd(Relative_abundance, na.rm = TRUE),
    .groups = "drop"
  )

microbe_sig_sub <- microbe_sub_long %>%
  dplyr::group_by(Microbe) %>%
  dplyr::summarise(
    y_max = max(Relative_abundance, na.rm = TRUE),
    y_min = min(Relative_abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(microbe_paired_stats, by = "Microbe") %>%
  dplyr::mutate(
    y_range = y_max - y_min,
    y_range = ifelse(is.na(y_range) | y_range == 0, abs(y_max), y_range),
    y_range = ifelse(is.na(y_range) | y_range == 0, 1, y_range),
    y_bracket = y_max + 0.15 * y_range,
    y_text = y_max + 0.25 * y_range,
    x_start = 1,
    x_end = 2,
    x_mid = 1.5
  )

# ------------------------------------------------------------
# Microbe paired bar plot
# ------------------------------------------------------------

p_microbe_sub <- ggplot() +
  geom_col(
    data = microbe_bar_sub,
    aes(x = environment, y = Mean, fill = environment),
    width = 0.6,
    alpha = 0.8,
    color = "black",
    linewidth = 0.3
  ) +
  geom_errorbar(
    data = microbe_bar_sub,
    aes(x = environment, ymin = Mean - SD, ymax = Mean + SD),
    width = 0.18,
    linewidth = 0.5
  ) +
  geom_line(
    data = microbe_sub_long,
    aes(x = environment, y = Relative_abundance, group = pair_id),
    color = "grey40",
    linewidth = 0.5,
    alpha = 0.8
  ) +
  geom_point(
    data = microbe_sub_long,
    aes(x = environment, y = Relative_abundance),
    color = "black",
    size = 2
  ) +
  geom_segment(
    data = microbe_sig_sub,
    aes(x = x_start, xend = x_end, y = y_bracket, yend = y_bracket),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_segment(
    data = microbe_sig_sub,
    aes(x = x_start, xend = x_start, y = y_bracket, yend = y_bracket - 0.04 * y_range),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_segment(
    data = microbe_sig_sub,
    aes(x = x_end, xend = x_end, y = y_bracket, yend = y_bracket - 0.04 * y_range),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_text(
    data = microbe_sig_sub,
    aes(x = x_mid, y = y_text, label = sig_label),
    inherit.aes = FALSE,
    size = 5
  ) +

  facet_wrap(~ Microbe, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("Bulk water" = "#F8766D", "Biofilm" = "#00BFC4")) +
  scale_y_continuous(

    labels = scales::label_number(accuracy = 0.01),
    expand = expansion(mult = c(0.05, 0.30))
  ) +
  theme_ms_plot() +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey85", color = "black"),

    panel.spacing = grid::unit(1.1, "lines"),
    axis.text.x = element_text(size = 13, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    title = "Selected microbial genera: Reservoir 1 paired comparison",
    subtitle = "sample 1/3/5 = Bulk water; sample 2/4/6 = Biofilm",
    x = "",

    y = ""
  )

print(p_microbe_sub)

ggsave(
  "Subset_sample1to6_Microbe_paired_barplot.png",
  p_microbe_sub,
  width = 12,
  height = 10,
  dpi = 300
)

# Replace the original cimDiablo() section with this
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)

# ------------------------------------------------------------
# 1. Decide sample order
# ------------------------------------------------------------

# Option A: paired order, recommended
sample_order <- paste("sample", 1:10)

# If you only want sample 1-6, use this instead:
# sample_order <- paste("sample", 1:6)

# ------------------------------------------------------------
# 2. Order selected DOM proxies and microbes by loading weight
# ------------------------------------------------------------

dom_order <- dom_selected[
  order(abs(diablo_pool_res$loadings$DOM_proxy[dom_selected, 1]),
        decreasing = TRUE)
]

microbe_order <- microbe_selected[
  order(abs(diablo_pool_res$loadings$Microbe[microbe_selected, 1]),
        decreasing = TRUE)
]

# ------------------------------------------------------------
# 3. Extract selected features
# ------------------------------------------------------------

dom_heat <- data_diablo_pool$DOM_proxy[
  sample_order,
  dom_order,
  drop = FALSE
]

microbe_heat <- data_diablo_pool$Microbe[
  sample_order,
  microbe_order,
  drop = FALSE
]

heat_mat <- cbind(dom_heat, microbe_heat)

# ------------------------------------------------------------
# 4. Z-score each feature for heatmap
# ------------------------------------------------------------

zscore <- function(x) {
  if (sd(x, na.rm = TRUE) == 0 || is.na(sd(x, na.rm = TRUE))) {
    return(rep(0, length(x)))
  } else {
    return(as.numeric(scale(x)))
  }
}

heat_mat_z <- apply(heat_mat, 2, zscore)
rownames(heat_mat_z) <- rownames(heat_mat)

# ------------------------------------------------------------
# 5. Build annotation information
# ------------------------------------------------------------

feature_info <- data.frame(
  Feature = colnames(heat_mat_z),
  Block = c(
    rep("DOM proxy", length(dom_order)),
    rep("Microbe", length(microbe_order))
  )
)

sample_info <- data.frame(
  MS_sample = rownames(heat_mat_z),
  environment = as.character(Y_diablo[rownames(heat_mat_z)])
)

# ------------------------------------------------------------
# 6. Convert to long format
# ------------------------------------------------------------

heat_df <- as.data.frame(heat_mat_z) %>%
  tibble::rownames_to_column("MS_sample") %>%
  tidyr::pivot_longer(
    cols = -MS_sample,
    names_to = "Feature",
    values_to = "Z_score"
  ) %>%
  dplyr::left_join(feature_info, by = "Feature") %>%
  dplyr::left_join(sample_info, by = "MS_sample") %>%
  dplyr::mutate(
    MS_sample = factor(MS_sample, levels = rev(sample_order)),
    Feature = factor(Feature, levels = c(dom_order, microbe_order)),
    environment = factor(environment, levels = c("Bulk water", "Biofilm")),
    
    # Make long labels easier to read
    Feature_label = stringr::str_replace_all(as.character(Feature), "_", "\n")
  )

# ------------------------------------------------------------
# 7. Plot custom heatmap
# ------------------------------------------------------------

p_custom_heatmap <- ggplot(
  heat_df,
  aes(x = Feature, y = MS_sample, fill = Z_score)
) +
  geom_tile(color = "white", linewidth = 0.25) +
  facet_grid(
    . ~ Block,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_fill_gradient2(
    low = "#3D5A80",
    mid = "#F7F7F2",
    high = "#C8553D",
    midpoint = 0,
    name = "Z-score"
  ) +
  scale_x_discrete(
    labels = function(x) stringr::str_replace_all(x, "_", " ")
  ) +
  theme_ms_plot() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey85", color = "black"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,

      size = 11,
      face = "bold"
    ),

    axis.text.y = element_text(size = 13, face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),

    panel.spacing.x = grid::unit(1.2, "lines"),
    legend.position = "right"
  ) +
  labs(
    title = "DIABLO-selected DOM proxies and microbial genera",
    x = "",
    y = "Sample"
  )

print(p_custom_heatmap)

ggsave(
  "DIABLO_custom_heatmap_better_arrangement.png",
  p_custom_heatmap,
  width = 15,
  height = 7,
  dpi = 300
)

# ============================================================
# 23. Proxy-based Spearman correlation analysis
# Filter: |rho| > 0.6 and p < 0.05
# ============================================================
# This section is added at the end of the script.
# It uses DOM proxy-level variables, not molecule/formula-level classes.
# DOM proxy matrix: ms_pool_mat
# Microbial matrix: otu_mat_rel

cat("\nStarting proxy-based Spearman correlation analysis...\n")

# ------------------------------------------------------------
# 23.1 User-adjustable parameters
# ------------------------------------------------------------

spearman_rho_cutoff <- 0.6
spearman_p_cutoff <- 0.05

# Use all 10 matched samples by default.
# If you only want Reservoir 1 paired samples, replace this with: paste("sample", 1:6)
samples_for_proxy_spearman <- rownames(ms_pool_mat)

# TRUE = only correlate DOM proxies with DIABLO-selected microbial genera.
# This is recommended for thesis figures because the network will be cleaner.
# FALSE = use top variable genera from the full relative abundance table.
use_diablo_selected_microbes_only <- TRUE
top_n_microbes_for_spearman <- 30

# ------------------------------------------------------------
# 23.2 Prepare DOM proxy and microbe matrices
# ------------------------------------------------------------

samples_for_proxy_spearman <- samples_for_proxy_spearman[
  samples_for_proxy_spearman %in% rownames(ms_pool_mat) &
    samples_for_proxy_spearman %in% rownames(otu_mat_rel) &
    samples_for_proxy_spearman %in% names(Y_diablo)
]

if (length(samples_for_proxy_spearman) < 4) {
  stop("Too few matched samples for Spearman correlation. Please check sample names.")
}

# DOM proxies: use all proxy-level variables already constructed in ms_pool_mat,
# including CRAM, BDOC_PCL, and RDOC_to_BDOC.
dom_proxy_features_for_cor <- colnames(ms_pool_mat)

dom_proxy_var <- apply(
  ms_pool_mat[samples_for_proxy_spearman, dom_proxy_features_for_cor, drop = FALSE],
  2,
  function(x) var(x, na.rm = TRUE)
)

dom_proxy_features_for_cor <- names(dom_proxy_var)[
  is.finite(dom_proxy_var) & !is.na(dom_proxy_var) & dom_proxy_var > 0
]

if (length(dom_proxy_features_for_cor) == 0) {
  stop("No DOM proxy feature has non-zero variance for correlation analysis.")
}

# Microbes: use DIABLO-selected genera by default; otherwise use top variable genera.
if (
  isTRUE(use_diablo_selected_microbes_only) &&
    exists("microbe_selected") &&
    length(microbe_selected) > 0
) {
  microbe_features_for_cor <- intersect(microbe_selected, colnames(otu_mat_rel))
} else {
  microbe_var_all <- apply(
    otu_mat_rel[samples_for_proxy_spearman, , drop = FALSE],
    2,
    function(x) var(x, na.rm = TRUE)
  )
  microbe_var_all <- microbe_var_all[
    is.finite(microbe_var_all) & !is.na(microbe_var_all) & microbe_var_all > 0
  ]
  microbe_features_for_cor <- names(sort(microbe_var_all, decreasing = TRUE))[
    1:min(top_n_microbes_for_spearman, length(microbe_var_all))
  ]
}

# Remove zero-variance microbial features after selection.
microbe_var_for_cor <- apply(
  otu_mat_rel[samples_for_proxy_spearman, microbe_features_for_cor, drop = FALSE],
  2,
  function(x) var(x, na.rm = TRUE)
)

microbe_features_for_cor <- names(microbe_var_for_cor)[
  is.finite(microbe_var_for_cor) & !is.na(microbe_var_for_cor) & microbe_var_for_cor > 0
]

if (length(microbe_features_for_cor) == 0) {
  stop("No microbial feature has non-zero variance for correlation analysis.")
}

dom_proxy_cor_mat <- ms_pool_mat[
  samples_for_proxy_spearman,
  dom_proxy_features_for_cor,
  drop = FALSE
]

microbe_cor_mat <- otu_mat_rel[
  samples_for_proxy_spearman,
  microbe_features_for_cor,
  drop = FALSE
]

cat("\nSamples used for proxy-based Spearman correlation:\n")
print(samples_for_proxy_spearman)

cat("\nDOM proxies used:\n")
print(dom_proxy_features_for_cor)

cat("\nMicrobial genera used:\n")
print(microbe_features_for_cor)

# ------------------------------------------------------------
# 23.3 Safe Spearman correlation function
# ------------------------------------------------------------

safe_spearman_test <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  n_used <- length(x)
  
  if (n_used < 4 || length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(list(rho = NA_real_, p_value = NA_real_, n = n_used))
  }
  
  test_res <- tryCatch(
    suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE)),
    error = function(e) NULL
  )
  
  if (is.null(test_res)) {
    return(list(rho = NA_real_, p_value = NA_real_, n = n_used))
  }
  
  return(list(
    rho = unname(test_res$estimate),
    p_value = test_res$p.value,
    n = n_used
  ))
}

# ------------------------------------------------------------
# 23.4 Calculate all pairwise proxy-microbe Spearman correlations
# ------------------------------------------------------------

cor_proxy_rows <- list()
row_i <- 1

for (dp in colnames(dom_proxy_cor_mat)) {
  for (mb in colnames(microbe_cor_mat)) {
    res <- safe_spearman_test(
      dom_proxy_cor_mat[, dp],
      microbe_cor_mat[, mb]
    )
    
    cor_proxy_rows[[row_i]] <- data.frame(
      DOM_proxy = dp,
      Microbe = mb,
      n = res$n,
      Spearman_rho = res$rho,
      p_value = res$p_value,
      stringsAsFactors = FALSE
    )
    row_i <- row_i + 1
  }
}

cor_proxy_df <- dplyr::bind_rows(cor_proxy_rows) %>%
  dplyr::mutate(
    abs_rho = abs(Spearman_rho),
    p_adj_BH = p.adjust(p_value, method = "BH"),
    Correlation_direction = dplyr::case_when(
      is.na(Spearman_rho) ~ NA_character_,
      Spearman_rho > 0 ~ "Positive",
      Spearman_rho < 0 ~ "Negative",
      TRUE ~ "Zero"
    )
  ) %>%
  dplyr::arrange(desc(abs_rho), p_value)

# ------------------------------------------------------------
# 23.5 Add group-level means to help biological interpretation
# ------------------------------------------------------------

dom_proxy_group_for_cor <- as.data.frame(dom_proxy_cor_mat) %>%
  tibble::rownames_to_column("MS_sample") %>%
  dplyr::mutate(
    environment = as.character(Y_diablo[MS_sample])
  ) %>%
  tidyr::pivot_longer(
    cols = -c(MS_sample, environment),
    names_to = "DOM_proxy",
    values_to = "DOM_proxy_value"
  ) %>%
  dplyr::group_by(DOM_proxy, environment) %>%
  dplyr::summarise(
    DOM_proxy_mean = mean(DOM_proxy_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = environment,
    values_from = DOM_proxy_mean,
    names_prefix = "DOM_mean_"
  ) %>%
  dplyr::mutate(
    DOM_proxy_higher_in = dplyr::case_when(
      `DOM_mean_Biofilm` > `DOM_mean_Bulk water` ~ "Biofilm",
      `DOM_mean_Biofilm` < `DOM_mean_Bulk water` ~ "Bulk water",
      TRUE ~ "similar"
    )
  )

microbe_group_for_cor <- as.data.frame(microbe_cor_mat) %>%
  tibble::rownames_to_column("MS_sample") %>%
  dplyr::mutate(
    environment = as.character(Y_diablo[MS_sample])
  ) %>%
  tidyr::pivot_longer(
    cols = -c(MS_sample, environment),
    names_to = "Microbe",
    values_to = "Relative_abundance"
  ) %>%
  dplyr::group_by(Microbe, environment) %>%
  dplyr::summarise(
    Microbe_mean = mean(Relative_abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = environment,
    values_from = Microbe_mean,
    names_prefix = "Microbe_mean_"
  ) %>%
  dplyr::mutate(
    Microbe_higher_in = dplyr::case_when(
      `Microbe_mean_Biofilm` > `Microbe_mean_Bulk water` ~ "Biofilm",
      `Microbe_mean_Biofilm` < `Microbe_mean_Bulk water` ~ "Bulk water",
      TRUE ~ "similar"
    )
  )

cor_proxy_df <- cor_proxy_df %>%
  dplyr::left_join(dom_proxy_group_for_cor, by = "DOM_proxy") %>%
  dplyr::left_join(microbe_group_for_cor, by = "Microbe")

# ------------------------------------------------------------
# 23.6 Strict filtering: |rho| > 0.6 and p < 0.05
# ------------------------------------------------------------

cor_proxy_filtered <- cor_proxy_df %>%
  dplyr::filter(
    !is.na(Spearman_rho),
    abs_rho > spearman_rho_cutoff,
    p_value < spearman_p_cutoff
  ) %>%
  dplyr::arrange(p_value, desc(abs_rho))

cat("\nTop proxy-based Spearman correlations before filtering:\n")
print(head(cor_proxy_df, 30))

cat("\nFiltered proxy-based Spearman correlations: |rho| > 0.6 and p < 0.05\n")
print(cor_proxy_filtered)

write.csv(
  cor_proxy_df,
  "Proxy_based_DOM_proxy_Microbe_Spearman_all_pairs.csv",
  row.names = FALSE
)

write.csv(
  cor_proxy_filtered,
  "Proxy_based_DOM_proxy_Microbe_Spearman_filtered_absrho_gt_0.6_p_lt_0.05.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 23.7 Heatmap for filtered significant correlations
# ------------------------------------------------------------

short_taxon_label <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace(x, "^.*\\|g__", "g__")
  x <- stringr::str_replace(x, "^.*;g__", "g__")
  x <- stringr::str_replace(x, "^.*g__", "g__")
  x <- stringr::str_replace_all(x, "_", " ")
  stringr::str_trunc(x, width = 45)
}

if (nrow(cor_proxy_filtered) > 0) {
  heatmap_proxy_df <- cor_proxy_filtered %>%
    dplyr::mutate(
      DOM_proxy = factor(
        DOM_proxy,
        levels = rev(unique(DOM_proxy[order(DOM_proxy)]))
      ),
      Microbe = factor(
        Microbe,
        levels = unique(Microbe[order(Microbe)])
      ),
      label_text = paste0(
        "rho=", round(Spearman_rho, 2),
        "\np=", signif(p_value, 2)
      )
    )
  
  p_proxy_spearman_heatmap <- ggplot(
    heatmap_proxy_df,
    aes(x = Microbe, y = DOM_proxy, fill = Spearman_rho)
  ) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = label_text), size = 3) +
    scale_fill_gradient2(
      low = "#3D5A80",
      mid = "#F7F7F2",
      high = "#C8553D",
      midpoint = 0,
      limits = c(-1, 1),
      name = "Spearman\nrho"
    ) +
    scale_x_discrete(labels = short_taxon_label) +
    scale_y_discrete(labels = function(x) stringr::str_replace_all(x, "_", "\n")) +

    theme_ms_plot(base_size = 15) +
    theme(
      panel.grid = element_blank(),

      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11, face = "bold"),

      axis.text.y = element_text(size = 12, face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(face = "bold", hjust = 0.5)
    ) +
    labs(
      title = "Proxy-based DOM-microbe Spearman correlations",
      subtitle = "Only correlations with |rho| > 0.6 and p < 0.05 are shown",
      x = "Microbial genus",
      y = "DOM proxy"
    )
  
  print(p_proxy_spearman_heatmap)
  
  ggsave(
    "Proxy_based_DOM_proxy_Microbe_Spearman_filtered_heatmap.png",
    p_proxy_spearman_heatmap,
    width = max(11, 0.65 * length(unique(cor_proxy_filtered$Microbe))),
    height = max(6, 0.60 * length(unique(cor_proxy_filtered$DOM_proxy))),
    dpi = 300
  )
} else {
  cat("\nNo correlations passed |rho| > 0.6 and p < 0.05. Filtered heatmap was not generated.\n")
}

# ------------------------------------------------------------
# 23.8 Simple bipartite network plot for filtered correlations
# ------------------------------------------------------------

if (nrow(cor_proxy_filtered) > 0) {
  dom_nodes <- unique(cor_proxy_filtered$DOM_proxy)
  microbe_nodes <- unique(cor_proxy_filtered$Microbe)
  
  dom_node_df <- data.frame(
    node = dom_nodes,
    node_label = stringr::str_replace_all(dom_nodes, "_", " "),
    x = 1,
    y = seq_along(dom_nodes),
    stringsAsFactors = FALSE
  )
  
  microbe_node_df <- data.frame(
    node = microbe_nodes,
    node_label = short_taxon_label(microbe_nodes),
    x = 2,
    y = seq_along(microbe_nodes),
    stringsAsFactors = FALSE
  )
  
  edge_plot_df <- cor_proxy_filtered %>%
    dplyr::left_join(
      dom_node_df %>% dplyr::select(DOM_proxy = node, x, y),
      by = "DOM_proxy"
    ) %>%
    dplyr::rename(x_start = x, y_start = y) %>%
    dplyr::left_join(
      microbe_node_df %>% dplyr::select(Microbe = node, x, y),
      by = "Microbe"
    ) %>%
    dplyr::rename(x_end = x, y_end = y)
  
  p_proxy_spearman_network <- ggplot() +
    geom_segment(
      data = edge_plot_df,
      aes(
        x = x_start,
        y = y_start,
        xend = x_end,
        yend = y_end,
        color = Correlation_direction,
        linewidth = abs_rho
      ),
      alpha = 0.75
    ) +
    geom_point(
      data = dom_node_df,
      aes(x = x, y = y),
      size = 4,
      color = "#A8D5BA"
    ) +
    geom_point(
      data = microbe_node_df,
      aes(x = x, y = y),
      size = 4,
      color = "#A9CCE3"
    ) +
    geom_text(
      data = dom_node_df,
      aes(x = x - 0.04, y = y, label = node_label),
      hjust = 1,
      size = 3.5
    ) +
    geom_text(
      data = microbe_node_df,
      aes(x = x + 0.04, y = y, label = node_label),
      hjust = 0,
      size = 3.2
    ) +
    scale_color_manual(
      values = c("Positive" = "#C8553D", "Negative" = "#3D5A80"),
      name = "Correlation"
    ) +
    scale_linewidth_continuous(range = c(0.4, 1.5), name = "|rho|") +
    coord_cartesian(xlim = c(0.5, 2.5), clip = "off") +

    theme_void(base_size = 15) +
    theme(
      legend.position = "bottom",
      plot.margin = margin(14, 130, 14, 130),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 13, face = "bold")
    ) +
    annotate("text", x = 1, y = max(dom_node_df$y) + 0.6, label = "DOM proxies", fontface = "bold") +
    annotate("text", x = 2, y = max(microbe_node_df$y) + 0.6, label = "Microbial genera", fontface = "bold") +
    labs(
      title = "Proxy-based DOM-microbe correlation network",
      subtitle = "Spearman |rho| > 0.6 and p < 0.05"
    )
  
  print(p_proxy_spearman_network)
  
  ggsave(
    "Proxy_based_DOM_proxy_Microbe_Spearman_filtered_network.png",
    p_proxy_spearman_network,
    width = 13,
    height = max(7, 0.45 * max(length(dom_nodes), length(microbe_nodes))),
    dpi = 300
  )
}

cat("\nProxy-based Spearman correlation analysis finished.\n")
cat("Output files saved:\n")
cat("1. Proxy_based_DOM_proxy_Microbe_Spearman_all_pairs.csv\n")
cat("2. Proxy_based_DOM_proxy_Microbe_Spearman_filtered_absrho_gt_0.6_p_lt_0.05.csv\n")
cat("3. Proxy_based_DOM_proxy_Microbe_Spearman_filtered_heatmap.png, if significant correlations exist\n")
cat("4. Proxy_based_DOM_proxy_Microbe_Spearman_filtered_network.png, if significant correlations exist\n")

cat("\nAll outputs were saved in:\n")
print(normalizePath(output_dir, winslash = "/", mustWork = FALSE))
