# Run through run_all.R; see README.md for inputs and figure mapping.
# Figures 2b and 2d use only the ten metadata-matched samples.
# Use sample-level PICRUSt2 input here, never the habitat-aggregated table.
# Heatmaps use clipped row z scores of log1p percentage abundance.
# Spearman stars use unadjusted p values; the CSV also reports BH-adjusted p values.
# Predicted pathways indicate functional potential, not measured carbon fluxes.
theme_corr_plot <- function(base_size = 16, base_family = "Arial") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(size = base_size + 3, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = base_size - 1, face = "bold", hjust = 0.5),
      axis.title = element_text(size = base_size + 1, face = "bold", colour = "black"),
      axis.text = element_text(size = base_size - 2, face = "bold", colour = "black"),
      legend.title = element_text(size = base_size - 1, face = "bold"),
      legend.text = element_text(size = base_size - 2, face = "bold"),
      panel.grid = element_blank()
    )
}
# 0. Load packages
# ============================================================

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)

# ============================================================
# 1. File paths
# ============================================================

otu_file <- data_file("total_featuretable.csv")

meta_file <- data_file("MataData.csv")

picrust2_sample_file <- data_file("picrust2.sample.pathways.relative.xls")

# ============================================================
# 2. Helper functions
# ============================================================

std_name <- function(x) {
  tolower(gsub("[^A-Za-z0-9]", "", x))
}

row_z_score <- function(x) {
  if (all(is.na(x)) || sd(x, na.rm = TRUE) == 0) {
    return(rep(0, length(x)))
  } else {
    return((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  }
}

clip_value <- function(x, low = -2, high = 2) {
  pmin(pmax(x, low), high)
}

shorten_label <- function(x, n = 50) {
  x <- gsub("&beta;", "\u03b2", x)
  x <- gsub("&alpha;", "\u03b1", x)
  x <- gsub("superpathway of ", "", x, ignore.case = TRUE)
  x <- gsub("biosynthesis", "biosyn.", x, ignore.case = TRUE)
  x <- gsub("degradation", "degrad.", x, ignore.case = TRUE)
  x <- gsub("fermentation", "ferm.", x, ignore.case = TRUE)
  x <- gsub("pathway", "path.", x, ignore.case = TRUE)
  x <- gsub("\\s+", " ", x)
  
  ifelse(
    nchar(x) > n,
    paste0(substr(x, 1, n - 3), "..."),
    x
  )
}

read_any_table <- function(file) {
  if (grepl("\\.csv$", file, ignore.case = TRUE)) {
    readr::read_csv(file, show_col_types = FALSE)
  } else {
    readr::read_tsv(file, show_col_types = FALSE)
  }
}

get_cluster_order <- function(data, row_col, sample_col, value_col, sample_order) {
  wide <- data %>%
    dplyr::select(tidyselect::all_of(c(row_col, sample_col, value_col))) %>%
    tidyr::pivot_wider(
      names_from = tidyselect::all_of(sample_col),
      values_from = tidyselect::all_of(value_col),
      values_fill = 0
    )
  
  rn <- wide[[row_col]]
  mat <- as.matrix(wide[, sample_order, drop = FALSE])
  rownames(mat) <- rn
  
  if (nrow(mat) >= 2) {
    rownames(mat)[hclust(dist(mat))$order]
  } else {
    rn
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
    MS_sample = trimws(`sample ID in MS`),
    MS_order = as.numeric(stringr::str_extract(MS_sample, "\\d+")),
    environment = tolower(trimws(as.character(environment))),
    environment = dplyr::recode(
      environment,
      "water" = "Bulk water",
      "biofilm" = "Biofilm"
    ),
    environment = factor(environment, levels = c("Bulk water", "Biofilm"))
  ) %>%
  dplyr::filter(!is.na(MS_sample)) %>%
  dplyr::arrange(MS_order)

sample_order <- matched_meta$MS_sample

print(matched_meta %>% dplyr::select(ID, MS_sample, MS_order, environment))
print(table(matched_meta$environment))

# ============================================================
# 4. Process sample-level 16S genus table
# ============================================================

otu_raw <- readr::read_csv(otu_file, show_col_types = FALSE)

colnames(otu_raw) <- trimws(colnames(otu_raw))
otu_raw$Taxonomy <- trimws(otu_raw$Taxonomy)

missing_otu_cols <- setdiff(matched_meta$ID, colnames(otu_raw))

if (length(missing_otu_cols) > 0) {
  stop(
    paste0(
      "Metadata sample IDs missing from the feature table: ",
      paste(missing_otu_cols, collapse = ", ")
    )
  )
}

otu_clean <- otu_raw %>%
  dplyr::filter(!is.na(Taxonomy)) %>%
  dplyr::filter(Taxonomy != "") %>%
  dplyr::filter(Taxonomy != "Others") %>%
  dplyr::filter(!stringr::str_detect(Taxonomy, stringr::regex("chloroplast|mitochondria", ignore_case = TRUE)))

taxa_long <- otu_clean %>%
  dplyr::select(Taxonomy, tidyselect::all_of(matched_meta$ID)) %>%
  tidyr::pivot_longer(
    cols = -Taxonomy,
    names_to = "ID",
    values_to = "Abundance"
  ) %>%
  dplyr::mutate(
    Abundance = as.numeric(Abundance)
  ) %>%
  dplyr::left_join(
    matched_meta %>% dplyr::select(ID, MS_sample, MS_order, environment),
    by = "ID"
  ) %>%
  dplyr::filter(!is.na(MS_sample)) %>%
  dplyr::group_by(MS_sample, MS_order, environment, Taxonomy) %>%
  dplyr::summarise(
    Abundance = sum(Abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(MS_sample) %>%
  dplyr::mutate(
    RelAbund = Abundance / sum(Abundance, na.rm = TRUE) * 100
  ) %>%
  dplyr::ungroup()

taxa_check <- taxa_long %>%
  dplyr::group_by(MS_sample) %>%
  dplyr::summarise(total = sum(RelAbund, na.rm = TRUE), .groups = "drop")

print(taxa_check)

# ============================================================
# 5. Process sample-level PICRUSt2 pathway table
# ============================================================

func_raw <- read_any_table(picrust2_sample_file)

colnames(func_raw) <- trimws(colnames(func_raw))
colnames(func_raw)[1] <- "Pathway"

desc_col <- names(func_raw)[tolower(names(func_raw)) %in% c("description", "desc", "pathway_description")]

if (length(desc_col) > 0 && desc_col[1] != "description") {
  colnames(func_raw)[colnames(func_raw) == desc_col[1]] <- "description"
}

if (!"description" %in% colnames(func_raw)) {
  func_raw$description <- func_raw$Pathway
}

func_raw <- func_raw %>%
  dplyr::mutate(
    Pathway = as.character(Pathway),
    description = as.character(description)
  )

non_sample_cols <- c("Pathway", "description", "classification", "taxon", "Taxon")
candidate_cols <- setdiff(colnames(func_raw), non_sample_cols)

sample_map_id <- matched_meta %>%
  dplyr::transmute(
    Func_col_std = std_name(ID),
    ID,
    MS_sample,
    MS_order,
    environment
  )

sample_map_ms <- matched_meta %>%
  dplyr::transmute(
    Func_col_std = std_name(MS_sample),
    ID,
    MS_sample,
    MS_order,
    environment
  )

func_col_map <- tibble(
  Func_col = candidate_cols,
  Func_col_std = std_name(candidate_cols)
) %>%
  dplyr::left_join(sample_map_id, by = "Func_col_std") %>%
  dplyr::bind_rows(
    tibble(
      Func_col = candidate_cols,
      Func_col_std = std_name(candidate_cols)
    ) %>%
      dplyr::left_join(sample_map_ms, by = "Func_col_std")
  ) %>%
  dplyr::filter(!is.na(MS_sample)) %>%
  dplyr::distinct(Func_col, .keep_all = TRUE)

if (nrow(func_col_map) < 3) {
  stop(
    paste0(
      "Too few sample columns were identified in the sample-level PICRUSt2 table.\n",
      "Check that PICRUSt2 column names match metadata ID or MS_sample.\n",
      "Available PICRUSt2 columns: ",
      paste(colnames(func_raw), collapse = ", ")
    )
  )
}

print(func_col_map)

func_long_all <- func_raw %>%
  dplyr::select(Pathway, description, tidyselect::all_of(func_col_map$Func_col)) %>%
  tidyr::pivot_longer(
    cols = tidyselect::all_of(func_col_map$Func_col),
    names_to = "Func_col",
    values_to = "AbsAbund"
  ) %>%
  dplyr::mutate(
    AbsAbund = as.numeric(AbsAbund)
  ) %>%
  dplyr::left_join(
    func_col_map %>% dplyr::select(Func_col, MS_sample, MS_order, environment),
    by = "Func_col"
  ) %>%
  dplyr::filter(!is.na(MS_sample)) %>%
  dplyr::group_by(MS_sample, MS_order, environment, Pathway, description) %>%
  dplyr::summarise(
    AbsAbund = sum(AbsAbund, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(MS_sample) %>%
  dplyr::mutate(
    RelAbund = AbsAbund / sum(AbsAbund, na.rm = TRUE) * 100
  ) %>%
  dplyr::ungroup()

func_check <- func_long_all %>%
  dplyr::group_by(MS_sample) %>%
  dplyr::summarise(total = sum(RelAbund, na.rm = TRUE), .groups = "drop")

print(func_check)

# ============================================================
# 6. Select carbon/DOC-related predicted pathways
# ============================================================

carbon_keywords <- c(
  "glycolysis",
  "tca",
  "tricarboxylic",
  "fermentation",
  "pyruvate",
  "acetate",
  "lactate",
  "formate",
  "succinate",
  "fumarate",
  "malate",
  "glyoxylate",
  "citrate",
  "methan",
  "methane",
  "methanol",
  "formaldehyde",
  "photosynth",
  "calvin",
  "carbon fixation",
  "wood",
  "ljungdahl",
  "reductive tca",
  "rtca",
  "hydroxypropionate",
  "aromatic",
  "benzoate",
  "phenyl",
  "catechol",
  "fatty acid",
  "sucrose",
  "starch",
  "cellulose",
  "xylan",
  "chitin",
  "glycogen",
  "respiration"
)

carbon_pattern <- paste(carbon_keywords, collapse = "|")

func_focus <- func_long_all %>%
  dplyr::filter(
    grepl(carbon_pattern, description, ignore.case = TRUE) |
      grepl(carbon_pattern, Pathway, ignore.case = TRUE)
  )

if (nrow(func_focus) == 0) {
  stop("No carbon-related pathways matched the specified keywords.")
}

label_df <- func_focus %>%
  dplyr::distinct(Pathway, description) %>%
  dplyr::mutate(Function_label = shorten_label(description, 55))

label_df$Function_label <- make.unique(label_df$Function_label)

func_focus <- func_focus %>%
  dplyr::left_join(label_df %>% dplyr::select(Pathway, Function_label), by = "Pathway")

# ============================================================
# Version 1:
# Sample 1-10 Microbial taxa + PICRUSt2 function joint heatmap
# ============================================================

top_taxa_n <- 20
top_func_n <- 25

# ---------- top taxa ----------
top_taxa <- taxa_long %>%
  dplyr::group_by(Taxonomy) %>%
  dplyr::summarise(
    MeanRel = mean(RelAbund, na.rm = TRUE),
    VarRel = var(log1p(RelAbund), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(MeanRel)) %>%
  dplyr::slice_head(n = top_taxa_n) %>%
  dplyr::pull(Taxonomy)

taxa_heat <- taxa_long %>%
  dplyr::filter(Taxonomy %in% top_taxa) %>%
  dplyr::mutate(LogRel = log1p(RelAbund)) %>%
  dplyr::group_by(Taxonomy) %>%
  dplyr::mutate(
    z_score = row_z_score(LogRel),
    z_score = clip_value(z_score, -2, 2)
  ) %>%
  dplyr::ungroup()

taxa_order <- get_cluster_order(
  data = taxa_heat,
  row_col = "Taxonomy",
  sample_col = "MS_sample",
  value_col = "z_score",
  sample_order = sample_order
)

taxa_heat <- taxa_heat %>%
  dplyr::mutate(
    MS_sample = factor(MS_sample, levels = sample_order),
    Taxonomy = factor(Taxonomy, levels = taxa_order),
    environment = factor(environment, levels = c("Bulk water", "Biofilm"))
  )

# ---------- top functions ----------
top_funcs <- func_focus %>%
  dplyr::group_by(Pathway, Function_label) %>%
  dplyr::summarise(
    MeanRel = mean(RelAbund, na.rm = TRUE),
    VarRel = var(log1p(RelAbund), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(MeanRel)) %>%
  dplyr::slice_head(n = top_func_n) %>%
  dplyr::pull(Pathway)

func_heat <- func_focus %>%
  dplyr::filter(Pathway %in% top_funcs) %>%
  dplyr::mutate(LogRel = log1p(RelAbund)) %>%
  dplyr::group_by(Function_label) %>%
  dplyr::mutate(
    z_score = row_z_score(LogRel),
    z_score = clip_value(z_score, -2, 2)
  ) %>%
  dplyr::ungroup()

func_order <- get_cluster_order(
  data = func_heat,
  row_col = "Function_label",
  sample_col = "MS_sample",
  value_col = "z_score",
  sample_order = sample_order
)

func_heat <- func_heat %>%
  dplyr::mutate(
    MS_sample = factor(MS_sample, levels = sample_order),
    Function_label = factor(Function_label, levels = func_order),
    environment = factor(environment, levels = c("Bulk water", "Biofilm"))
  )

# ---------- environment annotation ----------
env_anno <- matched_meta %>%
  dplyr::mutate(
    MS_sample = factor(MS_sample, levels = sample_order),
    environment = factor(environment, levels = c("Bulk water", "Biofilm"))
  )

env_colors <- c(
  "Bulk water" = "#4E79A7",
  "Biofilm" = "#C8553D"
)

p_env <- ggplot(env_anno, aes(x = MS_sample, y = "Environment", fill = environment)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(
    values = env_colors,
    labels = c("Bulk water", "Biofilm"),
    name = "Environment"
  ) +
  scale_x_discrete(labels = function(x) gsub("sample\\s*", "S", x, ignore.case = TRUE)) +
  labs(x = NULL, y = NULL) +

  theme_void(base_size = 15) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13, face = "bold"),
    axis.text.x = element_blank(),
    plot.margin = margin(0, 10, 0, 90)
  )

# ---------- taxa heatmap ----------
p_taxa <- ggplot(taxa_heat, aes(x = MS_sample, y = Taxonomy, fill = z_score)) +
  geom_tile(color = "white", linewidth = 0.35) +
  scale_fill_gradient2(
    low = "#3D5A80",
    mid = "#F7F7F2",
    high = "#C8553D",
    midpoint = 0,
    limits = c(-2, 2),
    name = "Taxa\nz-score"
  ) +
  scale_x_discrete(labels = function(x) gsub("sample\\s*", "S", x, ignore.case = TRUE)) +
  labs(
    title = "A. Microbial taxa",
    x = NULL,
    y = NULL
  ) +

  theme_corr_plot(base_size = 15) +
  theme(
    axis.text.x = element_blank(),

    axis.text.y = element_text(size = 12, face = "bold.italic", color = "gray15"),
    panel.grid = element_blank(),
    plot.margin = margin(5, 10, 5, 10)
  )

# ---------- function heatmap ----------
p_func <- ggplot(func_heat, aes(x = MS_sample, y = Function_label, fill = z_score)) +
  geom_tile(color = "white", linewidth = 0.35) +
  scale_fill_gradient2(
    low = "#3D5A80",
    mid = "#F7F7F2",
    high = "#C8553D",
    midpoint = 0,
    limits = c(-2, 2),
    name = "Function\nz-score"
  ) +
  scale_x_discrete(labels = function(x) gsub("sample\\s*", "S", x, ignore.case = TRUE)) +
  labs(
    title = "B. PICRUSt2-predicted carbon-cycling functions",
    x = NULL,
    y = NULL
  ) +

  theme_corr_plot(base_size = 15) +
  theme(

    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),

    axis.text.y = element_text(size = 11, face = "bold", color = "gray15"),
    panel.grid = element_blank(),
    plot.margin = margin(5, 10, 10, 10)
  )

p_joint_sample <- p_env / p_taxa / p_func +
  plot_layout(heights = c(0.18, 1.05, 1.45)) +
  plot_annotation(
    title = "Microbial Community and PICRUSt2-predicted Functional Potential",
    subtitle = "Sample-level comparison across Sample 1-10",
    theme = theme(

      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, face = "bold", hjust = 0.5, color = "gray35")
    )
  )

p_joint_sample

ggsave(
  filename = safe_out("V1_sample_level_taxa_function_joint_heatmap.png"),
  plot = p_joint_sample,

  width = 15,
  height = 12,
  dpi = 300,
  bg = "white"
)

# ============================================================
# Version 2:
# Taxa x predicted function Spearman correlation heatmap
# ============================================================

top_taxa_cor_n <- 20
top_func_cor_n <- 25

# ---------- select taxa by variance ----------
top_taxa_cor <- taxa_long %>%
  dplyr::mutate(LogRel = log1p(RelAbund)) %>%
  dplyr::group_by(Taxonomy) %>%
  dplyr::summarise(
    VarRel = var(LogRel, na.rm = TRUE),
    MeanRel = mean(RelAbund, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::filter(!is.na(VarRel), VarRel > 0) %>%
  dplyr::arrange(dplyr::desc(VarRel), dplyr::desc(MeanRel)) %>%
  dplyr::slice_head(n = top_taxa_cor_n) %>%
  dplyr::pull(Taxonomy)

# ---------- select functions by variance ----------
top_func_cor <- func_focus %>%
  dplyr::mutate(LogRel = log1p(RelAbund)) %>%
  dplyr::group_by(Function_label) %>%
  dplyr::summarise(
    VarRel = var(LogRel, na.rm = TRUE),
    MeanRel = mean(RelAbund, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::filter(!is.na(VarRel), VarRel > 0) %>%
  dplyr::arrange(dplyr::desc(VarRel), dplyr::desc(MeanRel)) %>%
  dplyr::slice_head(n = top_func_cor_n) %>%
  dplyr::pull(Function_label)

# ---------- taxa matrix: taxa x samples ----------
taxa_mat_df <- taxa_long %>%
  dplyr::filter(Taxonomy %in% top_taxa_cor) %>%
  dplyr::select(Taxonomy, MS_sample, RelAbund) %>%
  tidyr::pivot_wider(
    names_from = MS_sample,
    values_from = RelAbund,
    values_fill = 0
  )

taxa_mat <- as.data.frame(taxa_mat_df)
rownames(taxa_mat) <- taxa_mat$Taxonomy
taxa_mat$Taxonomy <- NULL

taxa_mat <- taxa_mat[, sample_order, drop = FALSE]

# ---------- function matrix: function x samples ----------
func_mat_df <- func_focus %>%
  dplyr::filter(Function_label %in% top_func_cor) %>%
  dplyr::group_by(Function_label, MS_sample) %>%
  dplyr::summarise(
    RelAbund = sum(RelAbund, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = MS_sample,
    values_from = RelAbund,
    values_fill = 0
  )

func_mat <- as.data.frame(func_mat_df)
rownames(func_mat) <- func_mat$Function_label
func_mat$Function_label <- NULL

func_mat <- func_mat[, sample_order, drop = FALSE]

# ---------- correlation ----------
cor_results <- expand.grid(
  Taxon = rownames(taxa_mat),
  Function = rownames(func_mat),
  stringsAsFactors = FALSE
)

cor_results$rho <- NA_real_
cor_results$p_value <- NA_real_

for (i in seq_len(nrow(cor_results))) {
  tx <- cor_results$Taxon[i]
  fn <- cor_results$Function[i]
  
  x <- as.numeric(taxa_mat[tx, sample_order])
  y <- as.numeric(func_mat[fn, sample_order])
  
  if (sd(x, na.rm = TRUE) == 0 || sd(y, na.rm = TRUE) == 0) {
    cor_results$rho[i] <- NA
    cor_results$p_value[i] <- NA
  } else {
    test <- suppressWarnings(
      cor.test(x, y, method = "spearman", exact = FALSE)
    )
    cor_results$rho[i] <- as.numeric(test$estimate)
    cor_results$p_value[i] <- test$p.value
  }
}

cor_results <- cor_results %>%
  dplyr::filter(!is.na(rho)) %>%
  dplyr::mutate(
    p_adj = p.adjust(p_value, method = "BH"),

    sig = dplyr::case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

write_csv(cor_results, safe_out("taxa_function_spearman_correlation_table.csv"))

# ---------- cluster order ----------
rho_wide <- cor_results %>%
  dplyr::select(Taxon, Function, rho) %>%
  tidyr::pivot_wider(
    names_from = Function,
    values_from = rho,
    values_fill = 0
  )

rho_mat <- as.data.frame(rho_wide)
rownames(rho_mat) <- rho_mat$Taxon
rho_mat$Taxon <- NULL
rho_mat <- as.matrix(rho_mat)

taxon_order_cor <- if (nrow(rho_mat) >= 2) {
  rownames(rho_mat)[hclust(dist(rho_mat))$order]
} else {
  rownames(rho_mat)
}

function_order_cor <- if (ncol(rho_mat) >= 2) {
  colnames(rho_mat)[hclust(dist(t(rho_mat)))$order]
} else {
  colnames(rho_mat)
}

cor_results <- cor_results %>%
  dplyr::mutate(
    Taxon = factor(Taxon, levels = taxon_order_cor),
    Function = factor(Function, levels = function_order_cor)
  )

# ---------- wrap function labels ----------
wrap_label <- function(x, width = 28) {
  sapply(x, function(y) paste(strwrap(y, width = width), collapse = "\n"))
}

# ---------- plot ----------
p_cor <- ggplot(cor_results, aes(x = Function, y = Taxon, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(
    aes(label = sig),

    size = 4,
    fontface = "bold",
    color = "black"
  ) +
  scale_fill_gradient2(
    low = "#3D5A80",
    mid = "#F7F7F2",
    high = "#C8553D",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman\nrho"
  ) +
  scale_x_discrete(labels = function(x) wrap_label(x, width = 26)) +
  labs(
    title = "Taxa-Function Association Heatmap",
    subtitle = "Spearman correlation between microbial taxa and PICRUSt2-predicted carbon-cycling pathways",
    x = "PICRUSt2-predicted carbon-cycling pathways",
    y = "Microbial taxa"
  ) +

  theme_corr_plot(base_size = 15) +
  theme(
    plot.subtitle = element_text(face = "bold", hjust = 0.5, color = "gray35"),

    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 11,
      face = "bold",
      color = "gray15"
    ),

    axis.text.y = element_text(
      size = 12,
      face = "bold.italic",
      color = "gray15"
    ),
    axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
    panel.grid = element_blank(),
    plot.margin = margin(18, 30, 18, 18)
  )

p_cor

ggsave(
  filename = safe_out("V2_sample_level_taxa_function_correlation_heatmap.png"),
  plot = p_cor,

  width = 16,
  height = 11,
  dpi = 300,
  bg = "white"
)

