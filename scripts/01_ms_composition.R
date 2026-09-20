# Run through run_all.R; see README.md for inputs and figure mapping.
# Figures 1a-d, S3 and S4. Pool rules are sequential and mutually exclusive.
# CRAM is a compositional proxy, not a direct measurement of persistence.
# RDOC = CRAM; BDOC = proteins + carbohydrates + lipids (relative signal).
library(readxl)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(scales)
library(patchwork)
library(readr)
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
# 1. Read Sample 1-10
# ============================================================

available_sheets <- excel_sheets(file)
wanted_sheets <- paste("sample", 1:10)

sheets <- available_sheets[tolower(available_sheets) %in% wanted_sheets]

if (!setequal(tolower(sheets), wanted_sheets)) {
  stop("Expected all ten sheets named sample 1 through sample 10.")
}

data_list <- map(
  sheets,
  ~ read_excel(file, sheet = .x) %>%
    mutate(Sample = tolower(.x))
)

df <- bind_rows(data_list)

sample_levels <- paste("sample", 1:10)

# ============================================================
# 2. Clean columns and calculate ratios / AImod
# ============================================================

for (nm in c("N", "S", "P")) {
  if (!nm %in% names(df)) {
    df[[nm]] <- 0
  }
}

numeric_cols <- intersect(
  c(
    "C", "H", "O", "N", "S", "P",
    "H_plus_1", "HC_ratio", "OC_ratio",
    "DBE", "AI_mod", "AI_denominator",
    "Observed Intens"
  ),
  names(df)
)

df <- df %>%
  mutate(across(all_of(numeric_cols), ~ suppressWarnings(as.numeric(.x)))) %>%
  mutate(
    N = ifelse(is.na(N), 0, N),
    S = ifelse(is.na(S), 0, S),
    P = ifelse(is.na(P), 0, P)
  )

if (!"H_plus_1" %in% names(df)) {
  df$H_plus_1 <- df$H + 1
}

if (!"HC_ratio" %in% names(df)) {
  df$HC_ratio <- df$H_plus_1 / df$C
}

if (!"OC_ratio" %in% names(df)) {
  df$OC_ratio <- df$O / df$C
}

if (!"DBE" %in% names(df)) {
  df$DBE <- 1 + (2 * df$C - df$H_plus_1 + df$N) / 2
}

df <- df %>%
  mutate(
    AI_denominator = C - 0.5 * O - S - N - P,
    AI_mod = ifelse(
      is.finite(AI_denominator) & AI_denominator > 0,
      (1 + C - 0.5 * O - S - 0.5 * (H_plus_1 + N + P)) / AI_denominator,
      NA_real_
    )
  )

# ============================================================
# 3. Add environment
# ============================================================

df <- df %>%
  mutate(
    Sample = factor(Sample, levels = sample_levels),
    environment = case_when(
      Sample %in% paste("sample", c(1, 3, 5, 7, 9)) ~ "Bulk water",
      Sample %in% paste("sample", c(2, 4, 6, 8, 10)) ~ "Biofilm",
      TRUE ~ NA_character_
    ),
    environment = factor(environment, levels = c("Bulk water", "Biofilm"))
  )

# ============================================================
# 4. Functional pool classification

# ============================================================

df_pool <- df %>%
  mutate(
    H_C = HC_ratio,
    O_C = OC_ratio,
    
    DBE_C = ifelse(C > 0, DBE / C, NA_real_),
    DBE_H = ifelse(H_plus_1 > 0, DBE / H_plus_1, NA_real_),
    DBE_O = ifelse(O > 0, DBE / O, NA_real_),
    N_C   = ifelse(C > 0, N / C, NA_real_),
    P_C   = ifelse(C > 0, P / C, NA_real_)
  ) %>%
  mutate(
    pool = case_when(
      is.finite(DBE_C) & is.finite(DBE_H) & is.finite(DBE_O) &
        DBE_C >= 0.30 & DBE_C <= 0.68 &
        DBE_H >= 0.20 & DBE_H <= 0.95 &
        DBE_O >= 0.77 & DBE_O <= 1.75 ~ "CRAM",

      is.finite(AI_mod) & AI_mod >= 0.67 ~ "Condensed Aromatics",
      
      is.finite(O_C) & is.finite(H_C) &
        O_C >= 0.8 & H_C >= 1.65 & H_C < 2.7 & N == 0 ~ "Carbohydrates",
      
      is.finite(O_C) & is.finite(H_C) & is.finite(N_C) & is.finite(P_C) &
        (
          (O_C > 0.12 & O_C <= 0.6 & H_C > 0.9 & H_C < 2.5 &
             N_C >= 0.126 & N_C <= 0.7 & P_C < 0.17) |
            (O_C > 0.6 & O_C <= 1.0 & H_C > 1.2 & H_C < 2.5 &
               N_C > 0.2 & N_C <= 0.7 & P_C < 0.17)
        ) ~ "Proteins",
      
      is.finite(O_C) & is.finite(H_C) & is.finite(N_C) & is.finite(P_C) &
        O_C <= 0.6 & H_C >= 1.32 & N_C <= 0.125 & P_C < 0.35 ~ "Lipids",
      
      is.finite(O_C) & is.finite(H_C) &
        O_C >= 0.67 & O_C <= 0.97 & H_C >= 0.53 & H_C <= 1.5 ~ "Tannins",
      
      is.finite(O_C) & is.finite(H_C) &
        O_C >= 0.25 & O_C <= 0.67 & H_C >= 0.75 & H_C <= 1.5 ~ "Lignins",
      
      is.finite(O_C) & is.finite(H_C) &
        O_C >= 0.00 & O_C <= 0.29 & H_C >= 1.0 & H_C <= 1.6 ~ "Unsaturated Hydrocarbons",
      
      TRUE ~ "Other"
    )
  )

table(df_pool$pool)

# ============================================================
# 5. Pool relative abundance
#    Fixed: normalize pool names and force stacked bars to 100%
# ============================================================

int_col <- "Observed Intens"

pool_levels <- c(
  "CRAM", 
  "Tannins", 
  "Condensed Aromatics", 
  "Lignins",
  "Unsaturated Hydrocarbons", 
  "Carbohydrates", 
  "Proteins", 
  "Lipids", 
  "Other"
)

df_pool <- df_pool %>%
  mutate(
    pool = as.character(pool),
    pool = case_when(
      pool %in% c("Condensed_Aromatics", "Condensed Aromatics") ~ "Condensed Aromatics",
      pool %in% c("Unsaturated_Hydrocarbons", "Unsaturated Hydrocarbons") ~ "Unsaturated Hydrocarbons",
      pool %in% pool_levels ~ pool,
      is.na(pool) ~ "Other",
      TRUE ~ "Other"
    ),
    pool = factor(pool, levels = pool_levels),
    Sample = factor(Sample, levels = sample_levels),
    environment = factor(environment, levels = c("Bulk water", "Biofilm"))
  )

df_comp <- df_pool %>%
  group_by(Sample, environment, pool) %>%
  summarise(
    I = sum(.data[[int_col]], na.rm = TRUE),
    count = n(),
    .groups = "drop"
  ) %>%
  group_by(Sample, environment) %>%
  mutate(
    total_I = sum(I, na.rm = TRUE),
    relI = ifelse(total_I > 0, I / total_I, 0)
  ) %>%
  ungroup()

check_sample_sum <- df_comp %>%
  group_by(Sample, environment) %>%
  summarise(sum_relI = sum(relI, na.rm = TRUE), .groups = "drop")

print(check_sample_sum)

write_csv(
  df_comp,
  file.path(out_dir, "Molecular_pool_relative_abundance_by_sample.csv")
)

# ============================================================
# 6. RDOC / BDOC
# ============================================================

df_rdoc_bdoc <- df_comp %>%
  mutate(
    RDOC = ifelse(pool == "CRAM", relI, 0),
    BDOC = ifelse(pool %in% c("Proteins", "Carbohydrates", "Lipids"), relI, 0)
  ) %>%
  group_by(Sample, environment) %>%
  summarise(
    RDOC = sum(RDOC, na.rm = TRUE),
    BDOC = sum(BDOC, na.rm = TRUE),
    RDOC_to_BDOC = ifelse(BDOC > 0, RDOC / BDOC, NA_real_),
    .groups = "drop"
  ) %>%
  mutate(Sample = factor(Sample, levels = sample_levels))

write_csv(df_rdoc_bdoc, file.path(out_dir, "RDOC_BDOC_bySample.csv"))

# RDOC vs BDOC barplot
df_long <- df_rdoc_bdoc %>%
  dplyr::select(Sample, RDOC, BDOC) %>%
  pivot_longer(cols = c(RDOC, BDOC), names_to = "Metric", values_to = "Value")

p_rdoc_bdoc <- ggplot(df_long, aes(x = Sample, y = Value, fill = Metric)) +
  geom_col(position = "dodge", width = 0.8) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_ms_plot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "RDOC (CRAM) vs BDOC (Protein + Carbohydrate + Lipid)",
    x = "Sample",
    y = "Relative abundance"
  )

p_rdoc_bdoc
ggsave(
  file.path(out_dir, "RDOC_vs_BDOC.png"),
  p_rdoc_bdoc,
  width = 12,
  height = 6,
  dpi = 300
)
# 9. VK density by environment
# ============================================================

df_vk <- df_pool %>%
  mutate(
    H_C = as.numeric(H_C),
    O_C = as.numeric(O_C),
    Intensity = as.numeric(.data[[int_col]])
  ) %>%
  filter(
    is.finite(H_C),
    is.finite(O_C),
    is.finite(Intensity),
    !is.na(environment),
    O_C >= 0,
    O_C <= 1.2,
    H_C >= 0,
    H_C <= 2.5
  )

p_vk_density <- ggplot(df_vk, aes(x = O_C, y = H_C)) +
  geom_density_2d_filled() +
  facet_wrap(~ environment) +
  coord_cartesian(xlim = c(0, 1.2), ylim = c(0, 2.5)) +

  scale_x_continuous(breaks = c(0, 0.4, 0.8, 1.2)) +
  theme_ms_plot() +
  theme(

    panel.spacing.x = grid::unit(1.2, "lines"),
    plot.margin = margin(10, 18, 14, 18)
  ) +
  labs(
    title = "Van Krevelen density distribution",
    x = "O/C",
    y = "H/C"
  )

p_vk_density
ggsave(
  file.path(out_dir, "VK_density_Biofilm_vs_Bulk.png"),
  p_vk_density,
  width = 12,
  height = 6,
  dpi = 300
)

# ============================================================
# 10. VK map by environment
# ============================================================

p_vkmap_env <- ggplot(df_vk, aes(x = O_C, y = H_C)) +
  geom_point(
    aes(size = Intensity, color = pool),
    alpha = 0.35
  ) +
  facet_wrap(~ environment) +
  coord_cartesian(xlim = c(0, 1.2), ylim = c(0, 2.5)) +

  scale_x_continuous(breaks = c(0, 0.4, 0.8, 1.2)) +
  scale_size_continuous(
    range = c(0.2, 2.5),
    labels = scientific
  ) +
  theme_ms_plot() +
  labs(
    title = "Van Krevelen Map of DOM Molecular Formulae",
    x = "O/C",
    y = "H/C",
    color = "Molecular pool",
    size = "Intensity"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank(),

    panel.spacing.x = grid::unit(1.2, "lines"),
    plot.margin = margin(10, 18, 14, 18)
  )

p_vkmap_env

ggsave(
  file.path(out_dir, "VK_map_by_environment.png"),
  p_vkmap_env,
  width = 13,
  height = 6.5,
  dpi = 300
)

# ============================================================
# 11. VK map for each sample
# ============================================================

df_vkmap_sample <- df_vk %>%
  mutate(
    Sample = factor(
      Sample,
      levels = c(
        "sample 1", "sample 3", "sample 5", "sample 7", "sample 9",
        "sample 2", "sample 4", "sample 6", "sample 8", "sample 10"
      )
    )
  )

p_vkmap_sample <- ggplot(df_vkmap_sample, aes(x = O_C, y = H_C)) +
  geom_point(
    aes(size = Intensity, color = pool),
    alpha = 0.35
  ) +

  facet_wrap(~ Sample, ncol = 5) +
  coord_cartesian(xlim = c(0, 1.2), ylim = c(0, 2.5)) +

  scale_x_continuous(breaks = c(0, 0.4, 0.8, 1.2)) +
  scale_size_continuous(
    range = c(0.2, 2.2),
    labels = scientific
  ) +
  theme_ms_plot() +
  labs(
    title = "Van Krevelen Maps of Individual Samples",
    x = "O/C",
    y = "H/C",
    color = "Molecular pool",
    size = "Intensity"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank(),

    axis.text.x = element_text(angle = 0, hjust = 0.5),

    panel.spacing = grid::unit(1.1, "lines"),
    plot.margin = margin(10, 18, 14, 18)
  )

p_vkmap_sample

ggsave(
  file.path(out_dir, "VK_map_each_sample.png"),
  p_vkmap_sample,
  width = 13,
  height = 14,
  dpi = 300
)
sample_env <- df_pool %>%
  distinct(Sample, environment) %>%
  filter(!is.na(Sample), !is.na(environment))

# Complete pools within observed sample-habitat pairs, never across habitats.
stopifnot(nrow(sample_env) == 10L, !anyDuplicated(sample_env$Sample))

df_proxy_grid <- sample_env %>%
  tidyr::expand_grid(pool = pool_levels)

df_proxy_rel_sample <- df_proxy_grid %>%
  left_join(
    df_comp %>%
      mutate(pool = as.character(pool)) %>%
      dplyr::select(Sample, environment, pool, I, count, relI),
    by = c("Sample", "environment", "pool")
  ) %>%
  mutate(
    I = replace_na(I, 0),
    count = replace_na(count, 0),
    relI = replace_na(relI, 0),
    Sample = factor(Sample, levels = sample_levels),
    environment = factor(environment, levels = c("Bulk water", "Biofilm")),
    pool = factor(pool, levels = pool_levels)
  ) %>%
  arrange(Sample, pool)

df_proxy_rel_sample %>%
  mutate(relI_percent = percent(relI, accuracy = 0.01)) %>%
  dplyr::select(Sample, environment, pool, relI, relI_percent) %>%
  print(n = Inf)

write_csv(
  df_proxy_rel_sample,
  file.path(out_dir, "Proxy_relative_abundance_by_sample.csv")
)

df_proxy_rel_env <- df_proxy_rel_sample %>%
  group_by(environment, pool) %>%
  summarise(
    mean_relI = mean(relI, na.rm = TRUE),
    sd_relI = sd(relI, na.rm = TRUE),
    se_relI = sd_relI / sqrt(n()),
    .groups = "drop"
  )

df_proxy_rel_env %>%
  mutate(mean_percent = percent(mean_relI, accuracy = 0.01)) %>%
  dplyr::select(environment, pool, mean_relI, mean_percent, sd_relI, se_relI) %>%
  print(n = Inf)

write_csv(
  df_proxy_rel_env,
  file.path(out_dir, "Proxy_relative_abundance_by_environment.csv")
)

p_proxy_sample <- ggplot(
  df_proxy_rel_sample,
  aes(x = Sample, y = relI, fill = pool)
) +
  geom_col(width = 0.85) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_ms_plot() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    title = "Relative abundance of molecular pools in each sample",
    x = "Sample",
    y = "Relative abundance",
    fill = "Molecular pool"
  )

p_proxy_sample

ggsave(
  file.path(out_dir, "Proxy_relative_abundance_by_sample.png"),
  p_proxy_sample,
  width = 12,
  height = 6,
  dpi = 300
)

df_proxy_rel_env_stacked <- df_proxy_rel_sample %>%
  group_by(environment, pool) %>%
  summarise(
    mean_relI = mean(relI, na.rm = TRUE),
    sd_relI = sd(relI, na.rm = TRUE),
    se_relI = sd_relI / sqrt(n()),
    .groups = "drop"
  ) %>%
  group_by(environment) %>%
  mutate(
    sum_mean_relI = sum(mean_relI, na.rm = TRUE),
    mean_relI_norm = ifelse(sum_mean_relI > 0, mean_relI / sum_mean_relI, 0)
  ) %>%
  ungroup()

check_env_sum <- df_proxy_rel_env_stacked %>%
  group_by(environment) %>%
  summarise(sum_relI = sum(mean_relI_norm, na.rm = TRUE), .groups = "drop")

print(check_env_sum)

write_csv(
  df_proxy_rel_env_stacked,
  file.path(out_dir, "Proxy_relative_abundance_Biofilm_vs_BulkWater_mean_stacked_FIXED.csv")
)

p_proxy_env_stacked <- ggplot(
  df_proxy_rel_env_stacked,
  aes(x = environment, y = mean_relI_norm, fill = pool)
) +
  geom_col(width = 0.75, color = "white", linewidth = 0.25) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  theme_ms_plot() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.6)
  ) +
  labs(
    title = "Overall molecular-pool relative abundance",
    subtitle = "Mean of sample-level relative abundance\nCondensed aromatics: AImod >= 0.67",
    x = "",
    y = "Relative abundance",
    fill = "Molecular pool"
  )

p_proxy_env_stacked

ggsave(
  file.path(out_dir, "Proxy_relative_abundance_Biofilm_vs_BulkWater_mean_stacked_FIXED.png"),
  p_proxy_env_stacked,
  width = 9,
  height = 6,
  dpi = 300
)
