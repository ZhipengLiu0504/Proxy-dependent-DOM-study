# Run through run_all.R; see README.md for inputs and figure mapping.
# Figures 2a, 2c and S5 use the supplied habitat-aggregated input tables.
# S5 retains Others during normalization; Fig. 2a excludes Others before normalization.
# Pathway predictions describe functional potential, not measured activity.
# The .xls pathway input is tab-delimited text, not an Excel workbook.
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
theme_genomic_plot <- function(base_size = 16, base_family = "Arial") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(size = base_size + 3, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = base_size - 1, face = "bold", hjust = 0.5),
      axis.title = element_text(size = base_size + 1, face = "bold", colour = "black"),
      axis.text = element_text(size = base_size - 2, face = "bold", colour = "black"),
      legend.title = element_text(size = base_size - 1, face = "bold"),
      legend.text = element_text(size = base_size - 2, face = "bold"),
      panel.grid.minor = element_blank()
    )
}

file <- data_file("WandB_featuretable.csv")

df <- read_csv(file)

print(head(df))
print(colnames(df))

# ============================================================

# Taxonomy = genus / taxa name
# Bulk = bulk water abundance
# Benthic = biofilm abundance
# Tax_detail = full taxonomy
# ============================================================

df_clean <- df %>%
  mutate(
    Taxonomy = as.character(Taxonomy),
    Tax_detail = as.character(Tax_detail),
    Bulk = as.numeric(Bulk),
    Benthic = as.numeric(Benthic)
  ) %>%
  filter(!is.na(Taxonomy)) %>%
  filter(Taxonomy != "") %>%
  filter(!grepl("chloroplast|mitochondria", Taxonomy, ignore.case = TRUE)) %>%
  filter(!grepl("chloroplast|mitochondria", Tax_detail, ignore.case = TRUE))

df_long <- df_clean %>%
  dplyr::select(Taxonomy, Bulk, Benthic) %>%
  pivot_longer(
    cols = c(Bulk, Benthic),
    names_to = "Environment",
    values_to = "Abundance"
  ) %>%
  mutate(
    Environment = case_when(
      Environment == "Bulk" ~ "Bulk water",
      Environment == "Benthic" ~ "Biofilm",
      TRUE ~ Environment
    )
  )

df_long <- df_long %>%
  group_by(Environment) %>%
  mutate(
    RelAbund = Abundance / sum(Abundance, na.rm = TRUE) * 100
  ) %>%
  ungroup()

check_sum <- df_long %>%
  group_by(Environment) %>%
  summarise(total_percent = sum(RelAbund, na.rm = TRUE))

print(check_sum)

top_n <- 10

top_taxa <- df_long %>%
  filter(Taxonomy != "Others") %>%
  group_by(Taxonomy) %>%
  summarise(
    mean_abund = mean(RelAbund, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abund)) %>%
  slice_head(n = top_n) %>%
  pull(Taxonomy)

print(top_taxa)

df_plot <- df_long %>%
  mutate(
    Taxa_plot = ifelse(Taxonomy %in% top_taxa, Taxonomy, "Others")
  ) %>%
  group_by(Environment, Taxa_plot) %>%
  summarise(
    RelAbund = sum(RelAbund, na.rm = TRUE),
    .groups = "drop"
  )

check_plot_sum <- df_plot %>%
  group_by(Environment) %>%
  summarise(total_percent = sum(RelAbund, na.rm = TRUE))

print(check_plot_sum)

taxa_levels <- unique(c(top_taxa, "Others"))

df_plot <- df_plot %>%
  mutate(
    Environment = factor(
      Environment,
      levels = c("Bulk water", "Biofilm")
    ),
    Taxa_plot = factor(
      Taxa_plot,
      levels = taxa_levels
    )
  )

taxa_levels_new <- c("Others", top_taxa)

df_plot <- df_plot %>%
  mutate(
    Taxa_plot = factor(Taxa_plot, levels = taxa_levels_new),
    Environment = factor(Environment, levels = c("Bulk water", "Biofilm"))
  )

main_colors <- c(
  "#4E79A7",  # muted blue
  "#59A14F",  # green
  "#F28E2B",  # orange
  "#76B7B2",  # teal
  "#E15759",  # coral red
  "#B07AA1",  # purple
  "#EDC948",  # soft yellow
  "#9C755F",  # brown
  "#BAB0AC",  # warm grey
  "#86BCB6",  # pale teal
  "#A0CBE8",  # light blue
  "#FFBE7D",  # peach
  "#8CD17D",  # light green
  "#D4A6C8",  # lavender pink
  "#B6992D"   # ochre
)

taxa_colors <- c(
  "Others" = "#ACD6FF",
  setNames(main_colors[1:length(top_taxa)], top_taxa)
)

legend_order <- c(top_taxa, "Others")

p <- ggplot(df_plot, aes(x = Environment, y = RelAbund, fill = Taxa_plot)) +
  geom_col(
    width = 0.55,
    color = "white",
    linewidth = 0.35
  ) +
  scale_fill_manual(
    values = taxa_colors,
    breaks = legend_order
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%"),
    expand = c(0, 0)
  ) +
  labs(
    title = "Microbial Community Composition",
    x = NULL,
    y = "Relative abundance",
    fill = "Taxa"
  ) +
  theme_genomic_plot() +
  theme(
    plot.title = element_text(margin = margin(b = 18)),
    
    axis.title.y = element_text(margin = margin(r = 12)),

    axis.text.x = element_text(
      size = 15,
      face = "bold",
      color = "gray20"
    ),

    axis.text.y = element_text(
      size = 12,
      face = "bold",
      color = "gray25"
    ),
    
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(
      color = "gray88",
      linewidth = 0.35
    ),
    
    axis.line.x = element_line(
      color = "gray35",
      linewidth = 0.5
    ),
    
    axis.line.y = element_line(
      color = "gray35",
      linewidth = 0.5
    ),

    legend.key.size = unit(0.55, "cm"),
    
    legend.position = "right",
    legend.background = element_blank(),
    legend.key = element_blank(),
    
    plot.margin = margin(15, 25, 15, 15)
  )

p

ggsave(
  filename = safe_out("microbial_community_composition_stacked_bar.png"),
  plot = p,
  width = 8,
  height = 6.5,
  dpi = 300,
  bg = "white"
)

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

file <- data_file("WandB_featuretable.csv")

df <- read_csv(file)

# ============================================================

# Taxonomy = genus name
# Bulk = bulk water abundance
# Benthic = Biofilm abundance
# Tax_detail = full taxonomy
# ============================================================

df_clean <- df %>%
  mutate(
    Taxonomy = as.character(Taxonomy),
    Tax_detail = as.character(Tax_detail),
    Bulk = as.numeric(Bulk),
    Benthic = as.numeric(Benthic)
  ) %>%
  filter(!is.na(Taxonomy)) %>%
  filter(Taxonomy != "") %>%
  filter(Taxonomy != "Others") %>%
  filter(!grepl("chloroplast|mitochondria", Taxonomy, ignore.case = TRUE)) %>%
  filter(!grepl("chloroplast|mitochondria", Tax_detail, ignore.case = TRUE))

df_long <- df_clean %>%
  dplyr::select(Taxonomy, Bulk, Benthic, Tax_detail) %>%
  pivot_longer(
    cols = c(Bulk, Benthic),
    names_to = "Environment",
    values_to = "Abundance"
  ) %>%
  mutate(
    Environment = case_when(
      Environment == "Bulk" ~ "Bulk water",
      Environment == "Benthic" ~ "Biofilm",
      TRUE ~ Environment
    )
  )

df_long <- df_long %>%
  group_by(Environment) %>%
  mutate(
    RelAbund = Abundance / sum(Abundance, na.rm = TRUE) * 100
  ) %>%
  ungroup()

top_n <- 20

top_taxa <- df_long %>%
  group_by(Taxonomy) %>%
  summarise(
    mean_abund = mean(RelAbund, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abund)) %>%
  slice_head(n = top_n) %>%
  pull(Taxonomy)

df_bubble <- df_long %>%
  filter(Taxonomy %in% top_taxa)

taxa_order <- df_bubble %>%
  group_by(Taxonomy) %>%
  summarise(mean_abund = mean(RelAbund, na.rm = TRUE), .groups = "drop") %>%
  arrange(mean_abund) %>%
  pull(Taxonomy)

df_bubble <- df_bubble %>%
  mutate(
    Taxonomy = factor(Taxonomy, levels = taxa_order),
    Environment = factor(Environment, levels = c("Bulk water", "Biofilm"))
  )

p_bubble <- ggplot(df_bubble, aes(x = Taxonomy, y = Environment)) +
  geom_point(
    aes(size = RelAbund, fill = Environment),
    shape = 21,
    color = "gray25",
    stroke = 0.25,
    alpha = 0.9
  ) +
  scale_size_continuous(
    range = c(1.5, 12),
    breaks = c(1, 2, 5, 10, 20),
    name = "Relative abundance (%)"
  ) +
  scale_fill_manual(
    values = c(
      "Bulk water" = "#4E79A7",
      "Biofilm" = "#E15759"
    )
  ) +
  labs(
    title = "Microbial Community Bubble Plot",
    subtitle = "Top genera in bulk water and Biofilm",
    x = "",
    y = "",
    fill = "Environment"
  ) +
  theme_genomic_plot() +
  theme(
    plot.subtitle = element_text(color = "gray35", margin = margin(b = 12)),

    axis.text.x = element_text(
      angle = 55,
      hjust = 1,
      vjust = 1,
      size = 12,
      face = "bold.italic",
      color = "gray20"
    ),

    axis.text.y = element_text(
      size = 14,
      face = "bold",
      color = "gray20"
    ),
    panel.grid.major = element_line(
      color = "gray88",
      linewidth = 0.35
    ),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      color = "gray35",
      linewidth = 0.6
    ),
    legend.position = "right"
  )

p_bubble

ggsave(
  filename = safe_out("microbial_community_bubble_plot.png"),
  plot = p_bubble,
  width = 11,
  height = 7,
  dpi = 300,
  bg = "white"
)

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

file <- data_file("picrust2.group.pathways.relative.xls")

df_raw <- read_tsv(file, show_col_types = FALSE)

colnames(df_raw)[1] <- "Pathway"

df <- df_raw %>%
  mutate(
    Pathway = as.character(Pathway),
    description = as.character(description),
    Bulk = as.numeric(Bulk),
    Benthic = as.numeric(Benthic)
  )

print(head(df))
print(colnames(df))

df_long_all <- df %>%
  dplyr::select(Pathway, description, Bulk, Benthic) %>%
  pivot_longer(
    cols = c(Bulk, Benthic),
    names_to = "Environment",
    values_to = "AbsAbund"
  ) %>%
  mutate(
    Environment = case_when(
      Environment == "Bulk" ~ "Bulk_water",
      Environment == "Benthic" ~ "Biofilm",
      TRUE ~ Environment
    )
  ) %>%
  group_by(Environment) %>%
  mutate(
    RelAbund = AbsAbund / sum(AbsAbund, na.rm = TRUE) * 100
  ) %>%
  ungroup()

check_sum <- df_long_all %>%
  group_by(Environment) %>%
  summarise(total_percent = sum(RelAbund, na.rm = TRUE), .groups = "drop")

print(check_sum)

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

df_focus <- df_long_all %>%
  filter(
    grepl(carbon_pattern, description, ignore.case = TRUE) |
      grepl(carbon_pattern, Pathway, ignore.case = TRUE)
  )

print(paste("Number of selected carbon-related pathways:", length(unique(df_focus$Pathway))))

shorten_label <- function(x, n = 58) {
  x <- gsub("&beta;", "\u03b2", x)
  x <- gsub("&alpha;", "\u03b1", x)
  x <- gsub("superpathway of ", "", x, ignore.case = TRUE)
  x <- gsub("biosynthesis", "biosyn.", x, ignore.case = TRUE)
  x <- gsub("degradation", "degrad.", x, ignore.case = TRUE)
  x <- gsub("fermentation", "ferm.", x, ignore.case = TRUE)
  x <- gsub("pathway", "path.", x, ignore.case = TRUE)
  x <- gsub("\\s+", " ", x)
  
  x <- ifelse(
    nchar(x) > n,
    paste0(substr(x, 1, n - 3), "..."),
    x
  )
  
  return(x)
}

label_df <- df_focus %>%
  distinct(Pathway, description) %>%
  mutate(Label = shorten_label(description, 58))

label_df$Label <- make.unique(label_df$Label)

df_focus <- df_focus %>%
  left_join(label_df %>% dplyr::select(Pathway, Label), by = "Pathway")
# B. Biofilm vs Bulk water log2FC enrichment bar plot

# ------------------------------------------------------------

pseudo <- 1e-6
# This pseudocount is applied to percentage abundance, as in the original code.
# The log2 ratio is descriptive; its sign does not establish significance.

df_fc <- df_focus %>%
  group_by(Pathway, Label, Environment) %>%
  summarise(
    MeanRel = mean(RelAbund, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Environment,
    values_from = MeanRel,
    values_fill = 0
  ) %>%
  mutate(
    log2FC = log2((Biofilm + pseudo) / (Bulk_water + pseudo)),
    Abs_log2FC = abs(log2FC),
    Direction = case_when(
      log2FC > 0 ~ "Biofilm enriched",
      log2FC < 0 ~ "Bulk water enriched",
      TRUE ~ "No difference"
    )
  )

print(head(df_fc))

top_n_fc <- 35

df_fc_plot <- df_fc %>%
  arrange(desc(Abs_log2FC)) %>%
  slice_head(n = top_n_fc) %>%
  arrange(log2FC)

wrap_label <- function(x, width = 38) {
  sapply(x, function(y) paste(strwrap(y, width = width), collapse = "\n"))
}

df_fc_plot$Label <- factor(df_fc_plot$Label, levels = df_fc_plot$Label)

max_abs_fc <- max(abs(df_fc_plot$log2FC), na.rm = TRUE)

p_B <- ggplot(df_fc_plot, aes(x = log2FC, y = Label, fill = Direction)) +
  geom_col(
    width = 0.72,
    color = "white",
    linewidth = 0.25
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "gray35",
    linewidth = 0.55
  ) +
  scale_fill_manual(
    values = c(
      "Biofilm enriched" = "#C8553D",
      "Bulk water enriched" = "#3D5A80",
      "No difference" = "gray70"
    )
  ) +
  scale_y_discrete(
    labels = function(x) wrap_label(x, width = 58)
  ) +
  scale_x_continuous(
    limits = c(-max_abs_fc * 1.12, max_abs_fc * 1.12),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    title = "Predicted Carbon-cycling Pathway Enrichment",
    subtitle = "Positive log2FC = enriched in Biofilm; negative log2FC = enriched in Bulk water",
    x = "log2FC  [Biofilm / Bulk water]",
    y = NULL,
    fill = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 15, base_family = "Arial") +
  theme(
    plot.title = element_text(
      size = 19,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 6)
    ),
    plot.subtitle = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5,
      color = "gray35",
      margin = margin(b = 10)
    ),
    axis.title.x = element_text(
      face = "bold",
      margin = margin(t = 8)
    ),

    axis.text.x = element_text(
      size = 12,
      face = "bold",
      color = "gray20"
    ),

    axis.text.y = element_text(
      size = 10,
      face = "bold",
      color = "gray15",
      lineheight = 0.9
    ),
    legend.position = "top",
    legend.justification = "center",
    legend.text = element_text(size = 13, face = "bold"),
    axis.line = element_line(
      color = "gray35",
      linewidth = 0.5
    ),
    axis.ticks = element_line(
      color = "gray35",
      linewidth = 0.4
    ),
    plot.margin = margin(
      t = 20,
      r = 40,
      b = 20,
      l = 20
    )
  )

p_B

ggsave(
  filename = safe_out("B_picrust2_carbon_pathway_log2FC.png"),
  plot = p_B,

  width = 13,
  height = 14,
  dpi = 300,
  bg = "white"
)
