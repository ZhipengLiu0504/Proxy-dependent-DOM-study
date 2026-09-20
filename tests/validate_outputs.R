# Run from the release root after run_all.R: Rscript tests/validate_outputs.R
library(readr)
library(dplyr)
root <- getwd()
read_stage <- function(name) readRDS(file.path(root, "outputs", name, "analysis_tables.rds"))
a <- read_stage("01_ms_composition")
b <- read_stage("03_proxy_diablo")
checks <- data.frame(check = character(), metric_value = numeric(), tolerance = numeric(), pass = logical())
add_check <- function(name, difference, tolerance = 1e-8) {
  checks <<- rbind(checks, data.frame(check = name, metric_value = difference, tolerance = tolerance,
                                   pass = is.finite(difference) && difference <= tolerance))
}
# Public checks do not require the unpublished SI reference-table files.
# They validate consistency of results regenerated from the shared inputs.
norm <- a$df_proxy_rel_sample %>% group_by(Sample) %>% summarise(total = sum(relI))
add_check("Sample pool fractions sum to one", max(abs(norm$total - 1)), 1e-12)
add_check("Exactly ten observed sample-habitat pairs", abs(nrow(a$sample_env) - 10), 0)
cross <- a$df_proxy_rel_sample %>%
  transmute(MS_sample = as.character(Sample), pool = gsub(" ", "_", as.character(pool)), rel1 = relI) %>%
  inner_join(b$df_comp %>% transmute(MS_sample, pool = as.character(pool), rel2 = relI),
             by = c("MS_sample", "pool"))
stopifnot(nrow(cross) == nrow(b$df_comp))
add_check("Pool classification agrees between MS and proxy modules", max(abs(cross$rel1-cross$rel2)), 1e-12)
mapping <- read.csv("figure_map.csv")
for (extension in c("png", "pdf")) {
  add_check(paste("All 22 mapped", extension, "panels exist"),
            sum(!file.exists(file.path("figures", paste0(mapping$figure, ".", extension)))), 0)
}
fit <- readRDS("outputs/03_proxy_diablo/proxy_diablo_model.rds")
labels <- setNames(as.character(b$matched_meta$environment), b$matched_meta$MS_sample)
score_summary <- lapply(names(fit$variates)[names(fit$variates) != "Y"], function(block) {
  scores <- fit$variates[[block]][, 1]
  habitat <- labels[rownames(fit$variates[[block]])]
  data.frame(block = block, bulk_water_mean = mean(scores[habitat == "Bulk water"]),
             biofilm_mean = mean(scores[habitat == "Biofilm"]))
}) %>% bind_rows()
write_csv(score_summary, "outputs/DIABLO_component1_orientation.csv")
add_check("Positive component 1 scores orient toward Biofilm", sum(score_summary$biofilm_mean <= score_summary$bulk_water_mean), 0)
write_csv(checks, "outputs/validation_checks.csv")
print(checks)
stopifnot(all(checks$pass))
cat("All output validation checks passed.\n")
