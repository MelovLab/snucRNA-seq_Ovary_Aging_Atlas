## Script name: Human Ovary Atlas - snucRNA-seq and GTEx deconvolved bulk compositional analysis
##
## Purpose: Run coarse cell-type-level compositional analysis on ovary snucRNA-seq and GTEx deconvolved bulk RNA-seq data
## using sccomp, generate model-based composition figures, and export supplementary result tables.
##
## Author: Josef Byrne
#########################

#### 1. Import libraries and set directories ####

# Import libraries
library(readr)
library(readxl)
library(naniar)
library(sccomp) # version 2.1.23 via Github
library(ggplot2)
library(zellkonverter)
library(SingleCellExperiment)
library(compositions)
library(Matrix)
library(dplyr)
library(tidyr)
library(stringr)
library(scales)
library(reshape2)
library(forcats)
library(purrr)
library(ggrepel)
library(patchwork)
library(tidybayes)
library(ggh4x)
library(openxlsx)

# Initialize paths
nostromo_datadir="/path/to/project_directory"
r_analyses_dir <- paste0(nostromo_datadir, "/R_analyses")
adata_objdir <- paste0(nostromo_datadir, "/python_analyses/adata_objs")
outs_dir <- paste0(r_analyses_dir, "/sccomp_results")



#### 2. Prepare snucRNA-seq and deconvolved bulk data for sccomp analysis ####

## Import AnnData object & prepare for sccomp
sce_obj <- readH5AD(paste0(adata_objdir, "/final_annotated_atlas.h5ad"), reader = "R", verbose = TRUE)

# Import snucRNA-seq pseudobulked metadata (from DGE preprocessing) that has per-sample aggregated StressScore values calculated from DGE & attach to sce_obj colData
meta_df <- readRDS(file.path(r_analyses_dir, "donor_pb_meta_df.rds"))
colData(sce_obj)$StressScore <- meta_df$StressScore[match(colData(sce_obj)$unique_sample_ID, rownames(meta_df))]

# Change colData name from "plex" to "Plex" and make it a factor
colData(sce_obj)$Plex <- factor(colData(sce_obj)$plex, levels = sort(unique(as.integer(colData(sce_obj)$plex))))
colData(sce_obj)$plex <- NULL

# Rename levels for readability
colData(sce_obj)$cell_type_fine <- forcats::fct_recode(
  colData(sce_obj)$cell_type_fine,
  "Activated Stromal" = "Activated (mito-high) Stromal",
  "IER Stromal" = "Immediate-early response Stromal",
  "ArtEnd" = "Arterial",
  "LympEnd" = "Lymphatic",
  "VenEnd" = "Venous",
  "CapEnd" = "Capillary",
  "nmSchwann" = "Nonmyelinating Schwann"
)

# Race counts - 19 White, 4 Black, 3 Hispanic, 2 Unknown
# Black has sufficient number for very large effect size inference, retain and validate if effects found
# Hispanic too rare for even large effect size inference - group with Unknown in "Other" to minimize factor bloat
# Make largest group (White) the reference level
colData(sce_obj)$Simplified_Race_Ethnicity <- fct_recode(
  colData(sce_obj)$Simplified_Race_Ethnicity,
  "Other" = "Hispanic",
  "Other" = "Unknown"
) %>%
  factor(levels = c("White", "Black", "Other"))

# Results suggest non-linear effect of Para on cell type proportions, bin to distinguish well powered groups (0, 1-3) from low powered groups (4+)
colData(sce_obj)$Para_bin <- factor(ifelse(colData(sce_obj)$Para == 0, "0", ifelse(colData(sce_obj)$Para <= 3, "1-3", "4+")), levels = c("0", "1-3", "4+"))

## Import deconvoluted bulk RNA-seq ovary data & format
bulk_prop_df <- read.csv(file.path(r_analyses_dir, "/deconv_results/deconv_results_with_metadata.csv"))

## Collapse sparse DTHHRDY level 1 into level 2 for GTEx metadata
bulk_prop_df$DTHHRDY <- ifelse(bulk_prop_df$DTHHRDY == 1, 2, bulk_prop_df$DTHHRDY)

# Change RACE to factor with levels "White" and "Black", where 2 = "Black" and 3 = "White"
bulk_prop_df$RACE <- factor(ifelse(bulk_prop_df$RACE == 2, "Black",
                             ifelse(bulk_prop_df$RACE == 3, "White", NA)),
                             levels = c("White", "Black"))


#### 3. Run sccomp analysis on snucRNA-seq data ####

## Remove rare cell types where compositional analysis is not reliable (requires >=100 cells and detection in >=10 donors)
# Compute counts per cell type and donors per cell type
ct_stats <- as.data.frame(colData(sce_obj)) %>%
  group_by(cell_type_fine) %>%
  dplyr::summarize(
    n_cells = n(),
    n_donors = n_distinct(SenNet_ID),
    .groups = "drop"
  )

min_donors <- 10
min_cells <- 100

# Determine eligibility for filtering
ct_stats <- ct_stats %>%
  mutate(eligible = n_cells >= min_cells & n_donors >= min_donors)

# Filter cell types: >= min_cells cells and in >= min_donors donors
eligible_cell_types <- ct_stats %>%
  filter(eligible) %>%
  pull(cell_type_fine)

# Copy sce_obj to sce_obj_sccomp for additional handling
sce_obj_sccomp <- sce_obj[, colData(sce_obj)$cell_type_fine %in% eligible_cell_types]

# Print counts for each cell_type_fine
ct_counts <- as.data.frame(colData(sce_obj_sccomp)) %>%
  group_by(cell_type_fine) %>%
  summarize(n_cells = n(), .groups = "drop") %>%
  arrange(desc(n_cells))
print(ct_counts, n = nrow(ct_counts))

# Print cell types removed, with n_cells and n_donors
removed_cts <- ct_stats %>% dplyr::filter(!eligible)
cat("Cell types removed due to too few cells or donors:\n")
removed_strings <- paste0(
  removed_cts$cell_type_fine,
  " (", removed_cts$n_cells, " cells, ", removed_cts$n_donors, " donors)"
)
result_string <- if (length(removed_strings) > 0) paste(removed_strings, collapse = ", ") else "None"
cat(result_string, "\n")
# Cell types removed due to too few cells or donors:
# CD16+ Monocyte (32 cells, 14 donors), Cyc Immune (29 cells, 12 donors), Lymphatic Valve Up (46 cells, 14 donors), Mast (19 cells, 5 donors),
# B Cell (11 cells, 10 donors), Myelinating Schwann (11 cells, 8 donors), NK (71 cells, 6 donors), Stress-response Secretory Epithelial (48 cells, 9 donors),
# cDC1 (45 cells, 19 donors)

# Ensure mathematical separability via VIF
cell_counts <- as.data.frame(table(
  sample = colData(sce_obj_sccomp)$unique_sample_ID,
  cell_group = colData(sce_obj_sccomp)$cell_type_fine
))
cell_totals <- as.data.frame(table(
  sample = colData(sce_obj_sccomp)$unique_sample_ID
))
colnames(cell_totals)[2] <- "n_cells_total"

cell_counts <- cell_counts %>%
  left_join(cell_totals, by = "sample") %>%
  mutate(obs_prop = Freq / n_cells_total)

# Get per-sample metadata needed for formula variables
meta_vars <- as.data.frame(colData(sce_obj_sccomp)) %>%
  dplyr::select(unique_sample_ID, Age, Para_bin, Para, BMI, Simplified_Race_Ethnicity, StressScore) %>% distinct()

# Merge observed proportions with sample metadata, expand to all cell_groups for each sample
obs_prop_df_snuc <- cell_counts %>% left_join(meta_vars, by = c("sample" = "unique_sample_ID"))

car::vif(lm(obs_prop ~ Age + Para_bin + BMI + Simplified_Race_Ethnicity, data = obs_prop_df_snuc))
#                               GVIF Df GVIF^(1/(2*Df))
# Age                       1.435726  1        1.198218
# Para_bin                  1.654405  2        1.134124
# BMI                       1.180125  1        1.086336
# Simplified_Race_Ethnicity 1.292532  2        1.066253

# Change "cell_type_fine" column/colData to "cell_group" to align with sccomp conventions
colData(sce_obj_sccomp)$cell_group <- colData(sce_obj_sccomp)$cell_type_fine
colData(sce_obj_sccomp)$cell_group <- droplevels(colData(sce_obj_sccomp)$cell_group)

form_composition_snuc <- ~ Age + Para_bin + BMI + Plex + Simplified_Race_Ethnicity + StressScore + (1 | SenNet_ID)
form_variability_snuc <- ~ 1

## Run sccomp compositional analysis on snuc
set.seed(42)
sccomp_res_snuc <- sce_obj_sccomp %>%
  sccomp_estimate(
    formula_composition = form_composition_snuc,
    formula_variability = form_variability_snuc,
    sample = "unique_sample_ID",
    cell_group = "cell_group",
    cores = 1, # run single-threaded for deterministic results
    mcmc_seed = 42,
    bimodal_mean_variability_association = TRUE,
    inference_method = "pathfinder", # default for sccomp_estimate
    cleanup_draw_files = FALSE,
    output_directory = paste0(r_analyses_dir, "/R_sccomp_cache")
  )

set.seed(42)
sccomp_res_snuc <- sccomp_remove_outliers(sccomp_res_snuc, mcmc_seed = 42, cores = 1, cleanup_draw_files = FALSE,
  output_directory = paste0(r_analyses_dir, "/R_sccomp_cache"),
  cache_stan_model = paste0(r_analyses_dir, "/R_sccomp_cache")) %>%
  sccomp_test()

# Extract final sccomp outliers (censored sample×cell_group pairs) with counts + metadata
extract_sccomp_final_outliers <- function(sccomp_res, meta = NULL, sample_col = "sample") {
  mod_input <- attr(sccomp_res, "model_input")
  if (is.null(mod_input)) stop("No `model_input` attribute found on sccomp_res.")

  # Find final outliers: truncation_down == -1 indicates censored cell
  out_idx <- which(mod_input$truncation_down == -1L, arr.ind = TRUE)

  # Handle case of no outliers
  if (nrow(out_idx) == 0) {
    out_tbl <- tibble::tibble(
      sample = character(0),
      cell_group = character(0),
      cell_count = integer(0)
    )
    return(out_tbl)
  }

  samples <- rownames(mod_input$X)
  cell_groups <- colnames(mod_input$y)
  if (is.null(samples) || is.null(cell_groups)) {
    stop("Could not find rownames(model_input$X) and/or colnames(model_input$y).")
  }

  # Build base outlier table
  out_tbl <- tibble::tibble(
    sample = samples[out_idx[, 1]],
    cell_group = cell_groups[out_idx[, 2]],
    cell_count = as.integer(mod_input$y[cbind(out_idx[, 1], out_idx[, 2])])
  ) %>%
    dplyr::arrange(sample, cell_group)

  # Attach metadata
  if (!is.null(meta)) {
    # user-supplied metadata
    if (!sample_col %in% colnames(meta)) {
      stop(sprintf("`meta` does not contain sample_col = '%s'.", sample_col))
    }
    out_tbl <- out_tbl %>%
      dplyr::left_join(dplyr::distinct(meta), by = c("sample" = sample_col))
  } else {
    # derive metadata from design matrix X (model matrix)
    X_df <- as.data.frame(mod_input$X)
    X_df$sample <- rownames(mod_input$X)

    out_tbl <- out_tbl %>%
      dplyr::left_join(X_df, by = "sample")
  }
  out_tbl
}

# Extract variable names from formula
form_vars <- setdiff(attr(terms(update(form_composition_snuc, . ~ . - (1 | SenNet_ID))), "term.labels"), ".")
form_vars <- unique(c(form_vars, "SenNet_ID"))

# Pull from colData(sce_obj_sccomp) including unique_sample_id which is present there
meta_df_used_snuc <- as.data.frame(SummarizedExperiment::colData(sce_obj_sccomp)) %>%
  dplyr::select(dplyr::all_of(c(form_vars, "unique_sample_ID"))) %>%
  distinct()

outliers_tbl_snuc <- extract_sccomp_final_outliers(sccomp_res_snuc, meta = meta_df_used_snuc, sample_col = "unique_sample_ID")
outliers_tbl_snuc

# Save sccomp results object & df and sccomp-censored outliers table to CSV
saveRDS(sccomp_res_snuc, file = paste0(outs_dir, "/sccomp_res_snuc.rds"))

write.csv(as.data.frame(sccomp_res_snuc), file = paste0(outs_dir, "/sccomp_res_snuc_df.csv"), row.names = FALSE)
write.csv(outliers_tbl_snuc, file = paste0(outs_dir, "/sccomp_snuc_outliers.csv"), row.names = FALSE)



#### 4. Run sccomp analysis on deconvoluted bulk RNA-seq data ####

# Coerce columns to factor or numeric as needed
factor_vars <- intersect(c("RACE", "DTHHRDY", "SMCENTER"), colnames(bulk_prop_df))
bulk_prop_df[factor_vars] <- lapply(bulk_prop_df[factor_vars], as.factor)

numeric_vars <- intersect(c("AGE", "BMI", "SMRIN", "SMTSISCH", "prop"), colnames(bulk_prop_df))
bulk_prop_df[numeric_vars] <- lapply(bulk_prop_df[numeric_vars], as.numeric)

# Change "Celltype" column to "cell_group" for consistency with downstream analysis
bulk_prop_df$cell_group <- bulk_prop_df$Celltype

# Find indexes (row numbers) of rows with NA in any of these numeric variables
# Check if there are any NA rows for any variables included in the formula
# Get the unique SUBJID that have an NA in any formula variable
formula_vars <- c("AGE", "BMI", "RACE", "DTHHRDY", "SMCENTER", "SMTSISCH", "SMRIN", "prop")
vars_in_df <- formula_vars[formula_vars %in% colnames(bulk_prop_df)]
na_rows_formula <- which(rowSums(is.na(bulk_prop_df[, vars_in_df, drop = FALSE])) > 0)
if (length(na_rows_formula) > 0) {
  # Get SUBJID for those rows
  subjids_with_na <- unique(bulk_prop_df$sample_name[na_rows_formula])
  cat("Unique SUBJID values with any NA in formula variables:\n")
  print(subjids_with_na)
} else {
  cat("No rows with NA in the formula variables.\n")
}

# Remove 2 donors with NA for Hardy score, this throws errors in sccomp_estimate
bulk_prop_df <- bulk_prop_df[!bulk_prop_df$sample_name %in% subjids_with_na, ]
bulk_prop_df <- droplevels(bulk_prop_df)
cat("Number of unique sample names after dropping NA donors:", length(unique(bulk_prop_df$sample_name)), "\n")
# Number of unique sample names after dropping NA donors: 104

# Multiply prop by 100k pseudo-count to allow sccomp modeling to include relatively rare cell types
# Note: Pseudocount multiple marginally impacts results for well represented cell types while avoiding clipping of rare cell types to zero.
bulk_100k_df <- bulk_prop_df %>% mutate(prop = round(prop * 100000))
bulk_100k_df$prop <- as.integer(bulk_100k_df$prop)

zero_prop_stats_100k <- bulk_100k_df %>%
  dplyr::group_by(cell_group) %>%
  dplyr::summarise(
    n = dplyr::n(),
    n_zero = sum(prop == 0, na.rm = TRUE),
    pct_zero = 100 * mean(prop == 0, na.rm = TRUE)
  ) %>%
  dplyr::arrange(desc(pct_zero))

print(zero_prop_stats_100k, n = 32)

# Remove cell types (cell_group) that have zeros in >50% samples in bulk_100k_df
celltypes_keep <- zero_prop_stats_100k %>%
  dplyr::filter(pct_zero <= 50) %>%
  dplyr::pull(cell_group)

bulk_100k_df <- bulk_100k_df %>%
  dplyr::filter(cell_group %in% celltypes_keep)

# Identify the cell types removed because they had zeros in >50% samples
celltypes_removed <- zero_prop_stats_100k %>%
  dplyr::filter(!(cell_group %in% celltypes_keep))

cat("Cell types removed due to >50% zero values:\n")
print(celltypes_removed[, c("cell_group", "pct_zero")], row.names = FALSE)
# Cell types removed due to >50% zero values:
# # A tibble: 1 × 2
#   cell_group pct_zero
#   <chr>         <dbl>
# 1 CD16- NK       51.0

# Ensure mathematical separability via VIF
obs_prop_df_bulk <- bulk_100k_df %>%
  group_by(sample_name) %>%
  mutate(n_cells_total = sum(prop),
         obs_prop = prop / n_cells_total) %>%
  ungroup()

# VIF (note: obs_prop is per sample×celltype; this is fine for a quick collinearity check)
car::vif(lm(obs_prop ~ AGE + BMI + RACE + DTHHRDY + SMCENTER + SMTSISCH + SMRIN,
            data = obs_prop_df_bulk))
#              GVIF Df GVIF^(1/(2*Df))
# AGE      1.199843  1        1.095373
# BMI      1.082820  1        1.040586
# RACE     1.089110  1        1.043604
# DTHHRDY  3.116857  3        1.208610
# SMCENTER 1.353442  1        1.163375
# SMTSISCH 2.642346  1        1.625529
# SMRIN    1.811254  1        1.345828

# Expected modest collinearity between technical variables (SMTSISCH, SMRIN)
# No concerning collinearity between variables of interest (AGE, RACE) and other variables

## Run sccomp compositional analysis on deconvolved bulk RNA-seq data with pseudocounts
# Sufficiently large sample size to include Age, BMI & QC/technical variables in the variability model without instability issues
form_composition_bulk <- ~ AGE + BMI + RACE + DTHHRDY + SMCENTER + SMTSISCH + SMRIN
form_variability_bulk <- ~ AGE + BMI + DTHHRDY + SMCENTER + SMTSISCH + SMRIN

set.seed(42)
sccomp_res_bulk <- sccomp::sccomp_estimate(
  bulk_100k_df,
  formula_composition = form_composition_bulk,
  formula_variability = form_variability_bulk,
  sample = "sample_name",
  cell_group = "cell_group",
  cores = 1, # run single-threaded for deterministic results
  mcmc_seed = 42,
  abundance = "prop",
  inference_method = "pathfinder", # default for sccomp_estimate
  cleanup_draw_files = FALSE,
  output_directory = paste0(r_analyses_dir, "/R_sccomp_cache")
)

set.seed(42)
sccomp_res_bulk <- sccomp_res_bulk %>%
  sccomp_remove_outliers(mcmc_seed = 42, cores = 1, cleanup_draw_files = FALSE, output_directory = paste0(r_analyses_dir, "/R_sccomp_cache")) %>%
  sccomp_test()

## Extract outliers detected and censored by sccomp analysis
form_vars <- setdiff(attr(terms(form_composition_bulk), "term.labels"), ".")
form_vars <- unique(c(form_vars, "sample_name"))

meta_df_used_bulk <- bulk_prop_df %>%
  dplyr::select(dplyr::all_of(form_vars)) %>%
  distinct()

outliers_tbl_bulk <- extract_sccomp_final_outliers(sccomp_res_bulk, meta = meta_df_used_bulk, sample_col = "sample_name")

# Save sccomp results object & df and sccomp-censored outliers table to CSV
saveRDS(sccomp_res_bulk, file = paste0(outs_dir, "/sccomp_res_bulk.rds"))

write.csv(as.data.frame(sccomp_res_bulk), file = paste0(outs_dir, "/sccomp_res_bulk_df.csv"), row.names = FALSE)
write.csv(outliers_tbl_bulk, file = paste0(outs_dir, "/sccomp_bulk_outliers.csv"), row.names = FALSE)


#### 5. Plot sccomp results plots and tables ####

# Function to plot sccomp results as a forest plot (individual)
plot_sccomp_forest_one <- function(
    df,
    parameter_keep = "Age",
    which = c("c", "v"), # "c" for composition, "v" for variability
    signif_stat = c("FDR", "pH0"), # which stat drives significance color
    signif_threshold = 0.05,
    effect_threshold = 0.1, # logit fold-change threshold (via Region of Practical Equivalence - ROPE)
    order_by = c("effect", "signif_then_effect", "custom"),
    custom_order = NULL, # character vector of cell_group top->bottom
    wrap_width = 28,
    base_size = 11,
    point_size = 2,
    xlab = "Composition effect\n(log-odds; 95% CrI)",
    title = NULL,
    x_breaks = waiver(), # Argument for custom x breaks/grid intervals, default waiver()
    digits = 2 # Number of digits after 0 to include
) {
  which <- match.arg(which)
  signif_stat <- match.arg(signif_stat)
  order_by <- match.arg(order_by)

  # pick columns based on which = "c" or "v"
  prefix <- paste0(which, "_")
  lower_col <- paste0(prefix, "lower")
  effect_col <- paste0(prefix, "effect")
  upper_col <- paste0(prefix, "upper")
  stat_col <- paste0(prefix, signif_stat)
  cell_group <- attr(df, ".cell_group")

  required <- c(cell_group, "parameter", lower_col, effect_col, upper_col, stat_col)
  if (!all(required %in% colnames(df))) {
    missing <- setdiff(required, colnames(df))
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }

  plot_df <- df %>%
    filter(parameter == parameter_keep) %>%
    transmute(
      !!cell_group,
      lower = .data[[lower_col]],
      effect = .data[[effect_col]],
      upper = .data[[upper_col]],
      stat = .data[[stat_col]],
      is_signif = !is.na(stat) & stat < signif_threshold,
      cell_group_wrapped = str_wrap(!!cell_group, width = wrap_width)
    ) %>%
    filter(!is.na(effect), !is.na(lower), !is.na(upper))

  # Order by effect size
  plot_df <- plot_df %>% mutate(cell_group_wrapped = fct_reorder(cell_group_wrapped, effect))

  # labels used both for plotting and for scale names
  lab_not <- paste0(signif_stat, " \u2265 ", signif_threshold) # ≥
  lab_up <- paste0(signif_stat, " < ", signif_threshold, " & up")
  lab_down <- paste0(signif_stat, " < ", signif_threshold, " & down")

  plot_df <- plot_df %>%
    mutate(color_group = dplyr::case_when(
      is_signif & effect >= 0 ~ lab_up,
      is_signif & effect < 0 ~ lab_down,
      TRUE ~ lab_not
    ))

  color_values <- c(
    setNames("grey55", lab_not),
    setNames("#b2182b", lab_up), # dark red for significant up
    setNames("#08519c", lab_down) # dark blue for significant down
  )

  # Dynamic accuracy for axis labels based on digits
  accuracy_val <- 10^(-digits)
  # Use element_text for axis.text.y since we don't do markdown/bold anymore
  ggplot(plot_df, aes(x = effect, y = cell_group_wrapped)) +
    geom_vline(xintercept = 0, linewidth = 0.2, color = "grey70") +
    # geom_vline(
    #   xintercept = c(-effect_threshold, effect_threshold), # Plot ROPE as dashed line
    #   linewidth = 0.2, color = "grey80", linetype = "dashed"
    # ) +
    # Add error bars with T-shaped ends (short "|" at each end)
    geom_errorbarh(aes(xmin = lower, xmax = upper, color = color_group), height = 0.5, linewidth = 0.35) +
    geom_point(aes(color = color_group), size = point_size) +
    scale_color_manual(values = color_values, name = NULL) +
    scale_x_continuous(labels = label_number(accuracy = accuracy_val), breaks = x_breaks) +  # user-controlled breaks and digits
    labs(
      title = title,
      x = xlab,
      y = NULL
    ) +
    theme_minimal(base_family = "DejaVu Sans", base_size = base_size) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(linewidth = 0.2), # Make major gridlines thinner
      axis.text.y = element_text(size = base_size, color = "black", margin = margin(r = 4)), # Padding between y ticks and plot, font black
      axis.text.x = element_text(size = base_size, color = "black", margin = margin(t = 2)), # Padding between x ticks and plot, font black
      axis.title.x = element_text(size = base_size + 2, color = "black", margin = margin(t = 4)), # Padding between x-axis label and plot, font black
      plot.title = element_text(size = base_size + 2, hjust = 0.5, color = "black", margin = margin(b = 4)), # font black
      legend.position = "none",
      legend.text = element_text(size = base_size, color = "black"), # font black
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2) # Adds box around main plot
    )
}

# Function to plot sccomp results as a forest plot in a faceted panel for easier comparison
plot_sccomp_forest_comparison <- function(
    df, # Now expects the combined dataframe
    parameter_keep = "Age",
    which = "c",
    signif_stat = "FDR",
    signif_threshold = 0.05,
    effect_threshold = 0.1,
    base_size = 11,
    point_size = 2,
    x_breaks = waiver(),
    digits = 2
) {


  prefix <- paste0(which, "_")
  lower_col <- paste0(prefix, "lower")
  effect_col <- paste0(prefix, "effect")
  upper_col <- paste0(prefix, "upper")
  stat_col <- paste0(prefix, signif_stat)
  cell_group <- attr(df, ".cell_group")

  plot_df <- df %>%
    filter(parameter == parameter_keep) %>%
    transmute(
      method,
      !!cell_group,
      lower = .data[[lower_col]],
      effect = .data[[effect_col]],
      upper = .data[[upper_col]],
      stat = .data[[stat_col]],
      is_signif = !is.na(stat) & stat < signif_threshold
    ) %>%
    filter(!is.na(effect))

  # Create a shared order based on effect
  order_ref <- plot_df %>%
    group_by(!!sym(cell_group)) %>%
    summarise(
      ordering_effect = dplyr::case_when(
        any(method == "snucRNA-seq") ~ effect[method == "snucRNA-seq"][1], # Use snucRNA-seq's effect if present
        TRUE ~ effect[1] # Else use whatever effect there is
      ),
      .groups = "drop"
    ) %>%
    arrange(ordering_effect)

  # Extract the levels as a character vector & apply
  ordered_levels <- order_ref %>% pull(cell_group)
  plot_df <- plot_df %>%
    mutate(cell_group = factor(cell_group, levels = ordered_levels))

  # Define colors
  lab_not <- "Non-significant"; lab_up <- "Significant Up"; lab_down <- "Significant Down"

  plot_df <- plot_df %>%
    mutate(color_group = case_when(
      is_signif & effect >= 0 ~ lab_up,
      is_signif & effect < 0 ~ lab_down,
      TRUE ~ lab_not
    ))

  color_values <- c("Non-significant" = "grey70",  "Significant Up" = "#b2182b", "Significant Down" = "#08519c")

  ggplot(plot_df, aes(x = effect, y = !!sym(cell_group))) +
    facet_wrap2(~method, axes = "y", remove_labels = "none", scales = "fixed") + # Keeps the x-axis scale the same for comparison
    geom_vline(xintercept = 0, linewidth = 0.2, color = "grey60") +
    geom_errorbarh(aes(xmin = lower, xmax = upper, color = color_group), height = 0.3, linewidth = 0.5) +
    geom_point(aes(color = color_group), size = point_size) +
    scale_color_manual(values = color_values, name = NULL) +
    scale_x_continuous(labels = scales::label_number(accuracy = 10^-digits), breaks = x_breaks) +
    labs(x = "Composition effect (log-odds; 95% CrI)", y = NULL) +
    theme_minimal(base_family = "DejaVu Sans", base_size = base_size) +
    theme(
      panel.grid.major.y = element_line(linewidth = 0.2, color = "gray90"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.y = element_text(size = base_size, color = "black", margin = margin(r = 4)), # Padding between y ticks and plot, font black
      axis.text.x = element_text(size = base_size, color = "black", margin = margin(t = 2)), # Padding between x ticks and plot, font black
      axis.title.x = element_text(size = base_size + 2, color = "black", margin = margin(t = 8)), # Padding between x-axis label and plot, font black
      plot.title = element_text(size = base_size + 2, hjust = 0.5, color = "black"),
      panel.spacing = unit(1, "lines"), # Adds space between the two plots
      strip.text = element_text(size = base_size + 2, margin = margin(b = 6)), # method labels
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2) # Adds box around main plot
    )

}

logit_fold_effect_threshold <- attr(sccomp_res_snuc, "test_composition_above_logit_fold_change") # aka delta, same for snuc and bulk

sccomp_res_snuc <- sccomp_res_snuc %>%
  mutate(parameter = recode(parameter, "Simplified_Race_EthnicityBlack" = "RaceBlack"), method = "snucRNA-seq")

sccomp_res_bulk <- sccomp_res_bulk %>%
  mutate(parameter = recode(parameter, "AGE" = "Age", "RACEBlack" = "RaceBlack"), method = "Deconvoluted Bulk (GTEx)")

# For deconvoluted bulk, remove cell type annotations that are ambiguous, sinks/combined types, or...
# cell types used to soak up variability where biological interpretation is unclear
bulk_types_exclude_plot <- c("Stress-signaling Mural", "Granulocyte-myeloid Sink", "Lymphatic Valve (ambiguous)",
  "Schwann/Nerve", "Erythrocyte", "Cyc Immune")

combined_res <- bind_rows(sccomp_res_snuc, sccomp_res_bulk) %>%
  filter(!cell_group %in% bulk_types_exclude_plot)
combined_res$method <- factor(combined_res$method, levels = c("snucRNA-seq", "Deconvoluted Bulk (GTEx)"))

# Make faceted plot for comparison of snuc and bulk sccomp Age results
p1 <- plot_sccomp_forest_comparison(combined_res, parameter_keep = "Age", which = "c", signif_stat = "FDR",
  signif_threshold = 0.05, effect_threshold = logit_fold_effect_threshold, base_size = 5, point_size = 0.8, x_breaks = waiver(), digits = 1)

ggsave(filename = paste0(outs_dir, "/sccomp_forest_snuc-bulk_comparison_age.pdf"), plot = p1, width = 180,
  height = 90, units = "mm", device = cairo_pdf, dpi = 500
)

# Make individual plot for snuc Parity results (no reproductive history metadata in GTEx to compare)
p2 <- plot_sccomp_forest_one(
  sccomp_res_snuc,
  parameter_keep = "Para_bin1-3",
  which = "c",
  signif_stat = "FDR",
  signif_threshold = 0.05,
  effect_threshold = logit_fold_effect_threshold,
  order_by = "effect",
  base_size = 5,
  title = "Parity (0 vs 1-3)",
  point_size = 0.8,
  x_breaks = seq(-2, 1, 1),
  digits = 1
)

ggsave(filename = paste0(outs_dir, "/sccomp_forest_snuc_parity_fin.pdf"), plot = p2, width = 55,
  height = 90, units = "mm", device = cairo_pdf, dpi = 500
)

# Make faceted plot for comparison of snuc and bulk sccomp Race (Black vs White) results
p3 <- plot_sccomp_forest_comparison(combined_res, parameter_keep = "RaceBlack", which = "c", signif_stat = "FDR",
  signif_threshold = 0.05, effect_threshold = logit_fold_effect_threshold, base_size = 5, point_size = 0.8, x_breaks = waiver(), digits = 1)

ggsave(filename = paste0(outs_dir, "/sccomp_forest_snuc-bulk_comparison_race.pdf"), plot = p3, width = 180,
  height = 90, units = "mm", device = cairo_pdf, dpi = 500
)



#### 6. Extract model predictions for snuc- and bulk- cell type changes for Age, Para_bin, and Race ####

# Function to extract predictions for a given focal variable (e.g. Age, Para_bin, Race) and marginalize over nuisance variables (e.g. donor+plex, donor+center).
# Works for continuous (e.g., Age grid) and categorical (e.g., Para_bin levels) focal variables.
sccomp_predict_marginal <- function(
  fit_sccomp,
  focal_var,                    # string, e.g. "Age" or "Para_bin"
  focal_values = NULL,          # numeric vector or character/factor levels; if NULL and numeric -> seq over range
  n_grid = 100,                 # only used when focal_values is NULL and focal_var is numeric
  marginal_vars = character(0), # strings, e.g. c("SenNet_ID","Plex") or c("DTHHRDY","SMCENTER")
  ref_levels = list(),          # named list of fixed reference levels for factors (and optionally numerics)
  mcmc_seed = 42,
  number_of_draws = 500,
  robust = TRUE
) {
  cd0 <- attr(fit_sccomp, "count_data")
  sample_col <- as.character(attr(fit_sccomp, ".sample"))
  cell_col   <- as.character(attr(fit_sccomp, ".cell_group"))
  count_col  <- as.character(attr(fit_sccomp, ".count"))
  form_comp  <- attr(fit_sccomp, "formula_composition")

  keep_cols <- setdiff(names(cd0), c(cell_col, count_col))
  covars <- cd0 %>%
    dplyr::select(dplyr::all_of(keep_cols)) %>%
    dplyr::distinct(.data[[sample_col]], .keep_all = TRUE)

  # Baseline/reference row: numeric -> median (unless provided in ref_levels), factor -> ref_levels if provided
  ref_covars <- covars[1, , drop = FALSE]
  for (nm in names(ref_covars)) {
    if (is.numeric(ref_covars[[nm]])) {
      ref_covars[[nm]] <- if (!is.null(ref_levels[[nm]])) ref_levels[[nm]] else median(covars[[nm]], na.rm = TRUE)
    } else if (is.factor(ref_covars[[nm]])) {
      if (!is.null(ref_levels[[nm]])) ref_covars[[nm]] <- factor(ref_levels[[nm]], levels = levels(covars[[nm]]))
    }
  }

  # Determine focal values
  if (is.null(focal_values)) {
    if (!is.numeric(covars[[focal_var]])) stop("focal_values must be supplied for non-numeric focal_var.")
    focal_values <- seq(min(covars[[focal_var]], na.rm = TRUE),
                        max(covars[[focal_var]], na.rm = TRUE),
                        length.out = n_grid)
  } else {
    # preserve factor levels if needed
    if (is.factor(covars[[focal_var]])) {
      focal_values <- factor(focal_values, levels = levels(covars[[focal_var]]))
    }
  }

  # Levels for marginal vars
  marg_levels <- lapply(marginal_vars, function(v) {
    x <- cd0[[v]]
    if (is.factor(x)) levels(x) else sort(unique(x))
  })
  names(marg_levels) <- marginal_vars

  # Build crossing grid: focal x marginal vars
  grid_args <- c(setNames(list(focal_values), focal_var), marg_levels)
  new_data <- do.call(tidyr::crossing, grid_args)

  # Make factor types consistent with training data for focal + marginal vars
  for (v in c(focal_var, marginal_vars)) {
    if (is.factor(covars[[v]])) new_data[[v]] <- factor(new_data[[v]], levels = levels(covars[[v]]))
  }

  # Fill remaining covariates with baseline refs
  fill_cols <- setdiff(names(ref_covars), c(sample_col, focal_var, marginal_vars))
  for (nm in fill_cols) new_data[[nm]] <- ref_covars[[nm]][[1]]

  # Dummy sample IDs
  new_data[[sample_col]] <- paste0("grid_", seq_len(nrow(new_data)))

  set.seed(mcmc_seed)
  pred <- sccomp_predict(
    fit_sccomp,
    formula_composition = form_comp,
    new_data = new_data,
    mcmc_seed = mcmc_seed,
    number_of_draws = number_of_draws,
    robust = robust
  )

  list(pred = pred, new_data = new_data, covars = covars)
}

## snuc Age, Para_bin (0 vs 1-3), and Race (White vs Black) predictions
# snuc Age curve (marginalize over donor+plex, fix race and Para_bin to most populous bins)
out_snuc_age <- sccomp_predict_marginal(
  fit_sccomp = sccomp_res_snuc,
  focal_var = "Age",
  focal_values = NULL, n_grid = 100,
  marginal_vars = c("SenNet_ID","Plex"),
  ref_levels = list(Simplified_Race_Ethnicity="White", Para_bin="1-3"),
  number_of_draws = 500, robust = TRUE
)
saveRDS(out_snuc_age, file = paste0(outs_dir, "/pred_snuc_age_sccomp_marginal.rds"))

# snuc Para categorical contrast (Para_bin = c("0","1-3") at fixed Age, marginalize donor+plex)
out_snuc_para <- sccomp_predict_marginal(
  fit_sccomp = sccomp_res_snuc,
  focal_var = "Para_bin",
  focal_values = c("0", "1-3"),
  marginal_vars = c("SenNet_ID", "Plex"),
  ref_levels = list(Simplified_Race_Ethnicity = "White", Age = 55),
  number_of_draws = 500, robust = TRUE
)
saveRDS(out_snuc_para, file = paste0(outs_dir, "/pred_snuc_para_sccomp_marginal.rds"))

# snuc Race categorical contrast (Race = c("White","Black") at fixed Age, marginalize donor+plex)
out_snuc_race <- sccomp_predict_marginal(
  fit_sccomp = sccomp_res_snuc,
  focal_var = "Simplified_Race_Ethnicity",
  focal_values = c("White","Black"),
  marginal_vars = c("SenNet_ID","Plex"),
  ref_levels = list(Para_bin="1-3", Age=55),
  number_of_draws = 500, robust = TRUE
)
saveRDS(out_snuc_race, file = paste0(outs_dir, "/pred_snuc_race_sccomp_marginal.rds"))

## bulk Age and Race (White vs Black) predictions
# bulk Age curve (marginalize over donor+center, fix race to most populous group)
out_bulk_age <- sccomp_predict_marginal(
  fit_sccomp = sccomp_res_bulk,
  focal_var = "AGE",
  focal_values = NULL, n_grid = 100,
  marginal_vars = c("DTHHRDY","SMCENTER"),
  ref_levels = list(RACE="White"),
  number_of_draws = 500, robust = TRUE
)
saveRDS(out_bulk_age, file = paste0(outs_dir, "/pred_bulk_age_sccomp_marginal.rds"))

# bulk Race categorical contrast (Race = c("White","Black") at fixed Age, marginalize donor+center)
out_bulk_race <- sccomp_predict_marginal(
  fit_sccomp = sccomp_res_bulk,
  focal_var = "RACE",
  focal_values = c("White","Black"),
  marginal_vars = c("DTHHRDY","SMCENTER"),
  ref_levels = list(AGE=55),
  number_of_draws = 500, robust = TRUE
)
saveRDS(out_bulk_race, file = paste0(outs_dir, "/pred_bulk_race_sccomp_marginal.rds"))


#### 7. Plot marginalized model prediction results for snuc- and bulk- cell type changes for Age, Para_bin, and Race ####

### Categorical: plot observed proportions + sccomp predictions for a given cell type and categorical variable (e.g. Para_bin, Race)
# Function to plot observed proportions + sccomp predictions for a given cell type and categorical variable
plot_cat_obs_plus_sccomp <- function(
  fit_sccomp,
  out_cat,
  celltype,
  cat_var = "Para_bin",
  levels_to_plot = NULL,
  age_var = "Age",
  ref_age = 55,
  y_labels = scales::label_percent(accuracy = 0.1),
  y_breaks = NULL,
  point_alpha = 0.25,
  outliers_tbl = NULL,
  base_size = 5,
  plot_title = NULL,
  plot_x_label = NULL
) {
  cd0 <- attr(fit_sccomp, "count_data")
  sample_col <- as.character(attr(fit_sccomp, ".sample"))
  cell_col   <- as.character(attr(fit_sccomp, ".cell_group"))
  count_col  <- as.character(attr(fit_sccomp, ".count"))

  # observed sample proportions for this cell type
  obs_df <- cd0 %>%
    dplyr::group_by(.data[[sample_col]]) %>%
    dplyr::mutate(obs_total = sum(.data[[count_col]], na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(obs_prop = .data[[count_col]] / obs_total) %>%
    dplyr::filter(.data[[cell_col]] == celltype) %>%
    dplyr::select(
      !!sample_col, !!cat_var, !!age_var, obs_prop
    )

  # Filter for chosen cat levels (if given)
  if (!is.null(levels_to_plot)) {
    obs_df <- obs_df %>% dplyr::filter(.data[[cat_var]] %in% levels_to_plot)
  }

  # If outliers_tbl is provided, annotate outliers
  obs_df$outlier <- FALSE
  obs_df$shape <- 16   # default: circle
  obs_df$color <- NA

  if (!is.null(outliers_tbl) &&
      all(c("sample", "cell_group") %in% names(outliers_tbl))) {

    # samples with the celltype that are outliers
    is_outlier <- paste0(obs_df[[sample_col]], "___", celltype) %in%
      paste0(outliers_tbl$sample, "___", outliers_tbl$cell_group)

    obs_df$outlier[is_outlier] <- TRUE
    obs_df$shape[is_outlier] <- 15   # square
    obs_df$color[is_outlier] <- "red"
    obs_df$color[!is_outlier] <- "blue"
  } else {
    obs_df$color <- "blue"
    obs_df$shape <- 16
    obs_df$outlier <- FALSE
  }

  # (optional) restrict observed points to same age reference window if you want:
  # obs_df <- obs_df %>% dplyr::filter(abs(.data[[age_var]] - ref_age) <= 5)

  # model-predicted group means ± CrI from out_cat$pred
  pred <- out_cat$pred

  pred_sum <- pred %>%
    dplyr::filter(.data[[cell_col]] == celltype) %>%
    {
      # If levels_to_plot is given, filter here too
      if (!is.null(levels_to_plot)) dplyr::filter(., .data[[cat_var]] %in% levels_to_plot) else .
    } %>%
    dplyr::group_by(.data[[cat_var]]) %>%
    dplyr::summarise(
      proportion_mean  = mean(.data$proportion_mean),
      proportion_lower = mean(.data$proportion_lower),
      proportion_upper = mean(.data$proportion_upper),
      .groups = "drop"
    )

  # Ensure cat ordering matches requested order (if factor)
  if (!is.null(levels_to_plot)) {
    pred_sum[[cat_var]] <- factor(pred_sum[[cat_var]], levels = levels_to_plot)
    if (is.factor(obs_df[[cat_var]])) {
      obs_df[[cat_var]] <- factor(obs_df[[cat_var]], levels = levels_to_plot)
    }
  } else if (is.factor(obs_df[[cat_var]])) {
    pred_sum[[cat_var]] <- factor(pred_sum[[cat_var]], levels = levels(obs_df[[cat_var]]))
  }

  # plot
  p <- ggplot() +
    geom_point(
      data = obs_df,
      aes(
        x = .data[[cat_var]], y = obs_prop,
        shape = factor(shape), color = color
      ),
      position = position_jitter(width = 0.12, height = 0),
      alpha = point_alpha,
      size = 0.2
    ) +
    scale_shape_manual(
      values = c(`15` = 15, `16` = 16),
      guide = "none"
    ) +
    scale_color_identity(guide = "none") +
    geom_errorbar(
      data = pred_sum,
      aes(x = .data[[cat_var]], ymin = proportion_lower, ymax = proportion_upper),
      width = 0.2, linewidth = 0.2
    ) +
    geom_point(
      data = pred_sum,
      aes(x = .data[[cat_var]], y = proportion_mean),
      size = 0.6
    ) +
    labs(
      title = plot_title,
      x = plot_x_label,
      y = "Proportion"
    ) +
    theme_classic(base_family = "DejaVu Sans", base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, colour = "black", size = base_size + 1),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
      axis.title.x = element_text(margin = margin(t = 4), colour = "black", size = base_size + 2),
      axis.title.y = element_text(margin = margin(r = 4), colour = "black", size = base_size + 2),
      axis.line = element_line(linewidth = 0, colour = "black"),
      axis.text.x = element_text(colour = "black", size = base_size, margin = margin(t = 2)),
      axis.text.y = element_text(colour = "black", size = base_size, margin = margin(r = 4)),
      axis.ticks.length.x = unit(0.1, "cm"),
      axis.ticks.length.y = unit(0.1, "cm"))

  # Add y scale transform with optional breaks
  if (is.null(y_breaks)) {
    p <- p + scale_y_sqrt(labels = y_labels)
  } else {
    p <- p + scale_y_sqrt(labels = y_labels, breaks = y_breaks)
  }

  p
}

## snuc cell type Para & Race predictions
# snuc CD4+ T Para predictions
p_cd4 <- plot_cat_obs_plus_sccomp(
  fit_sccomp = sccomp_res_snuc,
  out_cat   = out_snuc_para,
  celltype   = "CD4+ T",
  cat_var   = "Para_bin",
  levels_to_plot = c("0","1-3"),
  age_var    = "Age",
  ref_age    = 55,
  outliers_tbl = outliers_tbl_snuc,
  y_labels = scales::label_percent(accuracy = 0.01),
  plot_title = "CD4+ T",
  plot_x_label = "Parity",
  y_breaks = c(0, 0.0001, 0.001, 0.002, 0.003, 0.004)
)
p_cd4
ggsave(
  filename = paste0(outs_dir, "/sccomp_para_modelpredict_snuc_cd4.pdf"), plot = p_cd4,
  width = 30, height = 45, units = "mm", device = cairo_pdf, dpi = 500
)

# snuc TNF-activated Stromal Para predictions
p_tnfstromal <- plot_cat_obs_plus_sccomp(
  fit_sccomp = sccomp_res_snuc,
  out_cat   = out_snuc_para,
  celltype   = "TNF-activated Stromal",
  cat_var   = "Para_bin",
  levels_to_plot = c("0","1-3"),
  age_var    = "Age",
  ref_age    = 55,
  outliers_tbl = outliers_tbl_snuc,
  y_labels = scales::label_percent(accuracy = 1),
  plot_title = "TNF-activated\nStromal",
  plot_x_label = "Parity",
  y_breaks = c(0, 0.01, 0.05, 0.1, 0.2, 0.3)
)
p_tnfstromal
ggsave(
  filename = paste0(outs_dir, "/sccomp_para_modelpredict_snuc_tnfactivatedstromal.pdf"), plot = p_tnfstromal,
  width = 30, height = 45, units = "mm", device = cairo_pdf, dpi = 500
)

# snuc Steroidogenic Stromal Race predictions
race_snuc_steroidogenicstromal <- plot_cat_obs_plus_sccomp(
  fit_sccomp = sccomp_res_snuc,
  out_cat   = out_snuc_race,
  celltype   = "Steroidogenic Stromal",
  cat_var   = "Simplified_Race_Ethnicity",
  levels_to_plot = c("White","Black"),
  age_var    = "Age",
  ref_age    = 55,
  outliers_tbl = outliers_tbl_snuc,
  y_labels = scales::label_percent(accuracy = 0.1),
  plot_title = "Steroidogenic\nStromal",
  plot_x_label = "Race",
  y_breaks = c(0, 0.001, 0.01, 0.02)
)
race_snuc_steroidogenicstromal
ggsave(
  filename = paste0(outs_dir, "/sccomp_race_modelpredict_snuc_steroidogenicstromal.pdf"), plot = race_snuc_steroidogenicstromal,
  width = 30, height = 45, units = "mm", device = cairo_pdf, dpi = 500
)

# snuc Contractile Stromal Race predictions
race_snuc_contractilestromal <- plot_cat_obs_plus_sccomp(
  fit_sccomp = sccomp_res_snuc,
  out_cat   = out_snuc_race,
  celltype   = "Contractile Stromal",
  cat_var   = "Simplified_Race_Ethnicity",
  levels_to_plot = c("White","Black"),
  age_var    = "Age",
  ref_age    = 55,
  outliers_tbl = outliers_tbl_snuc,
  y_labels = scales::label_percent(accuracy = 1),
  plot_title = "Contractile\nStromal",
  plot_x_label = "Race",
  y_breaks = c(0, 0.01, 0.1, 0.2, 0.3)
)
race_snuc_contractilestromal
ggsave(
  filename = paste0(outs_dir, "/sccomp_race_modelpredict_snuc_contractilestromal.pdf"), plot = race_snuc_contractilestromal,
  width = 30, height = 45, units = "mm", device = cairo_pdf, dpi = 500
)

## bulk cell type Race predictions
# bulk Steroidogenic Stromal Race predictions
race_bulk_steroidogenicstromal <- plot_cat_obs_plus_sccomp(
  fit_sccomp = sccomp_res_bulk,
  out_cat   = out_bulk_race,
  celltype   = "Steroidogenic Stromal",
  cat_var   = "RACE",
  levels_to_plot = c("White","Black"),
  age_var    = "AGE",
  ref_age    = 55,
  outliers_tbl = outliers_tbl_bulk,
  y_labels = scales::label_percent(accuracy = 0.1),
  plot_title = "Steroidogenic\nStromal",
  plot_x_label = "Race",
  y_breaks = c(0, 0.001,0.01, 0.05, 0.1)
)
race_bulk_steroidogenicstromal
ggsave(
  filename = paste0(outs_dir, "/sccomp_race_modelpredict_bulk_steroidogenicstromal.pdf"), plot = race_bulk_steroidogenicstromal,
  width = 30, height = 45, units = "mm", device = cairo_pdf, dpi = 500
)

# bulk Contractile Stromal Race predictions
race_bulk_contractilestromal <- plot_cat_obs_plus_sccomp(
  fit_sccomp = sccomp_res_bulk,
  out_cat   = out_bulk_race,
  celltype   = "Contractile Stromal",
  cat_var   = "RACE",
  levels_to_plot = c("White","Black"),
  age_var    = "AGE",
  ref_age    = 55,
  outliers_tbl = outliers_tbl_bulk,
  y_labels = scales::label_percent(accuracy = 1),
  plot_title = "Contractile\nStromal",
  plot_x_label = "Race",
  y_breaks = c(0, 0.01, 0.1, 0.2, 0.4, 0.6)
)
race_bulk_contractilestromal
ggsave(
  filename = paste0(outs_dir, "/sccomp_race_modelpredict_bulk_contractilestromal.pdf"), plot = race_bulk_contractilestromal,
  width = 30, height = 45, units = "mm", device = cairo_pdf, dpi = 500
)

### Continuous: plot observed proportions + sccomp predictions for a given cell type and continuous variable (e.g. Age)
# Function to plot sccomp marginal age curve as fold-change vs a reference age
plot_sccomp_fc_curve <- function(
  pred_tbl,
  cell_group,
  plot_title,
  base_size = 8,
  ref_age = 55,
  eps = 1e-8
) {
  stopifnot(
    all(c("Age", "cell_group", "proportion_mean", "proportion_lower", "proportion_upper") %in% names(pred_tbl))
  )

  # 1) Marginalize over donor/plex/etc. already present in pred_tbl by averaging at each Age
  curve <- pred_tbl |>
    dplyr::filter(.data$cell_group == .env$cell_group) |>
    dplyr::group_by(.data$Age) |>
    dplyr::summarise(
      proportion_mean  = mean(.data$proportion_mean),
      proportion_lower = mean(.data$proportion_lower),
      proportion_upper = mean(.data$proportion_upper),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$Age)

  if (nrow(curve) < 2) stop("Not enough rows for this cell_group after filtering.")

  # 2) Lognormal / delta-style approximation for log2 fold-change vs ref_age
  log_mu <- log(curve$proportion_mean + eps)
  log_lo <- log(curve$proportion_lower + eps)
  log_hi <- log(curve$proportion_upper + eps)
  log_sd <- (log_hi - log_lo) / (2 * 1.96)

  ref_log_mu <- stats::approx(curve$Age, log_mu, xout = ref_age, rule = 2)$y
  ref_log_sd <- stats::approx(curve$Age, log_sd, xout = ref_age, rule = 2)$y

  fc_log_mu <- log_mu - ref_log_mu
  fc_log_sd <- sqrt(log_sd^2 + ref_log_sd^2)

  curve_fc <- curve |>
    dplyr::mutate(
      fc_mean  = fc_log_mu / log(2),
      fc_lower = (fc_log_mu - 1.96 * fc_log_sd) / log(2),
      fc_upper = (fc_log_mu + 1.96 * fc_log_sd) / log(2)
    )

  pick_fc_breaks <- function(yrng, min_breaks = 3) {
    stopifnot(length(yrng) == 2, is.finite(yrng[1]), is.finite(yrng[2]))
    y0 <- yrng[1]; y1 <- yrng[2]
    if (y0 > y1) { tmp <- y0; y0 <- y1; y1 <- tmp }

    # Always include 0 (1× reference)
    candidates <- c(-3,-2,-1,0,1,2,3)

    br <- candidates[candidates >= y0 & candidates <= y1]

    # If br contains 6 values, reduce to lowest, 0 and second highest
    if (length(br) == 6) {
      second_highest <- sort(br, decreasing = TRUE)[2]
      br <- sort(unique(c(min(br), 0, second_highest)))
    }
    # If still only 1-2 breaks because range is tiny, force a symmetric set around 0
    if (length(br) < min_breaks) {
      br <- c(-1, 0, 1)
    }

    br
  }

  fc_labels_x <- function(breaks) {
    paste0(format(round(2^breaks, 2), trim = TRUE), "×")
  }

  # 3) Choose nice log2 breaks and × labels
  yrng <- range(c(curve_fc$fc_lower, curve_fc$fc_upper), finite = TRUE)
  brks <- pick_fc_breaks(yrng, min_breaks = 3)

  # 4) Plot
  ggplot() +
    geom_ribbon(
      data = curve_fc,
      aes(x = .data$Age, ymin = .data$fc_lower, ymax = .data$fc_upper),
      alpha = 0.15
    ) +
    geom_hline(
      yintercept = 0,
      linewidth = 0.2,
      color = "grey80",
      linetype = "dotted"
    ) +
    geom_line(
      data = curve_fc,
      aes(x = .data$Age, y = .data$fc_mean),
      color = "black",
      linewidth = 0.5
    ) +
    scale_y_continuous(
      name = "Fold Change",
      breaks = brks,
      labels = fc_labels_x(brks),
      limits = yrng
    ) +
    labs(title = plot_title, x = "Age") +
    theme_classic(base_family = "DejaVu Sans", base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, colour = "black", size = base_size + 1),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
      axis.title.x = element_text(margin = margin(t = 4), colour = "black", size = base_size + 2),
      axis.title.y = element_text(margin = margin(r = 4), colour = "black", size = base_size + 2),
      axis.line = element_line(linewidth = 0, colour = "black"),
      axis.text.x = element_text(colour = "black", size = base_size),
      axis.text.y = element_text(colour = "black", size = base_size),
      axis.ticks.length.x = unit(0.1, "cm"),
      axis.ticks.length.y = unit(0.1, "cm")
    )
}

## snuc Age predictions
snuc_plotting_list <- list(
  "Steroidogenic Stromal" = "Steroidogenic\nStromal",
  "CD4+ T" = "CD4+ T",
  "LympEnd" = "LympEnd",
  "Homeostatic Fibroblast" = "Homeostatic\nFibroblast",
  "IER Stromal" = "IER Stromal",
  "TNF-activated Stromal" = "TNF-activated\nStromal"
)

for (cell_group in names(snuc_plotting_list)) {
  p <- plot_sccomp_fc_curve(
    pred_tbl   = out_snuc_age$pred,
    cell_group = cell_group,
    plot_title = snuc_plotting_list[[cell_group]],
    base_size  = 5,
    ref_age    = 55
  )
  ggsave(filename = paste0(outs_dir, "/sccomp_age_modelpredict_snuc_", gsub(" ", "_", cell_group), ".pdf"), plot = p,
  width = 30, height = 30, units = "mm", device = cairo_pdf, dpi = 500)
}

# Sanitize age column name
out_bulk_age$pred <- out_bulk_age$pred %>% dplyr::rename(Age = AGE)
bulk_plotting_list <- list(
  "Steroidogenic Stromal" = "Steroidogenic\nStromal",
  "CD4+ T" = "CD4+ T",
  "LympEnd" = "LympEnd"
)

for (cell_group in names(bulk_plotting_list)) {
  p <- plot_sccomp_fc_curve(
    pred_tbl   = out_bulk_age$pred,
    cell_group = cell_group,
    plot_title = bulk_plotting_list[[cell_group]],
    base_size  = 5,
    ref_age    = 55
  )
  ggsave(filename = paste0(outs_dir, "/sccomp_age_modelpredict_bulk_", gsub(" ", "_", cell_group), ".pdf"), plot = p,
  width = 30, height = 30, units = "mm", device = cairo_pdf, dpi = 500)
}



#### 8. Assess correlation between Steroidogenic Stromal and Contractile Stromal via CLR-transformed proportions ####

## Helpers for handling bulk and snuc data for CLR analysis
# Bulk prep:
prep_clr_input_from_prop_df <- function(
  prop_df,
  sample_col = "sample_name",
  celltype_col = "Celltype",
  prop_col = "prop"
) {
  stopifnot(all(c(sample_col, celltype_col, prop_col) %in% colnames(prop_df)))

  out <- prop_df %>%
    dplyr::select(
      sample_name = all_of(sample_col),
      Celltype = all_of(celltype_col),
      prop = all_of(prop_col)
    ) %>%
    dplyr::group_by(sample_name) %>%
    dplyr::mutate(prop = prop / sum(prop, na.rm = TRUE)) %>%
    dplyr::ungroup()

  return(out)
}

# snuc prep from SingleCellExperiment:
# Builds sample-level cell-type proportions from cell metadata
prep_clr_input_from_sce <- function(
  sce_obj,
  sample_col = "unique_sample_ID",
  celltype_col = "cell_type_fine"
) {
  meta <- as.data.frame(colData(sce_obj))

  stopifnot(all(c(sample_col, celltype_col) %in% colnames(meta)))

  out <- meta %>%
    dplyr::select(
      sample_name = all_of(sample_col),
      Celltype = all_of(celltype_col)
    ) %>%
    dplyr::count(sample_name, Celltype, name = "n_cells") %>%
    dplyr::group_by(sample_name) %>%
    dplyr::mutate(prop = n_cells / sum(n_cells)) %>%
    dplyr::ungroup() %>%
    dplyr::select(sample_name, Celltype, prop)

  return(out)
}

## Core CLR analysis
run_clr_pairwise_analysis <- function(
  clr_input_df,
  outs_dir = ".",
  out_prefix = "clr",
  pseudocount = 1e-6
) {
  stopifnot(all(c("sample_name", "Celltype", "prop") %in% colnames(clr_input_df)))

  # wide matrix: rows = samples, cols = cell types
  wide_df <- clr_input_df %>%
    tidyr::pivot_wider(
      names_from = Celltype,
      values_from = prop,
      values_fill = 0
    ) %>%
    as.data.frame()

  rownames(wide_df) <- wide_df$sample_name
  wide_df$sample_name <- NULL

  prop_mat <- as.matrix(wide_df)

  # pseudocount + renormalize if needed
  if (any(prop_mat == 0, na.rm = TRUE)) {
    prop_mat <- prop_mat + pseudocount
    prop_mat <- prop_mat / rowSums(prop_mat)
  }

  # CLR transform within sample
  clr_transform <- function(x) {
    lx <- log(x)
    lx - mean(lx)
  }

  clr_mat <- t(apply(prop_mat, 1, clr_transform))
  colnames(clr_mat) <- colnames(prop_mat)
  rownames(clr_mat) <- rownames(prop_mat)

  celltypes <- colnames(clr_mat)

  # unique unordered pairwise tests
  pairwise_unique <- lapply(seq_along(celltypes), function(i) {
    lapply(seq_along(celltypes), function(j) {
      if (j <= i) return(NULL)

      ct1 <- celltypes[i]
      ct2 <- celltypes[j]

      keep <- complete.cases(clr_mat[, ct1], clr_mat[, ct2])

      if (sum(keep) < 3) {
        return(data.frame(
          Celltype1 = ct1,
          Celltype2 = ct2,
          n = sum(keep),
          rho = NA_real_,
          p = NA_real_,
          stringsAsFactors = FALSE
        ))
      }

      tt <- suppressWarnings(
        cor.test(
          clr_mat[keep, ct1],
          clr_mat[keep, ct2],
          method = "spearman",
          exact = FALSE
        )
      )

      data.frame(
        Celltype1 = ct1,
        Celltype2 = ct2,
        n = sum(keep),
        rho = unname(tt$estimate),
        p = tt$p.value,
        stringsAsFactors = FALSE
      )
    })
  }) %>%
    unlist(recursive = FALSE) %>%
    bind_rows() %>%
    mutate(
      FDR = p.adjust(p, method = "BH"),
      abs_rho = abs(rho)
    ) %>%
    arrange(FDR, desc(abs_rho))

  # double-sided review table for Excel filtering
  pairwise_review <- bind_rows(
    pairwise_unique %>%
      transmute(
        Celltype = Celltype1,
        Partner_Celltype = Celltype2,
        n, rho, p, FDR, abs_rho
      ),
    pairwise_unique %>%
      transmute(
        Celltype = Celltype2,
        Partner_Celltype = Celltype1,
        n, rho, p, FDR, abs_rho
      )
  ) %>%
    arrange(Celltype, FDR, desc(abs_rho))

  # export CSVs
  write.csv(
    clr_input_df,
    file = file.path(outs_dir, paste0(out_prefix, "_clr_input_long.csv")),
    row.names = FALSE
  )

  write.csv(
    pairwise_unique,
    file = file.path(outs_dir, paste0(out_prefix, "_pairwise_clr_unique.csv")),
    row.names = FALSE
  )

  write.csv(
    pairwise_review,
    file = file.path(outs_dir, paste0(out_prefix, "_pairwise_clr_review_double_sided.csv")),
    row.names = FALSE
  )

  return(list(
    clr_input_df = clr_input_df,
    wide_df = wide_df,
    prop_mat = prop_mat,
    clr_mat = clr_mat,
    pairwise_unique = pairwise_unique,
    pairwise_review = pairwise_review
  ))
}

## Pair-specific extractor
get_clr_pair_result <- function(clr_result, celltype_x, celltype_y) {
  clr_mat <- clr_result$clr_mat
  pairwise_unique <- clr_result$pairwise_unique

  stopifnot(celltype_x %in% colnames(clr_mat))
  stopifnot(celltype_y %in% colnames(clr_mat))

  plot_df <- data.frame(
    sample_name = rownames(clr_mat),
    x = clr_mat[, celltype_x],
    y = clr_mat[, celltype_y],
    stringsAsFactors = FALSE
  )

  pair_row <- pairwise_unique %>%
    filter(
      (Celltype1 == celltype_x & Celltype2 == celltype_y) |
        (Celltype1 == celltype_y & Celltype2 == celltype_x)
    )

  if (nrow(pair_row) != 1) {
    stop("Could not uniquely identify requested pair in pairwise_unique table.")
  }

  list(
    plot_df = plot_df,
    stats = pair_row
  )
}

## Plotter for CLR pair
plot_clr_pair <- function(
  clr_result,
  celltype_x,
  celltype_y,
  title = "Cross-sample co-abundance (all samples)",
  point_alpha = 0.5,
  base_size = 5,
  point_size = 0.8,
  x_breaks = NULL,
  y_breaks = NULL
) {
  pair_obj <- get_clr_pair_result(clr_result, celltype_x, celltype_y)
  plot_df <- pair_obj$plot_df
  stats <- pair_obj$stats

  x_rng <- range(plot_df$x, na.rm = TRUE)
  y_rng <- range(plot_df$y, na.rm = TRUE)

  annot_size <- base_size / .pt

  x_annot <- x_rng[1] + 0.03 * diff(x_rng)
  y_annot_top <- y_rng[2] - 0.03 * diff(y_rng)
  y_annot_bottom <- y_rng[2] - 0.11 * diff(y_rng)

  # Set up x scale with optional breaks
  x_scale <- if (!is.null(x_breaks)) {
    scale_x_continuous(breaks = x_breaks)
  } else {
    scale_x_continuous()
  }

  # Set up y scale with optional breaks
  y_scale <- if (!is.null(y_breaks)) {
    scale_y_continuous(breaks = y_breaks)
  } else {
    scale_y_continuous()
  }

  p <- ggplot(plot_df, aes(x = x, y = y)) +
    geom_point(alpha = point_alpha, size = point_size, stroke = 0) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5) +
    annotate(
      "text",
      x = x_annot,
      y = y_annot_top,
      label = paste0("rho == ", sprintf("%.2f", stats$rho)),
      parse = TRUE,
      hjust = 0,
      vjust = 1,
      size = annot_size,
      family = "DejaVu Sans",
      colour = "black"
    ) +
    annotate(
      "text",
      x = x_annot,
      y = y_annot_bottom,
      label = paste0("italic(p)[adj] == '", sprintf("%.1e", stats$FDR), "'"),
      parse = TRUE,
      hjust = 0,
      vjust = 1,
      size = annot_size,
      family = "DejaVu Sans",
      colour = "black"
    ) +
    x_scale +
    y_scale +
    theme_classic(base_family = "DejaVu Sans", base_size = base_size) +
    labs(
      title = title,
      x = paste0("CLR ", celltype_x),
      y = paste0("CLR ", celltype_y)
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, colour = "black", size = base_size + 1),
      axis.title.x = element_text(margin = margin(t = 5), colour = "black", size = base_size+1),
      axis.title.y = element_text(margin = margin(r = 5), colour = "black", size = base_size+1),
      axis.text.x = element_text(colour = "black", size = base_size),
      axis.text.y = element_text(colour = "black", size = base_size),
      axis.ticks = element_line(colour = "black", linewidth = 0.3),
      axis.ticks.length.x = unit(0.1, "cm"),
      axis.ticks.length.y = unit(0.1, "cm"),
      axis.line = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4)
    )

  return(p)
}

## Running CLR analysis on bulk data
bulk_clr_input <- prep_clr_input_from_prop_df(bulk_prop_df, sample_col = "sample_name", celltype_col = "Celltype", prop_col = "prop")
bulk_clr_res <- run_clr_pairwise_analysis(clr_input_df = bulk_clr_input, outs_dir = outs_dir, out_prefix = "bulk")
bulk_target <- get_clr_pair_result(bulk_clr_res, celltype_x = "Steroidogenic Stromal", celltype_y = "Contractile Stromal")
bulk_target$stats
p_bulk <- plot_clr_pair(bulk_clr_res, celltype_x = "Steroidogenic Stromal", celltype_y = "Contractile Stromal",
  title = "Cross-sample co-abundance (all samples)", y_breaks = c(0, 4, 8))
ggsave(filename = paste0(outs_dir, "/bulk_steroidogenic_contractile_clr_plot.pdf"), plot = p_bulk, width = 30, height = 40, units = "mm",
    device = cairo_pdf, dpi = 500)

## Running CLR analysis on snuc data
snuc_clr_input <- prep_clr_input_from_sce(sce_obj, sample_col = "unique_sample_ID", celltype_col = "cell_type_fine")
snuc_clr_res <- run_clr_pairwise_analysis(clr_input_df = snuc_clr_input, outs_dir = outs_dir, out_prefix = "snuc")
snuc_target <- get_clr_pair_result(snuc_clr_res, celltype_x = "Steroidogenic Stromal", celltype_y = "Contractile Stromal")
snuc_target$stats
p_snuc <- plot_clr_pair(snuc_clr_res, celltype_x = "Steroidogenic Stromal", celltype_y = "Contractile Stromal",
  title = "Cross-sample co-abundance (all samples)", x_breaks = c(-6, -3, 0, 3))
ggsave(filename = paste0(outs_dir, "/snuc_steroidogenic_contractile_clr_plot.pdf"), plot = p_snuc, width = 30, height = 40, units = "mm",
  device = cairo_pdf, dpi = 500)



#### 9. Supplementary Figs showing observed props and model fit to observed props ####

# Patched local sccomp plotting helpers are used only for visualization so that
# sccomp-censored outlier sample × cell-group pairs are displayed as red squares.
# Statistical results are unchanged.
source(file.path(r_analyses_dir, "sccomp_outlier_plot_helpers.R"))

# Attach sccomp-censored outlier sample × cell-group pairs
sccomp_res_snuc_plot <- attach_sccomp_outliers(sccomp_res_snuc, outliers_tbl_snuc)
sccomp_res_bulk_plot <- attach_sccomp_outliers(sccomp_res_bulk, outliers_tbl_bulk)

# Generate patched sccomp plots with outliers shown
plots_snuc <- plot_sccomp_patched(sccomp_res_snuc_plot, significance_statistic = "FDR", significance_threshold = 0.05, show_fdr_message = FALSE,
  marginal_predictions = list(Age = out_snuc_age$pred),
  marginal_x_cols = list(Age = "Age"))

plots_bulk <- plot_sccomp_patched(sccomp_res_bulk_plot, significance_statistic = "FDR", significance_threshold = 0.05,  show_fdr_message = FALSE,
  marginal_predictions = list(AGE = out_bulk_age$pred),
  marginal_x_cols = list(AGE = "AGE"))

# Ensure plot list names match model factors
name_sccomp_boxplots <- function(fit_sccomp, plot_list) {
  factor_names <- fit_sccomp %>%
    dplyr::filter(!is.na(factor)) %>%
    dplyr::distinct(factor) %>%
    dplyr::pull(factor)

  names(plot_list$boxplot) <- factor_names
  plot_list
}

plots_snuc <- name_sccomp_boxplots(sccomp_res_snuc_plot, plots_snuc)
plots_bulk <- name_sccomp_boxplots(sccomp_res_bulk_plot, plots_bulk)

# Save larger raster plots for supplementary review and figure assembly
ggsave(filename = file.path(outs_dir, "sccomp_rawprops_model_fit_snuc_age.png"), plot = plots_snuc$boxplot[["Age"]],
  width = 280, height = 180, units = "mm", dpi = 500)

ggsave(filename = file.path(outs_dir, "sccomp_rawprops_model_fit_bulk_age.png"), plot = plots_bulk$boxplot[["AGE"]],
  width = 280, height = 180, units = "mm")

ggsave(filename = file.path(outs_dir, "sccomp_rawprops_model_fit_snuc_parity.png"), plot = plots_snuc$boxplot[["Para_bin"]],
  width = 280, height = 180, units = "mm", dpi = 500)




#### 10. HMC and pseudocount sensitivity analyses ####

set.seed(42)
sccomp_res_snuc_hmc <- sce_obj_sccomp %>%
  sccomp_estimate(
    formula_composition = ~ Age + Para_bin + BMI + Plex + Simplified_Race_Ethnicity + StressScore + (1 | SenNet_ID),
    formula_variability = ~ 1,
    sample = "unique_sample_ID",
    cell_group = "cell_group",
    cores = 100,
    mcmc_seed = 42,
    bimodal_mean_variability_association = TRUE,
    inference_method = "hmc",
    cleanup_draw_files = FALSE,
    output_directory = paste0(r_analyses_dir, "/R_sccomp_cache")
  ) %>%
  sccomp_remove_outliers(mcmc_seed = 42, cores = 100, cleanup_draw_files = FALSE, output_directory = paste0(r_analyses_dir, "/R_sccomp_cache")) %>%
  sccomp_test()

saveRDS(sccomp_res_snuc_hmc, file = paste0(outs_dir, "/sccomp_res_snuc_hmc.rds"))

## Extract raw draws for sccomp fits and plot as halfeye density plot
# Get parameter and cell group names from the design matrix
parameters_lookup <- colnames(attr(sccomp_res_snuc_hmc, "model_input")$X)
cell_groups_lookup <- colnames(attr(sccomp_res_snuc_hmc, "model_input")$y)

# Extract raw draws
age_raw_draws_for_plot <- sccomp_res_snuc_hmc %>%
  attr("fit") %>%
  gather_draws(beta[parameter_idx, cell_group_idx]) %>%
  mutate(
    cell_group = cell_groups_lookup[cell_group_idx],
    parameter = parameters_lookup[parameter_idx]
  ) %>%
  # Filter for your factor of interest and relevant cell groups
  filter(parameter == "Age")

# Extract raw draws
para_raw_draws_for_plot <- sccomp_res_snuc_hmc %>%
  attr("fit") %>%
  gather_draws(beta[parameter_idx, cell_group_idx]) %>%
  mutate(
    cell_group = cell_groups_lookup[cell_group_idx],
    parameter = parameters_lookup[parameter_idx]
  ) %>%
  # Filter for your factor of interest and relevant cell groups
  filter(parameter == "Para_bin1-3")

# Plot halfeye density plot for Age
age_hmc_halfeye <- ggplot(age_raw_draws_for_plot, aes(x = .value, y = reorder(cell_group, .value, median))) +
  # Shaded density + Point Interval
  stat_halfeye(
    aes(fill = after_stat(x < 0)),
    .width = c(0.66, 0.95),  # Shows 66% and 95% intervals
    point_interval = median_qi,
    alpha = 0.7
  ) +

  # Reference lines
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "red3", alpha = 0.4) +

  # Colors
  scale_fill_manual(values = c("TRUE" = "#377EB8", "FALSE" = "#E41A1C"),
                    labels = c("TRUE" = "Decrease", "FALSE" = "Increase"),
                    name = "Direction") +

  labs(x = "Age Effect Size (Logit Fold Change)", y = NULL) +
  theme_minimal()

# Plot halfeye density plot for Para_bin1-3
para_hmc_halfeye <- ggplot(para_raw_draws_for_plot, aes(x = .value, y = reorder(cell_group, .value, median))) +
  # Shaded density + Point Interval
  stat_halfeye(
    aes(fill = after_stat(x < 0)),
    .width = c(0.66, 0.95),  # Shows 66% and 95% intervals
    point_interval = median_qi,
    alpha = 0.7
  ) +

  # Reference lines
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "red3", alpha = 0.4) +

  # Colors
  scale_fill_manual(values = c("TRUE" = "#377EB8", "FALSE" = "#E41A1C"),
                    labels = c("TRUE" = "Decrease", "FALSE" = "Increase"),
                    name = "Direction") +

  labs(x = "Para_bin1-3 Effect Size Relative to Para_bin0 (Logit Fold Change)", y = NULL) +
  theme_minimal()

ggsave(filename = paste0(outs_dir, "/sccomp_age_hmc_halfeye.png"), plot = age_hmc_halfeye, width = 180,
  height = 160, units = "mm", dpi = 500)

ggsave(
  filename = paste0(outs_dir, "/sccomp_para_hmc_halfeye.png"), plot = para_hmc_halfeye, width = 180,
  height = 160, units = "mm", dpi = 500
)



#### 11. Bulk proportion multiplier/pseudocount sensitivity analysis ####

bulk_10k_df <- bulk_prop_df %>% mutate(prop = round(prop * 10000))
bulk_1M_df <- bulk_prop_df %>% mutate(prop = round(prop * 1000000))
bulk_10k_df$prop <- as.integer(bulk_10k_df$prop)
bulk_1M_df$prop <- as.integer(bulk_1M_df$prop)

zero_prop_stats_10k <- bulk_10k_df %>%
  dplyr::group_by(cell_group) %>%
  dplyr::summarise(
    n = dplyr::n(),
    n_zero = sum(prop == 0, na.rm = TRUE),
    pct_zero = 100 * mean(prop == 0, na.rm = TRUE)
  ) %>%
  dplyr::arrange(desc(pct_zero))

zero_prop_stats_1M <- bulk_1M_df %>%
  dplyr::group_by(cell_group) %>%
  dplyr::summarise(
    n = dplyr::n(),
    n_zero = sum(prop == 0, na.rm = TRUE),
    pct_zero = 100 * mean(prop == 0, na.rm = TRUE)
  ) %>%
  dplyr::arrange(desc(pct_zero))

print(zero_prop_stats_10k, n = 32)
print(zero_prop_stats_1M, n = 32)

# Remove cell types (cell_group) that have zeros in >50% samples in bulk_100k_df
celltypes_keep_10k <- zero_prop_stats_10k %>%
  dplyr::filter(pct_zero <= 50) %>%
  dplyr::pull(cell_group)

celltypes_keep_1M <- zero_prop_stats_1M %>%
  dplyr::filter(pct_zero <= 50) %>%
  dplyr::pull(cell_group)

bulk_10k_df <- bulk_10k_df %>%
  dplyr::filter(cell_group %in% celltypes_keep_10k)

bulk_1M_df <- bulk_1M_df %>%
  dplyr::filter(cell_group %in% celltypes_keep_1M)

# Identify the cell types removed because they had zeros in >50% samples
celltypes_removed_10k <- zero_prop_stats_10k %>%
  dplyr::filter(!(cell_group %in% celltypes_keep_10k))

celltypes_removed_1M <- zero_prop_stats_1M %>%
  dplyr::filter(!(cell_group %in% celltypes_keep_1M))

cat("Cell types removed due to >50% zero values:\n")
print(celltypes_removed_10k[, c("cell_group", "pct_zero")], row.names = FALSE)
# Cell types removed due to >50% zero values:
# A tibble: 3 × 2
#   cell_group               pct_zero
#   <chr>                       <dbl>
# 1 CD16- NK                     76.9
# 2 Mast                         64.4
# 3 Granulocyte-myeloid Sink     63.5

print(celltypes_removed_1M[, c("cell_group", "pct_zero")], row.names = FALSE)
# Cell types removed due to >50% zero values:
# # A tibble: 0 × 2
#   cell_group pct_zero
#   <chr>         <dbl>

set.seed(42)
sccomp_res_bulk_10k <- sccomp::sccomp_estimate(
  bulk_10k_df,
  formula_composition = form_composition_bulk,
  formula_variability = form_variability_bulk,
  sample = "sample_name",
  cell_group = "cell_group",
  cores = 1, # run single-threaded for deterministic results
  mcmc_seed = 42,
  abundance = "prop",
  inference_method = "pathfinder", # default for sccomp_estimate
  cleanup_draw_files = FALSE,
  output_directory = paste0(r_analyses_dir, "/R_sccomp_cache")
)

set.seed(42)
sccomp_res_bulk_10k <- sccomp_res_bulk_10k %>%
  sccomp_remove_outliers(mcmc_seed = 42, cores = 1, cleanup_draw_files = FALSE, output_directory = paste0(r_analyses_dir, "/R_sccomp_cache")) %>%
  sccomp_test()

## sccomp analysis with 1M pseudocounts
set.seed(42)
sccomp_res_bulk_1M <- sccomp::sccomp_estimate(
  bulk_1M_df,
  formula_composition = form_composition_bulk,
  formula_variability = form_variability_bulk,
  sample = "sample_name",
  cell_group = "cell_group",
  cores = 1, # run single-threaded for deterministic results
  mcmc_seed = 42,
  abundance = "prop",
  inference_method = "pathfinder", # default for sccomp_estimate
  cleanup_draw_files = FALSE,
  output_directory = paste0(r_analyses_dir, "/R_sccomp_cache")
)

set.seed(42)
sccomp_res_bulk_1M <- sccomp_res_bulk_1M %>%
  sccomp_remove_outliers(mcmc_seed = 42, cores = 1, cleanup_draw_files = FALSE, output_directory = paste0(r_analyses_dir, "/R_sccomp_cache")) %>%
  sccomp_test()

# Save sccomp results objects to RDS files
saveRDS(sccomp_res_bulk_10k, file = paste0(outs_dir, "/sccomp_res_bulk_10k.rds"))
saveRDS(sccomp_res_bulk_1M, file = paste0(outs_dir, "/sccomp_res_bulk_1M.rds"))

sccomp_res_bulk_10k <- sccomp_res_bulk_10k %>%
  filter(!cell_group %in% bulk_types_exclude_plot)
sccomp_res_bulk_1M <- sccomp_res_bulk_1M %>%
  filter(!cell_group %in% bulk_types_exclude_plot)

# Make individual plot for GTEx deconvolved bulk Age results using 10k pseudocounts
p_10k_age <- plot_sccomp_forest_one(
  sccomp_res_bulk_10k,
  parameter_keep = "AGE",
  which = "c",
  signif_stat = "FDR",
  signif_threshold = 0.05,
  effect_threshold = logit_fold_effect_threshold,
  order_by = "effect",
  base_size = 5,
  title = "Deconvoluted Bulk (GTEx) Aging\n(10k pseudocounts)",
  point_size = 0.8,
  x_breaks = seq(-0.4, 0.4, 0.4),
  digits = 1
)

ggsave(filename = paste0(outs_dir, "/sccomp_forest_10kage_bulk.pdf"), plot = p_10k_age, width = 75,
  height = 90, units = "mm", device = cairo_pdf, dpi = 500
)

# Make individual plot for GTEx deconvolved bulk Age results using 1M pseudocounts
p_1M_age <- plot_sccomp_forest_one(
  sccomp_res_bulk_1M,
  parameter_keep = "AGE",
  which = "c",
  signif_stat = "FDR",
  signif_threshold = 0.05,
  effect_threshold = logit_fold_effect_threshold,
  order_by = "effect",
  base_size = 5,
  title = "Deconvoluted Bulk (GTEx) Aging\n(1M pseudocounts)",
  point_size = 0.8,
  x_breaks = seq(-0.4, 0.4, 0.4),
  digits = 1
)

ggsave(filename = paste0(outs_dir, "/sccomp_forest_1Mage_bulk.pdf"), plot = p_1M_age, width = 75,
  height = 90, units = "mm", device = cairo_pdf, dpi = 500
)


#### 12. Export Supplementary Tables ####

# Helper for safe column selection
select_existing <- function(df, cols) {
  df %>% dplyr::select(dplyr::any_of(cols))
}

# 1. Cell-type inclusion QC for snuc sccomp
sccomp_celltype_inclusion_QC <- ct_stats %>%
  dplyr::mutate(
    eligible_for_snuc_sccomp = eligible,
    exclusion_reason = dplyr::case_when(
      eligible ~ "Included",
      n_cells < min_cells & n_donors < min_donors ~ paste0(
        "Excluded: fewer than ", min_cells, " nuclei and fewer than ", min_donors, " donors"
      ),
      n_cells < min_cells ~ paste0("Excluded: fewer than ", min_cells, " nuclei"),
      n_donors < min_donors ~ paste0("Excluded: fewer than ", min_donors, " donors"),
      TRUE ~ "Excluded"
    )
  ) %>%
  dplyr::rename(
    cell_type = cell_type_fine,
    total_nuclei = n_cells,
    n_donors_with_detected_nuclei = n_donors
  ) %>%
  dplyr::select(
    cell_type,
    total_nuclei,
    n_donors_with_detected_nuclei,
    eligible_for_snuc_sccomp,
    exclusion_reason
  ) %>%
  dplyr::arrange(desc(eligible_for_snuc_sccomp), desc(total_nuclei))

write.csv(
  sccomp_celltype_inclusion_QC,
  file.path(outs_dir, "sccomp_celltype_inclusion_QC.csv"),
  row.names = FALSE
)

# 2. Cell-type inclusion QC for deconvolved bulk sccomp
sccomp_GTEx_celltype_inclusion_QC <- zero_prop_stats_100k %>%
  dplyr::mutate(
    eligible_for_GTEx_sccomp = cell_group %in% celltypes_keep,
    exclusion_reason = dplyr::if_else(
      eligible_for_GTEx_sccomp,
      "Included",
      "Excluded: zero pseudo-count in >50% of samples"
    )
  ) %>%
  dplyr::rename(
    n_samples = n,
    n_zero_pseudocount_samples = n_zero,
    percent_zero_pseudocount_samples = pct_zero
  ) %>%
  dplyr::select(
    cell_group,
    n_samples,
    n_zero_pseudocount_samples,
    percent_zero_pseudocount_samples,
    eligible_for_GTEx_sccomp,
    exclusion_reason
  ) %>%
  dplyr::arrange(desc(eligible_for_GTEx_sccomp), desc(percent_zero_pseudocount_samples))

write.csv(sccomp_GTEx_celltype_inclusion_QC, file.path(outs_dir, "sccomp_GTEx_celltype_inclusion_QC.csv"), row.names = FALSE)

# 3. snucRNA-seq sccomp results
sccomp_snuc_results <- sccomp_res_snuc %>%
  dplyr::mutate(
    dataset = "snucRNA-seq",
    parameter = dplyr::recode(
      parameter,
      "Para_bin1-3" = "Parity_1-3_vs_0",
      "Para_bin4+" = "Parity_4plus_vs_0",
      "Simplified_Race_EthnicityBlack" = "Race_Black_vs_White",
      "Simplified_Race_EthnicityOther" = "Race_Other_vs_White",
      .default = parameter
    )
  ) %>%
  as.data.frame() %>%
  dplyr::filter(!stringr::str_detect(parameter, "___")) %>%
  dplyr::relocate(dataset, .before = 1) %>%
  dplyr::arrange(cell_group, parameter)

write.csv(sccomp_snuc_results, file.path(outs_dir, "sccomp_snuc_results.csv"), row.names = FALSE)

# 4. GTEx deconvolved bulk sccomp results
sccomp_GTEx_deconv_results <- sccomp_res_bulk %>%
  as.data.frame() %>%
  dplyr::mutate(
    dataset = "GTEx_deconvolved_bulk",
    parameter = dplyr::case_when(
      parameter == "AGE" ~ "Age",
      parameter == "RACEBlack" ~ "Race_Black_vs_White",
      parameter == "SMTSISCH" ~ "Ischemic_time",
      parameter == "SMRIN" ~ "RNA_integrity_number",
      stringr::str_detect(parameter, "^DTHHRDY") ~ parameter,
      stringr::str_detect(parameter, "^SMCENTER") ~ parameter,
      TRUE ~ parameter
    )
  ) %>%
  dplyr::filter(!stringr::str_detect(parameter, "___")) %>%
  dplyr::relocate(dataset, .after = parameter) %>%
  dplyr::arrange(cell_group, parameter)

write.csv(sccomp_GTEx_deconv_results, file.path(outs_dir, "sccomp_GTEx_deconv_results.csv"), row.names = FALSE)

# 5. snucRNA-seq sccomp outliers
sccomp_snuc_outliers <- outliers_tbl_snuc %>%
  as.data.frame() %>%
  dplyr::rename(
    sample_id = sample
  ) %>%
  dplyr::select(
    sample_id,
    cell_group,
    cell_count,
    dplyr::any_of(c(
      "SenNet_ID",
      "Age",
      "BMI",
      "Para_bin",
      "Plex",
      "Simplified_Race_Ethnicity",
      "StressScore"
    ))
  ) %>%
    dplyr::mutate(
    Para_bin = dplyr::recode(
      as.character(Para_bin),
      "0" = "0",
      "1-3" = "1_to_3",
      "4+" = "4plus"
    )
  ) %>%
  dplyr::arrange(sample_id, cell_group) %>%
  dplyr::mutate(dataset = "snucRNA-seq") %>%
  dplyr::relocate(dataset, .before = 1)

write.csv(sccomp_snuc_outliers, file.path(outs_dir, "sccomp_snuc_outliers.csv"), row.names = FALSE)

# 6. GTEx deconvolved bulk sccomp outliers
sccomp_GTEx_outliers <- outliers_tbl_bulk %>%
  as.data.frame() %>%
  dplyr::rename(
    gtex_sample_id = sample
  ) %>%
  dplyr::select(
    gtex_sample_id,
    cell_group,
    cell_count
  ) %>%
  dplyr::arrange(gtex_sample_id, cell_group) %>%
  dplyr::mutate(dataset = "GTEx_deconvolved_bulk") %>%
  dplyr::relocate(dataset, .before = 1)

write.csv(sccomp_GTEx_outliers, file.path(outs_dir, "sccomp_GTEx_outliers.csv"), row.names = FALSE)

# 7. Age sccomp prediction table
make_age_pred_tbl <- function(
  pred_obj,
  dataset_label,
  age_col = "Age",
  ref_age = 55,
  eps = 1e-8
) {
  pred_df <- pred_obj$pred %>%
    as.data.frame()

  # Harmonize GTEx AGE -> Age if needed
  if (age_col != "Age" && age_col %in% colnames(pred_df)) {
    pred_df <- pred_df %>%
      dplyr::rename(Age = dplyr::all_of(age_col))
  }

  pred_summary <- pred_df %>%
    dplyr::group_by(cell_group, Age) %>%
    dplyr::summarise(
      proportion_mean  = mean(proportion_mean, na.rm = TRUE),
      proportion_lower = mean(proportion_lower, na.rm = TRUE),
      proportion_upper = mean(proportion_upper, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(cell_group, Age)

  pred_summary %>%
    dplyr::group_by(cell_group) %>%
    dplyr::group_modify(~{
      curve <- .x %>% dplyr::arrange(Age)

      # Same lognormal / delta-style approximation used for plotting
      log_mu <- log(curve$proportion_mean + eps)
      log_lo <- log(curve$proportion_lower + eps)
      log_hi <- log(curve$proportion_upper + eps)
      log_sd <- (log_hi - log_lo) / (2 * 1.96)

      ref_log_mu <- stats::approx(curve$Age, log_mu, xout = ref_age, rule = 2)$y
      ref_log_sd <- stats::approx(curve$Age, log_sd, xout = ref_age, rule = 2)$y

      fc_log_mu <- log_mu - ref_log_mu
      fc_log_sd <- sqrt(log_sd^2 + ref_log_sd^2)

      curve %>%
        dplyr::mutate(
          ref_age = ref_age,
            fold_change_vs_age55_mean  = exp(fc_log_mu),
            fold_change_vs_age55_lower = exp(fc_log_mu - 1.96 * fc_log_sd),
            fold_change_vs_age55_upper = exp(fc_log_mu + 1.96 * fc_log_sd)
        )
    }) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      dataset = dataset_label,
      contrast = "Age"
    ) %>%
    dplyr::select(
      dataset,
      contrast,
      cell_group,
      Age,
      ref_age,
      proportion_mean,
      proportion_lower,
      proportion_upper,
      fold_change_vs_age55_mean,
      fold_change_vs_age55_lower,
      fold_change_vs_age55_upper
    ) %>%
    dplyr::arrange(cell_group, Age)
}

sccomp_age_predictions <- dplyr::bind_rows(
  make_age_pred_tbl(
    pred_obj = out_snuc_age,
    dataset_label = "snucRNA-seq",
    age_col = "Age",
    ref_age = 55
  ),
  make_age_pred_tbl(
    pred_obj = out_bulk_age,
    dataset_label = "GTEx_deconvolved_bulk",
    age_col = "AGE",
    ref_age = 55
  )
) %>%
  dplyr::arrange(factor(dataset, levels = c("snucRNA-seq", "GTEx_deconvolved_bulk")), cell_group, Age)

write.csv(
  sccomp_age_predictions,
  file = file.path(outs_dir, "sccomp_age_predictions.csv"),
  row.names = FALSE
)

# 8. Categorical sccomp prediction tables
make_cat_pred_tbl <- function(pred_obj, dataset_label, contrast_label, level_col) {
  pred_df <- pred_obj$pred %>%
    as.data.frame()

  if (!level_col %in% colnames(pred_df)) {
    stop("Column '", level_col, "' not found in prediction table.")
  }

  pred_df %>%
    dplyr::mutate(level = as.character(.data[[level_col]])) %>%
    dplyr::group_by(cell_group, level) %>%
    dplyr::summarise(
      proportion_mean  = mean(proportion_mean, na.rm = TRUE),
      proportion_lower = mean(proportion_lower, na.rm = TRUE),
      proportion_upper = mean(proportion_upper, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      dataset = dataset_label,
      contrast = contrast_label,
      level = dplyr::recode(
        level,
        "1-3" = "1_to_3",
        "4+" = "4plus",
        .default = level
      )
    ) %>%
    dplyr::select(
      dataset,
      contrast,
      cell_group,
      level,
      proportion_mean,
      proportion_lower,
      proportion_upper
    ) %>%
    dplyr::arrange(cell_group, level)
}

sccomp_categorical_predictions <- dplyr::bind_rows(
  make_cat_pred_tbl(
    pred_obj = out_snuc_para,
    dataset_label = "snucRNA-seq",
    contrast_label = "Parity_bin",
    level_col = "Para_bin"
  ),
  make_cat_pred_tbl(
    pred_obj = out_snuc_race,
    dataset_label = "snucRNA-seq",
    contrast_label = "Race_ethnicity",
    level_col = "Simplified_Race_Ethnicity"
  ),
  make_cat_pred_tbl(
    pred_obj = out_bulk_race,
    dataset_label = "GTEx_deconvolved_bulk",
    contrast_label = "Race_ethnicity",
    level_col = "RACE"
  )
) %>%
  dplyr::arrange(dataset, contrast, cell_group, level)

write.csv(
  sccomp_categorical_predictions,
  file = file.path(outs_dir, "sccomp_categorical_predictions.csv"),
  row.names = FALSE
)

# 9. Combined all-pairs CLR correlation table
make_clr_review_supp_tbl <- function(clr_res, dataset_label, pseudocount = 1e-6) {
  clr_res$pairwise_review %>%
    as.data.frame() %>%
    dplyr::mutate(
      dataset = dataset_label,
    ) %>%
    dplyr::rename(
      cell_type = Celltype,
      partner_cell_type = Partner_Celltype,
      spearman_rho = rho,
      p_value = p,
      bh_fdr = FDR,
      abs_spearman_rho = abs_rho
    ) %>%
    dplyr::select(
      dataset,
      cell_type,
      partner_cell_type,
      n,
      spearman_rho,
      p_value,
      bh_fdr,
      abs_spearman_rho
    ) %>%
    dplyr::arrange(
      dataset,
      cell_type,
      bh_fdr,
      dplyr::desc(abs_spearman_rho)
    )
}

sccomp_clr_pairwise_review_double_sided <- dplyr::bind_rows(
  make_clr_review_supp_tbl(
    clr_res = snuc_clr_res,
    dataset_label = "snucRNA-seq",
    pseudocount = 1e-6
  ),
  make_clr_review_supp_tbl(
    clr_res = bulk_clr_res,
    dataset_label = "GTEx_deconvolved_bulk",
    pseudocount = 1e-6
  )
)

write.csv(
  sccomp_clr_pairwise_review_double_sided,
  file.path(outs_dir, "sccomp_clr_pairwise_review_double_sided.csv"),
  row.names = FALSE
)



#### 13. Session Info ####
sessionInfo()

# R version 4.4.1 (2024-06-14)
# Platform: x86_64-conda-linux-gnu
# Running under: Ubuntu 22.04.5 LTS

# Matrix products: default
# BLAS/LAPACK: .../lib/libopenblasp-r0.3.30.so;  LAPACK version 3.12.0

# locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8
#  [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

# time zone: America/Los_Angeles
# tzcode source: system (glibc)

# attached base packages:
# [1] stats4    stats     graphics  grDevices utils     datasets  methods   base

# other attached packages:
#  [1] magrittr_2.0.4              openxlsx_4.2.8.1            ggh4x_0.3.1                 tidybayes_3.0.7             patchwork_1.3.1
#  [6] ggrepel_0.9.6               purrr_1.1.0                 forcats_1.0.0               reshape2_1.4.4              scales_1.4.0
# [11] stringr_1.5.1               tidyr_1.3.1                 dplyr_1.1.4                 Matrix_1.7-3                compositions_2.0-9
# [16] SingleCellExperiment_1.30.1 SummarizedExperiment_1.38.1 Biobase_2.68.0              GenomicRanges_1.60.0        GenomeInfoDb_1.44.0
# [21] IRanges_2.42.0              S4Vectors_0.48.0            BiocGenerics_0.54.0         generics_0.1.4              MatrixGenerics_1.20.0
# [26] matrixStats_1.5.0           zellkonverter_1.18.0        ggplot2_3.5.2               sccomp_2.1.23               instantiate_0.2.3
# [31] naniar_1.1.0                readxl_1.4.5                readr_2.1.5

# loaded via a namespace (and not attached):
#  [1] tidyselect_1.2.1        svUnit_1.0.8            farver_2.1.2            filelock_1.0.3          tensorA_0.36.2.1        lifecycle_1.0.5
#  [7] processx_3.8.6          posterior_1.6.1         compiler_4.4.1          rlang_1.1.7             tools_4.4.1             S4Arrays_1.8.1
# [13] reticulate_1.42.0       DelayedArray_0.34.1     plyr_1.8.9              RColorBrewer_1.1-3      abind_1.4-8             unigd_0.1.3
# [19] withr_3.0.2             grid_4.4.1              httpgd_2.0.4            MASS_7.3-65             cli_3.6.5               crayon_1.5.3
# [25] robustbase_0.99-7       httr_1.4.7              tzdb_0.5.0              bayesm_3.1-7            parallel_4.4.1          cellranger_1.1.0
# [31] XVector_0.48.0          basilisk_1.20.0         vctrs_0.7.1             jsonlite_2.0.0          dir.expiry_1.16.0       callr_3.7.6
# [37] hms_1.1.3               arrayhelpers_1.1-0      visdat_0.6.0            systemfonts_1.2.3       ggdist_3.3.3            glue_1.8.0
# [43] DEoptimR_1.1-4          ps_1.9.1                distributional_0.5.0    stringi_1.8.7           gtable_0.3.6            UCSC.utils_1.4.0
# [49] tibble_3.3.0            pillar_1.11.0           basilisk.utils_1.20.0   GenomeInfoDbData_1.2.14 R6_2.6.1                lattice_0.22-7
# [55] backports_1.5.0         png_0.1-8               zip_2.3.3               Rcpp_1.1.1              checkmate_2.3.3         coda_0.19-4.1
# [61] SparseArray_1.8.0       fs_1.6.6                pkgconfig_2.0.3



