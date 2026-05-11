## Script name: Human Ovary Atlas - cell-type pseudobulk differential expression and pathway analysis
##
## Purpose: Run cell-type pseudobulk differential expression analysis with dreamlet/variancePartition,
## derive sample-level StressScore covariates, summarize transcriptional perturbation by cell type, and
## perform Hallmark pathway testing with zenith.
##
## Author: Josef Byrne
#########################

#### 1. Import libraries and set directories ####

# Import libraries
box::use(
  readr[read_tsv],
  zellkonverter[readH5AD],
)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(patchwork)
library(grid)
library(BiocParallel)
library(RhpcBLASctl)

library(variancePartition)
library(SingleCellExperiment)
library(dreamlet)
library(zenith)
library(GSEABase)
library(muscat)
library(car) # For VIF check
library(forcats)
library(ggrepel)

# Enable parallel processing
# Don't allow oversubscription of threads which causes thrashing and slowdowns
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)

Sys.setenv(
  OPENBLAS_NUM_THREADS = "1",
  OMP_NUM_THREADS      = "1",
  MKL_NUM_THREADS      = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

RhpcBLASctl::blas_get_num_procs()
RhpcBLASctl::omp_get_num_procs()

parallel::detectCores()
bp <- MulticoreParam(workers = 64)

# Initialize paths
nostromo_datadir <- "/path/to/project_directory"
adata_objdir <- paste0(nostromo_datadir, "/python_analyses/adata_objs")
r_analyses_dir <- paste0(nostromo_datadir, "/R_analyses")
outs_dir <- paste0(r_analyses_dir, "/R_DGE")

#### 2. Import AnnData object and prepare annotations ####
# Import AnnData object
sce_obj <- readH5AD(paste0(adata_objdir, "/final_annotated_atlas.h5ad"), reader = "R", verbose = TRUE)

# Fix ENSG -> gene symbol missed naming from .h5ad 10x flex file (used older 10x Chromium Human Transcriptome Probe Set v1.0.1 GRCh38-2020 probe set file)
probe_df <- read_csv(paste0(nostromo_datadir, "/Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv"), comment = "#", show_col_types = FALSE)

probe_map <- probe_df %>%
  filter(included %in% TRUE) %>%
  mutate(
    gene_id = as.character(gene_id),
    probe_id = as.character(probe_id),
    gene_symbol_raw = str_split_fixed(probe_id, "\\|", 3)[, 2]
  ) %>%
  filter(
    !is.na(gene_id),
    gene_id != "",
    !startsWith(gene_id, "DEPRECATED_"),
    !is.na(gene_symbol_raw),
    gene_symbol_raw != ""
  ) %>%
  distinct(gene_id, gene_symbol_raw) %>%
  group_by(gene_id) %>%
  summarise(
    gene_symbol_raw = gene_symbol_raw[1],
    .groups = "drop"
  ) %>%
  mutate(
    # keep only symbols that look like actual gene symbols,
    # not generic clone / contig / accession-style aliases
    gene_symbol = ifelse(
      str_detect(gene_symbol_raw, "^(AC|AL|AP|BX|FP|RP[0-9]|CT[0-9]|CH[0-9])"),
      NA_character_,
      gene_symbol_raw
    )
  )

old_names <- rownames(sce_obj)
old_ids_clean <- sub("\\..*$", "", old_names)

conversion_df <- tibble(
  row_idx = seq_along(old_names),
  old_name = old_names,
  gene_id = old_ids_clean
) %>%
  left_join(probe_map, by = "gene_id") %>%
  mutate(
    was_ensg = startsWith(old_name, "ENSG"),
    new_name = ifelse(
      was_ensg & !is.na(gene_symbol) & gene_symbol != "",
      gene_symbol,
      old_name
    ),
    converted = old_name != new_name
  )

rowData(sce_obj)$original_name <- old_names
rowData(sce_obj)$gene_id <- old_ids_clean
rowData(sce_obj)$gene_symbol_raw_from_probe <- conversion_df$gene_symbol_raw
rowData(sce_obj)$gene_symbol_fixed <- conversion_df$new_name

new_names <- conversion_df$new_name

if (anyDuplicated(new_names)) {
  cat("\nDuplicate symbols detected after conversion. Making unique with make.unique().\n")
  print(table(new_names[new_names %in% new_names[duplicated(new_names)]]))
  new_names <- make.unique(new_names)
}

rownames(sce_obj) <- new_names

# Manual map remaining ENSGs -> Gene symbols where a clean annotation is available
manual_map <- c(
  "ENSG00000158483" = "FAM86C1",
  "ENSG00000182230" = "FAM153B",
  "ENSG00000182584" = "ACTL10",
  "ENSG00000189366" = "ALG1L",
  "ENSG00000196312" = "MFSD14C",
  "ENSG00000203690" = "TCP10",
  "ENSG00000205212" = "CCDC144NL",
  "ENSG00000214534" = "ZNF705E",
  "ENSG00000215790" = "SLC35E2A",
  "ENSG00000259511" = "UBE2Q2L",
  "ENSG00000274744" = "ELOA3D",
  "ENSG00000275993" = "SIK1B",
  "ENSG00000276076" = "CRYAA2",
  "ENSG00000276289" = "KCNE1B",
  "ENSG00000278674" = "ELOA3B",
  "ENSG00000280071" = "GATD3B"
)

old_names <- rownames(sce_obj)
new_names <- ifelse(old_names %in% names(manual_map), manual_map[old_names], old_names)

patch_df <- data.frame(
  old_name = old_names[old_names %in% names(manual_map)],
  new_name = unname(manual_map[old_names[old_names %in% names(manual_map)]]),
  row.names = NULL
)

rowData(sce_obj)$original_name_before_manual_map <- rownames(sce_obj)
rownames(sce_obj) <- make.unique(new_names)


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

# Bin parity to distinguish nulliparous donors, well-powered parous donors with 1-3 live births, and the lower-powered 4+ group.
colData(sce_obj)$Para_bin <- factor(ifelse(colData(sce_obj)$Para == 0, "0", ifelse(colData(sce_obj)$Para <= 3, "1-3", "4+")), levels = c("0", "1-3", "4+"))


#### 3. Compute StressScore and create pseudobulk ####

# For the previously determined IER scores (from Marsh 2022, Nat Neuroscience), z-score per cell type in preparation for pseudobulk
# Calculate modified (robust) z-score for IER score per cell type

# Define a function for modified (robust) z-scoring
mod_zscore <- function(x) {
  m  <- median(x, na.rm = TRUE)
  s  <- stats::mad(x, constant = 1.4826, na.rm = TRUE)  # consistent with SD under Normal
  return((x - m) / s)
}

# Add IER z-score to colData
colData(sce_obj)$IER_zscore <- as.data.frame(colData(sce_obj)) %>%
  group_by(cell_type_fine) %>%
  mutate(IER_zscore = mod_zscore(marsh_stress_UCell)) %>%
  pull(IER_zscore)

# Plot IER_zscore per cell_type_fine using ggplot2
df_ier <- as.data.frame(colData(sce_obj))

ggplot(df_ier, aes(x = cell_type_fine, y = IER_zscore, fill = cell_type_fine)) +
  geom_violin(trim = FALSE, scale = "width", color = NA, alpha = 0.7) +
  geom_boxplot(width = 0.2, outlier.shape = NA, fill = "white", color = "black", alpha = 0.8) +
  labs(
    title = "IER z-score per cell type",
    x = "Cell Type (Fine)",
    y = "IER z-score"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# Prep for Pseudobulk

# Create pseudobulk profiles per unique sample/tissue piece and cell type
pb_unfilt <- aggregateToPseudoBulk(
  sce_obj,
  assay = "raw_counts",
  cluster_id = "cell_type_fine",
  sample_id = "unique_sample_ID",
  verbose = TRUE
)
pb_unfilt

# Define a sample-level StressScore as the median cell-type-level IER z-score across available cell types per sample.
stress_score_tbl <- metadata(pb_unfilt)$aggr_means %>%
  group_by(unique_sample_ID) %>%
  summarise(
    StressScore = median(IER_zscore, na.rm = TRUE),
    n_types_used     = sum(!is.na(IER_zscore))
  )

colData(pb_unfilt)$StressScore<- stress_score_tbl$StressScore[match(rownames(colData(pb_unfilt)), stress_score_tbl$unique_sample_ID)]

# Export sample-level StressScore covariate
sample_stressscore_tbl <- stress_score_tbl %>%
  dplyr::rename(
    n_cell_types_used_for_StressScore = n_types_used
  ) %>%
  dplyr::left_join(
    as.data.frame(colData(pb_unfilt)) %>%
      tibble::rownames_to_column("unique_sample_ID") %>%
      dplyr::select(
        unique_sample_ID,
        SenNet_ID,
        Plex
      ) %>%
      dplyr::distinct(),
    by = "unique_sample_ID"
  ) %>%
  dplyr::select(
    unique_sample_ID,
    SenNet_ID,
    Plex,
    StressScore,
    n_cell_types_used_for_StressScore
  ) %>%
  dplyr::arrange(Plex, SenNet_ID, unique_sample_ID)

write.csv(
  sample_stressscore_tbl,
  file = file.path(outs_dir, "Sample_StressScore.csv"),
  row.names = FALSE
)


#### 4. Prepare metadata and covariate checks ####

# Remove incomplete Para-Gravida Ratio and handle NA for Uterine_Pathology
meta_df <- as.data.frame(colData(pb_unfilt)) %>%
  dplyr::select(-Para_Gravida_Ratio) %>%
  mutate(Uterine_Pathology = ifelse(is.na(Uterine_Pathology), "None", as.character(Uterine_Pathology)))

colData(pb_unfilt)$Uterine_Pathology <- ifelse(is.na(colData(pb_unfilt)$Uterine_Pathology), "None",
  as.character(colData(pb_unfilt)$Uterine_Pathology))

# Save meta_df and pb_unfilt to RDS files
saveRDS(meta_df, file = file.path(r_analyses_dir, "donor_pb_meta_df.rds")) # to be used by all other analyses scripts for sample-level StressScore
write.csv(meta_df, file = file.path(r_analyses_dir, "donor_pb_meta_df.csv"), row.names = TRUE)

saveRDS(pb_unfilt, file = file.path(outs_dir, "pb_unfilt.rds"))

## Check correlation between Age and StressScore
# Sample-level StressScore is used in models; donor medians are used only to test association with donor-level Age
sample_stress_df <- meta_df %>%
  as.data.frame() %>%
  tibble::rownames_to_column("unique_sample_ID") %>%
  dplyr::select(unique_sample_ID, SenNet_ID, Age, StressScore) %>%
  dplyr::filter(!is.na(Age), !is.na(StressScore), !is.na(SenNet_ID))

donor_stress_df <- sample_stress_df %>%
  dplyr::group_by(SenNet_ID) %>%
  dplyr::summarise(
    Age = dplyr::first(Age),
    donor_median_StressScore = median(StressScore, na.rm = TRUE),
    n_samples = dplyr::n(),
    .groups = "drop"
  )

spearman_test <- cor.test(
  donor_stress_df$Age,
  donor_stress_df$donor_median_StressScore,
  method = "spearman",
  exact = FALSE
)

rho <- unname(spearman_test$estimate)
pval <- spearman_test$p.value

corr_label <- sprintf(
  "Spearman \u03c1 = %.2f\np = %.2g\nn = %d donors",
  rho, pval, nrow(donor_stress_df)
)

p_stress_age <- ggplot() +
  geom_point(
    data = sample_stress_df,
    aes(x = Age, y = StressScore),
    color = "grey55",
    alpha = 0.35,
    size = 1.5,
    stroke = 0
  ) +
  geom_point(
    data = donor_stress_df,
    aes(x = Age, y = donor_median_StressScore),
    shape = 16,    # solid circle
    color = "black",
    size = 2.5,
    stroke = 0    # no border
  ) +
  annotate(
    "text",
    x = min(donor_stress_df$Age, na.rm = TRUE) + 3,
    y = max(sample_stress_df$StressScore, na.rm = TRUE) - 0.001 * diff(range(sample_stress_df$StressScore, na.rm = TRUE)),
    label = corr_label,
    hjust = 0,
    vjust = 1,
    size = 1.5
  ) +
  labs(
    x = "Donor Age",
    y = "StressScore",
  ) +
  theme_bw(base_size = 12, base_family = "DejaVu Sans") +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.title.x = element_text(margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 5)),
    axis.line = element_line(color = "black", linewidth = 0.25),
    axis.ticks = element_line(color = "black", linewidth = 0.25),
    plot.caption = element_text(size = 5.5, hjust = 0, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

p_stress_age

ggsave(filename = file.path(outs_dir, "StressScore_age_correlation.pdf"), plot = p_stress_age, width = 160, height = 90,
  units = "mm", device = cairo_pdf)

# VIF on fixed-effect design
tmp <- meta_df %>%
  dplyr::select(Age, BMI, Para_bin, StressScore, Plex) %>%
  na.omit()
set.seed(42)
tmp$y_dummy <- rnorm(nrow(tmp))
car::vif(lm(y_dummy ~ Age + BMI + Para_bin + StressScore + Plex, data = tmp))
# No evidence of multicollinearity in fixed-effect design
#                 GVIF Df GVIF^(1/(2*Df))
# Age         1.498466  1        1.224119
# BMI         1.157509  1        1.075876
# Para_bin    1.674591  2        1.137568
# StressScore 1.242309  1        1.114589
# Plex        1.199775  3        1.030821


#### 5. Run dreamlet differential expression ####

# Final primary model for Age and parity DGE analyses.
core_form <- as.formula(paste("~ Age + BMI + Para_bin + StressScore + Plex + (1|SenNet_ID)"))

res.proc.core <- processAssays(sceObj = pb_unfilt, formula = core_form, min.samples = 10, min.cells = 5, min.count = 5, BPPARAM = bp)
res.dl <- dreamlet(res.proc.core, core_form, use.eBayes = TRUE, robust = TRUE, BPPARAM = bp)

# Save the objects for later import
saveRDS(res.proc.core, file = file.path(outs_dir, "res_proc_core.rds"))
saveRDS(res.dl, file = file.path(outs_dir, "res_dl.rds"))

# dreamlet::topTable on the full result reports study-wide FDR. For cell-type-level DGE tables,
# extract topTable per assay to report per-cell-type FDR values.
get_tabs_merged <- function(res.dl, coef) {
  assays_used <- topTable(res.dl, coef = coef, number = Inf)$assay %>% unique()
  tabs_merged <- do.call(rbind, lapply(assays_used, function(a) {
    df <- dreamlet::topTable(assay(res.dl, a), coef = coef, number = Inf)
    df <- cbind(Celltype = a, Gene = rownames(df), df)
    rownames(df) <- NULL
    return(df)
  }))
  return(tabs_merged)
}

tabs_merged_age <- get_tabs_merged(res.dl, "Age")
tabs_merged_para <- get_tabs_merged(res.dl, "Para_bin1-3") # Primary parity contrast: 1-3 vs 0 live births

# Combine both tables and add a 'contrast' column with the coef name
tabs_merged_age$contrast <- "Age"
tabs_merged_para$contrast <- "Para_bin1-3"
tabs_merged_combined <- rbind(tabs_merged_age, tabs_merged_para)
tabs_merged_combined <- tabs_merged_combined[, c("Celltype", "contrast", setdiff(names(tabs_merged_combined), c("Celltype", "contrast")))]

# Save the merged tables to CSV files
write.csv(tabs_merged_combined, file = file.path(outs_dir, "dreamlet_dge_results.csv"), row.names = FALSE)


#### 6. Determine cell-type specific perturbation descriptively for ranking via root-mean-square (RMS) of moderated t-statistic ####

# Function to plot the top cell types for Age and Parity for Main Fig to highlight the most perturbed cell types
plot_top_celltype_perturbation <- function(
  tabs_age,
  tabs_para,
  stat_col = "t",
  n_top_genes = 500,
  top_n_celltypes = 6,
  base_size = 6
) {
  minmax01 <- function(x) {
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) == 0) return(rep(0.5, length(x)))
    (x - rng[1]) / diff(rng)
  }

  score_by_celltype_top <- function(df, source_name) {
    df %>%
      dplyr::group_by(Celltype) %>%
      dplyr::arrange(dplyr::desc(abs(.data[[stat_col]])), .by_group = TRUE) %>%
      dplyr::slice_head(n = n_top_genes) %>%
      dplyr::summarise(
        perturb_score = sqrt(mean((.data[[stat_col]])^2, na.rm = TRUE)),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        rel_score = minmax01(perturb_score),
        source = source_name
      ) %>%
      dplyr::arrange(dplyr::desc(perturb_score))
  }

  score_age  <- score_by_celltype_top(tabs_age, "Age")
  score_para <- score_by_celltype_top(tabs_para, "Parity")

  union_celltypes <- union(
    score_age %>% dplyr::slice_max(perturb_score, n = top_n_celltypes, with_ties = FALSE) %>% dplyr::pull(Celltype),
    score_para %>% dplyr::slice_max(perturb_score, n = top_n_celltypes, with_ties = FALSE) %>% dplyr::pull(Celltype)
  )

  plot_df <- dplyr::bind_rows(score_age, score_para) %>%
    dplyr::filter(Celltype %in% union_celltypes)

  cell_order <- plot_df %>%
    dplyr::filter(source == "Age") %>%
    dplyr::arrange(rel_score) %>%
    dplyr::pull(Celltype)

  plot_df <- plot_df %>%
    dplyr::mutate(Celltype = factor(Celltype, levels = cell_order))

  theme_perturb <- function() {
    theme_classic(base_size = base_size) +
      theme(
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = base_size, color = "black"),
        axis.text.x = element_text(size = base_size, color = "black"),
        axis.title.x = element_text(size = base_size, margin = margin(t = 4)),
        plot.title = element_text(face = "bold", hjust = 0.5, size = base_size + 1),
        plot.margin = margin(2, 6, 2, 2),
        panel.border = element_rect(color = "grey30", fill = NA, linewidth = 0.2),
        axis.line = element_line(linewidth = 0)
      )
  }

  p_age <- ggplot(
    dplyr::filter(plot_df, source == "Age"),
    aes(x = rel_score, y = Celltype)
  ) +
    geom_segment(
      aes(x = 0, xend = rel_score, yend = Celltype),
      linewidth = 0.5, lineend = "round"
    ) +
    geom_point(size = 2) +
    coord_cartesian(xlim = c(0, 1.02)) +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
    labs(title = "Age", x = "Relative Transcriptional Perturbation Score") +
    theme_perturb()

  p_para <- ggplot(
    dplyr::filter(plot_df, source == "Parity"),
    aes(x = rel_score, y = Celltype)
  ) +
    geom_segment(
      aes(x = 0, xend = rel_score, yend = Celltype),
      linewidth = 0.5, lineend = "round"
    ) +
    geom_point(size = 2) +
    coord_cartesian(xlim = c(0, 1.02)) +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
    labs(title = "Parity", x = "Relative Transcriptional Perturbation Score") +
    theme_perturb() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  p_age | p_para
}

# Plot top 6 cell types for Age and Parity for Main Fig to highlight the most perturbed cell types
perturb_plot_top6 <- plot_top_celltype_perturbation(
  tabs_age = tabs_merged_age,
  tabs_para = tabs_merged_para,
  stat_col = "t",
  n_top_genes = 500,
  top_n_celltypes = 6,
  base_size = 6
)

# Plot all cell types for Age and Parity for Extended Data Fig to highlight the most perturbed cell types
perturb_plot_all <- plot_top_celltype_perturbation(
  tabs_age = tabs_merged_age,
  tabs_para = tabs_merged_para,
  stat_col = "t",
  n_top_genes = 500,
  top_n_celltypes = 100,
  base_size = 6
)

ggsave(filename = paste0(outs_dir, "/perturbation_scores_age_parity_top6.pdf"), plot = perturb_plot_top6, width = 180,
  height = 40, units = "mm", device = cairo_pdf)
ggsave(filename = paste0(outs_dir, "/perturbation_scores_age_parity_all.pdf"), plot = perturb_plot_all, width = 180,
  height = 80, units = "mm", device = cairo_pdf)


#### 6b. Sensitivity of perturbation score rankings to number of genes retained ####

# Compute RMS perturbation score for one contrast across multiple top-N settings
compute_perturbation_scores_multiN <- function(
  df,
  contrast_name,
  stat_col = "t",
  n_top_genes_vec = c(250, 500, 1000)
) {
  minmax01 <- function(x) {
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) == 0) return(rep(0.5, length(x)))
    (x - rng[1]) / diff(rng)
  }

  purrr::map_dfr(n_top_genes_vec, function(n_top) {
    df %>%
      dplyr::group_by(Celltype) %>%
      dplyr::arrange(dplyr::desc(abs(.data[[stat_col]])), .by_group = TRUE) %>%
      dplyr::slice_head(n = n_top) %>%
      dplyr::summarise(
        perturb_score = sqrt(mean((.data[[stat_col]])^2, na.rm = TRUE)),
        n_genes_used = dplyr::n(),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        rel_score = minmax01(perturb_score),
        rank = dplyr::min_rank(dplyr::desc(rel_score)),
        contrast = contrast_name,
        n_top_genes = n_top
      )
  })
}

# Compute sensitivity table
perturb_sensitivity <- dplyr::bind_rows(
  compute_perturbation_scores_multiN(
    tabs_merged_age,
    contrast_name = "Age",
    stat_col = "t",
    n_top_genes_vec = c(250, 500, 1000)
  ),
  compute_perturbation_scores_multiN(
    tabs_merged_para,
    contrast_name = "Parity",
    stat_col = "t",
    n_top_genes_vec = c(250, 500, 1000)
  )
)

# Reference: reported top_n = 500
perturb_ref_500 <- perturb_sensitivity %>%
  dplyr::filter(n_top_genes == 500) %>%
  dplyr::select(
    Celltype,
    contrast,
    rel_score_500 = rel_score,
    rank_500 = rank,
    perturb_score_500 = perturb_score
  )

# Comparisons: 250 vs 500 and 1000 vs 500
perturb_compare <- perturb_sensitivity %>%
  dplyr::filter(n_top_genes %in% c(250, 1000)) %>%
  dplyr::left_join(
    perturb_ref_500,
    by = c("Celltype", "contrast")
  ) %>%
  dplyr::mutate(
    comparison = paste0("Top ", n_top_genes, " vs top 500"),
    comparison = factor(
      comparison,
      levels = c("Top 250 vs top 500", "Top 1000 vs top 500")
    ),
    contrast = factor(contrast, levels = c("Age", "Parity"))
  )

# Spearman correlations per panel
cor_df <- perturb_compare %>%
  dplyr::group_by(contrast, comparison) %>%
  dplyr::summarise(
    spearman_rho_score = cor(rel_score_500, rel_score, method = "spearman", use = "complete.obs"),
    spearman_rho_rank = cor(rank_500, rank, method = "spearman", use = "complete.obs"),
    n_cell_types = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    label = sprintf(
      "Spearman ρ = %.2f",
      spearman_rho_score
    ),
    x = 0.03,
    y = 0.97
  )

# Plot: score stability, with cell-type labels
p_perturb_sensitivity <- ggplot(
  perturb_compare,
  aes(x = rel_score_500, y = rel_score, label = Celltype)
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linewidth = 0.25,
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_point(size = 1.2, alpha = 0.9) +
  ggrepel::geom_text_repel(
    size = 1.8,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.linewidth = 0.15,
    box.padding = 0.15,
    point.padding = 0.1
  ) +
  geom_text(
    data = cor_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 2.0
  ) +
  facet_grid(contrast ~ comparison) +
  coord_cartesian(xlim = c(0, 1.02), ylim = c(0, 1.02), clip = "off") +
  scale_x_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1.0),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1.0),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    x = "Relative perturbation score using top 500 genes",
    y = "Relative perturbation score using alternate top-N genes"
  ) +
  theme_bw(base_size = 5) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(color = "black", size = 6),
    axis.text = element_text(color = "black", size = 6),
    axis.title = element_text(color = "black", size = 6),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.25),
    plot.margin = margin(4, 4, 4, 4)
  )

p_perturb_sensitivity

ggsave(filename = paste0(outs_dir, "/perturbation_score_topN_sensitivity_spearman.pdf"),  plot = p_perturb_sensitivity,
  width = 180, height = 110, units = "mm", device = cairo_pdf)

# Save values for Supplementary Table
write.csv(perturb_sensitivity, file = file.path(outs_dir, "RMS_topN_sensitivity_scores.csv"), row.names = FALSE)


#### 7. Run Hallmark pathway analysis with zenith ####

# Import GeneSetCollection objects
unique(msigdbr::msigdbr(species = "human", collection = "H")$db_version) #  "2025.1.Hs"
msigdb_h <- get_MSigDB(cat = 'H', to = 'SYMBOL', organism = 'HS') # Hallmarks

# Function to extract gene set metadata from a GeneSetCollection object
extract_gsc_info <- function(gsc) {
  data.frame(
    name = sapply(gsc@.Data, GSEABase::setName),
    size = sapply(gsc@.Data, function(gs) length(GSEABase::geneIds(gs))),
    category = sapply(gsc@.Data, function(gs) GSEABase::bcCategory(GSEABase::collectionType(gs))),
    subcategory = sapply(gsc@.Data, function(gs) GSEABase::bcSubCategory(GSEABase::collectionType(gs)))
  )
}

filter_msigdb <- function(msigdb_gsc, subcategory = NULL, min_size = 15, max_size = 500) {
  gsc_info <- extract_gsc_info(msigdb_gsc)
  cat(sprintf("Number of gene sets before filtering: %d\n", nrow(gsc_info)))

  filtered <- gsc_info %>%
    filter(size >= min_size, size <= max_size)
  if (!is.null(subcategory)) {
    filtered <- filtered %>% filter(subcategory == !!subcategory)
  }

  filtered_names <- pull(filtered, name)
  cat(sprintf("Number of gene sets after filtering: %d\n", length(filtered_names)))

  msigdb_gsc[gsc_info$name %in% filtered_names]
}

gsc_filt <- filter_msigdb(msigdb_h, min_size = 15, max_size = 500) # Hallmark gene sets were retained after standard size filtering.

# Run zenith_gsa per unique cell type (assay) and per desired coef, then concatenate results
genesets <- gsc_filt
cell_types <- assayNames(res.dl)
coefs <- c("Age", "Para_bin1-3")

# Run zenith_gsa for all combinations and concatenate results
res_zenith_all <- do.call(
  rbind,
  lapply(cell_types, function(cell_type) {
    lapply(coefs, function(coef) {
      message(sprintf("Running zenith_gsa for assay/cell type: %s, coef: %s", cell_type, coef))
      # Subset res.dl to just this cell_type (assay)
      res.dl.sub <- res.dl[assayNames(res.dl) %in% cell_type]
      res <- zenith_gsa(res.dl.sub, genesets, coefs = coef, use.ranks = FALSE, inter.gene.cor = 0.01, progress = TRUE)
      if (!"assay" %in% colnames(res)) res$assay <- cell_type
      res$coef <- coef
      res$t <- res$delta / res$se
      res
    })
  }) %>% unlist(recursive = FALSE)
)

# Geneset mapping (pretty labels) for all 50 Hallmark sets
geneset_map <- tribble(
  ~Geneset_clean, ~Geneset_pretty,
  "ALLOGRAFT_REJECTION", "Allograft rejection",
  "ANDROGEN_RESPONSE", "Androgen response",
  "ANGIOGENESIS", "Angiogenesis",
  "APICAL_JUNCTION", "Apical junction",
  "APICAL_SURFACE", "Apical surface",
  "APOPTOSIS", "Apoptosis",
  "ADIPOGENESIS", "Adipogenesis",
  "BILE_ACID_METABOLISM", "Bile acid metabolism",
  "CHOLESTEROL_HOMEOSTASIS", "Cholesterol homeostasis",
  "COAGULATION", "Coagulation",
  "COMPLEMENT", "Complement",
  "DNA_REPAIR", "DNA repair",
  "E2F_TARGETS", "E2F targets",
  "EPITHELIAL_MESENCHYMAL_TRANSITION", "Epithelial–mesenchymal transition",
  "ESTROGEN_RESPONSE_EARLY", "Estrogen response (early)",
  "ESTROGEN_RESPONSE_LATE", "Estrogen response (late)",
  "FATTY_ACID_METABOLISM", "Fatty acid metabolism",
  "G2M_CHECKPOINT", "G2/M checkpoint",
  "GLYCOLYSIS", "Glycolysis",
  "HEDGEHOG_SIGNALING", "Hedgehog signaling",
  "HEME_METABOLISM", "Heme metabolism",
  "HYPOXIA", "Hypoxia",
  "IL2_STAT5_SIGNALING", "IL-2/STAT5 signaling",
  "IL6_JAK_STAT3_SIGNALING", "IL-6/JAK/STAT3 signaling",
  "INFLAMMATORY_RESPONSE", "Inflammatory response",
  "INTERFERON_ALPHA_RESPONSE", "Interferon-α response",
  "INTERFERON_GAMMA_RESPONSE", "Interferon-γ response",
  "KRAS_SIGNALING_DN", "KRAS signaling (down)",
  "KRAS_SIGNALING_UP", "KRAS signaling (up)",
  "MITOTIC_SPINDLE", "Mitotic spindle",
  "MTORC1_SIGNALING", "mTORC1 signaling",
  "MYC_TARGETS_V1", "MYC targets v1",
  "MYC_TARGETS_V2", "MYC targets v2",
  "MYOGENESIS", "Myogenesis",
  "NOTCH_SIGNALING", "Notch signaling",
  "OXIDATIVE_PHOSPHORYLATION", "Oxidative phosphorylation",
  "PANCREAS_BETA_CELLS", "Pancreas beta cells",
  "P53_PATHWAY", "p53 pathway",
  "PEROXISOME", "Peroxisome",
  "PI3K_AKT_MTOR_SIGNALING", "PI3K/AKT/mTOR signaling",
  "PROTEIN_SECRETION", "Protein secretion",
  "REACTIVE_OXYGEN_SPECIES_PATHWAY", "Reactive oxygen species pathway",
  "SPERMATOGENESIS", "Spermatogenesis",
  "TGF_BETA_SIGNALING", "TGF-β signaling",
  "TNFA_SIGNALING_VIA_NFKB", "TNFα signaling via NF-κB",
  "UNFOLDED_PROTEIN_RESPONSE", "Unfolded protein response",
  "UV_RESPONSE_DN", "UV response (down)",
  "UV_RESPONSE_UP", "UV response (up)",
  "WNT_BETA_CATENIN_SIGNALING", "Wnt/β-catenin signaling",
  "XENOBIOTIC_METABOLISM", "Xenobiotic metabolism"
)

res_zenith_all_clean <- res_zenith_all %>%
  mutate(
    Geneset_clean = str_remove(Geneset, "^HALLMARK_")
  ) %>%
  left_join(geneset_map, by = c("Geneset_clean" = "Geneset_clean")) %>%
  mutate(
    Geneset = factor(Geneset_pretty, levels = geneset_map$Geneset_pretty)
  )

# Save res_zenith_all
saveRDS(res_zenith_all_clean, file = file.path(outs_dir, "res_zenith_hallmarks_all.rds"))
#res_zenith_all_clean <- readRDS(file.path(outs_dir, "res_zenith_hallmarks_all.rds"))
res_zenith_all_clean %>% head(15)

# Clean ASCII labels for pathways
clean_ascii_labels <- function(x) {
  x %>%
    stringr::str_replace_all("α", "alpha") %>%
    stringr::str_replace_all("β", "beta") %>%
    stringr::str_replace_all("κ", "kappa") %>%
    stringr::str_replace_all("γ", "gamma") %>%
    stringr::str_replace_all("–", "-") %>%
    stringr::str_replace_all("—", "-")
}

# Select main contrasts and clean ASCII labels for pathways
res_zenith_supp_table <- res_zenith_all_clean %>%
  mutate(pathway = clean_ascii_labels(Geneset)) %>%
  filter(coef %in% c("Age", "Para_bin1-3")) %>%
  dplyr::select(
    cell_type = assay,
    contrast = coef,
    pathway = pathway,
    n_genes = NGenes,
    inter_gene_correlation = Correlation,
    delta = delta,
    se = se,
    t = t,
    p_value = PValue,
    fdr = FDR,
    direction = Direction
  )

# Save to .csv for supplementary table
write.csv(res_zenith_supp_table, file = file.path(outs_dir, "res_zenith_main_contrasts.csv"), row.names = FALSE)


#### 8. Generate pathway dotplots ####

## Identify pathways significantly associated with Age and Parity across cell types
# Prioritize pathways that are FDR <= 0.05, t >= 3, and are present in many cell types

# Summarize top recurrent pathways for a given coefficient
get_top_recurrent_pathways <- function(df, coef_name, fdr_max = 0.05, t_min = 3, top_n = 10) {
  df %>%
    filter(coef == coef_name, FDR <= fdr_max, abs(t) >= t_min) %>%
    group_by(Geneset) %>%
    summarize(
      n_cell_types = n_distinct(assay),
      mean_t = mean(t, na.rm = TRUE),
      min_FDR = min(FDR, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n_cell_types), desc(abs(mean_t)), min_FDR) %>%
    slice_head(n = top_n)
}

# Get top recurrent pathways for multiple coefficients
get_recurrent_pathways_multi <- function(df, coef_names, fdr_max = 0.05, t_min = 3, top_n = 10) {
  lapply(setNames(coef_names, coef_names), function(x) {
    get_top_recurrent_pathways(
      df = df,
      coef_name = x,
      fdr_max = fdr_max,
      t_min = t_min,
      top_n = top_n
    )
  })
}

res_list <- get_recurrent_pathways_multi(res_zenith_all_clean, coef_names = c("Age", "Para_bin1-3"))

top_age_pathways  <- res_list[["Age"]]
top_para_pathways <- res_list[["Para_bin1-3"]]
common_pathways   <- intersect(top_age_pathways$Geneset, top_para_pathways$Geneset)

# Function to plot gene set enrichment results pathways as dotplots
plot_pathways <- function(
  df,
  pathways,
  coefs_plot = c("Age", "Para_bin1-3"),
  celltypes_plot = c("Homeostatic Fibroblast"),
  coef_labels = c("Age" = "↑Age", "Para_bin1-3" = "↑Parity"),
  assay_labels = NULL, #c("HOXD-patterned Myofibroblast" = "Myofibro")
  base_size = 6,
  sig_fdr_for_stroke = 0.05,
  dot_range = c(0.5, 4),
  panel_spacing = 1
) {

  # display labels for facets
  facet_labels <- setNames(celltypes_plot, celltypes_plot)
  facet_labels[names(assay_labels)] <- assay_labels[names(assay_labels)]

  assay_labeller <- function(x) {
    label_wrap_gen(width = 10)(unname(facet_labels[x]))
  }

  keep_df_long <- df %>%
    filter(
      Geneset %in% pathways,
      coef %in% coefs_plot,
      assay %in% celltypes_plot
    ) %>%
    mutate(
      Geneset = factor(Geneset, levels = pathways),
      coef_plot = dplyr::recode(coef, !!!coef_labels),
      assay = factor(assay, levels = celltypes_plot),
      assay_coef = factor(
        paste(assay, coef, sep = "_"),
        levels = as.vector(outer(celltypes_plot, coefs_plot, paste, sep = "_"))
      )
    )

  ggplot(keep_df_long, aes(x = coef_plot, y = Geneset, fill = t, size = -log10(FDR))) +
    geom_point(
      aes(stroke = ifelse(FDR <= sig_fdr_for_stroke, 0.5, 0.1)),
      shape = 21,
      alpha = 1
    ) +
    scale_fill_gradient2(
      name = "Enrichment\nScore",
      low = "blue", mid = "white", high = "red",
      guide = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        label.position = "right",
        barwidth = unit(0.3, "cm"),
        barheight = unit(1.8, "cm"),
        order = 1
      )
    ) +
    scale_size(
      name = expression(-log[10]~"FDR"),
      range = dot_range,
      guide = guide_legend(
        title.position = "top",
        title.hjust = 0.5,
        order = 2,
        keyheight = unit(0.4, "cm")
      )
    ) +
    labs(x = NULL, y = NULL) +
    facet_grid(
      . ~ assay,
      scales = "free_x",
      space = "free_x",
      switch = "x",
      labeller = labeller(assay = assay_labeller)
    ) +
    theme_bw(base_size = base_size, base_family = "DejaVu Sans") +
    theme(
      strip.placement = "outside",
      strip.background = element_blank(),
      axis.text.x = element_text(
        angle = 0, hjust = 0.5, vjust = 1,
        colour = "black", size = base_size,
        margin = margin(t = 4)
      ),
      axis.text.y = element_text(colour = "black", size = base_size),
      strip.text.x = element_text(colour = "black", size = base_size - 1),
      strip.text.y = element_text(colour = "black"),
      legend.title = element_text(colour = "black", hjust = 0.5, size = base_size),
      legend.text = element_text(colour = "black", size = base_size),
      legend.box.just = "center",
      legend.position = "right",
      legend.box = "vertical",
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(panel_spacing, "lines"),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

## Plotting the top pathways for Age and Parity
coefs_plot <- c("Age", "Para_bin1-3")

celltype_highlight <- c("Homeostatic Fibroblast", "IER Stromal", "HOXD-patterned Myofibroblast", "CapEnd", "LympEnd", "ArtSMC")
coef_labels <- c("Age" = "↑Age", "Para_bin1-3" = "↑Parity")

p_main_common_pathways <- plot_pathways(
  df = res_zenith_all_clean,
  pathways = common_pathways,
  celltypes_plot = celltype_highlight,
  coefs_plot = coefs_plot,
  coef_labels = coef_labels,
  base_size = 6
)
p_main_common_pathways

ggsave(file.path(outs_dir, "zenith_hallmarks_age_para_common_assays.pdf"), p_main_common_pathways,
  width = 180, height = 60, units = "mm", device = cairo_pdf)

# Plot all cell types & significant pathways for extended
all_celltypes <- res_zenith_all_clean$assay %>% unique()

# Get all significant pathways (FDR <0.05) across all cell types
all_sig_pathways <- unique(res_zenith_all_clean$Geneset[res_zenith_all_clean$FDR < 0.05])
# Count unique cell types for each pathway with FDR < 0.05
pathway_celltype_counts <- res_zenith_all_clean[res_zenith_all_clean$FDR < 0.05, ] %>%
  group_by(Geneset) %>%
  summarise(n_celltypes = n_distinct(assay)) %>%
  arrange(n_celltypes)  # sort from least to most significant cell types

# Sort all_sig_pathways by ascending number of unique cell types (FDR < 0.05)
all_sig_pathways <- factor(pathway_celltype_counts$Geneset, levels = pathway_celltype_counts$Geneset)
all_sig_pathways

# For readability and breadth, get top 3 pathways (by |t|) per cell type for either Age or Parity and create a non-redundant list
top_n <- 3
coefs_interest <- c("Age", "Para_bin1-3")

# Get top pathways per cell type for each coef, as before
top_pathways_per_celltype <- res_zenith_all_clean %>%
  filter(coef %in% coefs_interest, FDR < 0.05) %>%
  group_by(assay, coef) %>%
  arrange(desc(abs(t))) %>%
  slice_head(n = top_n) %>%
  ungroup()

# Count, for each pathway, in how many unique cell types (assay) is it significant (FDR < 0.05)
pathway_celltype_count <- top_pathways_per_celltype %>%
  group_by(Geneset) %>%
  summarise(n_celltypes = n_distinct(assay)) %>%
  arrange(n_celltypes)

# Sort unique list of pathways by most cell type presence (across top n and FDR < 0.05)
top_pathways_nonredundant <- factor(pathway_celltype_count$Geneset, levels = pathway_celltype_count$Geneset)

p_main_all_pathways_1 <- plot_pathways(
  df = res_zenith_all_clean,
  pathways = top_pathways_nonredundant,
  celltypes_plot = head(all_celltypes, floor(length(all_celltypes) / 2)),
  coefs_plot = coefs_plot,
  dot_range = c(0.2, 2),
  panel_spacing = 0.2,
  coef_labels = c("Age" = "↑Age", "Para_bin1-3" = "↑Para"),
  base_size = 5
)

p_main_all_pathways_2 <- plot_pathways(
  df = res_zenith_all_clean,
  pathways = top_pathways_nonredundant,
  celltypes_plot = tail(all_celltypes, floor(length(all_celltypes)/2)),
  coefs_plot = coefs_plot,
  dot_range = c(0.2, 2),
  panel_spacing = 0.2,
  coef_labels = c("Age" = "↑Age", "Para_bin1-3" = "↑Para"),
  base_size = 5
)

ggsave(file.path(outs_dir, "zenith_hallmarks_age_para_all_assays_1.pdf"), p_main_all_pathways_1,
  width = 180, height = 60, units = "mm", device = cairo_pdf)

ggsave(file.path(outs_dir, "zenith_hallmarks_age_para_all_assays_2.pdf"), p_main_all_pathways_2,
  width = 180, height = 60, units = "mm", device = cairo_pdf)




#### 9. Extract pathway driver genes and heatmaps ####

## Identify genes that are drivers for enrichment of a pathway across cell types
consensus_pathway_drivers <- function(
  res.dl,
  res_zenith_all,
  gsc_obj,
  geneset,
  coef_id = "Age",
  outs_dir = ".",
  save_csv = c("all", "concise", "none"),
  top_n_mainfig = 30,
  min_frac_present = 0.75,
  min_frac_pos = 0.80,
  min_mean_t = 1.0
) {
  save_csv <- match.arg(save_csv)

  # 1) Get genes in pathway
  genes_in_set <- geneIds(gsc_obj[[geneset]])

  # 2) Pick assays/cell types where pathway is positively enriched
  enriched_assays <- res_zenith_all %>%
    dplyr::filter(
      Geneset == !!geneset,
      FDR <= 0.05,
      delta > 0
    ) %>%
    dplyr::pull(assay) %>%
    unique()

  # Edge case: return if no assay is enriched
  if (length(enriched_assays) == 0) {
    warning(sprintf("No assays enriched for pathway '%s'. No files written.", geneset))
    return(invisible(NULL))
  }

  # 3) Build gene x assay matrix of t-statistics and logFC
  tstat_mat <- matrix(
    NA_real_,
    nrow = length(genes_in_set),
    ncol = length(enriched_assays),
    dimnames = list(genes_in_set, enriched_assays)
  )

  logfc_mat <- matrix(
    NA_real_,
    nrow = length(genes_in_set),
    ncol = length(enriched_assays),
    dimnames = list(genes_in_set, enriched_assays)
  )

  for (a in enriched_assays) {
    fit <- tryCatch(assay(res.dl, a), error = function(e) NULL)
    if (is.null(fit)) next

    tt <- tryCatch(
      topTable(fit, coef = coef_id, number = Inf, sort.by = "none") %>%
        tibble::rownames_to_column("gene") %>%
        dplyr::select(gene, logFC, t),
      error = function(e) NULL
    )
    if (is.null(tt)) next

    tt_sub <- tt %>%
      dplyr::filter(gene %in% genes_in_set)

    if (nrow(tt_sub) == 0) next

    tstat_mat[tt_sub$gene, a] <- tt_sub$t
    logfc_mat[tt_sub$gene, a] <- tt_sub$logFC
  }

  # Drop genes missing in all selected assays
  keep_gene <- rowSums(!is.na(tstat_mat)) > 0
  tstat_mat <- tstat_mat[keep_gene, , drop = FALSE]
  logfc_mat <- logfc_mat[keep_gene, , drop = FALSE]

  # 4) Rank genes within each assay (-t so higher t = better; lower rank = better)
  rank_mat <- apply(tstat_mat, 2, function(x) {
    if (all(is.na(x))) {
      rep(NA_real_, length(x))
    } else {
      rank(-x, ties.method = "average", na.last = "keep")
    }
  })
  # Restore matrix structure for n = 1
  if (is.null(dim(rank_mat))) {
    rank_mat <- matrix(rank_mat, ncol = 1, dimnames = dimnames(tstat_mat))
  }

  # 5) Summarize consensus per gene
  n_assays <- ncol(tstat_mat)

  gene_summary <- tibble(
    gene = rownames(tstat_mat),
    n_present = rowSums(!is.na(tstat_mat)),
    n_pos = rowSums(tstat_mat > 0, na.rm = TRUE),
    n_neg = rowSums(tstat_mat < 0, na.rm = TRUE),
    frac_present = rowMeans(!is.na(tstat_mat)),
    frac_pos = rowMeans(tstat_mat > 0, na.rm = TRUE),
    frac_neg = rowMeans(tstat_mat < 0, na.rm = TRUE),
    mean_t = rowMeans(tstat_mat, na.rm = TRUE),
    median_t = apply(tstat_mat, 1, median, na.rm = TRUE),
    min_t = apply(tstat_mat, 1, min, na.rm = TRUE),
    max_t = apply(tstat_mat, 1, max, na.rm = TRUE),
    mean_logFC = rowMeans(logfc_mat, na.rm = TRUE),
    median_logFC = apply(logfc_mat, 1, median, na.rm = TRUE),
    mean_rank = rowMeans(rank_mat, na.rm = TRUE),
    median_rank = apply(rank_mat, 1, median, na.rm = TRUE)
  ) %>%
    dplyr::mutate(
      consensus_score = frac_pos * mean_t
    ) %>%
    dplyr::arrange(desc(consensus_score), mean_rank)

  # 6) Define consensus drivers using provided thresholds
  consensus_drivers <- gene_summary %>%
    dplyr::filter(
      frac_present >= min_frac_present,
      frac_pos >= min_frac_pos,
      mean_t >= min_mean_t
    ) %>%
    dplyr::arrange(mean_rank, desc(mean_t), desc(frac_pos))

  consensus_drivers_top <- consensus_drivers %>%
    dplyr::slice_head(n = top_n_mainfig)

  # 7) Heatmap-ready matrices for top consensus genes
  top_genes <- consensus_drivers_top$gene
  tstat_top_mat <- tstat_mat[top_genes, , drop = FALSE]
  logfc_top_mat <- logfc_mat[top_genes, , drop = FALSE]

  # Order assays by pathway enrichment magnitude if available
  assay_order <- res_zenith_all %>%
    dplyr::filter(Geneset == !!geneset, coef == !!coef_id, assay %in% colnames(tstat_top_mat)) %>%
    dplyr::arrange(desc(delta)) %>%
    dplyr::pull(assay)

  assay_order <- assay_order[assay_order %in% colnames(tstat_top_mat)]
  if (length(assay_order) > 0) {
    tstat_top_mat <- tstat_top_mat[, assay_order, drop = FALSE]
    logfc_top_mat <- logfc_top_mat[, assay_order, drop = FALSE]
  }

  # Order genes by consensus ranking
  gene_order <- consensus_drivers_top$gene
  tstat_top_mat <- tstat_top_mat[gene_order, , drop = FALSE]
  logfc_top_mat <- logfc_top_mat[gene_order, , drop = FALSE]

  # 8) Write outputs
  if (save_csv != "none") {
    # Always write tstat_matrix & consensus_drivers_filtered
    write.csv(gene_summary, file = file.path(outs_dir, paste0("consensus_driver_summary_", geneset, "_", coef_id, ".csv")), row.names = FALSE)
    write.csv(tibble::rownames_to_column(as.data.frame(tstat_mat), "gene"), file = file.path(outs_dir, paste0("tstat_matrix_", geneset, "_", coef_id, ".csv")), row.names = FALSE)
  }

  if (save_csv == "all") {
    write.csv(tibble::rownames_to_column(as.data.frame(logfc_mat), "gene"), file = file.path(outs_dir, paste0("logfc_matrix_", geneset, "_", coef_id, ".csv")), row.names = FALSE)
    write.csv(consensus_drivers, file = file.path(outs_dir, paste0("consensus_drivers_filtered_", geneset, "_", coef_id, ".csv")), row.names = FALSE)
    write.csv(tibble::rownames_to_column(as.data.frame(tstat_top_mat), "gene"), file = file.path(outs_dir, paste0("tstat_matrix_topConsensus_", geneset, "_", coef_id, ".csv")), row.names = FALSE)
    write.csv(tibble::rownames_to_column(as.data.frame(logfc_top_mat), "gene"), file = file.path(outs_dir, paste0("logfc_matrix_topConsensus_", geneset, "_", coef_id, ".csv")), row.names = FALSE)
  }
  # else if 'none', write no tables

  invisible(list(
    tstat_mat = tstat_mat,
    logfc_mat = logfc_mat,
    rank_mat = rank_mat,
    gene_summary = gene_summary,
    consensus_drivers = consensus_drivers,
    consensus_drivers_top = consensus_drivers_top,
    tstat_top_mat = tstat_top_mat,
    logfc_top_mat = logfc_top_mat,
    enriched_assays = enriched_assays
  ))
}

# Function to make a heatmap of the t-statistics for the top consensus genes
make_tstat_heatmap <- function(
  tstat_top_mat,
  show_axis_text_x = FALSE,
  major_vlines = NULL,
  major_hlines = NULL,
  celltype_order = NULL,
  gene_order = NULL,
  base_size = 5,
  font_family = "DejaVu Sans",
  cbar_L = NULL
) {

  plot_df <- as.data.frame(tstat_top_mat) %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::pivot_longer(-gene, names_to = "cell_type", values_to = "t_stat")

  # If not provided, use the current order in matrix as default
  if (is.null(celltype_order)) {
    celltype_order <- colnames(tstat_top_mat)
  }

  plot_df$gene <- factor(plot_df$gene, levels = rev(gene_order))
  plot_df$cell_type <- factor(plot_df$cell_type, levels = celltype_order)

  L <- cbar_L

  axis_text_x_theme <- if (show_axis_text_x) {
    element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black")
  } else {
    element_blank()
  }

  p_heat <- ggplot(plot_df, aes(cell_type, gene, fill = t_stat)) +
    geom_tile(width = 1, height = 1, color = NA) +
    geom_vline(
      xintercept = major_vlines,
      linewidth = 0.1,
      color = "black"
    ) +
    geom_hline(
      yintercept = major_hlines,
      linewidth = 0.1,
      color = "black"
    ) +
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      limits = c(-L, L),
      oob = scales::squish,
      na.value = "grey90",
      name = "t-stat",
      guide = guide_colorbar(
        barwidth = unit(0.3, "cm"),
        barheight = unit(2.5, "cm"),
        label.theme = element_text(color = "black"),
        title.theme = element_text(color = "black")
      )) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = base_size, base_family = font_family) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = axis_text_x_theme,
      axis.text.y = element_text(size = base_size, color = "black")
    ) +
    labs(x = NULL, y = NULL)

  return(p_heat)
}

make_tstat_top_mat <- function(res.dl, coef_id, genes, cell_types) {
  mat <- matrix(NA_real_, nrow = length(genes), ncol = length(cell_types),
                dimnames = list(genes, cell_types))
  for (a in cell_types) {
    fit <- tryCatch(assay(res.dl, a), error = function(e) NULL)
    if (is.null(fit)) next
    tt <- tryCatch(
      topTable(fit, coef = coef_id, number = Inf, sort.by = "none") %>%
        tibble::rownames_to_column("gene") %>%
        dplyr::select(gene, t),
      error = function(e) NULL
    )
    if (is.null(tt)) next
    vals <- tt$t[match(genes, tt$gene)]
    mat[, a] <- vals
  }
  mat
}

## TNFA signaling via NF-κB pathway ##
# Determine & plot consensus gene drivers for TNFA signaling via NF-κB pathway
age_tnfa_consensus_pathway_drivers <- consensus_pathway_drivers(res.dl, res_zenith_all, gsc_filt, "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "Age",
  outs_dir, save_csv = "none")

para_tstat_top_mat_tnfa <- make_tstat_top_mat(
  res.dl, "Para_bin1-3", age_tnfa_consensus_pathway_drivers$consensus_drivers_top$gene,
  colnames(age_tnfa_consensus_pathway_drivers$tstat_top_mat)
)

cbar_L_tnfa <- ceiling(max(
  max(abs(para_tstat_top_mat_tnfa), na.rm = TRUE),
  max(abs(age_tnfa_consensus_pathway_drivers$tstat_top_mat), na.rm = TRUE)
))

tnfa_celltype_order <- c(
  "Homeostatic Fibroblast", "IER Stromal", "HOXD-patterned Myofibroblast", "Steroidogenic Stromal",
  "Contractile Stromal", "Perivascular-like Stromal", "TNF-activated Stromal", "Activated Stromal",
  "ArtEnd", "CapEnd", "VenEnd", "LympEnd",
  "ArtSMC", "VenSMC", "SpecSMC", "Stressed Perivascular",
  "CD4+ T", "CD8+ T", "Macrophage",
  "nmSchwann"
)
tnfa_major_vlines <- c(8.5, 12.5, 16.5, 19.5)

tnfa_gene_order <- c("FOSL2", "JUNB", "CEBPD", "LITAF", "NFKB2", "NFKBIA", "KDM6B", # Inflammatory TF / NFkB response
  "ZFP36", "SOD2", "SAT1", "BTG1", "BTG2", "CDKN1A", "RHOB", "PER1", # Stress-adaptive state
  "KLF2", "KLF4", "KLF6", "KLF9", "MCL1", "TRAF1", # State / survival regulators
  "TAP1", "LDLR", "B4GALT1", "CCNL1", "EIF1") # Interface / metabolic context
tnfa_major_hlines <- c(5.5, 11.5, 19.5)

tnfa_heat_age <- make_tstat_heatmap(age_tnfa_consensus_pathway_drivers$tstat_top_mat, show_axis_text_x = FALSE, major_vlines = tnfa_major_vlines,
  major_hlines = tnfa_major_hlines, celltype_order = tnfa_celltype_order, gene_order = tnfa_gene_order, cbar_L = cbar_L_tnfa)
ggsave(filename = paste0(outs_dir, "/tnfa_heatmap_age.pdf"), plot = tnfa_heat_age, width = 35, height = 90, units = "mm", device = cairo_pdf)

tnfa_heat_para <- make_tstat_heatmap(para_tstat_top_mat_tnfa, show_axis_text_x = FALSE, major_vlines = tnfa_major_vlines, major_hlines = tnfa_major_hlines,
  celltype_order = tnfa_celltype_order, gene_order = tnfa_gene_order, cbar_L = cbar_L_tnfa)
ggsave(filename = paste0(outs_dir, "/tnfa_heatmap_para.pdf"), plot = tnfa_heat_para, width = 35, height = 90, units = "mm", device = cairo_pdf)

## Androgen response pathway ##
# Determine & plot consensus gene drivers for Androgen response pathway
age_androgen_consensus_pathway_drivers <- consensus_pathway_drivers(res.dl, res_zenith_all, gsc_filt, "HALLMARK_ANDROGEN_RESPONSE", "Age",
  outs_dir, save_csv = "none")

para_tstat_top_mat_androgen <- make_tstat_top_mat(
  res.dl, "Para_bin1-3", age_androgen_consensus_pathway_drivers$consensus_drivers_top$gene,
  colnames(age_androgen_consensus_pathway_drivers$tstat_top_mat)
)

cbar_L_androgen <- ceiling(max(
  max(abs(para_tstat_top_mat_androgen), na.rm = TRUE),
  max(abs(age_androgen_consensus_pathway_drivers$tstat_top_mat), na.rm = TRUE)
))

androgen_celltype_order <- c(
  "Homeostatic Fibroblast", "IER Stromal", "HOXD-patterned Myofibroblast", "Steroidogenic Stromal",
  "Contractile Stromal", "Perivascular-like Stromal", "TNF-activated Stromal", "Activated Stromal",
  "CapEnd", "VenEnd", "LympEnd",
  "ArtSMC", "SpecSMC",
  "Macrophage"
)
androgen_major_vlines <- c(8.5, 11.5, 13.5)

androgen_gene_order <- c(
  "ELOVL5", "HMGCS1", "FADS1", "HSD17B14", "VAPA", # Lipid / sterol metabolism
  "FKBP5", "ELL2", "ARID5B", "PA2G4", "AKAP12", # Hormone-response state
  "ADAMTS1", "ACTN1", "MERTK", "AKT1", # Remodeling / structural state
  "SAT1", "B4GALT1", "UAP1", "SLC38A2" # Adaptive biosynthetic support
)
androgen_major_hlines <- c(4.5, 8.5, 13.5)

androgen_heat_age <- make_tstat_heatmap(age_androgen_consensus_pathway_drivers$tstat_top_mat, show_axis_text_x = FALSE, major_vlines = androgen_major_vlines,
  major_hlines = androgen_major_hlines, celltype_order = androgen_celltype_order, gene_order = androgen_gene_order, cbar_L = cbar_L_androgen)
ggsave(filename = paste0(outs_dir, "/androgen_heatmap_age.pdf"), plot = androgen_heat_age, width = 35, height = 90, units = "mm", device = cairo_pdf)

androgen_heat_para <- make_tstat_heatmap(para_tstat_top_mat_androgen, show_axis_text_x = FALSE, major_vlines = androgen_major_vlines, major_hlines = androgen_major_hlines,
  celltype_order = androgen_celltype_order, gene_order = androgen_gene_order, cbar_L = cbar_L_androgen)
ggsave(filename = paste0(outs_dir, "/androgen_heatmap_para.pdf"), plot = androgen_heat_para, width = 35, height = 90, units = "mm", device = cairo_pdf)

## EMT pathway ##
# Determine & plot consensus gene drivers for EMT pathway
age_emt_consensus_pathway_drivers <- consensus_pathway_drivers(res.dl, res_zenith_all, gsc_filt, "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "Age",
  outs_dir, save_csv = "none")

para_tstat_top_mat_emt <- make_tstat_top_mat(
  res.dl, "Para_bin1-3", age_emt_consensus_pathway_drivers$consensus_drivers_top$gene,
  colnames(age_emt_consensus_pathway_drivers$tstat_top_mat)
)

cbar_L_emt <- ceiling(max(
  max(abs(para_tstat_top_mat_emt), na.rm = TRUE),
  max(abs(age_emt_consensus_pathway_drivers$tstat_top_mat), na.rm = TRUE)
))

# Order cell types by compartment
emt_celltype_order <- c(
  "Homeostatic Fibroblast", "IER Stromal", "HOXD-patterned Myofibroblast",
  "Contractile Stromal", "Perivascular-like Stromal", "TNF-activated Stromal",
  "CapEnd", "VenEnd", "LympEnd",
  "ArtSMC", "SpecSMC",
  "nmSchwann"
)
emt_major_vlines <- c(6.5, 9.5, 11.5, 12.5)

emt_gene_order <- c(
  "TIMP1", "TIMP3", "THBS1", "COL6A2", "EFEMP2", "QSOX1", "P3H1", "PPIB", "CALU", "MGP", # ECM remodeling
  "TPM4", "TPM2", "TAGLN", "MYLK", "MYL9", "FLNA", "PDLIM4", "RHOB", # Contractile cytoskeleton
  "ITGA5", "FERMT2", "PVR", "EMP3", # Adhesion / interface
  "FSTL1", "FSTL3", "IGFBP2", "IGFBP4", "SAT1", "GADD45B", "NNMT") # Activated stromal state
emt_major_hlines <- c(7.5, 11.5, 19.5)

emt_heat_age <- make_tstat_heatmap(age_emt_consensus_pathway_drivers$tstat_top_mat, show_axis_text_x = FALSE, major_vlines = emt_major_vlines,
  major_hlines = emt_major_hlines, celltype_order = emt_celltype_order, gene_order = emt_gene_order, cbar_L = cbar_L_emt)
ggsave(filename = paste0(outs_dir, "/emt_heatmap_age.pdf"), plot = emt_heat_age, width = 35, height = 90, units = "mm", device = cairo_pdf)

emt_heat_para <- make_tstat_heatmap(para_tstat_top_mat_emt, show_axis_text_x = FALSE, major_vlines = emt_major_vlines, major_hlines = emt_major_hlines,
  celltype_order = emt_celltype_order, gene_order = emt_gene_order, cbar_L = cbar_L_emt)
ggsave(filename = paste0(outs_dir, "/emt_heatmap_para.pdf"), plot = emt_heat_para, width = 35, height = 90, units = "mm", device = cairo_pdf)







#### 10. Session Info ####
sessionInfo()

# R version 4.4.1 (2024-06-14)
# Platform: x86_64-conda-linux-gnu
# Running under: Ubuntu 22.04.5 LTS

# Matrix products: default
# BLAS/LAPACK: .../lib/libopenblasp-r0.3.30.so;  LAPACK version 3.12.0

# locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
#  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C
#  [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

# time zone: America/Los_Angeles
# tzcode source: system (glibc)

# attached base packages:
# [1] stats4    grid      stats     graphics  grDevices utils     datasets  methods   base

# other attached packages:
#  [1] ggrepel_0.9.6               car_3.1-3                   carData_3.0-6               muscat_1.22.0
#  [5] GSEABase_1.70.0             graph_1.86.0                annotate_1.86.1             XML_3.99-0.18
#  [9] AnnotationDbi_1.70.0        zenith_1.10.0               dreamlet_1.4.1              SingleCellExperiment_1.30.1
# [13] SummarizedExperiment_1.38.1 Biobase_2.68.0              GenomicRanges_1.60.0        GenomeInfoDb_1.44.0
# [17] IRanges_2.42.0              S4Vectors_0.48.0            BiocGenerics_0.54.0         generics_0.1.4
# [21] MatrixGenerics_1.20.0       matrixStats_1.5.0           variancePartition_1.38.0    limma_3.64.1
# [25] RhpcBLASctl_0.23-42         BiocParallel_1.42.1         patchwork_1.3.1             UCell_2.12.0
# [29] lubridate_1.9.4             forcats_1.0.0               stringr_1.5.1               dplyr_1.1.4
# [33] purrr_1.1.0                 readr_2.1.5                 tidyr_1.3.1                 tibble_3.3.0
# [37] ggplot2_3.5.2               tidyverse_2.0.0

# loaded via a namespace (and not attached):
#   [1] fs_1.6.6                  bitops_1.0-9              doParallel_1.0.17         httr_1.4.7
#   [5] RColorBrewer_1.1-3        Rgraphviz_2.52.0          numDeriv_2016.8-1.1       sctransform_0.4.2
#   [9] tools_4.4.1               backports_1.5.0           utf8_1.2.6                R6_2.6.1
#  [13] metafor_4.8-0             HDF5Array_1.34.0          mgcv_1.9-3                rhdf5filters_1.20.0
#  [17] GetoptLong_1.0.5          withr_3.0.2               prettyunits_1.2.0         gridExtra_2.3
#  [21] cli_3.6.5                 box_1.2.1                 labeling_0.4.3            KEGGgraph_1.68.0
#  [25] SQUAREM_2021.1            mvtnorm_1.3-3             blme_1.0-6                proxy_0.4-27
#  [29] mixsqp_0.3-54             systemfonts_1.2.3         scater_1.34.1             parallelly_1.46.1
#  [33] invgamma_1.2              readxl_1.4.5              rstudioapi_0.17.1         RSQLite_2.4.2
#  [37] shape_1.4.6.1             vroom_1.6.5               gtools_3.9.5              Matrix_1.7-3
#  [41] metadat_1.4-0             ggbeeswarm_0.7.2          DescTools_0.99.60         abind_1.4-8
#  [45] lifecycle_1.0.5           edgeR_4.6.3               mathjaxr_2.0-0            rhdf5_2.52.1
#  [49] gplots_3.3.0              SparseArray_1.8.0         blob_1.2.4                crayon_1.5.3
#  [53] dir.expiry_1.16.0         lattice_0.22-7            haven_2.5.5               beachmat_2.24.0
#  [57] msigdbr_25.1.1            KEGGREST_1.48.1           unigd_0.1.3               pillar_1.11.0
#  [61] knitr_1.50                ComplexHeatmap_2.24.1     rjson_0.2.23              boot_1.3-32
#  [65] gld_2.6.7                 estimability_1.5.1        corpcor_1.6.10            future.apply_1.20.0
#  [69] codetools_0.2-20          glue_1.8.0                data.table_1.18.2.1       vctrs_0.7.1
#  [73] png_0.1-8                 Rdpack_2.6.4              cellranger_1.1.0          gtable_0.3.6
#  [77] assertthat_0.2.1          cachem_1.1.0              zigg_0.0.2                xfun_0.56
#  [81] rbibutils_2.3             S4Arrays_1.8.1            Rfast_2.1.5.1             coda_0.19-4.1
#  [85] reformulas_0.4.1          pheatmap_1.0.13           iterators_1.0.14          statmod_1.5.0
#  [89] nlme_3.1-168              pbkrtest_0.5.5            bit64_4.6.0-1             progress_1.2.3
#  [93] EnvStats_3.1.0            filelock_1.0.3            TMB_1.9.17                irlba_2.3.5.1
#  [97] vipor_0.4.7               KernSmooth_2.23-26        colorspace_2.1-2          rmeta_3.0
# [101] DBI_1.2.3                 zellkonverter_1.18.0      DESeq2_1.48.1             Exact_3.3
# [105] tidyselect_1.2.1          emmeans_1.11.1            bit_4.6.0                 compiler_4.4.1
# [109] curl_7.0.0                BiocNeighbors_2.0.1       expm_1.0-0                basilisk.utils_1.20.0
# [113] DelayedArray_0.34.1       scales_1.4.0              caTools_1.18.3            remaCor_0.0.18
# [117] httpgd_2.0.4              digest_0.6.39             minqa_1.2.8               basilisk_1.20.0
# [121] aod_1.3.3                 XVector_0.48.0            pkgconfig_2.0.3           lme4_1.1-37
# [125] sparseMatrixStats_1.20.0  mashr_0.2.79              fastmap_1.2.0             GlobalOptions_0.1.3
# [129] rlang_1.1.7               UCSC.utils_1.4.0          DelayedMatrixStats_1.30.0 farver_2.1.2
# [133] jsonlite_2.0.0            BiocSingular_1.22.0       RCurl_1.98-1.17           magrittr_2.0.4
# [137] Formula_1.2-5             scuttle_1.16.0            GenomeInfoDbData_1.2.14   Rhdf5lib_1.30.0
# [141] Rcpp_1.1.1                babelgene_22.9            viridis_0.6.5             reticulate_1.42.0
# [145] EnrichmentBrowser_2.38.0  stringi_1.8.7             rootSolve_1.8.2.4         MASS_7.3-65
# [149] plyr_1.8.9                listenv_0.10.0            parallel_4.4.1            lmom_3.2
# [153] Biostrings_2.76.0         splines_4.4.1             circlize_0.4.16           hms_1.1.3
# [157] locfit_1.5-9.12           reshape2_1.4.4            ScaledMatrix_1.14.0       evaluate_1.0.5
# [161] RcppParallel_5.1.11-1     foreach_1.5.2             nloptr_2.2.1              tzdb_0.5.0
# [165] future_1.58.0             clue_0.3-66               scattermore_1.2           ashr_2.2-63
# [169] rsvd_1.0.5                broom_1.0.8               xtable_1.8-4              fANCOVA_0.6-1
# [173] e1071_1.7-16              viridisLite_0.4.2         class_7.3-23              truncnorm_1.0-9
# [177] glmmTMB_1.1.11            lmerTest_3.1-3            memoise_2.0.1             beeswarm_0.4.0
# [181] cluster_2.1.8.1           globals_0.19.0            timechange_0.3.0