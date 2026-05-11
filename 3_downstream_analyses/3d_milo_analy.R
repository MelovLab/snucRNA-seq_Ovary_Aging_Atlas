## Script name: Human Ovary Atlas - Milo differential abundance analysis on Immune compartment
##
## Purpose: Run Milo neighborhood-level differential abundance analysis in the immune compartment,
## generate DA beeswarm and UMAP projection plots, identify age-associated macrophage neighborhood markers,
## perform pathway analysis, and export supplementary result tables.
##
## Author: Josef Byrne
#########################

#### 1. Import libraries and set directories ####

# Import libraries
library(zellkonverter)
library(tidyverse)
library(ggplot2)
library(ggbeeswarm)
library(scales)
library(SingleCellExperiment)
library(miloR) # v2.6.0
library(scater)
library(patchwork)
library(dplyr)
library(tibble)
library(openxlsx)
library(stringr)
library(limma)
library(Matrix)
library(SummarizedExperiment)
library(msigdbr)

# Initialize paths
nostromo_datadir="/path/to/project_directory"
adata_objdir <- paste0(nostromo_datadir, "/python_analyses/adata_objs")
r_analyses_dir <- paste0(nostromo_datadir, "/R_analyses")
outs_dir <- paste0(r_analyses_dir, "/milo_objs")

# Import snucRNA-seq pseudobulked metadata (from DGE preprocessing) that has per-sample aggregated StressScore values calculated from DGE & attach to sce_obj colData
meta_df <- readRDS(file.path(r_analyses_dir, "donor_pb_meta_df.rds"))

# Import cell type compartment-specific data as sce_objs
immune_sce <- readH5AD(paste0(adata_objdir, "/immune_annot.h5ad"), reader = "R", verbose = TRUE)

# Ensure 'plex' column in colData(sce) is a factor (not integer)
# Add compartment label and format metadata for Milo modeling
immune_sce$celltype_compartment <- "immune"

sce_list <- list(immune_sce)
names(sce_list) <- c("immune")

# Ensure 'plex' column in colData(sce_list) is a factor (not integer)
for (sce_name in names(sce_list)) {
  # Convert 'plex' to factor, rename to 'Plex', drop old 'plex'
  sce <- sce_list[[sce_name]]
  colData(sce)$Plex <- as.factor(colData(sce)$plex)
  colData(sce)$plex <- NULL

  # Add Para_bin as a factor with 0 as the reference and 1-3 as the reported parity contrast
  para_vals <- as.numeric(colData(sce)$Para)
  para_bin <- cut(
    para_vals,
    breaks = c(-Inf, 0, 3, Inf),
    labels = c("0", "1-3", "4+"),
    right = TRUE
  )
  # Actually add Para_bin to the sce
  colData(sce)$Para_bin <- factor(as.character(para_bin), levels = c("0", "4+", "1-3"))

  # Make Age_per_dec variable: Age divided by 10
  if ("Age" %in% colnames(colData(sce))) {
    colData(sce)$Age_per_dec <- colData(sce)$Age / 10
  }

  # Add StressScore at the sample level (match unique_sample_ID to meta_df, from DGE preprocessing)
  if ("unique_sample_ID" %in% colnames(colData(sce))) {
    colData(sce)$StressScore <- meta_df$StressScore[
      match(colData(sce)$unique_sample_ID, rownames(meta_df))
    ]
  }

  sce_list[[sce_name]] <- sce
}




#### 2. Prep & run Milo workflow ####

## Milo workflow function
run_milo <- function(
  sce,               # SingleCellExperiment used to build Milo
  k,                 # Number of nearest neighbors for Milo
  prop,              # Proportion of neighbors to sample
  reduced.dim,       # Name of reduced dimension (e.g. "PCA", "UMAP", etc.)
  da_design_forms,   # List of formulas describing the design for DA testNhoods
  design_df,         # Data frame with experimental/metadata for model matrix
  donor_col,         # Name of donor column in colData
  sample_col,        # Name of sample column in colData
  seed = 42          # Random seed for reproducibility
) {
  set.seed(seed)

  num_dims <- dim(reducedDim(sce, reduced.dim))[2] # number of dimensions in the reduced dimension
  milo <- Milo(sce)
  milo <- buildGraph(milo, k = k, d = num_dims, reduced.dim = reduced.dim)
  milo <- makeNhoods(milo, k = k, prop = prop, d = num_dims, reduced_dim = reduced.dim, refined = TRUE)
  milo <- countCells(milo, meta.data = data.frame(colData(milo)), samples = sample_col)
  milo <- calcNhoodDistance(milo, d = num_dims, reduced.dim = reduced.dim)
  milo <- buildNhoodGraph(milo, overlap = k/2) # Only show connections between neighborhoods that are at least half # of cells

    # Count and filter nhoods so that
  # a) >=10 non-zero donors are present per nhood
  # b) 15 non-zero samples (tissue chunks, unit of analysis with donor random effect) are present per nhood
  # c) 3 donors don't make up more than 80% of nhood

  # Map sample IDs in milo_obj (colnames(cts)) to SenNet_ID for each sample
  cts <- nhoodCounts(milo)   # nhood x sample
  sample_ids <- colnames(cts)
  sample_metadata <- as.data.frame(milo@colData)
  sample_id_to_donor <- setNames(sample_metadata[[donor_col]], sample_metadata[[sample_col]])
  donor <- sample_id_to_donor[sample_ids]

  # a) sample breadth
  sample_nonzero <- rowSums(cts > 0)

  # donor-aggregated matrix (nhood x donor)
  donor_levels <- unique(donor)
  cts_donor <- sapply(donor_levels, function(d) Matrix::rowSums(cts[, donor == d, drop=FALSE]))
  cts_donor <- as.matrix(cts_donor)

  # b) donor breadth
  donor_nonzero <- rowSums(cts_donor > 0)

  # c) donor dominance (fractions from top donors)
  top_donors_num <- 3
  tot <- rowSums(cts_donor); tot[tot == 0] <- 1
  top_donors_cts <- apply(cts_donor, 1, function(x) sum(sort(x, decreasing=TRUE)[1:top_donors_num]))
  top_donors_frac <- top_donors_cts / tot

  # Filtering criteria
  keep <- (donor_nonzero >= 10) & (sample_nonzero >= 15) & (top_donors_frac <= 0.80)

  # Reorder design_df rownames to match columns of nhoodCounts(milo)
  design_df <- design_df[colnames(nhoodCounts(milo)), , drop = FALSE]

  # Run DA testNhoods for each design formula
  da_res_list <- list()
  for (da_name in names(da_design_forms)) {
    # don't use BPPARAM here, OpenMP parallelization will oversubscribe, use all cores and cause thrashing --> slow/ineffective
    # Fisher scoring initially tried, but near all produced negative variance estimates. To ensure numerically stable results, switched to HE-NNLS.
    da_res_list[[da_name]] <- testNhoods(milo, design = da_design_forms[[da_name]], design.df = design_df,
      reduced.dim = reduced.dim, fdr.weighting = "graph-overlap", glmm.solver = "HE-NNLS",
      subset.nhoods = keep, REML = TRUE, norm.method = "TMM", max.iters = 100, force = TRUE)
  }

  return(list(milo = milo, da_res = da_res_list, nhoods_kept = keep, k = k, prop = prop))
}

## Prepare design matrix and run Milo workflow
# Assign vars
dim_red_name <- "X_scVI"
vis_red_name <- "X_umap"
k_neighs <- 100 # Targeting rule of thumb of average neighborhood size ~ 5 x N biological replicates (28 here) or ~140 cells per neighborhood
donor_col <- "SenNet_ID"
sample_col <- "unique_sample_ID"
cell_type_col <- "cell_type"

form_para <- paste0("~ Age_per_dec + BMI + Plex + StressScore + Para_bin + (1|", donor_col, ")")
form_age <- paste0("~ Para_bin + BMI + Plex + StressScore + Age_per_dec + (1|", donor_col, ")")
da_design_forms <- list(age_per_dec = as.formula(form_age), para_bin = as.formula(form_para))

milo_runs_list <- list()
for (i in seq_along(sce_list)) {
  sce_i <- sce_list[[i]]
  sce_name <- names(sce_list)[i]

  # Determine dynamic prop and k-neighbors based on the number of cells in the compartment
  num_cells <- ncol(sce_i)
  prop_nhoods = 0.1 # May need to make variable based on compartment size

  # Prepare design matrix for each SCE
  design_mat_i <- colData(sce_i) %>%
    as.data.frame() %>%
    dplyr::select(all_of(c(donor_col, sample_col, "Age_per_dec", "BMI", "Para_bin", "Plex", "StressScore"))) %>%
    dplyr::distinct() %>%
    magrittr::set_rownames(.[[sample_col]]) %>%
    dplyr::mutate(
      Para_bin = factor(Para_bin, levels = c("0", "4+", "1-3")), # Set 0 as reference and 1-3 as the reported parity contrast; 4+ is retained as a nuisance level
      Plex = factor(Plex)
    )

  # Create output list to store results, warnings and messages
  outs <- list()
  warnings_collected <- character()
  messages_collected <- character()

  # Capture warnings and messages and save to file
  milo_run_i <- withCallingHandlers(
    {
      result <- run_milo(
        sce = sce_i, k = k_neighs, prop = prop_nhoods, reduced.dim = dim_red_name,
        da_design_forms = da_design_forms, design_df = design_mat_i, donor_col = donor_col, sample_col = sample_col, seed = 42
      )
      outs$result <- result
      result
    },
    warning = function(w) {
      warnings_collected <<- c(warnings_collected, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      messages_collected <<- c(messages_collected, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  outs$warnings <- warnings_collected
  outs$messages <- messages_collected

  prop_nhoods_rounded <- gsub("\\.", "_", as.character(round(prop_nhoods, 2)))
  warnings_file <- file.path(outs_dir, paste0("milo_run_", sce_name, "_k", k_neighs, "_prop", prop_nhoods_rounded, "_warnings_messages.txt"))
  writeLines(
    c(
      "WARNINGS:",
      if (length(warnings_collected) > 0) warnings_collected else "None",
      "",
      "MESSAGES:",
      if (length(messages_collected) > 0) messages_collected else "None"
    ),
    con = warnings_file
  )

  print(paste0("Warnings and messages saved to: ", warnings_file))
  saveRDS(outs, file = file.path(outs_dir, paste0("milo_run_", sce_name, "_k", k_neighs, "_prop", prop_nhoods_rounded, ".rds")))

  milo_runs_list[[sce_name]] <- outs
}



#### 3. Visualize results with beeswarm & UMAP projection plots ####

## 3.a Plot beeswarm plots for cell type level DA results
## Modified version of plotDAbeeswarm to allow for custom adjustments
plotDAbeeswarm_mod <- function(da.res, group.by = NULL, alpha = 0.1, subset.nhoods = NULL, base_font_size = 7, group.by.axis.title = NULL, point_size = 0.3, ylim = NULL) {
  if (!is.null(group.by)) {
    if (!group.by %in% colnames(da.res)) {
      stop(group.by, " is not a column in da.res. Have you forgot to run annotateNhoods(x, da.res, ", group.by,")?")
    }
    if (is.numeric(da.res[, group.by])) {
      # stop(group.by, " is a numeric variable. Please bin to use for grouping.")
    }
    da.res <- mutate(da.res, group_by = da.res[, group.by])
  } else {
    da.res <- mutate(da.res, group_by = "g1")
  }

  if (!is.factor(da.res[, "group_by"])) {
    message("Converting group_by to factor...")
    da.res <- mutate(da.res, group_by = factor(group_by, levels = unique(group_by)))
  }

  if (!is.null(subset.nhoods)) {
    da.res <- da.res[subset.nhoods, ]
  }

  # Set the x-axis label (either use provided title, or fall back to group.by)
  x_axis_label <- if(!is.null(group.by.axis.title)) group.by.axis.title else group.by

  # Get position with ggbeeswarm
  beeswarm_pos <- ggplot_build(
    da.res %>%
      mutate(is_signif = ifelse(SpatialFDR < alpha, 1, 0)) %>%
      arrange(group_by) %>%
      ggplot(aes(group_by, logFC)) +
      geom_quasirandom()
  )

  pos_x <- beeswarm_pos$data[[1]]$x
  pos_y <- beeswarm_pos$data[[1]]$y

  n_groups <- unique(da.res$group_by) %>% length()

  plot_obj <- da.res %>%
    mutate(is_signif = ifelse(SpatialFDR < alpha, 1, 0)) %>%
    mutate(logFC_color = ifelse(is_signif == 1, logFC, NA)) %>%
    arrange(group_by) %>%
    mutate(Nhood = factor(Nhood, levels = unique(Nhood))) %>%
    mutate(pos_x = pos_x, pos_y = pos_y) %>%
    ggplot(aes(pos_x, pos_y, color = logFC_color)) +
    scale_color_gradient2(low = "#0080ff", mid = "grey95", high = "#e70b25") +
    guides(color = "none") +
    xlab(x_axis_label) + ylab("Log Fold Change") +
    scale_x_continuous(
      breaks = seq(1, n_groups),
      labels = setNames(levels(da.res$group_by), seq(1, n_groups))
    ) +
    #scale_y_continuous(breaks = function(lims) seq(floor(min(lims, na.rm = TRUE)/0.5)*0.5, ceiling(max(lims, na.rm = TRUE)/0.5)*0.5, by = 0.5)) +
    geom_point(size = point_size) +
    theme_bw(base_size = base_font_size) +
    theme(
      strip.text.y = element_text(angle = 0, color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      axis.title.x = element_text(margin = margin(t = 6), color = "black"),  # more padding top
      axis.title.y = element_text(margin = margin(r = 6), color = "black"),  # more padding right
      plot.title = element_text(color = "black"),
      legend.text = element_text(color = "black"),
      legend.title = element_text(color = "black"),
      axis.ticks.length = unit(0.08, "cm"),
      axis.text.x = element_text(margin = margin(t = 2)), # increase tick label distancing x
      axis.text.y = element_text(margin = margin(r = 2))  # increase tick label distancing y
    )

  # If ylim provided, add it.
  if (!is.null(ylim)) {
    plot_obj <- plot_obj +
      scale_y_continuous(
        limits = ylim,
        expand = expansion(mult = c(0, 0))
      )
  }

  return(plot_obj)
}

## Function to plot beeswarm plots for cell type level DA results
plot_milo_beeswarm <- function(milo_runs, sce_name, coef, mixed_exclude = FALSE, min_neighs = 1, ylim = NULL) {
  milo_obj <- milo_runs[[sce_name]]$result$milo
  keep <- milo_runs[[sce_name]]$result$nhoods_kept
  da_results_var <- milo_runs[[sce_name]]$result$da_res[[coef]]
  # Skip if this result is missing for this compartment
  if (is.null(da_results_var)) {
    message(sprintf("Skipping %s - %s (no results)", sce_name, coef))
    return(invisible(NULL))
  }
  da_results_var_clean <- annotateNhoods(milo_obj, da_results_var, coldata_col = cell_type_col, subset.nhoods = keep)
  da_results_var_clean <- da_results_var_clean[!is.na(da_results_var_clean$logFC) & da_results_var_clean$Converged == TRUE, ]

  # Assign "Mixed" as before
  da_results_var_clean$cell_type <- ifelse(
    da_results_var_clean$cell_type_fraction < 0.70,
    "Mixed",
    da_results_var_clean$cell_type
  )

  # Exclude "Mixed" if mixed_exclude is TRUE
  if (mixed_exclude) {
    da_results_var_clean <- da_results_var_clean[da_results_var_clean$cell_type != "Mixed", ]
  }

  # Exclude cell_type if it has fewer than min_neighs neighborhoods
  if (!is.null(min_neighs) && min_neighs > 1) {
    nhood_counts <- table(da_results_var_clean$cell_type)
    valid_types <- names(nhood_counts)[nhood_counts >= min_neighs]
    da_results_var_clean <- da_results_var_clean[da_results_var_clean$cell_type %in% valid_types, ]
  }

  plot_obj <- plotDAbeeswarm_mod(
    da_results_var_clean,
    group.by = cell_type_col,
    alpha = 0.1,
    base_font_size = 5,
    group.by.axis.title = "Cell Type",
    point_size = 0.5,
    ylim = ylim
  ) + ggtitle(
    paste0("miloR DA: ", tools::toTitleCase(coef), " - ", sce_name)
  ) + theme(
    plot.title = element_text(hjust = 0.5)
  )
  return(plot_obj)
}

para_immune_bs <- plot_milo_beeswarm(milo_runs = milo_runs_list, sce_name = "immune", coef = "para_bin", mixed_exclude = TRUE, min_neighs = 10)
# Get the y axis limits for para_immune_bs to use for age_immune_bs for easier visual comparison of the two effects
para_y_lims <- ggplot_build(para_immune_bs)$layout$panel_params[[1]]$y.range

age_immune_bs <- plot_milo_beeswarm(milo_runs = milo_runs_list, sce_name = "immune", coef = "age_per_dec", mixed_exclude = TRUE, min_neighs = 10, ylim = para_y_lims)

ggsave(filename = paste0(outs_dir, "/Milo_Age_DA_immune_beeswarm.pdf"), plot = age_immune_bs, width = 80, height = 40, units = "mm", device = cairo_pdf, dpi = 500)
ggsave(filename = paste0(outs_dir, "/Milo_Para_DA_immune_beeswarm.pdf"), plot = para_immune_bs, width = 80, height = 40, units = "mm", device = cairo_pdf, dpi = 500)

### 3.b Plot UMAP projection plots for cell type level DA results
## Function to plot a subset of the nhood graph with DA results
plotNhoodGraphDA_subset <- function(
  x,                # Milo object
  da_res,           # Data frame returned from testNhoods(), containing DA results per neighborhood (incl. Nhood, logFC, PValue, etc.)
  nhoods_keep,      # Logical vector indicating which neighborhoods to retain for plotting (length = number of nhoods in Milo object)
  red_name,         # Name of the reduced dimensionality layout to use for plotting (string)
  alpha,            # Numeric [0,1]: transparency for plot points
  ...               # Additional arguments passed to plotNhoodGraphDA
) {

  #   Notes:
  #   - Current standard plotNhoodGraphDA() without modification is unable to handle subsets of neighborhoods.
  #   - Only neighborhoods present both in the Milo object's nhood graph and in 'nhoods_keep'
  #     are displayed.
  #   - Designed for fast visualization of DA results restricted to QC-filtered neighborhoods.
  # ------------------------------------------------------------------------------

  g <- nhoodGraph(x)

  # seed cell per nhood (this is the mapping between Nhood index and graph vertex name)
  seed_per_nhood <- vapply(nhoodIndex(x), `[`, numeric(1), 1)

  # graph vertex names are seed IDs as character
  graph_vids <- as.character(igraph::V(g)$name)

  # create a column that matches graph vertex names
  da_res$vid <- as.character(seed_per_nhood[da_res$Nhood])

  keep_nhood <- which(nhoods_keep) # indices 1..K where keep==TRUE
  keep_vid <- as.character(seed_per_nhood[keep_nhood])

  da_sub <- da_res %>%
    dplyr::mutate(vid = as.character(seed_per_nhood[Nhood])) %>%
    dplyr::filter(vid %in% graph_vids, vid %in% keep_vid)

  g_sub <- igraph::induced_subgraph(g, vids = unique(da_sub$vid))

  milo_tmp <- x
  milo_tmp@nhoodGraph$nhoodGraph <- g_sub

  return(plotNhoodGraphDA(milo_tmp, da_sub, layout = red_name, alpha = alpha, ...))
}

# Function to plot the nhood graph with DA results
plot_milo_da_umap <- function(
    milo_runs_list,
    compartment,
    var,
    base_size = 5,
    bold_sig = FALSE,
    red_name = vis_red_name
) {
  # Check that the requested compartment exists
  stopifnot(compartment %in% names(milo_runs_list))

  # Pull out the stored Milo results for this compartment
  res <- milo_runs_list[[compartment]]$result
  milo <- res$milo
  da_res <- res$da_res[[var]]
  nhoods_kept <- res$nhoods_kept

  # Rebuild graph with a stricter overlap threshold for clearer UMAP projection visualization
  milo <- buildNhoodGraph(milo, overlap = 100)

  # Compute a symmetric, rounded color limit for logFC
  # This keeps the blue/red scale balanced around 0
  get_pretty_limit <- function(x, q = 0.95) {
    #L <- quantile(abs(x), q, na.rm = TRUE)
    L <- max(abs(x), na.rm = TRUE)
    if (!is.finite(L) || L == 0) return(1)
    # exp <- floor(log10(abs(L)))
    # round_to <- 10^exp
    # ceiling(L / round_to) * round_to
    ceiling(L / 0.1) * 0.1
  }

  # Find point layers programmatically rather than hard-coding layer numbers
  get_point_layers <- function(p) {
    which(vapply(p$layers, function(x) inherits(x$geom, "GeomPoint"), logical(1)))
  }

  # Apply the shared final styling used for all plots
  style_plot <- function(p, L) {
    p <- p +
      ggraph::scale_edge_width(range = c(0.05, 1), guide = "none") +
      scale_fill_gradient2(
        low = "#0080ff",
        mid = "grey95",
        high = "#e70b25",
        midpoint = 0,
        limits = c(-L, L),
        oob = scales::squish,
        name = "logFC",
        breaks = c(-L, 0, L),
        labels = c(sprintf("-%g", L), "0", sprintf("%g", L))
      ) +
      scale_size_continuous(
        range = c(0.1, 1.5),
        name = "Nhood Size",
        guide = guide_legend(
          order = 2,
          keywidth = unit(6, "pt"),
          keyheight = unit(8, "pt"),
          title.theme = element_text(size = base_size),
          label.theme = element_text(size = base_size)
        )
      ) +
      guides(
        fill = guide_colorbar(
          order = 1,
          barwidth = unit(6, "pt"),
          barheight = unit(40, "pt"),
          title.theme = element_text(size = base_size, margin = margin(b = 6)),
          label.theme = element_text(size = base_size)
        )
      ) +
      labs(x = "UMAP 1", y = "UMAP 2") +
      theme_minimal(base_family = "DejaVu Sans", base_size = base_size) +
      theme(
        panel.border = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.2),
        axis.line.y.right = element_blank(),
        axis.line.x.top = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        legend.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_text(margin = margin(t = 4)),
        axis.title.y = element_text(margin = margin(r = 4))
      )

    # Force edge styling to match the original bold-outline version
    edge_idx <- which(vapply(p$layers, function(x) class(x$geom)[1] == "GeomEdgeSegment", logical(1)))[1]
    if (!is.na(edge_idx)) {
      p$layers[[edge_idx]]$aes_params$edge_colour <- "gray50"
      p$layers[[edge_idx]]$aes_params$edge_alpha <- 0.2
    }

    p
  }

  # Determine color scale limit from the observed logFC values
  L <- get_pretty_limit(da_res$logFC)

  # Base plot: draw all neighborhoods with the final shared styling
  p <- plotNhoodGraphDA_subset(
    milo,
    da_res,
    nhoods_kept,
    red_name = red_name,
    alpha = 1,
    highlight.da = FALSE
  )

  # Give all nodes a thin border
  point_layers <- get_point_layers(p)
  for (i in point_layers) {
    p$layers[[i]]$aes_params$stroke <- 0.1
  }

  # If requested, overlay a second layer containing only significant nodes
  # with a thicker border
  if (bold_sig) {
    p_sig <- plotNhoodGraphDA_subset(
      milo,
      da_res,
      nhoods_kept,
      red_name = red_name,
      alpha = 0.1,
      highlight.da = 0
    )

    # Style non-significant point layers in the base plot with light gray fill
    all_point_layers <- get_point_layers(p)
    if (length(all_point_layers) > 0) {
      for (i in all_point_layers) {
        # Set fill (color) to light gray for non-significant neighborhoods
        p$layers[[i]]$aes_params$fill <- "#e0e0e0"
        # If the layer also uses color aesthetic (not just fill), set that too
        if (!is.null(p$layers[[i]]$aes_params$colour)) {
          p$layers[[i]]$aes_params$colour <- "#e0e0e0"
        }
        p$layers[[i]]$aes_params$stroke <- 0.01
        # Lower alpha for the fill
        p$layers[[i]]$aes_params$alpha <- 0.3
        # Also decrease the alpha of the border (edge/line/colour/outline)
        if (!is.null(p$layers[[i]]$aes_params$stroke)) {
          # If the layer supports edge_alpha or linewidth alpha control, set a low value
          p$layers[[i]]$aes_params$edge_alpha <- 0.01
        }
      }
    }
    # Style significant point layers and add only those back (on top)
    sig_point_layers <- get_point_layers(p_sig)
    for (i in sig_point_layers) {
      p_sig$layers[[i]]$aes_params$stroke <- 0.3
      p <- p + p_sig$layers[[i]]
    }
  }

  # Apply shared styling last
  p <- style_plot(p, L)

  return(p)
}

age_immune_umap <- plot_milo_da_umap(milo_runs_list = milo_runs_list, compartment = "immune", var = "age_per_dec", base_size = 5, bold_sig = FALSE)
para_immune_umap <- plot_milo_da_umap(milo_runs_list = milo_runs_list, compartment = "immune", var = "para_bin", base_size = 5, bold_sig = FALSE)
ggsave(filename = paste0(outs_dir, "/Milo_Age_DA_immune_UMAP.pdf"), plot = age_immune_umap, width = 90, height = 45, units = "mm", device = cairo_pdf, dpi = 500)
ggsave(filename = paste0(outs_dir, "/Milo_Para_DA_immune_UMAP.pdf"), plot = para_immune_umap, width = 90, height = 45, units = "mm", device = cairo_pdf, dpi = 500)

age_plot_sig <- plot_milo_da_umap(milo_runs_list = milo_runs_list, compartment = "immune", var = "age_per_dec", base_size = 5, bold_sig = TRUE)
para_plot_sig <- plot_milo_da_umap(milo_runs_list = milo_runs_list, compartment = "immune", var = "para_bin", base_size = 5, bold_sig = TRUE)
ggsave(filename = paste0(outs_dir, "/Milo_Age_DA_immune_UMAP_sig.pdf"), plot = age_plot_sig, width = 90, height = 45, units = "mm", device = cairo_pdf, dpi = 500)
ggsave(filename = paste0(outs_dir, "/Milo_Para_DA_immune_UMAP_sig.pdf"), plot = para_plot_sig, width = 90, height = 45, units = "mm", device = cairo_pdf, dpi = 500)


#### 4. Identify age-associated macrophage neighborhood markers ####

# Function to perform donor-blocked DE across DA-defined Milo nhood groups
# Donor-blocked pos vs neg DE across DA-defined Milo nhood groups,
# using the SAME "unique cell assignment + sample×group aggregation" logic as findNhoodGroupMarkers(),
# but then fitting limma with donor blocking (duplicateCorrelation) and returning:
#   - topTable (gene-level DE)
#   - ranked statistic (moderated t)
milo_posneg_blocked_DE <- function(
  milo_obj,
  milo_da_res,
  subset.nhoods,                  # MUST be colnames(nhoods(milo_obj)) (seed IDs)
  assay_name = "log1p_norm",
  sample_col = "sample_ID",
  donor_col  = "SenNet_ID",
  groups_use = c("pos","neg"),
  gene.offset = FALSE,
  verbose = TRUE
) {
  stopifnot(inherits(milo_obj, "Milo"))
  if (!assay_name %in% assayNames(milo_obj)) stop("assay_name not found in milo_obj assays().")
  if (!all(c("Nhood","NhoodGroup") %in% colnames(milo_da_res))) stop("milo_da_res must contain Nhood and NhoodGroup.")
  if (!sample_col %in% colnames(colData(milo_obj))) stop("sample_col not found in milo_obj colData().")
  if (!donor_col  %in% colnames(colData(milo_obj))) stop("donor_col not found in milo_obj colData().")

  nhs_all <- nhoods(milo_obj)  # cells x nhoods
  if (!all(subset.nhoods %in% colnames(nhs_all))) {
    bad <- subset.nhoods[!subset.nhoods %in% colnames(nhs_all)]
    stop("subset.nhoods contains nhood names not found in nhoods(milo_obj). Example: ", paste(head(bad, 5), collapse=", "))
  }

  nhs <- nhs_all[, subset.nhoods, drop = FALSE]

  # Build vector of group labels for these nhood columns.
  # milo_da_res$Nhood is 1..N indexing columns of nhs_all; use that to align.
  # For each selected column name, find its column index in nhs_all, then map to milo_da_res row with Nhood==index.
  col_idx <- match(subset.nhoods, colnames(nhs_all))          # in 1..N
  da_idx  <- match(col_idx, milo_da_res$Nhood)                # rows in milo_da_res
  nhood_group <- as.character(milo_da_res$NhoodGroup[da_idx]) # "pos"/"neg"/"zero"
  names(nhood_group) <- subset.nhoods

  # Keep only requested groups at the nhood level
  keep_cols <- names(nhood_group)[nhood_group %in% groups_use & !is.na(nhood_group)]
  nhs <- nhs[, keep_cols, drop = FALSE]
  nhood_group <- nhood_group[keep_cols]

  if (verbose) {
    message("Selected nhoods provided: ", length(subset.nhoods))
    message("Selected nhoods retained (pos/neg): ", ncol(nhs))
    message("NhoodGroup counts: ", paste(names(table(nhood_group)), table(nhood_group), sep="=", collapse=", "))
  }
  if (ncol(nhs) == 0) stop("No nhoods left after filtering to groups_use.")

  # Keep cells that are in at least one retained nhood
  cells_keep <- rowSums(nhs) > 0
  nhs <- nhs[cells_keep, , drop = FALSE]
  cell_ids <- rownames(nhs)

  # Assign each cell to a group based on membership in any nhood of that group.
  fake_meta <- data.frame(CellID = cell_ids, Nhood.Group = NA_character_, stringsAsFactors = FALSE)
  rownames(fake_meta) <- fake_meta$CellID

  # Determine group membership per cell
  for (g in groups_use) {
    nhood_cols_g <- names(nhood_group)[nhood_group == g]
    if (!length(nhood_cols_g)) next
    in_g <- rowSums(nhs[, nhood_cols_g, drop = FALSE]) > 0

    # Drop cells represented in multiple comparison groups
    already <- !is.na(fake_meta[in_g, "Nhood.Group"])
    fake_meta[in_g, "Nhood.Group"] <- ifelse(already, NA_character_, g)
  }

  fake_meta <- fake_meta[!is.na(fake_meta$Nhood.Group), , drop = FALSE]
  if (verbose) {
    message("Cells retained after group assignment: ", nrow(fake_meta))
    message("Cell group counts: ", paste(names(table(fake_meta$Nhood.Group)), table(fake_meta$Nhood.Group), sep="=", collapse=", "))
  }

  if (length(unique(fake_meta$Nhood.Group)) < 2) {
    stop("Only one group present among cells (", unique(fake_meta$Nhood.Group),
         "). Relax cell type purity / DA thresholds to increase the number of cells in each group.")
  }

  # Expression genes x cells
  exprs <- SummarizedExperiment::assay(milo_obj, assay_name)[, fake_meta$CellID, drop = FALSE]

  # Add sample/donor IDs
  cd <- as.data.frame(colData(milo_obj)[fake_meta$CellID, , drop = FALSE])
  fake_meta$sample_id <- as.character(cd[[sample_col]])
  fake_meta$donor     <- as.character(cd[[donor_col]])
  fake_meta$BMI       <- as.numeric(cd$BMI)
  fake_meta$Plex      <- as.factor(cd$Plex)
  fake_meta$StressScore <- as.numeric(cd$StressScore)
  fake_meta$Para_bin  <- as.factor(cd$Para_bin)
  fake_meta$sample_group <- paste(fake_meta$sample_id, fake_meta$Nhood.Group, sep = "_")

  # Aggregate expression to sample_group (sample × group)
  sample_groups <- unique(fake_meta$sample_group)
  sg_mat <- matrix(0, nrow = nrow(fake_meta), ncol = length(sample_groups),
                   dimnames = list(rownames(fake_meta), sample_groups))
  for (s in sample_groups) sg_mat[fake_meta$sample_group == s, s] <- 1

  summFunc <- if (assay_name %in% c("counts","raw_counts")) rowSums else rowMeans
  exprs_smp <- matrix(0, nrow = nrow(exprs), ncol = ncol(sg_mat),
                      dimnames = list(rownames(exprs), colnames(sg_mat)))

  for (i in seq_len(ncol(sg_mat))) {
    idx <- which(sg_mat[, i] > 0)
    if (!length(idx)) next
    exprs_smp[, i] <- if (length(idx) > 1) summFunc(exprs[, idx, drop = FALSE]) else exprs[, idx]
  }

  # sample_group metadata
  smp_meta <- unique(fake_meta[, c("sample_group","sample_id","donor","Nhood.Group","BMI","Plex","StressScore","Para_bin")])
  rownames(smp_meta) <- smp_meta$sample_group
  smp_meta <- smp_meta[colnames(exprs_smp), , drop = FALSE]

  # limma pos vs neg with donor block
  smp_meta$Test <- factor(smp_meta$Nhood.Group, levels = c("neg","pos"))
  design <- model.matrix(~ BMI + Plex + StressScore + Para_bin + Test, data = smp_meta)
  colnames(design) <- gsub("^Test", "", colnames(design))
  colnames(design) <- make.names(colnames(design))
  print(colnames(design))

  if (isTRUE(gene.offset)) {
    n_genes <- apply(exprs_smp, 2, function(x) sum(x > 0))
    design <- cbind(design[, 1, drop=FALSE], NGenes = n_genes, design[, 2, drop=FALSE])
  }

  corfit <- duplicateCorrelation(exprs_smp, design, block = smp_meta$donor)
  fit <- lmFit(exprs_smp, design, block = smp_meta$donor, correlation = corfit$consensus.correlation)
  fit2 <- eBayes(fit, trend = TRUE)

  tt <- topTable(fit2, coef = "pos", number = Inf, sort.by = "P")
  rank_stat <- tt$t
  names(rank_stat) <- rownames(tt)
  rank_stat <- sort(rank_stat, decreasing = TRUE)

  if (verbose) {
    message("Consensus within-donor correlation: ", signif(corfit$consensus.correlation, 3))
    message("Aggregated columns (sample_group): ", ncol(exprs_smp),
            " | donors: ", length(unique(smp_meta$donor)),
            " | samples: ", length(unique(smp_meta$sample_id)))
  }

  list(
    fit = fit2,
    topTable = tt,
    rank_stat = rank_stat,
    exprs_smp = exprs_smp,
    smp_meta = smp_meta,
    cor = corfit$consensus.correlation
  )
}

# Prep for DA neighborhood markers
milo_obj <- milo_runs_list[["immune"]]$result$milo
milo_da_res <- milo_runs_list[["immune"]]$result$da_res$age_per_dec
milo_da_res <- annotateNhoods(milo_obj, milo_da_res, coldata_col = cell_type_col)

# Define age-depleted and age-enriched macrophage neighborhoods using Milo logFC thresholds.
milo_da_res$NhoodGroup <- ifelse(
  milo_da_res$logFC < -0.25, "neg",
  ifelse(
    milo_da_res$logFC > 0.25, "pos",
    "zero"
  )
)

nh_names <- colnames(nhoods(milo_obj))

keep_rows <- which(
  milo_da_res$cell_type == "Macrophage" &
  milo_da_res$cell_type_fraction > 0.99 &
  milo_da_res$NhoodGroup %in% c("pos", "neg")
)

subset.nhoods <- nh_names[milo_da_res$Nhood[keep_rows]]

table(milo_da_res$NhoodGroup[which(
  milo_da_res$cell_type == "Macrophage" & milo_da_res$cell_type_fraction > 0.99
)])
length(subset.nhoods)

res <- milo_posneg_blocked_DE(
  milo_obj = milo_obj,
  milo_da_res = milo_da_res,
  subset.nhoods = subset.nhoods,  # your 36 seed-ID nhood names
  assay_name = "log1p_norm",
  sample_col = "unique_sample_ID",
  donor_col  = "SenNet_ID",
  groups_use = c("pos","neg"),
  gene.offset = FALSE,
  verbose = TRUE
)

milo_tt_res <- cbind(Gene = rownames(res$topTable), res$topTable)
# Save top table to file
write.csv(milo_tt_res, file = file.path(outs_dir, "milo_macrophage_posneg_DE_topTable.csv"), row.names = FALSE)



#### 5. Correlate Milo and dreamlet DE results ####

# Import dreamlet DGE aging results & subset to Macrophage age contrast
dreamlet_dge_tt <- read.csv(paste0(r_analyses_dir, "/R_DGE/dreamlet_dge_results.csv"))
dreamlet_macrophage_tt <- subset(dreamlet_dge_tt, Celltype == "Macrophage" & contrast == "Age")

# Find common genes between milo_tt_res and dreamlet_dge_res_macrophage
common_genes <- intersect(milo_tt_res$Gene, dreamlet_macrophage_tt$Gene)

# Subset the tables to only these genes
milo_tt_sub <- milo_tt_res[milo_tt_res$Gene %in% common_genes, ]
dreamlet_tt_sub <- dreamlet_macrophage_tt[dreamlet_macrophage_tt$Gene %in% common_genes, ]

# For each gene, match the t statistics
# Ensure the same order of genes
milo_tt_sub <- milo_tt_sub[match(common_genes, milo_tt_sub$Gene), ]
dreamlet_tt_sub <- dreamlet_tt_sub[match(common_genes, dreamlet_tt_sub$Gene), ]

# Extract t statistics
t_milo <- milo_tt_sub$t
t_dreamlet <- dreamlet_tt_sub$t
all.equal(milo_tt_sub$Gene, dreamlet_tt_sub$Gene) # Should be TRUE

# Spearman correlation
spearman_test <- cor.test(t_milo, t_dreamlet, method = "spearman", use = "complete.obs")
spearman_cor <- spearman_test$estimate
print(paste("Spearman correlation (t stats):", round(spearman_cor, 3)))

# Scatter plot
df_corr <- data.frame(
  Gene = common_genes,
  t_milo = t_milo,
  t_dreamlet = t_dreamlet
)

base_size <- 7
p_corr <- ggplot(df_corr, aes(x = t_milo, y = t_dreamlet)) +
  geom_point(alpha = 0.2, color = "#2980b9", shape = 16, size = 1) +  # shape=16 is filled circle, no border
  geom_smooth(method = "lm", color = "black", se = FALSE, linetype = "dashed", linewidth = 0.75) +
  ggtitle("Macrophage Age t-statistics: dreamlet vs Milo") +
  xlab("t-stat (Milo + limma)") +
  ylab("t-stat (dreamlet)") +
  theme_minimal(base_family = "DejaVu Sans", base_size = base_size) +
  theme(
    panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.43),
    axis.title = element_text(colour = "black"),              # all font black
    axis.title.x = element_text(colour = "black", margin = margin(t = 8)),   # add margin to x axis title
    axis.title.y = element_text(colour = "black", margin = margin(r = 8)),   # add margin to y axis title
    axis.text = element_text(colour = "black", size = base_size-2),
    plot.title = element_text(colour = "black", hjust = 0.5, size = base_size, margin = margin(b = 8)),
    axis.ticks = element_line(colour = "black", linewidth = 0.3),                # add ticks
    axis.ticks.length.x = unit(0.15, "cm"),
    axis.ticks.length.y = unit(0.15, "cm")
  ) +
  annotate("text",
           x = min(df_corr$t_milo, na.rm=TRUE),
           y = max(df_corr$t_dreamlet, na.rm=TRUE),
           parse = TRUE,
           label = paste0("rho == ", sprintf("%.2f", spearman_cor)),
           hjust = 0, vjust = 1, size = 3)

ggsave(filename = paste0(outs_dir, "/Milo_Macrophage_dreamlet_correlation.pdf"), plot = p_corr, width = 90, height = 90,
  units = "mm", device = cairo_pdf, dpi = 500)




#### 6. Run pathway analysis using limma camera ####
msigdbr_version <- as.character(utils::packageVersion("msigdbr"))
cat("msigdbr version:", msigdbr_version, "\n")

msigdb_reactome <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
msigdb_list_reactome <- split(msigdb_reactome$gene_symbol, msigdb_reactome$gs_name)

# Choose GSEA gene sets
min_genes <- 5
max_genes <- 500
msigdb_list <- msigdb_list_reactome

# Map gene sets for this expression set & remove gene sets with too few or too many genes
idx <- ids2indices(gene.sets = msigdb_list, identifiers = milo_tt_res$Gene)
idx <- idx[vapply(idx, function(x) length(x) >= min_genes && length(x) <= max_genes, logical(1))]

# Use cameraPR to run camera analysis on the t statistics generated by donor blocked DE
cam_pr <- cameraPR(milo_tt_res$t, index = idx)
cam_pr <- cam_pr[order(cam_pr$FDR), ]
head(cam_pr, 30)

# Save cameraPR results to CSV
write.csv(cam_pr, file = file.path(outs_dir, "Milo_limma_cameraPR_results.csv"), row.names = TRUE)

## Plot top biologically informative pathways
make_panel_macrophage_pathways <- function(
    cam_pr,
    pathways_keep = rev(c(
      "REACTOME_THE_NLRP3_INFLAMMASOME",
      "REACTOME_INTERLEUKIN_1_PROCESSING",
      "REACTOME_INTERLEUKIN_10_SIGNALING",
      "REACTOME_CD163_MEDIATING_AN_ANTI_INFLAMMATORY_RESPONSE",
      "REACTOME_PTK6_PROMOTES_HIF1A_STABILIZATION",
      "REACTOME_ATF6_ATF6_ALPHA_ACTIVATES_CHAPERONE_GENES"
    )),
    pathway_labels = rev(c(
      REACTOME_THE_NLRP3_INFLAMMASOME = "NLRP3\ninflammasome",
      REACTOME_INTERLEUKIN_1_PROCESSING = "IL-1\nprocessing",
      REACTOME_INTERLEUKIN_10_SIGNALING = "IL-10\nsignaling",
      REACTOME_CD163_MEDIATING_AN_ANTI_INFLAMMATORY_RESPONSE = "CD163\nresponse",
      REACTOME_PTK6_PROMOTES_HIF1A_STABILIZATION = "HIF1A\nstabilization",
      REACTOME_ATF6_ATF6_ALPHA_ACTIVATES_CHAPERONE_GENES = "ATF6\nchaperones"
    )),
    base_size = 6,
    font_family = "DejaVu Sans") {
  df <- cam_pr |>
    rownames_to_column("Pathway") |>
    filter(Pathway %in% pathways_keep) |>
    mutate(
      Pathway_label = pathway_labels[Pathway],
      log_fdr = ifelse(Direction == "Up", 1, -1) * -log10(FDR),
      Pathway_label = factor(
        Pathway_label,
        levels = rev(pathway_labels[pathways_keep])
      )
    )

  ggplot(df, aes(y = log_fdr, x = Pathway_label)) +
    geom_hline(yintercept = 0, linewidth = 0.2, color = "black") +
    geom_segment(
      aes(y = 0, yend = log_fdr, x = Pathway_label, xend = Pathway_label),
      linewidth = 0.3, color = "grey60"
    ) +
    geom_point(
      aes(size = NGenes),
      shape = 21, fill = "#e70b25", color = "#222222", stroke = 0.2
    ) +
    scale_size_continuous(
      name = "# Pathway\nGenes",
      range = c(2, 4),
      guide = guide_legend(
        title.position = "top",
        title.hjust = 0.5,
        label.position = "right",
        keyheight = unit(0.6, "cm")
      )
    ) +
    labs(
      y = expression(paste(-log[10], "(FDR)")),
      x = NULL,
      fill = NULL
    ) +
    ylim(0, 3.1) +
    theme_classic(base_size = base_size, base_family = font_family) +
    theme(
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_text(color = "black", size = base_size),
      axis.text.y = element_text(color = "black", size = base_size),
      axis.title.x = element_text(color = "black", size = base_size, margin = margin(t = 8)),
      axis.title.y = element_text(color = "black", size = base_size, margin = margin(r = 8)),
      plot.margin = margin(8, 8, 8, 8),
      legend.title = element_text(color = "black"),
      legend.text = element_text(color = "black", size = base_size),
      axis.ticks = element_line(colour = "black", linewidth = 0.3),
      axis.ticks.length.y = unit(0.15, "cm")
    )
}

p_macrophage_pathways <- make_panel_macrophage_pathways(cam_pr = cam_pr)
ggsave(file.path(outs_dir, "macrophage_pathways_panel.pdf"), p_macrophage_pathways, width = 140, height = 50, units = "mm", device = cairo_pdf)

## Function to plot dotplot of top markers for each pathway
make_marker_theme_dotplot <- function(
    tt_df,
    marker_df = data.frame(
      gene = c("NLRP3","IL1B", "TNFRSF1B",
               "HIF1A", "HSPA5", "VEGFA",
               "IL10","CD163", "MCL1"),
      theme = c("Inflammasome","Inflammasome","Inflammasome",
                "Stress","Stress","Stress",
                "Regulatory","Regulatory","Regulatory"),
      stringsAsFactors = FALSE
    ),
    base_size = 5,
    font_family = "DejaVu Sans",
    point_size_range = c(0.5, 2)
) {

  tvec <- tt_df$t
  names(tvec) <- rownames(tt_df)

  df <- marker_df
  df$tstat <- unname(tvec[df$gene])
  df <- df[!is.na(df$tstat), , drop = FALSE]

  # Order genes within theme by t-stat
  df <- df |>
    dplyr::group_by(theme) |>
    dplyr::arrange(dplyr::desc(tstat), .by_group = TRUE) |>
    dplyr::ungroup()

  gene_levels <- rev(df$gene)
  theme_levels <- unique(df$theme)

  df$gene <- factor(df$gene, levels = gene_levels)
  df$theme <- factor(df$theme, levels = theme_levels)

  L <- max(abs(df$tstat), na.rm = TRUE)

  ggplot(df, aes(x = tstat, y = gene)) +
    geom_vline(xintercept = 0, linewidth = 0.2, color = "black") +
    geom_segment(
      aes(x = 0, xend = tstat, y = gene, yend = gene),
      linewidth = 0.3, color = "grey60"
    ) +
    geom_point(
      shape = 21, color = "black", fill = "#e70b25", stroke = 0.2, size = 1.5
    ) +
    facet_grid(theme ~ ., scales = "free_y", space = "free_y", switch = "y") +
    labs(x = "t-stat", y = NULL) +
    theme_classic(base_size = base_size+1, base_family = font_family) +
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_text(color = "black"),
      axis.text.x = element_text(color = "black"),
      strip.text.y.left = element_blank(),
      panel.spacing.y = grid::unit(0.5, "lines")
    ) +
    xlim(0, 8.2)
}

p_mac_markers <- make_marker_theme_dotplot(tt_df = milo_tt_res)
ggsave(file.path(outs_dir, "macrophage_markers_dotplot.pdf"), p_mac_markers, width = 40, height = 50, units = "mm", device = cairo_pdf)

extract_pathway_gene_tstats <- function(
    tstats,
    pathway_results,
    pathway_list,
    pathways = NULL,
    fdr_max = 0.1,
    keep_only_present = TRUE,
    sort_by = c("Pathway", "tstat_desc")
) {
  sort_by <- match.arg(sort_by)

  # Flatten 1-col matrix/data.frame to named vector
  if (!is.null(dim(tstats))) {
    if (ncol(tstats) != 1) stop("tstats must be a named vector or 1-column matrix/data.frame.")
    tstat_vec <- as.vector(tstats[, 1])
    names(tstat_vec) <- rownames(tstats)
  } else {
    tstat_vec <- tstats
  }

  if (is.null(names(tstat_vec))) stop("tstats must have gene names.")

  # Pick pathways
  pw_df <- pathway_results
  if (!is.null(pathways)) {
    pw_df <- pw_df[pw_df$Pathway %in% pathways, , drop = FALSE]
  } else {
    pw_df <- pw_df[pw_df$FDR <= fdr_max, , drop = FALSE]
  }

  if (nrow(pw_df) == 0) stop("No pathways selected.")

  out_list <- lapply(seq_len(nrow(pw_df)), function(i) {
    pw <- pw_df$Pathway[i]
    if (!pw %in% names(pathway_list)) return(NULL)

    genes <- unique(pathway_list[[pw]])

    df <- data.frame(
      Pathway = pw,
      Gene = genes,
      tstat = unname(tstat_vec[genes]),
      Direction = pw_df$Direction[i],
      FDR = pw_df$FDR[i],
      stringsAsFactors = FALSE
    )

    if (keep_only_present) {
      df <- df[!is.na(df$tstat), , drop = FALSE]
    }

    df
  })

  out <- dplyr::bind_rows(out_list)

  if (nrow(out) == 0) return(out)

  out$gene_direction <- ifelse(out$tstat > 0, "pos", ifelse(out$tstat < 0, "neg", "zero"))
  out$discordant <- ifelse(
    (out$Direction == "Up" & out$tstat < 0) |
      (out$Direction == "Down" & out$tstat > 0),
    TRUE, FALSE
  )

  if (sort_by == "Pathway") {
    out <- out |>
      dplyr::arrange(Pathway, dplyr::desc(abs(tstat)))
  } else {
    out <- out |>
      dplyr::arrange(Pathway, dplyr::desc(tstat))
  }

  out
}

cam_pr_df <- cam_pr |>
  tibble::rownames_to_column("Pathway")

key_pathways <- c(
  "REACTOME_INTERLEUKIN_10_SIGNALING",
  "REACTOME_INTERLEUKIN_4_AND_INTERLEUKIN_13_SIGNALING",
  "REACTOME_CD163_MEDIATING_AN_ANTI_INFLAMMATORY_RESPONSE",
  "REACTOME_THE_NLRP3_INFLAMMASOME",
  "REACTOME_INTERLEUKIN_1_PROCESSING",
  "REACTOME_PTK6_PROMOTES_HIF1A_STABILIZATION",
  "REACTOME_ATF6_ATF6_ALPHA_ACTIVATES_CHAPERONE_GENES"
)

pw_gene_df <- extract_pathway_gene_tstats(
  tstats = res$fit$t,
  pathway_results = cam_pr_df,
  pathway_list = msigdb_list_reactome,
  pathways = key_pathways
)

head(pw_gene_df, 50)
# Save pw_gene_df to a CSV file
write.csv(pw_gene_df, file = file.path(outs_dir, "pw_gene_df.csv"), row.names = FALSE)


#### 7. Export Supplementary Tables ####

# Helper: clean Milo DA results
make_milo_da_table <- function(
    milo_runs_list,
    compartment = "immune",
    contrast_key,
    contrast_label,
    cell_type_col = "cell_type",
    mixed_threshold = 0.70
) {
  milo_obj <- milo_runs_list[[compartment]]$result$milo
  keep <- milo_runs_list[[compartment]]$result$nhoods_kept
  da_res <- milo_runs_list[[compartment]]$result$da_res[[contrast_key]]

  if (is.null(da_res)) {
    stop("No DA results found for contrast_key: ", contrast_key)
  }

  da_annot <- annotateNhoods(
    milo_obj,
    da_res,
    coldata_col = cell_type_col,
    subset.nhoods = keep
  )

  # Add nhood size
  nhood_sizes <- Matrix::colSums(nhoods(milo_obj) > 0)
  da_annot$nhood_size <- as.numeric(nhood_sizes[da_annot$Nhood])

  da_annot |>
    dplyr::mutate(
      contrast = contrast_label,
      nhood = .data[["Nhood"]],
      kept_for_testing = keep[.data[["Nhood"]]],
      nhood_size = as.numeric(nhood_sizes[.data[["Nhood"]]]),
      cell_type_original = .data[[cell_type_col]],
      cell_type = ifelse(
        .data[["cell_type_fraction"]] < mixed_threshold,
        "Mixed",
        .data[["cell_type_original"]]
      ),
      significant_sFDR_0_1 = .data[["SpatialFDR"]] < 0.1
    ) |>
    dplyr::rename(
      p_value = PValue,
      spatial_fdr = SpatialFDR,
      converged = Converged,
      dominant_cell_type_fraction = cell_type_fraction
    ) |>
    dplyr::select(
      contrast,
      nhood,
      kept_for_testing,
      nhood_size,
      cell_type,
      cell_type_original,
      dominant_cell_type_fraction,
      logFC,
      logCPM,
      SE,
      tvalue,
      p_value,
      spatial_fdr,
      significant_sFDR_0_1,
      converged,
      -dplyr::any_of(c("Nhood", "SenNet_ID_variance", "Logliklihood"))
    )
}

# Helper: neighborhood QC table
make_milo_nhood_qc_table <- function(
    milo_runs_list,
    compartment = "immune",
    donor_col = "SenNet_ID",
    sample_col = "unique_sample_ID",
    cell_type_col = "cell_type"
) {
  milo_obj <- milo_runs_list[[compartment]]$result$milo
  keep <- milo_runs_list[[compartment]]$result$nhoods_kept

  # Neighborhood x sample counts
  cts <- nhoodCounts(milo_obj)
  sample_ids <- colnames(cts)
  n_total_samples <- length(sample_ids)

  cd <- as.data.frame(colData(milo_obj))

  sample_donor_df <- cd %>%
    select(all_of(c(sample_col, donor_col))) %>%
    distinct()

  donor <- sample_donor_df[[donor_col]][match(sample_ids, sample_donor_df[[sample_col]])]

  n_total_donors <- length(unique(donor))

  sample_nonzero <- Matrix::rowSums(cts > 0)

  donor_levels <- unique(donor)
  cts_donor <- sapply(donor_levels, function(d) {
    Matrix::rowSums(cts[, donor == d, drop = FALSE])
  })
  cts_donor <- as.matrix(cts_donor)

  donor_nonzero <- rowSums(cts_donor > 0)

  total_nhood_cells <- rowSums(cts_donor)
  total_nhood_cells_safe <- total_nhood_cells
  total_nhood_cells_safe[total_nhood_cells_safe == 0] <- 1

  top3_donor_cells <- apply(cts_donor, 1, function(x) {
    sum(sort(x, decreasing = TRUE)[seq_len(min(3, length(x)))])
  })

  top3_donor_fraction <- top3_donor_cells / total_nhood_cells_safe

  # Dominant cell type and purity for each neighborhood
  nh_mat <- nhoods(milo_obj) > 0
  cell_types <- as.character(colData(milo_obj)[[cell_type_col]])

  dominant_cell_type <- character(ncol(nh_mat))
  dominant_cell_type_fraction <- numeric(ncol(nh_mat))

  for (i in seq_len(ncol(nh_mat))) {
    idx <- which(nh_mat[, i])
    tab <- sort(table(cell_types[idx]), decreasing = TRUE)
    dominant_cell_type[i] <- names(tab)[1]
    dominant_cell_type_fraction[i] <- as.numeric(tab[1]) / sum(tab)
  }

  tibble::tibble(
    nhood = seq_len(ncol(nh_mat)),
    kept_for_testing = keep,
    total_nuclei = as.numeric(total_nhood_cells),
    n_samples_with_detected_nuclei = as.numeric(sample_nonzero),
    percent_samples_with_detected_nuclei = round(100 * as.numeric(sample_nonzero) / n_total_samples, 2),
    n_donors_with_detected_nuclei = as.numeric(donor_nonzero),
    percent_donors_with_detected_nuclei = round(100 * as.numeric(donor_nonzero) / n_total_donors, 2),
    top3_donor_fraction = round(top3_donor_fraction, 4),
    top3_donor_percent = round(100 * top3_donor_fraction, 2),
    dominant_cell_type = dominant_cell_type,
    dominant_cell_type_fraction = round(dominant_cell_type_fraction, 4)
  )
}

# 1. Neighborhood QC
milo_nhood_qc <- make_milo_nhood_qc_table(
  milo_runs_list = milo_runs_list,
  compartment = "immune",
  donor_col = donor_col,
  sample_col = sample_col,
  cell_type_col = cell_type_col
)

# 2. Combined Milo DA table
milo_da_age <- make_milo_da_table(
  milo_runs_list = milo_runs_list,
  compartment = "immune",
  contrast_key = "age_per_dec",
  contrast_label = "Age_per_decade",
  cell_type_col = cell_type_col
)

milo_da_parity <- make_milo_da_table(
  milo_runs_list = milo_runs_list,
  compartment = "immune",
  contrast_key = "para_bin",
  contrast_label = "Parity_1-3_vs_0",
  cell_type_col = cell_type_col
)

milo_da_results <- bind_rows(
  milo_da_age,
  milo_da_parity
) %>%
  arrange(contrast, nhood)

# 3. Cell type-level DA summary
milo_celltype_DA_summary <- milo_da_results %>%
  filter(
    kept_for_testing,
    converged,
    !is.na(logFC),
    cell_type != "Mixed"
  ) %>%
  group_by(contrast, cell_type) %>%
  summarise(
    n_neighborhoods = n(),
    n_significant_sFDR_0_1 = sum(spatial_fdr < 0.1, na.rm = TRUE),
    n_positive_sFDR_0_1 = sum(spatial_fdr < 0.1 & logFC > 0, na.rm = TRUE),
    n_negative_sFDR_0_1 = sum(spatial_fdr < 0.1 & logFC < 0, na.rm = TRUE),
    median_logFC = median(logFC, na.rm = TRUE),
    mean_logFC = mean(logFC, na.rm = TRUE),
    min_logFC = min(logFC, na.rm = TRUE),
    max_logFC = max(logFC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(contrast, desc(n_significant_sFDR_0_1), cell_type)

# 4. Macrophage pos/neg DE table
milo_macrophage_posneg_DE <- milo_tt_res %>%
  as.data.frame() %>%
  dplyr::rename(
    gene = Gene,
    p_value = P.Value,
    adj_p_value = adj.P.Val,
    B_statistic = B
  ) %>%
  select(
    gene,
    logFC,
    AveExpr,
    t,
    p_value,
    adj_p_value,
    B_statistic,
    everything()
  ) %>%
  arrange(p_value)

# 5. Macrophage pathway table
milo_macrophage_pathways <- cam_pr %>%
  as.data.frame() %>%
  rownames_to_column("pathway") %>%
  dplyr::rename(
    n_genes = NGenes,
    p_value = PValue,
    fdr = FDR,
    direction = Direction
  ) %>%
  select(
    pathway,
    n_genes,
    direction,
    p_value,
    fdr,
    everything()
  ) %>%
  arrange(fdr)

# 6. Selected macrophage pathway genes
milo_macrophage_pathway_genes <- pw_gene_df %>%
  as.data.frame() %>%
  dplyr::rename(
    pathway = Pathway,
    gene = Gene,
    t_statistic = tstat,
    pathway_direction = Direction,
    pathway_fdr = FDR
  ) %>%
  select(
    pathway,
    gene,
    t_statistic,
    pathway_direction,
    pathway_fdr,
    gene_direction,
    discordant
  ) %>%
  arrange(pathway, desc(abs(t_statistic)))

# Save individual CSVs
write.csv(milo_nhood_qc, file.path(outs_dir, "Milo_nhood_QC.csv"), row.names = FALSE)
write.csv(milo_da_results, file.path(outs_dir, "Milo_DA_results.csv"), row.names = FALSE)
write.csv(milo_celltype_DA_summary, file.path(outs_dir, "Milo_celltype_DA_summary.csv"), row.names = FALSE)
write.csv(milo_macrophage_posneg_DE, file.path(outs_dir, "Milo_macrophage_posneg_DE.csv"), row.names = FALSE)
write.csv(milo_macrophage_pathways, file.path(outs_dir, "Milo_macrophage_pathways.csv"), row.names = FALSE)
write.csv(milo_macrophage_pathway_genes, file.path(outs_dir, "Milo_macrophage_pathway_genes.csv"), row.names = FALSE)

# Save combined Supplementary Table workbook
milo_supp_tables <- list(
  Milo_nhood_QC = milo_nhood_qc,
  Milo_DA_results = milo_da_results,
  Milo_celltype_DA_summary = milo_celltype_DA_summary,
  Milo_macrophage_posneg_DE = milo_macrophage_posneg_DE,
  Milo_macrophage_pathways = milo_macrophage_pathways,
  Milo_pathway_gene_tstats = milo_macrophage_pathway_genes
)

openxlsx::write.xlsx(
  milo_supp_tables,
  file = file.path(outs_dir, "Supplementary_Table_X_Milo_immune_results.xlsx"),
  overwrite = TRUE
)


#### 8. Session Info ####
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
# [1] stats4    stats     graphics  grDevices utils     datasets  methods   base

# other attached packages:
#  [1] msigdbr_25.1.1              openxlsx_4.2.8.1            httpgd_2.0.4                vegan_2.7-1
#  [5] permute_0.9-8               ggeffects_2.3.1             patchwork_1.3.1             lmerTest_3.1-3
#  [9] lme4_1.1-37                 FNN_1.1.4.1                 scater_1.34.1               scuttle_1.16.0
# [13] miloR_2.6.0                 edgeR_4.6.3                 limma_3.64.1                SingleCellExperiment_1.30.1
# [17] SummarizedExperiment_1.38.1 Biobase_2.68.0              GenomicRanges_1.60.0        GenomeInfoDb_1.44.0
# [21] IRanges_2.42.0              S4Vectors_0.48.0            BiocGenerics_0.54.0         generics_0.1.4
# [25] MatrixGenerics_1.20.0       matrixStats_1.5.0           scales_1.4.0                ggbeeswarm_0.7.2
# [29] lubridate_1.9.4             forcats_1.0.0               stringr_1.5.1               dplyr_1.1.4
# [33] purrr_1.1.0                 readr_2.1.5                 tidyr_1.3.1                 tibble_3.3.0
# [37] ggplot2_3.5.2               tidyverse_2.0.0             Matrix_1.7-3                zellkonverter_1.18.0

# loaded via a namespace (and not attached):
#  [1] Rdpack_2.6.4            gridExtra_2.3           rlang_1.1.7             magrittr_2.0.4
#  [5] compiler_4.4.1          mgcv_1.9-3              dir.expiry_1.16.0       png_0.1-8
#  [9] systemfonts_1.2.3       vctrs_0.7.1             pkgconfig_2.0.3         crayon_1.5.3
# [13] fastmap_1.2.0           XVector_0.48.0          ggraph_2.2.2            tzdb_0.5.0
# [17] pracma_2.4.6            nloptr_2.2.1            UCSC.utils_1.4.0        cachem_1.1.0
# [21] beachmat_2.24.0         jsonlite_2.0.0          DelayedArray_0.34.1     BiocParallel_1.42.1
# [25] tweenr_2.0.3            cluster_2.1.8.1         irlba_2.3.5.1           parallel_4.4.1
# [29] R6_2.6.1                stringi_1.8.7           RColorBrewer_1.1-3      reticulate_1.42.0
# [33] boot_1.3-32             numDeriv_2016.8-1.1     assertthat_0.2.1        Rcpp_1.1.1
# [37] splines_4.4.1           igraph_2.1.4            timechange_0.3.0        tidyselect_1.2.1
# [41] abind_1.4-8             viridis_0.6.5           codetools_0.2-20        curl_7.0.0
# [45] lattice_0.22-7          basilisk.utils_1.20.0   withr_3.0.2             unigd_0.1.3
# [49] polyclip_1.10-7         zip_2.3.3               pillar_1.11.0           filelock_1.0.3
# [53] insight_1.4.5           reformulas_0.4.1        hms_1.1.3               minqa_1.2.8
# [57] gtools_3.9.5            glue_1.8.0              tools_4.4.1             BiocNeighbors_2.0.1
# [61] ScaledMatrix_1.14.0     locfit_1.5-9.12         babelgene_22.9          graphlayouts_1.2.2
# [65] tidygraph_1.3.1         cowplot_1.2.0           grid_4.4.1              rbibutils_2.3
# [69] nlme_3.1-168            GenomeInfoDbData_1.2.14 basilisk_1.20.0         beeswarm_0.4.0
# [73] BiocSingular_1.22.0     ggforce_0.5.0           vipor_0.4.7             cli_3.6.5
# [77] rsvd_1.0.5              S4Arrays_1.8.1          viridisLite_0.4.2       gtable_0.3.6
# [81] SparseArray_1.8.0       ggrepel_0.9.6           farver_2.1.2            memoise_2.0.1
# [85] lifecycle_1.0.5         httr_1.4.7              statmod_1.5.0           MASS_7.3-65





