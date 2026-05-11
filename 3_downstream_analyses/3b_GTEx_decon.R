## Script name: Human Ovary Atlas - GTEx ovary bulk RNA-seq deconvolution with InstaPrism
##
## Purpose: Deconvolve GTEx postmenopausal ovary bulk RNA-seq data using an ovary snucRNA-seq reference augmented with
## a Tabula Sapiens erythrocyte reference, then assess deconvolved cell-type profiles and export cell-type proportion summaries.
##
## Author: Josef Byrne
#########################


#### 1. Import libraries and set directories ####

# Import libraries
library(readr)
library(readxl)
library(naniar)
library(BayesPrism)
library(InstaPrism) # version v0.1.6 via Github
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
library(Polychrome)
library(UCell)
library(pheatmap)

# Initialize paths
nostromo_datadir <- "/path/to/project_directory"
r_analyses_dir <- paste0(nostromo_datadir, "/R_analyses")
adata_objdir <- paste0(nostromo_datadir, "/python_analyses/adata_objs")
gtex_path <- paste0(r_analyses_dir, "/GTEx_ovary")
tabula_sapiens_path <- "/path/to/tabula_sapiens_ref"
outs_dir <- paste0(r_analyses_dir, "/deconv_results")



#### 2. Prepare GTEx bulk RNA-seq data ####

## Import bulk RNA-seq ovary data & format
gtex_bulkrna <- read.delim(file = paste0(gtex_path, "/gene_reads_v10_ovary.gct.gz"), skip = 2, check.names = FALSE)
gene_id <- gtex_bulkrna$Description
gtex_bulkrna <- as.matrix(gtex_bulkrna[, -(1:2)])
rownames(gtex_bulkrna) <- gene_id
storage.mode(gtex_bulkrna) <- "numeric"
gtex_bulkrna <- t(gtex_bulkrna)

# Get current rownames (sample IDs)
rna_sample_ids <- rownames(gtex_bulkrna)

# Extract subject IDs from sample names from bulk RNA-seq data & make rownames
rownames(gtex_bulkrna) <- sub("^([^-]+)-([^-]+)-.*$", "\\1-\\2", rna_sample_ids)

## Import & align sample metadata
gtex_sample_dd <- as.data.frame(read_excel(
  paste0(gtex_path, "/GTEx_Analysis_2022-06-06_v10_Annotations_GTEx_Analysis_2022-06-06_v10_Annotations_SampleAttributesDD.xlsx")))
gtex_sample_ds <- read.csv(
  paste0(gtex_path, "/GTEx_Analysis_2022-06-06_v10_Annotations_GTEx_Analysis_2022-06-06_v10_Annotations_SampleAttributesDS.txt"), sep = '\t')

# Extract subject IDs from sample names from bulk RNA-seq data & make rownames
gtex_sample_ds$SUBJID <- sub("^([^-]+)-([^-]+)-.*$", "\\1-\\2", gtex_sample_ds$SAMPID)

## Import & align subject metadata
gtex_subject_dd <- as.data.frame(read_excel(
  paste0(gtex_path, "/GTEx_Analysis_2022-06-06_v10_Annotations_GTEx_Analysis_2022-06-06_v10_Annotations_SubjectPhenotypesDD.xlsx")))
gtex_subject_ds <- read.csv(
  paste0(gtex_path, "/GTEx_Analysis_2022-06-06_v10_Annotations_GTEx_Analysis_2022-06-06_v10_Annotations_SubjectPhenotypesDS.txt"), sep = '\t')

# Align by sample and subject data on SUBJID & subset to only include subjects with bulk RNA-seq data
sample_metadata <- merge(gtex_subject_ds, gtex_sample_ds, by = "SUBJID", all.x = TRUE, suffixes = c("_subject", "_sample"))
sample_metadata <- sample_metadata[sample_metadata$SAMPID %in% rna_sample_ids, ]

## Subset to include likely postmenopausal donors
# Subset to postmenopausal donors: either "post" and "meno" both in SMPTHNTS, OR AGE >= 55

# Identify subjects with explicit SMPTHNTS post-meno designation
has_post_meno <- !is.na(sample_metadata$SMPTHNTS) &
                 grepl("post", tolower(sample_metadata$SMPTHNTS), fixed = TRUE) &
                 grepl("meno", tolower(sample_metadata$SMPTHNTS), fixed = TRUE)
n_post_meno <- sum(has_post_meno, na.rm = TRUE)
cat("Number with explicit SMPTHNTS post-meno designation:", n_post_meno, "\n")
# Number with explicit SMPTHNTS post-meno designation: 52

# Identify subjects with AGE >= 55 and NOT explicit post-meno label, by 55 a significant majority of women in the USA are postmenopausal
age_55plus <- as.numeric(sample_metadata$AGE) >= 55
age_55plus_no_postmeno <- sum(age_55plus & !has_post_meno, na.rm = TRUE)
cat("Number with AGE >= 55 without explicit post-meno label:", age_55plus_no_postmeno, "\n")
# Number with AGE >= 55 without explicit post-meno label: 58

# Subset to postmenopausal: either explicit designation + age >=40    OR    age >= 55
sample_metadata <- sample_metadata[(has_post_meno & as.numeric(sample_metadata$AGE) >= 40) | age_55plus,]

## Remove sparse levels: remove n=1 of RACE=1, and remove n=1 of SMCENTER=D1, won't be included in downstream analyses
sample_metadata <- sample_metadata[sample_metadata$RACE != 1, ]
sample_metadata <- sample_metadata[sample_metadata$SMCENTER != "D1", ]

cat("Number of subjects after subsetting:", nrow(sample_metadata), "\n")
# Number of subjects after subsetting: 106

# Subset bulk RNA-seq data to only include subjects in sample_metadata
gtex_bulkrna <- gtex_bulkrna[rownames(gtex_bulkrna) %in% sample_metadata$SUBJID, ]




#### 3. Prepare ovary snucRNA-seq and Tabula Sapiens references ####

## Prepare ovary snucRNA-seq and Tabula Sapiens erythrocyte references
sce_obj_snuc <- readH5AD(paste0(adata_objdir, "/final_annotated_atlas.h5ad"), reader = "R", verbose = TRUE)
sce_obj_ts <- readH5AD(paste0(tabula_sapiens_path, "/tabula_sapiens_erythrocyte.h5ad"), reader = "R", verbose = TRUE)

# Rename snuc- levels for readability
colData(sce_obj_snuc)$cell_type_fine <- forcats::fct_recode(
  colData(sce_obj_snuc)$cell_type_fine,
  "Activated Stromal" = "Activated (mito-high) Stromal",
  "IER Stromal" = "Immediate-early response Stromal",
  "ArtEnd" = "Arterial",
  "LympEnd" = "Lymphatic",
  "VenEnd" = "Venous",
  "CapEnd" = "Capillary",
  "nmSchwann" = "Nonmyelinating Schwann"
)

## Refine cell type annotations to be forwarded to deconvolution & downsample to balance donor contribution
cd <- as.data.frame(colData(sce_obj_snuc))
cd$bp_cell_type  <- as.character(cd$cell_type_fine)

# Remove noisy or weakly distinguishable cell states that were poorly identifiable in bulk RNA-seq deconvolution
cd <- cd[cd$bp_cell_type != "Stress-response Secretory Epithelial", ]
cd <- cd[cd$bp_cell_type != "Activated Stromal", ]
cd <- cd[cd$bp_cell_type != "IER Stromal", ]
cd <- cd[cd$bp_cell_type != "TNF-activated Stromal", ]
cd <- cd[cd$bp_cell_type != "Stressed Perivascular", ]
cd <- cd[cd$bp_cell_type != "NK", ] # Removed due to donor-specific, weakly distinguishable profile in bulk deconvolution

# Rare cell types with cell counts between 10-20 (Mast, B Cell, Myelinating Schwann) are very transcriptionally distinct from
# other cell types and are not expected to impact quality of deconvolution and thus are retained.

## Downsample snuc- reference to balance donor contribution with preserving cell type representation
cap_per_donor <- 100
set.seed(42)
cd$cell_id <- rownames(cd)

keep_cells <- cd %>%
  group_by(SenNet_ID, bp_cell_type) %>%
  slice_sample(n = cap_per_donor) %>%
  ungroup() %>%
  pull(cell_id)

cd <- cd[cd$cell_id %in% keep_cells, ]

# Subset the sce_obj to only include the cells in cd
sce_obj_snuc_bp <- sce_obj_snuc[, rownames(cd)]
# Add the refined cell type annotations to the sce_obj
colData(sce_obj_snuc_bp)$bp_cell_type  <- factor(cd$bp_cell_type)

## Downsample Tabula Sapiens reference to balance donor contribution with preserving cell type representation
cd_ts <- as.data.frame(colData(sce_obj_ts))
cd_ts$cell_type <- as.character(cd_ts$cell_type)
cd_ts$cell_id <- rownames(cd_ts)

keep_cells_ts <- cd_ts %>%
  group_by(donor_id, cell_type) %>%
  slice_sample(n = 100) %>%
  ungroup() %>%
  pull(cell_id)

cd_ts <- cd_ts[cd_ts$cell_id %in% keep_cells_ts, ]

# Subset the Tabula Sapiens object to retained erythrocyte-reference cells.
sce_obj_ts_bp <- sce_obj_ts[, rownames(cd_ts)]
# Add the refined cell type annotations to the Tabula Sapiens reference object
colData(sce_obj_ts_bp)$bp_cell_type <- factor(cd_ts$cell_type)

# Make first char uppercase for bp_cell_type for Tabula Sapiens reference
colData(sce_obj_ts_bp)$bp_cell_type <- factor(stringr::str_replace(colData(sce_obj_ts_bp)$bp_cell_type, "^(\\w)", toupper))




#### 4. Merge references and align genes ####

assay_ts   <- "decontXcounts"
assay_snuc <- "raw_counts"

# helper: minimal SCE with just counts + bp_cell_type
min_sce <- function(X, cell_type) {
  SingleCellExperiment(
    assays  = list(counts = as(X, "dgCMatrix")),
    colData = DataFrame(bp_cell_type = factor(cell_type))
  )
}

# Tabula: ENSG -> symbol collapse (sum) into counts matrix
sym <- as.character(rowData(sce_obj_ts_bp)$feature_name)
ok  <- !is.na(sym) & sym != ""

ts_counts_sym <- rowsum(
  as(assay(sce_obj_ts_bp[ok, ], assay_ts), "dgCMatrix"),
  group = sym[ok],
  reorder = TRUE
)

ts_min <- min_sce(ts_counts_sym, colData(sce_obj_ts_bp[ok, ])$bp_cell_type)

# snuc: grab counts matrix
sn_counts <- as(assay(sce_obj_snuc_bp, assay_snuc), "dgCMatrix")
sn_min <- min_sce(sn_counts, colData(sce_obj_snuc_bp)$bp_cell_type)
rownames(sn_min) <- rownames(sce_obj_snuc_bp)  # keep your symbol rownames

# align genes + merge
g <- intersect(rownames(sn_min), rownames(ts_min))

colData(sn_min)$source <- "snuc"
colData(ts_min)$source <- "TabulaSapiens"
sce_ref_merged <- cbind(sn_min[g, ], ts_min[g, ])



#### 5. Build BayesPrism reference matrix and run InstaPrism deconvolution ####

# Note: only cell types are included for reference deconvolution, not cell states which are weakly identifiable in this context and result in
# states that serve as a "sink" for cell proportion attribution and yield inaccurate deconvolution results.
sc.data <- t(as.matrix(assay(sce_ref_merged, "counts")))
bk.data <- gtex_bulkrna
celltypes <- c(as.character(colData(sce_ref_merged)$bp_cell_type))

# Subset bk.data and sc.data to their common genes (colnames)
common_genes <- intersect(colnames(sc.data), colnames(bk.data))
sc.data <- sc.data[, common_genes, drop = FALSE]
bk.data <- bk.data[, common_genes, drop = FALSE]

# Filter noisy/low specificity genes & retain only protein coding genes
sc.data.filtered <- cleanup.genes(input = sc.data, input.type = "count.matrix", species = "hs",
 gene.group = c("Rb", "Mrp", "other_Rb", "chrM", "MALAT1", "chrX", "chrY"), exp.cells = 5)

sc.data.filtered.pc <- select.gene.type(sc.data.filtered, gene.type = "protein_coding")

# Create BayesPrism reference matrix and run InstaPrism deconvolution, using default parameters.
# Update reference + re-run deconvolution results using updated reference to get updated posterior with improved deconvolution accuracy.
prism_obj <- new.prism(reference = sc.data.filtered.pc, mixture = bk.data, input.type = "count.matrix",
  cell.type.labels = celltypes, cell.state.labels = NULL, key = NULL, outlier.cut = 0.01, outlier.fraction = 0.1)

set.seed(42)
deconv_res <- InstaPrism(input_type = "prism", prismObj = prism_obj, n.core = 100, verbose = TRUE) %>%
  InstaPrism_update(bulk_Expr = t(bk.data), cell.types.to.update = "all", key = NA, n.core = 100)

# Save full deconvolution results object to file
saveRDS(deconv_res, file = file.path(outs_dir, "deconv_res_full.rds"))

# Extract cell-type proportions: rows are GTEx samples, columns are cell types
celltype_theta <- deconv_res@theta



#### 6. Assess quality of deconvolved cell-type profiles ####

celltype_order_by_compartment <- c(
  "Homeostatic Fibroblast", "Steroidogenic Stromal", "Contractile Stromal", "Perivascular-like Stromal", "HOXD-patterned Myofibroblast", "Myelinating Schwann", "nmSchwann",
  "ArtEnd", "VenEnd", "CapEnd", "LympEnd", "Lymphatic Valve Up", "ArtSMC", "VenSMC", "SpecSMC", "Pericyte", "Mesothelium", "PAX2+ Mullerian Epithelial", "PAX2- Mullerian Epithelial",
  "CD4+ T", "CD8+ T", "CD16+ NK", "CD16- NK", "B Cell", "Plasma B", "Macrophage", "CD16+ Monocyte", "cDC1", "cDC2", "Mast", "Erythrocyte", "Cyc Immune"
)

## Extract top marker genes per deconvoluted cell types and take mean across samples
Z <- get_Z_array(deconv_res, resolution = 'ct', n.core = 100)
Z_state_gene_mean <- apply(Z, c(3, 2), mean)
write.csv(Z_state_gene_mean, file = paste0(outs_dir, "/Z_state_gene_mean_deconv_instaprism_raw.csv"), row.names = TRUE)

rank_deconv_profile_markers <- function(M_ct_gene, pseudo = 0.1, top_n = 25, min_mean = 0) {
  M <- as.matrix(M_ct_gene)

  out <- lapply(rownames(M), function(ct) {
    x <- M[ct, ]
    bg <- colMeans(M[setdiff(rownames(M), ct), , drop = FALSE])

    log2fc <- log2((x + pseudo) / (bg + pseudo))

    keep <- x >= min_mean
    log2fc2 <- log2fc[keep]

    ord <- order(log2fc2, decreasing = TRUE)
    genes <- names(log2fc2)[ord][seq_len(min(top_n, length(ord)))]

    data.frame(
      deconvolved_cell_type = ct,
      rank = seq_along(genes),
      gene = genes,
      log2FC_vs_other_deconv_profiles = as.numeric(log2fc[genes]),
      mean_expression_in_deconv_profile = as.numeric(x[genes]),
      mean_expression_in_other_profiles = as.numeric(bg[genes]),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  })

  do.call(rbind, out)
}

deconv_profile_top_markers <- rank_deconv_profile_markers(Z_state_gene_mean, pseudo = 0.1, top_n = 50, min_mean = 0)

deconv_profile_top_markers <- deconv_profile_top_markers %>%
  dplyr::arrange(
    deconvolved_cell_type,
    rank
  )

write.csv(deconv_profile_top_markers, file = file.path(outs_dir, "deconv_profile_top50_markers_instaprism.csv"),
  row.names = FALSE)

## Check for atlas annotation marker overlap via scoring to ensure consistency with annotated cell type marker genes
sig_list_atlas <- list(
  # ---- Stromal / fibro ----
  Fibro_ECM = c("DCN","PDGFRA"),
  Stromal1 = c("STAR","KLHDC8A","SIGLEC11"),
  Steroidogenic = c("CYP11A1","CYP17A1","FDX1","DHCR24","DLK1"),
  HOXD_Myofibro = c("TWIST2","HAND2","HOXD11","WNT2","DLK1"),

  # ---- Schwann / nerve ----
  Schwann = c("S100B","MPZ","PLP1","MBP"),

  # ---- Endothelium ----
  Endothelial_core = c("PECAM1","VWF","CDH5", "CLDN5", "ESAM"),
  Arterial_EC = c("SEMA3G","GJA5","HEY1", "SOX17", "NOTCH4"),
  Venous_EC = c("ACKR1","PLVAP","PLAT"),
  Capillary_EC = c("CA4","MFSD2A","SLC7A5"),
  Lymphatic_EC = c("PROX1","LYVE1","CCL21", "TFF3"),

  # ---- Mural / perivascular ----
  SMC_core = c("MYH11","ACTA2","TAGLN","CNN1","DES"),
  Art_SMC = c("RERGL","PDE4C","RGS6"),
  Ven_SMC = c("HMCN2","PDE4D","DES"),
  Spec_SMC = c("CNN1","COL16A1","SYTL4","GREB1","PRUNE2"),
  Pericyte = c("RGS5","NOTCH3","ADGRF5","BCAM"),

  # ---- Surface / epithelial ----
  General_epithelial = c('KRT18', 'MSLN', 'KRT8', 'KRT19'),
  Mesothelium = c('CALB2', 'BNC1', 'SFRP1', "WT1", "PDPN"),
  Mullerian_PAX2pos = c("PAX8","PAX2","EPCAM","KRT8","KRT18","KRT19","CLDN3","CLDN4","MUC1","KRT7"),
  Mullerian_PAX2neg = c("PAX8", "OCLN","FGFR3","EPHA7","ADGRG6","PLA2G4C"),

  # ---- Immune ----
  Tcell = c("IL7R","CD3D", "TRAC"),
  CD4_T = c("CD4", "CD40LG", "CCR6"),
  CD8_T = c("CD8A","KLRK1","GZMK","CCL5"),
  NK = c("NKG7","PRF1","GNLY","GZMB","NCAM1","KLRC3"),
  Bcell = c("MS4A1","BANK1","FCRL1", "FCRL2"),
  Plasma = c("JCHAIN", "XBP1","IGHG1", "MZB1", "PDK1"),
  Macrophage = c("CD163","CD14","C1QA","F13A1","LYVE1"),
  Monocyte_nonclassical = c("VCAN","FCGR3A","LILRA5"),
  cDC1 = c("FLT3","CLEC9A","XCR1","CADM1", "IRF8", "IDO1"),
  cDC2 = c("FLT3", "CD1C", "CD1D", 'CLEC10A', 'SIRPA', "CIITA", "ITGAX"),
  Mast = c("MS4A2","KIT", "IL1RL1"),

  # “Contaminant / QC”
  Erythrocyte = c("HBB","HBA1","HBA2","ALAS2","AHSP","SLC4A1","GYPA","CA1"),

  # Cycling state program
  Cycling = c("MKI67","TOP2A","CCNA2","FOXM1","DLGAP5","TPX2","CENPF","BUB1B")
)

# Get celltype x gene matrix & log-transform to stabilize scores
mat <- t(log1p(as.matrix(Z_state_gene_mean)))   # genes x celltypes

# Keep only genes present scoring list
sig_list2 <- lapply(sig_list_atlas, function(g) intersect(g, rownames(mat)))

# Score celltypes for each signature via UCell
scores <- ScoreSignatures_UCell(mat, features = sig_list2)
scores_mat <- as.matrix(scores)
scores_mat_ordered <- scores_mat[celltype_order_by_compartment, ]

supp_sig_list_atlas_plot <- pheatmap(t(scores_mat_ordered),
  cluster_rows = FALSE, cluster_cols = FALSE,
  color = colorRampPalette(c("white", "#FCA5A5", "#EF4444", "#DC2626"))(100), breaks = seq(0, 1, length.out = 101),
  border_color = "grey85", fontsize = 7, fontsize_row = 7, fontsize_col = 7, angle_col = 90,
  legend_breaks = c(0, 0.25, 0.5, 0.75, 1.0), legend_labels = c("0", "0.25", "0.5", "0.75", "1.0")
)

ggsave(filename = paste0(outs_dir, "/deconv_supp_atlas_sig_heatmap.pdf"), plot = supp_sig_list_atlas_plot, width = 180,
  height = 140, units = "mm", device = cairo_pdf, dpi = 500)

write.csv(scores_mat_ordered, file = paste0(outs_dir, "/deconv_supp_scores_atlas.csv"), row.names = TRUE)

## Update annotations for Z and theta based on marker review
updated_annotations <- list(
  SpecSMC = "Stress-signaling Mural",
  cDC2 = "Granulocyte-myeloid Sink",
  `Lymphatic Valve Up` = "Lymphatic Valve (ambiguous)",
  nmSchwann = "Schwann/Nerve",
  `Myelinating Schwann` = "Schwann/Nerve"
)

# Function to update rownames based on a mapping and combine specified duplicated rows by mean
update_and_combine_rows <- function(mat, annotation_map, combine_names = NULL, combine_fun = c("mean", "sum")) {
  combine_fun <- match.arg(combine_fun)

  # Ensure matrix
  mat <- as.matrix(mat)

  # Update rownames via mapping
  rn <- rownames(mat)
  idx <- match(rn, names(annotation_map))
  rn2 <- rn
  rn2[!is.na(idx)] <- unname(annotation_map[idx[!is.na(idx)]])
  rownames(mat) <- rn2

  # Combine duplicates for each target name
  for (nm in combine_names) {
    rows_to_combine <- which(rownames(mat) == nm)
    if (length(rows_to_combine) > 1) {

      if (combine_fun == "mean") {
        combined_row <- colMeans(mat[rows_to_combine, , drop = FALSE], na.rm = TRUE)
      } else {
        combined_row <- colSums(mat[rows_to_combine, , drop = FALSE], na.rm = TRUE)
      }

      # Remove old rows
      mat2 <- mat[-rows_to_combine, , drop = FALSE]

      # Make a 1xN matrix with rowname nm
      new_row <- matrix(combined_row, nrow = 1,
                        dimnames = list(nm, colnames(mat2)))

      mat <- rbind(mat2, new_row)
    }
  }
  mat
}

# Update annotations of Z_state_gene_mean and celltype_theta
Z_state_gene_mean <- update_and_combine_rows(Z_state_gene_mean, updated_annotations, combine_names = c("Schwann/Nerve"), combine_fun = "mean")
celltype_theta    <- update_and_combine_rows(celltype_theta,    updated_annotations, combine_names = c("Schwann/Nerve"), combine_fun = "sum")

celltype_order_by_compartment <- c(
  "Homeostatic Fibroblast", "Steroidogenic Stromal", "Contractile Stromal", "Perivascular-like Stromal", "HOXD-patterned Myofibroblast", "Schwann/Nerve",
  "ArtEnd", "VenEnd", "CapEnd", "LympEnd", "Lymphatic Valve (ambiguous)", "ArtSMC", "VenSMC", "Stress-signaling Mural", "Pericyte", "Mesothelium", "PAX2+ Mullerian Epithelial", "PAX2- Mullerian Epithelial",
  "CD4+ T", "CD8+ T", "CD16+ NK", "CD16- NK", "B Cell", "Plasma B", "Macrophage", "CD16+ Monocyte", "cDC1", "Granulocyte-myeloid Sink", "Mast", "Erythrocyte", "Cyc Immune"
)

## Compute & plot correlation matrix between deconvoluted cell types (across GTEx samples)
# Compute correlation between gene expression profiles of cell types
# Z_state_gene_mean: rows are cell types, columns are genes (from above context)
celltype_expr_cormat <- cor(t(Z_state_gene_mean), method = "pearson", use = "pairwise.complete.obs")

# Reorder the correlation matrix accordingly
celltype_expr_cormat_ordered <- celltype_expr_cormat[celltype_order_by_compartment, rev(celltype_order_by_compartment)]

# Melt for ggplot
celltype_expr_cormat_long <- melt(celltype_expr_cormat_ordered)
colnames(celltype_expr_cormat_long) <- c("Celltype1", "Celltype2", "Correlation")

base_size <- 5
# Adjust steepness of color change for ease of visualization
red_steepness <- 0.3
exp_values <- seq(0, 1, length.out = 10)^red_steepness

plt1 <- ggplot(celltype_expr_cormat_long, aes(x = Celltype1, y = Celltype2, fill = Correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(
    # Simple color ramp from white to red
    colours = c("white", "red"),
    # The exponential sequence handles the "speed" of the color change
    values = exp_values,
    limits = c(0, 1),
    name = "Pearson\nCorrelation",
    guide = guide_colorbar(barwidth = 0.4, barheight = 6)
  ) +
  theme_minimal(base_family = "DejaVu Sans", base_size = base_size) +
  labs(
    title = "Correlation Matrix of Deconvoluted Cell Type Gene Expression Profiles",
    x = "",
    y = ""
  ) +
  theme(
    axis.text.y = element_text(size = base_size, color = "black", margin = margin(r = 4)),
    axis.text.x = element_text(size = base_size, color = "black", margin = margin(t = 4), angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_text(size = base_size, color = "black", margin = margin(t = 10)),
    axis.title.y = element_text(size = base_size, color = "black", margin = margin(r = 10)),
    plot.title = element_text(size = base_size + 2, hjust = 0.5, color = "black", margin = margin(b = 4)),
    legend.text = element_text(size = base_size, color = "black"),
    legend.ticks = element_blank()
  )

ggsave(filename = paste0(outs_dir, "/celltype_expr_cormat_ordered.pdf"), plot = plt1, width = 180,
  height = 120, units = "mm", device = cairo_pdf, dpi = 500
)




#### 7. Generate deconvolution summary figures ####

# Theta: rows are cell types, columns are GTEx samples
theta_df <- as.data.frame(celltype_theta)
theta_df$Celltype <- rownames(celltype_theta)
theta_long_df <- as.data.frame(as.table(as.matrix(celltype_theta)))
colnames(theta_long_df) <- c("Celltype", "sample_name", "prop")
theta_long_df <- merge(theta_long_df, sample_metadata, by.x = "sample_name", by.y = "SUBJID", all.x = TRUE, sort = FALSE)

# Stacked bar plot: x = Sample (GTEx donor), y = theta (proportion), fill = Celltype
# Group fine cell types into coarse cell types for ease of visualization with plotting
theta_long_df$Celltype_coarse <- dplyr::case_when(
  theta_long_df$Celltype %in% c("ArtEnd", "LympEnd", "Lymphatic Valve (ambiguous)", "CapEnd", "VenEnd") ~ "Endothelial",
  theta_long_df$Celltype %in% c("ArtSMC", "Stress-signaling Mural", "VenSMC", "Pericyte") ~ "Perivascular",
  theta_long_df$Celltype %in% c("Mesothelium", "PAX2+ Mullerian Epithelial", "PAX2- Mullerian Epithelial") ~ "Epithelial",
  theta_long_df$Celltype %in% c("Plasma B", "CD16+ NK", "CD8+ T", "CD4+ T", "B Cell", "CD16- NK", "Cyc Immune") ~ "Lymphoid Immune",
  theta_long_df$Celltype %in% c("CD16+ Monocyte", "Granulocyte-myeloid Sink", "Macrophage", "cDC1", "Mast", "Erythrocyte") ~ "Myeloid Immune",
  TRUE ~ theta_long_df$Celltype
)

# Ensure Celltype_coarse is a factor
theta_long_df$Celltype_coarse <- factor(theta_long_df$Celltype_coarse)

# Determine average (mean) proportion per coarse cell type to find the largest cell type on average
avg_props <- theta_long_df %>%
  group_by(Celltype_coarse) %>%
  summarize(mean_prop = mean(prop, na.rm = TRUE)) %>%
  arrange(desc(mean_prop))

# Sort celltypes so the largest (on average) are at the bottom of the stack
celltype_levels <- avg_props$Celltype_coarse[order(avg_props$mean_prop)]
theta_long_df$Celltype_coarse <- factor(theta_long_df$Celltype_coarse, levels = celltype_levels)

# Order samples by the  proportion of the largest cell type, descending
sample_order <- theta_long_df %>%
  filter(Celltype_coarse == avg_props$Celltype_coarse[1]) %>%
  group_by(sample_name) %>%
  summarize(largest_ct_prop = mean(prop, na.rm = TRUE)) %>%
  arrange(desc(largest_ct_prop)) %>%
  pull(sample_name)

# Set Sample factor levels according to the new order
theta_long_df$sample_name <- factor(theta_long_df$sample_name, levels = sample_order)

# Make a high-contrast palette for many classes
n_ct <- nlevels(theta_long_df$Celltype_coarse)
base_cols <- Polychrome::createPalette(n_ct, seedcolors = c("#000000", "#FFFFFF"))

# Reorder colors to maximize contrast between adjacent stack levels:
idx <- as.vector(rbind(seq_len(ceiling(n_ct/2)), seq(from = n_ct, to = ceiling(n_ct/2)+1, by = -1)))
idx <- idx[!is.na(idx)][1:n_ct]

ct_cols <- base_cols[idx]
names(ct_cols) <- levels(theta_long_df$Celltype_coarse)

theta_plot <- theta_long_df %>%
  group_by(sample_name, Celltype_coarse) %>%
  summarize(prop = sum(prop, na.rm = TRUE), .groups = "drop")

plt2 <- ggplot(theta_plot, aes(x = sample_name, y = prop, fill = Celltype_coarse)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8, color = NA) +
  scale_fill_manual(values = ct_cols) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  labs(
    title = "Ovary Cell Type Proportions per GTEx Donor",
    x = "GTEx Donor Cell Type Composition",
    y = "Cell Type Proportion",
    fill = "Cell Type Composition"
  ) +
  theme_minimal(base_family = "DejaVu Sans", base_size = base_size+2) +
  theme(
    axis.text.y = element_text(size = base_size+2, color = "black", margin = margin(r = 4)),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = base_size+2, color = "black", margin = margin(t = 6)),
    axis.title.y = element_text(size = base_size+2, color = "black", margin = margin(r = 6)),
    plot.title = element_text(size = base_size + 2, hjust = 0.5, color = "black", margin = margin(b = 4)),
    legend.text = element_text(size = base_size+2, color = "black"),
  )
plt2

ggsave(filename = paste0(outs_dir, "/celltype_proportions_stacked_barplot.pdf"), plot = plt2, width = 180,
  height = 100, units = "mm", device = cairo_pdf, dpi = 500
)


#### 8. Export deconvolution tables ####

# Save deconvolution results along with attached metadata for downstream use
write.csv(theta_long_df, file = paste0(outs_dir, "/deconv_results_with_metadata.csv"), row.names = FALSE)

# Save deconvolution results for publication and supplementary tables
GTEx_deconv_props_pub_supp <- theta_long_df %>%
  dplyr::transmute(
    gtex_sample_id = sample_name,
    cell_type = Celltype,
    proportion = prop
  )

write.csv(GTEx_deconv_props_pub_supp, file = paste0(outs_dir, "/deconv_props_pub_supp.csv"), row.names = FALSE)

# Calculate summary statistics for each deconvolved cell type & save to file
celltype_summary_stats <- theta_long_df %>%
  dplyr::group_by(Celltype) %>%
  dplyr::summarize(
    n_samples = dplyr::n_distinct(sample_name),
    mean_percent = round(mean(prop, na.rm = TRUE) * 100, 4),
    median_percent = round(median(prop, na.rm = TRUE) * 100, 4),
    sd_percent = round(sd(prop, na.rm = TRUE) * 100, 4),
    iqr_percent = round(IQR(prop, na.rm = TRUE) * 100, 4),
    min_percent = round(min(prop, na.rm = TRUE) * 100, 4),
    max_percent = round(max(prop, na.rm = TRUE) * 100, 4),
    .groups = "drop"
  ) %>%
  dplyr::rename(cell_type = Celltype) %>%
  dplyr::arrange(dplyr::desc(mean_percent))

write.csv(celltype_summary_stats, file = paste0(outs_dir, "/deconv_celltype_summary_stats.csv"), row.names = FALSE)

#### 9. Session Info ####
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
#  [1] UCell_2.12.0                Polychrome_1.5.4            reshape2_1.4.4              scales_1.4.0                stringr_1.5.1
#  [6] compositions_2.0-9          SingleCellExperiment_1.30.1 SummarizedExperiment_1.38.1 GenomicRanges_1.60.0        GenomeInfoDb_1.44.0
# [11] IRanges_2.42.0              S4Vectors_0.48.0            MatrixGenerics_1.20.0       matrixStats_1.5.0           zellkonverter_1.18.0
# [16] InstaPrism_0.1.6            pbapply_1.7-4               abind_1.4-8                 pheatmap_1.0.13             caret_7.0-1
# [21] lattice_0.22-7              ggplot2_3.5.2               Matrix_1.7-3                tidyr_1.3.1                 dplyr_1.1.4
# [26] BayesPrism_2.2.2            NMF_0.28                    Biobase_2.68.0              BiocGenerics_0.54.0         generics_0.1.4
# [31] cluster_2.1.8.1             rngtools_1.5.2              registry_0.5-1              snowfall_1.84-6.3           snow_0.4-4
# [36] naniar_1.1.0                readxl_1.4.5                readr_2.1.5

# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3      tensorA_0.36.2.1        jsonlite_2.0.0          magrittr_2.0.4          farver_2.1.2            vctrs_0.7.1
#   [7] S4Arrays_1.8.1          BiocNeighbors_2.0.1     cellranger_1.1.0        pROC_1.19.0.1           SparseArray_1.8.0       parallelly_1.46.1
#  [13] KernSmooth_2.23-26      basilisk_1.20.0         plyr_1.8.9              lubridate_1.9.4         igraph_2.1.4            lifecycle_1.0.5
#  [19] iterators_1.0.14        pkgconfig_2.0.3         rsvd_1.0.5              R6_2.6.1                GenomeInfoDbData_1.2.14 future_1.58.0
#  [25] digest_0.6.39           colorspace_2.1-2        dqrng_0.4.1             irlba_2.3.5.1           beachmat_2.24.0         filelock_1.0.3
#  [31] timechange_0.3.0        httr_1.4.7              compiler_4.4.1          withr_3.0.2             doParallel_1.0.17       BiocParallel_1.42.1
#  [37] gplots_3.3.0            bayesm_3.1-7            MASS_7.3-65             lava_1.8.2              DelayedArray_0.34.1     scatterplot3d_0.3-44
#  [43] bluster_1.20.0          gtools_3.9.5            caTools_1.18.3          ModelMetrics_1.2.2.2    tools_4.4.1             visdat_0.6.0
#  [49] future.apply_1.20.0     nnet_7.3-20             glue_1.8.0              nlme_3.1-168            grid_4.4.1              gridBase_0.4-7
#  [55] recipes_1.3.1           gtable_0.3.6            tzdb_0.5.0              class_7.3-23            data.table_1.18.2.1     hms_1.1.3
#  [61] BiocSingular_1.22.0     ScaledMatrix_1.14.0     metapod_1.18.0          XVector_0.48.0          foreach_1.5.2           pillar_1.11.0
#  [67] limma_3.64.1            robustbase_0.99-7       splines_4.4.1           survival_3.8-3          tidyselect_1.2.1        httpgd_2.0.4
#  [73] locfit_1.5-9.12         scuttle_1.16.0          edgeR_4.6.3             statmod_1.5.0           hardhat_1.4.2           timeDate_4052.112
#  [79] DEoptimR_1.1-4          stringi_1.8.7           UCSC.utils_1.4.0        unigd_0.1.3             codetools_0.2-20        tibble_3.3.0
#  [85] BiocManager_1.30.27     cli_3.6.5               rpart_4.1.24            reticulate_1.42.0       systemfonts_1.2.3       Rcpp_1.1.1
#  [91] dir.expiry_1.16.0       globals_0.19.0          png_0.1-8               parallel_4.4.1          gower_1.0.2             basilisk.utils_1.20.0
#  [97] scran_1.38.0            bitops_1.0-9            listenv_0.10.0          ipred_0.9-15            prodlim_2025.04.28      purrr_1.1.0
# [103] crayon_1.5.3            rlang_1.1.7