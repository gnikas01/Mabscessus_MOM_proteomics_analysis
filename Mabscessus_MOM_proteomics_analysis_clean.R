# M. abscessus proteomics analysis
#
# This script contains the R-based analysis steps used for the manuscript workflow:
#   1. FragPipe combined_protein import and MaxLFQ processing
#   2. Replicate-level protein filtering
#   3. QC outputs, including feature counts, PCA, scree plot, and Pearson correlation matrix
#   4. limma differential abundance analysis
#   5. Volcano plots for P8/P100, mutant CF/WT CF, and CF/TL contrasts
#   6. Fig. 2-style P8/P100 marker volcano plot highlighting Antigen 85, MspA-like,
#      ATP synthase, and NDH-1 marker proteins
#
# Presence/absence heatmaps, hybrid heatmap assembly, and ESX-dependency
# classifications were performed in Microsoft Excel using FragPipe detection outputs
# and limma contrast tables. Those Excel-based steps are intentionally not reproduced
# in this script.

# ==============================================================================
# 0. Package setup
# ==============================================================================

install_missing_packages <- FALSE

cran_packages <- c(
  "tidyverse",
  "readxl",
  "writexl",
  "ggrepel"
)

bioc_packages <- c("limma")

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!install_missing_packages) {
      stop(
        "Package '", pkg, "' is required but is not installed. ",
        "Install it manually or set install_missing_packages <- TRUE."
      )
    }
    install.packages(pkg)
  }
}

for (pkg in cran_packages) install_if_missing(pkg)

if (!requireNamespace("BiocManager", quietly = TRUE) && install_missing_packages) {
  install.packages("BiocManager")
}

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!install_missing_packages) {
      stop(
        "Bioconductor package '", pkg, "' is required but is not installed. ",
        "Install it manually or set install_missing_packages <- TRUE."
      )
    }
    BiocManager::install(pkg)
  }
}

library(tidyverse)
library(readxl)
library(writexl)
library(ggrepel)
library(limma)

# ==============================================================================
# 1. User settings
# ==============================================================================

# Recommended repository structure:
#   data/combined_protein.tsv or data/combined_protein.xlsx
#   outputs/
#
# The script accepts TSV/TXT/CSV or XLSX input. If using XLSX, set the sheet name
# below to the sheet containing the FragPipe combined_protein table.

combined_protein_file <- file.path("data", "combined_protein.tsv")
combined_protein_sheet <- "combined_protein"

# Optional plain-text file containing one UniProt ID per line for all curated pMOMPs.
# If present, these IDs are highlighted in the candidate P8/P100 volcano plot.
pmomp_id_file <- file.path("data", "pMOMP_ids.txt")

output_dir <- "outputs"
qc_dir <- file.path(output_dir, "QC")
limma_dir <- file.path(output_dir, "limma_outputs")
volcano_dir <- file.path(output_dir, "volcano_plots")

for (d in c(output_dir, qc_dir, limma_dir, volcano_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Manuscript filtering: retain proteins detected in at least 2 of 3 biological
# replicates in at least one strain-fraction condition.
min_valid_per_group <- 2

# FragPipe MaxLFQ values of 0 are treated as missing. Non-zero values are
# log2-transformed without imputation for limma.

# PCA/QC settings. PCA is run on complete cases when enough proteins remain;
# otherwise missing values are imputed only for PCA/QC visualization.
use_pca_imputation_if_needed <- TRUE
min_complete_proteins_for_pca <- 300

# Volcano plot thresholds used for visualization/candidate prioritization.
p_cutoff_marker <- 0.05
logfc_cutoff_membrane <- 0.32
padj_cutoff_membrane <- 0.10

padj_cutoff_secretome <- 0.05
logfc_cutoff_secretome <- 1.00

# ==============================================================================
# 2. Protein groups highlighted in volcano plots
# ==============================================================================

mspa_like_ids <- c(
  "B1MI08", "B1MK78", "B1MAY9", "B1MCA8", "B1MDB0",
  "B1MLI2", "B1MCA7", "B1MNS3", "B1MK79"
)

antigen85_ids <- c("B1MEL1", "B1MEL2", "B1MEL3")

# Marker proteins used in the Fig. 2-style P8/P100 volcano plot.
atp_synthase_ids <- c("B1MLV6", "B1MLV7", "B1MLV8")
ndh1_ids <- c("B1MAF5", "B1MAF6", "B1MAF7")

esx_substrate_ids <- unique(c(
  "B1MG82", "B1MG83", "B1MAN8", "B1MAP0", "B1MI82", "B1MAI5",
  "B1MAR4", "B1MEI4", "B1MNS7", "B1ME49", "B1MM66", "B1MNL2",
  "B1MHU5", "B1ME83", "B1MKT1", "B1MKS6", "B1MG58", "B1MFE3",
  "B1ME55", "B1MG75", "B1MAP1", "B1MHU1", "B1ME54", "B1MHU0",
  "B1MAN9", "B1MG76", "B1MG59", "B1MAY6", "B1MD68", "B1MK25",
  "B1MD69", "B1ME85", "B1MKT2", "B1MFE2", "B1MFS5", "B1MNL3",
  "B1MAY7", "B1MJN1", "B1MK26", "B1MMB3", "B1MMB4", "B1MKS9",
  "B1ME84", "B1MEI5", "B1ME46", "B1MMB2", "B1MD70", "B1MK27",
  "B1MFZ3", "B1MIJ3", "B1MJ72", "B1MBY8", "B1MJB0", "B1MNX7",
  "B1MJA3", "B1MFM0", "B1MKA2", "B1ML06"
))

read_optional_id_file <- function(path) {
  if (!file.exists(path)) return(character())
  ids <- readLines(path, warn = FALSE)
  ids <- stringr::str_trim(ids)
  ids <- ids[ids != "" & !stringr::str_starts(ids, "#")]
  unique(ids)
}

pmomp_ids <- read_optional_id_file(pmomp_id_file)

# ==============================================================================
# 3. Helper functions
# ==============================================================================

find_column <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

read_combined_protein <- function(path, sheet = NULL) {
  if (!file.exists(path)) {
    stop(
      "Could not find combined_protein input file: ", path, "\n",
      "Place the file in data/ or update combined_protein_file at the top of this script."
    )
  }

  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    readxl::read_excel(path, sheet = sheet)
  } else if (ext %in% c("tsv", "txt")) {
    readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
  } else if (ext == "csv") {
    readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  } else {
    stop("Unsupported input file extension: ", ext)
  }
}

parse_sample_metadata <- function(maxlfq_cols) {
  tibble(sample_col = maxlfq_cols) %>%
    mutate(
      sample = str_remove(sample_col, "\\s*MaxLFQ Intensity$"),
      sample_type = case_when(
        str_detect(sample, "^BlankUltra_") ~ "BlankUltra",
        str_detect(sample, "^Blank_") ~ "Blank",
        TRUE ~ "Biological"
      ),
      genotype = case_when(
        str_detect(sample, "^WT_") ~ "WT",
        str_detect(sample, "^C3C4_") ~ "C3C4",
        str_detect(sample, "^C3_") ~ "C3",
        str_detect(sample, "^C4_") ~ "C4",
        TRUE ~ NA_character_
      ),
      fraction = case_when(
        str_detect(sample, "_CF_") ~ "CF",
        str_detect(sample, "_P100_") ~ "P100",
        str_detect(sample, "_P8_") ~ "P8",
        str_detect(sample, "_TL_") ~ "TL",
        TRUE ~ NA_character_
      ),
      replicate = str_extract(sample, "(?<=_)\\d+$"),
      group = if_else(sample_type == "Biological", paste(genotype, fraction, sep = "_"), sample)
    ) %>%
    mutate(
      genotype = factor(genotype, levels = c("WT", "C3", "C4", "C3C4")),
      fraction = factor(fraction, levels = c("TL", "CF", "P8", "P100"))
    )
}

impute_low_per_column <- function(mat) {
  mat_imp <- mat
  for (j in seq_len(ncol(mat_imp))) {
    observed <- mat_imp[, j][!is.na(mat_imp[, j])]
    if (length(observed) > 5) {
      imp_val <- as.numeric(stats::quantile(observed, probs = 0.01, na.rm = TRUE)) - 1
      mat_imp[is.na(mat_imp[, j]), j] <- imp_val
    }
  }
  mat_imp
}

safe_name <- function(x) {
  str_replace_all(x, "[^A-Za-z0-9_]+", "_")
}

contrast_if_groups_exist <- function(contrast_name, contrast_formula, design_colnames) {
  groups <- str_split(contrast_formula, " - ", simplify = TRUE)
  if (all(groups %in% design_colnames)) {
    tibble(contrast = contrast_name, formula = contrast_formula)
  } else {
    tibble(contrast = character(), formula = character())
  }
}

assign_highlight_group <- function(df, highlight_lists, highlight_priority) {
  df$highlight_group <- "Other"
  for (grp in rev(highlight_priority)) {
    ids <- highlight_lists[[grp]]
    if (!is.null(ids) && length(ids) > 0) {
      df$highlight_group[df$protein_id %in% ids] <- grp
    }
  }
  df
}

make_volcano_plot <- function(limma_df,
                              contrasts_to_plot,
                              facet_labels,
                              highlight_lists,
                              highlight_priority,
                              output_prefix,
                              plot_title,
                              p_cutoff,
                              logfc_cutoff,
                              output_subdir,
                              p_value_column = c("adj.P.Val", "P.Value"),
                              plot_width = 10,
                              plot_height = 8,
                              y_limit = NULL,
                              label_highlighted_if_significant = TRUE,
                              x_axis_label = "log2 fold change") {

  p_value_column <- match.arg(p_value_column)
  dir.create(output_subdir, recursive = TRUE, showWarnings = FALSE)

  plot_df <- limma_df %>%
    filter(contrast %in% contrasts_to_plot) %>%
    mutate(
      contrast = factor(contrast, levels = contrasts_to_plot),
      p_for_plot = .data[[p_value_column]],
      p_for_plot = if_else(is.na(p_for_plot), 1, p_for_plot),
      p_for_plot = pmax(p_for_plot, 1e-300),
      neg_log10_p = -log10(p_for_plot),
      significance = case_when(
        !is.na(.data[[p_value_column]]) & .data[[p_value_column]] < p_cutoff & logFC >= logfc_cutoff ~ "Significant_up",
        !is.na(.data[[p_value_column]]) & .data[[p_value_column]] < p_cutoff & logFC <= -logfc_cutoff ~ "Significant_down",
        TRUE ~ "Not_significant"
      )
    )

  plot_df <- assign_highlight_group(plot_df, highlight_lists, highlight_priority)
  highlight_groups_present <- names(highlight_lists)[vapply(highlight_lists, length, integer(1)) > 0]

  plot_df <- plot_df %>%
    mutate(
      plot_group = case_when(
        highlight_group %in% highlight_groups_present ~ highlight_group,
        significance == "Significant_up" ~ "Significant_up",
        significance == "Significant_down" ~ "Significant_down",
        TRUE ~ "Background"
      ),
      label = if_else(
        label_highlighted_if_significant & highlight_group %in% highlight_groups_present & significance != "Not_significant",
        protein_id,
        NA_character_
      )
    )

  x_max <- ceiling(max(abs(plot_df$logFC), na.rm = TRUE) + 0.25)

  if (is.null(y_limit)) {
    finite_y <- is.finite(plot_df$neg_log10_p)
    max_y <- if (any(finite_y)) max(plot_df$neg_log10_p[finite_y], na.rm = TRUE) * 1.05 else 1
  } else {
    max_y <- y_limit
  }

  point_colors <- c(
    pMOMPs = "#2C7FB8",
    Antigen85 = "#8C6D31",
    ESX_substrates = "#F28E2B",
    MspA_like = "#C0392B",
    ATP_synthase = "#FFD92F",
    NDH_1 = "#6A3D9A"
  )
  point_labels <- c(
    pMOMPs = "pMOMPs",
    Antigen85 = "Antigen 85 complex",
    ESX_substrates = "ESX substrates",
    MspA_like = "MspA homologs",
    ATP_synthase = "ATP synthase",
    NDH_1 = "NDH-1 (NuoL/M/N)"
  )
  point_colors <- point_colors[names(point_colors) %in% highlight_groups_present]
  point_labels <- point_labels[names(point_labels) %in% highlight_groups_present]

  p <- ggplot(plot_df, aes(x = logFC, y = neg_log10_p)) +
    geom_point(
      data = plot_df %>% filter(!(plot_group %in% highlight_groups_present)),
      color = "grey78", size = 2, alpha = 0.35, show.legend = FALSE
    ) +
    geom_point(
      data = plot_df %>% filter(plot_group %in% highlight_groups_present),
      aes(color = plot_group), size = 2, alpha = 0.90, show.legend = TRUE
    ) +
    geom_text_repel(
      data = plot_df %>% filter(!is.na(label)),
      aes(label = label),
      size = 3,
      fontface = "bold",
      max.overlaps = 100,
      box.padding = 0.6,
      point.padding = 0.2,
      segment.color = "grey40",
      segment.size = 0.3,
      min.segment.length = 0,
      show.legend = FALSE
    ) +
    geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dashed", linewidth = 0.4) +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", linewidth = 0.4) +
    facet_wrap(~ contrast, labeller = labeller(contrast = facet_labels)) +
    coord_cartesian(xlim = c(-x_max, x_max), ylim = c(0, max_y)) +
    labs(
      title = plot_title,
      subtitle = paste0("Cutoffs: ", p_value_column, " < ", p_cutoff, ", |log2FC| >= ", logfc_cutoff),
      x = x_axis_label,
      y = if_else(p_value_column == "adj.P.Val", "-log10 adjusted P value", "-log10 P value"),
      color = NULL
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 18),
      plot.subtitle = element_text(size = 12, margin = margin(b = 10)),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
      strip.text = element_text(size = 14, face = "bold"),
      legend.position = "bottom"
    )

  if (length(point_colors) > 0) {
    p <- p + scale_color_manual(values = point_colors, labels = point_labels, drop = FALSE)
  }

  print(p)

  write_csv(plot_df, file.path(output_subdir, paste0(output_prefix, "_plot_data.csv")))
  write_xlsx(plot_df, file.path(output_subdir, paste0(output_prefix, "_plot_data.xlsx")))
  ggsave(file.path(output_subdir, paste0(output_prefix, ".png")), p, width = plot_width, height = plot_height, dpi = 300)
  ggsave(file.path(output_subdir, paste0(output_prefix, ".pdf")), p, width = plot_width, height = plot_height)

  invisible(p)
}

# ==============================================================================
# 4. Import FragPipe combined_protein and prepare MaxLFQ matrix
# ==============================================================================

message("Reading combined_protein file: ", combined_protein_file)
protein_df <- read_combined_protein(combined_protein_file, sheet = combined_protein_sheet)

protein_id_col <- find_column(protein_df, c("Protein ID", "Protein", "UniProt ID", "Protein.ID"))
gene_col <- find_column(protein_df, c("Gene", "Gene Name", "gene_name"))
desc_col <- find_column(protein_df, c("Description", "Protein Description", "description"))

if (is.na(protein_id_col)) stop("Could not find a protein identifier column.")

protein_df <- protein_df %>%
  filter(!is.na(.data[[protein_id_col]]), .data[[protein_id_col]] != "") %>%
  distinct(.data[[protein_id_col]], .keep_all = TRUE)

maxlfq_cols <- grep("\\s*MaxLFQ Intensity$", colnames(protein_df), value = TRUE)
if (length(maxlfq_cols) == 0) stop("No MaxLFQ Intensity columns found.")

sample_info <- parse_sample_metadata(maxlfq_cols)
bio_info <- sample_info %>%
  filter(sample_type == "Biological", !is.na(genotype), !is.na(fraction)) %>%
  mutate(group = factor(group, levels = unique(group)))

bio_cols <- bio_info$sample_col
bio_samples <- bio_info$sample

feature_annot <- tibble(
  protein_id = as.character(protein_df[[protein_id_col]]),
  gene = if (!is.na(gene_col)) as.character(protein_df[[gene_col]]) else NA_character_,
  description = if (!is.na(desc_col)) as.character(protein_df[[desc_col]]) else NA_character_
) %>%
  mutate(gene = if_else(is.na(gene) | gene == "", protein_id, gene))

expr_df <- protein_df[, bio_cols, drop = FALSE]
expr_mat <- as.matrix(as.data.frame(lapply(expr_df, function(x) suppressWarnings(as.numeric(as.character(x)))), check.names = FALSE))
rownames(expr_mat) <- feature_annot$protein_id
colnames(expr_mat) <- bio_samples

expr_mat[expr_mat == 0] <- NA

write_csv(sample_info, file.path(qc_dir, "sample_metadata_all_columns.csv"))
write_csv(bio_info, file.path(qc_dir, "sample_metadata_biological.csv"))

# ==============================================================================
# 5. Protein filtering and log2 transformation
# ==============================================================================

group_factor <- bio_info$group
names(group_factor) <- bio_info$sample

valid_by_group <- sapply(levels(group_factor), function(grp) {
  cols <- bio_info$sample[bio_info$group == grp]
  rowSums(!is.na(expr_mat[, cols, drop = FALSE]))
})

if (is.null(dim(valid_by_group))) {
  valid_by_group <- matrix(valid_by_group, ncol = 1)
  colnames(valid_by_group) <- levels(group_factor)
}

keep_protein <- apply(valid_by_group, 1, function(x) any(x >= min_valid_per_group))
expr_keep_raw <- expr_mat[keep_protein, , drop = FALSE]
feature_annot_keep <- feature_annot[keep_protein, , drop = FALSE]
expr_keep_log2 <- log2(expr_keep_raw)

message("Proteins before filtering: ", nrow(expr_mat))
message("Proteins after filtering: ", nrow(expr_keep_log2))

filter_summary <- feature_annot %>%
  mutate(
    total_biological_detections = rowSums(!is.na(expr_mat)),
    retained_after_filtering = keep_protein
  ) %>%
  bind_cols(as_tibble(valid_by_group, .name_repair = "unique"))

write_csv(filter_summary, file.path(qc_dir, "protein_detection_filter_summary.csv"))
write_xlsx(filter_summary, file.path(qc_dir, "protein_detection_filter_summary.xlsx"))

filtered_matrix <- feature_annot_keep %>%
  bind_cols(as_tibble(expr_keep_log2, .name_repair = "unique"))
write_csv(filtered_matrix, file.path(qc_dir, "filtered_log2_maxlfq_matrix.csv"))
write_xlsx(filtered_matrix, file.path(qc_dir, "filtered_log2_maxlfq_matrix.xlsx"))

# ==============================================================================
# 6. QC: feature counts, PCA, scree plot, and Pearson correlation matrix
# ==============================================================================

feature_counts <- tibble(
  sample = colnames(expr_mat),
  n_quantified_proteins = colSums(!is.na(expr_mat))
) %>%
  left_join(bio_info, by = "sample") %>%
  arrange(fraction, genotype, replicate)

write_csv(feature_counts, file.path(qc_dir, "features_per_sample_counts.csv"))
write_xlsx(feature_counts, file.path(qc_dir, "features_per_sample_counts.xlsx"))

fraction_colors <- c(
  TL = "#4D4D4D",
  CF = "#377EB8",
  P8 = "#E41A1C",
  P100 = "#4DAF4A"
)

strain_shapes <- c(
  WT = 16,
  C3 = 17,
  C4 = 15,
  C3C4 = 18
)

p_features <- ggplot(feature_counts, aes(x = fraction, y = n_quantified_proteins, color = fraction, shape = genotype)) +
  geom_jitter(width = 0.15, height = 0, size = 3.5) +
  scale_color_manual(values = fraction_colors, drop = FALSE) +
  scale_shape_manual(values = strain_shapes, drop = FALSE) +
  labs(
    title = "Quantified proteins per sample",
    x = "Fraction",
    y = "Number of quantified proteins",
    color = "Fraction",
    shape = "Strain"
  ) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(file.path(qc_dir, "features_per_sample_grouped.png"), p_features, width = 7, height = 5, dpi = 300)
ggsave(file.path(qc_dir, "features_per_sample_grouped.pdf"), p_features, width = 7, height = 5)

log_expr_complete <- expr_keep_log2[complete.cases(expr_keep_log2), , drop = FALSE]

if (nrow(log_expr_complete) >= min_complete_proteins_for_pca) {
  pca_input <- log_expr_complete
  pca_mode <- "complete_cases"
} else if (use_pca_imputation_if_needed) {
  pca_input <- impute_low_per_column(expr_keep_log2)
  pca_mode <- "low_intensity_imputed_for_PCA_only"
} else {
  pca_input <- log_expr_complete
  pca_mode <- "complete_cases_low_N"
}

message("PCA mode: ", pca_mode)
message("Proteins used for PCA: ", nrow(pca_input))

pca_res <- prcomp(t(pca_input), center = TRUE, scale. = TRUE)
percent_var <- (pca_res$sdev^2 / sum(pca_res$sdev^2)) * 100

pca_df <- as.data.frame(pca_res$x[, seq_len(min(4, ncol(pca_res$x))), drop = FALSE]) %>%
  rownames_to_column("sample") %>%
  left_join(bio_info, by = "sample")

write_csv(pca_df, file.path(qc_dir, "PCA_coordinates.csv"))
write_xlsx(pca_df, file.path(qc_dir, "PCA_coordinates.xlsx"))

p_pca_12 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = fraction, shape = genotype)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = sample), size = 3, max.overlaps = 50) +
  scale_color_manual(values = fraction_colors, drop = FALSE) +
  scale_shape_manual(values = strain_shapes, drop = FALSE) +
  labs(
    title = paste0("PCA of MaxLFQ intensities (", pca_mode, ")"),
    x = paste0("PC1: ", round(percent_var[1], 1), "%"),
    y = paste0("PC2: ", round(percent_var[2], 1), "%"),
    color = "Fraction",
    shape = "Strain"
  ) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(file.path(qc_dir, "PCA_PC1_PC2.png"), p_pca_12, width = 10, height = 7, dpi = 300)
ggsave(file.path(qc_dir, "PCA_PC1_PC2.pdf"), p_pca_12, width = 10, height = 7)

if ("PC3" %in% colnames(pca_df)) {
  p_pca_13 <- ggplot(pca_df, aes(x = PC1, y = PC3, color = fraction, shape = genotype)) +
    geom_point(size = 4) +
    geom_text_repel(aes(label = sample), size = 3, max.overlaps = 50) +
    scale_color_manual(values = fraction_colors, drop = FALSE) +
    scale_shape_manual(values = strain_shapes, drop = FALSE) +
    labs(
      title = paste0("PCA of MaxLFQ intensities (", pca_mode, ")"),
      x = paste0("PC1: ", round(percent_var[1], 1), "%"),
      y = paste0("PC3: ", round(percent_var[3], 1), "%"),
      color = "Fraction",
      shape = "Strain"
    ) +
    theme_classic(base_size = 14) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))

  ggsave(file.path(qc_dir, "PCA_PC1_PC3.png"), p_pca_13, width = 10, height = 7, dpi = 300)
  ggsave(file.path(qc_dir, "PCA_PC1_PC3.pdf"), p_pca_13, width = 10, height = 7)
}

scree_df <- tibble(
  PC = paste0("PC", seq_along(percent_var)),
  variance_explained = percent_var
)
write_csv(scree_df, file.path(qc_dir, "PCA_scree_values.csv"))

p_scree <- ggplot(slice_head(scree_df, n = 10), aes(x = reorder(PC, variance_explained), y = variance_explained)) +
  geom_col(fill = "grey40") +
  coord_flip() +
  labs(title = "PCA scree plot", x = NULL, y = "% variance explained") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(file.path(qc_dir, "PCA_scree_plot.png"), p_scree, width = 6, height = 5, dpi = 300)
ggsave(file.path(qc_dir, "PCA_scree_plot.pdf"), p_scree, width = 6, height = 5)

cor_mat <- cor(pca_input, use = "pairwise.complete.obs", method = "pearson")
write_csv(as.data.frame(cor_mat) %>% rownames_to_column("sample"), file.path(qc_dir, "pearson_sample_correlation_matrix.csv"))
write_xlsx(as.data.frame(cor_mat) %>% rownames_to_column("sample"), file.path(qc_dir, "pearson_sample_correlation_matrix.xlsx"))

# ==============================================================================
# 7. limma differential abundance analysis
# ==============================================================================

design <- model.matrix(~ 0 + group, data = bio_info)
colnames(design) <- levels(bio_info$group)

fit <- lmFit(expr_keep_log2, design)

# Explicit contrasts used for manuscript figures and downstream Excel assembly.
contrast_rows <- list()

for (g in c("WT", "C3", "C4", "C3C4")) {
  contrast_rows[[paste0(g, "_P8_vs_", g, "_P100")]] <- contrast_if_groups_exist(
    paste0(g, "_P8_vs_", g, "_P100"),
    paste0(g, "_P8 - ", g, "_P100"),
    colnames(design)
  )
  contrast_rows[[paste0(g, "_CF_vs_", g, "_TL")]] <- contrast_if_groups_exist(
    paste0(g, "_CF_vs_", g, "_TL"),
    paste0(g, "_CF - ", g, "_TL"),
    colnames(design)
  )
}

for (fr in c("TL", "CF", "P8", "P100")) {
  for (g in c("C3", "C4", "C3C4")) {
    contrast_rows[[paste0(g, "_", fr, "_vs_WT_", fr)]] <- contrast_if_groups_exist(
      paste0(g, "_", fr, "_vs_WT_", fr),
      paste0(g, "_", fr, " - WT_", fr),
      colnames(design)
    )
  }
}

contrast_catalog <- bind_rows(contrast_rows)
if (nrow(contrast_catalog) == 0) stop("No valid contrasts were generated. Check sample names and metadata parsing.")

contrast_matrix <- makeContrasts(contrasts = contrast_catalog$formula, levels = design)
colnames(contrast_matrix) <- contrast_catalog$contrast

write_csv(contrast_catalog, file.path(limma_dir, "contrast_catalog.csv"))
write_xlsx(contrast_catalog, file.path(limma_dir, "contrast_catalog.xlsx"))
write.table(contrast_matrix, file = file.path(limma_dir, "contrast_matrix.tsv"), sep = "\t", quote = FALSE, col.names = NA)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

export_one_contrast <- function(contrast_name) {
  contrast_formula <- contrast_catalog$formula[match(contrast_name, contrast_catalog$contrast)]
  tt <- topTable(fit2, coef = contrast_name, number = Inf, sort.by = "P") %>%
    rownames_to_column("protein_id") %>%
    left_join(feature_annot_keep, by = "protein_id") %>%
    relocate(protein_id, gene, description) %>%
    mutate(
      contrast = contrast_name,
      contrast_formula = contrast_formula,
      significant_adjP_0.05 = adj.P.Val < 0.05,
      significant_adjP_0.10 = adj.P.Val < 0.10
    )

  groups_in_contrast <- str_split(contrast_formula, " - ", simplify = TRUE)
  group_a <- groups_in_contrast[1]
  group_b <- groups_in_contrast[2]

  samples_a <- bio_info %>% filter(group == group_a) %>% pull(sample)
  samples_b <- bio_info %>% filter(group == group_b) %>% pull(sample)

  mean_a <- rowMeans(expr_keep_log2[, samples_a, drop = FALSE], na.rm = TRUE)
  mean_b <- rowMeans(expr_keep_log2[, samples_b, drop = FALSE], na.rm = TRUE)

  tt[[paste0("mean_log2_", group_a)]] <- mean_a[match(tt$protein_id, names(mean_a))]
  tt[[paste0("mean_log2_", group_b)]] <- mean_b[match(tt$protein_id, names(mean_b))]

  write_csv(tt, file.path(limma_dir, paste0(contrast_name, "_limma.csv")))
  write_xlsx(tt, file.path(limma_dir, paste0(contrast_name, "_limma.xlsx")))

  tt
}

limma_results <- map_dfr(contrast_catalog$contrast, export_one_contrast)
write_csv(limma_results, file.path(limma_dir, "all_limma_contrasts_long.csv"))
write_xlsx(limma_results, file.path(limma_dir, "all_limma_contrasts_long.xlsx"))

limma_summary <- limma_results %>%
  group_by(contrast, contrast_formula) %>%
  summarise(
    n_tested = n(),
    n_adjP_lt_0.05 = sum(adj.P.Val < 0.05, na.rm = TRUE),
    n_adjP_lt_0.10 = sum(adj.P.Val < 0.10, na.rm = TRUE),
    top_hit = protein_id[which.min(P.Value)],
    top_hit_logFC = logFC[which.min(P.Value)],
    top_hit_adjP = adj.P.Val[which.min(P.Value)],
    .groups = "drop"
  )

write_csv(limma_summary, file.path(limma_dir, "limma_contrast_summary.csv"))
write_xlsx(limma_summary, file.path(limma_dir, "limma_contrast_summary.xlsx"))

# ==============================================================================
# 8. Volcano plots
# ==============================================================================

membrane_contrasts <- c(
  "WT_P8_vs_WT_P100",
  "C3_P8_vs_C3_P100",
  "C4_P8_vs_C4_P100",
  "C3C4_P8_vs_C3C4_P100"
)
membrane_contrasts <- membrane_contrasts[membrane_contrasts %in% unique(limma_results$contrast)]

membrane_labels <- c(
  WT_P8_vs_WT_P100 = "WT",
  C3_P8_vs_C3_P100 = "ΔeccC3",
  C4_P8_vs_C4_P100 = "ΔeccC4",
  C3C4_P8_vs_C3C4_P100 = "ΔeccC3/C4"
)

# 8A. Fig. 2-style marker volcano plot using nominal P.Value < 0.05 and |log2FC| > 0.32.
# This matches the marker volcano logic used for the manuscript Fig. 2 P8/P100 panel.
marker_highlights <- list(
  Antigen85 = antigen85_ids,
  MspA_like = mspa_like_ids,
  ATP_synthase = atp_synthase_ids,
  NDH_1 = ndh1_ids
)

if (length(membrane_contrasts) > 0) {
  make_volcano_plot(
    limma_df = limma_results,
    contrasts_to_plot = membrane_contrasts,
    facet_labels = membrane_labels,
    highlight_lists = marker_highlights,
    highlight_priority = c("MspA_like", "ATP_synthase", "Antigen85", "NDH_1"),
    output_prefix = "Fig2_P8_vs_P100_marker_volcano",
    plot_title = "P8 vs P100 marker volcano plots",
    p_cutoff = p_cutoff_marker,
    logfc_cutoff = logfc_cutoff_membrane,
    output_subdir = file.path(volcano_dir, "Fig2_P8_vs_P100_markers"),
    p_value_column = "P.Value",
    plot_width = 10.5,
    plot_height = 9,
    label_highlighted_if_significant = FALSE,
    x_axis_label = "log2 fold change (P8/P100)"
  )
}

# 8B. Candidate-focused P8/P100 volcano plot. Curated pMOMPs are highlighted if
# data/pMOMP_ids.txt is provided; otherwise this plot still highlights MspA-like,
# ESX substrate, and Antigen 85 groups.
membrane_candidate_highlights <- list(
  pMOMPs = pmomp_ids,
  Antigen85 = antigen85_ids,
  ESX_substrates = esx_substrate_ids,
  MspA_like = mspa_like_ids
)

if (length(membrane_contrasts) > 0) {
  make_volcano_plot(
    limma_df = limma_results,
    contrasts_to_plot = membrane_contrasts,
    facet_labels = membrane_labels,
    highlight_lists = membrane_candidate_highlights,
    highlight_priority = c("MspA_like", "ESX_substrates", "Antigen85", "pMOMPs"),
    output_prefix = "P8_vs_P100_candidate_volcano",
    plot_title = "P8 vs P100 candidate volcano plots",
    p_cutoff = padj_cutoff_membrane,
    logfc_cutoff = logfc_cutoff_membrane,
    output_subdir = file.path(volcano_dir, "P8_vs_P100_candidates"),
    p_value_column = "adj.P.Val",
    plot_width = 10.5,
    plot_height = 9,
    x_axis_label = "log2 fold change (P8/P100)"
  )
}

# 8C. Mutant CF vs WT CF volcano plots: ESX-associated culture filtrate recovery.
mutant_cf_contrasts <- c(
  "C3_CF_vs_WT_CF",
  "C4_CF_vs_WT_CF",
  "C3C4_CF_vs_WT_CF"
)
mutant_cf_contrasts <- mutant_cf_contrasts[mutant_cf_contrasts %in% unique(limma_results$contrast)]

mutant_cf_labels <- c(
  C3_CF_vs_WT_CF = "ΔeccC3 vs WT",
  C4_CF_vs_WT_CF = "ΔeccC4 vs WT",
  C3C4_CF_vs_WT_CF = "ΔeccC3/C4 vs WT"
)

esx_highlights <- list(ESX_substrates = esx_substrate_ids)

if (length(mutant_cf_contrasts) > 0) {
  make_volcano_plot(
    limma_df = limma_results,
    contrasts_to_plot = mutant_cf_contrasts,
    facet_labels = mutant_cf_labels,
    highlight_lists = esx_highlights,
    highlight_priority = c("ESX_substrates"),
    output_prefix = "Mutant_CF_vs_WT_CF_3panel_volcano",
    plot_title = "Mutant CF vs WT CF volcano plots",
    p_cutoff = padj_cutoff_secretome,
    logfc_cutoff = logfc_cutoff_secretome,
    output_subdir = file.path(volcano_dir, "Mutant_CF_vs_WT_CF"),
    p_value_column = "adj.P.Val",
    plot_width = 10,
    plot_height = 7.5,
    x_axis_label = "log2 fold change"
  )
}

# 8D. CF vs TL volcano plots: culture filtrate enrichment relative to total lysate.
cf_vs_tl_contrasts <- c(
  "WT_CF_vs_WT_TL",
  "C3_CF_vs_C3_TL",
  "C4_CF_vs_C4_TL",
  "C3C4_CF_vs_C3C4_TL"
)
cf_vs_tl_contrasts <- cf_vs_tl_contrasts[cf_vs_tl_contrasts %in% unique(limma_results$contrast)]

cf_vs_tl_labels <- c(
  WT_CF_vs_WT_TL = "WT",
  C3_CF_vs_C3_TL = "ΔeccC3",
  C4_CF_vs_C4_TL = "ΔeccC4",
  C3C4_CF_vs_C3C4_TL = "ΔeccC3/C4"
)

if (length(cf_vs_tl_contrasts) > 0) {
  make_volcano_plot(
    limma_df = limma_results,
    contrasts_to_plot = cf_vs_tl_contrasts,
    facet_labels = cf_vs_tl_labels,
    highlight_lists = esx_highlights,
    highlight_priority = c("ESX_substrates"),
    output_prefix = "CF_vs_TL_4panel_volcano",
    plot_title = "CF vs TL volcano plots",
    p_cutoff = padj_cutoff_secretome,
    logfc_cutoff = logfc_cutoff_secretome,
    output_subdir = file.path(volcano_dir, "CF_vs_TL"),
    p_value_column = "adj.P.Val",
    plot_width = 10.5,
    plot_height = 9,
    x_axis_label = "log2 fold change"
  )
}

# ==============================================================================
# 9. Session information
# ==============================================================================

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("Analysis complete. Outputs written to: ", normalizePath(output_dir))
