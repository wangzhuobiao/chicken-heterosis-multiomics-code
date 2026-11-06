#!/usr/bin/env Rscript

# ==============================================================================
# 05C_duodenum_F1_inheritance_pattern_complete_v4_all_samples.R
# Standalone F1 expression-inheritance analysis plus categorized-gene enrichment for duodenum RNA-seq
#
# Groups: YY, YX (F1), XX
#
# Required upstream files in RNA_CORE_OUTDIR:
#   Table3_VST_expression_matrix.csv
#   DEG_YY_vs_XX.csv
#   DEG_YX_vs_YY.csv
#   DEG_YX_vs_XX.csv
#
# Optional upstream file:
#   Table2_gene_annotation_used.csv
#
# Important:
#   1. This script does NOT rerun DESeq2.
#   3. Every sample present in both the upstream VST matrix and metadata is used.
#      A sample absent from the upstream VST matrix cannot be restored here and
#      will only be reported in the sample-QC table.
#   4. Percentages in the main pattern plot use the total number of classifiable
#      genes as denominator, independently for every pattern.
# ==============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

# ==============================================================================
# 0. User configuration through environment variables
# ==============================================================================

CORE_OUTDIR <- Sys.getenv(
  "RNA_CORE_OUTDIR",
  unset = ""
)

META_FILE <- Sys.getenv(
  "META_FILE",
  unset = ""
)

OUTDIR <- Sys.getenv(
  "RNA_F1_OUTDIR",
  unset = file.path(CORE_OUTDIR, "F1_inheritance_standalone2")
)

TISSUE_LABEL <- Sys.getenv("RNA_TISSUE_LABEL", unset = "duodenum")
TISSUE_TITLE <- Sys.getenv("RNA_TISSUE_TITLE", unset = "Duodenum transcriptome")
GENE_SET_CONTEXT <- tolower(Sys.getenv("RNA_F1_GENE_SET_CONTEXT", unset = "duodenum"))

# Significance column used for pattern classification: padj (recommended) or pvalue.
SIG_COLUMN <- tolower(Sys.getenv("RNA_PATTERN_SIGNIFICANCE_COLUMN", unset = "pvalue"))
PATTERN_ALPHA <- as.numeric(Sys.getenv("RNA_PATTERN_ALPHA", unset = "0.05"))
PAIRWISE_LFC_CUTOFF <- as.numeric(Sys.getenv("RNA_PATTERN_LFC_CUTOFF", unset = "0"))

# Numerical tolerance for checking whether YX lies within/outside the parental range.
POSITION_TOLERANCE <- as.numeric(Sys.getenv("RNA_PATTERN_POSITION_TOLERANCE", unset = "0.02"))

# Gene-label settings for the mid-parent-deviation plot.
LABEL_TOP_N_DEG <- as.integer(Sys.getenv("RNA_F1_LABEL_TOP_N_DEG", unset = "15"))
LABEL_TOP_N_DEVIATION <- as.integer(Sys.getenv("RNA_F1_LABEL_TOP_N_DEVIATION", unset = "10"))
LABEL_CANDIDATES <- toupper(Sys.getenv("RNA_F1_LABEL_CANDIDATES", unset = "TRUE")) == "TRUE"
MAX_TOTAL_LABELS <- as.integer(Sys.getenv("RNA_F1_MAX_TOTAL_LABELS", unset = "35"))

# Heatmap settings.
TOP_N_PER_PATTERN <- as.integer(Sys.getenv("RNA_F1_TOP_N_PER_PATTERN", unset = "10"))
HEATMAP_MAX_GENES <- as.integer(Sys.getenv("RNA_F1_HEATMAP_MAX_GENES", unset = "50"))

# Enrichment settings for genes after F1-pattern classification.
# This section is intentionally downstream of classification and does not change
# any pattern assignment. It only asks whether genes in each class are enriched
# for GO/KEGG terms or for curated duodenal lipid-handling modules.
RUN_PATTERN_ENRICHMENT <- toupper(Sys.getenv("RNA_F1_RUN_PATTERN_ENRICHMENT", unset = "TRUE")) == "TRUE"
ENRICH_MIN_SET_SIZE <- as.integer(Sys.getenv("RNA_F1_ENRICH_MIN_SET_SIZE", unset = "5"))
ENRICH_SHOW_CATEGORY <- as.integer(Sys.getenv("RNA_F1_ENRICH_SHOW_CATEGORY", unset = "12"))
ENRICH_PVALUE_CUTOFF <- as.numeric(Sys.getenv("RNA_F1_ENRICH_PVALUE_CUTOFF", unset = "0.05"))
ENRICH_QVALUE_CUTOFF <- as.numeric(Sys.getenv("RNA_F1_ENRICH_QVALUE_CUTOFF", unset = "0.20"))
ENRICH_REPORT_FDR <- as.numeric(Sys.getenv("RNA_F1_ENRICH_REPORT_FDR", unset = "0.25"))

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

if (!SIG_COLUMN %in% c("padj", "pvalue")) {
  stop("RNA_PATTERN_SIGNIFICANCE_COLUMN must be 'padj' or 'pvalue'.")
}
if (!is.finite(PATTERN_ALPHA) || PATTERN_ALPHA <= 0 || PATTERN_ALPHA >= 1) {
  stop("RNA_PATTERN_ALPHA must be between 0 and 1.")
}

# ==============================================================================
# 1. Packages
# ==============================================================================

required_packages <- c(
  "dplyr", "tidyr", "tibble", "readr", "ggplot2", "ggrepel",
  "pheatmap", "RColorBrewer", "scales", "stringr"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ", paste(missing_packages, collapse = ", "),
    "\nInstall them first, for example:\ninstall.packages(c(",
    paste(sprintf("'%s'", missing_packages), collapse = ", "), "))"
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
})

# Optional Bioconductor packages for GO/KEGG enrichment.
# If absent, the script will still complete and will write a clear skip message.
HAS_CLUSTERPROFILER <- requireNamespace("clusterProfiler", quietly = TRUE)
HAS_ORG_GG <- requireNamespace("org.Gg.eg.db", quietly = TRUE)
HAS_ANNOTATIONDBI <- requireNamespace("AnnotationDbi", quietly = TRUE)

# ==============================================================================
# 2. General helper functions
# ==============================================================================

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

first_existing_file <- function(directory, candidates, required = TRUE) {
  paths <- file.path(directory, candidates)
  hit <- paths[file.exists(paths)]
  if (length(hit) > 0) return(normalizePath(hit[[1]], winslash = "/", mustWork = TRUE))
  if (required) {
    stop(
      "None of the expected files was found in: ", directory,
      "\nExpected one of: ", paste(candidates, collapse = ", ")
    )
  }
  NA_character_
}

read_table_auto <- function(path) {
  ext <- tolower(tools::file_ext(path))

  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("Package 'readxl' is required to read: ", path)
    }
    return(as.data.frame(readxl::read_excel(path, check.names = FALSE)))
  }

  if (ext == "csv") {
    return(as.data.frame(readr::read_csv(path, show_col_types = FALSE, name_repair = "minimal")))
  }

  # readr::read_delim with delimiter guessing handles tab/space/semicolon metadata.
  guessed <- tryCatch(
    suppressMessages(readr::read_delim(path, delim = NULL, show_col_types = FALSE, name_repair = "minimal")),
    error = function(e) NULL
  )
  if (!is.null(guessed) && ncol(guessed) >= 2) return(as.data.frame(guessed))

  as.data.frame(
    utils::read.table(
      path, header = TRUE, sep = "", quote = "\"", comment.char = "",
      check.names = FALSE, fileEncoding = "UTF-8"
    )
  )
}

clean_gene_id <- function(x) {
  x <- trimws(as.character(x))
  sub("\\.[0-9]+$", "", x)
}

find_col <- function(dat, candidates, required = TRUE, label = "column") {
  nm <- names(dat)
  nm_low <- tolower(nm)
  idx <- match(tolower(candidates), nm_low)
  idx <- idx[!is.na(idx)]
  if (length(idx) > 0) return(nm[idx[[1]]])
  if (required) {
    stop("Could not identify ", label, ". Available columns: ", paste(nm, collapse = ", "))
  }
  NA_character_
}

safe_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

theme_rna <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, color = "grey90"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(color = "grey30"),
      strip.background = element_rect(fill = "grey95", color = "grey70"),
      strip.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
}

save_plot <- function(plot_object, stem, width, height, dpi = 400) {
  ggplot2::ggsave(
    filename = file.path(OUTDIR, paste0(stem, ".png")), plot = plot_object,
    width = width, height = height, dpi = dpi, bg = "white"
  )
  ggplot2::ggsave(
    filename = file.path(OUTDIR, paste0(stem, ".pdf")), plot = plot_object,
    width = width, height = height, device = grDevices::cairo_pdf, bg = "white"
  )
}

write_csv_safe <- function(x, filename) {
  readr::write_excel_csv(as.data.frame(x), file.path(OUTDIR, filename), na = "")
}

zscore_rows <- function(mat) {
  ans <- t(scale(t(mat)))
  ans[!is.finite(ans)] <- 0
  ans
}


# ==============================================================================
# 3. Locate and read inputs
# ==============================================================================

VST_FILE <- first_existing_file(
  CORE_OUTDIR,
  c(
    "Table3_VST_expression_matrix.csv",
    "VST_expression_matrix.csv",
    "Table3_vst_expression_matrix.csv"
  )
)

ANNOT_FILE <- first_existing_file(
  CORE_OUTDIR,
  c(
    "Table2_gene_annotation_used.csv",
    "gene_annotation_used.csv",
    "Table2_gene_annotation.csv"
  ),
  required = FALSE
)

DEG_FILES <- c(
  YY_vs_XX = first_existing_file(CORE_OUTDIR, c("DEG_YY_vs_XX.csv", "YY_vs_XX.csv")),
  YX_vs_YY = first_existing_file(CORE_OUTDIR, c("DEG_YX_vs_YY.csv", "YX_vs_YY.csv")),
  YX_vs_XX = first_existing_file(CORE_OUTDIR, c("DEG_YX_vs_XX.csv", "YX_vs_XX.csv"))
)

message("Reading VST matrix: ", VST_FILE)
vst_raw <- read_table_auto(VST_FILE)

gene_col_vst <- find_col(
  vst_raw,
  c("gene_id", "gene", "ensembl_gene_id", "Geneid", "X"),
  required = FALSE,
  label = "gene-ID column in VST matrix"
)
if (is.na(gene_col_vst)) gene_col_vst <- names(vst_raw)[[1]]

gene_ids <- clean_gene_id(vst_raw[[gene_col_vst]])
vst_numeric <- vst_raw[, setdiff(names(vst_raw), gene_col_vst), drop = FALSE]
vst_numeric[] <- lapply(vst_numeric, safe_numeric)
vst_mat <- as.matrix(vst_numeric)
rownames(vst_mat) <- gene_ids

# Remove invalid or duplicate gene rows. For duplicate IDs retain the row with
# the largest expression variance.
valid_gene <- !is.na(rownames(vst_mat)) & rownames(vst_mat) != ""
vst_mat <- vst_mat[valid_gene, , drop = FALSE]
if (anyDuplicated(rownames(vst_mat))) {
  gene_var <- apply(vst_mat, 1, stats::var, na.rm = TRUE)
  ord <- order(gene_var, decreasing = TRUE, na.last = TRUE)
  vst_mat <- vst_mat[ord, , drop = FALSE]
  vst_mat <- vst_mat[!duplicated(rownames(vst_mat)), , drop = FALSE]
}

if (anyNA(vst_mat)) {
  stop("The VST expression matrix contains non-numeric or missing expression values.")
}

message("Reading metadata: ", META_FILE)
meta_raw <- read_table_auto(META_FILE)

sample_col <- find_col(
  meta_raw,
  c("sample_id", "sample", "Sample", "sampleID", "id", "ID"),
  label = "sample-ID column in metadata"
)
group_col <- find_col(
  meta_raw,
  c("group", "Group", "treatment", "Treatment", "breed", "Breed"),
  label = "group column in metadata"
)

meta <- meta_raw %>%
  transmute(
    sample_id = trimws(as.character(.data[[sample_col]])),
    group = toupper(trimws(as.character(.data[[group_col]])))
  ) %>%
  filter(!is.na(sample_id), sample_id != "") %>%
  distinct(sample_id, .keep_all = TRUE)

meta <- meta %>% mutate(plot_label = sample_id)

unexpected_groups <- setdiff(unique(meta$group), c("YY", "YX", "XX"))
if (length(unexpected_groups) > 0) {
  stop("Unexpected metadata groups: ", paste(unexpected_groups, collapse = ", "),
       ". Expected exactly YY, YX and XX.")
}

# ==============================================================================
# 4. Match samples and retain all available samples
# ==============================================================================

matrix_samples_initial <- colnames(vst_mat)
metadata_samples_initial <- meta$sample_id

missing_in_metadata <- setdiff(matrix_samples_initial, metadata_samples_initial)
if (length(missing_in_metadata) > 0) {
  stop(
    "VST samples missing from metadata: ",
    paste(missing_in_metadata, collapse = ", ")
  )
}

missing_in_vst <- setdiff(metadata_samples_initial, matrix_samples_initial)
if (length(missing_in_vst) > 0) {
  message(
    "Metadata samples absent from the upstream VST matrix and unavailable for this analysis: ",
    paste(missing_in_vst, collapse = ", ")
  )
}

# Keep a complete metadata copy for QC reporting, then retain every sample that
# is available in the upstream VST matrix. No additional sample exclusion is
# performed in this script.
meta_all <- meta
meta <- meta_all %>% filter(sample_id %in% matrix_samples_initial)
message("No sample exclusion is applied; using all samples matched between metadata and the VST matrix.")

# Reorder VST columns to metadata order.
meta <- meta %>%
  mutate(group = factor(group, levels = c("YY", "YX", "XX"))) %>%
  arrange(group, sample_id)
vst_mat <- vst_mat[, meta$sample_id, drop = FALSE]

group_n <- table(meta$group)
if (any(group_n == 0)) stop("At least one required group has zero matched samples.")

sample_qc <- tibble(
  sample_id = meta$sample_id,
  plot_label = meta$plot_label,
  group = as.character(meta$group),
  present_in_upstream_vst = TRUE,
  used_in_this_analysis = TRUE
)

upstream_absent_qc <- meta_all %>%
  filter(sample_id %in% missing_in_vst) %>%
  transmute(
    sample_id,
    plot_label,
    group = as.character(group),
    present_in_upstream_vst = FALSE,
    used_in_this_analysis = FALSE
  )

write_csv_safe(bind_rows(sample_qc, upstream_absent_qc), "QC_samples_used_and_upstream_absent.csv")

message("Samples used: YY = ", group_n[["YY"]],
        ", YX = ", group_n[["YX"]],
        ", XX = ", group_n[["XX"]])

# ==============================================================================
# 5. Read gene annotation
# ==============================================================================

if (!is.na(ANNOT_FILE)) {
  message("Reading annotation: ", ANNOT_FILE)
  annot_raw <- read_table_auto(ANNOT_FILE)
  annot_gene_col <- find_col(
    annot_raw,
    c("gene_id", "gene", "ensembl_gene_id", "Geneid"),
    required = FALSE,
    label = "gene-ID column in annotation"
  )
  if (is.na(annot_gene_col)) annot_gene_col <- names(annot_raw)[[1]]

  symbol_col <- find_col(
    annot_raw,
    c("symbol", "gene_symbol", "external_gene_name", "SYMBOL"),
    required = FALSE
  )
  gene_name_col <- find_col(
    annot_raw,
    c("gene_name", "description", "gene_description", "GENENAME"),
    required = FALSE
  )
  entrez_col <- find_col(
    annot_raw,
    c("entrez", "entrez_id", "ENTREZID", "ncbi_gene_id"),
    required = FALSE
  )

  annot <- tibble(
    gene_id = clean_gene_id(annot_raw[[annot_gene_col]]),
    symbol = if (!is.na(symbol_col)) as.character(annot_raw[[symbol_col]]) else NA_character_,
    gene_name = if (!is.na(gene_name_col)) as.character(annot_raw[[gene_name_col]]) else NA_character_,
    entrez = if (!is.na(entrez_col)) as.character(annot_raw[[entrez_col]]) else NA_character_
  ) %>%
    filter(!is.na(gene_id), gene_id != "") %>%
    distinct(gene_id, .keep_all = TRUE)
} else {
  message("Annotation file not found; gene IDs will be used as plot labels.")
  annot <- tibble(
    gene_id = rownames(vst_mat), symbol = NA_character_,
    gene_name = NA_character_, entrez = NA_character_
  )
}

# ==============================================================================
# 6. Read and standardize pairwise DESeq2 result files
# ==============================================================================

standardize_deg <- function(path, contrast_name) {
  dat <- read_table_auto(path)
  gene_col <- find_col(
    dat,
    c("gene_id", "gene", "ensembl_gene_id", "Geneid", "X"),
    required = FALSE,
    label = paste0("gene-ID column in ", basename(path))
  )
  if (is.na(gene_col)) gene_col <- names(dat)[[1]]

  lfc_col <- find_col(dat, c("log2FoldChange", "log2FC", "logFC"), label = "log2 fold-change column")
  p_col <- find_col(dat, c("pvalue", "PValue", "p_value", "P.Value"), required = FALSE)
  padj_col <- find_col(dat, c("padj", "FDR", "qvalue", "adj.P.Val"), required = FALSE)

  if (SIG_COLUMN == "pvalue" && is.na(p_col)) {
    stop("pvalue column is absent from: ", path)
  }
  if (SIG_COLUMN == "padj" && is.na(padj_col)) {
    stop("padj/FDR column is absent from: ", path)
  }

  out <- tibble(
    gene_id = clean_gene_id(dat[[gene_col]]),
    log2FoldChange = safe_numeric(dat[[lfc_col]]),
    pvalue = if (!is.na(p_col)) safe_numeric(dat[[p_col]]) else NA_real_,
    padj = if (!is.na(padj_col)) safe_numeric(dat[[padj_col]]) else NA_real_
  ) %>%
    filter(!is.na(gene_id), gene_id != "") %>%
    distinct(gene_id, .keep_all = TRUE)

  sig_value <- if (SIG_COLUMN == "padj") out$padj else out$pvalue
  out$significant <- !is.na(sig_value) & sig_value < PATTERN_ALPHA &
    !is.na(out$log2FoldChange) & abs(out$log2FoldChange) >= PAIRWISE_LFC_CUTOFF

  names(out)[-1] <- paste0(names(out)[-1], "_", contrast_name)
  out
}

message("Reading pairwise DESeq2 results...")
deg_yy_xx <- standardize_deg(DEG_FILES[["YY_vs_XX"]], "YY_vs_XX")
deg_yx_yy <- standardize_deg(DEG_FILES[["YX_vs_YY"]], "YX_vs_YY")
deg_yx_xx <- standardize_deg(DEG_FILES[["YX_vs_XX"]], "YX_vs_XX")

# ==============================================================================
# 7. Calculate group means, F1 position and expression-inheritance patterns
# ==============================================================================

group_mean <- function(group_name) {
  samples <- meta$sample_id[meta$group == group_name]
  rowMeans(vst_mat[, samples, drop = FALSE], na.rm = TRUE)
}

mean_YY <- group_mean("YY")
mean_YX <- group_mean("YX")
mean_XX <- group_mean("XX")

pattern_df <- tibble(
  gene_id = rownames(vst_mat),
  mean_YY = mean_YY[rownames(vst_mat)],
  mean_YX = mean_YX[rownames(vst_mat)],
  mean_XX = mean_XX[rownames(vst_mat)]
) %>%
  mutate(
    mid_parent = (mean_YY + mean_XX) / 2,
    F1_deviation_from_mid = mean_YX - mid_parent,
    abs_mid_deviation = abs(F1_deviation_from_mid),
    parental_difference = mean_YY - mean_XX,
    parental_range = abs(mean_YY - mean_XX),
    F1_position = if_else(
      abs(parental_difference) > 1e-8,
      (mean_YX - mean_XX) / parental_difference,
      NA_real_
    ),
    lower_parent = pmin(mean_YY, mean_XX),
    upper_parent = pmax(mean_YY, mean_XX)
  ) %>%
  left_join(annot, by = "gene_id") %>%
  left_join(deg_yy_xx, by = "gene_id") %>%
  left_join(deg_yx_yy, by = "gene_id") %>%
  left_join(deg_yx_xx, by = "gene_id")

# Missing DE-result rows are treated as not significant, while their statistics
# remain NA in the output table.
pattern_df <- pattern_df %>%
  mutate(
    sig_YY_vs_XX = coalesce(significant_YY_vs_XX, FALSE),
    sig_YX_vs_YY = coalesce(significant_YX_vs_YY, FALSE),
    sig_YX_vs_XX = coalesce(significant_YX_vs_XX, FALSE),
    DEG_any = sig_YY_vs_XX | sig_YX_vs_YY | sig_YX_vs_XX,
    DEG_F1_vs_any_parent = sig_YX_vs_YY | sig_YX_vs_XX,
    within_parental_range = mean_YX >= lower_parent - POSITION_TOLERANCE &
      mean_YX <= upper_parent + POSITION_TOLERANCE,
    above_parental_range = mean_YX > upper_parent + POSITION_TOLERANCE,
    below_parental_range = mean_YX < lower_parent - POSITION_TOLERANCE,

    # Classification rules:
    # - Transgressive: outside both parents and significantly different from both.
    # - Dominant: parents differ; YX resembles one parent (not significant) but
    #   differs from the other parent.
    # - Intermediate: parents differ; YX is inside their range and differs from both.
    # - The rest remain unresolved rather than being forced into a class.
    pattern = case_when(
      above_parental_range & sig_YX_vs_YY & sig_YX_vs_XX ~ "YX_high_transgressive",
      below_parental_range & sig_YX_vs_YY & sig_YX_vs_XX ~ "YX_low_transgressive",
      sig_YY_vs_XX & !sig_YX_vs_YY & sig_YX_vs_XX ~ "YY_dominant",
      sig_YY_vs_XX & sig_YX_vs_YY & !sig_YX_vs_XX ~ "XX_dominant",
      sig_YY_vs_XX & within_parental_range & sig_YX_vs_YY & sig_YX_vs_XX ~ "YX_intermediate",
      !sig_YY_vs_XX ~ "parental_similar_or_unresolved",
      TRUE ~ "unclassified"
    ),
    inheritance_class = case_when(
      pattern == "YX_intermediate" ~ "additive_like",
      pattern %in% c("YY_dominant", "XX_dominant") ~ "dominance_like",
      pattern %in% c("YX_high_transgressive", "YX_low_transgressive") ~ "overdominance_like",
      TRUE ~ "unresolved"
    ),
    gene_label = case_when(
      !is.na(symbol) & trimws(symbol) != "" & symbol != gene_id ~ symbol,
      TRUE ~ gene_id
    )
  )

pattern_levels_all <- c(
  "YX_intermediate", "YY_dominant", "XX_dominant",
  "YX_low_transgressive", "YX_high_transgressive",
  "parental_similar_or_unresolved", "unclassified"
)

pattern_levels_main <- c(
  "YX_intermediate", "YY_dominant", "XX_dominant",
  "YX_low_transgressive", "YX_high_transgressive"
)

pattern_label_map <- c(
  "YX_intermediate" = "Intermediate",
  "YY_dominant" = "YY-dominant",
  "XX_dominant" = "XX-dominant",
  "YX_low_transgressive" = "Low transgressive",
  "YX_high_transgressive" = "High transgressive",
  "parental_similar_or_unresolved" = "Parental similar / unresolved",
  "unclassified" = "Unclassified"
)

PATTERN_COLORS <- c(
  "YX_intermediate" = "#4DBBD5",
  "YY_dominant" = "#E64B35",
  "XX_dominant" = "#3C5488",
  "YX_low_transgressive" = "#00A087",
  "YX_high_transgressive" = "#F39B7F",
  "parental_similar_or_unresolved" = "#BDBDBD",
  "unclassified" = "#737373"
)

# Defensive conversion: avoid failures if an upstream object was accidentally stored as a list.
if (is.list(pattern_df) && !inherits(pattern_df, "data.frame")) {
  if (length(pattern_df) > 0 && all(vapply(pattern_df, function(z) inherits(z, "data.frame"), logical(1)))) {
    pattern_df <- dplyr::bind_rows(pattern_df)
  } else {
    stop("pattern_df is a list, not a data.frame/tibble. Please check the object construction before pattern summary.")
  }
}
pattern_df <- tibble::as_tibble(pattern_df)

write_csv_safe(pattern_df, "Table_F1_expression_pattern_all_genes.csv")

# ==============================================================================
# 8. Correct pattern summaries and main count plot
# ==============================================================================

pattern_summary_all <- pattern_df %>%
  mutate(pattern = factor(pattern, levels = pattern_levels_all)) %>%
  dplyr::count(pattern, .drop = FALSE, name = "N_genes") %>%
  mutate(
    pattern = as.character(pattern),
    pattern_label = unname(pattern_label_map[pattern]),
    N_total_genes = sum(N_genes),
    Percent_all_genes = if_else(N_total_genes > 0, N_genes / N_total_genes * 100, 0)
  )

n_total_genes <- sum(pattern_summary_all$N_genes)
n_classified_total <- pattern_summary_all %>%
  filter(pattern %in% pattern_levels_main) %>%
  summarise(n = sum(N_genes)) %>%
  pull(n)
n_unresolved <- n_total_genes - n_classified_total

pattern_summary_main <- pattern_summary_all %>%
  filter(pattern %in% pattern_levels_main) %>%
  mutate(
    pattern = factor(pattern, levels = pattern_levels_main),
    N_classified_total = n_classified_total,
    # Critical bug fix: calculate each row's own percentage.
    Percent_classified = if_else(
      N_classified_total > 0,
      N_genes / N_classified_total * 100,
      0
    ),
    label_text = sprintf("%s\n%.1f%%", scales::comma(N_genes), Percent_classified)
  ) %>%
  arrange(pattern)

# Internal validation prevents recurrence of the repeated-percentage plotting bug.
if (n_classified_total > 0) {
  if (abs(sum(pattern_summary_main$Percent_classified) - 100) > 0.05) {
    stop("Internal error: classified pattern percentages do not sum to 100%.")
  }
  expected_percent <- pattern_summary_main$N_genes / n_classified_total * 100
  if (!isTRUE(all.equal(pattern_summary_main$Percent_classified, expected_percent, tolerance = 1e-10))) {
    stop("Internal error: Percent_classified does not match N_genes / classified total.")
  }
}

write_csv_safe(pattern_summary_all, "Table_F1_expression_pattern_summary_all.csv")
write_csv_safe(pattern_summary_main, "Table_F1_expression_pattern_summary_classified_only.csv")

message("\nPattern percentage check (denominator = ", n_classified_total, " classifiable genes):")
print(
  pattern_summary_main %>%
    dplyr::select(pattern, N_genes, N_classified_total, Percent_classified, Percent_all_genes)
)

unresolved_percent <- if (n_total_genes > 0) n_unresolved / n_total_genes * 100 else 0

p_pattern_main <- ggplot(
  pattern_summary_main,
  aes(x = reorder(pattern_label, N_genes), y = N_genes, fill = pattern)
) +
  geom_col(width = 0.70, show.legend = FALSE) +
  geom_text(aes(label = label_text), hjust = -0.06, size = 3.9, lineheight = 0.9) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = PATTERN_COLORS, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.22))) +
  labs(
    title = paste0("F1 expression inheritance patterns: ", TISSUE_LABEL),
    subtitle = paste0(
      "Percentages use classifiable genes only (n = ", scales::comma(n_classified_total),
      "); unresolved/unclassified = ", scales::comma(n_unresolved),
      " (", sprintf("%.1f", unresolved_percent), "% of all genes)"
    ),
    x = NULL,
    y = "Number of genes"
  ) +
  theme_rna(11) +
  theme(plot.margin = margin(8, 38, 8, 8))

save_plot(
  p_pattern_main,
  "Fig_F1_expression_pattern_counts_for_manuscript",
  width = 8.4,
  height = 5.0
)

# ==============================================================================
# 9. Inheritance-class summary
# ==============================================================================

inheritance_levels <- c("additive_like", "dominance_like", "overdominance_like", "unresolved")
inheritance_colors <- c(
  "additive_like" = "#4DBBD5",
  "dominance_like" = "#7E6148",
  "overdominance_like" = "#E64B35",
  "unresolved" = "#BDBDBD"
)

inheritance_summary <- pattern_df %>%
  mutate(inheritance_class = factor(inheritance_class, levels = inheritance_levels)) %>%
  dplyr::count(inheritance_class, .drop = FALSE, name = "N_genes") %>%
  mutate(
    Percent_all_genes = N_genes / sum(N_genes) * 100,
    label_text = sprintf("%s (%.1f%%)", scales::comma(N_genes), Percent_all_genes)
  )

write_csv_safe(inheritance_summary, "Table_expression_inheritance_class_summary.csv")

p_inheritance <- ggplot(
  inheritance_summary,
  aes(x = inheritance_class, y = N_genes, fill = inheritance_class)
) +
  geom_col(width = 0.68, show.legend = FALSE) +
  geom_text(aes(label = label_text), vjust = -0.35, size = 3.5) +
  scale_fill_manual(values = inheritance_colors, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.14))) +
  labs(
    title = paste0("Expression inheritance classes: ", TISSUE_LABEL),
    x = NULL,
    y = "Number of genes"
  ) +
  theme_rna(11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

save_plot(p_inheritance, "Fig_expression_inheritance_class_summary", 7.3, 5.0)

# ==============================================================================
# 10. F1 position distribution
# ==============================================================================

position_plot_df <- pattern_df %>%
  filter(pattern %in% pattern_levels_main, is.finite(F1_position)) %>%
  mutate(
    pattern = factor(pattern, levels = pattern_levels_main),
    F1_position_plot = pmax(pmin(F1_position, 2.5), -1.5)
  )

p_position <- ggplot(position_plot_df, aes(x = F1_position_plot, fill = pattern, color = pattern)) +
  geom_density(alpha = 0.18, linewidth = 0.8, adjust = 1.1) +
  geom_vline(xintercept = c(0, 0.5, 1), linetype = c("dotted", "dashed", "dotted"),
             color = c("grey35", "black", "grey35"), linewidth = c(0.5, 0.6, 0.5)) +
  scale_fill_manual(values = PATTERN_COLORS, drop = FALSE) +
  scale_color_manual(values = PATTERN_COLORS, drop = FALSE) +
  scale_x_continuous(breaks = seq(-1.5, 2.5, 0.5)) +
  labs(
    title = "F1 expression position relative to parental means",
    subtitle = "XX = 0, mid-parent = 0.5, YY = 1; extreme values truncated only for display",
    x = "F1 position: (YX - XX) / (YY - XX)",
    y = "Density",
    color = "Pattern",
    fill = "Pattern"
  ) +
  theme_rna(10) +
  theme(legend.position = "right")

save_plot(p_position, "Fig_F1_expression_position_distribution", 8.3, 5.6)

# ==============================================================================
# 11. Mid-parent-deviation plot with selected gene labels
# ==============================================================================

duodenum_candidate_genes <- c(
  # Lipid digestion, absorption and lipoprotein assembly
  "APOA1", "APOA4", "APOB", "APOC3", "MTTP", "SAR1B", "MOGAT2", "DGAT1", "DGAT2",
  "CD36", "SCARB1", "NPC1L1", "FABP1", "FABP2", "FABP3", "FABP4", "PLIN2",
  # Fatty-acid uptake, activation and oxidation
  "SLC27A1", "SLC27A2", "SLC27A4", "ACSL1", "ACSL3", "ACSL5", "CPT1A", "CPT2",
  "ACOX1", "ACOX2", "HADHA", "HADHB", "PPARA", "PPARD", "PPARG",
  # Lipogenesis and lipid transport
  "FASN", "ACACA", "SCD", "LPL", "ABCA1", "ABCG5", "ABCG8",
  # Bile-acid transport and intestinal barrier
  "SLC10A2", "SLC51A", "SLC51B", "FABP6", "FXR", "NR1H4", "TJP1", "OCLN", "CLDN1",
  "CLDN2", "CLDN3", "CLDN4", "MUC2"
)

plot_mid_df <- pattern_df %>%
  filter(pattern %in% pattern_levels_main) %>%
  mutate(pattern = factor(pattern, levels = pattern_levels_main))

# 1) Pairwise DEGs with the largest absolute mid-parent deviation.
label_deg <- plot_mid_df %>%
  filter(DEG_any) %>%
  arrange(desc(abs_mid_deviation)) %>%
  slice_head(n = max(0, LABEL_TOP_N_DEG)) %>%
  mutate(label_reason = "pairwise_DEG")

# 2) Genes with the largest absolute mid-parent deviation.
label_deviation <- plot_mid_df %>%
  arrange(desc(abs_mid_deviation)) %>%
  slice_head(n = max(0, LABEL_TOP_N_DEVIATION)) %>%
  mutate(label_reason = "largest_mid_parent_deviation")

# 3) Prespecified duodenum candidates, prioritizing large deviations.
label_candidate <- plot_mid_df[0, ] %>% mutate(label_reason = character())
if (LABEL_CANDIDATES) {
  label_candidate <- plot_mid_df %>%
    filter(toupper(gene_label) %in% duodenum_candidate_genes) %>%
    arrange(desc(DEG_any), desc(abs_mid_deviation)) %>%
    slice_head(n = 12) %>%
    mutate(label_reason = "duodenum_candidate")
}

label_df <- bind_rows(label_deg, label_deviation, label_candidate) %>%
  arrange(desc(DEG_any), desc(abs_mid_deviation)) %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  slice_head(n = max(1, MAX_TOTAL_LABELS))

# Make duplicate symbols unique in labels without changing the annotation table.
label_df <- label_df %>%
  group_by(gene_label) %>%
  mutate(plot_gene_label = if_else(n() > 1, paste0(gene_label, " [", gene_id, "]"), gene_label)) %>%
  ungroup()

write_csv_safe(
  label_df %>%
    dplyr::select(
      gene_id, symbol, gene_name, plot_gene_label, label_reason, pattern,
      mean_YY, mean_YX, mean_XX, parental_difference,
      F1_deviation_from_mid, F1_position,
      DEG_any, DEG_F1_vs_any_parent,
      sig_YY_vs_XX, sig_YX_vs_YY, sig_YX_vs_XX,
      padj_YY_vs_XX, padj_YX_vs_YY, padj_YX_vs_XX,
      pvalue_YY_vs_XX, pvalue_YX_vs_YY, pvalue_YX_vs_XX
    ),
  "Table_F1_mid_parent_deviation_labeled_genes.csv"
)

p_mid_dev <- ggplot(
  plot_mid_df,
  aes(x = parental_difference, y = F1_deviation_from_mid, color = pattern)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.45) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.45) +
  geom_point(size = 1.15, alpha = 0.48) +
  ggrepel::geom_text_repel(
    data = label_df,
    aes(label = plot_gene_label),
    size = 2.8,
    max.overlaps = Inf,
    box.padding = 0.35,
    point.padding = 0.18,
    min.segment.length = 0,
    segment.size = 0.25,
    segment.alpha = 0.65,
    seed = 20260626,
    show.legend = FALSE
  ) +
  scale_color_manual(values = PATTERN_COLORS, drop = FALSE) +
  labs(
    title = "F1 deviation from mid-parent expression",
    subtitle = paste0(
      "Labels: pairwise DEGs and genes with large deviation",
      if (LABEL_CANDIDATES) ", plus selected duodenum candidates" else "",
      "; significance = ", SIG_COLUMN, " < ", PATTERN_ALPHA
    ),
    x = "Parental expression difference: YY - XX (VST)",
    y = "F1 deviation: YX - (YY + XX)/2 (VST)",
    color = "Pattern"
  ) +
  theme_rna(10) +
  theme(legend.position = "right")

save_plot(p_mid_dev, "Fig_F1_mid_parent_deviation_labeled", 9.3, 6.7)

# Also save an unlabelled version for dense datasets.
p_mid_dev_plain <- p_mid_dev +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2.2)))
p_mid_dev_plain$layers <- p_mid_dev_plain$layers[1:3]
save_plot(p_mid_dev_plain, "Fig_F1_mid_parent_deviation_plain", 8.0, 5.8)

# ==============================================================================
# 12. Select top genes by pattern
# ==============================================================================

# Ranking reflects strength and interpretability of each pattern.
top_pattern_genes <- pattern_df %>%
  filter(pattern %in% pattern_levels_main) %>%
  mutate(
    pattern_rank_score = case_when(
      pattern == "YX_intermediate" ~ abs(parental_difference),
      pattern %in% c("YY_dominant", "XX_dominant") ~ abs(parental_difference),
      pattern %in% c("YX_high_transgressive", "YX_low_transgressive") ~ abs_mid_deviation,
      TRUE ~ abs_mid_deviation
    )
  ) %>%
  group_by(pattern) %>%
  arrange(desc(pattern_rank_score), .by_group = TRUE) %>%
  slice_head(n = TOP_N_PER_PATTERN) %>%
  ungroup() %>%
  arrange(factor(pattern, levels = pattern_levels_main), desc(pattern_rank_score))

write_csv_safe(top_pattern_genes, "Table_F1_top_genes_by_expression_pattern.csv")

nonadditive_genes <- pattern_df %>%
  filter(inheritance_class %in% c("dominance_like", "overdominance_like")) %>%
  arrange(factor(inheritance_class, levels = c("overdominance_like", "dominance_like")),
          desc(abs_mid_deviation))
write_csv_safe(nonadditive_genes, "Table_nonadditive_like_expression_genes.csv")

# ==============================================================================
# 13. Heatmaps of selected pattern genes
# ==============================================================================

heatmap_genes <- unique(top_pattern_genes$gene_id)
if (length(heatmap_genes) > HEATMAP_MAX_GENES) {
  heatmap_genes <- heatmap_genes[seq_len(HEATMAP_MAX_GENES)]
}

if (length(heatmap_genes) >= 2) {
  heatmap_mat <- vst_mat[heatmap_genes, , drop = FALSE]
  heatmap_mat_z <- zscore_rows(heatmap_mat)

  heatmap_annot <- top_pattern_genes %>%
    filter(gene_id %in% heatmap_genes) %>%
    distinct(gene_id, .keep_all = TRUE)

  row_labels <- heatmap_annot$gene_label[match(rownames(heatmap_mat_z), heatmap_annot$gene_id)]
  row_labels[is.na(row_labels) | row_labels == ""] <- rownames(heatmap_mat_z)[is.na(row_labels) | row_labels == ""]
  rownames(heatmap_mat_z) <- make.unique(row_labels)

  annotation_col <- data.frame(Group = meta$group)
  rownames(annotation_col) <- meta$sample_id

  annotation_row <- data.frame(
    Pattern = heatmap_annot$pattern[match(heatmap_genes, heatmap_annot$gene_id)]
  )
  rownames(annotation_row) <- make.unique(row_labels)

  ann_colors <- list(
    Group = c(YY = "#E64B35", YX = "#4DBBD5", XX = "#3C5488"),
    Pattern = PATTERN_COLORS[pattern_levels_main]
  )

  png(file.path(OUTDIR, "Fig_F1_top_pattern_genes_heatmap_samples.png"),
      width = 2600, height = max(2200, 48 * nrow(heatmap_mat_z)), res = 300)
  pheatmap::pheatmap(
    heatmap_mat_z,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = annotation_col,
    annotation_row = annotation_row,
    annotation_colors = ann_colors,
    color = colorRampPalette(c("#3C5488", "white", "#E64B35"))(101),
    breaks = seq(-2.5, 2.5, length.out = 102),
    fontsize_row = 7,
    fontsize_col = 8,
    border_color = NA,
    main = paste0("Top F1 expression-pattern genes: ", TISSUE_LABEL)
  )
  dev.off()

  pdf(file.path(OUTDIR, "Fig_F1_top_pattern_genes_heatmap_samples.pdf"),
      width = 9.0, height = max(7.0, 0.16 * nrow(heatmap_mat_z) + 2.5), useDingbats = FALSE)
  pheatmap::pheatmap(
    heatmap_mat_z,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = annotation_col,
    annotation_row = annotation_row,
    annotation_colors = ann_colors,
    color = colorRampPalette(c("#3C5488", "white", "#E64B35"))(101),
    breaks = seq(-2.5, 2.5, length.out = 102),
    fontsize_row = 7,
    fontsize_col = 8,
    border_color = NA,
    main = paste0("Top F1 expression-pattern genes: ", TISSUE_LABEL)
  )
  dev.off()

  group_mean_mat <- cbind(
    YY = pattern_df$mean_YY[match(heatmap_genes, pattern_df$gene_id)],
    YX = pattern_df$mean_YX[match(heatmap_genes, pattern_df$gene_id)],
    XX = pattern_df$mean_XX[match(heatmap_genes, pattern_df$gene_id)]
  )
  rownames(group_mean_mat) <- make.unique(row_labels)
  group_mean_z <- zscore_rows(group_mean_mat)

  png(file.path(OUTDIR, "Fig_F1_top_pattern_genes_heatmap_group_mean.png"),
      width = 1700, height = max(2200, 48 * nrow(group_mean_z)), res = 300)
  pheatmap::pheatmap(
    group_mean_z,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    annotation_row = annotation_row,
    annotation_colors = ann_colors,
    color = colorRampPalette(c("#3C5488", "white", "#E64B35"))(101),
    breaks = seq(-2.5, 2.5, length.out = 102),
    fontsize_row = 7,
    fontsize_col = 10,
    border_color = NA,
    main = "Group-mean expression of top pattern genes"
  )
  dev.off()

  pdf(file.path(OUTDIR, "Fig_F1_top_pattern_genes_heatmap_group_mean.pdf"),
      width = 6.2, height = max(7.0, 0.16 * nrow(group_mean_z) + 2.5), useDingbats = FALSE)
  pheatmap::pheatmap(
    group_mean_z,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    annotation_row = annotation_row,
    annotation_colors = ann_colors,
    color = colorRampPalette(c("#3C5488", "white", "#E64B35"))(101),
    breaks = seq(-2.5, 2.5, length.out = 102),
    fontsize_row = 7,
    fontsize_col = 10,
    border_color = NA,
    main = "Group-mean expression of top pattern genes"
  )
  dev.off()
} else {
  message("Fewer than two top pattern genes; heatmaps were skipped.")
}

# ==============================================================================
# 14. Candidate-gene summary for duodenum
# ==============================================================================

candidate_summary <- pattern_df %>%
  filter(toupper(gene_label) %in% duodenum_candidate_genes) %>%
  arrange(factor(pattern, levels = pattern_levels_all), desc(abs_mid_deviation))

write_csv_safe(candidate_summary, "Table_duodenum_candidate_genes_in_F1_patterns.csv")

candidate_pattern_summary <- candidate_summary %>%
  mutate(pattern = factor(pattern, levels = pattern_levels_all)) %>%
  dplyr::count(pattern, .drop = FALSE, name = "N_candidate_genes") %>%
  mutate(pattern_label = unname(pattern_label_map[as.character(pattern)]))

write_csv_safe(candidate_pattern_summary, "Table_duodenum_candidate_gene_pattern_summary.csv")

if (sum(candidate_pattern_summary$N_candidate_genes) > 0) {
  p_candidate <- candidate_pattern_summary %>%
    filter(N_candidate_genes > 0) %>%
    ggplot(aes(x = reorder(pattern_label, N_candidate_genes), y = N_candidate_genes,
               fill = as.character(pattern))) +
    geom_col(width = 0.7, show.legend = FALSE) +
    geom_text(aes(label = N_candidate_genes), hjust = -0.15, size = 3.7) +
    coord_flip(clip = "off") +
    scale_fill_manual(values = PATTERN_COLORS, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.18))) +
    labs(
      title = "Duodenum candidate genes across F1 expression patterns",
      x = NULL,
      y = "Number of matched candidate genes"
    ) +
    theme_rna(10) +
    theme(plot.margin = margin(8, 25, 8, 8))

  save_plot(p_candidate, "Fig_duodenum_candidate_genes_in_F1_patterns", 7.5, 5.0)
}


# ==============================================================================
# 15. Enrichment analysis for classified F1-pattern gene sets
# ==============================================================================

# Biological purpose:
#   The previous sections assign each gene to an expression-inheritance pattern.
#   This section treats those categories as gene sets and tests whether each set
#   is enriched for GO biological process, GO cellular component, GO molecular
#   function, KEGG pathways, or curated duodenum lipid-handling modules.
#
# Statistical note:
#   - GO/KEGG ORA uses genes with Entrez IDs as the testable universe.
#   - Curated-module enrichment uses all classified/tested genes with gene labels
#     as the universe and therefore remains useful even when Entrez annotation is
#     incomplete.
#   - Results from small gene sets, especially transgressive classes, should be
#     treated as exploratory and interpreted together with the actual overlapping
#     genes, not only by p.adjust.

parse_gene_ratio <- function(x) {
  x <- as.character(x)
  out <- rep(NA_real_, length(x))
  has_slash <- grepl("/", x)
  if (any(has_slash)) {
    sp <- strsplit(x[has_slash], "/", fixed = TRUE)
    out[has_slash] <- vapply(sp, function(z) {
      if (length(z) != 2) return(NA_real_)
      suppressWarnings(as.numeric(z[1]) / as.numeric(z[2]))
    }, numeric(1))
  }
  if (any(!has_slash)) {
    out[!has_slash] <- suppressWarnings(as.numeric(x[!has_slash]))
  }
  out
}

shorten_text <- function(x, width = 70) {
  x <- as.character(x)
  ifelse(nchar(x) > width, paste0(substr(x, 1, width - 3), "..."), x)
}

get_first_nonempty <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) NA_character_ else as.character(x[[1]])
}

# Initialize objects so the Excel workbook can be written even if enrichment is skipped.
ora_pattern_kegg <- tibble()
ora_pattern_go_bp <- tibble()
ora_pattern_go_cc <- tibble()
ora_pattern_go_mf <- tibble()
ora_class_kegg <- tibble()
ora_class_go_bp <- tibble()
ora_class_go_cc <- tibble()
ora_class_go_mf <- tibble()
module_pattern_enrichment <- tibble()
module_class_enrichment <- tibble()
enrichment_input_summary <- tibble()

prepare_entrez_for_enrichment <- function(pattern_table) {
  out <- pattern_table %>%
    dplyr::mutate(
      entrez_enrich = as.character(entrez),
      entrez_enrich = dplyr::if_else(is.na(entrez_enrich) | entrez_enrich == "", NA_character_, entrez_enrich)
    )

  # Fill missing Entrez IDs from org.Gg.eg.db when possible.
  if (HAS_ORG_GG && HAS_ANNOTATIONDBI) {
    suppressPackageStartupMessages({
      library(org.Gg.eg.db)
      library(AnnotationDbi)
    })

    # 1) Try Ensembl gene ID mapping.
    missing_entrez <- out %>% dplyr::filter(is.na(entrez_enrich)) %>% dplyr::pull(gene_id) %>% unique()
    keytypes_available <- tryCatch(AnnotationDbi::keytypes(org.Gg.eg.db::org.Gg.eg.db), error = function(e) character())
    columns_available <- tryCatch(AnnotationDbi::columns(org.Gg.eg.db::org.Gg.eg.db), error = function(e) character())

    if (length(missing_entrez) > 0 && "ENTREZID" %in% columns_available) {
      for (kt in intersect(c("ENSEMBL", "ENSEMBLGENE"), keytypes_available)) {
        mp <- tryCatch(
          AnnotationDbi::select(
            org.Gg.eg.db::org.Gg.eg.db,
            keys = missing_entrez,
            keytype = kt,
            columns = c("ENTREZID")
          ),
          error = function(e) NULL
        )
        if (!is.null(mp) && nrow(mp) > 0) {
          names(mp)[names(mp) == kt] <- "gene_id"
          mp <- mp %>%
            dplyr::filter(!is.na(ENTREZID), ENTREZID != "") %>%
            dplyr::group_by(gene_id) %>%
            dplyr::summarise(mapped_entrez = get_first_nonempty(ENTREZID), .groups = "drop")
          out <- out %>%
            dplyr::left_join(mp, by = "gene_id") %>%
            dplyr::mutate(entrez_enrich = dplyr::coalesce(entrez_enrich, mapped_entrez)) %>%
            dplyr::select(-mapped_entrez)
        }
        missing_entrez <- out %>% dplyr::filter(is.na(entrez_enrich)) %>% dplyr::pull(gene_id) %>% unique()
        if (length(missing_entrez) == 0) break
      }
    }

    # 2) Try symbol mapping for still-unmapped genes.
    if ("SYMBOL" %in% keytypes_available && "ENTREZID" %in% columns_available) {
      missing_symbol <- out %>%
        dplyr::filter(is.na(entrez_enrich), !is.na(symbol), symbol != "") %>%
        dplyr::pull(symbol) %>% unique()
      if (length(missing_symbol) > 0) {
        smp <- tryCatch(
          AnnotationDbi::select(
            org.Gg.eg.db::org.Gg.eg.db,
            keys = missing_symbol,
            keytype = "SYMBOL",
            columns = c("ENTREZID")
          ),
          error = function(e) NULL
        )
        if (!is.null(smp) && nrow(smp) > 0) {
          smp <- smp %>%
            dplyr::filter(!is.na(ENTREZID), ENTREZID != "") %>%
            dplyr::group_by(SYMBOL) %>%
            dplyr::summarise(mapped_entrez_symbol = get_first_nonempty(ENTREZID), .groups = "drop")
          out <- out %>%
            dplyr::left_join(smp, by = c("symbol" = "SYMBOL")) %>%
            dplyr::mutate(entrez_enrich = dplyr::coalesce(entrez_enrich, mapped_entrez_symbol)) %>%
            dplyr::select(-mapped_entrez_symbol)
        }
      }
    }
  }

  out %>%
    dplyr::mutate(
      entrez_enrich = stringr::str_extract(as.character(entrez_enrich), "[0-9]+"),
      entrez_enrich = dplyr::if_else(is.na(entrez_enrich) | entrez_enrich == "", NA_character_, entrez_enrich)
    )
}

plot_ora_dot <- function(ora_table, group_col, title, stem, max_terms = ENRICH_SHOW_CATEGORY) {
  if (is.null(ora_table) || nrow(ora_table) == 0) return(invisible(NULL))
  if (!all(c(group_col, "Description", "p.adjust") %in% names(ora_table))) return(invisible(NULL))

  plot_df <- ora_table %>%
    dplyr::mutate(
      Group = as.character(.data[[group_col]]),
      GeneRatio_num = parse_gene_ratio(GeneRatio),
      Description_short = shorten_text(Description, 72)
    ) %>%
    dplyr::filter(!is.na(GeneRatio_num), !is.na(p.adjust)) %>%
    dplyr::group_by(Group) %>%
    dplyr::arrange(p.adjust, dplyr::desc(GeneRatio_num), .by_group = TRUE) %>%
    dplyr::slice_head(n = max_terms) %>%
    dplyr::ungroup()

  if (nrow(plot_df) == 0) return(invisible(NULL))

  p <- ggplot(plot_df, aes(x = Group, y = reorder(Description_short, GeneRatio_num))) +
    geom_point(aes(size = GeneRatio_num, color = p.adjust), alpha = 0.9) +
    scale_color_gradient(low = "#D73027", high = "#4575B4", trans = "reverse", name = "FDR") +
    scale_size_continuous(name = "Gene ratio") +
    labs(title = title, x = NULL, y = NULL) +
    theme_rna(9) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.text.y = element_text(size = 7))

  save_plot(p, stem, width = 9.8, height = max(4.8, 0.23 * dplyr::n_distinct(plot_df$Description_short) + 2.5))
  invisible(plot_df)
}

run_compare_kegg <- function(input_df, group_col, universe_entrez, out_prefix, title) {
  if (!HAS_CLUSTERPROFILER) return(tibble())
  suppressPackageStartupMessages(library(clusterProfiler))

  if (nrow(input_df) == 0) return(tibble())
  input_df <- input_df %>%
    dplyr::filter(!is.na(entrez_enrich), entrez_enrich != "") %>%
    dplyr::distinct(.data[[group_col]], entrez_enrich)

  group_sizes <- input_df %>% dplyr::count(.data[[group_col]], name = "N_entrez")
  valid_groups <- group_sizes %>% dplyr::filter(N_entrez >= ENRICH_MIN_SET_SIZE) %>% dplyr::pull(.data[[group_col]])
  input_df <- input_df %>% dplyr::filter(.data[[group_col]] %in% valid_groups)
  if (nrow(input_df) == 0) return(tibble())

  names(input_df)[names(input_df) == group_col] <- "Group"
  res <- tryCatch(
    clusterProfiler::compareCluster(
      entrez_enrich ~ Group,
      data = input_df,
      fun = "enrichKEGG",
      organism = "gga",
      universe = universe_entrez,
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1
    ),
    error = function(e) {
      message(out_prefix, " KEGG enrichment failed: ", e$message)
      NULL
    }
  )
  if (is.null(res)) return(tibble())

  df <- as.data.frame(res)
  if (nrow(df) == 0) return(tibble())
  df <- df %>%
    dplyr::mutate(
      Analysis = out_prefix,
      Reported = !is.na(p.adjust) & p.adjust <= ENRICH_REPORT_FDR
    ) %>%
    dplyr::arrange(p.adjust, pvalue)
  write_csv_safe(df, paste0("Enrichment_", out_prefix, "_KEGG_all.csv"))
  write_csv_safe(df %>% dplyr::filter(Reported), paste0("Enrichment_", out_prefix, "_KEGG_reported_FDR", ENRICH_REPORT_FDR, ".csv"))
  plot_ora_dot(df %>% dplyr::filter(Reported), "Cluster", title, paste0("Fig_Enrichment_", out_prefix, "_KEGG_dotplot"))
  df
}

run_compare_go <- function(input_df, group_col, universe_entrez, out_prefix, title, ont = "BP") {
  if (!HAS_CLUSTERPROFILER || !HAS_ORG_GG) return(tibble())
  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Gg.eg.db)
  })

  if (nrow(input_df) == 0) return(tibble())
  input_df <- input_df %>%
    dplyr::filter(!is.na(entrez_enrich), entrez_enrich != "") %>%
    dplyr::distinct(.data[[group_col]], entrez_enrich)

  group_sizes <- input_df %>% dplyr::count(.data[[group_col]], name = "N_entrez")
  valid_groups <- group_sizes %>% dplyr::filter(N_entrez >= ENRICH_MIN_SET_SIZE) %>% dplyr::pull(.data[[group_col]])
  input_df <- input_df %>% dplyr::filter(.data[[group_col]] %in% valid_groups)
  if (nrow(input_df) == 0) return(tibble())

  names(input_df)[names(input_df) == group_col] <- "Group"
  res <- tryCatch(
    clusterProfiler::compareCluster(
      entrez_enrich ~ Group,
      data = input_df,
      fun = "enrichGO",
      OrgDb = org.Gg.eg.db::org.Gg.eg.db,
      keyType = "ENTREZID",
      ont = ont,
      universe = universe_entrez,
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      readable = TRUE
    ),
    error = function(e) {
      message(out_prefix, " GO-", ont, " enrichment failed: ", e$message)
      NULL
    }
  )
  if (is.null(res)) return(tibble())

  df <- as.data.frame(res)
  if (nrow(df) == 0) return(tibble())
  df <- df %>%
    dplyr::mutate(
      Analysis = out_prefix,
      Ontology = ont,
      Reported = !is.na(p.adjust) & p.adjust <= ENRICH_REPORT_FDR
    ) %>%
    dplyr::arrange(p.adjust, pvalue)
  write_csv_safe(df, paste0("Enrichment_", out_prefix, "_GO_", ont, "_all.csv"))
  write_csv_safe(df %>% dplyr::filter(Reported), paste0("Enrichment_", out_prefix, "_GO_", ont, "_reported_FDR", ENRICH_REPORT_FDR, ".csv"))
  plot_ora_dot(df %>% dplyr::filter(Reported), "Cluster", paste0(title, " GO ", ont), paste0("Fig_Enrichment_", out_prefix, "_GO_", ont, "_dotplot"))
  df
}

run_module_fisher <- function(pattern_table, category_col, category_levels, modules, out_prefix) {
  universe_df <- pattern_table %>%
    dplyr::mutate(
      category = as.character(.data[[category_col]]),
      gene_symbol_upper = toupper(dplyr::coalesce(symbol, gene_label, gene_id))
    ) %>%
    dplyr::filter(category %in% category_levels) %>%
    dplyr::distinct(gene_id, category, gene_symbol_upper)

  universe_genes <- unique(universe_df$gene_id)
  if (length(universe_genes) == 0) return(tibble())

  res <- list()
  k <- 1
  for (cat in category_levels) {
    cat_genes <- universe_df$gene_id[universe_df$category == cat]
    for (module_name in names(modules)) {
      module_symbols <- toupper(unique(modules[[module_name]]))
      module_genes <- universe_df$gene_id[universe_df$gene_symbol_upper %in% module_symbols]
      overlap <- intersect(cat_genes, module_genes)

      a <- length(overlap)
      b <- length(setdiff(cat_genes, module_genes))
      c <- length(setdiff(module_genes, cat_genes))
      d <- length(setdiff(universe_genes, union(cat_genes, module_genes)))

      ft <- tryCatch(stats::fisher.test(matrix(c(a, b, c, d), nrow = 2)), error = function(e) NULL)
      overlap_labels <- universe_df %>%
        dplyr::filter(gene_id %in% overlap) %>%
        dplyr::mutate(label = dplyr::if_else(!is.na(gene_symbol_upper) & gene_symbol_upper != "", gene_symbol_upper, gene_id)) %>%
        dplyr::pull(label) %>% unique()

      res[[k]] <- tibble(
        category = cat,
        module = module_name,
        N_category_genes = length(unique(cat_genes)),
        N_module_genes_in_universe = length(unique(module_genes)),
        N_overlap = length(unique(overlap)),
        odds_ratio = if (!is.null(ft)) unname(ft$estimate) else NA_real_,
        pvalue = if (!is.null(ft)) ft$p.value else NA_real_,
        overlap_genes = paste(sort(overlap_labels), collapse = "; ")
      )
      k <- k + 1
    }
  }

  out <- dplyr::bind_rows(res) %>%
    dplyr::mutate(p.adjust = stats::p.adjust(pvalue, method = "BH")) %>%
    dplyr::arrange(p.adjust, pvalue, dplyr::desc(N_overlap))

  write_csv_safe(out, paste0("Enrichment_", out_prefix, "_curated_duodenum_modules.csv"))

  plot_df <- out %>%
    dplyr::filter(N_overlap > 0) %>%
    dplyr::mutate(module = gsub("_", " ", module))

  if (nrow(plot_df) > 0) {
    p <- ggplot(plot_df, aes(x = category, y = module)) +
      geom_point(aes(size = N_overlap, color = p.adjust), alpha = 0.9) +
      scale_color_gradient(low = "#D73027", high = "#4575B4", trans = "reverse", name = "FDR") +
      scale_size_continuous(name = "Overlap genes") +
      labs(title = paste0("Curated duodenum-module enrichment: ", out_prefix), x = NULL, y = NULL) +
      theme_rna(9) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
    save_plot(p, paste0("Fig_Enrichment_", out_prefix, "_curated_module_dotplot"), 8.8, 5.2)
  }

  out
}

if (RUN_PATTERN_ENRICHMENT) {
  message("\nRunning enrichment analyses for F1 expression-pattern gene sets...")

  pattern_df <- prepare_entrez_for_enrichment(pattern_df)
  write_csv_safe(
    pattern_df %>% dplyr::select(gene_id, symbol, gene_label, entrez, entrez_enrich, pattern, inheritance_class),
    "QC_enrichment_entrez_mapping_used.csv"
  )

  universe_entrez <- pattern_df %>%
    dplyr::filter(!is.na(entrez_enrich), entrez_enrich != "") %>%
    dplyr::pull(entrez_enrich) %>% unique()

  pattern_enrich_input <- pattern_df %>%
    dplyr::filter(pattern %in% pattern_levels_main) %>%
    dplyr::filter(!is.na(entrez_enrich), entrez_enrich != "") %>%
    dplyr::distinct(pattern, entrez_enrich)

  class_enrich_input <- pattern_df %>%
    dplyr::filter(inheritance_class %in% c("additive_like", "dominance_like", "overdominance_like")) %>%
    dplyr::filter(!is.na(entrez_enrich), entrez_enrich != "") %>%
    dplyr::distinct(inheritance_class, entrez_enrich)

  enrichment_input_summary <- dplyr::bind_rows(
    pattern_enrich_input %>% dplyr::count(pattern, name = "N_entrez") %>% dplyr::rename(Set = pattern) %>% dplyr::mutate(Analysis = "F1_pattern"),
    class_enrich_input %>% dplyr::count(inheritance_class, name = "N_entrez") %>% dplyr::rename(Set = inheritance_class) %>% dplyr::mutate(Analysis = "inheritance_class")
  ) %>% dplyr::select(Analysis, Set, N_entrez)
  write_csv_safe(enrichment_input_summary, "Enrichment_input_gene_set_sizes.csv")

  if (length(universe_entrez) >= ENRICH_MIN_SET_SIZE) {
    ora_pattern_kegg <- run_compare_kegg(pattern_enrich_input, "pattern", universe_entrez, "F1_pattern", "KEGG enrichment by F1 expression pattern")
    ora_pattern_go_bp <- run_compare_go(pattern_enrich_input, "pattern", universe_entrez, "F1_pattern", "GO enrichment by F1 expression pattern", ont = "BP")
    ora_pattern_go_cc <- run_compare_go(pattern_enrich_input, "pattern", universe_entrez, "F1_pattern", "GO enrichment by F1 expression pattern", ont = "CC")
    ora_pattern_go_mf <- run_compare_go(pattern_enrich_input, "pattern", universe_entrez, "F1_pattern", "GO enrichment by F1 expression pattern", ont = "MF")

    ora_class_kegg <- run_compare_kegg(class_enrich_input, "inheritance_class", universe_entrez, "inheritance_class", "KEGG enrichment by inheritance class")
    ora_class_go_bp <- run_compare_go(class_enrich_input, "inheritance_class", universe_entrez, "inheritance_class", "GO enrichment by inheritance class", ont = "BP")
    ora_class_go_cc <- run_compare_go(class_enrich_input, "inheritance_class", universe_entrez, "inheritance_class", "GO enrichment by inheritance class", ont = "CC")
    ora_class_go_mf <- run_compare_go(class_enrich_input, "inheritance_class", universe_entrez, "inheritance_class", "GO enrichment by inheritance class", ont = "MF")
  } else {
    message("Too few Entrez-mapped genes for GO/KEGG ORA. Curated module enrichment will still be attempted.")
  }

  duodenum_modules <- list(
    lipid_absorption_lipoprotein = c("APOA1", "APOA4", "APOB", "APOC3", "MTTP", "SAR1B", "MOGAT2", "DGAT1", "DGAT2", "CD36", "SCARB1", "NPC1L1", "FABP1", "FABP2", "PLIN2"),
    fatty_acid_transport_oxidation = c("SLC27A1", "SLC27A2", "SLC27A4", "ACSL1", "ACSL3", "ACSL5", "CPT1A", "CPT2", "ACOX1", "ACOX2", "HADHA", "HADHB", "PPARA", "PPARD", "PPARG"),
    bile_acid_sterol_transport = c("SLC10A2", "SLC51A", "SLC51B", "FABP6", "NR1H4", "ABCA1", "ABCG5", "ABCG8", "NPC1L1"),
    intestinal_barrier_mucin = c("TJP1", "OCLN", "CLDN1", "CLDN2", "CLDN3", "CLDN4", "CDH1", "MUC2", "VIL1"),
    carbohydrate_digestion_absorption = c("MGAM", "SI", "SLC5A1", "SLC2A2", "SLC2A5", "AMY2A", "LCT"),
    immune_inflammatory_response = c("IL1R1", "IL1R2", "IL6", "TNF", "TLR2", "TLR4", "MYD88", "NFKB1", "CXCL8", "CCL20")
  )

  module_pattern_enrichment <- run_module_fisher(pattern_df, "pattern", pattern_levels_main, duodenum_modules, "F1_pattern")
  module_class_enrichment <- run_module_fisher(
    pattern_df,
    "inheritance_class",
    c("additive_like", "dominance_like", "overdominance_like"),
    duodenum_modules,
    "inheritance_class"
  )
} else {
  message("RNA_F1_RUN_PATTERN_ENRICHMENT=FALSE; enrichment analyses skipped.")
}

# ==============================================================================
# 15. Configuration and optional Excel workbook
# ==============================================================================

config_table <- tibble(
  parameter = c(
    "CORE_OUTDIR", "META_FILE", "OUTDIR", "TISSUE_LABEL", "GENE_SET_CONTEXT",
    "SIG_COLUMN", "PATTERN_ALPHA", "PAIRWISE_LFC_CUTOFF", "POSITION_TOLERANCE",
    "LABEL_TOP_N_DEG", "LABEL_TOP_N_DEVIATION",
    "LABEL_CANDIDATES", "MAX_TOTAL_LABELS", "TOP_N_PER_PATTERN",
    "RNA_F1_RUN_PATTERN_ENRICHMENT", "RNA_F1_ENRICH_MIN_SET_SIZE",
    "RNA_F1_ENRICH_REPORT_FDR",
    "N_YY", "N_YX", "N_XX", "N_total_genes", "N_classifiable_genes"
  ),
  value = as.character(c(
    CORE_OUTDIR, META_FILE, OUTDIR, TISSUE_LABEL, GENE_SET_CONTEXT,
    SIG_COLUMN, PATTERN_ALPHA, PAIRWISE_LFC_CUTOFF, POSITION_TOLERANCE,
    LABEL_TOP_N_DEG, LABEL_TOP_N_DEVIATION,
    LABEL_CANDIDATES, MAX_TOTAL_LABELS, TOP_N_PER_PATTERN,
    RUN_PATTERN_ENRICHMENT, ENRICH_MIN_SET_SIZE, ENRICH_REPORT_FDR,
    group_n[["YY"]], group_n[["YX"]], group_n[["XX"]],
    n_total_genes, n_classified_total
  ))
)
write_csv_safe(config_table, "Analysis_configuration.csv")

if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(
    list(
      configuration = config_table,
      samples_used = bind_rows(sample_qc, upstream_absent_qc),
      pattern_summary_all = pattern_summary_all,
      pattern_summary_classified = pattern_summary_main,
      inheritance_summary = inheritance_summary,
      all_gene_patterns = pattern_df,
      top_pattern_genes = top_pattern_genes,
      nonadditive_like_genes = nonadditive_genes,
      labeled_deviation_genes = label_df,
      duodenum_candidates = candidate_summary,
      enrichment_input_summary = enrichment_input_summary,
      ORA_pattern_KEGG = ora_pattern_kegg,
      ORA_pattern_GO_BP = ora_pattern_go_bp,
      ORA_pattern_GO_CC = ora_pattern_go_cc,
      ORA_pattern_GO_MF = ora_pattern_go_mf,
      ORA_class_KEGG = ora_class_kegg,
      ORA_class_GO_BP = ora_class_go_bp,
      ORA_class_GO_CC = ora_class_go_cc,
      ORA_class_GO_MF = ora_class_go_mf,
      curated_module_pattern = module_pattern_enrichment,
      curated_module_class = module_class_enrichment
    ),
    path = file.path(OUTDIR, "duodenum_F1_expression_inheritance_results.xlsx")
  )
} else {
  message("Optional package 'writexl' is not installed; CSV outputs were created, Excel workbook skipped.")
}

# ==============================================================================
# 16. Console summary
# ==============================================================================

message("\n============================================================")
message("Duodenum F1 expression-inheritance analysis complete.")
message("Core input:  ", CORE_OUTDIR)
message("Output:      ", OUTDIR)
message("Significance: ", SIG_COLUMN, " < ", PATTERN_ALPHA,
        ", |log2FC| >= ", PAIRWISE_LFC_CUTOFF)
message("Samples: YY = ", group_n[["YY"]],
        ", YX = ", group_n[["YX"]],
        ", XX = ", group_n[["XX"]])
message("All genes: ", n_total_genes)
message("Classifiable genes: ", n_classified_total)
message("Main pattern percentages sum to: ",
        sprintf("%.2f%%", sum(pattern_summary_main$Percent_classified)))
message("\nKey tables:")
message("  QC_samples_used_and_upstream_absent.csv")
message("  Table_F1_expression_pattern_all_genes.csv")
message("  Table_F1_expression_pattern_summary_classified_only.csv")
message("  Table_F1_mid_parent_deviation_labeled_genes.csv")
message("  Table_F1_top_genes_by_expression_pattern.csv")
message("  Table_nonadditive_like_expression_genes.csv")
message("\nKey figures:")
message("  Fig_F1_expression_pattern_counts_for_manuscript.png/pdf")
message("  Fig_F1_expression_position_distribution.png/pdf")
message("  Fig_F1_mid_parent_deviation_labeled.png/pdf")
message("  Fig_F1_top_pattern_genes_heatmap_samples.png/pdf")
message("  Fig_F1_top_pattern_genes_heatmap_group_mean.png/pdf")
message("============================================================\n")

