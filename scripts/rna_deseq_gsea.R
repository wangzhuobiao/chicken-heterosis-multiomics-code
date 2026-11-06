###############################################################################
# Duodenum transcriptome core analysis for YY, YX F1 and XX chickens
# v1.0: DESeq2 + DEG + GSEA + duodenum lipid-absorption candidate modules
#
# This is the CORE script. It intentionally removes F1 inheritance/heterosis
# pattern classification. Run the standalone script below after this script:
#   05C_transcriptome_F1_inheritance_pattern_standalone.R
###############################################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(20260623)

# =============================================================================
# 0. User parameters
# =============================================================================

TISSUE_LABEL <- Sys.getenv("RNA_TISSUE_LABEL", unset = "duodenum")
TISSUE_TITLE <- Sys.getenv("RNA_TISSUE_TITLE", unset = "Duodenum transcriptome")

COUNT_FILE <- Sys.getenv(
  "COUNT_FILE",
  unset = ""
)
META_FILE <- Sys.getenv(
  "META_FILE",
  unset = ""
)
ANNOT_FILE <- Sys.getenv(
  "ANNOT_FILE",
  unset = ""
)
OUTDIR <- Sys.getenv(
  "RNA_OUTDIR",
  unset = ""
)
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

GROUP_LEVELS <- c("YY", "YX", "XX")
GROUP_COLORS <- c("YY" = "#E64B35", "YX" = "#4DBBD5", "XX" = "#00A087")
DESIGN_FORMULA_TEXT <- Sys.getenv("RNA_DESIGN", unset = "~ group")

# DEG filtering modes:
# strict / conventional / balanced / exploratory / fdr_only / pvalue / custom
DEG_FILTER_MODE <- tolower(Sys.getenv("RNA_DEG_FILTER_MODE", unset = "pvalue"))
get_deg_preset <- function(mode) {
  presets <- list(
    strict       = list(alpha = 0.01, lfc = 1.00, min_count = 10L, min_samples = 3L, min_basemean = 10, sig_col = "padj"),
    conventional = list(alpha = 0.05, lfc = 1.00, min_count = 10L, min_samples = 3L, min_basemean = 10, sig_col = "padj"),
    balanced     = list(alpha = 0.05, lfc = 0.58, min_count = 5L,  min_samples = 3L, min_basemean = 5,  sig_col = "padj"),
    exploratory  = list(alpha = 0.10, lfc = 0.58, min_count = 3L,  min_samples = 2L, min_basemean = 3,  sig_col = "padj"),
    fdr_only     = list(alpha = 0.05, lfc = 0.00, min_count = 5L,  min_samples = 3L, min_basemean = 5,  sig_col = "padj"),
    pvalue       = list(alpha = 0.05, lfc = 0.58, min_count = 5L,  min_samples = 3L, min_basemean = 5,  sig_col = "pvalue"),
    custom       = list(alpha = 0.05, lfc = 0.58, min_count = 5L,  min_samples = 3L, min_basemean = 5,  sig_col = "padj")
  )
  if (!mode %in% names(presets)) mode <- "balanced"
  presets[[mode]]
}
DEG_PRESET <- get_deg_preset(DEG_FILTER_MODE)
ALPHA <- as.numeric(Sys.getenv("RNA_ALPHA", unset = as.character(DEG_PRESET$alpha)))
LFC_CUTOFF <- as.numeric(Sys.getenv("RNA_LFC_CUTOFF", unset = as.character(DEG_PRESET$lfc)))
MIN_COUNT <- as.integer(Sys.getenv("RNA_MIN_COUNT", unset = as.character(DEG_PRESET$min_count)))
MIN_SAMPLES <- as.integer(Sys.getenv("RNA_MIN_SAMPLES", unset = as.character(DEG_PRESET$min_samples)))
MIN_BASEMEAN_FOR_DEG_CALL <- as.numeric(Sys.getenv("RNA_MIN_BASEMEAN_FOR_DEG_CALL", unset = as.character(DEG_PRESET$min_basemean)))
SIGNIFICANCE_COLUMN <- Sys.getenv("RNA_SIGNIFICANCE_COLUMN", unset = DEG_PRESET$sig_col)
if (!SIGNIFICANCE_COLUMN %in% c("padj", "pvalue")) SIGNIFICANCE_COLUMN <- "padj"
SIGNIFICANCE_LABEL <- ifelse(SIGNIFICANCE_COLUMN == "pvalue", "nominal P value", "FDR")
SIGNIFICANCE_AXIS_LABEL <- ifelse(SIGNIFICANCE_COLUMN == "pvalue", "-log10 nominal P value", "-log10 adjusted P value")

DEG_STANDARD_GRID <- data.frame(
  standard = c("strict", "conventional", "balanced", "exploratory", "fdr_only", "pvalue"),
  sig_col = c("padj", "padj", "padj", "padj", "padj", "pvalue"),
  alpha = c(0.01, 0.05, 0.05, 0.10, 0.05, 0.05),
  lfc = c(1.00, 1.00, 0.58, 0.58, 0.00, 0.58),
  min_baseMean = c(10, 10, 5, 3, 5, 5),
  stringsAsFactors = FALSE
)

VOLCANO_Y_CAP <- as.numeric(Sys.getenv("RNA_VOLCANO_Y_CAP", unset = "50"))
VOLCANO_X_CAP <- as.numeric(Sys.getenv("RNA_VOLCANO_X_CAP", unset = "8"))
RUN_DEG_ENRICHMENT <- Sys.getenv("RNA_RUN_DEG_ENRICHMENT", unset = "TRUE") == "TRUE"
RUN_GSEA <- Sys.getenv("RNA_RUN_GSEA", unset = "TRUE") == "TRUE"
GSEA_RANK_METRIC <- tolower(Sys.getenv("RNA_GSEA_RANK_METRIC", unset = "stat"))
GSEA_PADJ_REPORT_CUTOFF <- as.numeric(Sys.getenv("RNA_GSEA_PADJ_REPORT_CUTOFF", unset = "0.25"))
GSEA_MIN_GSSIZE <- as.integer(Sys.getenv("RNA_GSEA_MIN_GSSIZE", unset = "10"))
GSEA_MAX_GSSIZE <- as.integer(Sys.getenv("RNA_GSEA_MAX_GSSIZE", unset = "500"))
ENRICH_MIN_GENESET_SIZE <- as.integer(Sys.getenv("RNA_ENRICH_MIN_GENESET_SIZE", unset = "3"))
ENRICH_SHOW_CATEGORY <- as.integer(Sys.getenv("RNA_ENRICH_SHOW_CATEGORY", unset = "8"))
TOP_N_HEATMAP_PER_SET <- as.integer(Sys.getenv("RNA_DUODENUM_HEATMAP_TOP_N_PER_SET", unset = "5"))
MAX_HEATMAP_GENES <- as.integer(Sys.getenv("RNA_DUODENUM_HEATMAP_MAX_GENES", unset = "50"))
USE_ORGDB_ANNOT <- Sys.getenv("RNA_USE_ORGDB_ANNOT", unset = "TRUE") == "TRUE"
EXCLUDE_SAMPLES_RAW <- Sys.getenv("DUODENUM_EXCLUDE_SAMPLES", unset = "")
parse_exclude <- function(x) {
  if (is.na(x) || trimws(x) == "") return(character(0))
  x <- gsub("[,;]+", " ", x)
  z <- unlist(strsplit(x, "[[:space:]]+"))
  trimws(z[z != ""])
}
EXCLUDE_SAMPLES <- parse_exclude(EXCLUDE_SAMPLES_RAW)

# =============================================================================
# 1. Packages
# =============================================================================

cran_pkgs <- c("tidyverse", "data.table", "readxl", "writexl", "ggpubr", "pheatmap", "ggrepel", "patchwork")
bioc_pkgs <- c("DESeq2", "clusterProfiler", "AnnotationDbi", "org.Gg.eg.db", "enrichplot")
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
if (!requireNamespace("BiocManager", quietly = TRUE)) stop(paste("Install required R package:", "BiocManager"))
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) stop(paste("Install required R package:", p))
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}
for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    tryCatch(stop(paste("Install required Bioconductor package:", p)), error = function(e) message("Install failed: ", p, ": ", e$message))
  }
}
suppressPackageStartupMessages(library(DESeq2))
HAS_CLUSTERPROFILER <- requireNamespace("clusterProfiler", quietly = TRUE)
HAS_ORG_GG <- requireNamespace("org.Gg.eg.db", quietly = TRUE)
HAS_ENRICHPLOT <- requireNamespace("enrichplot", quietly = TRUE)

# =============================================================================
# 2. Helper functions
# =============================================================================

`%||%` <- function(a, b) if (length(a) == 0 || is.null(a) || all(is.na(a))) b else a[1]

safe_read_table <- function(file) {
  if (!file.exists(file)) stop("Input file does not exist: ", file)
  ext <- tolower(tools::file_ext(file))
  if (ext %in% c("xlsx", "xls")) as.data.frame(readxl::read_excel(file, .name_repair = "unique"))
  else data.table::fread(file, data.table = FALSE)
}

safe_numeric <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, ",|，", "")
  x <- stringr::str_trim(x)
  x[x %in% c("", "NA", "NaN", "null", "NULL", "-", "--")] <- NA_character_
  suppressWarnings(as.numeric(x))
}

clean_gene_id <- function(x) stringr::str_replace(stringr::str_trim(as.character(x)), "\\.\\d+$", "")

normalize_group <- function(x, sample_id = NULL) {
  normalize_one <- function(z) {
    if (is.na(z) || stringr::str_trim(z) == "") return(NA_character_)
    u <- stringr::str_to_upper(z)
    u <- stringr::str_replace_all(u, "×|\\u00d7", "X")
    u <- stringr::str_replace_all(u, "\\s+|[-_\\.]", "")
    if (stringr::str_detect(u, "F1|杂交|YX")) return("YX")
    if (stringr::str_detect(u, "YY|^Y$|Y组|Y系|YLINE")) return("YY")
    if (stringr::str_detect(u, "XX|^X$|X组|X系|XLINE")) return("XX")
    NA_character_
  }
  out <- vapply(as.character(x), normalize_one, character(1))
  if (!is.null(sample_id)) {
    sid <- vapply(as.character(sample_id), normalize_one, character(1))
    out[is.na(out) & !is.na(sid)] <- sid[is.na(out) & !is.na(sid)]
  }
  factor(out, levels = GROUP_LEVELS)
}

theme_rna <- function(base_size = 11) {
  ggpubr::theme_pubr(base_size = base_size) +
    theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(color = "grey35"),
          axis.title = element_text(face = "bold"), axis.text = element_text(color = "black"),
          strip.background = element_rect(fill = "grey95", color = "grey70"), strip.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold"))
}

save_plot <- function(p, name, width, height, dpi = 300) {
  ggsave(file.path(OUTDIR, paste0(name, ".png")), p, width = width, height = height, dpi = dpi)
  ggsave(file.path(OUTDIR, paste0(name, ".pdf")), p, width = width, height = height)
}

label_compact_count <- function(x) {
  x <- as.numeric(x)
  ifelse(abs(x) >= 1e9, paste0(round(x / 1e9, 1), "B"),
         ifelse(abs(x) >= 1e6, paste0(round(x / 1e6, 1), "M"),
                ifelse(abs(x) >= 1e3, paste0(round(x / 1e3, 1), "K"), format(round(x), big.mark = ",", scientific = FALSE))))
}

sig_to_deg <- function(sig_value, lfc, baseMean = NULL, alpha = ALPHA, lfc_cutoff = LFC_CUTOFF, min_baseMean = MIN_BASEMEAN_FOR_DEG_CALL) {
  ok <- !is.na(sig_value) & sig_value < alpha & !is.na(lfc) & abs(lfc) >= lfc_cutoff
  if (!is.null(baseMean)) ok <- ok & !is.na(baseMean) & baseMean >= min_baseMean
  ok
}

sig_by_standard <- function(sig_value, lfc, baseMean, alpha, lfc_cutoff, min_baseMean) {
  !is.na(sig_value) & !is.na(lfc) & !is.na(baseMean) & baseMean >= min_baseMean & sig_value < alpha & abs(lfc) >= lfc_cutoff
}

parse_gene_ratio <- function(x) {
  x <- as.character(x); out <- rep(NA_real_, length(x)); has <- grepl("/", x)
  if (any(has)) out[has] <- vapply(strsplit(x[has], "/"), function(z) if (length(z) == 2) as.numeric(z[1]) / as.numeric(z[2]) else NA_real_, numeric(1))
  suppressWarnings(out[!has] <- as.numeric(x[!has])); out
}

shorten_term <- function(x, max_len = 72) ifelse(nchar(as.character(x)) > max_len, paste0(substr(as.character(x), 1, max_len - 3), "..."), as.character(x))

extract_gtf_attribute <- function(attr, key) {
  attr <- as.character(attr)
  pat1 <- paste0("(?:^|;|\\s)", key, "\\s+\\\"([^\\\"]+)\\\"")
  out <- stringr::str_match(attr, pat1)[, 2]
  miss <- is.na(out)
  if (any(miss)) {
    pat2 <- paste0("(?:^|;)", key, "=([^;]+)")
    out[miss] <- stringr::str_match(attr[miss], pat2)[, 2]
  }
  stringr::str_replace_all(out, "%20", " ")
}

parse_gtf_annotation <- function(file) {
  message("Parsing annotation: ", file)
  gtf <- data.table::fread(file, sep = "\t", header = FALSE, quote = "", comment.char = "#", data.table = FALSE,
                           col.names = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"))
  if (any(gtf$feature == "gene")) gtf <- gtf %>% dplyr::filter(feature == "gene") else gtf <- gtf %>% dplyr::filter(feature %in% c("transcript", "mRNA", "exon", "CDS"))
  gene_id <- dplyr::coalesce(extract_gtf_attribute(gtf$attribute, "gene_id"), extract_gtf_attribute(gtf$attribute, "gene"), extract_gtf_attribute(gtf$attribute, "ID"))
  gene_id <- stringr::str_replace(gene_id, "^gene:", "")
  tibble::tibble(
    gene_id_clean = clean_gene_id(gene_id),
    symbol_file = dplyr::coalesce(extract_gtf_attribute(gtf$attribute, "gene_name"), extract_gtf_attribute(gtf$attribute, "gene_symbol"), extract_gtf_attribute(gtf$attribute, "external_gene_name"), extract_gtf_attribute(gtf$attribute, "Name")),
    gene_name_file = dplyr::coalesce(extract_gtf_attribute(gtf$attribute, "description"), extract_gtf_attribute(gtf$attribute, "gene_description"), extract_gtf_attribute(gtf$attribute, "product")),
    gene_biotype_file = dplyr::coalesce(extract_gtf_attribute(gtf$attribute, "gene_biotype"), extract_gtf_attribute(gtf$attribute, "gene_type"), extract_gtf_attribute(gtf$attribute, "biotype")),
    entrez_file = dplyr::coalesce(extract_gtf_attribute(gtf$attribute, "entrez_id"), extract_gtf_attribute(gtf$attribute, "EntrezGene"), extract_gtf_attribute(gtf$attribute, "Dbxref"), extract_gtf_attribute(gtf$attribute, "db_xref"))
  ) %>%
    dplyr::filter(!is.na(gene_id_clean), gene_id_clean != "") %>%
    dplyr::mutate(entrez_file = stringr::str_extract(as.character(entrez_file), "\\d+"), dplyr::across(dplyr::everything(), ~ dplyr::na_if(.x, ""))) %>%
    dplyr::distinct(gene_id_clean, .keep_all = TRUE)
}

standardize_annotation_table <- function(x) {
  x <- as.data.frame(x); names(x) <- stringr::str_trim(names(x))
  id_col <- intersect(c("gene_id", "Geneid", "gene", "Gene stable ID", "ensembl_gene_id", "ENSEMBL", "ID", "gene_id_clean"), names(x))[1] %||% names(x)[1]
  sym_col <- intersect(c("symbol", "gene_symbol", "external_gene_name", "Gene name", "SYMBOL"), names(x))[1]
  name_col <- intersect(c("gene_name", "description", "Gene description", "GENENAME", "Description"), names(x))[1]
  biotype_col <- intersect(c("gene_biotype", "gene_type", "Gene type", "biotype"), names(x))[1]
  entrez_col <- intersect(c("entrez", "entrez_id", "EntrezGene ID", "NCBI gene ID", "ENTREZID"), names(x))[1]
  tibble::tibble(
    gene_id_clean = clean_gene_id(x[[id_col]]),
    symbol_file = if (!is.na(sym_col)) as.character(x[[sym_col]]) else NA_character_,
    gene_name_file = if (!is.na(name_col)) as.character(x[[name_col]]) else NA_character_,
    gene_biotype_file = if (!is.na(biotype_col)) as.character(x[[biotype_col]]) else NA_character_,
    entrez_file = if (!is.na(entrez_col)) as.character(x[[entrez_col]]) else NA_character_
  ) %>% dplyr::filter(!is.na(gene_id_clean), gene_id_clean != "") %>% dplyr::mutate(dplyr::across(dplyr::everything(), ~ dplyr::na_if(.x, ""))) %>% dplyr::distinct(gene_id_clean, .keep_all = TRUE)
}

annotate_from_orgdb <- function(gene_ids) {
  base <- tibble::tibble(gene_id_clean = unique(clean_gene_id(gene_ids)))
  empty <- base %>% dplyr::mutate(symbol_orgdb = NA_character_, gene_name_orgdb = NA_character_, entrez_orgdb = NA_character_)
  if (!USE_ORGDB_ANNOT || !requireNamespace("org.Gg.eg.db", quietly = TRUE) || !requireNamespace("AnnotationDbi", quietly = TRUE)) return(empty)
  keytypes_avail <- tryCatch(AnnotationDbi::keytypes(org.Gg.eg.db::org.Gg.eg.db), error = function(e) character(0))
  columns_avail <- tryCatch(AnnotationDbi::columns(org.Gg.eg.db::org.Gg.eg.db), error = function(e) character(0))
  try_keytypes <- intersect(c("ENSEMBL", "ENSEMBLGENE", "ALIAS", "SYMBOL"), keytypes_avail)
  cols <- intersect(c("SYMBOL", "GENENAME", "ENTREZID"), columns_avail)
  if (length(try_keytypes) == 0 || length(cols) == 0) return(empty)
  best <- NULL; best_rate <- -1
  for (kt in try_keytypes) {
    tmp <- tryCatch(AnnotationDbi::select(org.Gg.eg.db::org.Gg.eg.db, keys = unique(base$gene_id_clean), keytype = kt, columns = cols), error = function(e) NULL)
    if (is.null(tmp) || nrow(tmp) == 0) next
    tmp2 <- tmp %>% as.data.frame() %>% dplyr::rename(gene_id_clean = dplyr::all_of(kt)) %>% dplyr::mutate(gene_id_clean = clean_gene_id(gene_id_clean)) %>%
      dplyr::transmute(gene_id_clean, symbol_orgdb = if ("SYMBOL" %in% names(.)) as.character(SYMBOL) else NA_character_, gene_name_orgdb = if ("GENENAME" %in% names(.)) as.character(GENENAME) else NA_character_, entrez_orgdb = if ("ENTREZID" %in% names(.)) as.character(ENTREZID) else NA_character_) %>%
      dplyr::group_by(gene_id_clean) %>% dplyr::summarise(symbol_orgdb = dplyr::first(stats::na.omit(symbol_orgdb)) %||% NA_character_, gene_name_orgdb = dplyr::first(stats::na.omit(gene_name_orgdb)) %||% NA_character_, entrez_orgdb = dplyr::first(stats::na.omit(entrez_orgdb)) %||% NA_character_, .groups = "drop")
    joined <- base %>% dplyr::left_join(tmp2, by = "gene_id_clean")
    rate <- mean(!is.na(joined$symbol_orgdb) & joined$symbol_orgdb != "")
    if (rate > best_rate) { best_rate <- rate; best <- joined }
  }
  if (is.null(best)) empty else best
}

build_gene_annotation <- function(gene_ids, annot_file = "") {
  input <- tibble::tibble(gene_id = gene_ids, gene_id_clean = clean_gene_id(gene_ids)) %>% dplyr::distinct(gene_id, .keep_all = TRUE)
  file_annot <- if (!is.na(annot_file) && annot_file != "" && file.exists(annot_file)) {
    ext <- tolower(tools::file_ext(annot_file))
    if (ext %in% c("gtf", "gff", "gff3")) parse_gtf_annotation(annot_file) else standardize_annotation_table(safe_read_table(annot_file))
  } else tibble::tibble(gene_id_clean = character(), symbol_file = character(), gene_name_file = character(), gene_biotype_file = character(), entrez_file = character())
  org_annot <- annotate_from_orgdb(gene_ids)
  # Safety: make sure OrgDb annotation always has these columns.
  # This prevents downstream coalesce() errors when OrgDb returns no usable mapping
  # or when a helper returns only gene_id_clean.
  for (cc in c("symbol_orgdb", "gene_name_orgdb", "entrez_orgdb")) {
    if (!cc %in% names(org_annot)) org_annot[[cc]] <- NA_character_
  }
  org_annot <- org_annot %>% dplyr::select(gene_id_clean, symbol_orgdb, gene_name_orgdb, entrez_orgdb) %>% dplyr::distinct(gene_id_clean, .keep_all = TRUE)
  annot <- input %>%
    dplyr::left_join(file_annot, by = "gene_id_clean") %>%
    dplyr::left_join(org_annot, by = "gene_id_clean") %>%
    dplyr::mutate(
      symbol = dplyr::coalesce(symbol_file, symbol_orgdb, gene_id),
      gene_name = dplyr::coalesce(gene_name_file, gene_name_orgdb),
      gene_biotype = gene_biotype_file,
      entrez = dplyr::coalesce(entrez_file, entrez_orgdb),
      annotation_source = dplyr::case_when(!is.na(symbol_file) & symbol_file != "" ~ "ANNOT_FILE", !is.na(symbol_orgdb) & symbol_orgdb != "" ~ "org.Gg.eg.db", TRUE ~ "gene_id_only"),
      display_gene = dplyr::if_else(!is.na(symbol) & symbol != "" & symbol != gene_id, symbol, gene_id)
    ) %>% dplyr::select(gene_id, gene_id_clean, symbol, display_gene, gene_name, gene_biotype, entrez, annotation_source, dplyr::everything())
  message(sprintf("Final gene-symbol annotation rate: %.2f%%", mean(!is.na(annot$symbol) & annot$symbol != "" & annot$symbol != annot$gene_id) * 100))
  annot
}

# =============================================================================
# 3. Duodenum candidate gene sets
# =============================================================================

DUODENUM_GENE_SETS <- list(
  intestinal_lipid_absorption_transport = c("CD36", "FABP1", "FABP2", "FABP6", "SLC27A1", "SLC27A2", "SLC27A4", "SLC27A5", "APOA1", "APOA4", "APOB", "MTTP", "NPC1L1", "SCARB1", "LPL", "ABCA1", "ABCG5", "ABCG8", "SOAT2", "DGAT1", "DGAT2", "MGAT2", "MOGAT1", "MOGAT2"),
  bile_acid_transport_metabolism = c("SLC10A2", "SLC51A", "SLC51B", "ABCB11", "ABCC2", "ABCC3", "ABCG5", "ABCG8", "CYP7A1", "CYP8B1", "CYP27A1", "CYP7B1", "NR1H4", "GPBAR1", "SLC27A5"),
  ppar_fatty_acid_metabolism = c("PPARA", "PPARD", "PPARG", "RXRA", "CPT1A", "CPT1B", "CPT2", "ACOX1", "ACADL", "ACADM", "ACADS", "ACSL1", "ACSL5", "EHHADH", "HADHA", "HADHB"),
  abc_transporters_xenobiotic = c("ABCA1", "ABCA2", "ABCA3", "ABCA5", "ABCB1", "ABCB4", "ABCB11", "ABCC2", "ABCC3", "ABCC4", "ABCC5", "ABCG2", "ABCG5", "ABCG8", "AHR", "CYP1A1", "CYP3A4", "GSTA1"),
  intestinal_barrier_tight_junction = c("TJP1", "TJP2", "TJP3", "OCLN", "CLDN1", "CLDN2", "CLDN3", "CLDN4", "CLDN5", "MUC2", "MUC13", "VIL1", "CDH1", "EPCAM", "ITGB1", "ACTB"),
  immune_inflammation = c("TLR2", "TLR4", "MYD88", "NFKB1", "NFKBIA", "IL1B", "IL6", "IL10", "TNF", "CXCL8", "CCL20", "PTGS2", "STAT1", "STAT3", "IFNG", "MUC2"),
  oxidative_stress_glutathione = c("GPX1", "GPX2", "GPX3", "GPX4", "GSR", "GSS", "GCLC", "GCLM", "SOD1", "SOD2", "CAT", "NQO1", "TXN", "PRDX1", "PRDX3"),
  carbohydrate_energy_absorption = c("SLC2A2", "SLC2A5", "SLC5A1", "SLC5A2", "HK1", "HK2", "PFKM", "PKM", "LDHA", "PDHA1", "CS", "IDH3A", "SDHA", "UQCRC1", "ATP5F1A")
)
write.csv(tibble::tibble(gene_set = rep(names(DUODENUM_GENE_SETS), lengths(DUODENUM_GENE_SETS)), symbol_query = unlist(DUODENUM_GENE_SETS, use.names = FALSE)), file.path(OUTDIR, "Table_Duodenum_candidate_gene_sets_defined.csv"), row.names = FALSE)

# =============================================================================
# 4. Read count matrix, metadata, annotation
# =============================================================================

message("Reading count matrix: ", COUNT_FILE)
count_df <- safe_read_table(COUNT_FILE)
if (!"gene_id" %in% names(count_df)) {
  id_candidates <- intersect(names(count_df), c("Geneid", "gene", "Gene stable ID", "ensembl_gene_id", "ENSEMBL"))
  if (length(id_candidates) > 0) names(count_df)[names(count_df) == id_candidates[1]] <- "gene_id" else names(count_df)[1] <- "gene_id"
}
count_df <- count_df %>% dplyr::mutate(gene_id = clean_gene_id(gene_id)) %>% dplyr::filter(!is.na(gene_id), gene_id != "", !stringr::str_starts(gene_id, "__"))
sample_cols <- setdiff(names(count_df), "gene_id")
count_df[sample_cols] <- lapply(count_df[sample_cols], safe_numeric)
if (any(duplicated(count_df$gene_id))) count_df <- count_df %>% dplyr::group_by(gene_id) %>% dplyr::summarise(dplyr::across(dplyr::all_of(sample_cols), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
count_mat <- as.matrix(count_df[, sample_cols, drop = FALSE]); mode(count_mat) <- "numeric"; rownames(count_mat) <- count_df$gene_id
if (any(is.na(count_mat))) stop("Count matrix contains NA values.")
if (any(count_mat < 0)) stop("Count matrix contains negative values.")
if (any(abs(count_mat - round(count_mat)) > 1e-6)) warning("Count matrix contains non-integers; rounded for DESeq2. Please confirm this is not TPM/FPKM.")
count_mat <- round(count_mat)

message("Reading metadata: ", META_FILE)
meta <- safe_read_table(META_FILE)
if (!"sample_id" %in% names(meta)) {
  candidate <- intersect(names(meta), c("Sample", "sample", "SampleID", "sample_name", "library_id", "RNA_id"))
  if (length(candidate) == 0) stop("Metadata must contain sample_id.")
  names(meta)[names(meta) == candidate[1]] <- "sample_id"
}
if (!"group" %in% names(meta)) {
  candidate <- intersect(names(meta), c("Group", "group_raw", "line", "Line", "population"))
  if (length(candidate) == 0) stop("Metadata must contain group.")
  names(meta)[names(meta) == candidate[1]] <- "group"
}
meta <- meta %>% dplyr::mutate(sample_id = as.character(sample_id), group_raw_before_normalization = as.character(group), group = normalize_group(group, sample_id = sample_id)) %>% dplyr::distinct(sample_id, .keep_all = TRUE)
write.csv(meta %>% dplyr::select(sample_id, group_raw_before_normalization, group), file.path(OUTDIR, "QC_metadata_group_normalization.csv"), row.names = FALSE)
invalid_group_meta <- meta %>% dplyr::filter(is.na(group))
if (nrow(invalid_group_meta) > 0) write.csv(invalid_group_meta, file.path(OUTDIR, "QC_metadata_invalid_group_labels.csv"), row.names = FALSE)
meta <- meta %>% dplyr::filter(!is.na(group))
if (length(EXCLUDE_SAMPLES) > 0) {
  meta <- meta %>% dplyr::filter(!sample_id %in% EXCLUDE_SAMPLES)
  message("Excluded samples: ", paste(EXCLUDE_SAMPLES, collapse = ", "))
}
shared_samples <- intersect(colnames(count_mat), meta$sample_id)
if (length(shared_samples) < 6) stop("Too few shared samples between count matrix and metadata: ", length(shared_samples))
write.csv(tibble::tibble(sample_id = setdiff(colnames(count_mat), meta$sample_id)), file.path(OUTDIR, "QC_samples_in_counts_not_metadata.csv"), row.names = FALSE)
write.csv(tibble::tibble(sample_id = setdiff(meta$sample_id, colnames(count_mat))), file.path(OUTDIR, "QC_samples_in_metadata_not_counts.csv"), row.names = FALSE)
count_mat <- count_mat[, shared_samples, drop = FALSE]
meta <- meta %>% dplyr::filter(sample_id %in% shared_samples) %>% dplyr::arrange(match(sample_id, shared_samples))
rownames(meta) <- meta$sample_id
stopifnot(all(colnames(count_mat) == rownames(meta)))
meta$group <- factor(meta$group, levels = GROUP_LEVELS)
sample_summary <- meta %>% as.data.frame() %>% dplyr::count(group, name = "N_samples") %>% tidyr::complete(group = factor(GROUP_LEVELS, levels = GROUP_LEVELS), fill = list(N_samples = 0))
write.csv(sample_summary, file.path(OUTDIR, "Table1_sample_counts_by_group.csv"), row.names = FALSE)
if (any(sample_summary$N_samples < 2)) stop("Each group needs at least 2 samples. Check Table1_sample_counts_by_group.csv.")

annot <- build_gene_annotation(rownames(count_mat), ANNOT_FILE)
write.csv(annot, file.path(OUTDIR, "Table2_gene_annotation_used.csv"), row.names = FALSE)
annotation_summary <- annot %>% dplyr::summarise(N_genes = dplyr::n(), N_symbol_mapped = sum(!is.na(symbol) & symbol != "" & symbol != gene_id), Symbol_mapping_rate_percent = N_symbol_mapped / N_genes * 100, N_entrez_mapped = sum(!is.na(entrez) & entrez != ""), Entrez_mapping_rate_percent = N_entrez_mapped / N_genes * 100, N_gene_id_only = sum(annotation_source == "gene_id_only"))
write.csv(annotation_summary, file.path(OUTDIR, "Table2A_gene_annotation_summary.csv"), row.names = FALSE)

# =============================================================================
# 5. DESeq2 normalization and QC
# =============================================================================

raw_lib_size <- colSums(count_mat)
lib_df <- tibble::tibble(sample_id = names(raw_lib_size), library_size = raw_lib_size) %>% dplyr::left_join(meta %>% dplyr::select(sample_id, group), by = "sample_id")
write.csv(lib_df, file.path(OUTDIR, "QC_library_size_raw_counts.csv"), row.names = FALSE)
p_lib <- ggplot(lib_df, aes(x = sample_id, y = library_size, fill = group)) + geom_col(width = 0.75) + scale_fill_manual(values = GROUP_COLORS, drop = FALSE) + scale_y_continuous(labels = label_compact_count) + labs(title = paste0(TISSUE_TITLE, " library sizes"), x = NULL, y = "Total assigned gene counts", fill = "Group") + theme_rna(10) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
save_plot(p_lib, "QC_library_size_raw_counts", 8.5, 4.5)

keep_gene <- rowSums(count_mat >= MIN_COUNT) >= MIN_SAMPLES
filter_table <- tibble::tibble(Total_genes = nrow(count_mat), Genes_after_filter = sum(keep_gene), Removed_low_count_genes = sum(!keep_gene), MIN_COUNT = MIN_COUNT, MIN_SAMPLES = MIN_SAMPLES)
write.csv(filter_table, file.path(OUTDIR, "QC_gene_filter_summary.csv"), row.names = FALSE)
count_mat_filt <- count_mat[keep_gene, , drop = FALSE]
DESIGN_FORMULA <- as.formula(DESIGN_FORMULA_TEXT)
missing_design_vars <- setdiff(all.vars(DESIGN_FORMULA), names(meta))
if (length(missing_design_vars) > 0) stop("RNA_DESIGN variable(s) not found in metadata: ", paste(missing_design_vars, collapse = ", "))

dds <- DESeqDataSetFromMatrix(countData = count_mat_filt, colData = meta, design = DESIGN_FORMULA)
dds$group <- relevel(dds$group, ref = "XX")
dds <- DESeq(dds)

norm_counts <- counts(dds, normalized = TRUE) %>% as.data.frame() %>% tibble::rownames_to_column("gene_id") %>% dplyr::left_join(annot, by = "gene_id")
write.csv(norm_counts, file.path(OUTDIR, "Table2_DESeq2_normalized_counts.csv"), row.names = FALSE)
vsd <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)
write.csv(vst_mat %>% as.data.frame() %>% tibble::rownames_to_column("gene_id"), file.path(OUTDIR, "Table3_VST_expression_matrix.csv"), row.names = FALSE)

pca_df <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pca_df, "percentVar"))
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group, label = name)) + geom_point(size = 3, alpha = 0.9) + ggrepel::geom_text_repel(size = 2.6, show.legend = FALSE, max.overlaps = Inf) + stat_ellipse(aes(group = group), linewidth = 0.6, alpha = 0.6, show.legend = FALSE) + scale_color_manual(values = GROUP_COLORS, drop = FALSE) + labs(title = paste0("PCA of ", TISSUE_TITLE), x = paste0("PC1: ", percentVar[1], "% variance"), y = paste0("PC2: ", percentVar[2], "% variance"), color = "Group") + theme_rna(11)
save_plot(p_pca, "FigRNA1A_PCA_VST", 6.5, 5.2)

sample_dist_mat <- as.matrix(dist(t(vst_mat))); rownames(sample_dist_mat) <- colnames(vst_mat); colnames(sample_dist_mat) <- colnames(vst_mat)
ann_col <- meta %>% dplyr::select(group) %>% as.data.frame(); rownames(ann_col) <- meta$sample_id
png(file.path(OUTDIR, "FigRNA1B_sample_distance_heatmap.png"), width = 2200, height = 2000, res = 300)
pheatmap::pheatmap(sample_dist_mat, annotation_col = ann_col, annotation_row = ann_col, annotation_colors = list(group = GROUP_COLORS), main = paste0("Sample distance: ", TISSUE_TITLE), fontsize = 8)
dev.off()
pdf(file.path(OUTDIR, "FigRNA1B_sample_distance_heatmap.pdf"), width = 7.5, height = 7)
pheatmap::pheatmap(sample_dist_mat, annotation_col = ann_col, annotation_row = ann_col, annotation_colors = list(group = GROUP_COLORS), main = paste0("Sample distance: ", TISSUE_TITLE), fontsize = 8)
dev.off()

# =============================================================================
# 6. Differential expression
# =============================================================================

prepare_volcano_data <- function(res_df) {
  res_df %>% dplyr::mutate(raw_lfc = log2FoldChange, raw_sig = .data[[SIGNIFICANCE_COLUMN]], sig_for_plot = dplyr::case_when(is.na(raw_sig) ~ NA_real_, raw_sig <= 0 ~ 10^(-VOLCANO_Y_CAP), TRUE ~ raw_sig), neg_log10_sig = pmin(-log10(sig_for_plot), VOLCANO_Y_CAP), log2FC_plot = pmax(pmin(raw_lfc, VOLCANO_X_CAP), -VOLCANO_X_CAP), DEG_status = dplyr::case_when(is.na(raw_sig) ~ "Not tested", raw_sig < ALPHA & raw_lfc >= LFC_CUTOFF & !is.na(baseMean) & baseMean >= MIN_BASEMEAN_FOR_DEG_CALL ~ "Up", raw_sig < ALPHA & raw_lfc <= -LFC_CUTOFF & !is.na(baseMean) & baseMean >= MIN_BASEMEAN_FOR_DEG_CALL ~ "Down", TRUE ~ "Not significant"), label_gene = dplyr::coalesce(symbol, gene_id)) %>% dplyr::filter(DEG_status != "Not tested")
}
make_volcano <- function(res_df, contrast_label, label_n = 10) {
  df <- prepare_volcano_data(res_df)
  lab_df <- df %>% dplyr::filter(DEG_status %in% c("Up", "Down")) %>% dplyr::arrange(raw_sig, dplyr::desc(abs(raw_lfc))) %>% dplyr::slice_head(n = label_n)
  ggplot(df, aes(x = log2FC_plot, y = neg_log10_sig)) + geom_point(aes(color = DEG_status), size = 0.8, alpha = 0.65) + geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed", color = "grey45", linewidth = 0.35) + geom_hline(yintercept = -log10(ALPHA), linetype = "dashed", color = "grey45", linewidth = 0.35) + ggrepel::geom_text_repel(data = lab_df, aes(label = label_gene), size = 2.5, max.overlaps = Inf, box.padding = 0.25, min.segment.length = 0) + scale_color_manual(values = c("Up" = "#D55E00", "Down" = "#0072B2", "Not significant" = "grey70"), drop = FALSE) + labs(title = paste0("Volcano plot: ", contrast_label), subtitle = paste0(SIGNIFICANCE_LABEL, " < ", ALPHA, ", |log2FC| >= ", LFC_CUTOFF), x = "log2 fold change", y = SIGNIFICANCE_AXIS_LABEL, color = "DEG status") + theme_rna(11)
}

contrast_list <- list(YY_vs_XX = c("group", "YY", "XX"), YX_vs_YY = c("group", "YX", "YY"), YX_vs_XX = c("group", "YX", "XX"))
res_list <- list()
for (nm in names(contrast_list)) {
  message("Running contrast: ", nm)
  res <- results(dds, contrast = contrast_list[[nm]], alpha = ALPHA)
  res_df <- as.data.frame(res) %>% tibble::rownames_to_column("gene_id") %>% dplyr::left_join(annot, by = "gene_id") %>% dplyr::mutate(contrast = nm, significance_column = SIGNIFICANCE_COLUMN, significance_value = .data[[SIGNIFICANCE_COLUMN]], DEG = sig_to_deg(significance_value, log2FoldChange, baseMean), direction = dplyr::case_when(DEG & log2FoldChange > 0 ~ "Up", DEG & log2FoldChange < 0 ~ "Down", TRUE ~ "Not_DEG")) %>% dplyr::arrange(significance_value, padj)
  res_list[[nm]] <- res_df
  write.csv(res_df, file.path(OUTDIR, paste0("DEG_", nm, ".csv")), row.names = FALSE)
  save_plot(make_volcano(res_df, nm), paste0("FigRNA2_volcano_", nm), 6.8, 5.6)
  png(file.path(OUTDIR, paste0("FigRNA2_MA_", nm, ".png")), width = 1800, height = 1500, res = 300); plotMA(res, ylim = c(-6, 6), main = paste0("MA plot: ", nm)); dev.off()
  pdf(file.path(OUTDIR, paste0("FigRNA2_MA_", nm, ".pdf")), width = 6, height = 5); plotMA(res, ylim = c(-6, 6), main = paste0("MA plot: ", nm)); dev.off()
}
all_deg <- dplyr::bind_rows(res_list)
write.csv(all_deg, file.path(OUTDIR, "Table4_all_DEG_results_long.csv"), row.names = FALSE)

deg_count <- all_deg %>% dplyr::filter(direction %in% c("Up", "Down")) %>% dplyr::count(contrast, direction, name = "N") %>% tidyr::complete(contrast = names(contrast_list), direction = c("Up", "Down"), fill = list(N = 0)) %>% dplyr::mutate(filter_mode = DEG_FILTER_MODE, significance_column = SIGNIFICANCE_COLUMN, alpha = ALPHA, lfc_cutoff = LFC_CUTOFF, min_baseMean = MIN_BASEMEAN_FOR_DEG_CALL)
write.csv(deg_count, file.path(OUTDIR, "Table5_DEG_counts_by_contrast.csv"), row.names = FALSE)
p_deg <- ggplot(deg_count, aes(x = contrast, y = N, fill = direction)) + geom_col(position = position_dodge(width = 0.7), width = 0.62) + geom_text(aes(label = N, group = direction), position = position_dodge(width = 0.7), vjust = -0.35, size = 3.2) + scale_fill_manual(values = c("Up" = "#D55E00", "Down" = "#0072B2")) + labs(title = "Differentially expressed genes across contrasts", subtitle = paste0("Main standard: ", DEG_FILTER_MODE, "; ", SIGNIFICANCE_LABEL, " < ", ALPHA, ", |log2FC| >= ", LFC_CUTOFF), x = NULL, y = "Number of DEGs", fill = "Direction") + theme_rna(11) + theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_plot(p_deg, "FigRNA2A_DEG_counts", 7.2, 4.8)

deg_count_by_standard <- purrr::map_dfr(seq_len(nrow(DEG_STANDARD_GRID)), function(i) {
  st <- DEG_STANDARD_GRID[i, ]
  all_deg %>% dplyr::mutate(DEG_standard = sig_by_standard(.data[[st$sig_col]], log2FoldChange, baseMean, st$alpha, st$lfc, st$min_baseMean), direction_standard = dplyr::case_when(DEG_standard & log2FoldChange > 0 ~ "Up", DEG_standard & log2FoldChange < 0 ~ "Down", TRUE ~ "Not_DEG")) %>% dplyr::filter(direction_standard %in% c("Up", "Down")) %>% dplyr::count(contrast, direction = direction_standard, name = "N") %>% tidyr::complete(contrast = names(contrast_list), direction = c("Up", "Down"), fill = list(N = 0)) %>% dplyr::mutate(standard = st$standard, significance_column = st$sig_col, alpha = st$alpha, lfc_cutoff = st$lfc, min_baseMean = st$min_baseMean)
})
write.csv(deg_count_by_standard, file.path(OUTDIR, "Table5B_DEG_counts_by_filter_standard.csv"), row.names = FALSE)
p_deg_standard <- ggplot(deg_count_by_standard, aes(x = standard, y = N, fill = direction)) + geom_col(position = position_dodge(width = 0.72), width = 0.64) + geom_text(aes(label = N, group = direction), position = position_dodge(width = 0.72), vjust = -0.25, size = 2.8) + facet_wrap(~ contrast, scales = "free_y", ncol = 3) + scale_fill_manual(values = c("Up" = "#D55E00", "Down" = "#0072B2")) + labs(title = "DEG counts under alternative filtering standards", subtitle = "FDR-based standards plus nominal-P exploratory option", x = "Filtering standard", y = "Number of DEGs", fill = "Direction") + theme_rna(10) + theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_plot(p_deg_standard, "FigRNA2B_DEG_counts_by_filter_standard", 10.5, 4.6)

# =============================================================================
# 7. Duodenum candidate module scores and compact heatmap
# =============================================================================

symbol_to_gene <- annot %>% dplyr::mutate(symbol_upper = toupper(symbol)) %>% dplyr::filter(!is.na(symbol_upper), symbol_upper != "") %>% dplyr::select(gene_id, symbol, display_gene, gene_name, entrez, symbol_upper) %>% dplyr::distinct(gene_id, .keep_all = TRUE)
candidate_matches <- purrr::imap_dfr(DUODENUM_GENE_SETS, function(symbols, set_name) {
  symbol_to_gene %>% dplyr::filter(symbol_upper %in% toupper(symbols)) %>% dplyr::mutate(gene_set = set_name)
}) %>% dplyr::distinct(gene_set, gene_id, .keep_all = TRUE)
write.csv(candidate_matches, file.path(OUTDIR, "Table_Duodenum_candidate_module_gene_matches.csv"), row.names = FALSE)

gene_z <- t(scale(t(vst_mat))); gene_z[is.na(gene_z)] <- 0
module_scores <- purrr::map_dfr(unique(candidate_matches$gene_set), function(gs) {
  genes <- intersect(candidate_matches$gene_id[candidate_matches$gene_set == gs], rownames(gene_z))
  if (length(genes) < 2) return(NULL)
  score <- colMeans(gene_z[genes, , drop = FALSE], na.rm = TRUE)
  tibble::tibble(sample_id = names(score), module_score = as.numeric(score), gene_set = gs, n_genes = length(genes))
}) %>% dplyr::left_join(meta %>% dplyr::select(sample_id, group), by = "sample_id")
write.csv(module_scores, file.path(OUTDIR, "Table_Duodenum_candidate_module_scores.csv"), row.names = FALSE)
if (nrow(module_scores) > 0) {
  p_module <- ggplot(module_scores, aes(x = group, y = module_score, fill = group)) + geom_boxplot(alpha = 0.35, outlier.shape = NA, width = 0.65) + geom_jitter(aes(color = group), width = 0.12, size = 2, alpha = 0.85, show.legend = FALSE) + facet_wrap(~ gene_set, scales = "free_y", ncol = 3) + scale_fill_manual(values = GROUP_COLORS, drop = FALSE, guide = "none") + scale_color_manual(values = GROUP_COLORS, drop = FALSE, guide = "none") + labs(title = "Targeted pathway scores in duodenum", subtitle = "Mean Z-scored VST expression of matched genes in each candidate gene set", x = NULL, y = "Module score") + theme_rna(10) + theme(strip.text = element_text(size = 8.5))
  save_plot(p_module, "Fig_Duodenum_targeted_pathway_scores", 10.5, 7.2)
}

if (nrow(candidate_matches) > 0) {
  group_means_for_heat <- vst_mat[intersect(candidate_matches$gene_id, rownames(vst_mat)), , drop = FALSE] %>% as.data.frame() %>% tibble::rownames_to_column("gene_id") %>% tidyr::pivot_longer(-gene_id, names_to = "sample_id", values_to = "vst") %>% dplyr::left_join(meta %>% dplyr::select(sample_id, group), by = "sample_id") %>% dplyr::group_by(gene_id, group) %>% dplyr::summarise(mean_vst = mean(vst, na.rm = TRUE), .groups = "drop") %>% tidyr::pivot_wider(names_from = group, values_from = mean_vst, names_prefix = "mean_") %>% dplyr::left_join(candidate_matches %>% dplyr::select(gene_set, gene_id, display_gene, symbol), by = "gene_id") %>% dplyr::mutate(row_var = apply(dplyr::select(., dplyr::starts_with("mean_")), 1, stats::var, na.rm = TRUE))
  heat_select <- group_means_for_heat %>% dplyr::group_by(gene_set) %>% dplyr::arrange(dplyr::desc(row_var), .by_group = TRUE) %>% dplyr::slice_head(n = TOP_N_HEATMAP_PER_SET) %>% dplyr::ungroup() %>% dplyr::arrange(gene_set, dplyr::desc(row_var)) %>% dplyr::slice_head(n = MAX_HEATMAP_GENES)
  write.csv(heat_select, file.path(OUTDIR, "Table_Duodenum_candidate_heatmap_gene_selection.csv"), row.names = FALSE)
  if (nrow(heat_select) >= 2) {
    heat_mat <- heat_select %>% dplyr::select(mean_YY, mean_YX, mean_XX) %>% as.matrix()
    rownames(heat_mat) <- make.unique(heat_select$display_gene)
    heat_z <- t(scale(t(heat_mat))); heat_z[is.na(heat_z)] <- 0
    row_annot <- data.frame(gene_set = heat_select$gene_set); rownames(row_annot) <- rownames(heat_z)
    pdf(file.path(OUTDIR, "Fig_Duodenum_targeted_candidate_gene_heatmap_compact.pdf"), width = 7.2, height = max(5.5, min(12, nrow(heat_z) * 0.18 + 2)))
    pheatmap::pheatmap(heat_z, annotation_row = row_annot, color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101), cluster_cols = FALSE, cluster_rows = FALSE, show_rownames = TRUE, fontsize_row = 6, fontsize_col = 10, border_color = NA, main = "Duodenum lipid absorption and transport candidate genes")
    dev.off()
    png(file.path(OUTDIR, "Fig_Duodenum_targeted_candidate_gene_heatmap_compact.png"), width = 2200, height = max(1800, min(3600, nrow(heat_z) * 55 + 650)), res = 300)
    pheatmap::pheatmap(heat_z, annotation_row = row_annot, color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101), cluster_cols = FALSE, cluster_rows = FALSE, show_rownames = TRUE, fontsize_row = 6, fontsize_col = 10, border_color = NA, main = "Duodenum lipid absorption and transport candidate genes")
    dev.off()
  }
}

# =============================================================================
# 8. ORA and GSEA
# =============================================================================

plot_enrich_dot <- function(enrich_obj, title, filename, show_category = ENRICH_SHOW_CATEGORY, width = 9, height = 6) {
  if (is.null(enrich_obj)) return(NULL)
  df <- as.data.frame(enrich_obj)
  if (nrow(df) == 0) return(NULL)
  df <- df %>% dplyr::mutate(GeneRatio_num = parse_gene_ratio(GeneRatio), Description_short = shorten_term(Description), Cluster = as.character(Cluster)) %>% dplyr::group_by(Cluster) %>% dplyr::arrange(p.adjust, dplyr::desc(GeneRatio_num), .by_group = TRUE) %>% dplyr::slice_head(n = show_category) %>% dplyr::ungroup()
  p <- ggplot(df, aes(x = Cluster, y = reorder(Description_short, p.adjust))) + geom_point(aes(size = GeneRatio_num, color = p.adjust), alpha = 0.9) + scale_color_gradient(low = "#D73027", high = "#4575B4", trans = "reverse", name = "FDR") + scale_size_continuous(name = "Gene ratio") + labs(title = title, x = NULL, y = NULL) + theme_rna(10) + theme(axis.text.x = element_text(angle = 35, hjust = 1), axis.text.y = element_text(size = 8))
  save_plot(p, filename, width, height)
}

run_kegg_compare <- function(gene_sets, out_prefix, title) {
  if (!HAS_CLUSTERPROFILER || length(gene_sets) == 0) return(NULL)
  suppressPackageStartupMessages(library(clusterProfiler))
  res <- tryCatch(clusterProfiler::compareCluster(geneCluster = gene_sets, fun = "enrichKEGG", organism = "gga", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.20), error = function(e) { message(out_prefix, " KEGG failed: ", e$message); NULL })
  if (!is.null(res)) { write.csv(as.data.frame(res), file.path(OUTDIR, paste0(out_prefix, "_KEGG.csv")), row.names = FALSE); plot_enrich_dot(res, title, paste0("Fig", out_prefix, "_KEGG_dotplot")) }
  res
}

run_go_compare <- function(gene_sets, out_prefix, title, ont = "BP") {
  if (!HAS_CLUSTERPROFILER || !HAS_ORG_GG || length(gene_sets) == 0) return(NULL)
  suppressPackageStartupMessages(library(clusterProfiler)); suppressPackageStartupMessages(library(org.Gg.eg.db))
  res <- tryCatch(clusterProfiler::compareCluster(geneCluster = gene_sets, fun = "enrichGO", OrgDb = org.Gg.eg.db, keyType = "ENTREZID", ont = ont, pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.20, readable = TRUE), error = function(e) { message(out_prefix, " GO failed: ", e$message); NULL })
  if (!is.null(res)) { write.csv(as.data.frame(res), file.path(OUTDIR, paste0(out_prefix, "_GO_", ont, ".csv")), row.names = FALSE); plot_enrich_dot(res, paste0(title, " GO ", ont), paste0("Fig", out_prefix, "_GO_", ont, "_dotplot")) }
  res
}

if (RUN_DEG_ENRICHMENT && "entrez" %in% names(all_deg)) {
  deg_sets_df <- all_deg %>% dplyr::filter(DEG, direction %in% c("Up", "Down"), !is.na(entrez), entrez != "") %>% dplyr::mutate(Set = paste(contrast, direction, sep = "_")) %>% dplyr::distinct(Set, entrez)
  deg_gene_sets <- split(as.character(deg_sets_df$entrez), deg_sets_df$Set); deg_gene_sets <- lapply(deg_gene_sets, unique); deg_gene_sets <- deg_gene_sets[lengths(deg_gene_sets) >= ENRICH_MIN_GENESET_SIZE]
  write.csv(tibble::tibble(Set = names(deg_gene_sets), N_entrez = lengths(deg_gene_sets)), file.path(OUTDIR, "Enrichment_input_DEG_sets_by_contrast_direction.csv"), row.names = FALSE)
  if (length(deg_gene_sets) > 0) { run_kegg_compare(deg_gene_sets, "RNA4A_DEG_by_contrast_direction", "KEGG enrichment of DEGs by contrast and direction"); run_go_compare(deg_gene_sets, "RNA4B_DEG_by_contrast_direction", "DEG enrichment by contrast and direction", ont = "BP"); run_go_compare(deg_gene_sets, "RNA4C_DEG_by_contrast_direction", "DEG enrichment by contrast and direction", ont = "CC") }
}

make_ranked_genelist <- function(res_df, entrez_col = "entrez", metric = GSEA_RANK_METRIC) {
  df <- res_df %>% dplyr::filter(!is.na(.data[[entrez_col]]), .data[[entrez_col]] != "")
  df$rank_metric <- switch(metric, stat = df$stat, log2fc = df$log2FoldChange, signed_log10p = sign(df$log2FoldChange) * -log10(pmax(df$pvalue, 1e-300)), df$stat)
  df <- df %>% dplyr::filter(is.finite(rank_metric))
  if (nrow(df) == 0) return(numeric(0))
  df <- df %>% dplyr::mutate(entrez_tmp = as.character(.data[[entrez_col]])) %>% dplyr::group_by(entrez_tmp) %>% dplyr::slice_max(order_by = abs(rank_metric), n = 1, with_ties = FALSE) %>% dplyr::ungroup()
  gl <- df$rank_metric; names(gl) <- df$entrez_tmp; sort(gl, decreasing = TRUE)
}

run_gsea_kegg_one <- function(ranked_list, out_prefix) {
  if (!HAS_CLUSTERPROFILER || length(ranked_list) < GSEA_MIN_GSSIZE) return(NULL)
  suppressPackageStartupMessages(library(clusterProfiler))
  res <- tryCatch(clusterProfiler::gseKEGG(geneList = ranked_list, organism = "gga", minGSSize = GSEA_MIN_GSSIZE, maxGSSize = GSEA_MAX_GSSIZE, pvalueCutoff = 1, pAdjustMethod = "BH", seed = TRUE, eps = 0), error = function(e) { message(out_prefix, " GSEA-KEGG failed: ", e$message); NULL })
  if (!is.null(res) && nrow(as.data.frame(res)) > 0) write.csv(as.data.frame(res), file.path(OUTDIR, paste0(out_prefix, "_GSEA_KEGG.csv")), row.names = FALSE)
  res
}
run_gsea_go_one <- function(ranked_list, out_prefix, ont = "BP") {
  if (!HAS_CLUSTERPROFILER || !HAS_ORG_GG || length(ranked_list) < GSEA_MIN_GSSIZE) return(NULL)
  suppressPackageStartupMessages(library(clusterProfiler)); suppressPackageStartupMessages(library(org.Gg.eg.db))
  res <- tryCatch(clusterProfiler::gseGO(geneList = ranked_list, OrgDb = org.Gg.eg.db, keyType = "ENTREZID", ont = ont, minGSSize = GSEA_MIN_GSSIZE, maxGSSize = GSEA_MAX_GSSIZE, pvalueCutoff = 1, pAdjustMethod = "BH", seed = TRUE, eps = 0), error = function(e) { message(out_prefix, " GSEA-GO failed: ", e$message); NULL })
  if (!is.null(res) && nrow(as.data.frame(res)) > 0) { res <- tryCatch(clusterProfiler::setReadable(res, OrgDb = org.Gg.eg.db, keyType = "ENTREZID"), error = function(e) res); write.csv(as.data.frame(res), file.path(OUTDIR, paste0(out_prefix, "_GSEA_GO_", ont, ".csv")), row.names = FALSE) }
  res
}

gsea_all_long <- data.frame(); gsea_sig <- data.frame()
if (RUN_GSEA && "entrez" %in% names(all_deg)) {
  gsea_list <- list()
  for (nm in names(res_list)) {
    ranked <- make_ranked_genelist(res_list[[nm]])
    write.csv(tibble::tibble(entrez = names(ranked), rank_metric = as.numeric(ranked)), file.path(OUTDIR, paste0("Enrichment_input_GSEA_ranked_list_", nm, ".csv")), row.names = FALSE)
    if (length(ranked) < GSEA_MIN_GSSIZE) next
    gsea_list[[paste0(nm, "_KEGG")]] <- run_gsea_kegg_one(ranked, paste0("RNA6A_GSEA_", nm))
    gsea_list[[paste0(nm, "_GO_BP")]] <- run_gsea_go_one(ranked, paste0("RNA6B_GSEA_", nm), ont = "BP")
  }
  gsea_all_long <- dplyr::bind_rows(lapply(names(gsea_list), function(nm) {
    x <- gsea_list[[nm]]; if (is.null(x)) return(NULL); df <- as.data.frame(x); if (nrow(df) == 0) return(NULL)
    parts <- stringr::str_match(nm, "(.+)_(KEGG|GO_BP)$"); df$contrast <- parts[, 2]; df$ontology <- parts[, 3]; df
  }))
  if (nrow(gsea_all_long) > 0) {
    write.csv(gsea_all_long, file.path(OUTDIR, "Table14_GSEA_all_contrasts_long.csv"), row.names = FALSE)
    gsea_sig <- gsea_all_long %>% dplyr::filter(p.adjust < GSEA_PADJ_REPORT_CUTOFF)
    write.csv(gsea_sig, file.path(OUTDIR, "Table15_GSEA_suggestive_q25_pathways.csv"), row.names = FALSE)
    if (nrow(gsea_sig) > 0) {
      plot_df <- gsea_sig %>% dplyr::mutate(Description_short = shorten_term(Description)) %>% dplyr::group_by(contrast, ontology) %>% dplyr::slice_min(order_by = p.adjust, n = 8, with_ties = FALSE) %>% dplyr::ungroup()
      p_gsea <- ggplot(plot_df, aes(x = contrast, y = reorder(Description_short, NES))) + geom_point(aes(size = -log10(p.adjust), color = NES)) + scale_color_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0, name = "NES") + scale_size_continuous(name = "-log10(FDR q)") + facet_wrap(~ ontology, scales = "free_y") + labs(title = "GSEA: suggestive pathways (FDR q < 0.25)", x = NULL, y = NULL) + theme_rna(9) + theme(axis.text.y = element_text(size = 7), axis.text.x = element_text(angle = 30, hjust = 1))
      save_plot(p_gsea, "FigRNA6C_GSEA_comparison_bubble", 10, 7.5)
    }
  }
}

# =============================================================================
# 9. Export workbook
# =============================================================================

writexl::write_xlsx(list(sample_counts = sample_summary, gene_filter_summary = filter_table, annotation_summary = annotation_summary, DEG_counts = deg_count, DEG_counts_by_standard = deg_count_by_standard, duodenum_candidate_gene_matches = candidate_matches, duodenum_candidate_module_scores = module_scores, YY_vs_XX = res_list$YY_vs_XX, YX_vs_YY = res_list$YX_vs_YY, YX_vs_XX = res_list$YX_vs_XX, GSEA_all_contrasts_long = gsea_all_long, GSEA_suggestive_q25 = gsea_sig), path = file.path(OUTDIR, "duodenum_transcriptome_core_results.xlsx"))

message("\n============================================================")
message("Duodenum transcriptome CORE analysis complete.")
message("Output: ", OUTDIR)
message("Design: ", DESIGN_FORMULA_TEXT)
message("Main DEG threshold: ", SIGNIFICANCE_LABEL, " < ", ALPHA, ", |log2FC| >= ", LFC_CUTOFF)
message("Next: run 05C_transcriptome_F1_inheritance_pattern_standalone.R if you need F1 inheritance/nonadditive analysis.")
message("============================================================\n")
