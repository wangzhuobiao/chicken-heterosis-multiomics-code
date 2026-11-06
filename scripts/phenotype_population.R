###############################################################################
# Phenotype story analysis for YY, YX F1 and XX chickens - Article 1 core v4.1
# Article 1 publication boundary:
#   - Main: longitudinal growth, overall carcass-weight traits and abdominal fat.
#   - Supplementary: interval ADG and body measurements.
#   - Reserved for Article 2: breast/leg yield and meat-quality traits.
# Revised focus:
#   1. QC raw traits before deriving ADG and slaughter ratios.
#   2. Missing values are handled trait-by-trait.
#   3. Negative ADG intervals are exported and optionally set to NA.
#   4. Main paper figures focus on growth, core carcass traits and F1 relative
#      position instead of over-fragmented ADG panels and one-trait-per-bar pattern.
#   5. Body measurements and full ADG interval plots are moved to supplementary use.
#
# Input default:
###############################################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

# =============================================================================
# 0. User parameters
# =============================================================================

COMPLETE_FILE <- Sys.getenv(
  "PHENO_FILE",
  unset = ""
)
COMPLETE_SHEET <- as.integer(Sys.getenv("PHENO_SHEET", unset = "1"))

OUTDIR <- Sys.getenv(
  "PHENO_OUTDIR",
  unset = ""
)
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

GROUP_LEVELS <- c("YY", "YX", "XX")
GROUP_COLORS <- c(
  "YY" = "#E64B35",
  "YX" = "#4DBBD5",
  "XX" = "#00A087"
)

APPLY_HARD_OUTLIERS_TO_NA <- TRUE
APPLY_ROBUST_OUTLIERS_TO_NA <- FALSE
FILTER_NEGATIVE_GROWTH_INDIVIDUALS <- TRUE
SET_NEGATIVE_ADG_TO_NA <- TRUE
INCLUDE_BODY_MEASUREMENTS_IN_MAIN <- FALSE
# The previous broad phenotype-correlation heatmap mixed Article 1 and Article 2 traits.
# It is not needed for the Article 1 main story, so it is OFF by default.
MAKE_CORE_CORRELATION_HEATMAP <- FALSE
BOOT_N <- as.integer(Sys.getenv("PHENO_BOOT_N", unset = "2000"))
set.seed(as.integer(Sys.getenv("PHENO_SEED", unset = "20260610")))

ARTICLE1_CORE_TRAITS <- c(
  "BW_wk15", "ADG_0_15",
  "Live_wt", "Carcass_wt", "Evis_wt",
  "AF_wt", "AF_live_pct"
)

ARTICLE1_RESERVED_FOR_ARTICLE2 <- c(
  "Breast_wt", "Breast_live_pct", "Breast_evis_pct",
  "Leg_wt", "Leg_live_pct",
  "PH_mean", "Drip_24h", "L_mean", "a_mean", "b_mean", "Shear_mean"
)

LOWER_IS_BETTER_TRAITS <- c("AF_wt", "AF_live_pct")

# =============================================================================
# 1. Packages
# =============================================================================

cran_pkgs <- c("tidyverse", "readxl", "writexl", "patchwork", "scales")
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

for (pkg in cran_pkgs) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    stop(paste("Install required R package:", pkg))
    library(pkg, character.only = TRUE)
  }
}

# =============================================================================
# 2. Helper functions
# =============================================================================

is_nonempty <- function(x) {
  x <- as.character(x)
  !is.na(x) &
    stringr::str_trim(x) != "" &
    !stringr::str_to_lower(stringr::str_trim(x)) %in% c("na", "nan", "null", "none")
}

safe_numeric <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, ",", "")
  x <- stringr::str_replace_all(x, "，", "")
  x <- stringr::str_replace_all(x, "%", "")
  x <- stringr::str_trim(x)
  x[x %in% c("", "NA", "NaN", "null", "NULL", "-", "--")] <- NA_character_
  suppressWarnings(as.numeric(x))
}

normalize_group <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x <- stringr::str_replace_all(x, "F1", "YX")
  x <- stringr::str_replace_all(x, "Y×X", "YX")
  x <- stringr::str_replace_all(x, "Yx", "YX")
  x <- stringr::str_replace_all(x, "yx", "YX")
  x <- stringr::str_replace_all(x, "杂交", "YX")
  factor(x, levels = GROUP_LEVELS)
}

normalize_id <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x[x %in% c("", "NA", "NaN", "NULL", "null")] <- NA_character_
  x
}

make_numeric_if_exists <- function(df, vars) {
  vars <- vars[vars %in% names(df)]
  if (length(vars) == 0) return(df)
  df %>% mutate(across(all_of(vars), safe_numeric))
}

ensure_cols <- function(df, vars, value = NA_real_) {
  for (v in vars) {
    if (!v %in% names(df)) df[[v]] <- value
  }
  df
}

mean_if_available <- function(df, cols, out_name) {
  cols <- cols[cols %in% names(df)]
  if (length(cols) == 0) {
    df[[out_name]] <- NA_real_
  } else {
    df[[out_name]] <- rowMeans(df[, cols, drop = FALSE], na.rm = TRUE)
    df[[out_name]][is.nan(df[[out_name]])] <- NA_real_
  }
  df
}

safe_ratio <- function(num, den, multiplier = 100) {
  out <- ifelse(!is.na(num) & !is.na(den) & den != 0, num / den * multiplier, NA_real_)
  as.numeric(out)
}

se <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) <= 1) return(NA_real_)
  sd(x) / sqrt(length(x))
}

save_plot <- function(plot_obj, filename, width, height, dpi = 300) {
  ggsave(file.path(OUTDIR, paste0(filename, ".png")), plot_obj, width = width, height = height, dpi = dpi)
  ggsave(file.path(OUTDIR, paste0(filename, ".pdf")), plot_obj, width = width, height = height)
}

theme_manuscript <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2),
      plot.subtitle = element_text(size = base_size),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      strip.background = element_rect(fill = "grey95", color = "grey70"),
      strip.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.25),
      panel.grid.major.x = element_blank()
    )
}

label_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

# =============================================================================
# 3. Read and standardize data
# =============================================================================

if (!file.exists(COMPLETE_FILE)) stop("Input file does not exist: ", COMPLETE_FILE)

raw <- readxl::read_excel(COMPLETE_FILE, sheet = COMPLETE_SHEET, .name_repair = "unique")
raw <- raw %>% select(where(~ !all(is.na(.x))))
message("Input rows: ", nrow(raw))
message("Input columns: ", ncol(raw))

rename_map <- c(
  "组别" = "group_raw",
  "翅号" = "wing_id",
  "脚号" = "leg_id",
  "铝制脚号" = "aluminum_leg_id",
  "胸肌转录组编号" = "muscle_rna_id",
  "盲肠宏基因组" = "cecum_meta_id",
  "十二指肠转录组送样编号" = "duodenum_rna_id",
  "性别（矫正）" = "sex_corrected",
  "初生重" = "BW_birth",
  "第一周体重" = "BW_wk1",
  "第二周体重" = "BW_wk2",
  "第三周体重" = "BW_wk3",
  "第五周体重" = "BW_wk5",
  "第七周体重" = "BW_wk7",
  "第九周体重" = "BW_wk9",
  "第十一周体重" = "BW_wk11",
  "第十三周体重" = "BW_wk13",
  "第十五周体重" = "BW_wk15",
  "第9周髋宽(mm)" = "HW9",
  "第9周胫长(mm)" = "SL9",
  "第9周胫围(cm)" = "SC9",
  "第9周体斜长(cm)" = "BSL9",
  "第9周龙骨长(cm)" = "KL9",
  "第11周髋宽(mm)" = "HW11",
  "第11周胫长(mm)" = "SL11",
  "第11周胫围(cm)" = "SC11",
  "第11周体斜长(cm)" = "BSL11",
  "第11周龙骨长(cm)" = "KL11",
  "第13周髋宽(mm)" = "HW13",
  "第13周胫长(mm)" = "SL13",
  "第13周胫围(cm)" = "SC13",
  "第13周体斜长(cm)" = "BSL13",
  "第13周龙骨长(cm)" = "KL13",
  "第15周髋宽(mm)" = "HW15",
  "第15周胫长(mm)" = "SL15",
  "第15周胫围(cm)" = "SC15",
  "第15周体斜长(cm)" = "BSL15",
  "第15周龙骨长(cm)" = "KL15",
  "屠体重" = "Carcass_wt",
  "半净膛重" = "Semi_evis_wt",
  "全净堂重" = "Evis_wt",
  "腿重(单侧)" = "Leg_wt",
  "胸肌重(单侧)" = "Breast_wt",
  "翅膀重(双侧）" = "Wing_wt",
  "头重" = "Head_wt",
  "腹脂重" = "AF_wt",
  "爪子重" = "Paw_wt",
  "活重" = "Live_wt",
  "髋宽_屠宰(mm)" = "HW_slaughter",
  "胫长_屠宰(mm)" = "SL_slaughter",
  "胫围_屠宰(cm)" = "SC_slaughter",
  "体斜长_屠宰(cm)" = "BSL_slaughter",
  "龙骨长_屠宰(cm)" = "KL_slaughter",
  "滴水损失(挂前)" = "Drip_pre",
  "滴水损失(24h)" = "Drip_24h",
  "剪切力1" = "Shear1",
  "剪切力2" = "Shear2",
  "剪切力3" = "Shear3"
)

for (old in names(rename_map)) {
  if (old %in% names(raw)) names(raw)[names(raw) == old] <- rename_map[[old]]
}

color_map <- c(
  "L" = "L1", "L.1" = "L2", "L.2" = "L3",
  "a" = "a1", "a.1" = "a2", "a.2" = "a3",
  "b" = "b1", "b.1" = "b2", "b.2" = "b3"
)
for (old in names(color_map)) {
  if (old %in% names(raw)) names(raw)[names(raw) == old] <- color_map[[old]]
}

id_vars <- c("group_raw", "wing_id", "leg_id", "aluminum_leg_id", "muscle_rna_id", "cecum_meta_id", "duodenum_rna_id", "sex_corrected")
for (v in id_vars) if (!v %in% names(raw)) raw[[v]] <- NA_character_

bw_vars <- c("BW_birth", "BW_wk1", "BW_wk2", "BW_wk3", "BW_wk5", "BW_wk7", "BW_wk9", "BW_wk11", "BW_wk13", "BW_wk15")
# Canonical body-weight order used for longitudinal QC and deletion logs.
# Keeping this global avoids object-not-found errors when running QC chunks interactively.
bw_order <- bw_vars
body_vars <- c(
  "HW9", "SL9", "SC9", "BSL9", "KL9",
  "HW11", "SL11", "SC11", "BSL11", "KL11",
  "HW13", "SL13", "SC13", "BSL13", "KL13",
  "HW15", "SL15", "SC15", "BSL15", "KL15",
  "HW_slaughter", "SL_slaughter", "SC_slaughter", "BSL_slaughter", "KL_slaughter"
)
slaughter_vars <- c("Live_wt", "Carcass_wt", "Semi_evis_wt", "Evis_wt", "Leg_wt", "Breast_wt", "Wing_wt", "Head_wt", "AF_wt", "Paw_wt")
meat_raw_vars <- c("PH1", "PH2", "PH3", "Drip_pre", "Drip_24h", "L1", "L2", "L3", "a1", "a2", "a3", "b1", "b2", "b3", "Shear1", "Shear2", "Shear3")
all_raw_numeric <- unique(c(bw_vars, body_vars, slaughter_vars, meat_raw_vars))

raw <- ensure_cols(raw, all_raw_numeric, NA_real_)

dat_raw <- raw %>%
  mutate(
    group = normalize_group(group_raw),
    wing_id = normalize_id(wing_id),
    leg_id = normalize_id(leg_id),
    aluminum_leg_id = normalize_id(aluminum_leg_id),
    muscle_rna_id = normalize_id(muscle_rna_id),
    cecum_meta_id = normalize_id(cecum_meta_id),
    duodenum_rna_id = normalize_id(duodenum_rna_id),
    is_muscle_rna = is_nonempty(muscle_rna_id),
    is_cecum_meta = is_nonempty(cecum_meta_id),
    is_duodenum_rna = is_nonempty(duodenum_rna_id),
    # Article 1 only: breast-muscle RNA is reserved for Article 2.
    # Plasma sample IDs are not stored in this phenotype sheet; matched three-omics
    # counts are therefore handled by the dedicated integration/Figure-1 scripts.
    is_any_omics = is_cecum_meta | is_duodenum_rna,
    # Legacy audit flag only; never use this as the Article-1 three-omics count.
    is_all_three_omics = is_muscle_rna & is_cecum_meta & is_duodenum_rna
  ) %>%
  make_numeric_if_exists(all_raw_numeric) %>%
  filter(!is.na(group), !is.na(wing_id))

# =============================================================================
# 4. Trait labels
# =============================================================================

trait_labels <- c(
  BW_birth = "Hatch weight", BW_wk1 = "BW wk1", BW_wk2 = "BW wk2", BW_wk3 = "BW wk3",
  BW_wk5 = "BW wk5", BW_wk7 = "BW wk7", BW_wk9 = "BW wk9", BW_wk11 = "BW wk11",
  BW_wk13 = "BW wk13", BW_wk15 = "BW wk15",
  ADG_0_1 = "ADG 0-1 wk", ADG_1_2 = "ADG 1-2 wk", ADG_2_3 = "ADG 2-3 wk",
  ADG_3_5 = "ADG 3-5 wk", ADG_5_7 = "ADG 5-7 wk", ADG_7_9 = "ADG 7-9 wk",
  ADG_9_11 = "ADG 9-11 wk", ADG_11_13 = "ADG 11-13 wk", ADG_13_15 = "ADG 13-15 wk",
  ADG_0_15 = "ADG 0-15 wk",
  HW15 = "Hip width wk15", SL15 = "Shank length wk15", SC15 = "Shank circ. wk15",
  BSL15 = "Body slant length wk15", KL15 = "Keel length wk15",
  Live_wt = "Live weight", Carcass_wt = "Carcass weight", Semi_evis_wt = "Semi-evis. weight",
  Evis_wt = "Eviscerated weight", Leg_wt = "Leg weight", Breast_wt = "Breast weight",
  Wing_wt = "Wing weight", Head_wt = "Head weight", AF_wt = "Abdominal fat", Paw_wt = "Paw weight",
  Carcass_live_pct = "Carcass/live (%)", Semi_evis_live_pct = "Semi-evis./live (%)",
  Evis_live_pct = "Evis./live (%)", Leg_live_pct = "Leg/live (%)",
  Breast_live_pct = "Breast/live (%)", Breast_evis_pct = "Breast/evis. (%)",
  Wing_live_pct = "Wing/live (%)", Head_live_pct = "Head/live (%)",
  AF_live_pct = "Abdominal fat/live (%)", Paw_live_pct = "Paw/live (%)",
  PH_mean = "pH mean", Drip_24h = "Drip loss 24 h", L_mean = "L*", a_mean = "a*", b_mean = "b*",
  Shear_mean = "Shear force"
)

label_trait <- function(x) {
  x <- as.character(x)
  out <- unname(trait_labels[x])
  out[is.na(out)] <- x[is.na(out)]
  out
}

# =============================================================================
# 5. QC before derived traits
# =============================================================================

missing_summary_raw <- dat_raw %>%
  summarise(across(any_of(all_raw_numeric), ~ sum(is.na(.x)))) %>%
  pivot_longer(everything(), names_to = "Trait", values_to = "N_missing") %>%
  mutate(Trait_label = label_trait(Trait), N_total = nrow(dat_raw), Missing_percent = N_missing / N_total * 100) %>%
  arrange(desc(Missing_percent), Trait)
write.csv(missing_summary_raw, file.path(OUTDIR, "QC_missing_summary_raw_traits.csv"), row.names = FALSE, fileEncoding = "UTF-8")

hard_outlier_rules <- list(
  BW_birth = c(min = 20, max = 80),
  BW_wk1 = c(min = 30, max = 250),
  BW_wk2 = c(min = 50, max = 500),
  BW_wk3 = c(min = 80, max = 800),
  BW_wk5 = c(min = 150, max = 1400),
  BW_wk7 = c(min = 250, max = 2000),
  BW_wk9 = c(min = 350, max = 2600),
  BW_wk11 = c(min = 450, max = 3200),
  BW_wk13 = c(min = 500, max = 3600),
  BW_wk15 = c(min = 600, max = 4200),
  SL15 = c(min = 45, max = 95),
  HW15 = c(min = 35, max = 100),
  SC15 = c(min = 2.0, max = 6.5),
  BSL15 = c(min = 10, max = 40),
  KL15 = c(min = 5, max = 25),
  Live_wt = c(min = 500, max = 4500),
  Carcass_wt = c(min = 300, max = 4000),
  Evis_wt = c(min = 200, max = 3500),
  AF_wt = c(min = 0, max = 300)
)

hard_outlier_table <- purrr::map_dfr(names(hard_outlier_rules), function(v) {
  if (!v %in% names(dat_raw)) return(NULL)
  min_v <- hard_outlier_rules[[v]][["min"]]
  max_v <- hard_outlier_rules[[v]][["max"]]
  dat_raw %>%
    filter(!is.na(.data[[v]])) %>%
    filter(.data[[v]] < min_v | .data[[v]] > max_v) %>%
    transmute(
      Outlier_type = "hard_range",
      group, wing_id, Trait = v, Trait_label = label_trait(v), Value = .data[[v]],
      Lower_threshold = min_v, Upper_threshold = max_v,
      is_any_omics, is_muscle_rna, is_cecum_meta, is_duodenum_rna
    )
})
write.csv(hard_outlier_table, file.path(OUTDIR, "QC_hard_outlier_table.csv"), row.names = FALSE, fileEncoding = "UTF-8")

flag_robust_outliers <- function(df, traits, z_cut = 3.5) {
  traits <- traits[traits %in% names(df)]
  if (length(traits) == 0) return(tibble())
  df %>%
    select(group, wing_id, is_any_omics, all_of(traits)) %>%
    pivot_longer(cols = all_of(traits), names_to = "Trait", values_to = "Value") %>%
    filter(!is.na(Value)) %>%
    group_by(group, Trait) %>%
    mutate(
      N_group_trait = n(),
      med = median(Value, na.rm = TRUE),
      mad_value = mad(Value, constant = 1.4826, na.rm = TRUE),
      robust_z = ifelse(!is.na(mad_value) & mad_value > 0, (Value - med) / mad_value, NA_real_),
      is_robust_outlier = !is.na(robust_z) & abs(robust_z) > z_cut
    ) %>%
    ungroup() %>%
    filter(is_robust_outlier) %>%
    transmute(
      Outlier_type = "group_MAD",
      group, wing_id, Trait, Trait_label = label_trait(Trait), Value,
      Group_median = med, Group_MAD = mad_value, Robust_z = robust_z,
      N_group_trait, is_any_omics
    )
}
robust_outlier_table <- flag_robust_outliers(dat_raw, all_raw_numeric)
write.csv(robust_outlier_table, file.path(OUTDIR, "QC_robust_outlier_table.csv"), row.names = FALSE, fileEncoding = "UTF-8")

dat_qc <- dat_raw
if (APPLY_HARD_OUTLIERS_TO_NA && nrow(hard_outlier_table) > 0) {
  for (i in seq_len(nrow(hard_outlier_table))) {
    idx <- which(dat_qc$group == hard_outlier_table$group[i] & dat_qc$wing_id == hard_outlier_table$wing_id[i])
    tr <- hard_outlier_table$Trait[i]
    if (length(idx) > 0 && tr %in% names(dat_qc)) dat_qc[[tr]][idx] <- NA_real_
  }
}
if (APPLY_ROBUST_OUTLIERS_TO_NA && nrow(robust_outlier_table) > 0) {
  for (i in seq_len(nrow(robust_outlier_table))) {
    idx <- which(dat_qc$group == robust_outlier_table$group[i] & dat_qc$wing_id == robust_outlier_table$wing_id[i])
    tr <- robust_outlier_table$Trait[i]
    if (length(idx) > 0 && tr %in% names(dat_qc)) dat_qc[[tr]][idx] <- NA_real_
  }
}

# Use only body-weight columns present after standardization/QC.
bw_order <- bw_order[bw_order %in% names(dat_qc)]
if (length(bw_order) < 2) {
  stop("Fewer than two body-weight columns are available for longitudinal growth QC.")
}

bw_long_qc <- dat_qc %>%
  select(group, wing_id, is_any_omics, all_of(bw_order)) %>%
  pivot_longer(cols = all_of(bw_order), names_to = "Trait", values_to = "BW") %>%
  mutate(Week = case_when(Trait == "BW_birth" ~ 0, TRUE ~ as.numeric(stringr::str_extract(Trait, "\\d+")))) %>%
  filter(!is.na(BW)) %>%
  arrange(group, wing_id, Week)

bw_anomaly_table <- bw_long_qc %>%
  group_by(group, wing_id) %>%
  arrange(Week, .by_group = TRUE) %>%
  mutate(
    Prev_week = lag(Week),
    Prev_BW = lag(BW),
    Delta_BW = BW - Prev_BW,
    Interval_weeks = Week - Prev_week,
    ADG_interval = Delta_BW / (Interval_weeks * 7),
    is_negative_growth = !is.na(Delta_BW) & Delta_BW < 0
  ) %>%
  ungroup() %>%
  filter(is_negative_growth) %>%
  transmute(group, wing_id, From_week = Prev_week, To_week = Week, Prev_BW, BW, Delta_BW, ADG_interval, is_any_omics)
write.csv(bw_anomaly_table, file.path(OUTDIR, "QC_negative_growth_intervals.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# Remove individuals with longitudinally inconsistent body weight records.
# The rule is individual-level rather than interval-level: if any available later
# body-weight time point is lower than its previous available time point, that bird
# is excluded from downstream calculations. This is conservative and avoids using
# biologically implausible growth trajectories in heterosis estimates.
negative_growth_individuals_removed <- bw_anomaly_table %>%
  group_by(group, wing_id) %>%
  summarise(
    N_negative_intervals = n(),
    Negative_intervals = paste0(From_week, "-", To_week, " wk", collapse = "; "),
    Delta_BW_values = paste(round(Delta_BW, 3), collapse = "; "),
    is_any_omics = any(is_any_omics, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    dat_qc %>%
      select(
        group, wing_id,
        any_of(c("leg_id", "aluminum_leg_id", "muscle_rna_id", "cecum_meta_id", "duodenum_rna_id")),
        any_of(bw_order)
      ),
    by = c("group", "wing_id")
  ) %>%
  arrange(group, wing_id)

write.csv(
  negative_growth_individuals_removed,
  file.path(OUTDIR, "QC_negative_growth_individuals_removed.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

negative_growth_removal_summary <- negative_growth_individuals_removed %>%
  count(group, name = "N_removed") %>%
  complete(group = factor(GROUP_LEVELS, levels = GROUP_LEVELS), fill = list(N_removed = 0)) %>%
  mutate(group = as.character(group))
write.csv(
  negative_growth_removal_summary,
  file.path(OUTDIR, "QC_negative_growth_individuals_removed_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

if (FILTER_NEGATIVE_GROWTH_INDIVIDUALS && nrow(negative_growth_individuals_removed) > 0) {
  dat_qc <- dat_qc %>%
    anti_join(
      negative_growth_individuals_removed %>% select(group, wing_id),
      by = c("group", "wing_id")
    )
}


body_jump_table <- dat_qc %>%
  mutate(
    d_HW_13_15 = ifelse(!is.na(HW13) & !is.na(HW15), HW15 - HW13, NA_real_),
    d_SL_13_15 = ifelse(!is.na(SL13) & !is.na(SL15), SL15 - SL13, NA_real_),
    d_SC_13_15 = ifelse(!is.na(SC13) & !is.na(SC15), SC15 - SC13, NA_real_),
    d_BSL_13_15 = ifelse(!is.na(BSL13) & !is.na(BSL15), BSL15 - BSL13, NA_real_),
    d_KL_13_15 = ifelse(!is.na(KL13) & !is.na(KL15), KL15 - KL13, NA_real_)
  ) %>%
  select(group, wing_id, is_any_omics, starts_with("d_"), HW13, HW15, SL13, SL15, SC13, SC15, BSL13, BSL15, KL13, KL15) %>%
  pivot_longer(cols = starts_with("d_"), names_to = "Jump_trait", values_to = "Jump_value") %>%
  filter(!is.na(Jump_value)) %>%
  group_by(group, Jump_trait) %>%
  mutate(
    med = median(Jump_value, na.rm = TRUE),
    mad_value = mad(Jump_value, constant = 1.4826, na.rm = TRUE),
    robust_z = ifelse(!is.na(mad_value) & mad_value > 0, (Jump_value - med) / mad_value, NA_real_)
  ) %>%
  ungroup() %>%
  filter(!is.na(robust_z) & abs(robust_z) > 3.5) %>%
  arrange(desc(abs(robust_z)))
write.csv(body_jump_table, file.path(OUTDIR, "QC_body_measurement_jump_outliers.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# =============================================================================
# 6. Derive traits after QC
# =============================================================================

derive_traits <- function(df) {
  df %>%
    mean_if_available(c("PH1", "PH2", "PH3"), "PH_mean") %>%
    mean_if_available(c("L1", "L2", "L3"), "L_mean") %>%
    mean_if_available(c("a1", "a2", "a3"), "a_mean") %>%
    mean_if_available(c("b1", "b2", "b3"), "b_mean") %>%
    mean_if_available(c("Shear1", "Shear2", "Shear3"), "Shear_mean") %>%
    mutate(
      ADG_0_1 = safe_ratio(BW_wk1 - BW_birth, 7, 1),
      ADG_1_2 = safe_ratio(BW_wk2 - BW_wk1, 7, 1),
      ADG_2_3 = safe_ratio(BW_wk3 - BW_wk2, 7, 1),
      ADG_3_5 = safe_ratio(BW_wk5 - BW_wk3, 14, 1),
      ADG_5_7 = safe_ratio(BW_wk7 - BW_wk5, 14, 1),
      ADG_7_9 = safe_ratio(BW_wk9 - BW_wk7, 14, 1),
      ADG_9_11 = safe_ratio(BW_wk11 - BW_wk9, 14, 1),
      ADG_11_13 = safe_ratio(BW_wk13 - BW_wk11, 14, 1),
      ADG_13_15 = safe_ratio(BW_wk15 - BW_wk13, 14, 1),
      ADG_0_15 = safe_ratio(BW_wk15 - BW_birth, 105, 1),
      Carcass_live_pct = safe_ratio(Carcass_wt, Live_wt),
      Semi_evis_live_pct = safe_ratio(Semi_evis_wt, Live_wt),
      Evis_live_pct = safe_ratio(Evis_wt, Live_wt),
      Leg_live_pct = safe_ratio(Leg_wt, Live_wt),
      Breast_live_pct = safe_ratio(Breast_wt, Live_wt),
      Breast_evis_pct = safe_ratio(Breast_wt, Evis_wt),
      Wing_live_pct = safe_ratio(Wing_wt, Live_wt),
      Head_live_pct = safe_ratio(Head_wt, Live_wt),
      AF_live_pct = safe_ratio(AF_wt, Live_wt),
      Paw_live_pct = safe_ratio(Paw_wt, Live_wt),
      is_slaughter = if_any(any_of(c("Live_wt", "Carcass_wt", "Semi_evis_wt", "Evis_wt", "AF_wt")), ~ !is.na(.x))
    )
}

dat <- derive_traits(dat_qc)

adg_vars <- c("ADG_0_1", "ADG_1_2", "ADG_2_3", "ADG_3_5", "ADG_5_7", "ADG_7_9", "ADG_9_11", "ADG_11_13", "ADG_13_15", "ADG_0_15")
if (SET_NEGATIVE_ADG_TO_NA) {
  adg_negative_table <- dat %>%
    select(group, wing_id, is_any_omics, all_of(adg_vars), all_of(bw_vars)) %>%
    pivot_longer(cols = all_of(adg_vars), names_to = "Trait", values_to = "ADG") %>%
    filter(!is.na(ADG), ADG < 0) %>%
    arrange(group, wing_id, Trait)
  write.csv(adg_negative_table, file.path(OUTDIR, "QC_negative_ADG_set_to_NA.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  for (v in adg_vars) dat[[v]][!is.na(dat[[v]]) & dat[[v]] < 0] <- NA_real_
}

# =============================================================================
# 7. Statistics and heterosis
# =============================================================================

# Article 1 analysis scope.
# Reserved Article-2 traits are still readable from the source workbook, but they
# are excluded from Article-1 statistical families, heterosis summaries and figures.
article1_slaughter_vars <- c(
  "Live_wt", "Carcass_wt", "Semi_evis_wt", "Evis_wt", "AF_wt"
)
article1_slaughter_ratio_vars <- c(
  "Carcass_live_pct", "Semi_evis_live_pct", "Evis_live_pct", "AF_live_pct"
)
body15_vars <- c("HW15", "SL15", "SC15", "BSL15", "KL15")

analysis_traits <- unique(c(
  bw_vars, adg_vars, body15_vars,
  article1_slaughter_vars, article1_slaughter_ratio_vars
))
analysis_traits <- analysis_traits[
  analysis_traits %in% names(dat) &
    !analysis_traits %in% ARTICLE1_RESERVED_FOR_ARTICLE2
]

main_traits <- ARTICLE1_CORE_TRAITS[ARTICLE1_CORE_TRAITS %in% analysis_traits]

sample_counts <- dat %>%
  group_by(group) %>%
  summarise(
    N_total = n(),
    N_slaughter = sum(is_slaughter, na.rm = TRUE),
    N_cecum_metagenome = sum(is_cecum_meta, na.rm = TRUE),
    N_duodenum_RNAseq = sum(is_duodenum_rna, na.rm = TRUE),
    N_any_article1_omics_ID = sum(is_any_omics, na.rm = TRUE),
    .groups = "drop"
  )

trait_summary <- purrr::map_dfr(analysis_traits, function(v) {
  dat %>%
    group_by(group) %>%
    summarise(
      Trait = v,
      Trait_label = label_trait(v),
      N = sum(!is.na(.data[[v]])),
      Mean = mean(.data[[v]], na.rm = TRUE),
      SD = sd(.data[[v]], na.rm = TRUE),
      SE = se(.data[[v]]),
      Median = median(.data[[v]], na.rm = TRUE),
      Min = min(.data[[v]], na.rm = TRUE),
      Max = max(.data[[v]], na.rm = TRUE),
      .groups = "drop"
    )
}) %>% mutate(across(c(Mean, SD, SE, Median, Min, Max), ~ ifelse(is.infinite(.x), NA_real_, .x)))

global_trait_test <- function(df, trait) {
  sub <- df %>% select(group, value = all_of(trait)) %>% filter(!is.na(group), !is.na(value))
  if (nrow(sub) < 3 || n_distinct(sub$group) < 2) return(tibble(Trait = trait, ANOVA_P = NA_real_, Kruskal_P = NA_real_, Eta2 = NA_real_))
  aov_p <- tryCatch(summary(aov(value ~ group, data = sub))[[1]][["Pr(>F)"]][1], error = function(e) NA_real_)
  kw_p <- tryCatch(kruskal.test(value ~ group, data = sub)$p.value, error = function(e) NA_real_)
  eta2 <- tryCatch({
    fit <- aov(value ~ group, data = sub)
    ss <- summary(fit)[[1]][["Sum Sq"]]
    ss[1] / sum(ss, na.rm = TRUE)
  }, error = function(e) NA_real_)
  tibble(Trait = trait, ANOVA_P = aov_p, Kruskal_P = kw_p, Eta2 = eta2)
}

pairwise_trait_test <- function(df, trait) {
  sub <- df %>% select(group, value = all_of(trait)) %>% filter(!is.na(group), !is.na(value))
  if (nrow(sub) < 3 || n_distinct(sub$group) < 2) return(NULL)
  pw <- tryCatch(pairwise.wilcox.test(sub$value, sub$group, p.adjust.method = "BH", exact = FALSE)$p.value, error = function(e) NULL)
  if (is.null(pw)) return(NULL)
  as.data.frame(as.table(pw), stringsAsFactors = FALSE) %>%
    filter(!is.na(Freq)) %>%
    transmute(Trait = trait, Group1 = Var1, Group2 = Var2, Wilcoxon_BH_P = Freq)
}

global_tests <- purrr::map_dfr(analysis_traits, ~ global_trait_test(dat, .x)) %>%
  mutate(
    ANOVA_FDR = p.adjust(ANOVA_P, method = "BH"),
    Kruskal_FDR = p.adjust(Kruskal_P, method = "BH"),
    Trait_label = label_trait(Trait)
  )
pairwise_tests <- purrr::map_dfr(analysis_traits, ~ pairwise_trait_test(dat, .x)) %>%
  mutate(Trait_label = label_trait(Trait))

classify_position <- function(yy, yx, xx, tol = 0.15) {
  if (any(is.na(c(yy, yx, xx)))) return(NA_character_)
  parent_diff <- abs(yy - xx)
  if (parent_diff == 0) {
    if (abs(yx - yy) <= abs(yy) * 0.05 + 1e-8) return("mid-parent-like")
    if (yx > yy) return("above_parent")
    if (yx < yy) return("below_parent")
    return("mid-parent-like")
  }
  lo <- min(yy, xx)
  hi <- max(yy, xx)
  pos <- (yx - lo) / (hi - lo)
  if (pos < -tol) return("below_lower_parent")
  if (pos > 1 + tol) return("above_higher_parent")
  if (abs(pos - 0.5) <= tol) return("mid-parent-like")
  if (pos <= tol) return(ifelse(yy < xx, "YY-like", "XX-like"))
  if (pos >= 1 - tol) return(ifelse(yy > xx, "YY-like", "XX-like"))
  return("intermediate")
}

calc_metrics_from_means <- function(yy, yx, xx, lower_is_better = FALSE) {
  mp <- (yy + xx) / 2
  high_parent <- max(yy, xx)
  low_parent <- min(yy, xx)
  high_parent_name <- ifelse(yy >= xx, "YY", "XX")
  low_parent_name <- ifelse(yy <= xx, "YY", "XX")
  favorable_parent <- if (lower_is_better) low_parent else high_parent
  favorable_parent_name <- if (lower_is_better) low_parent_name else high_parent_name
  mph <- ifelse(!is.na(mp) && mp != 0, (yx - mp) / mp * 100, NA_real_)
  bph_high <- ifelse(!is.na(high_parent) && high_parent != 0, (yx - high_parent) / high_parent * 100, NA_real_)
  directional <- if (lower_is_better) {
    ifelse(!is.na(favorable_parent) && favorable_parent != 0, (favorable_parent - yx) / favorable_parent * 100, NA_real_)
  } else {
    ifelse(!is.na(favorable_parent) && favorable_parent != 0, (yx - favorable_parent) / favorable_parent * 100, NA_real_)
  }
  f1_position <- ifelse(!is.na(high_parent - low_parent) && (high_parent - low_parent) != 0,
                        (yx - low_parent) / (high_parent - low_parent), NA_real_)
  tibble(
    MP = mp,
    High_parent = high_parent_name,
    Low_parent = low_parent_name,
    Favorable_parent = favorable_parent_name,
    MPH_percent = mph,
    BPH_high_parent_percent = bph_high,
    Directional_BPH_percent = directional,
    F1_position = f1_position,
    F1_pattern = classify_position(yy, yx, xx)
  )
}

heterosis_one_trait <- function(df, trait, boot_n = BOOT_N) {
  sub <- df %>% select(group, value = all_of(trait)) %>% filter(!is.na(group), !is.na(value))
  yy <- sub$value[sub$group == "YY"]
  yx <- sub$value[sub$group == "YX"]
  xx <- sub$value[sub$group == "XX"]
  lower_is_better <- trait %in% LOWER_IS_BETTER_TRAITS

  if (length(yy) < 2 || length(yx) < 2 || length(xx) < 2) {
    return(tibble(
      Trait = trait, Trait_label = label_trait(trait),
      N_YY = length(yy), N_YX = length(yx), N_XX = length(xx),
      Mean_YY = mean(yy, na.rm = TRUE), Mean_YX = mean(yx, na.rm = TRUE), Mean_XX = mean(xx, na.rm = TRUE),
      MP = NA_real_, High_parent = NA_character_, Low_parent = NA_character_, Favorable_parent = NA_character_,
      MPH_percent = NA_real_, MPH_CI_low = NA_real_, MPH_CI_high = NA_real_,
      BPH_high_parent_percent = NA_real_, BPH_high_CI_low = NA_real_, BPH_high_CI_high = NA_real_,
      Directional_BPH_percent = NA_real_, Directional_BPH_CI_low = NA_real_, Directional_BPH_CI_high = NA_real_,
      F1_position = NA_real_, F1_position_CI_low = NA_real_, F1_position_CI_high = NA_real_,
      F1_pattern = NA_character_, Lower_is_better = lower_is_better
    ))
  }

  means <- c(YY = mean(yy), YX = mean(yx), XX = mean(xx))
  metrics <- calc_metrics_from_means(means[["YY"]], means[["YX"]], means[["XX"]], lower_is_better)

  boot_df <- replicate(boot_n, {
    bm <- c(
      YY = mean(sample(yy, length(yy), replace = TRUE), na.rm = TRUE),
      YX = mean(sample(yx, length(yx), replace = TRUE), na.rm = TRUE),
      XX = mean(sample(xx, length(xx), replace = TRUE), na.rm = TRUE)
    )
    as.numeric(calc_metrics_from_means(bm[["YY"]], bm[["YX"]], bm[["XX"]], lower_is_better) %>%
                 select(MPH_percent, BPH_high_parent_percent, Directional_BPH_percent, F1_position))
  })
  rownames(boot_df) <- c("MPH", "BPH_high", "Directional_BPH", "F1_position")

  tibble(
    Trait = trait,
    Trait_label = label_trait(trait),
    N_YY = length(yy), N_YX = length(yx), N_XX = length(xx),
    Mean_YY = means[["YY"]], Mean_YX = means[["YX"]], Mean_XX = means[["XX"]],
    MP = metrics$MP,
    High_parent = metrics$High_parent,
    Low_parent = metrics$Low_parent,
    Favorable_parent = metrics$Favorable_parent,
    MPH_percent = metrics$MPH_percent,
    MPH_CI_low = quantile(boot_df["MPH", ], 0.025, na.rm = TRUE),
    MPH_CI_high = quantile(boot_df["MPH", ], 0.975, na.rm = TRUE),
    BPH_high_parent_percent = metrics$BPH_high_parent_percent,
    BPH_high_CI_low = quantile(boot_df["BPH_high", ], 0.025, na.rm = TRUE),
    BPH_high_CI_high = quantile(boot_df["BPH_high", ], 0.975, na.rm = TRUE),
    Directional_BPH_percent = metrics$Directional_BPH_percent,
    Directional_BPH_CI_low = quantile(boot_df["Directional_BPH", ], 0.025, na.rm = TRUE),
    Directional_BPH_CI_high = quantile(boot_df["Directional_BPH", ], 0.975, na.rm = TRUE),
    F1_position = metrics$F1_position,
    F1_position_CI_low = quantile(boot_df["F1_position", ], 0.025, na.rm = TRUE),
    F1_position_CI_high = quantile(boot_df["F1_position", ], 0.975, na.rm = TRUE),
    F1_pattern = metrics$F1_pattern,
    Lower_is_better = lower_is_better
  )
}

heterosis_table <- purrr::map_dfr(analysis_traits, ~ heterosis_one_trait(dat, .x)) %>%
  left_join(global_tests %>% select(Trait, ANOVA_P, ANOVA_FDR, Kruskal_P, Kruskal_FDR, Eta2), by = "Trait") %>%
  mutate(Trait_label = label_trait(Trait))

# Omics-sample representativeness.
representativeness_one <- function(df, platform_col, platform_label, trait) {
  purrr::map_dfr(GROUP_LEVELS, function(g) {
    sub <- df %>% filter(group == g) %>% select(group, wing_id, is_seq = all_of(platform_col), value = all_of(trait)) %>% filter(!is.na(value))
    n_total <- nrow(sub)
    n_seq <- sum(sub$is_seq, na.rm = TRUE)
    n_nonseq <- sum(!sub$is_seq, na.rm = TRUE)
    mean_total <- ifelse(n_total > 0, mean(sub$value, na.rm = TRUE), NA_real_)
    mean_seq <- ifelse(n_seq > 0, mean(sub$value[sub$is_seq], na.rm = TRUE), NA_real_)
    mean_nonseq <- ifelse(n_nonseq > 0, mean(sub$value[!sub$is_seq], na.rm = TRUE), NA_real_)
    sd_total <- ifelse(n_total > 1, sd(sub$value, na.rm = TRUE), NA_real_)
    rel_dev <- ifelse(!is.na(mean_total) && mean_total != 0, (mean_seq - mean_total) / mean_total * 100, NA_real_)
    z_dev <- ifelse(!is.na(sd_total) && sd_total != 0, (mean_seq - mean_total) / sd_total, NA_real_)
    wilcox_p <- if (n_seq >= 2 && n_nonseq >= 3) {
      tryCatch(wilcox.test(sub$value[sub$is_seq], sub$value[!sub$is_seq], exact = FALSE)$p.value, error = function(e) NA_real_)
    } else NA_real_
    result <- case_when(
      is.na(rel_dev) | is.na(z_dev) ~ "Insufficient data",
      abs(rel_dev) <= 5 & abs(z_dev) <= 0.5 ~ "Highly representative",
      abs(rel_dev) <= 10 & abs(z_dev) <= 1.0 ~ "Acceptable",
      TRUE ~ "Potentially biased"
    )
    tibble(
      Platform = platform_label, Trait = trait, Trait_label = label_trait(trait), Group = g,
      N_total = n_total, N_seq = n_seq, N_nonseq = n_nonseq,
      Mean_total = mean_total, Mean_seq = mean_seq, Mean_nonseq = mean_nonseq,
      Relative_deviation_percent = rel_dev, Z_from_group_mean = z_dev,
      Wilcoxon_P = wilcox_p, Result = result
    )
  })
}

platforms <- tribble(
  ~platform_col, ~platform_label,
  "is_any_omics", "Any Article-1 intestinal omics ID",
  "is_cecum_meta", "Cecal metagenome",
  "is_duodenum_rna", "Duodenum transcriptome"
)
rep_traits <- c(
  "BW_wk15", "ADG_0_15",
  "Live_wt", "Carcass_wt", "Evis_wt",
  "AF_wt", "AF_live_pct"
)
rep_traits <- rep_traits[rep_traits %in% names(dat)]
representativeness_table <- purrr::map_dfr(seq_len(nrow(platforms)), function(i) {
  purrr::map_dfr(rep_traits, ~ representativeness_one(dat, platforms$platform_col[i], platforms$platform_label[i], .x))
})

# =============================================================================
# 8. Export tables
# =============================================================================

article1_export_cols <- unique(c(
  "group", "group_raw", "wing_id", "leg_id", "aluminum_leg_id", "sex_corrected",
  "cecum_meta_id", "duodenum_rna_id",
  "is_cecum_meta", "is_duodenum_rna", "is_any_omics", "is_slaughter",
  bw_vars, adg_vars, body_vars,
  "Live_wt", "Carcass_wt", "Semi_evis_wt", "Evis_wt",
  "Carcass_live_pct", "Semi_evis_live_pct", "Evis_live_pct",
  "AF_wt", "AF_live_pct"
))
article1_export_cols <- article1_export_cols[article1_export_cols %in% names(dat)]
dat_article1 <- dat %>% select(all_of(article1_export_cols))

omics_sample_key <- dat_article1 %>%
  filter(is_any_omics) %>%
  select(
    group, wing_id, leg_id, aluminum_leg_id, sex_corrected,
    cecum_meta_id, duodenum_rna_id,
    is_cecum_meta, is_duodenum_rna,
    all_of(main_traits)
  )

# Canonical Article-1 downstream phenotype table. Breast/leg yield and meat
# quality are intentionally excluded so downstream scripts cannot select them.
write.csv(dat_article1, file.path(OUTDIR, "Table0_clean_phenotype_with_derived_traits.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(omics_sample_key, file.path(OUTDIR, "Table0_omics_sample_key_for_integration.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(sample_counts, file.path(OUTDIR, "Table1_sample_counts_by_group.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(trait_summary, file.path(OUTDIR, "Table2_trait_summary_by_group.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(global_tests, file.path(OUTDIR, "Table3_global_group_tests.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(pairwise_tests, file.path(OUTDIR, "Table4_pairwise_wilcoxon_tests.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(heterosis_table, file.path(OUTDIR, "Table5_heterosis_and_F1_relative_position.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(representativeness_table, file.path(OUTDIR, "Table6_omics_sample_representativeness.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# =============================================================================
# 9. Figures for manuscript and supplement
# =============================================================================

# Fig1A: growth trajectory.
growth_long <- dat %>%
  select(group, wing_id, is_any_omics, all_of(bw_vars)) %>%
  pivot_longer(cols = all_of(bw_vars), names_to = "Trait", values_to = "BW") %>%
  mutate(Week = case_when(Trait == "BW_birth" ~ 0, TRUE ~ as.numeric(stringr::str_extract(Trait, "\\d+")))) %>%
  filter(!is.na(BW), !is.na(Week))

growth_summary <- growth_long %>%
  group_by(group, Week) %>%
  summarise(Mean = mean(BW), SE = se(BW), N = n(), .groups = "drop")

p_growth <- ggplot(growth_long, aes(x = Week, y = BW, group = wing_id, color = group)) +
  geom_line(alpha = 0.10, linewidth = 0.32) +
  geom_point(alpha = 0.16, size = 0.7) +
  geom_line(data = growth_summary, aes(x = Week, y = Mean, group = group, color = group), linewidth = 1.25, inherit.aes = FALSE) +
  geom_point(data = growth_summary, aes(x = Week, y = Mean, color = group), size = 2.3, inherit.aes = FALSE) +
  geom_errorbar(data = growth_summary, aes(x = Week, ymin = Mean - SE, ymax = Mean + SE, color = group), width = 0.25, linewidth = 0.55, inherit.aes = FALSE) +
  scale_color_manual(values = GROUP_COLORS) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 5, 7, 9, 11, 13, 15)) +
  labs(
    title = "Growth trajectories of YY, YX F1 and XX hens",
    x = "Age (weeks)",
    y = "Body weight (g)",
    color = "Group"
  ) +
  theme_manuscript(12)
save_plot(p_growth, "Fig1A_growth_trajectory", 8.2, 5.4)

# Fig1B: core growth traits instead of all ADG intervals.
growth_box_traits <- c("BW_wk15", "ADG_0_15")
growth_box_traits <- growth_box_traits[growth_box_traits %in% names(dat)]
growth_box_long <- dat %>%
  select(group, wing_id, is_any_omics, all_of(growth_box_traits)) %>%
  pivot_longer(cols = all_of(growth_box_traits), names_to = "Trait", values_to = "Value") %>%
  filter(!is.na(Value)) %>%
  mutate(Trait_label = factor(label_trait(Trait), levels = label_trait(growth_box_traits)))

p_growth_box <- ggplot(growth_box_long, aes(x = group, y = Value, fill = group)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.72) +
  geom_point(
    aes(color = group, size = is_any_omics),
    position = position_jitter(width = 0.12, height = 0),
    shape = 21,
    fill = "white",
    stroke = 0.9,
    alpha = 0.92
  ) +
  facet_wrap(~ Trait_label, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = GROUP_COLORS, guide = "none") +
  scale_color_manual(values = GROUP_COLORS, guide = "none") +
  scale_size_manual(values = c("FALSE" = 1.35, "TRUE" = 2.5), name = "Any omics") +
  labs(title = "Core growth phenotypes", x = NULL, y = "Trait value") +
  theme_manuscript(11)
save_plot(p_growth_box, "Fig1B_core_growth_boxplots", 6.8, 4.3)

# Fig1C: core slaughter traits and carcass composition.
core_carcass_traits <- c("Live_wt", "Carcass_wt", "Evis_wt")
core_carcass_traits <- core_carcass_traits[core_carcass_traits %in% names(dat)]
carcass_long <- dat %>%
  select(group, wing_id, is_any_omics, all_of(core_carcass_traits)) %>%
  pivot_longer(cols = all_of(core_carcass_traits), names_to = "Trait", values_to = "Value") %>%
  filter(!is.na(Value)) %>%
  mutate(Trait_label = factor(label_trait(Trait), levels = label_trait(core_carcass_traits)))

p_carcass <- ggplot(carcass_long, aes(x = group, y = Value, fill = group)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.72) +
  geom_point(
    aes(color = group, size = is_any_omics),
    position = position_jitter(width = 0.12, height = 0),
    shape = 21,
    fill = "white",
    stroke = 0.9,
    alpha = 0.92
  ) +
  facet_wrap(~ Trait_label, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = GROUP_COLORS, guide = "none") +
  scale_color_manual(values = GROUP_COLORS, guide = "none") +
  scale_size_manual(values = c("FALSE" = 1.35, "TRUE" = 2.5), name = "Any omics") +
  labs(
    title = "Core overall carcass-weight traits",
    subtitle = NULL,
    x = NULL,
    y = "Trait value"
  ) +
  theme_manuscript(10)
save_plot(p_carcass, "Fig1C_core_carcass_boxplots", 11, 5.4)

# Fig1D: F1 relative position, replacing one-trait-per-bar pattern plot.
position_traits <- c("BW_wk15", "ADG_0_15", "Live_wt", "Carcass_wt", "Evis_wt")
if (INCLUDE_BODY_MEASUREMENTS_IN_MAIN) {
  position_traits <- c(position_traits, "HW15", "SL15", "SC15", "BSL15", "KL15")
}
position_traits <- position_traits[position_traits %in% heterosis_table$Trait]
pos_df <- heterosis_table %>%
  filter(Trait %in% position_traits, !is.na(F1_position)) %>%
  mutate(Trait_label = factor(Trait_label, levels = rev(label_trait(position_traits))))

if (nrow(pos_df) > 0) {
  p_pos <- ggplot(pos_df, aes(x = F1_position, y = Trait_label, color = F1_pattern)) +
    geom_vline(xintercept = c(0, 0.5, 1), linetype = c("dashed", "solid", "dashed"), color = "grey55", linewidth = c(0.45, 0.55, 0.45)) +
    geom_errorbarh(aes(xmin = F1_position_CI_low, xmax = F1_position_CI_high), height = 0.18, alpha = 0.65, linewidth = 0.45) +
    geom_point(size = 2.5, alpha = 0.9) +
    scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"), limits = c(-0.15, 1.15), oob = scales::squish) +
    labs(
      title = "Relative position of the YX F1 mean between parental lines",
      subtitle = "0 = lower parental mean; 0.5 = mid-parent; 1 = higher parental mean; wide CIs are truncated to improve readability",
      x = "F1 relative position",
      y = NULL,
      color = "Pattern"
    ) +
    theme_manuscript(10)
  save_plot(p_pos, "Fig1D_F1_relative_position", 8.2, 6.6)
}

# Supplement: full ADG intervals.
adg_long <- dat %>%
  select(group, wing_id, is_any_omics, all_of(adg_vars)) %>%
  pivot_longer(cols = all_of(adg_vars), names_to = "Trait", values_to = "ADG") %>%
  mutate(Trait_label = factor(label_trait(Trait), levels = label_trait(adg_vars))) %>%
  filter(!is.na(ADG))

if (nrow(adg_long) > 0) {
  p_adg <- ggplot(adg_long, aes(x = group, y = ADG, fill = group)) +
    geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.72) +
    geom_point(
    aes(color = group, size = is_any_omics),
    position = position_jitter(width = 0.12, height = 0),
    shape = 21,
    fill = "white",
    stroke = 0.9,
    alpha = 0.92
  ) +
    facet_wrap(~ Trait_label, scales = "free_y", ncol = 5) +
    scale_fill_manual(values = GROUP_COLORS, guide = "none") +
    scale_color_manual(values = GROUP_COLORS, guide = "none") +
    scale_size_manual(values = c("FALSE" = 1.35, "TRUE" = 2.5), name = "Any omics") +
    labs(title = "Average daily gain across growth intervals", x = NULL, y = "ADG (g/day)") +
    theme_manuscript(10)
  save_plot(p_adg, "FigS1_ADG_interval_boxplots", 13, 7)
}

# Supplement: body measurement trajectories, with units in facet labels.
body_long <- dat %>%
  select(group, wing_id, is_any_omics, all_of(body_vars)) %>%
  pivot_longer(cols = all_of(body_vars), names_to = "Trait", values_to = "Value") %>%
  mutate(
    Week = as.numeric(stringr::str_extract(Trait, "\\d+")),
    Measure = case_when(
      stringr::str_detect(Trait, "^HW") ~ "Hip width (mm)",
      stringr::str_detect(Trait, "^SL") ~ "Shank length (mm)",
      stringr::str_detect(Trait, "^SC") ~ "Shank circumference (cm)",
      stringr::str_detect(Trait, "^BSL") ~ "Body slant length (cm)",
      stringr::str_detect(Trait, "^KL") ~ "Keel length (cm)",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Value), !is.na(Week), Week %in% c(9, 11, 13, 15), !is.na(Measure))

body_summary <- body_long %>%
  group_by(group, Week, Measure) %>%
  summarise(Mean = mean(Value), SE = se(Value), N = n(), .groups = "drop")

if (nrow(body_summary) > 0) {
  p_body <- ggplot(body_summary, aes(x = Week, y = Mean, color = group, group = group)) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 2.2) +
    geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.18, linewidth = 0.5) +
    facet_wrap(~ Measure, scales = "free_y", ncol = 3) +
    scale_color_manual(values = GROUP_COLORS) +
    scale_x_continuous(breaks = c(9, 11, 13, 15)) +
    labs(title = "Body measurement trajectories from week 9 to week 15", x = "Age (weeks)", y = "Mean ± SE", color = "Group") +
    theme_manuscript(11)
  save_plot(p_body, "FigS2_body_measurement_trajectories", 10.5, 6.4)
}

# Supplement: core phenotype PCA.
pca_vars <- c("BW_wk15", "ADG_0_15", "Live_wt", "Carcass_wt", "Evis_wt")
pca_vars <- pca_vars[pca_vars %in% names(dat)]
pca_df <- dat %>% select(group, wing_id, is_any_omics, all_of(pca_vars)) %>% filter(if_all(all_of(pca_vars), ~ !is.na(.x)))

if (nrow(pca_df) >= 6 && length(pca_vars) >= 3) {
  pca_fit <- prcomp(pca_df[, pca_vars], center = TRUE, scale. = TRUE)
  pca_scores <- as.data.frame(pca_fit$x[, 1:2]) %>% bind_cols(pca_df %>% select(group, wing_id, is_any_omics))
  var_exp <- round(summary(pca_fit)$importance[2, 1:2] * 100, 1)
  p_pca <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = group)) +
    stat_ellipse(aes(group = group), linewidth = 0.75, alpha = 0.65, show.legend = FALSE) +
    geom_point(aes(shape = is_any_omics), size = 2.6, alpha = 0.86) +
    scale_color_manual(values = GROUP_COLORS) +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17), name = "Any omics") +
    labs(title = "PCA of core growth and overall carcass-weight phenotypes", x = paste0("PC1 (", var_exp[1], "%)"), y = paste0("PC2 (", var_exp[2], "%)"), color = "Group") +
    theme_manuscript(12)
  save_plot(p_pca, "FigS3_core_phenotype_PCA", 7, 5.8)
}

# Supplement: heterosis heatmap.
het_heat <- heterosis_table %>%
  filter(Trait %in% position_traits) %>%
  select(Trait, Trait_label, MPH_percent, Directional_BPH_percent) %>%
  pivot_longer(cols = c(MPH_percent, Directional_BPH_percent), names_to = "Metric", values_to = "Value") %>%
  mutate(
    Metric = recode(Metric, MPH_percent = "MPH (%)", Directional_BPH_percent = "Directional BPH (%)"),
    Trait_label = factor(Trait_label, levels = rev(label_trait(position_traits)))
  )

if (nrow(het_heat) > 0) {
  p_het <- ggplot(het_heat, aes(x = Metric, y = Trait_label, fill = Value)) +
    geom_tile(color = "white", linewidth = 0.45) +
    geom_text(aes(label = ifelse(is.na(Value), "", sprintf("%.1f", Value))), size = 3.1) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, na.value = "grey90", name = "%") +
    labs(title = "Mid-parent and direction-aware better-parent heterosis", x = NULL, y = NULL) +
    theme_manuscript(10)
  save_plot(p_het, "FigS4_MPH_directional_BPH_heatmap", 6.5, 6.8)
}

# Supplement (optional): Article 1 core-trait correlation.
# Disabled by default. If enabled, it is strictly restricted to Article 1 traits.
if (MAKE_CORE_CORRELATION_HEATMAP) {
  corr_vars <- c(
    "BW_wk15", "ADG_0_15",
    "Live_wt", "Carcass_wt", "Evis_wt",
    "AF_wt", "AF_live_pct"
  )
  corr_vars <- corr_vars[corr_vars %in% names(dat)]

  if (length(corr_vars) >= 3) {
    corr_mat <- cor(
      dat[, corr_vars, drop = FALSE],
      use = "pairwise.complete.obs",
      method = "spearman"
    )

    corr_long <- as.data.frame(as.table(corr_mat), stringsAsFactors = FALSE) %>%
      transmute(
        Trait1 = Var1,
        Trait2 = Var2,
        R = Freq,
        Trait1_label = factor(
          label_trait(Trait1),
          levels = label_trait(corr_vars)
        ),
        Trait2_label = factor(
          label_trait(Trait2),
          levels = rev(label_trait(corr_vars))
        )
      )

    p_corr <- ggplot(
      corr_long,
      aes(x = Trait1_label, y = Trait2_label, fill = R)
    ) +
      geom_tile(color = "white", linewidth = 0.2) +
      scale_fill_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#B2182B",
        midpoint = 0,
        limits = c(-1, 1),
        name = "Spearman r"
      ) +
      labs(
        title = "Spearman correlation among Article 1 core phenotypes",
        x = NULL,
        y = NULL
      ) +
      theme_manuscript(9) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    save_plot(
      p_corr,
      "FigS5_article1_core_trait_correlation_heatmap",
      7.8,
      6.6
    )
  }
}

# Supplement: omics-sample representativeness.
rep_plot <- representativeness_table %>%
  filter(
    Platform %in% c(
      "Any Article-1 intestinal omics ID",
      "Cecal metagenome",
      "Duodenum transcriptome"
    )
  ) %>%
  filter(!is.na(Relative_deviation_percent)) %>%
  mutate(
    Trait_label = factor(Trait_label, levels = label_trait(rep_traits)),
    Platform = factor(
      Platform,
      levels = c(
        "Any Article-1 intestinal omics ID",
        "Cecal metagenome",
        "Duodenum transcriptome"
      )
    )
  )

if (nrow(rep_plot) > 0) {
  p_rep <- ggplot(rep_plot, aes(x = Trait_label, y = Relative_deviation_percent, color = Group, group = Group)) +
    geom_hline(yintercept = 0, color = "grey30", linewidth = 0.5) +
    geom_hline(yintercept = c(-5, 5), color = "grey55", linetype = "dashed", linewidth = 0.35) +
    geom_hline(yintercept = c(-10, 10), color = "grey70", linetype = "dotted", linewidth = 0.35) +
    geom_point(size = 2.1) +
    geom_line(linewidth = 0.65, alpha = 0.75) +
    facet_grid(Platform ~ Group) +
    scale_color_manual(values = GROUP_COLORS, guide = "none") +
    labs(title = "Representativeness of omics individuals within phenotype distributions", x = NULL, y = "Relative deviation (%)") +
    theme_manuscript(9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p_rep, "FigS6_omics_sample_representativeness", 13.5, 8.5)
}

# Supplement: QC summary.
qc_counts <- bind_rows(
  hard_outlier_table %>% count(Outlier_type, name = "N_records"),
  robust_outlier_table %>% count(Outlier_type, name = "N_records"),
  bw_anomaly_table %>% summarise(Outlier_type = "negative_growth_interval", N_records = n()),
  negative_growth_individuals_removed %>% summarise(Outlier_type = "negative_growth_individual_removed", N_records = n()),
  body_jump_table %>% summarise(Outlier_type = "body_measurement_jump", N_records = n())
) %>% filter(!is.na(Outlier_type))

if (nrow(qc_counts) > 0) {
  p_qc <- ggplot(qc_counts, aes(x = reorder(Outlier_type, N_records), y = N_records)) +
    geom_col(width = 0.68) +
    coord_flip() +
    labs(title = "QC flags generated before phenotype story analysis", x = NULL, y = "Number of flagged records") +
    theme_manuscript(10)
  save_plot(p_qc, "FigS7_QC_flag_summary", 6.6, 4.2)
}

# =============================================================================
# 10. Combined Excel workbook
# =============================================================================

trait_dictionary <- tribble(
  ~Trait, ~Original_or_derived, ~Description, ~Unit_or_note,
  "BW_birth", "初生重", "Body weight at hatch", "g",
  "BW_wk1-BW_wk15", "第一周体重-第十五周体重", "Longitudinal body weight", "g",
  "ADG_*", "derived", "Average daily gain between age intervals", "g/day; negative intervals optionally set to NA",
  "HW/SL/SC/BSL/KL", "第9/11/13/15周体尺", "Hip width, shank length, shank circumference, body slant length and keel length", "mm or cm according to original header",
  "Carcass_wt", "屠体重", "Carcass weight", "g",
  "Evis_wt", "全净堂重", "Eviscerated weight", "g",
  "AF_wt", "腹脂重", "Abdominal fat weight", "g",
  "AF_live_pct", "derived", "Abdominal fat weight / live weight", "%",
  "F1_position", "derived", "Relative position of YX between parental means", "0 = lower parent, 0.5 = mid-parent, 1 = higher parent"
)

writexl::write_xlsx(
  list(
    clean_phenotype = dat_article1,
    omics_sample_key = omics_sample_key,
    sample_counts = sample_counts,
    missing_summary_raw = missing_summary_raw,
    hard_outliers = hard_outlier_table,
    robust_outliers = robust_outlier_table,
    negative_growth = bw_anomaly_table,
    negative_growth_removed = negative_growth_individuals_removed,
    negative_growth_removed_summary = negative_growth_removal_summary,
    body_measurement_jumps = body_jump_table,
    trait_summary = trait_summary,
    global_tests = global_tests,
    pairwise_tests = pairwise_tests,
    heterosis = heterosis_table,
    representativeness = representativeness_table,
    trait_dictionary = trait_dictionary
  ),
  path = file.path(OUTDIR, "phenotype_story_tables_v2.xlsx")
)

message("\n============================================================")
message("Phenotype story analysis v2 complete.")
message("Input:  ", COMPLETE_FILE)
message("Output: ", OUTDIR)
message("Missing values: handled trait-by-trait; derived traits become NA when source traits are missing.")
message("Hard outliers applied to NA: ", APPLY_HARD_OUTLIERS_TO_NA)
message("Robust outliers applied to NA: ", APPLY_ROBUST_OUTLIERS_TO_NA)
message("Negative-growth individuals removed: ", FILTER_NEGATIVE_GROWTH_INDIVIDUALS, " (N=", nrow(negative_growth_individuals_removed), ")")
message("Negative interval ADG set to NA after filtering: ", SET_NEGATIVE_ADG_TO_NA)
message("Main figures: Fig1A_growth_trajectory, Fig1B_core_growth_boxplots, Fig1C_core_carcass_boxplots, Fig1D_F1_relative_position")
message("Supplementary figures: FigS1-FigS7")
message("============================================================\n")
