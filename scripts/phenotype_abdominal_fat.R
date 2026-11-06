###############################################################################
# Abdominal fat-focused phenotype analysis for Article 1 - v4.0
#
# Purpose:
#   Generate a standalone abdominal-fat phenotype figure set and tables:
#     A. Abdominal fat weight boxplot
#     B. Abdominal fat percentage boxplot
#     C. Proportion of individuals with measurable abdominal fat
#     D. Abdominal fat weight vs live weight scatter plot + regression line
#     E. Live-weight-adjusted abdominal fat residual plot
#     F. Favorable heterosis of abdominal fat-related traits
#   Additional robustness/sensitivity outputs added in v3:
#     G. Abdominal fat/carcass and abdominal fat/eviscerated percentages
#     H. Carcass-weight- and eviscerated-weight-adjusted abdominal fat residuals
#     I. Positive-only analysis excluding 0-g records
#     J. Sensitivity analysis treating 0-g AF records as below-measurable values
#
# Design notes:
#   YY and XX are contemporaneous purebred progeny; YX is the Y male x X female F1.
#   Abdominal fat is treated as a lower-is-favorable phenotype.
#   AF_wt == 0 is retained as a biologically meaningful "not measurable" record,
#   not as missing, unless the user changes the QC options below.
#
# Input:
#
# Output:
#   Tables and figures in AF_OUTDIR
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
COMPLETE_SHEET <- 1

AF_OUTDIR <- Sys.getenv(
  "AF_OUTDIR",
  unset = ""
)
dir.create(AF_OUTDIR, showWarnings = FALSE, recursive = TRUE)

GROUP_LEVELS <- c("YY", "YX", "XX")
GROUP_COLORS <- c(
  "YY" = "#E64B35",
  "YX" = "#4DBBD5",
  "XX" = "#00A087"
)


# QC options
FILTER_NEGATIVE_GROWTH_INDIVIDUALS <- TRUE
SET_NEGATIVE_AF_TO_NA <- TRUE
SET_NONPOSITIVE_LIVEWT_TO_NA <- TRUE

# Regression/residual options
# TRUE means AF_wt == 0 records are included in the live-weight adjustment model.
# This is recommended if zero means "no measurable abdominal fat" rather than missing.
INCLUDE_ZERO_AF_IN_REGRESSION <- TRUE

# Bootstrap settings for favorable heterosis estimates
BOOT_N <- as.integer(Sys.getenv("AF_BOOT_N", unset = "2000"))
set.seed(20260610)

# =============================================================================
# 1. Packages
# =============================================================================

cran_pkgs <- c(
  "tidyverse", "readxl", "writexl", "ggpubr", "patchwork", "scales"
)

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

`%||%` <- function(a, b) if (!is.null(a)) a else b

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

is_nonempty <- function(x) {
  x <- as.character(x)
  !is.na(x) &
    stringr::str_trim(x) != "" &
    !stringr::str_to_lower(stringr::str_trim(x)) %in% c("na", "nan", "null", "none")
}

make_numeric_if_exists <- function(df, vars) {
  vars <- vars[vars %in% names(df)]
  if (length(vars) == 0) return(df)
  df %>% mutate(across(all_of(vars), safe_numeric))
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
  ggsave(
    file.path(AF_OUTDIR, paste0(filename, ".png")),
    plot_obj,
    width = width,
    height = height,
    dpi = dpi
  )
  ggsave(
    file.path(AF_OUTDIR, paste0(filename, ".pdf")),
    plot_obj,
    width = width,
    height = height
  )
}

# Boxplot style requested by the user:
# - colored transparent box
# - small transparent points for full population
# - white-filled larger points for sequenced/omics individuals
plot_single_box <- function(df, value_col, trait_title, y_lab, subtitle = NULL) {
  box_data <- df %>%
    select(group, wing_id_std, is_sequenced, seq_label, value = all_of(value_col)) %>%
    filter(!is.na(group), !is.na(value)) %>%
    mutate(trait = factor(trait_title, levels = trait_title))

  ggplot(box_data, aes(x = group, y = value)) +
    geom_boxplot(aes(fill = group), alpha = 0.35, outlier.shape = NA, width = 0.65) +
    geom_jitter(
      data = box_data %>% filter(!is_sequenced),
      aes(color = group),
      alpha = 0.16,
      width = 0.18,
      size = 0.8,
      show.legend = FALSE
    ) +
    geom_jitter(
      data = box_data %>% filter(is_sequenced),
      aes(color = group),
      width = 0.11,
      size = 2.6,
      shape = 21,
      fill = "white",
      stroke = 1.1,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = GROUP_COLORS, drop = FALSE, guide = "none") +
    scale_color_manual(values = GROUP_COLORS, drop = FALSE, guide = "none") +
    facet_wrap(~ trait, scales = "free_y", ncol = 1) +
    labs(
      title = trait_title,
      subtitle = subtitle %||% "Small transparent points = population individuals; white-filled points = sequenced individuals",
      x = "Group",
      y = y_lab
    ) +
    ggpubr::theme_pubr(base_size = 11) +
    theme(
      strip.text = element_text(size = 9.5, face = "bold"),
      plot.title = element_text(face = "bold")
    )
}

pairwise_trait_test <- function(df, trait) {
  sub <- df %>%
    select(group, value = all_of(trait)) %>%
    filter(!is.na(group), !is.na(value))

  if (nrow(sub) < 3 || n_distinct(sub$group) < 2) return(NULL)

  pw <- tryCatch(
    pairwise.wilcox.test(sub$value, sub$group, p.adjust.method = "BH", exact = FALSE)$p.value,
    error = function(e) NULL
  )
  if (is.null(pw)) return(NULL)

  as.data.frame(as.table(pw), stringsAsFactors = FALSE) %>%
    filter(!is.na(Freq)) %>%
    transmute(Trait = trait, Group1 = Var1, Group2 = Var2, Wilcoxon_BH_P = Freq)
}

summarise_trait_by_group <- function(df, trait, label) {
  df %>%
    group_by(group) %>%
    summarise(
      Trait = trait,
      Trait_label = label,
      N = sum(!is.na(.data[[trait]])),
      Mean = ifelse(N > 0, mean(.data[[trait]], na.rm = TRUE), NA_real_),
      SD = ifelse(N > 1, sd(.data[[trait]], na.rm = TRUE), NA_real_),
      SE = se(.data[[trait]]),
      Median = ifelse(N > 0, median(.data[[trait]], na.rm = TRUE), NA_real_),
      Min = ifelse(N > 0, min(.data[[trait]], na.rm = TRUE), NA_real_),
      Max = ifelse(N > 0, max(.data[[trait]], na.rm = TRUE), NA_real_),
      .groups = "drop"
    )
}

favorable_heterosis_one_trait <- function(df, trait, label, unit = c("%", "g"), boot_n = BOOT_N) {
  unit <- match.arg(unit)

  sub <- df %>%
    select(group, value = all_of(trait)) %>%
    filter(!is.na(group), !is.na(value))

  vals <- split(sub$value, sub$group)
  yy <- vals[["YY"]]
  yx <- vals[["YX"]]
  xx <- vals[["XX"]]

  empty_result <- tibble(
    Trait = trait,
    Trait_label = label,
    Unit = unit,
    N_YY = length(yy), N_YX = length(yx), N_XX = length(xx),
    Mean_YY = mean(yy, na.rm = TRUE),
    Mean_YX = mean(yx, na.rm = TRUE),
    Mean_XX = mean(xx, na.rm = TRUE),
    Favorable_parent = NA_character_,
    Favorable_MPH = NA_real_, Favorable_MPH_CI_low = NA_real_, Favorable_MPH_CI_high = NA_real_,
    Favorable_BPH = NA_real_, Favorable_BPH_CI_low = NA_real_, Favorable_BPH_CI_high = NA_real_
  )

  if (length(yy) < 2 || length(yx) < 2 || length(xx) < 2) return(empty_result)

  calc <- function(a, b, c) {
    yy_m <- mean(a)
    yx_m <- mean(b)
    xx_m <- mean(c)
    mp <- (yy_m + xx_m) / 2
    low_parent <- min(yy_m, xx_m)

    if (unit == "%") {
      # For lower-is-better positive-scale traits: positive values mean F1 is lower/favorable.
      fav_mph <- ifelse(!is.na(mp) && mp != 0, (mp - yx_m) / mp * 100, NA_real_)
      fav_bph <- ifelse(!is.na(low_parent) && low_parent != 0, (low_parent - yx_m) / low_parent * 100, NA_real_)
    } else {
      # For residuals that may cross zero, use raw units instead of percentage.
      fav_mph <- mp - yx_m
      fav_bph <- low_parent - yx_m
    }
    c(Favorable_MPH = fav_mph, Favorable_BPH = fav_bph)
  }

  obs <- calc(yy, yx, xx)
  boot <- replicate(boot_n, {
    calc(
      sample(yy, length(yy), replace = TRUE),
      sample(yx, length(yx), replace = TRUE),
      sample(xx, length(xx), replace = TRUE)
    )
  })

  yy_mean <- mean(yy)
  yx_mean <- mean(yx)
  xx_mean <- mean(xx)
  fav_parent <- ifelse(yy_mean <= xx_mean, "YY", "XX")

  tibble(
    Trait = trait,
    Trait_label = label,
    Unit = unit,
    N_YY = length(yy), N_YX = length(yx), N_XX = length(xx),
    Mean_YY = yy_mean, Mean_YX = yx_mean, Mean_XX = xx_mean,
    Favorable_parent = fav_parent,
    Favorable_MPH = obs[["Favorable_MPH"]],
    Favorable_MPH_CI_low = quantile(boot["Favorable_MPH", ], 0.025, na.rm = TRUE),
    Favorable_MPH_CI_high = quantile(boot["Favorable_MPH", ], 0.975, na.rm = TRUE),
    Favorable_BPH = obs[["Favorable_BPH"]],
    Favorable_BPH_CI_low = quantile(boot["Favorable_BPH", ], 0.025, na.rm = TRUE),
    Favorable_BPH_CI_high = quantile(boot["Favorable_BPH", ], 0.975, na.rm = TRUE)
  )
}

# =============================================================================
# 3. Read and standardize data
# =============================================================================

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
  "活重" = "Live_wt",
  "腹脂重" = "AF_wt",
  "屠体重" = "Carcass_wt",
  "胴体重" = "Carcass_wt",
  "半净膛重" = "Semi_evis_wt",
  "半净堂重" = "Semi_evis_wt",
  "全净膛重" = "Evis_wt",
  "全净堂重" = "Evis_wt"
)

for (old in names(rename_map)) {
  if (old %in% names(raw)) names(raw)[names(raw) == old] <- rename_map[[old]]
}

required_cols <- c("group_raw", "wing_id", "AF_wt", "Live_wt")
missing_required <- setdiff(required_cols, names(raw))
if (length(missing_required) > 0) {
  stop("Missing required columns after renaming: ", paste(missing_required, collapse = ", "))
}

bw_vars <- c(
  "BW_birth", "BW_wk1", "BW_wk2", "BW_wk3", "BW_wk5",
  "BW_wk7", "BW_wk9", "BW_wk11", "BW_wk13", "BW_wk15"
)
bw_vars <- bw_vars[bw_vars %in% names(raw)]

numeric_vars <- unique(c(bw_vars, "Live_wt", "AF_wt", "Carcass_wt", "Semi_evis_wt", "Evis_wt"))

optional_id_cols <- c("leg_id", "aluminum_leg_id", "muscle_rna_id", "cecum_meta_id", "duodenum_rna_id")
for (v in optional_id_cols) {
  if (!v %in% names(raw)) raw[[v]] <- NA_character_
}

dat0 <- raw %>%
  mutate(
    group = normalize_group(group_raw),
    wing_id = normalize_id(wing_id),
    wing_id_std = wing_id,
    leg_id = normalize_id(leg_id),
    aluminum_leg_id = normalize_id(aluminum_leg_id),
    muscle_rna_id = normalize_id(muscle_rna_id),
    cecum_meta_id = normalize_id(cecum_meta_id),
    duodenum_rna_id = normalize_id(duodenum_rna_id),
    is_muscle_rna = is_nonempty(muscle_rna_id),
    is_cecum_meta = is_nonempty(cecum_meta_id),
    is_duodenum_rna = is_nonempty(duodenum_rna_id),
    is_any_omics = is_muscle_rna | is_cecum_meta | is_duodenum_rna,
    is_sequenced = is_any_omics,
    seq_label = ifelse(is_sequenced, "Sequenced/omics", "Population")
  ) %>%
  make_numeric_if_exists(numeric_vars) %>%
  filter(!is.na(group), !is.na(wing_id))

# =============================================================================
# 4. Quality control
# =============================================================================

qc_missing <- dat0 %>%
  summarise(
    N_total = n(),
    Missing_AF_wt = sum(is.na(AF_wt)),
    Missing_Live_wt = sum(is.na(Live_wt)),
    Zero_AF_wt = sum(!is.na(AF_wt) & AF_wt == 0),
    Negative_AF_wt = sum(!is.na(AF_wt) & AF_wt < 0),
    Nonpositive_Live_wt = sum(!is.na(Live_wt) & Live_wt <= 0)
  )
write.csv(qc_missing, file.path(AF_OUTDIR, "QC_missing_and_basic_AF_records.csv"), row.names = FALSE, fileEncoding = "UTF-8")

negative_af_table <- dat0 %>%
  filter(!is.na(AF_wt), AF_wt < 0) %>%
  select(group, wing_id, leg_id, aluminum_leg_id, AF_wt, Live_wt, is_any_omics)
write.csv(negative_af_table, file.path(AF_OUTDIR, "QC_negative_AF_values.csv"), row.names = FALSE, fileEncoding = "UTF-8")

dat_qc <- dat0
if (SET_NEGATIVE_AF_TO_NA) {
  dat_qc <- dat_qc %>% mutate(AF_wt = ifelse(!is.na(AF_wt) & AF_wt < 0, NA_real_, AF_wt))
}
if (SET_NONPOSITIVE_LIVEWT_TO_NA) {
  dat_qc <- dat_qc %>% mutate(Live_wt = ifelse(!is.na(Live_wt) & Live_wt <= 0, NA_real_, Live_wt))
}

# Longitudinal body-weight consistency check.
bw_anomaly_table <- tibble()
negative_growth_individuals_removed <- tibble()

if (length(bw_vars) >= 2) {
  bw_long_qc <- dat_qc %>%
    select(group, wing_id, is_any_omics, all_of(bw_vars)) %>%
    pivot_longer(cols = all_of(bw_vars), names_to = "Trait", values_to = "BW") %>%
    mutate(
      Week = case_when(
        Trait == "BW_birth" ~ 0,
        TRUE ~ as.numeric(stringr::str_extract(Trait, "\\d+"))
      )
    ) %>%
    filter(!is.na(BW), !is.na(Week)) %>%
    arrange(group, wing_id, Week)

  bw_anomaly_table <- bw_long_qc %>%
    group_by(group, wing_id) %>%
    mutate(
      Previous_week = lag(Week),
      Previous_BW = lag(BW),
      Delta_BW = BW - Previous_BW
    ) %>%
    ungroup() %>%
    filter(!is.na(Delta_BW), Delta_BW < 0) %>%
    transmute(
      group, wing_id,
      From_week = Previous_week,
      To_week = Week,
      Previous_BW,
      Current_BW = BW,
      Delta_BW,
      is_any_omics
    )

  write.csv(bw_anomaly_table, file.path(AF_OUTDIR, "QC_negative_growth_intervals.csv"), row.names = FALSE, fileEncoding = "UTF-8")

  if (nrow(bw_anomaly_table) > 0) {
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
            any_of(bw_vars),
            AF_wt, Live_wt, any_of(c("Carcass_wt", "Semi_evis_wt", "Evis_wt"))
          ),
        by = c("group", "wing_id")
      ) %>%
      arrange(group, wing_id)

    write.csv(
      negative_growth_individuals_removed,
      file.path(AF_OUTDIR, "QC_negative_growth_individuals_removed.csv"),
      row.names = FALSE,
      fileEncoding = "UTF-8"
    )

    write.csv(
      negative_growth_individuals_removed %>% count(group, name = "N_removed"),
      file.path(AF_OUTDIR, "QC_negative_growth_individuals_removed_summary.csv"),
      row.names = FALSE,
      fileEncoding = "UTF-8"
    )

    if (FILTER_NEGATIVE_GROWTH_INDIVIDUALS) {
      dat_qc <- dat_qc %>%
        anti_join(
          negative_growth_individuals_removed %>% select(group, wing_id),
          by = c("group", "wing_id")
        )
    }
  } else {
    write.csv(tibble(), file.path(AF_OUTDIR, "QC_negative_growth_individuals_removed.csv"), row.names = FALSE)
    write.csv(tibble(), file.path(AF_OUTDIR, "QC_negative_growth_individuals_removed_summary.csv"), row.names = FALSE)
  }
} else {
  write.csv(tibble(), file.path(AF_OUTDIR, "QC_negative_growth_intervals.csv"), row.names = FALSE)
  write.csv(tibble(), file.path(AF_OUTDIR, "QC_negative_growth_individuals_removed.csv"), row.names = FALSE)
  write.csv(tibble(), file.path(AF_OUTDIR, "QC_negative_growth_individuals_removed_summary.csv"), row.names = FALSE)
}

# Ensure optional carcass-composition denominators exist even if absent from the workbook.
for (v in c("Carcass_wt", "Semi_evis_wt", "Evis_wt")) {
  if (!v %in% names(dat_qc)) dat_qc[[v]] <- NA_real_
}

# Derived abdominal fat traits.
dat_af <- dat_qc %>%
  mutate(
    AF_live_pct = safe_ratio(AF_wt, Live_wt, multiplier = 100),
    AF_carcass_pct = safe_ratio(AF_wt, Carcass_wt, multiplier = 100),
    AF_evis_pct = safe_ratio(AF_wt, Evis_wt, multiplier = 100),
    AF_measurable = case_when(
      is.na(AF_wt) ~ NA,
      AF_wt > 0 ~ TRUE,
      AF_wt == 0 ~ FALSE,
      TRUE ~ NA
    ),
    AF_measurable_num = case_when(
      is.na(AF_measurable) ~ NA_real_,
      AF_measurable ~ 100,
      !AF_measurable ~ 0
    )
  )

# Live-weight-adjusted residuals.
resid_fit_data <- dat_af %>%
  filter(!is.na(AF_wt), !is.na(Live_wt))
if (!INCLUDE_ZERO_AF_IN_REGRESSION) {
  resid_fit_data <- resid_fit_data %>% filter(AF_wt > 0)
}

if (nrow(resid_fit_data) >= 6 && n_distinct(resid_fit_data$Live_wt) >= 3) {
  af_resid_fit <- lm(AF_wt ~ Live_wt, data = resid_fit_data)
  af_resid_df <- resid_fit_data %>%
    mutate(
      AF_expected_by_Live_wt = fitted(af_resid_fit),
      AF_resid_livewt = residuals(af_resid_fit)
    ) %>%
    select(group, wing_id, AF_expected_by_Live_wt, AF_resid_livewt)

  dat_af <- dat_af %>%
    left_join(af_resid_df, by = c("group", "wing_id"))

  sink(file.path(AF_OUTDIR, "Table4A_lm_AF_wt_by_Live_wt_model_summary.txt"))
  cat("Model for residual calculation: AF_wt ~ Live_wt\n")
  cat("INCLUDE_ZERO_AF_IN_REGRESSION = ", INCLUDE_ZERO_AF_IN_REGRESSION, "\n\n", sep = "")
  print(summary(af_resid_fit))
  sink()
} else {
  dat_af$AF_expected_by_Live_wt <- NA_real_
  dat_af$AF_resid_livewt <- NA_real_
  writeLines(
    "Insufficient non-missing AF_wt and Live_wt data for AF_wt ~ Live_wt residual model.",
    con = file.path(AF_OUTDIR, "Table4A_lm_AF_wt_by_Live_wt_model_summary.txt")
  )
}


# Additional size-adjusted residuals using carcass and eviscerated weight.
# These are useful sensitivity analyses because carcass/eviscerated weights reduce
# the influence of live gut fill and other non-carcass components.
add_size_adjusted_residual <- function(df, size_col, expected_col, resid_col, label) {
  if (!size_col %in% names(df)) {
    df[[expected_col]] <- NA_real_
    df[[resid_col]] <- NA_real_
    return(df)
  }

  fit_data <- df %>% filter(!is.na(AF_wt), !is.na(.data[[size_col]]))
  if (!INCLUDE_ZERO_AF_IN_REGRESSION) {
    fit_data <- fit_data %>% filter(AF_wt > 0)
  }

  if (nrow(fit_data) >= 6 && n_distinct(fit_data[[size_col]]) >= 3) {
    form <- as.formula(paste("AF_wt ~", size_col))
    fit <- lm(form, data = fit_data)
    resid_df <- fit_data %>%
      mutate(
        !!expected_col := fitted(fit),
        !!resid_col := residuals(fit)
      ) %>%
      select(group, wing_id, all_of(c(expected_col, resid_col)))

    df <- df %>% left_join(resid_df, by = c("group", "wing_id"))

    sink(file.path(AF_OUTDIR, paste0("Table4A_lm_AF_wt_by_", size_col, "_model_summary.txt")))
    cat("Model for residual calculation: AF_wt ~ ", size_col, "\n", sep = "")
    cat("Adjustment label: ", label, "\n", sep = "")
    cat("INCLUDE_ZERO_AF_IN_REGRESSION = ", INCLUDE_ZERO_AF_IN_REGRESSION, "\n\n", sep = "")
    print(summary(fit))
    sink()
  } else {
    df[[expected_col]] <- NA_real_
    df[[resid_col]] <- NA_real_
    writeLines(
      paste0("Insufficient data for AF_wt ~ ", size_col, " residual model."),
      con = file.path(AF_OUTDIR, paste0("Table4A_lm_AF_wt_by_", size_col, "_model_summary.txt"))
    )
  }
  df
}

dat_af <- add_size_adjusted_residual(dat_af, "Carcass_wt", "AF_expected_by_Carcass_wt", "AF_resid_carcasswt", "carcass-weight-adjusted")
dat_af <- add_size_adjusted_residual(dat_af, "Evis_wt", "AF_expected_by_Evis_wt", "AF_resid_eviswt", "eviscerated-weight-adjusted")

# ANCOVA model for group effect after live-weight adjustment.
ancova_data <- dat_af %>% filter(!is.na(AF_wt), !is.na(Live_wt), !is.na(group))
if (!INCLUDE_ZERO_AF_IN_REGRESSION) {
  ancova_data <- ancova_data %>% filter(AF_wt > 0)
}

ancova_coeff <- tibble()
ancova_anova <- tibble()
if (nrow(ancova_data) >= 9 && n_distinct(ancova_data$group) >= 2 && n_distinct(ancova_data$Live_wt) >= 3) {
  af_ancova_fit <- lm(AF_wt ~ Live_wt + group, data = ancova_data)
  ancova_coeff <- as.data.frame(summary(af_ancova_fit)$coefficients) %>%
    rownames_to_column("Term") %>%
    as_tibble() %>%
    rename(SE = `Std. Error`, T_value = `t value`, P_value = `Pr(>|t|)`)
  ancova_anova <- as.data.frame(anova(af_ancova_fit)) %>%
    rownames_to_column("Term") %>%
    as_tibble()
  write.csv(ancova_coeff, file.path(AF_OUTDIR, "Table4B_ANCOVA_AF_wt_coefficients.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(ancova_anova, file.path(AF_OUTDIR, "Table4C_ANCOVA_AF_wt_anova.csv"), row.names = FALSE, fileEncoding = "UTF-8")
} else {
  write.csv(ancova_coeff, file.path(AF_OUTDIR, "Table4B_ANCOVA_AF_wt_coefficients.csv"), row.names = FALSE)
  write.csv(ancova_anova, file.path(AF_OUTDIR, "Table4C_ANCOVA_AF_wt_anova.csv"), row.names = FALSE)
}


# Additional ANCOVA models using alternative denominators/adjusters.
run_size_adjusted_ancova <- function(df, size_col, model_label, positive_only = FALSE) {
  if (!size_col %in% names(df)) return(list(coeff = tibble(), anova = tibble()))
  model_data <- df %>% filter(!is.na(AF_wt), !is.na(.data[[size_col]]), !is.na(group))
  if (positive_only) model_data <- model_data %>% filter(AF_wt > 0)
  if (nrow(model_data) < 9 || n_distinct(model_data$group) < 2 || n_distinct(model_data[[size_col]]) < 3) {
    return(list(coeff = tibble(), anova = tibble()))
  }
  fit <- lm(as.formula(paste("AF_wt ~", size_col, "+ group")), data = model_data)
  coeff <- as.data.frame(summary(fit)$coefficients) %>%
    rownames_to_column("Term") %>%
    as_tibble() %>%
    rename(SE = `Std. Error`, T_value = `t value`, P_value = `Pr(>|t|)`) %>%
    mutate(Model = model_label, Adjuster = size_col, Positive_only = positive_only, .before = 1)
  anova_tbl <- as.data.frame(anova(fit)) %>%
    rownames_to_column("Term") %>%
    as_tibble() %>%
    mutate(Model = model_label, Adjuster = size_col, Positive_only = positive_only, .before = 1)
  list(coeff = coeff, anova = anova_tbl)
}

all_size_models <- list(
  Live_all = run_size_adjusted_ancova(dat_af, "Live_wt", "AF_wt ~ Live_wt + group", positive_only = FALSE),
  Carcass_all = run_size_adjusted_ancova(dat_af, "Carcass_wt", "AF_wt ~ Carcass_wt + group", positive_only = FALSE),
  Evis_all = run_size_adjusted_ancova(dat_af, "Evis_wt", "AF_wt ~ Evis_wt + group", positive_only = FALSE),
  Live_positive = run_size_adjusted_ancova(dat_af, "Live_wt", "AF_wt ~ Live_wt + group, AF_wt > 0", positive_only = TRUE),
  Carcass_positive = run_size_adjusted_ancova(dat_af, "Carcass_wt", "AF_wt ~ Carcass_wt + group, AF_wt > 0", positive_only = TRUE),
  Evis_positive = run_size_adjusted_ancova(dat_af, "Evis_wt", "AF_wt ~ Evis_wt + group, AF_wt > 0", positive_only = TRUE)
)
size_adjusted_coeff <- purrr::map_dfr(all_size_models, "coeff")
size_adjusted_anova <- purrr::map_dfr(all_size_models, "anova")
write.csv(size_adjusted_coeff, file.path(AF_OUTDIR, "Table4E_size_adjusted_ANCOVA_coefficients_all_models.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(size_adjusted_anova, file.path(AF_OUTDIR, "Table4F_size_adjusted_ANCOVA_anova_all_models.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# Positive-only and below-measurable sensitivity analyses for 0-g AF records.
# These analyses address the possibility that 0-g records mean "too small to be
# reliably separated/weighed" rather than true biological absence.
min_positive_af <- suppressWarnings(min(dat_af$AF_wt[!is.na(dat_af$AF_wt) & dat_af$AF_wt > 0], na.rm = TRUE))
if (!is.finite(min_positive_af)) min_positive_af <- NA_real_
zero_replacement_scenarios <- tibble(
  Scenario = c("Retain zeros", "Exclude zeros", "Replace zeros with 0.5 g", "Replace zeros with 1 g", "Replace zeros with half minimum positive"),
  Replacement = c(0, NA_real_, 0.5, 1, ifelse(is.na(min_positive_af), NA_real_, min_positive_af / 2)),
  Exclude_zero = c(FALSE, TRUE, FALSE, FALSE, FALSE)
)

make_sensitivity_dataset <- function(df, replacement, exclude_zero = FALSE) {
  out <- df
  if (exclude_zero) {
    out <- out %>% filter(is.na(AF_wt) | AF_wt > 0)
    out$AF_wt_sens <- out$AF_wt
  } else {
    out$AF_wt_sens <- ifelse(!is.na(out$AF_wt) & out$AF_wt == 0, replacement, out$AF_wt)
  }
  out <- out %>%
    mutate(
      AF_live_pct_sens = safe_ratio(AF_wt_sens, Live_wt, multiplier = 100),
      AF_carcass_pct_sens = safe_ratio(AF_wt_sens, Carcass_wt, multiplier = 100),
      AF_evis_pct_sens = safe_ratio(AF_wt_sens, Evis_wt, multiplier = 100)
    )
  out
}

sensitivity_heterosis_one <- function(df, scenario, trait, label, unit = "%") {
  if (!trait %in% names(df) || sum(!is.na(df[[trait]])) < 6) return(tibble())
  het <- favorable_heterosis_one_trait(df, trait, label, unit, boot_n = BOOT_N) %>%
    mutate(Scenario = scenario, .before = 1)
  test_data <- df %>% filter(!is.na(group), !is.na(.data[[trait]]))
  kw_p <- tryCatch(kruskal.test(test_data[[trait]] ~ test_data$group)$p.value, error = function(e) NA_real_)
  het %>% mutate(Kruskal_P = kw_p)
}

sensitivity_heterosis <- purrr::pmap_dfr(
  zero_replacement_scenarios,
  function(Scenario, Replacement, Exclude_zero) {
    tmp <- make_sensitivity_dataset(dat_af, Replacement, Exclude_zero)
    bind_rows(
      sensitivity_heterosis_one(tmp, Scenario, "AF_wt_sens", "Abdominal fat weight", "%"),
      sensitivity_heterosis_one(tmp, Scenario, "AF_live_pct_sens", "Abdominal fat/live weight (%)", "%"),
      sensitivity_heterosis_one(tmp, Scenario, "AF_carcass_pct_sens", "Abdominal fat/carcass weight (%)", "%"),
      sensitivity_heterosis_one(tmp, Scenario, "AF_evis_pct_sens", "Abdominal fat/eviscerated weight (%)", "%")
    )
  }
)
write.csv(sensitivity_heterosis, file.path(AF_OUTDIR, "Table6_zero_AF_sensitivity_heterosis.csv"), row.names = FALSE, fileEncoding = "UTF-8")

positive_af_data <- dat_af %>% filter(!is.na(AF_wt), AF_wt > 0)
positive_af_traits <- c("AF_wt", "AF_live_pct", "AF_carcass_pct", "AF_evis_pct")
positive_af_traits <- positive_af_traits[positive_af_traits %in% names(positive_af_data)]
positive_af_trait_summary <- if (nrow(positive_af_data) > 0 && length(positive_af_traits) > 0) {
  purrr::map_dfr(positive_af_traits, function(v) summarise_trait_by_group(positive_af_data, v, v))
} else {
  tibble()
}
write.csv(positive_af_trait_summary, file.path(AF_OUTDIR, "Table7_positive_only_AF_trait_summary.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# =============================================================================
# 5. Tables: summaries, proportions, associations and favorable heterosis
# =============================================================================

sample_counts <- dat_af %>%
  group_by(group) %>%
  summarise(
    N_total_after_QC = n(),
    N_any_omics = sum(is_any_omics, na.rm = TRUE),
    N_AF_nonmissing = sum(!is.na(AF_wt)),
    N_AF_measurable = sum(AF_measurable, na.rm = TRUE),
    N_AF_zero = sum(!is.na(AF_wt) & AF_wt == 0),
    N_Live_wt_nonmissing = sum(!is.na(Live_wt)),
    N_Carcass_wt_nonmissing = sum(!is.na(Carcass_wt)),
    N_Evis_wt_nonmissing = sum(!is.na(Evis_wt)),
    .groups = "drop"
  )
write.csv(sample_counts, file.path(AF_OUTDIR, "Table1_sample_counts_AF_QC.csv"), row.names = FALSE, fileEncoding = "UTF-8")

summary_traits <- c("Live_wt", "Carcass_wt", "Evis_wt", "AF_wt", "AF_live_pct", "AF_carcass_pct", "AF_evis_pct", "AF_resid_livewt", "AF_resid_carcasswt", "AF_resid_eviswt")
summary_labels <- c(
  Live_wt = "Live weight",
  Carcass_wt = "Carcass weight",
  Evis_wt = "Eviscerated weight",
  AF_wt = "Abdominal fat weight",
  AF_live_pct = "Abdominal fat/live weight (%)",
  AF_carcass_pct = "Abdominal fat/carcass weight (%)",
  AF_evis_pct = "Abdominal fat/eviscerated weight (%)",
  AF_resid_livewt = "Live-weight-adjusted abdominal fat residual",
  AF_resid_carcasswt = "Carcass-weight-adjusted abdominal fat residual",
  AF_resid_eviswt = "Eviscerated-weight-adjusted abdominal fat residual"
)
summary_traits <- summary_traits[summary_traits %in% names(dat_af)]

trait_summary <- purrr::map_dfr(summary_traits, function(v) {
  summarise_trait_by_group(dat_af, v, summary_labels[[v]])
})
write.csv(trait_summary, file.path(AF_OUTDIR, "Table2_abdominal_fat_trait_summary_by_group.csv"), row.names = FALSE, fileEncoding = "UTF-8")

pairwise_tests <- purrr::map_dfr(summary_traits, ~ pairwise_trait_test(dat_af, .x))
write.csv(pairwise_tests, file.path(AF_OUTDIR, "Table2B_pairwise_wilcoxon_AF_traits.csv"), row.names = FALSE, fileEncoding = "UTF-8")

measurable_prop <- dat_af %>%
  filter(!is.na(AF_measurable), !is.na(group)) %>%
  group_by(group) %>%
  summarise(
    N_nonmissing_AF = n(),
    N_measurable_AF = sum(AF_measurable),
    N_not_measurable_AF = sum(!AF_measurable),
    Proportion_measurable = N_measurable_AF / N_nonmissing_AF,
    Percent_measurable = Proportion_measurable * 100,
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    CI_low = ifelse(N_nonmissing_AF > 0, binom.test(N_measurable_AF, N_nonmissing_AF)$conf.int[1] * 100, NA_real_),
    CI_high = ifelse(N_nonmissing_AF > 0, binom.test(N_measurable_AF, N_nonmissing_AF)$conf.int[2] * 100, NA_real_),
    Label = paste0(N_measurable_AF, "/", N_nonmissing_AF)
  ) %>%
  ungroup()
write.csv(measurable_prop, file.path(AF_OUTDIR, "Table3_measurable_abdominal_fat_proportion.csv"), row.names = FALSE, fileEncoding = "UTF-8")

fisher_table <- dat_af %>%
  filter(!is.na(AF_measurable), !is.na(group)) %>%
  count(group, AF_measurable) %>%
  tidyr::pivot_wider(names_from = AF_measurable, values_from = n, values_fill = 0)

sink(file.path(AF_OUTDIR, "Table3B_measurable_abdominal_fat_fisher_test.txt"))
cat("Fisher's exact test for measurable abdominal fat status across groups\n\n")
print(fisher_table)
if (nrow(fisher_table) >= 2) {
  tab <- as.matrix(fisher_table %>% select(-group))
  rownames(tab) <- fisher_table$group
  if (ncol(tab) >= 2) {
    print(fisher.test(tab))
  } else {
    cat("Fisher test was not run because only one measurable-status category was present.\n")
  }
}
sink()

association_by_group_one <- function(df, size_col) {
  if (!size_col %in% names(df)) return(tibble())
  sub <- df %>% filter(!is.na(AF_wt), !is.na(.data[[size_col]]), !is.na(group))
  if (nrow(sub) == 0) return(tibble())
  split(sub, sub$group) %>%
    purrr::imap_dfr(function(d, g) {
      n <- nrow(d)
      x <- d[[size_col]]
      slope <- if (n >= 3 && n_distinct(x) >= 2) {
        coef(lm(AF_wt ~ adjuster, data = tibble(AF_wt = d$AF_wt, adjuster = x)))[2]
      } else {
        NA_real_
      }
      tibble(
        group = g,
        Adjuster = size_col,
        N = n,
        Pearson_r = ifelse(n >= 3, cor(x, d$AF_wt, method = "pearson"), NA_real_),
        Spearman_r = ifelse(n >= 3, cor(x, d$AF_wt, method = "spearman"), NA_real_),
        Slope_AF_per_g_adjuster = as.numeric(slope)
      )
    })
}
association_by_group <- purrr::map_dfr(c("Live_wt", "Carcass_wt", "Evis_wt"), ~ association_by_group_one(dat_af, .x))
write.csv(association_by_group, file.path(AF_OUTDIR, "Table4D_AF_wt_size_weight_association_by_group.csv"), row.names = FALSE, fileEncoding = "UTF-8")

heterosis_trait_specs <- tribble(
  ~Trait, ~Trait_label, ~Unit,
  "AF_wt", "Abdominal fat weight", "%",
  "AF_live_pct", "Abdominal fat/live (%)", "%",
  "AF_carcass_pct", "Abdominal fat/carcass (%)", "%",
  "AF_evis_pct", "Abdominal fat/eviscerated (%)", "%",
  "AF_measurable_num", "Measurable abdominal fat (%)", "%"
)
for (resid_spec in list(
  c("AF_resid_livewt", "AF residual adjusted for live weight"),
  c("AF_resid_carcasswt", "AF residual adjusted for carcass weight"),
  c("AF_resid_eviswt", "AF residual adjusted for eviscerated weight")
)) {
  resid_trait <- resid_spec[[1]]
  resid_label <- resid_spec[[2]]
  if (resid_trait %in% names(dat_af) && sum(!is.na(dat_af[[resid_trait]])) > 0) {
    heterosis_trait_specs <- bind_rows(
      heterosis_trait_specs,
      tribble(~Trait, ~Trait_label, ~Unit, resid_trait, resid_label, "g")
    )
  }
}

favorable_heterosis <- purrr::pmap_dfr(
  heterosis_trait_specs,
  function(Trait, Trait_label, Unit) {
    favorable_heterosis_one_trait(dat_af, Trait, Trait_label, Unit, boot_n = BOOT_N)
  }
)
write.csv(favorable_heterosis, file.path(AF_OUTDIR, "Table5_favorable_heterosis_abdominal_fat_traits.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# Export cleaned data for integration/checking.
dat_af_export <- dat_af %>%
  select(
    group, wing_id, wing_id_std, leg_id, aluminum_leg_id,
    muscle_rna_id, cecum_meta_id, duodenum_rna_id,
    is_any_omics, is_sequenced, seq_label,
    any_of(bw_vars),
    Live_wt, Carcass_wt, Semi_evis_wt, Evis_wt,
    AF_wt, AF_live_pct, AF_carcass_pct, AF_evis_pct, AF_measurable, AF_measurable_num,
    AF_expected_by_Live_wt, AF_resid_livewt,
    AF_expected_by_Carcass_wt, AF_resid_carcasswt,
    AF_expected_by_Evis_wt, AF_resid_eviswt
  )
write.csv(dat_af_export, file.path(AF_OUTDIR, "Table0_abdominal_fat_clean_data.csv"), row.names = FALSE, fileEncoding = "UTF-8")

writexl::write_xlsx(
  list(
    clean_data = dat_af_export,
    sample_counts = sample_counts,
    trait_summary = trait_summary,
    pairwise_tests = pairwise_tests,
    measurable_prop = measurable_prop,
    association_by_group = association_by_group,
    ANCOVA_coefficients = ancova_coeff,
    ANCOVA_anova = ancova_anova,
    size_adjusted_ANCOVA_coefficients = size_adjusted_coeff,
    size_adjusted_ANCOVA_anova = size_adjusted_anova,
    favorable_heterosis = favorable_heterosis,
    zero_AF_sensitivity = sensitivity_heterosis,
    positive_only_AF_summary = positive_af_trait_summary,
    negative_growth_removed = negative_growth_individuals_removed,
    negative_growth_intervals = bw_anomaly_table,
    negative_AF_values = negative_af_table
  ),
  path = file.path(AF_OUTDIR, "abdominal_fat_analysis_tables.xlsx")
)

# =============================================================================
# 6. Figures A-F: optimized manuscript plotting section
# =============================================================================
# This block replaces the original Figure A-F plotting section.
# It removes redundant single-trait facet strips, shortens subplot titles,
# moves long explanations into a single figure caption, and gives panel F
# more space in the combined layout.

plot_single_box_clean <- function(df, value_col, y_lab, title_text = NULL, y_expand_mult = c(0.03, 0.08)) {
  box_data <- df %>%
    select(group, wing_id_std, is_sequenced, seq_label, value = all_of(value_col)) %>%
    filter(!is.na(value), !is.na(group))

  ggplot(box_data, aes(x = group, y = value)) +
    geom_boxplot(aes(fill = group), alpha = 0.35, outlier.shape = NA, width = 0.65, linewidth = 0.45) +
    geom_jitter(
      data = box_data %>% filter(!is_sequenced),
      aes(color = group),
      alpha = 0.16,
      width = 0.18,
      size = 0.8,
      show.legend = FALSE
    ) +
    geom_jitter(
      data = box_data %>% filter(is_sequenced),
      aes(color = group),
      width = 0.11,
      size = 2.6,
      shape = 21,
      fill = "white",
      stroke = 1.1,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = GROUP_COLORS, drop = FALSE, guide = "none") +
    scale_color_manual(values = GROUP_COLORS, drop = FALSE, guide = "none") +
    scale_y_continuous(expand = expansion(mult = y_expand_mult)) +
    labs(title = title_text, x = "Group", y = y_lab) +
    ggpubr::theme_pubr(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 11.5, margin = margin(b = 3)),
      axis.title = element_text(face = "bold"),
      plot.margin = margin(5, 5, 5, 5)
    )
}

# A. Abdominal fat weight
pA <- plot_single_box_clean(
  dat_af,
  value_col = "AF_wt",
  y_lab = "Abdominal fat weight (g)",
  title_text = "Abdominal fat weight"
)

# B. Abdominal fat percentage
pB <- plot_single_box_clean(
  dat_af,
  value_col = "AF_live_pct",
  y_lab = "Abdominal fat/live weight (%)",
  title_text = "Abdominal fat percentage"
)

# C. Proportion with measurable abdominal fat
pC <- ggplot(measurable_prop, aes(x = group, y = Percent_measurable, fill = group)) +
  geom_col(width = 0.65, alpha = 0.55) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.16, linewidth = 0.55) +
  geom_text(aes(label = Label), vjust = -0.45, size = 3.3) +
  scale_fill_manual(values = GROUP_COLORS, drop = FALSE, guide = "none") +
  scale_y_continuous(
    limits = c(0, 115),
    breaks = seq(0, 100, 25),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Measurable abdominal fat",
    x = "Group",
    y = "Individuals with measurable AF (%)"
  ) +
  ggpubr::theme_pubr(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 11.5, margin = margin(b = 3)),
    axis.title = element_text(face = "bold"),
    plot.margin = margin(5, 5, 5, 5)
  )

# D. Abdominal fat weight versus live weight
scatter_data <- dat_af %>% filter(!is.na(AF_wt), !is.na(Live_wt), !is.na(group))
pD <- ggplot(scatter_data, aes(x = Live_wt, y = AF_wt, color = group)) +
  geom_point(
    data = scatter_data %>% filter(!is_sequenced),
    alpha = 0.20,
    size = 1.0,
    show.legend = FALSE
  ) +
  geom_point(
    data = scatter_data %>% filter(is_sequenced),
    aes(color = group),
    shape = 21,
    fill = "white",
    stroke = 1.1,
    size = 2.7,
    show.legend = FALSE
  ) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.75, alpha = 0.14) +
  scale_color_manual(values = GROUP_COLORS, drop = FALSE) +
  labs(
    title = "AF weight versus live weight",
    x = "Live weight (g)",
    y = "Abdominal fat weight (g)",
    color = "Group"
  ) +
  ggpubr::theme_pubr(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 11.5, margin = margin(b = 3)),
    axis.title = element_text(face = "bold"),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(5, 5, 5, 5)
  )

# E. Live-weight-adjusted abdominal fat residual
pE <- plot_single_box_clean(
  dat_af,
  value_col = "AF_resid_livewt",
  y_lab = "Residual AF weight (g)",
  title_text = "Live-weight-adjusted AF residual",
  y_expand_mult = c(0.08, 0.10)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.4)

# F. Favorable heterosis in abdominal fat-related traits
het_long <- favorable_heterosis %>%
  select(
    Trait, Trait_label, Unit,
    Favorable_MPH, Favorable_MPH_CI_low, Favorable_MPH_CI_high,
    Favorable_BPH, Favorable_BPH_CI_low, Favorable_BPH_CI_high
  ) %>%
  pivot_longer(
    cols = c(Favorable_MPH, Favorable_BPH),
    names_to = "Metric",
    values_to = "Estimate"
  ) %>%
  mutate(
    CI_low = ifelse(Metric == "Favorable_MPH", Favorable_MPH_CI_low, Favorable_BPH_CI_low),
    CI_high = ifelse(Metric == "Favorable_MPH", Favorable_MPH_CI_high, Favorable_BPH_CI_high),
    Metric = recode(Metric, Favorable_MPH = "Favorable MPH", Favorable_BPH = "Favorable BPH"),
    Trait_label_plot = case_when(
      Trait == "AF_wt" ~ "AF weight (%)",
      Trait == "AF_live_pct" ~ "AF/live weight (%)",
      Trait == "AF_carcass_pct" ~ "AF/carcass weight (%)",
      Trait == "AF_evis_pct" ~ "AF/eviscerated weight (%)",
      Trait == "AF_measurable_num" ~ "Measurable AF (%)",
      Trait == "AF_resid_livewt" ~ "Live-weight-adjusted AF residual (g)",
      Trait == "AF_resid_carcasswt" ~ "Carcass-weight-adjusted AF residual (g)",
      Trait == "AF_resid_eviswt" ~ "Eviscerated-weight-adjusted AF residual (g)",
      TRUE ~ paste0(Trait_label, " (", Unit, ")")
    ),
    Trait_label_plot = factor(
      Trait_label_plot,
      levels = rev(c(
        "AF weight (%)",
        "AF/live weight (%)",
        "AF/carcass weight (%)",
        "AF/eviscerated weight (%)",
        "Measurable AF (%)",
        "Live-weight-adjusted AF residual (g)",
        "Carcass-weight-adjusted AF residual (g)",
        "Eviscerated-weight-adjusted AF residual (g)"
      ))
    ),
    Metric = factor(Metric, levels = c("Favorable BPH", "Favorable MPH"))
  )

pF <- ggplot(het_long, aes(x = Estimate, y = Trait_label_plot)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.45) +
  geom_segment(
    aes(x = CI_low, xend = CI_high, y = Trait_label_plot, yend = Trait_label_plot),
    color = "grey45",
    linewidth = 0.65,
    na.rm = TRUE
  ) +
  geom_point(aes(fill = Metric), shape = 21, size = 2.8, color = "black", stroke = 0.25, na.rm = TRUE) +
  facet_wrap(~ Metric, nrow = 1, scales = "free_x") +
  labs(
    title = "Favorable heterosis/effect",
    x = "Effect size",
    y = NULL,
    fill = "Metric"
  ) +
  ggpubr::theme_pubr(base_size = 10.2) +
  theme(
    plot.title = element_text(face = "bold", size = 11.5, margin = margin(b = 3)),
    strip.text = element_text(face = "bold", size = 9.2),
    axis.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 8.7),
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5)
  )

save_plot(pA, "FigAF_A_abdominal_fat_weight_boxplot", 4.3, 4.3)
save_plot(pB, "FigAF_B_abdominal_fat_percentage_boxplot", 4.3, 4.3)
save_plot(pC, "FigAF_C_measurable_abdominal_fat_proportion", 4.6, 4.2)
save_plot(pD, "FigAF_D_abdominal_fat_vs_live_weight_regression", 5.4, 4.5)
save_plot(pE, "FigAF_E_live_weight_adjusted_AF_residual_boxplot", 4.5, 4.3)
save_plot(pF, "FigAF_F_favorable_heterosis_abdominal_fat", 7.4, 4.6)

# Combined figure.
# Widths give panel F extra space and avoid title/subtitle overlap.
p_combined <- (pA | pB | pC) / (pD | pE | pF) +
  patchwork::plot_layout(widths = c(1, 1, 1.35), heights = c(1, 1.08), guides = "collect") +
  patchwork::plot_annotation(
    title = "Abdominal fat phenotype and favorable heterosis in YY, YX and XX chickens",
    caption = paste(
      "Small transparent points indicate population individuals; white-filled points indicate sequenced individuals.",
      "AF = abdominal fat. Error bars in panel C indicate exact binomial 95% confidence intervals.",
      "Positive values in panel F indicate a favorable reduction in YX F1 relative to the mid-parent or favorable parent.",
      sep = " "
    ),
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, margin = margin(b = 6)),
      plot.caption = element_text(size = 9, hjust = 0, margin = margin(t = 6)),
      plot.margin = margin(8, 8, 8, 8)
    )
  )

save_plot(p_combined, "FigAF_abdominal_fat_multiplot_A_to_F", 16.5, 9.7)

# Additional stacked proportion plot for measured vs non-measured status.
prop_stacked <- dat_af %>%
  filter(!is.na(AF_measurable), !is.na(group)) %>%
  mutate(Status = ifelse(AF_measurable, "Measurable AF", "No measurable AF")) %>%
  count(group, Status) %>%
  group_by(group) %>%
  mutate(Percent = n / sum(n) * 100) %>%
  ungroup()

pC_stacked <- ggplot(prop_stacked, aes(x = group, y = Percent, fill = Status)) +
  geom_col(width = 0.65, alpha = 0.75) +
  geom_text(
    aes(label = paste0(n, "\n", sprintf("%.1f%%", Percent))),
    position = position_stack(vjust = 0.5),
    size = 3.1
  ) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.03))) +
  labs(
    title = "Measurable versus non-measurable AF status",
    x = "Group",
    y = "Proportion",
    fill = "Status"
  ) +
  ggpubr::theme_pubr(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))
save_plot(pC_stacked, "FigAF_C2_measurable_status_stacked_bar", 4.8, 4.3)


# Supplementary robustness figure: carcass/eviscerated adjustment and zero-record sensitivity.
if ("AF_carcass_pct" %in% names(dat_af) && sum(!is.na(dat_af$AF_carcass_pct)) > 0) {
  pG <- plot_single_box_clean(
    dat_af,
    value_col = "AF_carcass_pct",
    y_lab = "Abdominal fat/carcass weight (%)",
    title_text = "AF/carcass percentage"
  )
  save_plot(pG, "FigAF_G_abdominal_fat_carcass_percentage_boxplot", 4.3, 4.3)
} else {
  pG <- ggplot() + theme_void() + labs(title = "AF/carcass percentage\nnot available")
}

if ("AF_evis_pct" %in% names(dat_af) && sum(!is.na(dat_af$AF_evis_pct)) > 0) {
  pH <- plot_single_box_clean(
    dat_af,
    value_col = "AF_evis_pct",
    y_lab = "Abdominal fat/eviscerated weight (%)",
    title_text = "AF/eviscerated percentage"
  )
  save_plot(pH, "FigAF_H_abdominal_fat_eviscerated_percentage_boxplot", 4.3, 4.3)
} else {
  pH <- ggplot() + theme_void() + labs(title = "AF/eviscerated percentage\nnot available")
}

if ("AF_resid_carcasswt" %in% names(dat_af) && sum(!is.na(dat_af$AF_resid_carcasswt)) > 0) {
  pI <- plot_single_box_clean(
    dat_af,
    value_col = "AF_resid_carcasswt",
    y_lab = "Residual AF weight (g)",
    title_text = "Carcass-weight-adjusted AF residual",
    y_expand_mult = c(0.08, 0.10)
  ) + geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.4)
  save_plot(pI, "FigAF_I_carcass_weight_adjusted_AF_residual_boxplot", 4.5, 4.3)
} else {
  pI <- ggplot() + theme_void() + labs(title = "Carcass-adjusted residual\nnot available")
}

positive_plot_data <- dat_af %>% filter(!is.na(AF_wt), AF_wt > 0)
if (nrow(positive_plot_data) > 0) {
  pJ <- plot_single_box_clean(
    positive_plot_data,
    value_col = "AF_wt",
    y_lab = "Abdominal fat weight (g)",
    title_text = "Positive-only AF weight"
  )
  save_plot(pJ, "FigAF_J_positive_only_AF_weight_boxplot", 4.3, 4.3)
} else {
  pJ <- ggplot() + theme_void() + labs(title = "Positive-only AF\nnot available")
}

sens_plot <- sensitivity_heterosis %>%
  filter(Trait %in% c("AF_wt_sens", "AF_live_pct_sens", "AF_carcass_pct_sens")) %>%
  mutate(
    Trait_label_plot = case_when(
      Trait == "AF_wt_sens" ~ "AF weight",
      Trait == "AF_live_pct_sens" ~ "AF/live %",
      Trait == "AF_carcass_pct_sens" ~ "AF/carcass %",
      TRUE ~ Trait_label
    ),
    Trait_label_plot = factor(Trait_label_plot, levels = rev(c("AF weight", "AF/live %", "AF/carcass %"))),
    Scenario = factor(Scenario, levels = zero_replacement_scenarios$Scenario)
  )

if (nrow(sens_plot) > 0) {
  pK <- ggplot(sens_plot, aes(x = Favorable_MPH, y = Trait_label_plot)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.45) +
    geom_segment(aes(x = Favorable_MPH_CI_low, xend = Favorable_MPH_CI_high, yend = Trait_label_plot), color = "grey45", linewidth = 0.55, na.rm = TRUE) +
    geom_point(shape = 21, fill = "white", color = "black", size = 2.1, stroke = 0.35, na.rm = TRUE) +
    facet_wrap(~ Scenario, nrow = 1) +
    labs(title = "Zero-record sensitivity: favorable MPH", x = "Favorable MPH (%)", y = NULL) +
    ggpubr::theme_pubr(base_size = 9.2) +
    theme(
      plot.title = element_text(face = "bold", size = 11.2),
      strip.text = element_text(face = "bold", size = 7.6),
      axis.text.y = element_text(size = 8.2)
    )
  save_plot(pK, "FigAF_K_zero_record_sensitivity_favorable_MPH", 9.8, 3.2)
} else {
  pK <- ggplot() + theme_void() + labs(title = "Sensitivity analysis\nnot available")
}

p_robust <- (pG | pH | pI | pJ) / pK +
  patchwork::plot_layout(heights = c(1, 0.9)) +
  patchwork::plot_annotation(
    title = "Robustness checks for abdominal fat interpretation",
    caption = paste(
      "AF/carcass and AF/eviscerated percentages test whether the abdominal-fat pattern persists after alternative carcass-composition normalization.",
      "Positive-only analysis excludes 0-g records.",
      "Zero-record sensitivity treats 0-g values as below-measurable records rather than definitive biological absence.",
      sep = " "
    ),
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, margin = margin(b = 6)),
      plot.caption = element_text(size = 8.6, hjust = 0, margin = margin(t = 6))
    )
  )
save_plot(p_robust, "FigAF_robustness_carcass_adjustment_and_zero_sensitivity", 15.5, 8.8)

# =============================================================================
# 7. Console summary
# =============================================================================

message("\n============================================================")
message("Abdominal fat-focused analysis complete.")
message("Input:  ", COMPLETE_FILE)
message("Output: ", AF_OUTDIR)
message("\nKey QC logs:")
message("  QC_negative_growth_intervals.csv")
message("  QC_negative_growth_individuals_removed.csv")
message("  QC_negative_AF_values.csv")
message("\nKey tables:")
message("  Table0_abdominal_fat_clean_data.csv")
message("  Table2_abdominal_fat_trait_summary_by_group.csv")
message("  Table3_measurable_abdominal_fat_proportion.csv")
message("  Table4B_ANCOVA_AF_wt_coefficients.csv")
message("  Table4E_size_adjusted_ANCOVA_coefficients_all_models.csv")
message("  Table5_favorable_heterosis_abdominal_fat_traits.csv")
message("  Table6_zero_AF_sensitivity_heterosis.csv")
message("  Table7_positive_only_AF_trait_summary.csv")
message("  abdominal_fat_analysis_tables.xlsx")
message("\nKey figures:")
message("  FigAF_abdominal_fat_multiplot_A_to_F.png/pdf")
message("  FigAF_A_abdominal_fat_weight_boxplot.png/pdf")
message("  FigAF_B_abdominal_fat_percentage_boxplot.png/pdf")
message("  FigAF_C_measurable_abdominal_fat_proportion.png/pdf")
message("  FigAF_D_abdominal_fat_vs_live_weight_regression.png/pdf")
message("  FigAF_E_live_weight_adjusted_AF_residual_boxplot.png/pdf")
message("  FigAF_F_favorable_heterosis_abdominal_fat.png/pdf")
message("  FigAF_robustness_carcass_adjustment_and_zero_sensitivity.png/pdf")
message("============================================================\n")
