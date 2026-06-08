rm(list = ls())

suppressPackageStartupMessages({
  library(readxl)
  library(tidyverse)
  library(stringr)
  library(ggpubr)
  library(patchwork)
  library(scales)
  library(survival)
  library(survminer)
  library(broom)
  library(grid)
  library(ggrepel)
})

# ============================================================
# 0. paths and analysis settings
# ============================================================

parse_cli_args <- function(defaults) {
  args <- commandArgs(trailingOnly = TRUE)
  config <- defaults

  if (length(args) == 0) {
    return(config)
  }

  for (arg in args) {
    if (!str_starts(arg, "--") || !str_detect(arg, "=")) {
      stop(
        "Arguments must use --name=value format. Example: ",
        "--clinical-xlsx=data/ecDNA-survival.xlsx",
        call. = FALSE
      )
    }

    key_value <- str_split_fixed(str_remove(arg, "^--"), "=", 2)
    key <- str_replace_all(key_value[1], "-", "_")
    value <- key_value[2]

    if (!key %in% names(config)) {
      stop("Unknown argument: --", str_replace_all(key, "_", "-"), call. = FALSE)
    }

    config[[key]] <- value
  }

  config
}

CONFIG <- parse_cli_args(list(
  clinical_xlsx = file.path("data", "ecDNA-survival.xlsx"),
  ecdna_gene_rds = file.path("data", "ecDNA_genelist_all.rds"),
  output_dir = file.path("results", "figure5_survival"),
  clinical_sheet = "临床生存信息汇总表整理-latest",
  wgs_sample_sheet = "ecDNA WGS样本测序统计-refined",
  included_disease_status = "Brain metastases"
))

PATH_CLINICAL_XLSX <- CONFIG$clinical_xlsx
PATH_ECDNA_GENE <- CONFIG$ecdna_gene_rds
SHEET_CLINICAL <- CONFIG$clinical_sheet
SHEET_WGS_SAMPLE <- CONFIG$wgs_sample_sheet
INCLUDED_DISEASE_STATUS <- CONFIG$included_disease_status

DIR_OUT <- CONFIG$output_dir
DIR_RDS <- file.path(DIR_OUT, "rds")
DIR_TABLE <- file.path(DIR_OUT, "tables")
DIR_PLOT <- file.path(DIR_OUT, "plots")

dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_RDS, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_TABLE, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_PLOT, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(PATH_CLINICAL_XLSX)) {
  stop("Clinical workbook not found: ", PATH_CLINICAL_XLSX, call. = FALSE)
}

if (!file.exists(PATH_ECDNA_GENE)) {
  stop("ecDNA gene annotation RDS not found: ", PATH_ECDNA_GENE, call. = FALSE)
}

KM_XLIM <- c(0, 60)
KM_BREAK <- 12

MIN_N_PER_CANCER_FOR_KM <- 8
MIN_N_PER_GROUP_FOR_KM <- 2
MIN_EVENTS_FOR_COX <- 3
MIN_ECDNA_POS_FOR_COX <- 2

ECDNA_2GROUP_COLS <- c(
  "ecDNA-" = "#E0E0E0",
  "ecDNA+" = "#8B0000"
)

ECDNA_2GROUP_KM_COLS <- c(
  "ecDNA-" = "blue",
  "ecDNA+" = "#8B0000"
)

ECDNA_3GROUP_KM_COLS <- c(
  "No amplicon" = "blue",
  "ecDNA- (amplicon)" = "#DE9B13",
  "ecDNA+" = "#8B0000"
)

ONCOGENE_3GROUP_COLS <- c(
  "No amplicon/oncogene" = "#4D4D4D",
  "ecDNA- (oncogene)" = "#D39B14",
  "ecDNA+" = "#8B1A16"
)

CANONICAL_ONCOGENE_3GROUP_COLS <- c(
  "No canonical oncogene / undetected" = "#4D4D4D",
  "ecDNA- / canonical oncogene amplicon" = "#D99A00",
  "ecDNA+" = "#8B0000"
)

theme_figure5 <- theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(face = "bold", size = 13)
  )

# ============================================================
# 1. small helpers
# ============================================================

to_numeric <- function(x) {
  x_chr <- str_trim(as.character(x))
  x_chr[x_chr %in% c("", "/", "NA", "NaN", "Lost")] <- NA_character_
  suppressWarnings(as.numeric(x_chr))
}

to_os_event <- function(x) {
  x_chr <- str_trim(as.character(x))
  case_when(
    x_chr == "1" ~ 1,
    x_chr == "0" ~ 0,
    x_chr == "Lost" ~ 0,
    x_chr %in% c("", "/", "NA", "NaN") ~ NA_real_,
    TRUE ~ NA_real_
  )
}

clean_text <- function(x) {
  x_chr <- str_trim(as.character(x))
  x_chr[x_chr %in% c("", "/", "NA", "NaN")] <- NA_character_
  x_chr
}

standardize_sample_id <- function(x) {
  clean_text(x) %>%
    str_remove("-.*$")
}

safe_median <- function(x) {
  x <- x[is.finite(x) & !is.na(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  median(x)
}

format_p_value <- function(p_value) {
  case_when(
    is.na(p_value) ~ "",
    p_value < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3f", p_value)
  )
}

format_plot_p_value <- function(p_value) {
  case_when(
    is.na(p_value) ~ "p = NA",
    p_value < 0.001 ~ "p < 0.001",
    TRUE ~ paste0("p = ", sprintf("%.3f", p_value))
  )
}

format_facet_p_value <- function(p_value) {
  case_when(
    is.na(p_value) ~ "p = NA",
    p_value < 0.001 ~ "p < 0.001",
    TRUE ~ paste0("p = ", sprintf("%.2f", p_value))
  )
}

wilcox_ecdna_negative_greater <- function(data, y_col) {
  ecDNA_negative_values <- data %>%
    filter(ecDNA_Status == "ecDNA-") %>%
    pull(.data[[y_col]])

  ecDNA_positive_values <- data %>%
    filter(ecDNA_Status == "ecDNA+") %>%
    pull(.data[[y_col]])

  ecDNA_negative_values <- ecDNA_negative_values[!is.na(ecDNA_negative_values)]
  ecDNA_positive_values <- ecDNA_positive_values[!is.na(ecDNA_positive_values)]

  if (length(ecDNA_negative_values) == 0 || length(ecDNA_positive_values) == 0) {
    return(NA_real_)
  }

  test_out <- wilcox.test(
    x = ecDNA_negative_values,
    y = ecDNA_positive_values,
    alternative = "greater",
    exact = FALSE
  )

  test_out$p.value
}

truncate_survival_time <- function(data, time_col, event_col, max_time) {
  data_out <- data
  original_time <- data_out[[time_col]]
  over_limit <- !is.na(original_time) & original_time > max_time
  data_out[[time_col]] <- pmin(original_time, max_time)
  data_out[[event_col]] <- ifelse(over_limit, 0, data_out[[event_col]])
  data_out
}

keep_cancer_types_with_both_groups <- function(data, cancer_col, status_col) {
  group_counts <- data %>%
    filter(!is.na(.data[[cancer_col]]), !is.na(.data[[status_col]])) %>%
    distinct(.data[[cancer_col]], .data[[status_col]]) %>%
    count(.data[[cancer_col]], name = "n_status")

  cancer_types_to_keep <- group_counts %>%
    filter(n_status >= 2) %>%
    pull(.data[[cancer_col]])

  data %>%
    filter(.data[[cancer_col]] %in% cancer_types_to_keep)
}

save_plot_pdf_png <- function(plot_obj, file_stub, width, height, dpi = 300) {
  ggsave(file.path(DIR_PLOT, paste0(file_stub, ".pdf")), plot_obj, width = width, height = height)
  ggsave(file.path(DIR_PLOT, paste0(file_stub, ".png")), plot_obj, width = width, height = height, dpi = dpi)
}

column_or_na <- function(data, column_name) {
  if (column_name %in% names(data)) {
    return(data[[column_name]])
  }
  rep(NA_character_, nrow(data))
}

as_logical_flag <- function(x) {
  x_chr <- str_to_lower(str_trim(as.character(x)))
  x_chr %in% c("true", "t", "1", "yes", "y")
}

# ============================================================
# 2. read and clean the single clinical survival sheet
# ============================================================

clinical_raw <- read_excel(
  path = PATH_CLINICAL_XLSX,
  sheet = SHEET_CLINICAL,
  col_types = "text",
  .name_repair = "unique"
)

wgs_sample_raw <- read_excel(
  path = PATH_CLINICAL_XLSX,
  sheet = SHEET_WGS_SAMPLE,
  .name_repair = "unique"
)

wgs_sample_info <- wgs_sample_raw %>%
  filter(used == "yes", Pri_Met == "Met") %>%
  transmute(
    SampleID = standardize_sample_id(Sample),
    WGS_CancerType = clean_text(Cancer_type),
    WGS_CancerSubtype = clean_text(Sub_Cancer_type),
    WGS_SeqType = clean_text(Seq_type),
    WGS_ecDNA_count = to_numeric(ecDNA_count),
    WGS_ecDNA_oncogene_count = to_numeric(`ecDNA count (with oncogene)`)
  ) %>%
  distinct(SampleID, .keep_all = TRUE)

gene_sample_summary <- if (file.exists(PATH_ECDNA_GENE)) {
  readRDS(PATH_ECDNA_GENE) %>%
    mutate(
      SampleID = standardize_sample_id(Sample),
      is_ecDNA_flag = as_logical_flag(is_ecDNA),
      is_canonical_oncogene_flag = as_logical_flag(is_canonical_oncogene)
    ) %>%
    group_by(SampleID) %>%
    summarise(
      has_ecDNA_gene = any(is_ecDNA_flag, na.rm = TRUE),
      has_canonical_oncogene = any(is_canonical_oncogene_flag, na.rm = TRUE),
      .groups = "drop"
    )
} else {
  tibble(
    SampleID = character(),
    has_ecDNA_gene = logical(),
    has_canonical_oncogene = logical()
  )
}

clinical_analysis <- clinical_raw %>%
  transmute(
    Seq = to_numeric(Seq),
    Uploaded = clean_text(Uploaded),
    PathologyNumberFinal = clean_text(`Pathology number (Final)`),
    PairIDUpdated = clean_text(`Paired ID （Updated）`),
    PathologyNumber = clean_text(column_or_na(clinical_raw, "Pathology number...5")),
    SampleID = standardize_sample_id(`Sample ID`),
    T0_initial_diagnosis = clean_text(`T0 initial_diagnosis`),
    T1_brain_met_date = clean_text(`T1 brain_met_date`),
    T2_followup_date = clean_text(`T2 followup_date`),
    Latency_Months = to_numeric(`Brain_met_time (T1-T0)`),
    OS_Months = to_numeric(`OS_months (T2-T0)`),
    OS_status_raw = clean_text(`OS_status manual`),
    OS_event = to_os_event(`OS_status manual`),
    CancerType = clean_text(`Cancer type`),
    CancerSubtype = clean_text(`Cancer subtype`),
    DiseaseStatus = clean_text(`Disease status`),
    PathologyNumberOriginal = clean_text(column_or_na(clinical_raw, "Pathology number...16")),
    PairID = clean_text(`Paired ID`),
    SequencingType = clean_text(`Sequencing type`),
    Age = to_numeric(Age),
    Sex = clean_text(Sex),
    TreatmentClass = clean_text(`Treatment class`),
    ecDNA_status_raw = clean_text(`ecDNA status`)
  ) %>%
  mutate(
    Sex = str_to_lower(Sex),
    Sex = factor(Sex, levels = c("male", "female")),
    ecDNA_Status_3group = case_when(
      ecDNA_status_raw == "Negative" ~ "No amplicon",
      ecDNA_status_raw == "Amplicon (Linear)" ~ "ecDNA- (amplicon)",
      ecDNA_status_raw == "ecDNA+" ~ "ecDNA+",
      TRUE ~ NA_character_
    ),
    ecDNA_Status_3group = factor(
      ecDNA_Status_3group,
      levels = c("No amplicon", "ecDNA- (amplicon)", "ecDNA+")
    ),
    ecDNA_Status = case_when(
      ecDNA_Status_3group == "ecDNA+" ~ "ecDNA+",
      ecDNA_Status_3group %in% c("No amplicon", "ecDNA- (amplicon)") ~ "ecDNA-",
      TRUE ~ NA_character_
    ),
    ecDNA_Status = factor(ecDNA_Status, levels = c("ecDNA-", "ecDNA+")),
    CancerTypeFacet = case_when(
      CancerType == "Breast cancer" ~ "BRCA",
      CancerType == "Colorectal cancer" ~ "COAD",
      CancerType == "Esophageal cancer" ~ "ESCC",
      CancerType == "Gastric cancer" ~ "STAD",
      CancerType == "Ovarian cancer" ~ "OV",
      CancerType == "Melanoma" ~ "SKCM",
      CancerType == "Lung cancer" & !is.na(CancerSubtype) ~ CancerSubtype,
      TRUE ~ CancerType
    ),
    Age_z = ifelse(sum(!is.na(Age)) >= 2, as.numeric(scale(Age)), NA_real_),
    Cohort = "Internal"
  ) %>%
  filter(
    DiseaseStatus == INCLUDED_DISEASE_STATUS,
    !is.na(SampleID),
    !is.na(ecDNA_Status)
  ) %>%
  left_join(wgs_sample_info, by = "SampleID") %>%
  left_join(gene_sample_summary, by = "SampleID") %>%
  mutate(
    has_ecDNA_gene = replace_na(has_ecDNA_gene, FALSE),
    has_canonical_oncogene = replace_na(has_canonical_oncogene, FALSE),
    ecDNA_gene_binary = ecDNA_status_raw == "ecDNA+" | has_ecDNA_gene,
    OncogeneAmplicon_Group = case_when(
      ecDNA_gene_binary ~ "ecDNA+",
      !ecDNA_gene_binary & has_canonical_oncogene ~ "ecDNA- (oncogene)",
      !ecDNA_gene_binary & !has_canonical_oncogene ~ "No amplicon/oncogene",
      TRUE ~ NA_character_
    ),
    OncogeneAmplicon_Group = factor(
      OncogeneAmplicon_Group,
      levels = c("No amplicon/oncogene", "ecDNA- (oncogene)", "ecDNA+")
    ),
    CanonicalOncogene_Group = case_when(
      ecDNA_gene_binary ~ "ecDNA+",
      !ecDNA_gene_binary & has_canonical_oncogene ~ "ecDNA- / canonical oncogene amplicon",
      !ecDNA_gene_binary & !has_canonical_oncogene ~ "No canonical oncogene / undetected",
      TRUE ~ NA_character_
    ),
    CanonicalOncogene_Group = factor(
      CanonicalOncogene_Group,
      levels = c(
        "No canonical oncogene / undetected",
        "ecDNA- / canonical oncogene amplicon",
        "ecDNA+"
      )
    ),
    CancerSubtypeFigure = case_when(
      CancerType == "Breast cancer" ~ "BRCA",
      CancerType == "Ovarian cancer" ~ "OV",
      !is.na(CancerSubtype) ~ CancerSubtype,
      !is.na(WGS_CancerSubtype) ~ WGS_CancerSubtype,
      TRUE ~ CancerType
    ),
    CancerTypeFacet = factor(
      CancerTypeFacet,
      levels = c("BRCA", "COAD", "ESCC", "LUAD", "LUSC", "Lung-PD", "OV", "SKCM", "STAD")
    ),
    CancerTypeCox = case_when(
      CancerType == "Breast cancer" ~ "BRCA",
      CancerType == "Cervical cancer" ~ "CESC",
      CancerType == "Colorectal cancer" ~ "COAD",
      CancerType == "Esophageal cancer" ~ "ESCC",
      CancerType == "Kidney cancer" ~ "KIRC",
      CancerType == "Ovarian cancer" ~ "OV",
      CancerType == "Melanoma" ~ "SKCM",
      CancerType == "Gastric cancer" ~ "STAD",
      CancerType == "Lung cancer" & !is.na(CancerSubtype) ~ CancerSubtype,
      TRUE ~ NA_character_
    ),
    CancerTypeCox = factor(
      CancerTypeCox,
      levels = c("BRCA", "CESC", "COAD", "ESCC", "KIRC", "LUAD", "Lung-PD", "LUSC", "OV", "SCLC", "SKCM", "STAD")
    )
  ) %>%
  arrange(SampleID, desc(!is.na(OS_Months)), desc(!is.na(Latency_Months))) %>%
  distinct(SampleID, .keep_all = TRUE)

saveRDS(clinical_analysis, file.path(DIR_RDS, "figure5_analysis_table.rds"))
write.csv(clinical_analysis, file.path(DIR_TABLE, "figure5_analysis_table.csv"), row.names = FALSE)

# ============================================================
# 3. descriptive tables
# ============================================================

summary_by_ecdna <- clinical_analysis %>%
  group_by(ecDNA_Status) %>%
  summarise(
    n_total = n(),
    n_latency = sum(!is.na(Latency_Months)),
    median_latency_months = safe_median(Latency_Months),
    n_os = sum(!is.na(OS_Months)),
    median_os_months = safe_median(OS_Months),
    n_os_event = sum(OS_event == 1, na.rm = TRUE),
    n_os_censored = sum(OS_event == 0, na.rm = TRUE),
    n_age = sum(!is.na(Age)),
    n_sex = sum(!is.na(Sex)),
    .groups = "drop"
  )

summary_by_ecdna_3group <- clinical_analysis %>%
  group_by(ecDNA_Status_3group) %>%
  summarise(
    n_total = n(),
    n_latency = sum(!is.na(Latency_Months)),
    median_latency_months = safe_median(Latency_Months),
    n_os = sum(!is.na(OS_Months)),
    median_os_months = safe_median(OS_Months),
    n_os_event = sum(OS_event == 1, na.rm = TRUE),
    .groups = "drop"
  )

cancer_type_counts <- clinical_analysis %>%
  count(CancerType, ecDNA_Status, name = "n") %>%
  group_by(CancerType) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup()

write.csv(summary_by_ecdna, file.path(DIR_TABLE, "figure5_summary_by_ecDNA_2group.csv"), row.names = FALSE)
write.csv(summary_by_ecdna_3group, file.path(DIR_TABLE, "figure5_summary_by_ecDNA_3group.csv"), row.names = FALSE)
write.csv(cancer_type_counts, file.path(DIR_TABLE, "figure5_ecDNA_by_cancer_type_summary.csv"), row.names = FALSE)

# ============================================================
# 4. boxplots for latency and OS
# ============================================================

plot_box_by_ecdna <- function(data, y_col, y_label, title_label, label_y) {
  p_value <- wilcox_ecdna_negative_greater(data, y_col)
  p_label <- paste("Wilcoxon,", format_plot_p_value(p_value))

  ggplot(data, aes(x = ecDNA_Status, y = .data[[y_col]], fill = ecDNA_Status)) +
    stat_boxplot(geom = "errorbar", width = 0.25, linewidth = 0.6, color = "black") +
    geom_boxplot(width = 0.5, linewidth = 0.5, outlier.shape = NA, alpha = 0.9, color = "black") +
    geom_jitter(width = 0.15, height = 0, size = 1.5, alpha = 0.9, color = "#2c3e50") +
    annotate("text", x = 1.5, y = label_y, label = p_label, size = 4.5) +
    scale_fill_manual(values = ECDNA_2GROUP_COLS, drop = FALSE) +
    labs(y = y_label, x = NULL, title = title_label) +
    theme_figure5 +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
      aspect.ratio = 1.1
    )
}

plot_facet_box_by_cancer <- function(data, y_col, y_label, tag_label = NULL) {
  p_labels <- data %>%
    group_by(CancerTypeFacet) %>%
    group_modify(~ {
      tibble(
        p_value = wilcox_ecdna_negative_greater(.x, y_col),
        y_position = max(.x[[y_col]], na.rm = TRUE) * 0.92
      )
    }) %>%
    ungroup() %>%
    mutate(
      x_position = 1.5,
      label = format_facet_p_value(p_value)
    )

  ggplot(data, aes(x = ecDNA_Status, y = .data[[y_col]], fill = ecDNA_Status)) +
    stat_boxplot(geom = "errorbar", width = 0.25, linewidth = 0.5, color = "black") +
    geom_boxplot(width = 0.5, linewidth = 0.5, outlier.shape = NA, alpha = 0.9, color = "black") +
    geom_jitter(width = 0.15, height = 0, size = 1, alpha = 0.85, color = "#2c3e50") +
    geom_text(
      data = p_labels,
      aes(x = x_position, y = y_position, label = label),
      inherit.aes = FALSE,
      size = 3.4
    ) +
    scale_fill_manual(values = ECDNA_2GROUP_COLS, drop = FALSE) +
    labs(y = y_label, x = NULL, tag = tag_label) +
    facet_wrap(~ CancerTypeFacet, nrow = 1, scales = "fixed", strip.position = "top", drop = TRUE) +
    theme_classic(base_size = 13) +
    theme(
      axis.text = element_text(color = "black", size = 12),
      axis.text.x = element_text(color = "black", size = 13, angle = 40, hjust = 1, vjust = 1),
      axis.title = element_text(color = "black", size = 15),
      axis.line = element_line(linewidth = 0.6, lineend = "square"),
      axis.ticks = element_line(linewidth = 0.6, color = "black"),
      legend.position = "none",
      strip.background = element_rect(color = "black", fill = "white", linewidth = 0.6),
      strip.text = element_text(color = "black", size = 15),
      panel.spacing.x = unit(0.12, "lines"),
      panel.spacing.y = unit(0, "lines"),
      plot.tag = element_text(face = "plain", size = 18),
      plot.tag.position = c(-0.035, 1.02),
      plot.margin = margin(4, 6, 4, 22)
    )
}

latency_data <- clinical_analysis %>%
  filter(!is.na(Latency_Months), !is.na(ecDNA_Status))

os_data <- clinical_analysis %>%
  filter(!is.na(OS_Months), !is.na(ecDNA_Status))

latency_box <- plot_box_by_ecdna(
  data = latency_data,
  y_col = "Latency_Months",
  y_label = "Time to brain metastasis (months)",
  title_label = "Latency",
  label_y = max(latency_data$Latency_Months, na.rm = TRUE) * 1.05
)

os_box <- plot_box_by_ecdna(
  data = os_data,
  y_col = "OS_Months",
  y_label = "OS (months)",
  title_label = "OS",
  label_y = max(os_data$OS_Months, na.rm = TRUE) * 1.05
)

latency_facet_data <- keep_cancer_types_with_both_groups(latency_data, "CancerTypeFacet", "ecDNA_Status")
os_facet_data <- keep_cancer_types_with_both_groups(os_data, "CancerTypeFacet", "ecDNA_Status")

latency_facet_box <- plot_facet_box_by_cancer(
  data = latency_facet_data,
  y_col = "Latency_Months",
  y_label = "Time to BrM (months)",
  tag_label = "C"
)

os_facet_box <- plot_facet_box_by_cancer(
  data = os_facet_data,
  y_col = "OS_Months",
  y_label = "OS (months)",
  tag_label = "D"
)

save_plot_pdf_png(latency_box, "figure5_latency_box_2group", width = 4.2, height = 4.8)
save_plot_pdf_png(os_box, "figure5_OS_box_2group", width = 4.2, height = 4.8)
save_plot_pdf_png(latency_facet_box, "figure5_latency_facet_2group_dualStatusOnly", width = 18, height = 3.6)
save_plot_pdf_png(os_facet_box, "figure5_OS_facet_2group_dualStatusOnly", width = 18, height = 3.6)

# ============================================================
# 5. Kaplan-Meier analyses
# ============================================================

make_km_plot <- function(data, time_col, event_col, group_col, palette, x_label, y_label, title_label) {
  km_data <- data %>%
    mutate(
      km_time = .data[[time_col]],
      km_event = .data[[event_col]],
      km_group = .data[[group_col]]
    )

  km_group_levels <- levels(km_data[[group_col]])
  km_data$km_group <- factor(km_data$km_group, levels = km_group_levels)

  fit <- survfit(Surv(km_time, km_event) ~ km_group, data = km_data)

  ggsurvplot(
    fit,
    data = km_data,
    risk.table = TRUE,
    pval = TRUE,
    conf.int = FALSE,
    palette = unname(palette),
    xlab = x_label,
    ylab = y_label,
    title = title_label,
    legend.title = NULL,
    legend.labs = km_group_levels,
    risk.table.height = 0.23,
    xlim = KM_XLIM,
    break.time.by = KM_BREAK,
    ggtheme = theme_classic(base_size = 12)
  )
}

km_os_data <- clinical_analysis %>%
  filter(!is.na(OS_Months), !is.na(OS_event), !is.na(ecDNA_Status))

km_os_plot_data <- truncate_survival_time(
  data = km_os_data,
  time_col = "OS_Months",
  event_col = "OS_event",
  max_time = KM_XLIM[2]
)

km_os_2group <- make_km_plot(
  data = km_os_plot_data,
  time_col = "OS_Months",
  event_col = "OS_event",
  group_col = "ecDNA_Status",
  palette = ECDNA_2GROUP_KM_COLS,
  x_label = "OS (months)",
  y_label = "Survival probability",
  title_label = "OS by ecDNA status"
)

logrank_os_2group <- survdiff(Surv(OS_Months, OS_event) ~ ecDNA_Status, data = km_os_data)

km_latency_data <- clinical_analysis %>%
  filter(!is.na(Latency_Months), !is.na(ecDNA_Status)) %>%
  mutate(Latency_event = 1)

km_latency_plot_data <- truncate_survival_time(
  data = km_latency_data,
  time_col = "Latency_Months",
  event_col = "Latency_event",
  max_time = KM_XLIM[2]
)

km_latency_2group <- make_km_plot(
  data = km_latency_plot_data,
  time_col = "Latency_Months",
  event_col = "Latency_event",
  group_col = "ecDNA_Status",
  palette = ECDNA_2GROUP_KM_COLS,
  x_label = "Time to brain metastasis (months)",
  y_label = "Event-free proportion",
  title_label = "Latency by ecDNA status"
)

logrank_latency_2group <- survdiff(Surv(Latency_Months, Latency_event) ~ ecDNA_Status, data = km_latency_data)

km_os_3group_data <- clinical_analysis %>%
  filter(!is.na(OS_Months), !is.na(OS_event), !is.na(ecDNA_Status_3group))

km_os_3group_plot_data <- truncate_survival_time(
  data = km_os_3group_data,
  time_col = "OS_Months",
  event_col = "OS_event",
  max_time = KM_XLIM[2]
)

km_os_3group <- make_km_plot(
  data = km_os_3group_plot_data,
  time_col = "OS_Months",
  event_col = "OS_event",
  group_col = "ecDNA_Status_3group",
  palette = ECDNA_3GROUP_KM_COLS,
  x_label = "OS (months)",
  y_label = "Survival probability",
  title_label = "OS by ecDNA category"
)

logrank_os_3group <- survdiff(Surv(OS_Months, OS_event) ~ ecDNA_Status_3group, data = km_os_3group_data)

km_latency_3group_data <- clinical_analysis %>%
  filter(!is.na(Latency_Months), !is.na(ecDNA_Status_3group)) %>%
  mutate(Latency_event = 1)

km_latency_3group_plot_data <- truncate_survival_time(
  data = km_latency_3group_data,
  time_col = "Latency_Months",
  event_col = "Latency_event",
  max_time = KM_XLIM[2]
)

km_latency_3group <- make_km_plot(
  data = km_latency_3group_plot_data,
  time_col = "Latency_Months",
  event_col = "Latency_event",
  group_col = "ecDNA_Status_3group",
  palette = ECDNA_3GROUP_KM_COLS,
  x_label = "Time to brain metastasis (months)",
  y_label = "Event-free proportion",
  title_label = "Latency by ecDNA category"
)

logrank_latency_3group <- survdiff(
  Surv(Latency_Months, Latency_event) ~ ecDNA_Status_3group,
  data = km_latency_3group_data
)

sink(file.path(DIR_TABLE, "figure5_logrank_tests.txt"))
cat("=== OS: ecDNA 2-group ===\n")
print(logrank_os_2group)
cat("\n=== Latency: ecDNA 2-group ===\n")
print(logrank_latency_2group)
cat("\n=== OS: ecDNA 3-group ===\n")
print(logrank_os_3group)
cat("\n=== Latency: ecDNA 3-group ===\n")
print(logrank_latency_3group)
sink()

ggsave(file.path(DIR_PLOT, "figure5_KM_OS_2group_curve.pdf"), km_os_2group$plot, width = 6.2, height = 5.2)
ggsave(file.path(DIR_PLOT, "figure5_KM_OS_2group_curve.png"), km_os_2group$plot, width = 6.2, height = 5.2, dpi = 300)
ggsave(file.path(DIR_PLOT, "figure5_KM_OS_2group_risktable.pdf"), km_os_2group$table, width = 6.2, height = 2.2)

ggsave(file.path(DIR_PLOT, "figure5_KM_latency_2group_curve.pdf"), km_latency_2group$plot, width = 6.2, height = 5.2)
ggsave(file.path(DIR_PLOT, "figure5_KM_latency_2group_curve.png"), km_latency_2group$plot, width = 6.2, height = 5.2, dpi = 300)
ggsave(file.path(DIR_PLOT, "figure5_KM_latency_2group_risktable.pdf"), km_latency_2group$table, width = 6.2, height = 2.2)

ggsave(file.path(DIR_PLOT, "figure5_KM_OS_3group_curve.pdf"), km_os_3group$plot, width = 6.8, height = 5.2)
ggsave(file.path(DIR_PLOT, "figure5_KM_OS_3group_curve.png"), km_os_3group$plot, width = 6.8, height = 5.2, dpi = 300)
ggsave(file.path(DIR_PLOT, "figure5_KM_OS_3group_risktable.pdf"), km_os_3group$table, width = 6.8, height = 2.4)

ggsave(file.path(DIR_PLOT, "figure5_KM_latency_3group_curve.pdf"), km_latency_3group$plot, width = 6.8, height = 5.2)
ggsave(file.path(DIR_PLOT, "figure5_KM_latency_3group_curve.png"), km_latency_3group$plot, width = 6.8, height = 5.2, dpi = 300)
ggsave(file.path(DIR_PLOT, "figure5_KM_latency_3group_risktable.pdf"), km_latency_3group$table, width = 6.8, height = 2.4)

# ============================================================
# 6. Cox models
# ============================================================

tidy_cox_model <- function(model, model_name) {
  broom::tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(model = model_name) %>%
    rename(
      HR = estimate,
      CI_low = conf.low,
      CI_high = conf.high,
      P_value = p.value
    )
}

os_cox_data <- clinical_analysis %>%
  filter(!is.na(OS_Months), !is.na(OS_event), !is.na(ecDNA_Status), !is.na(CancerTypeCox)) %>%
  mutate(
    ecDNA_Status = factor(ecDNA_Status, levels = c("ecDNA-", "ecDNA+")),
    CancerType = CancerTypeCox,
    Sex = factor(Sex, levels = c("male", "female"))
  )

latency_cox_data <- clinical_analysis %>%
  filter(!is.na(Latency_Months), !is.na(ecDNA_Status), !is.na(CancerTypeCox)) %>%
  mutate(
    Latency_event = 1,
    ecDNA_Status = factor(ecDNA_Status, levels = c("ecDNA-", "ecDNA+")),
    CancerType = CancerTypeCox,
    Sex = factor(Sex, levels = c("male", "female"))
  )

os_clinical_cox_data <- os_cox_data %>%
  filter(!is.na(Sex))

latency_clinical_cox_data <- latency_cox_data %>%
  filter(!is.na(Sex))

cox_os_ecDNA <- coxph(Surv(OS_Months, OS_event) ~ ecDNA_Status, data = os_cox_data)
cox_os_cancer <- coxph(Surv(OS_Months, OS_event) ~ CancerType, data = os_cox_data)
cox_os_ecDNA_cancer <- coxph(Surv(OS_Months, OS_event) ~ ecDNA_Status + CancerType, data = os_cox_data)
cox_os_ecDNA_cancer_sex <- coxph(
  Surv(OS_Months, OS_event) ~ ecDNA_Status + CancerType + Sex,
  data = os_clinical_cox_data
)

cox_latency_ecDNA <- coxph(Surv(Latency_Months, Latency_event) ~ ecDNA_Status, data = latency_cox_data)
cox_latency_cancer <- coxph(Surv(Latency_Months, Latency_event) ~ CancerType, data = latency_cox_data)
cox_latency_ecDNA_cancer <- coxph(
  Surv(Latency_Months, Latency_event) ~ ecDNA_Status + CancerType,
  data = latency_cox_data
)
cox_latency_ecDNA_cancer_sex <- coxph(
  Surv(Latency_Months, Latency_event) ~ ecDNA_Status + CancerType + Sex,
  data = latency_clinical_cox_data
)

cox_results <- bind_rows(
  tidy_cox_model(cox_os_ecDNA, "OS: univariable ecDNA"),
  tidy_cox_model(cox_os_cancer, "OS: univariable cancer type"),
  tidy_cox_model(cox_os_ecDNA_cancer, "OS: ecDNA + cancer type"),
  tidy_cox_model(cox_os_ecDNA_cancer_sex, "OS: ecDNA + cancer type + sex"),
  tidy_cox_model(cox_latency_ecDNA, "Latency: univariable ecDNA"),
  tidy_cox_model(cox_latency_cancer, "Latency: univariable cancer type"),
  tidy_cox_model(cox_latency_ecDNA_cancer, "Latency: ecDNA + cancer type"),
  tidy_cox_model(cox_latency_ecDNA_cancer_sex, "Latency: ecDNA + cancer type + sex")
)

cox_ecdna_results <- cox_results %>%
  filter(str_detect(term, "ecDNA_Status"))

cox_model_order <- c(
  "OS: univariable ecDNA",
  "OS: univariable cancer type",
  "OS: ecDNA + cancer type",
  "OS: ecDNA + cancer type + sex",
  "Latency: univariable ecDNA",
  "Latency: univariable cancer type",
  "Latency: ecDNA + cancer type",
  "Latency: ecDNA + cancer type + sex"
)

cox_term_order <- c(
  "ecDNA_StatusecDNA+",
  paste0("CancerType", c("CESC", "COAD", "ESCC", "KIRC", "LUAD", "Lung-PD", "LUSC", "OV", "SCLC", "SKCM", "STAD")),
  "Sexfemale"
)

format_cox_results_for_export <- function(cox_table) {
  cox_table %>%
    mutate(
      model = factor(model, levels = cox_model_order),
      term = factor(term, levels = cox_term_order),
      term = as.character(term),
      term = case_when(
        term == "ecDNA_StatusecDNA+" ~ "ecDNA Status ecDNA+",
        str_starts(term, "CancerType") ~ paste("Cancer type", str_remove(term, "^CancerType")),
        term == "Sexfemale" ~ "Sexfemale",
        TRUE ~ term
      )
    ) %>%
    arrange(model, match(term, c(
      "ecDNA Status ecDNA+",
      paste("Cancer type", c("CESC", "COAD", "ESCC", "KIRC", "LUAD", "Lung-PD", "LUSC", "OV", "SCLC", "SKCM", "STAD")),
      "Sexfemale"
    ))) %>%
    select(term, HR, std.error, statistic, P_value, CI_low, CI_high, model) %>%
    mutate(model = as.character(model))
}

cox_results_for_export <- format_cox_results_for_export(cox_results)
cox_ecdna_results_for_export <- format_cox_results_for_export(cox_ecdna_results)

write.csv(os_cox_data, file.path(DIR_TABLE, "figure5_OS_cox_input.csv"), row.names = FALSE)
write.csv(latency_cox_data, file.path(DIR_TABLE, "figure5_Latency_cox_input.csv"), row.names = FALSE)
write.csv(cox_results_for_export, file.path(DIR_TABLE, "figure5_cox_results_all.csv"), row.names = FALSE)
write.csv(cox_ecdna_results_for_export, file.path(DIR_TABLE, "figure5_cox_results_ecDNA_only.csv"), row.names = FALSE)

if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(
    list(figure5_cox_results_all = cox_results_for_export),
    file.path(DIR_TABLE, "figure5_cox_results_all.xlsx")
  )
}

sink(file.path(DIR_TABLE, "figure5_cox_summary.txt"))
cat("=== OS: ecDNA ===\n")
print(summary(cox_os_ecDNA))
cat("\n=== OS: cancer type ===\n")
print(summary(cox_os_cancer))
cat("\n=== OS: ecDNA + cancer type ===\n")
print(summary(cox_os_ecDNA_cancer))
cat("\n=== OS: ecDNA + cancer type + sex ===\n")
print(summary(cox_os_ecDNA_cancer_sex))
cat("\n=== Latency: ecDNA ===\n")
print(summary(cox_latency_ecDNA))
cat("\n=== Latency: cancer type ===\n")
print(summary(cox_latency_cancer))
cat("\n=== Latency: ecDNA + cancer type ===\n")
print(summary(cox_latency_ecDNA_cancer))
cat("\n=== Latency: ecDNA + cancer type + sex ===\n")
print(summary(cox_latency_ecDNA_cancer_sex))
cat("\n=== PH assumption: OS ecDNA ===\n")
print(cox.zph(cox_os_ecDNA))
cat("\n=== PH assumption: OS ecDNA + cancer type ===\n")
print(cox.zph(cox_os_ecDNA_cancer))
cat("\n=== PH assumption: Latency ecDNA ===\n")
print(cox.zph(cox_latency_ecDNA))
cat("\n=== PH assumption: Latency ecDNA + cancer type ===\n")
print(cox.zph(cox_latency_ecDNA_cancer))
sink()

# ============================================================
# 7. forest plots
# ============================================================

prepare_ecDNA_forest_table <- function(cox_table) {
  cox_table %>%
    mutate(
      model_label = model,
      hr_text = sprintf("%.2f (%.2f-%.2f)", HR, CI_low, CI_high),
      p_text = format_p_value(P_value),
      row_id = rev(seq_len(n()))
    )
}

make_ecDNA_forest <- function(cox_table, title_label) {
  forest_data <- prepare_ecDNA_forest_table(cox_table)

  table_panel <- ggplot(forest_data) +
    geom_text(aes(x = 0.0, y = row_id, label = model_label), hjust = 0, size = 4) +
    geom_text(aes(x = 6.0, y = row_id, label = hr_text), hjust = 0, size = 4) +
    geom_text(aes(x = 9.3, y = row_id, label = p_text), hjust = 0, size = 4) +
    annotate("text", x = 0.0, y = max(forest_data$row_id) + 0.8, label = title_label, fontface = "bold", hjust = 0, size = 6) +
    annotate("text", x = 0.4, y = max(forest_data$row_id) + 0.8, label = "Model", fontface = "bold", hjust = 0, size = 4.2) +
    annotate("text", x = 6.0, y = max(forest_data$row_id) + 0.8, label = "HR (95% CI)", fontface = "bold", hjust = 0, size = 4.2) +
    annotate("text", x = 9.3, y = max(forest_data$row_id) + 0.8, label = "P", fontface = "bold", hjust = 0, size = 4.2) +
    xlim(0, 10) +
    ylim(0.5, max(forest_data$row_id) + 1.2) +
    theme_void()

  forest_panel <- ggplot(forest_data, aes(x = HR, y = row_id)) +
    geom_vline(xintercept = 1, linetype = 2, color = "grey50") +
    geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.15, color = "#8B0000") +
    geom_point(size = 3, color = "#8B0000") +
    scale_x_log10() +
    labs(x = "Hazard ratio for ecDNA+ (log scale)", y = NULL) +
    ylim(0.5, max(forest_data$row_id) + 1.2) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank()
    )

  table_panel + forest_panel + plot_layout(widths = c(2.4, 1.4))
}

ecdna_forest <- make_ecDNA_forest(cox_ecdna_results, title_label = "ecDNA Cox models")
save_plot_pdf_png(ecdna_forest, "figure5_ecDNA_cox_forest", width = 13, height = 5.2)

build_full_forest_table <- function(model_table, model_data, reference_cancer_type) {
  ecdna_counts <- model_data %>%
    count(ecDNA_Status, name = "n")

  sex_counts <- model_data %>%
    count(Sex, name = "n")

  cancer_counts <- model_data %>%
    count(CancerType, name = "n")

  ecDNA_rows <- tibble(
    section = "ecDNA status",
    level = c("ecDNA-", "ecDNA+"),
    n = ecdna_counts$n[match(c("ecDNA-", "ecDNA+"), ecdna_counts$ecDNA_Status)],
    term = c(NA_character_, "ecDNA_StatusecDNA+"),
    is_reference = c(TRUE, FALSE)
  )

  sex_rows <- tibble(
    section = "Sex",
    level = c("male", "female"),
    n = sex_counts$n[match(c("male", "female"), sex_counts$Sex)],
    term = c(NA_character_, "Sexfemale"),
    is_reference = c(TRUE, FALSE)
  )

  age_row <- tibble(
    section = "Age",
    level = "per 1 SD increase",
    n = sum(!is.na(model_data$Age_z)),
    term = "Age_z",
    is_reference = FALSE
  )

  cancer_term_rows <- model_table %>%
    filter(str_detect(term, "^CancerType")) %>%
    transmute(
      section = "Cancer type",
      level = str_remove(term, "^CancerType"),
      n = cancer_counts$n[match(level, cancer_counts$CancerType)],
      term = term,
      is_reference = FALSE
    )

  cancer_reference_row <- tibble(
    section = "Cancer type",
    level = reference_cancer_type,
    n = cancer_counts$n[match(reference_cancer_type, cancer_counts$CancerType)],
    term = NA_character_,
    is_reference = TRUE
  )

  forest_table <- bind_rows(ecDNA_rows, sex_rows, age_row, cancer_reference_row, cancer_term_rows)

  forest_table %>%
    left_join(model_table, by = "term") %>%
    mutate(
      hr_text = case_when(
        is_reference ~ "Reference",
        is.na(HR) ~ "NA",
        TRUE ~ sprintf("%.2f (%.2f-%.2f)", HR, CI_low, CI_high)
      ),
      p_text = case_when(
        is_reference ~ "",
        is.na(P_value) ~ "NA",
        TRUE ~ format_p_value(P_value)
      ),
      row_id = rev(seq_len(n()))
    )
}

make_full_forest_plot <- function(forest_table, title_label, footer_text) {
  plot_values <- c(forest_table$CI_low, forest_table$CI_high)
  plot_values <- plot_values[is.finite(plot_values) & !is.na(plot_values) & plot_values > 0]

  x_limits <- c(0.05, 20)
  if (length(plot_values) > 0) {
    x_limits[1] <- max(0.02, 10^floor(log10(min(plot_values) / 1.2)))
    x_limits[2] <- min(50, 10^ceiling(log10(max(plot_values) * 1.2)))
  }

  x_breaks <- c(0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50)
  x_breaks <- x_breaks[x_breaks >= x_limits[1] & x_breaks <= x_limits[2]]

  background_data <- forest_table %>%
    mutate(
      ymin = row_id - 0.5,
      ymax = row_id + 0.5,
      fill_color = rep(c("#D9D9D9", "#BFBFBF"), length.out = n())
    )

  table_panel <- ggplot() +
    geom_rect(
      data = background_data,
      aes(xmin = -0.2, xmax = 11.6, ymin = ymin, ymax = ymax, fill = fill_color),
      color = NA
    ) +
    scale_fill_identity() +
    geom_text(data = forest_table, aes(x = 0.00, y = row_id, label = section), hjust = 0, size = 4) +
    geom_text(data = forest_table, aes(x = 1.75, y = row_id, label = level), hjust = 0, size = 4) +
    geom_text(data = forest_table, aes(x = 4.55, y = row_id, label = ifelse(is.na(n), "", paste0("(n = ", n, ")"))), hjust = 0, size = 4) +
    geom_text(data = forest_table, aes(x = 6.25, y = row_id, label = hr_text), hjust = 0, size = 4) +
    geom_text(data = forest_table, aes(x = 10.15, y = row_id, label = p_text), hjust = 0, size = 4) +
    annotate("text", x = -0.1, y = max(forest_table$row_id) + 0.9, label = title_label, fontface = "bold", hjust = 0, size = 6) +
    annotate("text", x = 0, y = 0.05, label = footer_text, hjust = 0, vjust = 0, size = 4, fontface = "italic") +
    xlim(-0.2, 11.6) +
    ylim(0.5, max(forest_table$row_id) + 1.2) +
    theme_void()

  plot_data <- forest_table %>%
    filter(!is_reference, !is.na(HR), !is.na(CI_low), !is.na(CI_high)) %>%
    mutate(
      HR_plot = pmin(pmax(HR, x_limits[1]), x_limits[2]),
      CI_low_plot = pmax(CI_low, x_limits[1]),
      CI_high_plot = pmin(CI_high, x_limits[2])
    )

  forest_panel <- ggplot() +
    geom_rect(
      data = background_data,
      aes(xmin = x_limits[1], xmax = x_limits[2], ymin = ymin, ymax = ymax, fill = fill_color),
      color = NA
    ) +
    scale_fill_identity() +
    geom_vline(xintercept = 1, linetype = 2, linewidth = 0.7) +
    geom_errorbarh(
      data = plot_data,
      aes(y = row_id, xmin = CI_low_plot, xmax = CI_high_plot),
      height = 0.15,
      linewidth = 1,
      color = "black"
    ) +
    geom_point(data = plot_data, aes(x = HR_plot, y = row_id), shape = 15, size = 3, color = "black") +
    scale_x_log10(breaks = x_breaks, labels = x_breaks) +
    coord_cartesian(xlim = x_limits, clip = "off") +
    ylim(0.5, max(forest_table$row_id) + 1.2) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      plot.margin = margin(5.5, 15, 5.5, 0)
    )

  table_panel + forest_panel + plot_layout(widths = c(2.6, 1.6))
}

reference_cancer_type <- os_clinical_cox_data %>%
  count(CancerType, sort = TRUE) %>%
  slice(1) %>%
  pull(CancerType) %>%
  as.character()

os_full_forest_data <- os_clinical_cox_data %>%
  mutate(CancerType = relevel(CancerType, ref = reference_cancer_type))

cox_os_full_forest <- coxph(
  Surv(OS_Months, OS_event) ~ ecDNA_Status + Sex + Age_z + CancerType,
  data = os_full_forest_data
)

os_full_model_table <- tidy_cox_model(cox_os_full_forest, "OS: ecDNA + sex + age + cancer type")

os_full_forest_table <- build_full_forest_table(
  model_table = os_full_model_table,
  model_data = os_full_forest_data,
  reference_cancer_type = reference_cancer_type
)

os_full_forest_footer <- paste0(
  "Events: ", sum(os_full_forest_data$OS_event, na.rm = TRUE),
  "; N = ", nrow(os_full_forest_data),
  "; Reference cancer type: ", reference_cancer_type
)

os_full_forest <- make_full_forest_plot(
  forest_table = os_full_forest_table,
  title_label = "OS full model",
  footer_text = os_full_forest_footer
)

write.csv(os_full_model_table, file.path(DIR_TABLE, "figure5_OS_full_model.csv"), row.names = FALSE)
write.csv(os_full_forest_table, file.path(DIR_TABLE, "figure5_OS_full_forest_table.csv"), row.names = FALSE)
save_plot_pdf_png(os_full_forest, "figure5_OS_full_model_forest", width = 16, height = max(6, 0.48 * nrow(os_full_forest_table)))

# ============================================================
# 8. cancer-type-specific Cox analyses
# ============================================================

run_within_cancer_cox <- function(data, time_col, event_col) {
  cancer_types <- sort(unique(data$CancerType))
  result_list <- vector("list", length(cancer_types))

  for (i in seq_along(cancer_types)) {
    current_cancer_type <- cancer_types[i]

    model_data <- data %>%
      filter(CancerType == current_cancer_type) %>%
      filter(!is.na(.data[[time_col]]), !is.na(.data[[event_col]]), !is.na(ecDNA_Status)) %>%
      mutate(ecDNA_Status = factor(ecDNA_Status, levels = c("ecDNA-", "ecDNA+")))

    n_total <- nrow(model_data)
    n_pos <- sum(model_data$ecDNA_Status == "ecDNA+", na.rm = TRUE)
    n_neg <- sum(model_data$ecDNA_Status == "ecDNA-", na.rm = TRUE)
    n_event <- sum(model_data[[event_col]] == 1, na.rm = TRUE)

    base_row <- tibble(
      CancerType = current_cancer_type,
      n = n_total,
      n_ecDNA_pos = n_pos,
      n_ecDNA_neg = n_neg,
      n_event = n_event,
      HR = NA_real_,
      CI_low = NA_real_,
      CI_high = NA_real_,
      P_value = NA_real_,
      note = "insufficient data"
    )

    has_enough_data <- n_total >= MIN_N_PER_CANCER_FOR_KM &&
      n_pos >= MIN_ECDNA_POS_FOR_COX &&
      n_neg >= MIN_ECDNA_POS_FOR_COX &&
      n_event >= MIN_EVENTS_FOR_COX

    if (!has_enough_data) {
      result_list[[i]] <- base_row
      next
    }

    model_formula <- as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ ecDNA_Status"))
    model_fit <- tryCatch(coxph(model_formula, data = model_data), error = function(e) NULL)

    if (is.null(model_fit)) {
      base_row$note <- "model failed"
      result_list[[i]] <- base_row
      next
    }

    model_tidy <- broom::tidy(model_fit, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(str_detect(term, "ecDNA_Status"))

    if (nrow(model_tidy) == 0) {
      base_row$note <- "no ecDNA term"
      result_list[[i]] <- base_row
      next
    }

    result_list[[i]] <- tibble(
      CancerType = current_cancer_type,
      n = n_total,
      n_ecDNA_pos = n_pos,
      n_ecDNA_neg = n_neg,
      n_event = n_event,
      HR = model_tidy$estimate[1],
      CI_low = model_tidy$conf.low[1],
      CI_high = model_tidy$conf.high[1],
      P_value = model_tidy$p.value[1],
      note = "ok"
    )
  }

  bind_rows(result_list) %>%
    arrange(P_value, desc(HR), CancerType)
}

within_cancer_os <- run_within_cancer_cox(os_cox_data, "OS_Months", "OS_event")
within_cancer_latency <- run_within_cancer_cox(latency_cox_data, "Latency_Months", "Latency_event")

write.csv(within_cancer_os, file.path(DIR_TABLE, "figure5_OS_ecDNA_cox_within_cancer_type.csv"), row.names = FALSE)
write.csv(within_cancer_latency, file.path(DIR_TABLE, "figure5_Latency_ecDNA_cox_within_cancer_type.csv"), row.names = FALSE)

# ============================================================
# 9. ecDNA status by cancer type
# ============================================================

ecdna_cancer_table <- table(clinical_analysis$CancerType, clinical_analysis$ecDNA_Status)
chisq_ecdna_cancer <- suppressWarnings(chisq.test(ecdna_cancer_table))
fisher_ecdna_cancer <- tryCatch(fisher.test(ecdna_cancer_table), error = function(e) NULL)

write.csv(
  as.data.frame.matrix(ecdna_cancer_table),
  file.path(DIR_TABLE, "figure5_ecDNA_by_cancer_type_count_table.csv")
)

sink(file.path(DIR_TABLE, "figure5_ecDNA_vs_cancer_type_test.txt"))
cat("=== CancerType x ecDNA_Status count table ===\n")
print(ecdna_cancer_table)
cat("\n=== Chi-square test ===\n")
print(chisq_ecdna_cancer)
cat("\n=== Fisher exact test ===\n")
print(fisher_ecdna_cancer)
sink()

# ============================================================
# 10. cancer-level ecDNA prevalence vs survival summaries
# ============================================================

safe_spearman_test <- function(x, y) {
  keep <- is.finite(x) & is.finite(y) & !is.na(x) & !is.na(y)
  x_keep <- x[keep]
  y_keep <- y[keep]

  if (length(x_keep) < 3) {
    return(tibble(method = "spearman", estimate = NA_real_, p.value = NA_real_, n = length(x_keep)))
  }

  test_out <- suppressWarnings(cor.test(x_keep, y_keep, method = "spearman", exact = FALSE))

  tibble(
    method = "spearman",
    estimate = unname(test_out$estimate),
    p.value = test_out$p.value,
    n = length(x_keep)
  )
}

cor_label <- function(cor_table) {
  if (nrow(cor_table) == 0 || is.na(cor_table$estimate[1])) {
    return("Spearman rho = NA, p = NA")
  }

  paste0(
    "Spearman rho = ", sprintf("%.2f", cor_table$estimate[1]),
    ", p = ", format_p_value(cor_table$p.value[1]),
    ", n = ", cor_table$n[1]
  )
}

cancer_level_summary <- clinical_analysis %>%
  filter(!is.na(CancerType), !is.na(ecDNA_Status)) %>%
  group_by(CancerType) %>%
  summarise(
    n_total = n(),
    n_ecDNA_pos = sum(ecDNA_Status == "ecDNA+", na.rm = TRUE),
    ecDNA_pct = 100 * n_ecDNA_pos / n_total,
    median_OS = safe_median(OS_Months),
    median_Latency = safe_median(Latency_Months),
    n_os = sum(!is.na(OS_Months)),
    n_latency = sum(!is.na(Latency_Months)),
    .groups = "drop"
  ) %>%
  filter(n_total >= 5)

cor_ecDNA_pct_os <- safe_spearman_test(cancer_level_summary$ecDNA_pct, cancer_level_summary$median_OS)
cor_ecDNA_pct_latency <- safe_spearman_test(cancer_level_summary$ecDNA_pct, cancer_level_summary$median_Latency)

write.csv(cancer_level_summary, file.path(DIR_TABLE, "figure5_cancer_type_ecDNA_prevalence_vs_survival.csv"), row.names = FALSE)
write.csv(cor_ecDNA_pct_os, file.path(DIR_TABLE, "figure5_cor_ecDNApct_vs_medianOS.csv"), row.names = FALSE)
write.csv(cor_ecDNA_pct_latency, file.path(DIR_TABLE, "figure5_cor_ecDNApct_vs_medianLatency.csv"), row.names = FALSE)

ecDNA_pct_os_plot <- ggplot(cancer_level_summary, aes(x = ecDNA_pct, y = median_OS, label = CancerType)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7, linetype = 2) +
  geom_text_repel(size = 3.5, max.overlaps = 30) +
  labs(
    title = "Cancer-type ecDNA prevalence vs median OS",
    x = "ecDNA-positive cases by cancer type (%)",
    y = "Median OS by cancer type (months)",
    subtitle = cor_label(cor_ecDNA_pct_os)
  ) +
  theme_classic(base_size = 12)

ecDNA_pct_latency_plot <- ggplot(cancer_level_summary, aes(x = ecDNA_pct, y = median_Latency, label = CancerType)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7, linetype = 2) +
  geom_text_repel(size = 3.5, max.overlaps = 30) +
  labs(
    title = "Cancer-type ecDNA prevalence vs median latency",
    x = "ecDNA-positive cases by cancer type (%)",
    y = "Median latency by cancer type (months)",
    subtitle = cor_label(cor_ecDNA_pct_latency)
  ) +
  theme_classic(base_size = 12)

save_plot_pdf_png(ecDNA_pct_os_plot, "figure5_OS_ecDNA_spearman_by_cancer", width = 6.5, height = 5.2)
save_plot_pdf_png(ecDNA_pct_latency_plot, "figure5_latency_ecDNA_spearman_by_cancer", width = 6.5, height = 5.2)

# ============================================================
# 11. WGS oncogene-group figure
# ============================================================

make_oncogene_km_panel <- function(data, time_col, event_col, x_label, y_label, p_coord) {
  plot_data <- data %>%
    filter(!is.na(.data[[time_col]]), !is.na(.data[[event_col]]), !is.na(OncogeneAmplicon_Group)) %>%
    mutate(
      time_plot = .data[[time_col]],
      event_plot = .data[[event_col]],
      group_plot = OncogeneAmplicon_Group
    )

  plot_data <- truncate_survival_time(
    data = plot_data,
    time_col = "time_plot",
    event_col = "event_plot",
    max_time = KM_XLIM[2]
  )

  fit <- survfit(Surv(time_plot, event_plot) ~ group_plot, data = plot_data)
  logrank <- survdiff(Surv(time_plot, event_plot) ~ group_plot, data = plot_data)
  p_value <- 1 - pchisq(logrank$chisq, df = length(logrank$n) - 1)

  km <- ggsurvplot(
    fit,
    data = plot_data,
    risk.table = TRUE,
    pval = format_plot_p_value(p_value),
    pval.coord = p_coord,
    conf.int = FALSE,
    palette = unname(ONCOGENE_3GROUP_COLS),
    xlab = x_label,
    ylab = y_label,
    legend.title = NULL,
    legend.labs = levels(plot_data$group_plot),
    xlim = KM_XLIM,
    break.time.by = KM_BREAK,
    risk.table.height = 0.28,
    risk.table.y.text.col = TRUE,
    ggtheme = theme_classic(base_size = 15),
    tables.theme = theme_cleantable(base_size = 14)
  )

  km$plot <- km$plot +
    guides(color = guide_legend(title = NULL, ncol = 1)) +
    theme(
      aspect.ratio = 1,
      legend.position = c(0.28, 0.18),
      legend.background = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black", size = 16),
      axis.line = element_line(linewidth = 0.8, color = "black"),
      axis.ticks = element_line(linewidth = 0.7, color = "black"),
      plot.margin = margin(4, 8, 0, 8)
    )

  km$table <- km$table +
    labs(x = x_label) +
    theme(
      plot.title = element_text(size = 14, hjust = 0),
      axis.title.x = element_text(size = 14, color = "black"),
      axis.text.x = element_text(size = 13, color = "black"),
      plot.margin = margin(0, 8, 4, 8)
    )

  list(
    panel = ggpubr::ggarrange(km$plot, km$table, ncol = 1, heights = c(3.2, 1.12), align = "v"),
    plot = km$plot,
    table = km$table,
    p_value = p_value,
    input = plot_data
  )
}

make_subtype_box_panel <- function(data, y_col, y_label) {
  subtype_order <- data %>%
    count(CancerSubtypeFigure, name = "n") %>%
    arrange(CancerSubtypeFigure) %>%
    pull(CancerSubtypeFigure)

  plot_data <- data %>%
    filter(!is.na(.data[[y_col]]), !is.na(ecDNA_Status), !is.na(CancerSubtypeFigure)) %>%
    group_by(CancerSubtypeFigure) %>%
    filter(n_distinct(ecDNA_Status) == 2) %>%
    ungroup() %>%
    mutate(
      CancerSubtypeFigure = factor(CancerSubtypeFigure, levels = subtype_order),
      ecDNA_Status = factor(ecDNA_Status, levels = c("ecDNA-", "ecDNA+"))
    )

  p_labels <- plot_data %>%
    group_by(CancerSubtypeFigure) %>%
    group_modify(~ {
      p_value <- wilcox_ecdna_negative_greater(.x, y_col)
      tibble(
        p_value = p_value,
        label = paste0("p = ", ifelse(is.na(p_value), "NA", sprintf("%.2f", p_value))),
        y_position = max(.x[[y_col]], na.rm = TRUE) * 0.96
      )
    }) %>%
    ungroup()

  ggplot(plot_data, aes(x = ecDNA_Status, y = .data[[y_col]], fill = ecDNA_Status)) +
    stat_boxplot(geom = "errorbar", width = 0.22, linewidth = 0.55, color = "black") +
    geom_boxplot(width = 0.52, linewidth = 0.55, outlier.shape = NA, alpha = 0.95, color = "black") +
    geom_jitter(width = 0.12, height = 0, size = 0.9, alpha = 0.9, color = "#1B2A3A") +
    geom_text(
      data = p_labels,
      aes(x = 1.5, y = y_position, label = label),
      inherit.aes = FALSE,
      size = 3.5
    ) +
    scale_fill_manual(values = ECDNA_2GROUP_COLS, drop = FALSE) +
    facet_wrap(~ CancerSubtypeFigure, nrow = 1, scales = "free_x", strip.position = "top") +
    labs(x = NULL, y = y_label) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 12),
      axis.title.y = element_text(color = "black", size = 15),
      axis.line = element_line(linewidth = 0.7, color = "black"),
      strip.background = element_rect(color = "black", fill = "white", linewidth = 0.7),
      strip.text = element_text(color = "black", size = 14),
      panel.spacing.x = unit(0.12, "lines"),
      plot.margin = margin(4, 8, 4, 8)
    )
}

oncogene_latency_data <- clinical_analysis %>%
  filter(!is.na(Latency_Months), !is.na(OncogeneAmplicon_Group)) %>%
  mutate(Latency_event = 1)

oncogene_os_data <- clinical_analysis %>%
  filter(!is.na(OS_Months), !is.na(OS_event), !is.na(OncogeneAmplicon_Group))

oncogene_latency_km <- make_oncogene_km_panel(
  data = oncogene_latency_data,
  time_col = "Latency_Months",
  event_col = "Latency_event",
  x_label = "Time to brain metastasis (Months)",
  y_label = "Event-free proportion",
  p_coord = c(36, 0.44)
)

oncogene_os_km <- make_oncogene_km_panel(
  data = oncogene_os_data,
  time_col = "OS_Months",
  event_col = "OS_event",
  x_label = "OS (Months)",
  y_label = "Overall survival probability",
  p_coord = c(22, 0.46)
)

subtype_latency_box <- make_subtype_box_panel(
  data = clinical_analysis,
  y_col = "Latency_Months",
  y_label = "Time to brain metastasis (Months)"
)

subtype_os_box <- make_subtype_box_panel(
  data = clinical_analysis,
  y_col = "OS_Months",
  y_label = "OS (Months)"
)

figure5_wgs_oncogene_panels <- (
  (oncogene_latency_km$panel | oncogene_os_km$panel) /
    subtype_latency_box /
    subtype_os_box
) +
  plot_layout(heights = c(1.25, 0.78, 0.78)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "plain", size = 20))

write.csv(
  clinical_analysis %>%
    count(OncogeneAmplicon_Group, ecDNA_status_raw, name = "n"),
  file.path(DIR_TABLE, "figure5_WGS_oncogene_group_definition_counts.csv"),
  row.names = FALSE
)

write.csv(
  bind_rows(
    tibble(endpoint = "Latency", p_value = oncogene_latency_km$p_value),
    tibble(endpoint = "OS", p_value = oncogene_os_km$p_value)
  ),
  file.path(DIR_TABLE, "figure5_WGS_oncogene_group_logrank.csv"),
  row.names = FALSE
)

save_plot_pdf_png(
  figure5_wgs_oncogene_panels,
  "figure5_WGS_oncogene_panels",
  width = 15.5,
  height = 13.2
)

make_canonical_oncogene_km_panel <- function(data, time_col, event_col, x_label, y_label, title_label, p_coord) {
  plot_data <- data %>%
    filter(!is.na(.data[[time_col]]), !is.na(.data[[event_col]]), !is.na(CanonicalOncogene_Group)) %>%
    mutate(
      time_plot = .data[[time_col]],
      event_plot = .data[[event_col]],
      group_plot = CanonicalOncogene_Group
    )

  plot_data <- truncate_survival_time(
    data = plot_data,
    time_col = "time_plot",
    event_col = "event_plot",
    max_time = KM_XLIM[2]
  )

  fit <- survfit(Surv(time_plot, event_plot) ~ group_plot, data = plot_data)
  logrank <- survdiff(Surv(time_plot, event_plot) ~ group_plot, data = plot_data)
  p_value <- 1 - pchisq(logrank$chisq, df = length(logrank$n) - 1)

  km <- ggsurvplot(
    fit,
    data = plot_data,
    risk.table = TRUE,
    pval = format_plot_p_value(p_value),
    pval.coord = p_coord,
    conf.int = FALSE,
    palette = unname(CANONICAL_ONCOGENE_3GROUP_COLS),
    xlab = x_label,
    ylab = y_label,
    title = title_label,
    legend.title = "strata",
    legend.labs = levels(plot_data$group_plot),
    xlim = KM_XLIM,
    break.time.by = KM_BREAK,
    risk.table.height = 0.28,
    risk.table.y.text.col = TRUE,
    ggtheme = theme_classic(base_size = 14),
    tables.theme = theme_cleantable(base_size = 13)
  )

  km$plot <- km$plot +
    guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
      legend.position = "top",
      legend.justification = "center",
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 11),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black", size = 14),
      axis.line = element_line(linewidth = 0.7, color = "black"),
      axis.ticks = element_line(linewidth = 0.6, color = "black"),
      plot.margin = margin(4, 8, 0, 8)
    )

  km$table <- km$table +
    labs(x = x_label) +
    theme(
      plot.title = element_text(size = 14, hjust = 0),
      axis.title.x = element_text(size = 13, color = "black"),
      axis.text.x = element_text(size = 12, color = "black"),
      plot.margin = margin(0, 8, 4, 8)
    )

  list(
    panel = ggpubr::ggarrange(km$plot, km$table, ncol = 1, heights = c(3.4, 1.05), align = "v"),
    p_value = p_value,
    input = plot_data
  )
}

canonical_os_km <- make_canonical_oncogene_km_panel(
  data = clinical_analysis,
  time_col = "OS_Months",
  event_col = "OS_event",
  x_label = "OS (Months)",
  y_label = "Survival probability",
  title_label = "OS by ecDNA / oncogene amplicon status",
  p_coord = c(1.5, 0.18)
)

canonical_latency_km <- make_canonical_oncogene_km_panel(
  data = clinical_analysis %>% mutate(Latency_event = 1),
  time_col = "Latency_Months",
  event_col = "Latency_event",
  x_label = "Time to brain metastasis (Months)",
  y_label = "Event-free proportion",
  title_label = "Latency by ecDNA / oncogene amplicon status",
  p_coord = c(1.5, 0.18)
)

figure5_canonical_oncogene_km <- ggpubr::ggarrange(
  canonical_os_km$panel,
  canonical_latency_km$panel,
  ncol = 2,
  widths = c(1, 1),
  align = "hv"
)

write.csv(
  clinical_analysis %>%
    count(CanonicalOncogene_Group, ecDNA_status_raw, name = "n"),
  file.path(DIR_TABLE, "figure5_canonical_oncogene_group_definition_counts.csv"),
  row.names = FALSE
)

write.csv(
  bind_rows(
    tibble(endpoint = "OS", p_value = canonical_os_km$p_value),
    tibble(endpoint = "Latency", p_value = canonical_latency_km$p_value)
  ),
  file.path(DIR_TABLE, "figure5_canonical_oncogene_group_logrank.csv"),
  row.names = FALSE
)

save_plot_pdf_png(
  figure5_canonical_oncogene_km,
  "figure5_canonical_oncogene_KM_with_risktables",
  width = 13.5,
  height = 6.4
)

# ============================================================
# 12. publication-style main figure
# ============================================================

publication_box_plot <- function(data, y_col, y_label, tag_label) {
  p_value <- wilcox_ecdna_negative_greater(data, y_col)
  p_label <- paste("Wilcoxon,", format_plot_p_value(p_value))
  y_max <- max(data[[y_col]], na.rm = TRUE)

  ggplot(data, aes(x = ecDNA_Status, y = .data[[y_col]], fill = ecDNA_Status)) +
    stat_boxplot(geom = "errorbar", width = 0.32, linewidth = 0.7, color = "black") +
    geom_boxplot(width = 0.52, linewidth = 0.7, outlier.shape = NA, alpha = 0.92, color = "black") +
    geom_jitter(
      width = 0.12,
      height = 0,
      size = 1.25,
      shape = 21,
      stroke = 0.35,
      alpha = 0.9,
      color = "#0F1B2B",
      fill = "#8EA0B8"
    ) +
    annotate("text", x = 1.82, y = y_max * 1.02, label = p_label, size = 5.0, hjust = 0.5) +
    scale_fill_manual(values = ECDNA_2GROUP_COLS, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.14))) +
    labs(x = NULL, y = y_label, tag = tag_label) +
    theme_classic(base_size = 16) +
    theme(
      legend.position = "none",
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(size = 15),
      axis.title.y = element_text(color = "black", size = 16),
      axis.line = element_line(linewidth = 0.8, color = "black"),
      axis.ticks = element_line(linewidth = 0.7, color = "black"),
      plot.tag = element_text(face = "plain", size = 20),
      plot.tag.position = c(-0.16, 1.04),
      plot.margin = margin(8, 8, 8, 8)
    )
}

publication_km_object <- function(data, time_col, event_col, group_col, x_label, y_label, p_value, tag_label) {
  km_data <- data %>%
    transmute(
      km_time = .data[[time_col]],
      km_event = .data[[event_col]],
      km_group = factor(.data[[group_col]], levels = c("ecDNA-", "ecDNA+"))
    )

  km_fit <- survfit(Surv(km_time, km_event) ~ km_group, data = km_data)
  p_label <- format_plot_p_value(p_value)

  km_plot <- ggsurvplot(
    km_fit,
    data = km_data,
    risk.table = TRUE,
    pval = p_label,
    pval.coord = c(1.5, 0.18),
    conf.int = FALSE,
    palette = unname(ECDNA_2GROUP_KM_COLS),
    xlab = x_label,
    ylab = y_label,
    legend.title = NULL,
    legend.labs = c("ecDNA-", "ecDNA+"),
    xlim = KM_XLIM,
    break.time.by = KM_BREAK,
    risk.table.height = 0.24,
    risk.table.y.text.col = TRUE,
    ggtheme = theme_classic(base_size = 15),
    tables.theme = theme_cleantable(base_size = 14)
  )

  km_plot$plot <- km_plot$plot +
    labs(color = NULL, tag = tag_label) +
    guides(color = guide_legend(title = NULL)) +
    theme(
      legend.position = "top",
      legend.justification = "center",
      legend.text = element_text(size = 13),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black", size = 16),
      axis.line = element_line(linewidth = 0.8, color = "black"),
      axis.ticks = element_line(linewidth = 0.7, color = "black"),
      plot.tag = element_text(face = "plain", size = 20),
      plot.tag.position = c(-0.16, 1.04),
      plot.margin = margin(4, 8, 0, 8)
    )

  km_plot$table <- km_plot$table +
    labs(x = x_label) +
    theme(
      plot.title = element_text(size = 14, hjust = 0),
      axis.title.x = element_text(size = 14, color = "black"),
      axis.text.x = element_text(size = 13, color = "black"),
      plot.margin = margin(0, 8, 4, 8)
    )

  km_plot$plot / km_plot$table + plot_layout(heights = c(3.2, 1.0))
}

format_forest_p_value <- function(p_value) {
  case_when(
    is.na(p_value) ~ "NA",
    p_value < 0.001 ~ "<0.001***",
    p_value < 0.01 ~ paste0(sprintf("%.3f", p_value), "**"),
    p_value < 0.05 ~ paste0(sprintf("%.3f", p_value), "*"),
    TRUE ~ sprintf("%.3f", p_value)
  )
}

make_publication_ecDNA_forest <- function(cox_table) {
  row_order <- c(
    "Latency: univariable ecDNA",
    "Latency: ecDNA + cancer type",
    "Latency: ecDNA + cancer type + sex",
    "OS: univariable ecDNA",
    "OS: ecDNA + cancer type",
    "OS: ecDNA + cancer type + sex"
  )

  forest_data <- cox_table %>%
    filter(model %in% row_order) %>%
    mutate(
      model = factor(model, levels = row_order),
      model_label = recode(
        as.character(model),
        "Latency: univariable ecDNA" = "Latency: univariable ecDNA",
        "Latency: ecDNA + cancer type" = "Latency: ecDNA + cancer type",
        "Latency: ecDNA + cancer type + sex" = "Latency: ecDNA + cancer type + sex",
        "OS: univariable ecDNA" = "OS: univariable ecDNA",
        "OS: ecDNA + cancer type" = "OS: ecDNA + cancer type",
        "OS: ecDNA + cancer type + sex" = "OS: ecDNA + cancer type + sex"
      ),
      hr_text = sprintf("%.2f (%.2f-%.2f)", HR, CI_low, CI_high),
      p_text = format_forest_p_value(P_value)
    ) %>%
    arrange(model) %>%
    mutate(row_id = rev(seq_len(n())))

  row_background <- forest_data %>%
    mutate(
      fill_color = rep(c("#D7D7D7", "white"), length.out = n()),
      ymin = row_id - 0.5,
      ymax = row_id + 0.5
    )

  table_panel <- ggplot() +
    geom_rect(
      data = row_background,
      aes(xmin = 0, xmax = 6.4, ymin = ymin, ymax = ymax, fill = fill_color),
      color = NA
    ) +
    scale_fill_identity() +
    annotate("rect", xmin = 0, xmax = 6.4, ymin = max(forest_data$row_id) + 0.45, ymax = max(forest_data$row_id) + 1.45, fill = "#D7D7D7") +
    annotate("text", x = 0.25, y = max(forest_data$row_id) + 0.95, label = "Model", fontface = "bold", hjust = 0, size = 5) +
    annotate("text", x = 4.15, y = max(forest_data$row_id) + 0.95, label = "HR (95% CI)", fontface = "bold", hjust = 0, size = 5) +
    geom_text(data = forest_data, aes(x = 0.25, y = row_id, label = model_label), hjust = 0, size = 4.1) +
    geom_text(data = forest_data, aes(x = 4.15, y = row_id, label = hr_text), hjust = 0, size = 4.1) +
    xlim(0, 6.4) +
    ylim(0.5, max(forest_data$row_id) + 1.45) +
    labs(tag = "E") +
    theme_void() +
    theme(plot.tag = element_text(face = "plain", size = 20), plot.tag.position = c(0, 1))

  forest_panel <- ggplot() +
    geom_rect(
      data = row_background,
      aes(xmin = 0.55, xmax = 7.0, ymin = ymin, ymax = ymax, fill = fill_color),
      color = NA
    ) +
    scale_fill_identity() +
    geom_vline(xintercept = 1, linetype = 2, linewidth = 0.7, color = "grey45") +
    geom_errorbarh(
      data = forest_data,
      aes(y = row_id, xmin = CI_low, xmax = CI_high),
      height = 0.12,
      linewidth = 0.75,
      color = "#8B2A25"
    ) +
    geom_point(data = forest_data, aes(x = HR, y = row_id), size = 3.2, color = "#8B2A25") +
    scale_x_log10(breaks = c(1, 3, 5), limits = c(0.55, 7.0), labels = c("1", "3", "5")) +
    labs(x = "Hazard ratio for ecDNA+ (log scale)", y = NULL) +
    ylim(0.5, max(forest_data$row_id) + 1.45) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      axis.title.x = element_text(size = 14),
      plot.margin = margin(0, 0, 0, 0)
    )

  p_panel <- ggplot() +
    geom_rect(
      data = row_background,
      aes(xmin = 0, xmax = 1, ymin = ymin, ymax = ymax, fill = fill_color),
      color = NA
    ) +
    scale_fill_identity() +
    annotate("rect", xmin = 0, xmax = 1, ymin = max(forest_data$row_id) + 0.45, ymax = max(forest_data$row_id) + 1.45, fill = "#D7D7D7") +
    annotate("text", x = 0.1, y = max(forest_data$row_id) + 0.95, label = "p value", fontface = "bold", hjust = 0, size = 5) +
    geom_text(data = forest_data, aes(x = 0.1, y = row_id, label = p_text), hjust = 0, size = 4.4) +
    xlim(0, 1) +
    ylim(0.5, max(forest_data$row_id) + 1.45) +
    theme_void()

  table_panel + forest_panel + p_panel + plot_layout(widths = c(2.25, 1.65, 0.65))
}

latency_box_pub <- publication_box_plot(
  data = latency_data,
  y_col = "Latency_Months",
  y_label = "Time to brain metastasis (months)",
  tag_label = "A"
)

os_box_pub <- publication_box_plot(
  data = os_data,
  y_col = "OS_Months",
  y_label = "OS (months)",
  tag_label = "B"
)

latency_km_p <- 1 - pchisq(logrank_latency_2group$chisq, df = length(logrank_latency_2group$n) - 1)
os_km_p <- 1 - pchisq(logrank_os_2group$chisq, df = length(logrank_os_2group$n) - 1)

latency_km_pub <- publication_km_object(
  data = km_latency_plot_data,
  time_col = "Latency_Months",
  event_col = "Latency_event",
  group_col = "ecDNA_Status",
  x_label = "Time to brain metastasis (Months)",
  y_label = "Event-free proportion",
  p_value = latency_km_p,
  tag_label = "C"
)

os_km_pub <- publication_km_object(
  data = km_os_plot_data,
  time_col = "OS_Months",
  event_col = "OS_event",
  group_col = "ecDNA_Status",
  x_label = "OS (Months)",
  y_label = "Overall survival probability",
  p_value = os_km_p,
  tag_label = "D"
)

publication_forest <- make_publication_ecDNA_forest(cox_ecdna_results)

figure5_main_publication_style <- (
  (latency_box_pub | os_box_pub) /
    (latency_km_pub | os_km_pub) /
    publication_forest
) +
  plot_layout(heights = c(1.05, 1.25, 1.25))

save_plot_pdf_png(
  figure5_main_publication_style,
  "figure5_main_publication_style",
  width = 12.5,
  height = 12.8
)

# ============================================================
# 12. assembled figures and final run summary
# ============================================================

figure5_core <- (latency_box | km_os_2group$plot) / (km_latency_2group$plot | ecdna_forest)
figure5_boxplots <- latency_box | os_box
figure5_km_curves <- km_os_2group$plot | km_latency_2group$plot
figure5_km_with_tables <- (km_os_2group$plot / km_os_2group$table) | (km_latency_2group$plot / km_latency_2group$table)
figure5_facets <- (
  (latency_facet_box + labs(tag = NULL)) /
    (os_facet_box + labs(tag = NULL))
) +
  plot_annotation(tag_levels = list(c("C", "D"))) &
  theme(
    plot.tag = element_text(face = "plain", size = 18),
    plot.tag.position = c(0.005, 1)
  )
figure5_three_group_km <- km_os_3group$plot | km_latency_3group$plot

make_square_km_with_risktable <- function(km_object, panel_title = NULL) {
  curve_plot <- km_object$plot +
    labs(title = panel_title) +
    theme(
      aspect.ratio = 1,
      legend.position = "top",
      legend.title = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.margin = margin(4, 8, 2, 8)
    )

  risk_table <- km_object$table +
    theme(
      plot.title = element_text(size = 12, hjust = 0),
      axis.title.x = element_text(size = 12, color = "black"),
      axis.text.x = element_text(size = 11, color = "black"),
      plot.margin = margin(0, 8, 4, 8)
    )

  ggpubr::ggarrange(
    curve_plot,
    risk_table,
    ncol = 1,
    heights = c(3.2, 1.15),
    align = "v"
  )
}

km_os_square_with_table <- make_square_km_with_risktable(
  km_object = km_os_2group,
  panel_title = "OS by ecDNA status"
)

km_latency_square_with_table <- make_square_km_with_risktable(
  km_object = km_latency_2group,
  panel_title = "Latency by ecDNA status"
)

figure5_km_with_risktables_square_vertical <- ggpubr::ggarrange(
  km_os_square_with_table,
  km_latency_square_with_table,
  ncol = 1,
  heights = c(1, 1),
  align = "v"
)

save_plot_pdf_png(km_os_square_with_table, "figure5_KM_OS_square_with_risktable", width = 5.2, height = 6.4)
save_plot_pdf_png(km_latency_square_with_table, "figure5_KM_latency_square_with_risktable", width = 5.2, height = 6.4)
save_plot_pdf_png(figure5_km_with_risktables_square_vertical, "figure5_KM_with_risktables_square_vertical", width = 5.4, height = 12.8)

save_plot_pdf_png(figure5_core, "figure5_core_panels", width = 12.5, height = 9.5)
save_plot_pdf_png(figure5_boxplots, "figure5_boxplots", width = 9, height = 4.8)
save_plot_pdf_png(figure5_km_curves, "figure5_KM_curves", width = 11, height = 6)
save_plot_pdf_png(figure5_km_with_tables, "figure5_KM_with_risktables", width = 12.5, height = 8.5)
save_plot_pdf_png(figure5_facets, "figure5_all_facets_dualStatusOnly", width = 18, height = 7.2)
save_plot_pdf_png(figure5_three_group_km, "figure5_KM_3group_curves", width = 13, height = 6)

capture.output(
  list(
    source_xlsx = PATH_CLINICAL_XLSX,
    source_sheet = SHEET_CLINICAL,
    included_disease_status = INCLUDED_DISEASE_STATUS,
    n_total_analysis = nrow(clinical_analysis),
    summary_by_ecDNA = summary_by_ecdna,
    summary_by_ecDNA_3group = summary_by_ecdna_3group,
    latency_group_table = table(latency_data$ecDNA_Status),
    os_group_table = table(os_data$ecDNA_Status),
    km_os_group_table = table(km_os_data$ecDNA_Status),
    age_missing = sum(is.na(clinical_analysis$Age)),
    sex_missing = sum(is.na(clinical_analysis$Sex)),
    cancer_level_summary = cancer_level_summary,
    within_cancer_os = within_cancer_os,
    within_cancer_latency = within_cancer_latency,
    dual_latency_cancer_types = sort(unique(as.character(latency_facet_data$CancerTypeFacet))),
    dual_os_cancer_types = sort(unique(as.character(os_facet_data$CancerTypeFacet)))
  ),
  file = file.path(DIR_TABLE, "figure5_run_summary.txt")
)