# ============================================================
# 0) Packages
# ============================================================
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

# ============================================================
# 1) Input configuration and output directories
# ============================================================
input_dir <- "Inputs"
output_dir <- "Outputs"

DIR_OUT   <- file.path(output_dir, "Figure5_survival")
DIR_RDS   <- file.path(DIR_OUT, "rds")
DIR_TABLE <- file.path(DIR_OUT, "tables")
DIR_PLOT  <- file.path(DIR_OUT, "plots")
DIR_SUPP  <- file.path(DIR_OUT, "supplementary")

purrr::walk(c(DIR_OUT, DIR_RDS, DIR_TABLE, DIR_PLOT, DIR_SUPP), ~{
  dir.create(.x, showWarnings = FALSE, recursive = TRUE)
})

PATH_CLINICAL_XLSX <- file.path(input_dir, "metadata", "sample_metadata.xlsx")
SHEET_CLINICAL     <- "Clinical survival information summary table"
SHEET_SAMPLE_INFO  <- "ecDNA WGS sample sequencing statistics-refined"

PATH_FIG2_ECDNA_SAMPLE <- file.path(output_dir, "Figure2_genomics", "rds", "ecDNA_sample_summary.rds")
PATH_FIG2_SAMPLE_INFO  <- file.path(output_dir, "Figure2_genomics", "rds", "sample_info_summary.rds")

# Optional public cohort data.
PATH_PUBLIC_SURVIVAL <- NA_character_

ECDNA2_COLS <- c(
  "ecDNA-" = "#E0E0E0",
  "ecDNA+" = "#8B0000"
)

ECDNA2_COLS_KM <- c(
  "ecDNA-" = "blue",
  "ecDNA+" = "#8B0000"
)

ECDNA3_COLS_KM <- c(
  "No amplicon" = "blue",
  "ecDNA- (amplicon)" = "#DE9B13",
  "ecDNA+" = "#8B0000"
)


KM_XLIM  <- c(0, 60)
KM_BREAK <- 12

MIN_N_PER_CANCER_FOR_KM <- 8
MIN_N_PER_GROUP_FOR_KM  <- 2
MIN_EVENTS_FOR_COX      <- 3
MIN_ECDNA_POS_FOR_COX   <- 2

theme_fig5 <- theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", size = 13),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )

# ============================================================
# 2) Helper functions
# ============================================================
truncate_km_data <- function(df, time_col, event_col = NULL, max_time = 60) {
  df2 <- df
  time_raw <- df2[[time_col]]
  over_max <- !is.na(time_raw) & time_raw > max_time
  df2[[time_col]] <- pmin(time_raw, max_time)
  if (!is.null(event_col) && event_col %in% names(df2)) {
    df2[[event_col]] <- ifelse(over_max, 0, df2[[event_col]])
  }
  df2
}

read_clinical_sheet <- function(path_xlsx, sheet) {
  temp_data <- read_excel(path = path_xlsx, sheet = sheet, n_max = 1)
  total_cols <- ncol(temp_data)
  
  col_types_vec <- rep("guess", total_cols)
  col_types_vec[6] <- "date"
  col_types_vec[7] <- "date"
  col_types_vec[8] <- "date"
  
  read_excel(
    path = path_xlsx,
    sheet = sheet,
    col_types = col_types_vec
  )
}

safe_as_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

safe_median <- function(x) {
  x <- x[is.finite(x) & !is.na(x)]
  if (length(x) == 0) return(NA_real_)
  median(x)
}

standardize_sample_id <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    str_remove("-.*$")
}

clean_sample_id <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    str_replace("-\\d+$", "")
}

plot_box_2group <- function(df, yvar, ylab, title = NULL, label_y = NULL) {
  p <- ggplot(df, aes(x = ecDNA_Status, y = .data[[yvar]], fill = ecDNA_Status)) +
    stat_boxplot(geom = "errorbar", width = 0.25, linewidth = 0.6, color = "black") +
    geom_boxplot(width = 0.5, linewidth = 0.5, outlier.shape = NA, alpha = 0.9, color = "black") +
    geom_jitter(width = 0.15, height = 0, size = 1.5, alpha = 0.9, color = "#2c3e50") +
    scale_fill_manual(values = ECDNA2_COLS, drop = FALSE) +
    labs(y = ylab, x = NULL, title = title) +
    theme_fig5 +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
      aspect.ratio = 1.1
    )
  
  if (!is.null(label_y)) {
    p <- p +
      stat_compare_means(
        method = "wilcox.test",
        method.args = list(alternative = "greater"),
        label.y = label_y,
        size = 4.5
      )
  }
  
  p
}


format_ecDNA_forest_table <- function(tbl) {
  tbl %>%
    mutate(
      model_label = as.character(model),
      hr_text = sprintf("%.2f (%.2f-%.2f)", HR, CI_low, CI_high),
      p_text = case_when(
        is.na(P_value) ~ "",
        P_value < 0.001 ~ "<0.001",
        TRUE ~ sprintf("%.3f", P_value)
      )
    ) %>%
    mutate(
      row_id = rev(seq_len(n()))
    )
}

make_ecDNA_forest_with_text <- function(tbl, title_letter = "d") {
  df_plot <- format_ecDNA_forest_table(tbl)
  
  p_table <- ggplot(df_plot) +
    geom_text(aes(x = 0.0, y = row_id, label = model_label), hjust = 0, size = 4) +
    geom_text(aes(x = 5.8, y = row_id, label = hr_text), hjust = 0, size = 4) +
    geom_text(aes(x = 9.0, y = row_id, label = p_text), hjust = 0, size = 4) +
    annotate("text", x = 0.0, y = max(df_plot$row_id) + 0.8,
             label = title_letter, fontface = "bold", hjust = 0, size = 6) +
    annotate("text", x = 0.4, y = max(df_plot$row_id) + 0.8,
             label = "Model", fontface = "bold", hjust = 0, size = 4.2) +
    annotate("text", x = 5.8, y = max(df_plot$row_id) + 0.8,
             label = "HR (95% CI)", fontface = "bold", hjust = 0, size = 4.2) +
    annotate("text", x = 9.0, y = max(df_plot$row_id) + 0.8,
             label = "P", fontface = "bold", hjust = 0, size = 4.2) +
    xlim(0, 10) +
    ylim(0.5, max(df_plot$row_id) + 1.2) +
    theme_void()
  
  p_forest <- ggplot(df_plot, aes(x = HR, y = row_id)) +
    geom_vline(xintercept = 1, linetype = 2, color = "grey50") +
    geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.15, color = "#8B0000") +
    geom_point(size = 3, color = "#8B0000") +
    scale_x_log10() +
    labs(
      x = "Hazard ratio for ecDNA+ (log scale)",
      y = NULL
    ) +
    ylim(0.5, max(df_plot$row_id) + 1.2) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank()
    )
  
  p_table + p_forest + patchwork::plot_layout(widths = c(2.3, 1.4))
}

auto_forest_xlim <- function(df, min_lower = 0.02, max_upper = 20, expand_mult = 1.15) {
  vals <- c(df$CI_low, df$CI_high)
  vals <- vals[is.finite(vals) & !is.na(vals) & vals > 0]
  
  raw_min <- min(vals)
  raw_max <- max(vals)
  
  lower <- 10^floor(log10(raw_min / expand_mult))
  upper <- 10^ceiling(log10(raw_max * expand_mult))
  
  lower <- max(min_lower, lower)
  upper <- min(max_upper, upper)
  
  c(lower, upper)
}

make_publication_forest <- function(forest_df,
                                    title_letter = "e",
                                    xbreaks = NULL,
                                    xlimits = NULL,
                                    footer_text = NULL) {
  
  if (is.null(xlimits)) {
    xlimits <- auto_forest_xlim(forest_df)
  }
  
  if (is.null(xbreaks)) {
    candidate_breaks <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100)
    xbreaks <- candidate_breaks[candidate_breaks >= xlimits[1] & candidate_breaks <= xlimits[2]]
  }
  
  bg_df <- forest_df %>%
    mutate(
      ymin = row_id - 0.5,
      ymax = row_id + 0.5,
      bg_fill = rep(c("#D9D9D9", "#BFBFBF"), length.out = n())
    )
  
  p_table <- ggplot() +
    geom_rect(
      data = bg_df,
      aes(xmin = -0.2, xmax = 11.6, ymin = ymin, ymax = ymax, fill = bg_fill),
      color = NA
    ) +
    scale_fill_identity() +
    geom_text(
      data = forest_df,
      aes(x = 0.00, y = row_id, label = section,
          fontface = ifelse(duplicated(section), "plain", "bold")),
      hjust = 0, size = 4
    ) +
    geom_text(
      data = forest_df,
      aes(x = 1.75, y = row_id, label = level),
      hjust = 0, size = 4
    ) +
    geom_text(
      data = forest_df,
      aes(x = 4.45, y = row_id,
          label = ifelse(is.na(n), "", paste0("(n = ", n, ")"))),
      hjust = 0, size = 4
    ) +
    geom_text(
      data = forest_df,
      aes(x = 6.25, y = row_id, label = hr_text),
      hjust = 0, size = 4
    ) +
    geom_text(
      data = forest_df,
      aes(x = 10.15, y = row_id, label = p_text),
      hjust = 0, size = 4
    ) +
    annotate("text", x = -0.1, y = max(forest_df$row_id) + 0.9,
             label = title_letter, fontface = "bold", hjust = 0, size = 6) +
    xlim(-0.2, 11.6) +
    ylim(0.5, max(forest_df$row_id) + 1.2) +
    theme_void()
  
  plot_df <- forest_df %>%
    filter(!is_reference, !is.na(HR)) %>%
    mutate(
      CI_low_plot = pmax(CI_low, xlimits[1]),
      CI_high_plot = pmin(CI_high, xlimits[2]),
      left_trunc = CI_low < xlimits[1],
      right_trunc = CI_high > xlimits[2],
      HR_plot = pmin(pmax(HR, xlimits[1]), xlimits[2])
    )
  
  p_forest <- ggplot() +
    geom_rect(
      data = bg_df,
      aes(xmin = xlimits[1], xmax = xlimits[2], ymin = ymin, ymax = ymax, fill = bg_fill),
      color = NA
    ) +
    scale_fill_identity() +
    geom_vline(xintercept = 1, linetype = 2, linewidth = 0.7) +
    geom_errorbarh(
      data = plot_df,
      aes(y = row_id, xmin = CI_low_plot, xmax = CI_high_plot),
      height = 0.15, linewidth = 1, color = "black"
    ) +
    geom_point(
      data = plot_df,
      aes(x = HR_plot, y = row_id),
      shape = 15, size = 3, color = "black"
    ) +
    geom_segment(
      data = plot_df %>% filter(left_trunc),
      aes(
        x = xlimits[1] * 1.02, xend = xlimits[1],
        y = row_id, yend = row_id
      ),
      arrow = arrow(length = unit(0.12, "cm"), ends = "first", type = "closed"),
      linewidth = 0.8, color = "black"
    ) +
    geom_segment(
      data = plot_df %>% filter(right_trunc),
      aes(
        x = xlimits[2] / 1.02, xend = xlimits[2],
        y = row_id, yend = row_id
      ),
      arrow = arrow(length = unit(0.12, "cm"), ends = "last", type = "closed"),
      linewidth = 0.8, color = "black"
    ) +
    scale_x_log10(
      breaks = xbreaks,
      labels = xbreaks
    ) +
    coord_cartesian(xlim = xlimits, clip = "off") +
    ylim(0.5, max(forest_df$row_id) + 1.2) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      plot.margin = margin(5.5, 15, 5.5, 0)
    )
  
  if (!is.null(footer_text)) {
    p_table <- p_table +
      annotate(
        "text",
        x = 0,
        y = 0.05,
        label = footer_text,
        hjust = 0,
        vjust = 0,
        size = 4,
        fontface = "italic"
      )
  }
  
  p_table + p_forest + patchwork::plot_layout(widths = c(2.6, 1.6))
}

safe_cor_test <- function(x, y, method = "spearman") {
  keep <- is.finite(x) & is.finite(y) & !is.na(x) & !is.na(y)
  x2 <- x[keep]
  y2 <- y[keep]
  if (length(x2) < 3) {
    return(tibble(
      method = method,
      estimate = NA_real_,
      p.value = NA_real_,
      n = length(x2)
    ))
  }
  ct <- suppressWarnings(cor.test(x2, y2, method = method, exact = FALSE))
  tibble(
    method = method,
    estimate = unname(ct$estimate),
    p.value = ct$p.value,
    n = length(x2)
  )
}

label_cor_text <- function(ctbl) {
  if (nrow(ctbl) == 0 || is.na(ctbl$estimate[1])) return("Spearman rho = NA, p = NA")
  paste0(
    "Spearman rho = ", sprintf("%.2f", ctbl$estimate[1]),
    ", p = ", ifelse(ctbl$p.value[1] < 0.001, "<0.001", sprintf("%.3f", ctbl$p.value[1])),
    ", n = ", ctbl$n[1]
  )
}

harmonize_external_survival_data <- function(df) {
  df %>%
    transmute(
      SampleID = as.character(SampleID),
      CancerType = str_trim(as.character(CancerType)),
      ecDNA_Status = factor(as.character(ecDNA_Status), levels = c("ecDNA-", "ecDNA+")),
      OS_Months = suppressWarnings(as.numeric(OS_Months)),
      OS_event = suppressWarnings(as.numeric(OS_event)),
      Latency_Months = suppressWarnings(as.numeric(Latency_Months)),
      Age = suppressWarnings(as.numeric(Age)),
      Sex = factor(tolower(as.character(Sex)), levels = c("male", "female")),
      Cohort = "Public"
    ) %>%
    mutate(
      Age_z = ifelse(sum(!is.na(Age)) >= 2, as.numeric(scale(Age)), NA_real_)
    )
}

run_ecDNA_cox_within_cancer <- function(df, time_col, event_col, cancer_col = "CancerType",
                                        min_n = 8, min_events = 3, min_pos = 2) {
  split_df <- split(df, df[[cancer_col]])
  
  purrr::map_dfr(names(split_df), function(ct) {
    d <- split_df[[ct]] %>%
      filter(
        !is.na(.data[[time_col]]),
        !is.na(.data[[event_col]]),
        !is.na(ecDNA_Status)
      ) %>%
      mutate(ecDNA_Status = factor(ecDNA_Status, levels = c("ecDNA-", "ecDNA+")))
    
    n_total <- nrow(d)
    n_pos <- sum(d$ecDNA_Status == "ecDNA+", na.rm = TRUE)
    n_neg <- sum(d$ecDNA_Status == "ecDNA-", na.rm = TRUE)
    n_event <- sum(d[[event_col]] == 1, na.rm = TRUE)
    
    if (n_total < min_n || n_pos < min_pos || n_neg < min_pos || n_event < min_events) {
      return(tibble(
        CancerType = ct,
        n = n_total,
        n_ecDNA_pos = n_pos,
        n_ecDNA_neg = n_neg,
        n_event = n_event,
        HR = NA_real_,
        CI_low = NA_real_,
        CI_high = NA_real_,
        P_value = NA_real_,
        note = "insufficient data"
      ))
    }
    
    fit <- tryCatch(
      coxph(as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ ecDNA_Status")), data = d),
      error = function(e) NULL
    )
    
    if (is.null(fit)) {
      return(tibble(
        CancerType = ct,
        n = n_total,
        n_ecDNA_pos = n_pos,
        n_ecDNA_neg = n_neg,
        n_event = n_event,
        HR = NA_real_,
        CI_low = NA_real_,
        CI_high = NA_real_,
        P_value = NA_real_,
        note = "model failed"
      ))
    }
    
    tt <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(str_detect(term, "ecDNA_Status"))
    
    if (nrow(tt) == 0) {
      return(tibble(
        CancerType = ct,
        n = n_total,
        n_ecDNA_pos = n_pos,
        n_ecDNA_neg = n_neg,
        n_event = n_event,
        HR = NA_real_,
        CI_low = NA_real_,
        CI_high = NA_real_,
        P_value = NA_real_,
        note = "no ecDNA term"
      ))
    }
    
    tibble(
      CancerType = ct,
      n = n_total,
      n_ecDNA_pos = n_pos,
      n_ecDNA_neg = n_neg,
      n_event = n_event,
      HR = tt$estimate[1],
      CI_low = tt$conf.low[1],
      CI_high = tt$conf.high[1],
      P_value = tt$p.value[1],
      note = "ok"
    )
  })
}

make_surv_facet_data <- function(df, time_col, event_col, max_time = 60) {
  df %>%
    filter(
      !is.na(.data[[time_col]]),
      !is.na(.data[[event_col]]),
      !is.na(ecDNA_Status),
      !is.na(CancerType)
    ) %>%
    mutate(
      time_plot = pmin(.data[[time_col]], max_time),
      event_plot = ifelse(.data[[time_col]] > max_time, 0, .data[[event_col]]),
      CancerType = str_trim(as.character(CancerType)),
      ecDNA_Status = factor(ecDNA_Status, levels = c("ecDNA-", "ecDNA+"))
    ) %>%
    group_by(CancerType) %>%
    filter(
      n() >= MIN_N_PER_CANCER_FOR_KM,
      sum(ecDNA_Status == "ecDNA+", na.rm = TRUE) >= MIN_N_PER_GROUP_FOR_KM,
      sum(ecDNA_Status == "ecDNA-", na.rm = TRUE) >= MIN_N_PER_GROUP_FOR_KM
    ) %>%
    ungroup()
}

plot_facet_km_manual <- function(df, time_col, event_col, xlab, file_stub) {
  cancer_types_use <- unique(df$CancerType)
  
  plist <- purrr::map(cancer_types_use, function(ct) {
    d <- df %>% filter(CancerType == ct)
    
    surv_formula <- as.formula(
      paste0("Surv(", time_col, ", ", event_col, ") ~ ecDNA_Status")
    )
    
    out <- tryCatch({
      fit <- survfit(surv_formula, data = d)
      
      g <- ggsurvplot(
        fit,
        data = d,
        risk.table = FALSE,
        pval = TRUE,
        conf.int = FALSE,
        palette = unname(ECDNA2_COLS_KM),
        xlab = xlab,
        ylab = "Survival probability",
        title = ct,
        legend.title = NULL,
        legend.labs = levels(d$ecDNA_Status),
        xlim = KM_XLIM,
        break.time.by = KM_BREAK,
        ggtheme = theme_classic(base_size = 10)
      )
      
      g$plot +
        theme(
          plot.title = element_text(face = "bold", size = 10),
          legend.position = "none"
        )
    }, error = function(e) {
      ggplot() +
        annotate("text", x = 1, y = 1, label = paste0(ct, "\nKM failed")) +
        theme_void()
    })
    
    out
  })
  
  ncol_use <- min(4, length(plist))
  p <- wrap_plots(plist, ncol = ncol_use) +
    plot_annotation(title = file_stub)
  
  ggsave(
    file.path(DIR_PLOT, paste0(file_stub, ".pdf")),
    p,
    width = 4 * ncol_use,
    height = 3.6 * ceiling(length(plist) / ncol_use)
  )
  ggsave(
    file.path(DIR_PLOT, paste0(file_stub, ".png")),
    p,
    width = 4 * ncol_use,
    height = 3.6 * ceiling(length(plist) / ncol_use),
    dpi = 300
  )
  
  p
}

keep_dual_ecdna_cancer_types <- function(df, cancer_col = "CancerType", status_col = "ecDNA_Status") {
  keep_types <- df %>%
    filter(!is.na(.data[[cancer_col]]), !is.na(.data[[status_col]])) %>%
    distinct(.data[[cancer_col]], .data[[status_col]]) %>%
    count(.data[[cancer_col]], name = "n_status") %>%
    filter(n_status >= 2) %>%
    pull(.data[[cancer_col]])
  
  df %>%
    filter(.data[[cancer_col]] %in% keep_types)
}

# ============================================================
# 3) Read sample information with age and sex
# ============================================================
sampleinfo_detail <- read_excel(PATH_CLINICAL_XLSX, sheet = SHEET_SAMPLE_INFO) %>%
  filter(used == "yes") %>%
  transmute(
    Sample      = as.character(Sample),
    SampleID    = clean_sample_id(Sample),
    Cancer_type = str_trim(as.character(Cancer_type)),
    Pri_Met     = as.character(Pri_Met),
    Sex         = factor(tolower(as.character(Sex)), levels = c("male", "female")),
    Age         = suppressWarnings(as.numeric(Age))
  ) %>%
  distinct(SampleID, .keep_all = TRUE)

sampleInfo <- readRDS(PATH_FIG2_SAMPLE_INFO)
sampleInfo$SampleID <- sub("-.*", "", sampleInfo$Sample)
Met_sample <- sampleInfo$SampleID[sampleInfo$Pri_Met == "Met"]

# ============================================================
# 4) Read clinical data
# ============================================================
tbl_clinical_raw <- read_clinical_sheet(PATH_CLINICAL_XLSX, SHEET_CLINICAL)
tbl_clinical_raw <- tbl_clinical_raw[tbl_clinical_raw$SampleID %in% Met_sample, ]

tbl_fig2_ecsample <- if (file.exists(PATH_FIG2_ECDNA_SAMPLE)) {
  readRDS(PATH_FIG2_ECDNA_SAMPLE) %>%
    select(Sample, ecDNA_status) %>%
    distinct() %>%
    mutate(Sample = standardize_sample_id(Sample))
} else {
  tibble(Sample = character(), ecDNA_status = character())
}

# ============================================================
# 5) Clean clinical table and recode ecDNA
# ============================================================
tbl_clinical <- tbl_clinical_raw %>%
  rename(
    SampleID         = SampleID,
    CancerType       = CancerType,
    Pair             = Pair,
    T0_Date          = T0_Date,
    T1_Date          = T1_Date,
    T2_Date          = T2_Date,
    Latency_Months   = Latency_Months,
    OS_Months        = OS_Months,
    Status           = Status,
    ecDNA_status_raw = `ecDNA.status`,
    LastFollowUp     = `最近随访日期`,
    Notes            = Notes
  ) %>%
  mutate(
    SampleID         = standardize_sample_id(SampleID),
    CancerType       = str_trim(as.character(CancerType)),
    Pair             = str_trim(as.character(Pair)),
    ecDNA_status_raw = str_trim(as.character(ecDNA_status_raw)),
    Latency_Months   = safe_as_numeric(Latency_Months),
    OS_Months        = safe_as_numeric(OS_Months),
    Status           = safe_as_numeric(Status)
  ) %>%
  left_join(
    tbl_fig2_ecsample %>% rename(ecDNA_status_fig2 = ecDNA_status),
    by = c("SampleID" = "Sample")
  ) %>%
  mutate(
    ecDNA_Status_3group = case_when(
      ecDNA_status_raw == "ecDNA+" ~ "ecDNA+",
      ecDNA_status_raw == "Amplicon (Linear)" ~ "ecDNA- (amplicon)",
      ecDNA_status_raw == "Negative" ~ "No amplicon",
      ecDNA_status_fig2 == "ecDNA+" ~ "ecDNA+",
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
    OS_event = case_when(
      Status == 1 ~ 1,
      Status == 0 ~ 0,
      TRUE ~ NA_real_
    )
  )
# ============================================================
# 6) Merge age and sex and build analysis table
# ============================================================
tbl_analysis <- tbl_clinical %>%
  filter(!is.na(ecDNA_Status)) %>%
  filter(!(str_starts(SampleID, "A") & is.na(OS_Months) & is.na(Latency_Months))) %>%
  filter(!(is.na(OS_Months) & is.na(Latency_Months) & is.na(T1_Date) & is.na(T2_Date))) %>%
  left_join(
    sampleinfo_detail %>% select(SampleID, Age, Sex, Cancer_type),
    by = "SampleID"
  ) %>%
  mutate(
    CancerType = coalesce(CancerType, Cancer_type),
    Age = suppressWarnings(as.numeric(Age)),
    Sex = factor(Sex, levels = c("male", "female"))
  ) %>%
  mutate(
    info_score =
      1 * !is.na(OS_Months) +
      1 * !is.na(Status) +
      1 * !is.na(Latency_Months) +
      1 * !is.na(T2_Date) +
      1 * !is.na(Age) +
      1 * !is.na(Sex)
  ) %>%
  arrange(SampleID, desc(info_score), desc(T2_Date)) %>%
  distinct(SampleID, .keep_all = TRUE) %>%
  select(-info_score, -Cancer_type) %>%
  mutate(
    Age_z = ifelse(sum(!is.na(Age)) >= 2, as.numeric(scale(Age)), NA_real_),
    Cohort = "Internal"
  )
saveRDS(tbl_analysis, file.path(DIR_RDS, "figure5_analysis_table_2group.rds"))
write.csv(tbl_analysis, file.path(DIR_TABLE, "figure5_analysis_table_2group.csv"), row.names = FALSE)

# ============================================================
# 6.1) Optional public data integration
# ============================================================
tbl_public_analysis <- tibble()

if (!is.na(PATH_PUBLIC_SURVIVAL) && file.exists(PATH_PUBLIC_SURVIVAL)) {
  if (grepl("\\.rds$", PATH_PUBLIC_SURVIVAL, ignore.case = TRUE)) {
    tbl_public_raw <- readRDS(PATH_PUBLIC_SURVIVAL)
  } else if (grepl("\\.csv$", PATH_PUBLIC_SURVIVAL, ignore.case = TRUE)) {
    tbl_public_raw <- read.csv(PATH_PUBLIC_SURVIVAL, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    stop("Unsupported public survival file format. Please use .csv or .rds")
  }
  
  tbl_public_analysis <- harmonize_external_survival_data(tbl_public_raw)
}

tbl_analysis_combined <- bind_rows(
  tbl_analysis %>% mutate(Cohort = "Internal"),
  tbl_public_analysis
) %>%
  mutate(
    ecDNA_Status = factor(ecDNA_Status, levels = c("ecDNA-", "ecDNA+")),
    Sex = factor(Sex, levels = c("male", "female"))
  )

saveRDS(tbl_analysis_combined, file.path(DIR_RDS, "figure5_analysis_table_combined.rds"))
write.csv(tbl_analysis_combined, file.path(DIR_TABLE, "figure5_analysis_table_combined.csv"), row.names = FALSE)

# ============================================================
# 7) Descriptive summary
# ============================================================
tbl_summary <- tbl_analysis_combined %>%
  group_by(ecDNA_Status) %>%
  summarise(
    n_total = n(),
    n_latency = sum(!is.na(Latency_Months)),
    median_latency = median(Latency_Months, na.rm = TRUE),
    n_os = sum(!is.na(OS_Months)),
    median_os = median(OS_Months, na.rm = TRUE),
    n_event = sum(OS_event == 1, na.rm = TRUE),
    n_age = sum(!is.na(Age)),
    n_sex = sum(!is.na(Sex)),
    .groups = "drop"
  )

write.csv(tbl_summary, file.path(DIR_TABLE, "figure5_group_summary_2group.csv"), row.names = FALSE)

# ============================================================
# 8) Boxplots
# ============================================================
tbl_latency <- tbl_analysis_combined %>%
  filter(!is.na(Latency_Months), !is.na(ecDNA_Status))

plot_latency <- plot_box_2group(
  df = tbl_latency,
  yvar = "Latency_Months",
  ylab = "Time to brain metastasis (Months)",
  title = "a",
  label_y = max(tbl_latency$Latency_Months, na.rm = TRUE) * 1.05
)

plot_latency_facet <- ggplot(
  tbl_latency,
  aes(x = ecDNA_Status, y = Latency_Months, fill = ecDNA_Status)
) +
  stat_boxplot(geom = "errorbar", width = 0.25, linewidth = 0.5, color = "black") +
  geom_boxplot(width = 0.5, linewidth = 0.5, outlier.shape = NA, alpha = 0.9, color = "black") +
  geom_jitter(width = 0.15, height = 0, size = 1, alpha = 0.85, color = "#2c3e50") +
  scale_fill_manual(values = ECDNA2_COLS, drop = FALSE) +
  labs(y = "Time to brain metastasis (Months)", x = NULL) +
  facet_wrap(~ CancerType, ncol = 5, scales = "fixed", strip.position = "top") +
  theme_classic(base_size = 10) +
  theme(
    axis.text = element_text(color = "black", size = 10),
    axis.text.x = element_text(color = "black", size = 10, angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(color = "black", size = 12),
    axis.line = element_line(linewidth = 0.6, lineend = "square"),
    axis.ticks = element_line(linewidth = 0.6, color = "black"),
    legend.position = "none",
    strip.background = element_rect(color = "black", fill = "white", linewidth = 0.6),
    strip.text = element_text(color = "black", size = 10),
    panel.spacing = unit(1, "lines")
  ) +
  stat_compare_means(
    method = "wilcox.test",
    method.args = list(alternative = "greater"),
    label.x = 1.2,
    label.y.npc = 0.9,
    size = 3.2
  )

tbl_os <- tbl_analysis_combined %>%
  filter(!is.na(OS_Months), !is.na(ecDNA_Status))

plot_os <- plot_box_2group(
  df = tbl_os,
  yvar = "OS_Months",
  ylab = "OS (Months)",
  title = "Supplementary",
  label_y = max(tbl_os$OS_Months, na.rm = TRUE) * 1.05
)

plot_os_facet <- ggplot(
  tbl_os,
  aes(x = ecDNA_Status, y = OS_Months, fill = ecDNA_Status)
) +
  stat_boxplot(geom = "errorbar", width = 0.25, linewidth = 0.5, color = "black") +
  geom_boxplot(width = 0.5, linewidth = 0.5, outlier.shape = NA, alpha = 0.9, color = "black") +
  geom_jitter(width = 0.15, height = 0, size = 2, alpha = 0.85, color = "#2c3e50") +
  scale_fill_manual(values = ECDNA2_COLS, drop = FALSE) +
  labs(y = "OS (Months)", x = NULL) +
  facet_wrap(~ CancerType, ncol = 5, scales = "fixed", strip.position = "top") +
  theme_classic(base_size = 10) +
  theme(
    axis.text = element_text(color = "black", size = 10),
    axis.text.x = element_text(color = "black", size = 10, angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(color = "black", size = 12),
    axis.line = element_line(linewidth = 0.6, lineend = "square"),
    axis.ticks = element_line(linewidth = 0.6, color = "black"),
    legend.position = "none",
    strip.background = element_rect(color = "black", fill = "white", linewidth = 0.6),
    strip.text = element_text(color = "black", size = 10),
    panel.spacing = unit(1, "lines")
  ) +
  stat_compare_means(
    method = "wilcox.test",
    method.args = list(alternative = "greater"),
    label.x = 1.2,
    label.y.npc = 0.9,
    size = 3.2
  )

# dual-status-only facet plots
tbl_latency_facet_dual <- tbl_latency %>%
  keep_dual_ecdna_cancer_types(cancer_col = "CancerType", status_col = "ecDNA_Status")

tbl_latency_facet_dual$CancerType <- factor(
  tbl_latency_facet_dual$CancerType,
  levels = sort(unique(tbl_latency_facet_dual$CancerType))
)

plot_latency_facet_dual <- ggplot(
  tbl_latency_facet_dual,
  aes(x = ecDNA_Status, y = Latency_Months, fill = ecDNA_Status)
) +
  stat_boxplot(geom = "errorbar", width = 0.25, linewidth = 0.5, color = "black") +
  geom_boxplot(width = 0.5, linewidth = 0.5, outlier.shape = NA, alpha = 0.9, color = "black") +
  geom_jitter(width = 0.15, height = 0, size = 1, alpha = 0.85, color = "#2c3e50") +
  scale_fill_manual(values = ECDNA2_COLS, drop = FALSE) +
  labs(
    y = "Time to brain metastasis (Months)",
    x = NULL,
    title = "Supplementary"
  ) +
  facet_wrap(~ CancerType, ncol = 5, scales = "fixed", strip.position = "top") +
  theme_classic(base_size = 10) +
  theme(
    axis.text = element_text(color = "black", size = 10),
    axis.text.x = element_text(color = "black", size = 10, angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(color = "black", size = 12),
    axis.line = element_line(linewidth = 0.6, lineend = "square"),
    axis.ticks = element_line(linewidth = 0.6, color = "black"),
    legend.position = "none",
    strip.background = element_rect(color = "black", fill = "white", linewidth = 0.6),
    strip.text = element_text(color = "black", size = 10),
    panel.spacing = unit(1, "lines")
  ) +
  stat_compare_means(
    method = "wilcox.test",
    method.args = list(alternative = "greater"),
    label.x = 1.2,
    label.y.npc = 0.9,
    size = 3.2
  )

tbl_os_facet_dual <- tbl_os %>%
  keep_dual_ecdna_cancer_types(cancer_col = "CancerType", status_col = "ecDNA_Status")

tbl_os_facet_dual$CancerType <- factor(
  tbl_os_facet_dual$CancerType,
  levels = sort(unique(tbl_os_facet_dual$CancerType))
)

plot_os_facet_dual <- ggplot(
  tbl_os_facet_dual,
  aes(x = ecDNA_Status, y = OS_Months, fill = ecDNA_Status)
) +
  stat_boxplot(geom = "errorbar", width = 0.25, linewidth = 0.5, color = "black") +
  geom_boxplot(width = 0.5, linewidth = 0.5, outlier.shape = NA, alpha = 0.9, color = "black") +
  geom_jitter(width = 0.15, height = 0, size = 1, alpha = 0.85, color = "#2c3e50") +
  scale_fill_manual(values = ECDNA2_COLS, drop = FALSE) +
  labs(
    y = "OS (Months)",
    x = NULL,
    title = "Supplementary"
  ) +
  facet_wrap(~ CancerType, ncol = 5, scales = "fixed", strip.position = "top") +
  theme_classic(base_size = 10) +
  theme(
    axis.text = element_text(color = "black", size = 10),
    axis.text.x = element_text(color = "black", size = 10, angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(color = "black", size = 12),
    axis.line = element_line(linewidth = 0.6, lineend = "square"),
    axis.ticks = element_line(linewidth = 0.6, color = "black"),
    legend.position = "none",
    strip.background = element_rect(color = "black", fill = "white", linewidth = 0.6),
    strip.text = element_text(color = "black", size = 10),
    panel.spacing = unit(1, "lines")
  ) +
  stat_compare_means(
    method = "wilcox.test",
    method.args = list(alternative = "greater"),
    label.x = 1.2,
    label.y.npc = 0.9,
    size = 3.2
  )

# ============================================================
# 9) Survival curves
# ============================================================
tbl_km_os <- tbl_analysis_combined %>%
  filter(!is.na(OS_Months), !is.na(OS_event), !is.na(ecDNA_Status))

tbl_km_os_plot <- truncate_km_data(
  df = tbl_km_os,
  time_col = "OS_Months",
  event_col = "OS_event",
  max_time = KM_XLIM[2]
)

fit_km_os <- survfit(
  Surv(OS_Months, OS_event) ~ ecDNA_Status,
  data = tbl_km_os_plot
)

km_os <- ggsurvplot(
  fit_km_os,
  data = tbl_km_os_plot,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  palette = unname(ECDNA2_COLS_KM),
  xlab = "OS (Months)",
  ylab = "Survival probability",
  title = "b",
  legend.title = NULL,
  legend.labs = levels(tbl_km_os_plot$ecDNA_Status),
  risk.table.height = 0.23,
  xlim = KM_XLIM,
  break.time.by = KM_BREAK,
  ggtheme = theme_classic(base_size = 12)
)

logrank_os <- survdiff(Surv(OS_Months, OS_event) ~ ecDNA_Status, data = tbl_km_os)
sink(file.path(DIR_TABLE, "figure5_logrank_OS_2group.txt"))
print(logrank_os)
sink()

tbl_km_latency <- tbl_analysis_combined %>%
  filter(!is.na(Latency_Months), !is.na(ecDNA_Status)) %>%
  mutate(Latency_event = 1)

tbl_km_latency_plot <- truncate_km_data(
  df = tbl_km_latency,
  time_col = "Latency_Months",
  event_col = "Latency_event",
  max_time = KM_XLIM[2]
)

fit_km_latency <- survfit(
  Surv(Latency_Months, Latency_event) ~ ecDNA_Status,
  data = tbl_km_latency_plot
)

km_latency <- ggsurvplot(
  fit_km_latency,
  data = tbl_km_latency_plot,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  palette = unname(ECDNA2_COLS_KM),
  xlab = "Time to brain metastasis (Months)",
  ylab = "Event-free proportion",
  title = "Supplementary",
  legend.title = NULL,
  legend.labs = levels(tbl_km_latency_plot$ecDNA_Status),
  xlim = KM_XLIM,
  break.time.by = KM_BREAK,
  ggtheme = theme_classic(base_size = 12)
)

# ============================================================
# 9.1) Three-group KM curves
# Goal:
#   3-category KM for:
#   1) OS
#   2) Time to brain metastasis (Latency)
# Groups:
#   No amplicon / ecDNA- (amplicon) / ecDNA+
# ============================================================

plot_km_3group <- function(df, time_col, event_col, group_col,
                           xlab, ylab, title,
                           legend_labs = c("No amplicon", "ecDNA- (amplicon)", "ecDNA+")) {
  
  df_use <- df %>%
    filter(
      !is.na(.data[[time_col]]),
      !is.na(.data[[event_col]]),
      !is.na(.data[[group_col]])
    ) %>%
    mutate(
      .time  = .data[[time_col]],
      .event = .data[[event_col]],
      .group = factor(.data[[group_col]],
                      levels = c("No amplicon", "ecDNA- (amplicon)", "ecDNA+"))
    )
  
  fit <- survfit(
    survival::Surv(.time, .event) ~ .group,
    data = df_use
  )
  
  ggsurvplot(
    fit,
    data = df_use,
    risk.table = TRUE,
    pval = TRUE,
    conf.int = FALSE,
    palette = unname(ECDNA3_COLS_KM),
    xlab = xlab,
    ylab = ylab,
    title = title,
    legend.title = NULL,
    legend.labs = legend_labs,
    risk.table.height = 0.23,
    xlim = KM_XLIM,
    break.time.by = KM_BREAK,
    ggtheme = theme_classic(base_size = 12)
  )
}

# ----------------------------
# 9.1.1) Three-group OS KM
# ----------------------------
tbl_km_os_3group <- tbl_analysis_combined %>%
  filter(!is.na(OS_Months), !is.na(OS_event), !is.na(ecDNA_Status_3group))

tbl_km_os_3group_plot <- truncate_km_data(
  df = tbl_km_os_3group,
  time_col = "OS_Months",
  event_col = "OS_event",
  max_time = KM_XLIM[2]
)

km_os_3group <- plot_km_3group(
  df = tbl_km_os_3group_plot,
  time_col = "OS_Months",
  event_col = "OS_event",
  group_col = "ecDNA_Status_3group",
  xlab = "OS (Months)",
  ylab = "Survival probability",
  title = "Supplementary"
)

logrank_os_3group <- survdiff(
  Surv(OS_Months, OS_event) ~ ecDNA_Status_3group,
  data = tbl_km_os_3group
)

sink(file.path(DIR_TABLE, "figure5_logrank_OS_3group.txt"))
print(logrank_os_3group)
sink()

# ----------------------------
# 9.1.2) Three-group latency KM
# ----------------------------
tbl_km_latency_3group <- tbl_analysis_combined %>%
  filter(!is.na(Latency_Months), !is.na(ecDNA_Status_3group)) %>%
  mutate(Latency_event = 1)

tbl_km_latency_3group_plot <- truncate_km_data(
  df = tbl_km_latency_3group,
  time_col = "Latency_Months",
  event_col = "Latency_event",
  max_time = KM_XLIM[2]
)

km_latency_3group <- plot_km_3group(
  df = tbl_km_latency_3group_plot,
  time_col = "Latency_Months",
  event_col = "Latency_event",
  group_col = "ecDNA_Status_3group",
  xlab = "Time to brain metastasis (Months)",
  ylab = "Event-free proportion",
  title = "Supplementary"
)

logrank_latency_3group <- survdiff(
  Surv(Latency_Months, Latency_event) ~ ecDNA_Status_3group,
  data = tbl_km_latency_3group
)

sink(file.path(DIR_TABLE, "figure5_logrank_Latency_3group.txt"))
print(logrank_latency_3group)
sink()


# ============================================================
km_os_3group$plot      <- km_os_3group$plot + coord_cartesian(xlim = KM_XLIM, clip = "on")
km_latency_3group$plot <- km_latency_3group$plot + coord_cartesian(xlim = KM_XLIM, clip = "on")

ggsave(file.path(DIR_PLOT, "figure5_KM_OS_3group_curve.pdf"), km_os_3group$plot, width = 6.5, height = 5.2)
ggsave(file.path(DIR_PLOT, "figure5_KM_OS_3group_curve.png"), km_os_3group$plot, width = 6.5, height = 5.2, dpi = 300)
ggsave(file.path(DIR_PLOT, "figure5_KM_OS_3group_risktable.pdf"), km_os_3group$table, width = 6.5, height = 2.4)

ggsave(file.path(DIR_PLOT, "figure5_KM_latency_3group_curve.pdf"), km_latency_3group$plot, width = 6.5, height = 5.2)
ggsave(file.path(DIR_PLOT, "figure5_KM_latency_3group_curve.png"), km_latency_3group$plot, width = 6.5, height = 5.2, dpi = 300)
ggsave(file.path(DIR_PLOT, "figure5_KM_latency_3group_risktable.pdf"), km_latency_3group$table, width = 6.5, height = 2.4)




# ============================================================
# 10) Main inferential analyses
# Goal:
#   1) ecDNA+ vs OS
#   2) ecDNA+ vs latency
#   3) ecDNA_Status vs CancerType
# ============================================================

# ----------------------------
# 10.0) Helper
# ----------------------------
tidy_cox_model2 <- function(model, model_name) {
  broom::tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(model = model_name) %>%
    rename(
      HR = estimate,
      CI_low = conf.low,
      CI_high = conf.high,
      P_value = p.value
    )
}


# ----------------------------
# 10.1) Helper: build full forest rows from a Cox model
# ----------------------------
build_full_forest_from_model <- function(model_tbl, df_model, ref_cancer_level) {
  
  # ecDNA
  tab_ecdna <- df_model %>%
    count(ecDNA_Status, name = "n")
  
  ref_ecdna <- tibble(
    section = "ecDNA status",
    level = "ecDNA-",
    n = tab_ecdna$n[match("ecDNA-", tab_ecdna$ecDNA_Status)],
    HR = NA_real_,
    CI_low = NA_real_,
    CI_high = NA_real_,
    P_value = NA_real_,
    is_reference = TRUE
  )
  
  est_ecdna <- tibble(
    section = "ecDNA status",
    level = "ecDNA+",
    n = tab_ecdna$n[match("ecDNA+", tab_ecdna$ecDNA_Status)],
    HR = model_tbl$HR[match("ecDNA_StatusecDNA+", model_tbl$term)],
    CI_low = model_tbl$CI_low[match("ecDNA_StatusecDNA+", model_tbl$term)],
    CI_high = model_tbl$CI_high[match("ecDNA_StatusecDNA+", model_tbl$term)],
    P_value = model_tbl$P_value[match("ecDNA_StatusecDNA+", model_tbl$term)],
    is_reference = FALSE
  )
  
  # Sex
  tab_sex <- df_model %>%
    count(Sex, name = "n")
  
  ref_sex <- tibble(
    section = "Sex",
    level = "male",
    n = tab_sex$n[match("male", tab_sex$Sex)],
    HR = NA_real_,
    CI_low = NA_real_,
    CI_high = NA_real_,
    P_value = NA_real_,
    is_reference = TRUE
  )
  
  est_sex <- tibble(
    section = "Sex",
    level = "female",
    n = tab_sex$n[match("female", tab_sex$Sex)],
    HR = model_tbl$HR[match("Sexfemale", model_tbl$term)],
    CI_low = model_tbl$CI_low[match("Sexfemale", model_tbl$term)],
    CI_high = model_tbl$CI_high[match("Sexfemale", model_tbl$term)],
    P_value = model_tbl$P_value[match("Sexfemale", model_tbl$term)],
    is_reference = FALSE
  )
  
  # Cancer type
  tab_cancer <- df_model %>%
    count(CancerType, name = "n")
  
  ref_cancer <- tibble(
    section = "Cancer type",
    level = ref_cancer_level,
    n = tab_cancer$n[match(ref_cancer_level, tab_cancer$CancerType)],
    HR = NA_real_,
    CI_low = NA_real_,
    CI_high = NA_real_,
    P_value = NA_real_,
    is_reference = TRUE
  )
  
  est_cancer <- model_tbl %>%
    filter(str_detect(term, "^CancerType")) %>%
    mutate(
      level = str_remove(term, "^CancerType"),
      section = "Cancer type",
      n = tab_cancer$n[match(level, tab_cancer$CancerType)],
      is_reference = FALSE
    ) %>%
    select(section, level, n, HR, CI_low, CI_high, P_value, is_reference)
  
  bind_rows(
    ref_ecdna, est_ecdna,
    ref_sex, est_sex,
    ref_cancer, est_cancer
  ) %>%
    mutate(
      hr_text = case_when(
        is_reference ~ "Reference",
        is.na(HR) ~ "",
        TRUE ~ sprintf("%.2f (%.2f-%.2f)", HR, CI_low, CI_high)
      ),
      p_text = case_when(
        is_reference ~ "",
        is.na(P_value) ~ "",
        P_value < 0.001 ~ "<0.001",
        TRUE ~ sprintf("%.3f", P_value)
      )
    ) %>%
    mutate(
      row_id = rev(seq_len(n()))
    )
}

# ============================================================
# 10.2) Pooled OS analysis
# ============================================================
tbl_os_analysis <- tbl_analysis_combined %>%
  filter(
    !is.na(OS_Months),
    !is.na(OS_event),
    !is.na(ecDNA_Status),
    !is.na(CancerType)
  ) %>%
  mutate(
    ecDNA_Status = factor(ecDNA_Status, levels = c("ecDNA-", "ecDNA+")),
    CancerType = factor(str_trim(as.character(CancerType))),
    Sex = factor(Sex, levels = c("male", "female"))
  )

write.csv(
  tbl_os_analysis,
  file.path(DIR_TABLE, "figure5_OS_analysis_input.csv"),
  row.names = FALSE
)

tbl_os_analysis_clinical <- tbl_os_analysis %>%
  filter(!is.na(Sex))

# ----------------------------
# 10.2.1) OS models
# ----------------------------
cox_os_uni_ecDNA <- coxph(
  Surv(OS_Months, OS_event) ~ ecDNA_Status,
  data = tbl_os_analysis
)

cox_os_uni_cancertype <- coxph(
  Surv(OS_Months, OS_event) ~ CancerType,
  data = tbl_os_analysis
)

cox_os_multi_cancertype <- coxph(
  Surv(OS_Months, OS_event) ~ ecDNA_Status + CancerType,
  data = tbl_os_analysis
)


cox_os_multi_cancertype_clinical <- coxph(
  Surv(OS_Months, OS_event) ~ ecDNA_Status + CancerType + Sex,
  data = tbl_os_analysis_clinical
)

cox_zph_os_uni_ecDNA <- cox.zph(cox_os_uni_ecDNA)
cox_zph_os_multi_cancertype <- cox.zph(cox_os_multi_cancertype)
cox_zph_os_multi_cancertype_clinical <- cox.zph(cox_os_multi_cancertype_clinical)

# ============================================================
# 10.3) Pooled latency analysis
# Assumption:
#   All included latency cases are treated as observed events
#   in the latency time-to-event analysis.
# ============================================================
tbl_latency_analysis <- tbl_analysis_combined %>%
  filter(
    !is.na(Latency_Months),
    !is.na(ecDNA_Status),
    !is.na(CancerType)
  ) %>%
  mutate(
    Latency_event = 1,
    ecDNA_Status = factor(ecDNA_Status, levels = c("ecDNA-", "ecDNA+")),
    CancerType = factor(str_trim(as.character(CancerType))),
    Sex = factor(Sex, levels = c("male", "female"))
  )

write.csv(
  tbl_latency_analysis,
  file.path(DIR_TABLE, "figure5_Latency_analysis_input.csv"),
  row.names = FALSE
)

tbl_latency_analysis_clinical <- tbl_latency_analysis %>%
  filter(!is.na(Sex))

# ----------------------------
# 10.3.1) Latency models
# ----------------------------
cox_latency_uni_ecDNA <- coxph(
  Surv(Latency_Months, Latency_event) ~ ecDNA_Status,
  data = tbl_latency_analysis
)

cox_latency_uni_cancertype <- coxph(
  Surv(Latency_Months, Latency_event) ~ CancerType,
  data = tbl_latency_analysis
)

cox_latency_multi_cancertype <- coxph(
  Surv(Latency_Months, Latency_event) ~ ecDNA_Status + CancerType,
  data = tbl_latency_analysis
)

cox_latency_multi_cancertype_clinical <- coxph(
  Surv(Latency_Months, Latency_event) ~ ecDNA_Status + CancerType + Sex,
  data = tbl_latency_analysis_clinical
)

cox_zph_latency_uni_ecDNA <- cox.zph(cox_latency_uni_ecDNA)
cox_zph_latency_multi_cancertype <- cox.zph(cox_latency_multi_cancertype)
cox_zph_latency_multi_cancertype_clinical <- cox.zph(cox_latency_multi_cancertype_clinical)

# ============================================================
# 10.4) ecDNA status versus cancer type association
# This is NOT Cox; use contingency analysis
# ============================================================
tbl_ecdna_cancertype <- tbl_analysis_combined %>%
  filter(
    !is.na(ecDNA_Status),
    !is.na(CancerType)
  ) %>%
  mutate(
    ecDNA_Status = factor(ecDNA_Status, levels = c("ecDNA-", "ecDNA+")),
    CancerType = factor(str_trim(as.character(CancerType)))
  )

tab_ecdna_cancertype <- table(
  tbl_ecdna_cancertype$CancerType,
  tbl_ecdna_cancertype$ecDNA_Status
)

tbl_ecdna_cancertype_summary <- tbl_ecdna_cancertype %>%
  count(CancerType, ecDNA_Status, name = "n") %>%
  group_by(CancerType) %>%
  mutate(
    prop = n / sum(n),
    pct = 100 * prop
  ) %>%
  ungroup()

write.csv(
  as.data.frame.matrix(tab_ecdna_cancertype),
  file.path(DIR_TABLE, "figure5_ecDNA_by_cancer_type_count_table.csv")
)

write.csv(
  tbl_ecdna_cancertype_summary,
  file.path(DIR_TABLE, "figure5_ecDNA_by_cancer_type_summary.csv"),
  row.names = FALSE
)

chisq_out <- suppressWarnings(chisq.test(tab_ecdna_cancertype))
fisher_out <- tryCatch(
  fisher.test(tab_ecdna_cancertype),
  error = function(e) NULL
)

sink(file.path(DIR_TABLE, "figure5_ecDNA_vs_CancerType_test.txt"))
cat("=== CancerType x ecDNA_Status count table ===\n")
print(tab_ecdna_cancertype)

cat("\n=== Chi-square test ===\n")
print(chisq_out)

cat("\n=== Fisher's exact test ===\n")
print(fisher_out)
sink()

# ============================================================
# 10.5) Export pooled Cox model results
# ============================================================
tbl_cox_all <- bind_rows(
  tidy_cox_model2(cox_os_uni_ecDNA, "OS: univariable ecDNA"),
  tidy_cox_model2(cox_os_uni_cancertype, "OS: univariable cancer type"),
  tidy_cox_model2(cox_os_multi_cancertype, "OS: ecDNA + cancer type"),
  tidy_cox_model2(cox_os_multi_cancertype_clinical, "OS: ecDNA + cancer type + sex"),
  tidy_cox_model2(cox_latency_uni_ecDNA, "Latency: univariable ecDNA"),
  tidy_cox_model2(cox_latency_uni_cancertype, "Latency: univariable cancer type"),
  tidy_cox_model2(cox_latency_multi_cancertype, "Latency: ecDNA + cancer type"),
  tidy_cox_model2(cox_latency_multi_cancertype_clinical, "Latency: ecDNA + cancer type + sex")
)

write.csv(
  tbl_cox_all,
  file.path(DIR_TABLE, "figure5_cox_results_all.csv"),
  row.names = FALSE
)

tbl_cox_ecDNA_only <- tbl_cox_all %>%
  filter(str_detect(term, "ecDNA_Status"))

write.csv(
  tbl_cox_ecDNA_only,
  file.path(DIR_TABLE, "figure5_cox_results_ecDNA_only.csv"),
  row.names = FALSE
)

# ============================================================
# 10.6) Export pooled Cox model summaries
# ============================================================
sink(file.path(DIR_TABLE, "figure5_cox_summary.txt"))

cat("=== OS: univariable ecDNA ===\n")
print(summary(cox_os_uni_ecDNA))

cat("\n=== OS: univariable cancer type ===\n")
print(summary(cox_os_uni_cancertype))

cat("\n=== OS: ecDNA + cancer type ===\n")
print(summary(cox_os_multi_cancertype))

cat("\n=== OS: ecDNA + cancer type + sex ===\n")
print(summary(cox_os_multi_cancertype_clinical))

cat("\n=== Latency: univariable ecDNA ===\n")
print(summary(cox_latency_uni_ecDNA))

cat("\n=== Latency: univariable cancer type ===\n")
print(summary(cox_latency_uni_cancertype))

cat("\n=== Latency: ecDNA + cancer type ===\n")
print(summary(cox_latency_multi_cancertype))

cat("\n=== Latency: ecDNA + cancer type + sex ===\n")
print(summary(cox_latency_multi_cancertype_clinical))

cat("\n=== PH assumption: OS univariable ecDNA ===\n")
print(cox_zph_os_uni_ecDNA)

cat("\n=== PH assumption: OS ecDNA + cancer type ===\n")
print(cox_zph_os_multi_cancertype)

cat("\n=== PH assumption: Latency univariable ecDNA ===\n")
print(cox_zph_latency_uni_ecDNA)

cat("\n=== PH assumption: Latency ecDNA + cancer type ===\n")
print(cox_zph_latency_multi_cancertype)

sink()

# ============================================================
# 10.7) ecDNA-only summary forest across pooled models
# ============================================================
plot_ecDNA_forest_detail <- make_ecDNA_forest_with_text(
  tbl = tbl_cox_ecDNA_only,
  title_letter = "Sup"
)

ggsave(
  file.path(DIR_PLOT, "figure5dsup_cox_forest.pdf"),
  plot_ecDNA_forest_detail,
  width = 13,
  height = 5.0
)

ggsave(
  file.path(DIR_PLOT, "figure5dsup_cox_forest.png"),
  plot_ecDNA_forest_detail,
  width = 13,
  height = 5.0,
  dpi = 300
)


format_ecDNA_forest_table.1 <- function(tbl) {
  tbl %>%
    mutate(
      model_label = as.character(term),
      hr_text = sprintf("%.2f (%.2f-%.2f)", HR, CI_low, CI_high),
      p_text = case_when(
        is.na(P_value) ~ "",
        P_value < 0.001 ~ "<0.001",
        TRUE ~ sprintf("%.3f", P_value)
      )
    ) %>%
    mutate(
      row_id = rev(seq_len(n()))
    )
}

make_ecDNA_forest_with_text.1 <- function(tbl, title_letter = "d") {
  df_plot <- format_ecDNA_forest_table.1(tbl)
  
  p_table <- ggplot(df_plot) +
    geom_text(aes(x = 0.0, y = row_id, label = model_label), hjust = 0, size = 4) +
    geom_text(aes(x = 5.8, y = row_id, label = hr_text), hjust = 0, size = 4) +
    geom_text(aes(x = 9.0, y = row_id, label = p_text), hjust = 0, size = 4) +
    annotate("text", x = 0.0, y = max(df_plot$row_id) + 0.8,
             label = title_letter, fontface = "bold", hjust = 0, size = 6) +
    annotate("text", x = 0.4, y = max(df_plot$row_id) + 0.8,
             label = "Model", fontface = "bold", hjust = 0, size = 4.2) +
    annotate("text", x = 5.8, y = max(df_plot$row_id) + 0.8,
             label = "HR (95% CI)", fontface = "bold", hjust = 0, size = 4.2) +
    annotate("text", x = 9.0, y = max(df_plot$row_id) + 0.8,
             label = "P", fontface = "bold", hjust = 0, size = 4.2) +
    xlim(0, 10) +
    ylim(0.5, max(df_plot$row_id) + 1.2) +
    theme_void()
  
  p_forest <- ggplot(df_plot, aes(x = HR, y = row_id)) +
    geom_vline(xintercept = 1, linetype = 2, color = "grey50") +
    geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.15, color = "#8B0000") +
    geom_point(size = 3, color = "#8B0000") +
    scale_x_log10() +
    labs(
      x = "Hazard ratio for ecDNA+ (log scale)",
      y = NULL
    ) +
    ylim(0.5, max(df_plot$row_id) + 1.2) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank()
    )
  
  p_table + p_forest + patchwork::plot_layout(widths = c(2.3, 1.4))
}

tbl_cox_ecDNA_only.1 <- tbl_cox_all %>%
  filter(model %in% c("OS: univariable ecDNA",
                      "OS: ecDNA + cancer type + sex",
                      "Latency: univariable ecDNA",
                      "Latency: ecDNA + cancer type + sex"))

plot_ecDNA_forest_detail.1 <- make_ecDNA_forest_with_text.1(
  tbl = tbl_cox_ecDNA_only.1,
  title_letter = "Sup"
)


ggsave(
  file.path(DIR_PLOT, "figure5dsup_cox_forest_detail.pdf"),
  plot_ecDNA_forest_detail.1,
  width = 13,
  height =10
)

ggsave(
  file.path(DIR_PLOT, "figure5dsup_cox_forest_detail.png"),
  plot_ecDNA_forest_detail.1,
  width = 13,
  height = 10,
  dpi = 300
)


# ============================================================
# 10.8) Publication-style full forest for OS
# Goal:
#   Draw the full multivariable OS forest using ecDNA status,
#   sex, and all cancer types.
# ============================================================

# ----------------------------
# 10.8.1) Choose reference cancer type explicitly
# IMPORTANT:
#   The displayed reference level must match the Cox model reference.
# ----------------------------
REF_CANCER_TYPE_OS <- "LUAD"

# Use the same analysis dataset as the model
df_model_os <- tbl_os_analysis_clinical %>%
  mutate(
    ecDNA_Status = factor(ecDNA_Status, levels = c("ecDNA-", "ecDNA+")),
    Sex = factor(Sex, levels = c("male", "female")),
    CancerType = str_trim(as.character(CancerType)),
    CancerType = factor(CancerType)
  ) %>%
  filter(!is.na(ecDNA_Status), !is.na(Sex), !is.na(CancerType)) %>%
  mutate(
    CancerType = relevel(CancerType, ref = REF_CANCER_TYPE_OS)
  )

# Refit the model so that the forest reference matches the actual model reference
cox_os_full_for_forest <- coxph(
  Surv(OS_Months, OS_event) ~ ecDNA_Status + Sex + CancerType,
  data = df_model_os
)

tbl_model_os <- tidy_cox_model2(
  cox_os_full_for_forest,
  model_name = "OS: ecDNA + sex + cancer type (all cancer types)"
)

write.csv(
  tbl_model_os,
  file.path(DIR_TABLE, "figure5_OS_full_model_all_cancer_types.csv"),
  row.names = FALSE
)



# ============================================================
# 10.9) Cancer-type-specific OS analysis
# Goal:
#   Evaluate the ecDNA effect on OS separately within each cancer type.
# ============================================================
# ============================================================
# 10.9b) Cancer-type-specific latency analysis
# Goal:
#   Evaluate the ecDNA effect on latency separately within each cancer type.
# ============================================================
tbl_latency_by_cancer <- tbl_analysis_combined %>%
  filter(
    !is.na(Latency_Months),
    !is.na(ecDNA_Status),
    !is.na(CancerType)
  ) %>%
  mutate(
    Latency_event = 1,
    CancerType = str_trim(as.character(CancerType)),
    ecDNA_Status = factor(ecDNA_Status, levels = c("ecDNA-", "ecDNA+"))
  )

tbl_stratified_latency <- run_ecDNA_cox_within_cancer(
  df = tbl_latency_by_cancer,
  time_col = "Latency_Months",
  event_col = "Latency_event",
  cancer_col = "CancerType",
  min_n = MIN_N_PER_CANCER_FOR_KM,
  min_events = MIN_EVENTS_FOR_COX,
  min_pos = MIN_ECDNA_POS_FOR_COX
) %>%
  arrange(P_value, desc(HR))

write.csv(
  tbl_stratified_latency,
  file.path(DIR_TABLE, "figure5_Latency_ecDNA_cox_within_cancer_type.csv"),
  row.names = FALSE
)


# ============================================================
# 10.10) Supplementary faceted KM curves by cancer type
# Goal:
#   Visualize OS differences by ecDNA status within eligible cancer types.
# ============================================================
tbl_os_facet_km <- make_surv_facet_data(tbl_analysis_combined, "OS_Months", "OS_event", max_time = 60)

plot_os_km_facet <- plot_facet_km_manual(
  df = tbl_os_facet_km,
  time_col = "time_plot",
  event_col = "event_plot",
  xlab = "OS (Months)",
  file_stub = "figure5_OS_KM_by_cancer_type"
)

# ============================================================
# 10.11) Cancer-type-level ecDNA prevalence versus survival
# Goal:
#   Explore whether cancer types with higher ecDNA prevalence also show
#   shorter median OS or latency.
# ============================================================
tbl_cancer_level <- tbl_analysis_combined %>%
  filter(!is.na(CancerType), !is.na(ecDNA_Status)) %>%
  mutate(CancerType = str_trim(as.character(CancerType))) %>%
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

write.csv(
  tbl_cancer_level,
  file.path(DIR_TABLE, "figure5_cancer_type_ecDNA_prevalence_vs_survival.csv"),
  row.names = FALSE
)

cor_os <- safe_cor_test(tbl_cancer_level$ecDNA_pct, tbl_cancer_level$median_OS, method = "spearman")
cor_latency <- safe_cor_test(tbl_cancer_level$ecDNA_pct, tbl_cancer_level$median_Latency, method = "spearman")
write.csv(cor_os, file.path(DIR_TABLE, "figure5_cor_ecDNApct_vs_medianOS.csv"), row.names = FALSE)
write.csv(cor_latency, file.path(DIR_TABLE, "figure5_cor_ecDNApct_vs_medianLatency.csv"), row.names = FALSE)

# Sensitivity analysis:
# repeat the cancer-type-level correlation after removing OV and SKCM.
tbl_cancer_level_rm_OV_SKCM <- tbl_cancer_level %>% filter(CancerType != "OV" & CancerType != "SKCM")
cor_os.rmCT <- safe_cor_test(tbl_cancer_level_rm_OV_SKCM$ecDNA_pct, tbl_cancer_level_rm_OV_SKCM$median_OS, method = "spearman")

cor_os_spr <- cor.test(
  tbl_cancer_level$ecDNA_pct,
  tbl_cancer_level$median_OS,alternative = 'greater',
  method = "spearman"
)

plot_ecDNA_vs_medianOS.p1 <- ggplot(
  tbl_cancer_level,
  aes(x = ecDNA_pct, y = median_OS, label = CancerType)
) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7, linetype = 2) +
  geom_text_repel(size = 3.5, max.overlaps = 30) +
  labs(
    title = "Supplementary",
    x = "ecDNA-positive cases by cancer type (%)",
    y = "Median OS by cancer type (months)",
    subtitle = label_cor_text(cor_os)
  ) +
  theme_classic(base_size = 12)

plot_ecDNA_vs_medianOS.p2 <- ggplot(
  tbl_cancer_level_rm_OV_SKCM,
  aes(x = ecDNA_pct, y = median_OS, label = CancerType)
) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7, linetype = 2) +
  geom_text_repel(size = 3.5, max.overlaps = 30) +
  labs(
    title = "Supplementary",
    x = "ecDNA-positive cases by cancer type (%)",
    y = "Median OS by cancer type (months)",
    subtitle = label_cor_text(cor_os.rmCT)
  ) +
  theme_classic(base_size = 12)
plot_ecDNA_vs_medianOS <- plot_ecDNA_vs_medianOS.p1 + plot_ecDNA_vs_medianOS.p2
ggsave((file.path(DIR_PLOT, "figure5_OS_ecDNA_spearman_by_cancer.pdf")), plot_ecDNA_vs_medianOS, width = 10, height = 5.2)
ggsave((file.path(DIR_PLOT, "figure5_OS_ecDNA_spearman_by_cancer.png")), plot_ecDNA_vs_medianOS, width = 10, height = 5.2, dpi = 300)

cor_latency.rmCT <- safe_cor_test(tbl_cancer_level_rm_OV_SKCM$ecDNA_pct, tbl_cancer_level_rm_OV_SKCM$median_Latency, method = "spearman")
cor_latency_spr <- cor.test(
  tbl_cancer_level_rm_OV_SKCM$ecDNA_pct,
  tbl_cancer_level_rm_OV_SKCM$median_Latency,alternative = 'less',
  method = "spearman"
)

plot_ecDNA_vs_medianLatency.p1 <- ggplot(
  tbl_cancer_level,
  aes(x = ecDNA_pct, y = median_Latency, label = CancerType)
) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7, linetype = 2) +
  geom_text_repel(size = 3.5, max.overlaps = 30) +
  labs(
    title = "Supplementary",
    x = "ecDNA-positive cases by cancer type (%)",
    y = "Median latency by cancer type (months)",
    subtitle = label_cor_text(cor_latency)
  ) +
  theme_classic(base_size = 12)

plot_ecDNA_vs_medianLatency.p2 <- ggplot(
  tbl_cancer_level_rm_OV_SKCM,
  aes(x = ecDNA_pct, y = median_Latency, label = CancerType)
) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7, linetype = 2) +
  geom_text_repel(size = 3.5, max.overlaps = 30) +
  labs(
    title = "Supplementary",
    x = "ecDNA-positive cases by cancer type (%)",
    y = "Median latency by cancer type (months)",
    subtitle = label_cor_text(cor_latency.rmCT)
  ) +
  theme_classic(base_size = 12)

plot_ecDNA_vs_medianLatency <- plot_ecDNA_vs_medianLatency.p1 + plot_ecDNA_vs_medianLatency.p2

ggsave(file.path(DIR_PLOT, "figure5_ecDNApct_vs_medianLatency_by_cancer.pdf"),
       plot_ecDNA_vs_medianLatency, width = 10, height = 5.5)
ggsave(file.path(DIR_PLOT, "figure5_ecDNApct_vs_medianLatency_by_cancer.png"),
       plot_ecDNA_vs_medianLatency, width = 10, height = 5.5, dpi = 300)



# ============================================================
# 11) Assemble and export figure panels
# ============================================================
km_os$plot      <- km_os$plot + coord_cartesian(xlim = KM_XLIM, clip = "on")
km_latency$plot <- km_latency$plot + coord_cartesian(xlim = KM_XLIM, clip = "on")

fig5_main <- (plot_latency | km_os$plot) / (km_latency$plot | plot_ecDNA_forest_detail)

ggsave(file.path(DIR_PLOT, "figure5_main.pdf"), fig5_main, width = 12.5, height = 9.5)
ggsave(file.path(DIR_PLOT, "figure5_main.png"), fig5_main, width = 12.5, height = 9.5, dpi = 300)

fig5_supp <- (plot_os | plot_os_facet) / (plot_latency_facet | plot_ecDNA_forest_detail.1)

ggsave(file.path(DIR_PLOT, "figure5_supplementary.pdf"), fig5_supp, width = 12.5, height = 10)
ggsave(file.path(DIR_PLOT, "figure5_supplementary.png"), fig5_supp, width = 12.5, height = 10, dpi = 300)

ggsave(file.path(DIR_PLOT, "figure5_KM_OS_curve.pdf"), km_os$plot, width = 6.2, height = 5.2)
ggsave(file.path(DIR_PLOT, "figure5_KM_OS_risktable.pdf"), km_os$table, width = 6.2, height = 2.2)

ggsave(file.path(DIR_PLOT, "figure5_latency_facet_2group.pdf"), plot_latency_facet, width = 12, height = 7)
ggsave(file.path(DIR_PLOT, "figure5_latency_facet_2group.png"), plot_latency_facet, width = 12, height = 7, dpi = 300)

ggsave(file.path(DIR_PLOT, "figure5_OS_facet_2group.pdf"), plot_os_facet, width = 12, height = 7)
ggsave(file.path(DIR_PLOT, "figure5_OS_facet_2group.png"), plot_os_facet, width = 12, height = 7, dpi = 300)

ggsave(file.path(DIR_PLOT, "figure5_latency_facet_2group_dualStatusOnly.pdf"), plot_latency_facet_dual, width = 12, height = 7)
ggsave(file.path(DIR_PLOT, "figure5_latency_facet_2group_dualStatusOnly.png"), plot_latency_facet_dual, width = 12, height = 7, dpi = 300)

ggsave(file.path(DIR_PLOT, "figure5_OS_facet_2group_dualStatusOnly.pdf"), plot_os_facet_dual, width = 12, height = 7)
ggsave(file.path(DIR_PLOT, "figure5_OS_facet_2group_dualStatusOnly.png"), plot_os_facet_dual, width = 12, height = 7, dpi = 300)

ggsave(file.path(DIR_PLOT, "figure5_latency_KM.pdf"), km_latency$plot, width = 6, height = 5)
ggsave(file.path(DIR_PLOT, "figure5_latency_risktable.pdf"), km_latency$table, width = 6, height = 3)

ggsave(file.path(DIR_PLOT, "figure5a_latency_box_2group.pdf"), plot_latency, width = 4, height = 4.5)
ggsave(file.path(DIR_PLOT, "figure5b_OS_KM_2group.pdf"), km_os$plot, width = 6.2, height = 5.2)