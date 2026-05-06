# ============================================================
# 0) Packages
# ============================================================
library(tidyverse)
library(data.table)
library(readxl)
library(ggpubr)
library(ggbeeswarm)
library(patchwork)
library(CINmetrics)
library(scales)
library(ggrepel)

# ============================================================
# 1) Input configuration and output directories
# ============================================================
input_dir <- "Inputs"
output_dir <- "Outputs"

DIR_OUT <- file.path(output_dir, "Figure2_genomics")
DIR_RDS <- file.path(DIR_OUT, "rds")
DIR_TABLE <- file.path(DIR_OUT, "tables")
DIR_PLOT <- file.path(DIR_OUT, "plots")
DIR_PUB_LOCAL <- file.path(DIR_PLOT, "pub_local")

purrr::walk(c(DIR_OUT, DIR_RDS, DIR_TABLE, DIR_PLOT, DIR_PUB_LOCAL), ~{
  dir.create(.x, showWarnings = FALSE, recursive = TRUE)
})

PATH_SAMPLE_INFO <- file.path(input_dir, "metadata", "sample_metadata.xlsx")
SHEET_SAMPLE_INFO <- "ecDNA WGS Sample Sequencing Statistics-refined"

PATH_ECDNA_GENE <- file.path(input_dir, "database", "ecDNA_genelist_all.rds")
PATH_ONCOGENE <- file.path(input_dir, "database", "combined_oncogene_list.txt")
PATH_PUBLIC_MUTATION <- file.path(input_dir, "database", "pancan_pcawg_2020_clinical_data.tsv")

PATH_WES_DIR <- file.path(input_dir, "mutation", "vcf_output_WES")
PATH_WGS_GLOB <- file.path(input_dir, "mutation", "WGS_batch?_CNVkit_analysis_20260*", "*.maf")
PATH_CNVKIT_DIR <- file.path(input_dir, "cnv_kit")
PATH_FEATURE_DIR <- file.path(input_dir, "feature_basic_properties_copy")
PATH_WGS_GENELIST_DIR <- file.path(input_dir, "WGS_genelist")

# ============================================================
# 2) Constants, colors, and plotting theme
# ============================================================
KEEP_CLASSES <- c(
  "Missense_Mutation",
  "Nonsense_Mutation",
  "Frame_Shift_Del",
  "Frame_Shift_Ins",
  "In_Frame_Del",
  "In_Frame_Ins",
  "Splice_Site",
  "Nonstop_Mutation",
  "DE_NOVO_START_OUT_FRAME",
  "DE_NOVO_START_IN_FRAME",
  "START_CODON_SNP",
  "Translation_Start_Site",
  "START_CODON_INS"
)

KEEP_NO_PROTEIN_CHANGE_CLASSES <- c(
  "Splice_Site",
  "DE_NOVO_START_IN_FRAME",
  "DE_NOVO_START_OUT_FRAME",
  "START_CODON_SNP",
  "Translation_Start_Site",
  "START_CODON_INS"
)

ECDNA_SAMPLE_COLS <- c(
  "ecDNA-" = "grey70",
  "ecDNA+" = "#D62728"
)

ECDNA_SEGMENT_COLS <- c(
  "Other focal amplification feature" = "grey70",
  "ecDNA feature" = "#D62728"
)

ECDNA_GENE_COLS <- c(
  "Other focal CNA" = "grey70",
  "ecDNA" = "#D62728"
)

ECDNA_AMP_COLS <- c(
  "Other focal amplicon" = "grey70",
  "ecDNA amplicon" = "#D62728"
)

CANCER_COLS <- c(
  "Lung Cancer"          = "#A38A77",
  "Breast Cancer"        = "#FFB6C1",
  "Gastric Cancer"       = "#984EA3",
  "Renal Cell Carcinoma" = "#FF7F00",
  "Melanoma"             = "#A65628",
  "Colorectal Cancer"    = "#009FE8",
  "Cervical Cancer"      = "#542788",
  "Ovarian Cancer"       = "#C51B7D",
  "Esophageal Cancer"    = "#008B8B"
)

theme_fig2 <- theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", size = 13)
  )

format_fisher_label <- function(ft) {
  paste0(
    "Fisher's exact p = ", signif(ft$p.value, 3),
    "\nOR = ", format(unname(ft$estimate), digits = 3),
    " (95% CI ",
    format(ft$conf.int[1], digits = 3), "-",
    format(ft$conf.int[2], digits = 3), ")"
  )
}

# ============================================================
# 3) Helper functions
# ============================================================
clean_missing_strings <- function(df) {
  df[df == ""] <- NA
  df[df == "__UNKNOWN__"] <- NA
  df[df == "."] <- NA
  df
}

as_logical_ecDNA <- function(x) {
  x %in% c(TRUE, "TRUE", "True", "T", "1", 1)
}

standardize_sample_id <- function(x) {
  x %>%
    basename() %>%
    str_remove("\\.cnr$") %>%
    str_remove("-.*$")
}

read_sample_info <- function(path_xlsx, sheet) {
  readxl::read_excel(path_xlsx, sheet = sheet) %>%
    filter(used == "yes") %>%
    transmute(
      Sample = as.character(Sample),
      Cancer_type = str_to_title(as.character(Cancer_type)),
      Pri_Met = as.character(Pri_Met),
      Seq_type = as.character(Seq_type),
      Pair = na_if(trimws(as.character(Pair)), "")
    ) %>%
    mutate(
      pair_type = case_when(
        is.na(Pair) ~ "Unpaired",
        str_detect(Pair, "Primary/Metastatic/Recurrent") ~ "Triad",
        TRUE ~ "Paired"
      )
    ) %>%
    distinct()
}


read_oncogene_list <- function(path_txt) {
  read.table(
    path_txt,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE
  ) %>%
    pull(1) %>%
    unique()
}

read_maf_files <- function(maf_files, seq_type) {
  if (length(maf_files) == 0) return(tibble())

  maf_list <- lapply(maf_files, function(f) {
    sample_id <- basename(f) %>%
      str_remove("_somatic_PASS_ONLY\\.maf$")

    dat <- fread(
      f,
      sep = "\t",
      header = TRUE,
      data.table = FALSE,
      fill = TRUE,
      quote = "",
      skip = "Hugo_Symbol"
    )

    dat[] <- lapply(dat, as.character)

    dat %>%
      mutate(
        Sample = sample_id,
        Seq_type = seq_type
      )
  })

  bind_rows(maf_list) %>%
    clean_missing_strings()
}

filter_maf_tumor_only <- function(
    maf_all,
    keep_classes,
    keep_no_protein_change_classes,
    min_alt = 5,
    min_depth = 15,
    min_vaf = 0.10,
    max_pop_af = 0.001
) {
  maf_all %>%
    filter(Variant_Classification %in% keep_classes) %>%
    mutate(
      Protein_Change = na_if(trimws(as.character(Protein_Change)), ""),
      Protein_Change = na_if(Protein_Change, "__UNKNOWN__"),
      dbSNP_COMMON = trimws(as.character(dbSNP_COMMON)),
      dbSNP_COMMON = na_if(dbSNP_COMMON, ""),
      dbSNP_COMMON = na_if(dbSNP_COMMON, "NA"),
      dbSNP_COMMON = na_if(dbSNP_COMMON, "__UNKNOWN__"),
      t_alt_count = suppressWarnings(as.numeric(t_alt_count)),
      t_ref_count = suppressWarnings(as.numeric(t_ref_count)),
      tumor_f = suppressWarnings(as.numeric(tumor_f)),
      ClinVar_VCF_AF_TGP = suppressWarnings(as.numeric(ClinVar_VCF_AF_TGP)),
      ClinVar_VCF_AF_EXAC = suppressWarnings(as.numeric(ClinVar_VCF_AF_EXAC)),
      ClinVar_VCF_AF_ESP = suppressWarnings(as.numeric(ClinVar_VCF_AF_ESP))
    ) %>%
    mutate(
      tumor_depth = t_alt_count + t_ref_count,
      tumor_vaf_calc = if_else(
        !is.na(t_alt_count) & !is.na(t_ref_count) & tumor_depth > 0,
        t_alt_count / tumor_depth,
        NA_real_
      ),
      tumor_vaf = coalesce(tumor_vaf_calc, tumor_f)
    ) %>%
    filter(
      !is.na(t_alt_count),
      !is.na(t_ref_count),
      t_alt_count >= min_alt,
      tumor_depth >= min_depth,
      !is.na(tumor_vaf),
      tumor_vaf >= min_vaf
    ) %>%
    filter(is.na(dbSNP_COMMON) | dbSNP_COMMON == "0") %>%
    filter(
      is.na(ClinVar_VCF_AF_TGP)  | ClinVar_VCF_AF_TGP  < max_pop_af,
      is.na(ClinVar_VCF_AF_EXAC) | ClinVar_VCF_AF_EXAC < max_pop_af,
      is.na(ClinVar_VCF_AF_ESP)  | ClinVar_VCF_AF_ESP  < max_pop_af
    ) %>%
    filter(
      (!is.na(Protein_Change) & Protein_Change != "") |
        Variant_Classification %in% keep_no_protein_change_classes
    )
}

choose_one_seq_per_sample <- function(df, priority = c("WES", "WGS")) {
  priority_tbl <- tibble(
    Seq_type = priority,
    seq_priority = seq_along(priority)
  )

  df %>%
    left_join(priority_tbl, by = "Seq_type") %>%
    mutate(seq_priority = replace_na(seq_priority, 999L)) %>%
    group_by(Sample) %>%
    mutate(min_priority = min(seq_priority, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(seq_priority == min_priority) %>%
    select(-seq_priority, -min_priority)
}

make_ec_sample_table <- function(ec_raw, sample_info) {
  ec_summary <- ec_raw %>%
    mutate(
      Sample = as.character(Sample),
      Seq_type = as.character(Seq_type),
      gene = as.character(gene),
      gene_cn = as.numeric(gene_cn),
      is_ecDNA = as_logical_ecDNA(is_ecDNA)
    ) %>%
    group_by(Sample, Seq_type) %>%
    summarise(
      n_ecDNA_genes = sum(is_ecDNA, na.rm = TRUE),
      max_ecDNA_gene_cn = if_else(any(is_ecDNA), max(gene_cn[is_ecDNA], na.rm = TRUE), NA_real_),
      n_gene_events = n(),
      max_gene_cn_all = max(gene_cn, na.rm = TRUE),
      ecDNA_status = if_else(any(is_ecDNA), "ecDNA+", "ecDNA-"),
      ecDNA_binary = as.integer(any(is_ecDNA)),
      .groups = "drop"
    )

  sample_info %>%
    mutate(
      Sample = as.character(Sample),
      Seq_type = as.character(Seq_type)
    ) %>%
    full_join(ec_summary, by = c("Sample", "Seq_type"))
}

make_ec_event_table <- function(ec_raw, sample_info, oncogenes) {
  ec_filtered <- ec_raw %>%
    mutate(
      Sample = as.character(Sample),
      Seq_type = as.character(Seq_type),
      gene = as.character(gene),
      gene_cn = as.numeric(gene_cn),
      is_ecDNA = as_logical_ecDNA(is_ecDNA)
    ) %>%
    inner_join(sample_info %>% select(Sample, Seq_type), by = c("Sample", "Seq_type"))

  ec_event <- ec_filtered %>%
    mutate(
      ecDNA_segment = if_else(is_ecDNA, "ecDNA", "Other focal CNA"),
      ecDNA_segment = factor(ecDNA_segment, levels = c("Other focal CNA", "ecDNA"))
    )

  ec_event_onco <- ec_event %>%
    filter(gene %in% oncogenes)

  list(all = ec_event, oncogene = ec_event_onco)
}

make_mutation_sample_table <- function(maf_func_all, oncogenes, ec_sample) {
  tbl_mut <- maf_func_all %>%
    group_by(Sample, Seq_type) %>%
    summarise(
      n_func_mut = n(),
      n_mut_gene = n_distinct(Hugo_Symbol),
      .groups = "drop"
    ) %>%
    rename(mutseq = Seq_type)

  tbl_onco_mut <- maf_func_all %>%
    filter(Hugo_Symbol %in% oncogenes) %>%
    group_by(Sample, Seq_type) %>%
    summarise(
      n_oncogene_mut = n(),
      n_oncogene_mut_gene = n_distinct(Hugo_Symbol),
      .groups = "drop"
    ) %>%
    rename(mutseq = Seq_type)

  tbl_ec_status <- ec_sample %>%
    distinct(Sample, ecDNA_status)

  tbl_mut %>%
    left_join(tbl_onco_mut, by = c("Sample", "mutseq")) %>%
    mutate(
      n_oncogene_mut = replace_na(n_oncogene_mut, 0),
      n_oncogene_mut_gene = replace_na(n_oncogene_mut_gene, 0),
      TMB_func_mut = n_func_mut / 45,
      TMB_oncogene_mut = n_oncogene_mut / 45
    ) %>%
    left_join(tbl_ec_status, by = "Sample") %>%
    mutate(ecDNA_status = replace_na(ecDNA_status, "ecDNA-"))
}

plot_box_compare <- function(df, xvar, yvar, ylab, fill_cols,
                             label_y = NULL, alternative = "greater") {
  p <- ggplot(df, aes(x = .data[[xvar]], y = .data[[yvar]], fill = .data[[xvar]])) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 1, staplewidth = 0.5) +
    geom_jitter(width = 0.15, height = 0, alpha = 0.7, size = 1.8) +
    scale_fill_manual(values = fill_cols) +
    labs(x = NULL, y = ylab) +
    theme_fig2

  if (!is.null(label_y)) {
    p <- p + stat_compare_means(
      method = "wilcox.test",
      label = "p.format",
      label.y = label_y,
      method.args = list(alternative = alternative)
    )
  }

  p
}

plot_gene_cn_compare <- function(df, xvar, yvar = "gene_cn", ylab,
                                 fill_cols, ymax = NULL) {
  if (is.null(ymax)) {
    ymax <- quantile(df[[yvar]], 0.99, na.rm = TRUE) * 1.05
  }

  ggplot(df, aes(x = .data[[xvar]], y = log10(.data[[yvar]]), fill = .data[[xvar]])) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 1, staplewidth = 0.5) +
    geom_jitter(width = 0.15, height = 0, size = 0.5, alpha = 0.2) +
    stat_compare_means(
      method = "wilcox.test",
      label = "p.format",
      label.y = ymax
    ) +
    scale_fill_manual(values = fill_cols) +
    coord_cartesian(ylim = c(0, ymax)) +
    labs(x = NULL, y = ylab) +
    theme_fig2
}

read_props_safe <- function(f) {
  if (!file.exists(f) || file.size(f) == 0) return(NULL)
  df <- tryCatch(readr::read_tsv(f, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df %>% mutate(source_file = basename(f))
}

read_genelist_safe <- function(x) {
  if (!file.exists(x) || file.size(x) == 0) return(NULL)

  tmp <- tryCatch(
    read.table(
      x,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE,
      quote = "",
      comment.char = "",
      fill = TRUE
    ),
    error = function(e) NULL
  )

  if (is.null(tmp) || nrow(tmp) == 0) return(NULL)

  tmp[] <- lapply(tmp, as.character)
  tmp %>% mutate(file_source = basename(x))
}

# ============================================================
# 4) Read core inputs
# ============================================================
tbl_sample_info <- read_sample_info(PATH_SAMPLE_INFO, SHEET_SAMPLE_INFO)
stopifnot(nrow(tbl_sample_info) > 0)
saveRDS(tbl_sample_info, file.path(DIR_RDS, "sample_info_summary.rds"))

vec_oncogenes <- read_oncogene_list(PATH_ONCOGENE)
obj_ec_raw <- readRDS(PATH_ECDNA_GENE)

# ============================================================
# 5) Build analysis tables
# ============================================================

# 5.1 Mutation tables
wes_files <- list.files(
  path = PATH_WES_DIR,
  pattern = "_somatic_PASS_ONLY\\.maf$",
  full.names = TRUE
)

wgs_files <- Sys.glob(PATH_WGS_GLOB)

tbl_wes_maf_all <- read_maf_files(wes_files, seq_type = "WES")
tbl_wgs_maf_all <- read_maf_files(wgs_files, seq_type = "WGS")

tbl_wes_maf_func <- filter_maf_tumor_only(
  maf_all = tbl_wes_maf_all,
  keep_classes = KEEP_CLASSES,
  keep_no_protein_change_classes = KEEP_NO_PROTEIN_CHANGE_CLASSES,
  min_alt = 5,
  min_depth = 15,
  min_vaf = 0.10,
  max_pop_af = 0.001
)

tbl_wgs_maf_func <- filter_maf_tumor_only(
  maf_all = tbl_wgs_maf_all,
  keep_classes = KEEP_CLASSES,
  keep_no_protein_change_classes = KEEP_NO_PROTEIN_CHANGE_CLASSES,
  min_alt = 5,
  min_depth = 12,
  min_vaf = 0.10,
  max_pop_af = 0.001
)

saveRDS(tbl_wes_maf_func, file.path(DIR_RDS, "WES_maf_func_filtered.rds"))
saveRDS(tbl_wgs_maf_func, file.path(DIR_RDS, "WGS_maf_func_filtered.rds"))

tbl_maf_func_all <- bind_rows(tbl_wes_maf_func, tbl_wgs_maf_func) %>%
  mutate(
    Sample = as.character(Sample),
    Hugo_Symbol = as.character(Hugo_Symbol),
    Protein_Change = as.character(Protein_Change),
    Seq_type = as.character(Seq_type)
  ) %>%
  choose_one_seq_per_sample(priority = c("WES", "WGS")) %>%
  filter(Sample %in% tbl_sample_info$Sample)

# 5.2 ecDNA event and sample tables
tbl_ec_sample <- make_ec_sample_table(obj_ec_raw, tbl_sample_info)
tbl_ec_sample$ecDNA_status = ifelse(is.na(tbl_ec_sample$ecDNA_status), "undetected",tbl_ec_sample$ecDNA_status)
tbl_ec_sample$ecDNA_binary = ifelse(is.na(tbl_ec_sample$ecDNA_binary), 2,tbl_ec_sample$ecDNA_binary)

saveRDS(tbl_ec_sample, file.path(DIR_RDS, "ecDNA_sample_summary.rds"))

lst_ec_event <- make_ec_event_table(obj_ec_raw, tbl_sample_info, vec_oncogenes)
tbl_ec_event <- lst_ec_event$all
tbl_ec_event_onco <- lst_ec_event$oncogene

saveRDS(tbl_ec_event, file.path(DIR_RDS, "fig2_event_level_table.rds"))
saveRDS(tbl_ec_event_onco, file.path(DIR_RDS, "fig2_event_level_oncogene_table.rds"))

tbl_mut_sample <- make_mutation_sample_table(tbl_maf_func_all, vec_oncogenes, tbl_ec_sample)
saveRDS(tbl_mut_sample, file.path(DIR_RDS, "fig2_sample_level_table_Mutation.rds"))

# 5.3 FGA table
infer_seq_type_from_cnv_path <- function(x) {
  x_low <- tolower(x)
  case_when(
    str_detect(x_low, "wes") ~ "WES",
    str_detect(x_low, "wgs") ~ "WGS",
    TRUE ~ NA_character_
  )
}
cnv_files <- list.files(PATH_CNVKIT_DIR, recursive = TRUE, full.names = TRUE)

tbl_cnv <- bind_rows(lapply(cnv_files, function(f) {
  read.table(f, header = TRUE, stringsAsFactors = FALSE, sep = "\t") %>%
    mutate(
      Sample = standardize_sample_id(f),
      Seq_type = infer_seq_type_from_cnv_path(f),
      source_file = basename(f)
    )
}))

tbl_cnv <- tbl_cnv %>%
  choose_one_seq_per_sample(priority = c("WES", "WGS"))

tbl_cnv_gainloss <- tbl_cnv %>%
  filter(cn != 2)

tbl_fga_input <- tbl_cnv_gainloss %>%
  transmute(
    Sample = Sample,
    Chromosome = gsub("chr", "", chromosome),
    Start = start,
    End = end,
    Num_Probes = probes,
    Segment_Mean = log2
  ) %>%
  filter(Chromosome %in% as.character(1:22))

tbl_fga <- fga(
  tbl_fga_input,
  segmentMean = 0.2,
  numProbes = NA,
  genomeSize = 2875001522
)

write.csv(tbl_fga, file.path(DIR_TABLE, "CNV-index.csv"), row.names = FALSE)

tbl_fga_plot <- tbl_fga %>%
  left_join(
    tbl_ec_sample %>%
      mutate(SampleID = standardize_sample_id(Sample)) %>%
      distinct(SampleID, ecDNA_status),
    by = c("sample_id" = "SampleID")
  ) %>%
  mutate(ecDNA_status = replace_na(ecDNA_status, "ecDNA-"))

# 5.4 Feature-level CN table
feature_files <- list.files(
  PATH_FEATURE_DIR,
  pattern = "_feature_basic_properties\\.tsv$",
  full.names = TRUE
)

tbl_props <- purrr::map(feature_files, read_props_safe) %>%
  purrr::compact() %>%
  bind_rows()

if (nrow(tbl_props) == 0) {
  stop("No non-empty feature_basic_properties.tsv files were loaded.")
}

tbl_feature <- tbl_props %>%
  mutate(
    Sample = str_extract(feature_ID, "^[^_]+"),
    is_ecDNA = str_detect(feature_ID, "_ecDNA_")
  ) %>%
  left_join(
    tbl_sample_info %>%
      filter(Seq_type == "WGS") %>%
      select(Sample, Cancer_type, Pri_Met),
    by = "Sample"
  ) %>%
  filter(!is.na(Pri_Met), !is.na(Cancer_type)) %>%
  mutate(
    Pri_Met = factor(Pri_Met, levels = c("Pri", "Met")),
    Cancer_type_plot = factor(Cancer_type, levels = names(CANCER_COLS)),
    ecDNA_group = if_else(is_ecDNA, "ecDNA", "Other focal CNA"),
    ecDNA_group = factor(ecDNA_group, levels = c("Other focal CNA", "ecDNA"))
  )

saveRDS(tbl_feature, file.path(DIR_RDS, "feature_df.rds"))

# 5.5 WGS amplicon oncogene table
wgs_genelist_files <- list.files(PATH_WGS_GENELIST_DIR, pattern = "\\.tsv$", full.names = TRUE)
if (length(wgs_genelist_files) == 0) stop("No WGS genelist files found.")

df_list <- lapply(wgs_genelist_files, read_genelist_safe)
df_list <- df_list[!vapply(df_list, is.null, logical(1))]
if (length(df_list) == 0) stop("No valid genelist tables were loaded.")

tbl_wgs_gene_raw <- bind_rows(df_list)

tbl_wgs_gene <- tbl_wgs_gene_raw %>%
  mutate(
    sample_name = str_trim(as.character(sample_name)),
    amplicon_number = as.character(amplicon_number),
    feature = as.character(feature),
    gene = str_trim(as.character(gene)),
    gene_cn = suppressWarnings(as.numeric(gene_cn)),
    is_canonical_oncogene = recode(
      as.character(is_canonical_oncogene),
      "True" = "TRUE",
      "False" = "FALSE"
    )
  ) %>%
  filter(sample_name %in% tbl_sample_info$Sample[tbl_sample_info$Seq_type == "WGS"]) %>%
  filter(!is.na(gene), gene != "") %>%
  filter(!is.na(amplicon_number), amplicon_number != "")

write.csv(tbl_wgs_gene, file.path(DIR_TABLE, "WGS_genelist_cleaned.csv"), row.names = FALSE)

tbl_feature_count <- tbl_wgs_gene %>%
  count(feature, sort = TRUE)

write.csv(tbl_feature_count, file.path(DIR_TABLE, "WGS_genelist_feature_counts.csv"), row.names = FALSE)

tbl_wgs_gene2 <- tbl_wgs_gene %>%
  mutate(
    amplicon_id = paste(sample_name, amplicon_number, sep = "__"),
    is_oncogene = is_canonical_oncogene == "TRUE"
  )

tbl_amplicon_onco <- tbl_wgs_gene2 %>%
  group_by(sample_name, amplicon_id, amplicon_number) %>%
  summarise(
    amplicon_type = if_else(
      any(str_detect(feature, regex("ecDNA", ignore_case = TRUE))),
      "ecDNA amplicon",
      "Other focal amplicon"
    ),
    n_genes = n_distinct(gene),
    n_oncogene = n_distinct(gene[is_oncogene]),
    has_oncogene = as.integer(any(is_oncogene)),
    .groups = "drop"
  ) %>%
  mutate(
    amplicon_type = factor(
      amplicon_type,
      levels = c("Other focal amplicon", "ecDNA amplicon")
    )
  )

write.csv(tbl_amplicon_onco, file.path(DIR_TABLE, "amplicon_oncogene_flag_WGS.csv"), row.names = FALSE)

tbl_amplicon_tab <- tbl_amplicon_onco %>%
  count(amplicon_type, has_oncogene) %>%
  tidyr::complete(
    amplicon_type = factor(
      c("Other focal amplicon", "ecDNA amplicon"),
      levels = c("Other focal amplicon", "ecDNA amplicon")
    ),
    has_oncogene = c(0, 1),
    fill = list(n = 0)
  ) %>%
  mutate(
    oncogene_status = if_else(has_oncogene == 1, ">=1 oncogene", "0 oncogene")
  )

write.csv(tbl_amplicon_tab, file.path(DIR_TABLE, "amplicon_oncogene_2x2_table_WGS_long.csv"), row.names = FALSE)

mat_amplicon_onco <- tbl_amplicon_tab %>%
  mutate(amplicon_type = as.character(amplicon_type)) %>%
  select(-has_oncogene) %>%
  pivot_wider(names_from = oncogene_status, values_from = n) %>%
  tibble::column_to_rownames("amplicon_type") %>%
  as.matrix()

write.csv(as.data.frame(mat_amplicon_onco), file.path(DIR_TABLE, "amplicon_oncogene_2x2_table_WGS_matrix.csv"), row.names = TRUE)

res_amplicon_fisher <- fisher.test(mat_amplicon_onco)

tbl_amplicon_fisher <- tibble(
  comparison = "ecDNA amplicon vs Other focal amplicon",
  p_value = res_amplicon_fisher$p.value,
  odds_ratio = unname(res_amplicon_fisher$estimate),
  conf_low = res_amplicon_fisher$conf.int[1],
  conf_high = res_amplicon_fisher$conf.int[2]
)

write.csv(tbl_amplicon_fisher, file.path(DIR_TABLE, "amplicon_oncogene_fisher_result_WGS.csv"), row.names = FALSE)

tbl_prop_amplicon_onco <- tbl_amplicon_onco %>%
  group_by(amplicon_type) %>%
  summarise(
    n_amplicons = n(),
    n_with_oncogene = sum(has_oncogene),
    prop_with_oncogene = n_with_oncogene / n_amplicons,
    .groups = "drop"
  )

write.csv(tbl_prop_amplicon_onco, file.path(DIR_TABLE, "amplicon_oncogene_proportion_WGS.csv"), row.names = FALSE)

# 5.6 Sample-level sensitivity tables
tbl_wgs_sample_universe <- tbl_sample_info %>%
  filter(Seq_type == "WGS") %>%
  select(Sample, Cancer_type) %>%
  distinct() %>%
  left_join(
    tbl_ec_sample %>% select(Sample, ecDNA_status, ecDNA_binary) %>% distinct(),
    by = "Sample"
  ) %>%
  mutate(
    ecDNA_status = replace_na(ecDNA_status, "ecDNA-"),
    ecDNA_binary = replace_na(ecDNA_binary, 0L)
  )

tbl_sample_max_segment_cn <- tbl_feature %>%
  group_by(Sample) %>%
  summarise(
    n_features = n(),
    n_ecDNA_features = sum(is_ecDNA, na.rm = TRUE),
    max_segment_cn_raw = max(max_feature_CN, na.rm = TRUE),
    max_segment_cn_log2 = max(log2(max_feature_CN), na.rm = TRUE),
    max_segment_cn_log10 = max(log10(max_feature_CN), na.rm = TRUE),
    .groups = "drop"
  )

tbl_sample_max_oncogene_cn <- tbl_wgs_gene2 %>%
  filter(is_oncogene) %>%
  group_by(sample_name) %>%
  summarise(
    n_oncogene_instances = n(),
    n_distinct_oncogenes = n_distinct(gene),
    max_oncogene_cn_raw = max(gene_cn, na.rm = TRUE),
    max_oncogene_cn_log2 = max(log2(gene_cn), na.rm = TRUE),
    max_oncogene_cn_log10 = max(log10(gene_cn), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(Sample = sample_name)

tbl_sample_wgs_summary <- tbl_wgs_sample_universe %>%
  left_join(tbl_sample_max_segment_cn, by = "Sample") %>%
  left_join(tbl_sample_max_oncogene_cn, by = "Sample") %>%
  mutate(
    n_features = replace_na(n_features, 0L),
    n_ecDNA_features = replace_na(n_ecDNA_features, 0L),
    n_oncogene_instances = replace_na(n_oncogene_instances, 0L),
    n_distinct_oncogenes = replace_na(n_distinct_oncogenes, 0L)
  )

saveRDS(tbl_sample_wgs_summary, file.path(DIR_RDS, "sample_level_wgs_summary_for_sensitivity.rds"))
write.csv(tbl_sample_wgs_summary, file.path(DIR_TABLE, "sample_level_wgs_summary_for_sensitivity.csv"), row.names = FALSE)

# feature plot table for main panel d (2-group, feature-level)
tbl_feature_plot <- tbl_feature %>%
  mutate(
    feature_group2 = if_else(
      is_ecDNA,
      "ecDNA feature",
      "Other focal amplification feature"
    ),
    feature_group2 = factor(
      feature_group2,
      levels = c("Other focal amplification feature", "ecDNA feature")
    )
  )

# explicit join fix for gene-level CN
tbl_gene_cn <- obj_ec_raw %>%
  mutate(
    Sample = as.character(Sample),
    Seq_type = as.character(Seq_type),
    gene = as.character(gene),
    gene_cn = as.numeric(gene_cn),
    is_ecDNA = as_logical_ecDNA(is_ecDNA),
    ecDNA_group = if_else(is_ecDNA, "ecDNA", "Other focal CNA")
  ) %>%
  left_join(
    tbl_sample_info %>% select(Sample, Seq_type, Cancer_type, Pri_Met),
    by = c("Sample", "Seq_type")
  ) %>%
  filter(!is.na(gene_cn), !is.na(Cancer_type), !is.na(Pri_Met)) %>%
  mutate(
    ecDNA_group = factor(ecDNA_group, levels = c("Other focal CNA", "ecDNA")),
    Cancer_type = factor(Cancer_type, levels = names(CANCER_COLS)),
    Pri_Met = factor(Pri_Met, levels = c("Pri", "Met"))
  )

tbl_gene_cn_onco <- tbl_gene_cn %>%
  filter(is_canonical_oncogene %in% c(TRUE, "TRUE"))

# ============================================================
# 6) Main figure plots
# ============================================================
tbl_mut_sample <- tbl_mut_sample %>%
  mutate(ecDNA_status = if_else(ecDNA_status == "undetected", "ecDNA-", ecDNA_status))

plot_main_func_mut <- plot_box_compare(
  df = tbl_mut_sample,
  xvar = "ecDNA_status",
  yvar = "TMB_func_mut",
  ylab = "Burden of putative functional alterations",
  fill_cols = ECDNA_SAMPLE_COLS,
  label_y = 30
) + ggtitle("a") + ylim(0,30)

plot_main_oncogene_mut <- plot_box_compare(
  df = tbl_mut_sample,
  xvar = "ecDNA_status",
  yvar = "TMB_oncogene_mut",
  ylab = "Burden of candidate oncogene mutations",
  fill_cols = ECDNA_SAMPLE_COLS,
  label_y = 2
) + ggtitle("b") + ylim(0,2)

tbl_fga_plot$ecDNA_status = ifelse(tbl_fga_plot$ecDNA_status=='undetected', 'ecDNA-', 
                                   tbl_fga_plot$ecDNA_status)
tbl_fga_plot$ecDNA_status = factor(tbl_fga_plot$ecDNA_status, levels = c("ecDNA-", "ecDNA+"))

plot_main_gii <- plot_box_compare(
  df = tbl_fga_plot,
  xvar = "ecDNA_status",
  yvar = "fga",
  ylab = "Genomic instability index",
  fill_cols = ECDNA_SAMPLE_COLS,
  label_y = 1
) +
  ggtitle("c")

plot_main_segment_cn_2group <- ggplot(
  tbl_feature_plot,
  aes(x = feature_group2, y = log2(max_feature_CN), fill = feature_group2)
) +
  geom_boxplot(
    width = 0.5,
    outlier.shape = NA,
    alpha = 1,
    staplewidth = 0.5,
    color = "black",
    linewidth = 0.6
  ) +
  ggbeeswarm::geom_quasirandom(
    aes(color = Cancer_type_plot),
    width = 0.2,
    size = 1,
    alpha = 0.8
  ) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    label.y = 10,
    method.args = list(alternative = "greater")
  ) +
  scale_fill_manual(values = ECDNA_SEGMENT_COLS) +
  scale_color_manual(values = CANCER_COLS, drop = TRUE) +
  coord_cartesian(ylim = c(2, 10)) +
  labs(
    x = NULL,
    y = "Maximum segment copy number"
  ) +
  theme_fig2 +
  ggtitle("d")

plot_main_oncogene_cn <- plot_gene_cn_compare(
  df = tbl_gene_cn_onco,
  xvar = "ecDNA_group",
  ylab = "Copy number of oncogenes in amplified regions (log10)",
  fill_cols = ECDNA_GENE_COLS,
  ymax = 3
) +
  ggtitle("e")

plot_main_amplicon_oncogene <- ggplot(
  tbl_prop_amplicon_onco,
  aes(x = amplicon_type, y = prop_with_oncogene, fill = amplicon_type)
) +
  geom_col(width = 0.65, color = "black") +
  geom_text(
    aes(label = paste0(n_with_oncogene, "/", n_amplicons)),
    vjust = -0.3,
    size = 3.8
  ) +
  annotate(
    "text",
    x = 1.5,
    y = 1.02,
    label = format_fisher_label(res_amplicon_fisher),
    size = 3.5
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.08)) +
  scale_fill_manual(values = ECDNA_AMP_COLS) +
  labs(
    x = NULL,
    y = "Amplicons with >=1 oncogene (%)"
  ) +
  theme_fig2 +
  ggtitle("f")

# ============================================================
# 7) Supplementary figure plots
# ============================================================
ec_pos_samples <- tbl_ec_sample %>%
  filter(ecDNA_status == "ecDNA+") %>%
  pull(Sample) %>%
  unique()

tbl_feature_plot_supp <- tbl_feature %>%
  mutate(
    ecDNA_group1 = case_when(
      Sample %in% ec_pos_samples & is_ecDNA  ~ "ecDNA+ feature",
      Sample %in% ec_pos_samples & !is_ecDNA ~ "ecDNA+ non-ecDNA feature",
      !Sample %in% ec_pos_samples            ~ "ecDNA- feature",
      TRUE ~ NA_character_
    ),
    ecDNA_group1 = factor(
      ecDNA_group1,
      levels = c("ecDNA- feature", "ecDNA+ non-ecDNA feature", "ecDNA+ feature")
    )
  )

plot_supp_segment_cn_3group <- ggplot(
  tbl_feature_plot_supp,
  aes(x = ecDNA_group1, y = log2(max_feature_CN), fill = ecDNA_group1)
) +
  geom_boxplot(
    width = 0.5, outlier.shape = NA, alpha = 1,
    staplewidth = 0.5, color = "black", linewidth = 0.6
  ) +
  ggbeeswarm::geom_quasirandom(
    aes(color = Cancer_type_plot),
    width = 0.2, size = 1, alpha = 0.8
  ) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(c(1, 2), c(1, 3), c(2, 3)),
    label = "p.format",
    tip.length = 0.01,
    size = 3.5
  ) +
  scale_fill_manual(values = c(
    "ecDNA- feature" = "grey70",
    "ecDNA+ non-ecDNA feature" = "#F8766D",
    "ecDNA+ feature" = "#D62728"
  )) +
  scale_color_manual(values = CANCER_COLS, drop = TRUE) +
  labs(
    x = NULL,
    y = "Maximum segment copy number"
  ) +
  theme_fig2 +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) +
  ggtitle("a")

plot_supp_gene_cn_all <- plot_gene_cn_compare(
  df = tbl_gene_cn,
  xvar = "ecDNA_group",
  ylab = "Copy number of genes in amplified regions(log10)",
  fill_cols = ECDNA_GENE_COLS,
  ymax = 3
) + ggtitle("b")

tbl_sample_wgs_summary$ecDNA_status = ifelse(tbl_sample_wgs_summary$ecDNA_status=='undetected', 'ecDNA-', tbl_sample_wgs_summary$ecDNA_status)
tbl_sample_wgs_summary$ecDNA_status = factor(tbl_sample_wgs_summary$ecDNA_status, levels = c("ecDNA-", "ecDNA+"))
plot_supp_sample_max_segment_cn <- plot_box_compare(
  df = tbl_sample_wgs_summary %>% filter(!is.na(max_segment_cn_log2)),
  xvar = "ecDNA_status",
  yvar = "max_segment_cn_log10",
  ylab = "Per-sample maximum segment copy number(log10)",
  fill_cols = ECDNA_SAMPLE_COLS,
  label_y = 3
) + ggtitle("c")

plot_supp_sample_max_oncogene_cn <- plot_box_compare(
  df = tbl_sample_wgs_summary %>% filter(!is.na(max_oncogene_cn_log2)),
  xvar = "ecDNA_status",
  yvar = "max_oncogene_cn_log10",
  ylab = "Per-sample maximum oncogene copy number(log10)",
  fill_cols = ECDNA_SAMPLE_COLS,
  label_y = 3
) + ggtitle("d")

plot_supp_amplicon_oncogene_or <- ggplot(
  tbl_amplicon_fisher,
  aes(x = odds_ratio, y = comparison)
) +
  geom_point(size = 3, color = "#D62728") +
  geom_errorbarh(aes(xmin = conf_low, xmax = conf_high), height = 0.15) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey50") +
  scale_x_log10() +
  labs(
    x = "Odds ratio (log scale)",
    y = NULL
  ) +
  theme_fig2

save(tbl_mut_sample,
     tbl_fga_plot,
     tbl_feature_plot,
     tbl_gene_cn_onco,
     tbl_prop_amplicon_onco,
     tbl_ec_sample,
     tbl_feature_plot_supp,
     tbl_gene_cn,
     tbl_sample_wgs_summary, file = file.path(DIR_RDS, "data_for_draw_plot.rds"))

# ============================================================
# 8) Assemble and save main / supplementary figures
# ============================================================
fig2_main <- (
  plot_main_func_mut |
    plot_main_oncogene_mut |
    plot_main_gii
) / (
  plot_main_segment_cn_2group |
    plot_main_oncogene_cn |
    plot_main_amplicon_oncogene
)

fig2_supp <- (
  plot_supp_gene_cn_all +
    plot_supp_sample_max_segment_cn +
    plot_supp_sample_max_oncogene_cn
)

ggsave(file.path(DIR_PLOT, "Figure2_main.pdf"), fig2_main, width = 13, height = 9)
ggsave(file.path(DIR_PLOT, "Figure2_main.png"), fig2_main, width = 13, height = 9, dpi = 300)
saveRDS(fig2_main, file.path(DIR_RDS, "Figure2_main.rds"))

ggsave(file.path(DIR_PLOT, "Figure2_supp.pdf"), fig2_supp, width = 13, height = 5)
ggsave(file.path(DIR_PLOT, "Figure2_supp.png"), fig2_supp, width = 13, height = 5, dpi = 300)
saveRDS(fig2_supp, file.path(DIR_RDS, "Figure2_supp.rds"))

# ============================================================
# 9) Local versus public mutation burden comparison
# ============================================================

tbl_ec_sample <- tbl_ec_sample %>%
  select(Sample, Pri_Met, Cancer_type)

tbl_mut_sample <- tbl_mut_sample %>%
  select(Sample, mutseq, n_func_mut, TMB_func_mut) %>%
  left_join(tbl_ec_sample, by = "Sample") %>%
  filter(Pri_Met == "Met") %>%
  select(-Pri_Met)

pub <- readr::read_tsv(PATH_PUBLIC_MUTATION, show_col_types = FALSE)

pub_sample_col <- "Sample ID"
pub_cancer_col <- "Cancer Type"
pub_mut_col <- "Mutation Count"

# Define cancer type name mapping.
normalize_cancer_type <- function(x) {
  case_when(
    is.na(x) ~ NA_character_,
    x %in% c("Lung Cancer", "Non-Small Cell Lung Cancer") ~ "Lung Cancer",
    x %in% c("Esophageal Cancer", "Gastric Cancer", "Esophagogastric Cancer") ~ "Esophagogastric cancer",
    x %in% c("Renal Cell Carcinoma") ~ "Kidney cancer",
    x %in% c("Breast Cancer") ~ "Breast cancer",
    x %in% c("Colorectal Cancer") ~ "Colorectal cancer",
    x %in% c("Cervical Cancer") ~ "Cervical cancer",
    x %in% c("Melanoma") ~ "Melanoma",
    x %in% c("Ovarian Cancer") ~ "Ovarian cancer",
    TRUE ~ x
  )
}

local_df <- tbl_mut_sample %>%
  transmute(
    Sample = as.character(Sample),
    source = "local",
    mutseq = as.character(mutseq),
    n_func_mut = as.numeric(n_func_mut),
    TMB_func_mut = as.numeric(TMB_func_mut),
    Cancer_type_raw = as.character(Cancer_type),
    cancer_type = normalize_cancer_type(Cancer_type)
  ) %>%
  filter(!is.na(Sample), !is.na(cancer_type), !is.na(n_func_mut))

public_df <- pub %>%
  transmute(
    Sample = as.character(.data[[pub_sample_col]]),
    source = "public",
    mutseq = NA_character_,
    n_func_mut = suppressWarnings(as.numeric(.data[[pub_mut_col]])),
    TMB_func_mut = `TMB (nonsynonymous)`,
    Cancer_type_raw = as.character(.data[[pub_cancer_col]]),
    cancer_type = normalize_cancer_type(.data[[pub_cancer_col]])
  ) %>%
  filter(!is.na(Sample), !is.na(cancer_type), !is.na(n_func_mut))

# Preserve both local and public shared cancer types.
common_cancers <- intersect(
  unique(local_df$cancer_type),
  unique(public_df$cancer_type)
)

all_df0 <- bind_rows(local_df, public_df) %>%
  filter(cancer_type %in% common_cancers)

min_n_each_group <- 7

cancer_n_table <- all_df0 %>%
  count(cancer_type, source) %>%
  pivot_wider(names_from = source, values_from = n, values_fill = 0)

keep_cancers <- cancer_n_table %>%
  filter(local >= min_n_each_group, public >= min_n_each_group) %>%
  pull(cancer_type)

all_df <- all_df0 %>%
  filter(cancer_type %in% keep_cancers)

cancer_order <- all_df %>%
  group_by(cancer_type) %>%
  summarise(overall_median = median(n_func_mut, na.rm = TRUE), .groups = "drop") %>%
  arrange(overall_median) %>%
  pull(cancer_type)

all_df <- all_df %>%
  mutate(cancer_type = factor(cancer_type, levels = cancer_order))

plot_df <- all_df %>%
  group_by(cancer_type, source) %>%
  arrange(n_func_mut, .by_group = TRUE) %>%
  mutate(rank_in_group = row_number()) %>%
  ungroup() %>%
  mutate(
    x_plot = if_else(source == "local", rank_in_group - 0.12, rank_in_group + 0.12)
  )

# ============================================================
# 10) Summary statistics of cancer types.
# ============================================================
cancer_summary <- all_df %>%
  group_by(cancer_type, source) %>%
  summarise(
    n = n(),
    mean_mut = mean(n_func_mut, na.rm = TRUE),
    median_mut = median(n_func_mut, na.rm = TRUE),
    sd_mut = sd(n_func_mut, na.rm = TRUE),
    iqr_mut = IQR(n_func_mut, na.rm = TRUE),
    q1_mut = quantile(n_func_mut, 0.25, na.rm = TRUE),
    q3_mut = quantile(n_func_mut, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_tsv(cancer_summary, file.path(DIR_TABLE, "public_local_cancer_summary.tsv"))

cancer_compare <- cancer_summary %>%
  select(cancer_type, source, n, mean_mut, median_mut, q1_mut, q3_mut) %>%
  pivot_wider(
    names_from = source,
    values_from = c(n, mean_mut, median_mut, q1_mut, q3_mut),
    names_sep = "_"
  ) %>%
  drop_na()

readr::write_tsv(cancer_compare, file.path(DIR_TABLE, "public_local_cancer_compare_wide.tsv"))

# =========================================================
# 11) Correlation analysis
# =========================================================
cancer_compare.rmsp <- cancer_compare %>%
  filter(cancer_type %in% cancer_order)

cor_median <- cor.test(
  cancer_compare.rmsp$median_mut_public,
  cancer_compare.rmsp$median_mut_local,
  method = "spearman",
  exact = FALSE,
  alternative = "greater"
)

p_public_local_median <- ggplot(
  cancer_compare.rmsp,
  aes(x = median_mut_public, y = median_mut_local, label = cancer_type)
) +
  geom_point(size = 2.8, alpha = 0.85) +
  geom_smooth(method = "glm", se = FALSE, color = "#e63946", fill = "#f1aeb5", linewidth = 1.2) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = Inf) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  labs(
    x = "Public median functional mutation count",
    y = "Local median functional mutation count",
    title = "Cancer-level comparison of median mutation burden",
    subtitle = paste0(
      "Spearman rho = ", round(unname(cor_median$estimate), 3),
      ", p = ", signif(cor_median$p.value, 3)
    )
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(DIR_PUB_LOCAL, "01_cancer_level_median_correlation.pdf"),
  p_public_local_median,
  width = 6.2,
  height = 6
)

# =========================================================
# 12) Differences in cancer type distribution (Wilcoxon)
# =========================================================
per_cancer_test <- all_df %>%
  group_by(cancer_type) %>%
  group_modify(~{
    x <- .x %>% filter(source == "local") %>% pull(n_func_mut)
    y <- .x %>% filter(source == "public") %>% pull(n_func_mut)
    
    tibble(
      n_local = length(x),
      n_public = length(y),
      median_local = median(x, na.rm = TRUE),
      median_public = median(y, na.rm = TRUE),
      mean_local = mean(x, na.rm = TRUE),
      mean_public = mean(y, na.rm = TRUE),
      wilcox_p = tryCatch(
        wilcox.test(x, y)$p.value,
        error = function(e) NA_real_
      )
    )
  }) %>%
  ungroup() %>%
  mutate(wilcox_fdr = p.adjust(wilcox_p, method = "BH")) %>%
  arrange(wilcox_fdr)

plot_tcga <- all_df %>%
  filter(!is.na(cancer_type), !is.na(n_func_mut), !is.na(source)) %>%
  mutate(y = log10(n_func_mut + 1)) %>%
  group_by(cancer_type, source) %>%
  arrange(y, .by_group = TRUE) %>%
  mutate(
    rank_in_group = row_number(),
    n_in_group = n(),
    frac = if_else(n_in_group == 1, 0.5, (rank_in_group - 1) / (n_in_group - 1))
  ) %>%
  ungroup()

# Cancer type order: from highest to lowest public median.
cancer_order <- plot_tcga %>%
  filter(source == "public") %>%
  group_by(cancer_type) %>%
  summarise(pub_median = median(y, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(pub_median)) %>%
  pull(cancer_type)

cancer_order <- c(cancer_order, setdiff(unique(plot_tcga$cancer_type), cancer_order))

plot_tcga <- plot_tcga %>%
  mutate(
    cancer_type = factor(cancer_type, levels = cancer_order),
    cancer_id = as.numeric(cancer_type)
  )

spread_public <- 0.34
spread_local <- 0.22

plot_tcga <- plot_tcga %>%
  mutate(
    x = case_when(
      source == "public" ~ cancer_id - spread_public + frac * spread_public,
      source == "local" ~ cancer_id - spread_local + frac * spread_local,
      TRUE ~ as.numeric(cancer_id)
    )
  )

median_df <- plot_tcga %>%
  group_by(cancer_type, cancer_id, source) %>%
  summarise(
    med_y = median(y, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    xstart = case_when(
      source == "public" ~ cancer_id - 0.46,
      source == "local" ~ cancer_id - 0.18,
      TRUE ~ cancer_id - 0.2
    ),
    xend = case_when(
      source == "public" ~ cancer_id + 0.46,
      source == "local" ~ cancer_id + 0.18,
      TRUE ~ cancer_id + 0.2
    )
  )

median_pub <- median_df %>%
  filter(source == "public")

median_local <- median_df %>%
  filter(source == "local")

# Significance tests were performed on each cancer type.
stat_df <- all_df %>%
  filter(cancer_type %in% cancer_order) %>%
  mutate(cancer_type = factor(cancer_type, levels = cancer_order)) %>%
  group_by(cancer_type) %>%
  group_modify(~{
    x <- .x %>% filter(source == "local") %>% pull(TMB_func_mut)
    y <- .x %>% filter(source == "public") %>% pull(TMB_func_mut)

    p_val <- tryCatch(
      wilcox.test(x, y)$p.value,
      error = function(e) NA_real_
    )

    tibble(
      n_local = length(x),
      n_public = length(y),
      p_value = p_val
    )
  }) %>%
  ungroup() %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    signif_label = case_when(
      is.na(p_adj) ~ "NA",
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

y_pos_df <- plot_tcga %>%
  group_by(cancer_type, cancer_id) %>%
  summarise(
    y_max = max(y, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(y_sig = y_max + 0.12)

stat_df <- stat_df %>%
  left_join(y_pos_df, by = "cancer_type")

bg_df <- tibble(
  cancer_id = seq_along(cancer_order),
  xmin = cancer_id - 0.5,
  xmax = cancer_id + 0.5,
  fill_group = ifelse(cancer_id %% 2 == 0, "even", "odd")
)

p_public_local_distribution <- ggplot() +
  geom_rect(
    data = bg_df,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill_group),
    inherit.aes = FALSE,
    alpha = 0.22
  ) +
  scale_fill_manual(values = c("odd" = "#f3f6df", "even" = "#dbe8f4"), guide = "none") +
  geom_path(
    data = plot_tcga %>% filter(source == "public"),
    aes(x = x, y = y, group = cancer_type),
    color = "grey70",
    linewidth = 1.0,
    lineend = "round"
  ) +
  geom_point(
    data = plot_tcga %>% filter(source == "public"),
    aes(x = x, y = y),
    color = "grey65",
    size = 0.9,
    alpha = 0.7
  ) +
  geom_path(
    data = plot_tcga %>% filter(source == "local"),
    aes(x = x, y = y, group = cancer_type),
    color = "blue",
    linewidth = 1.0,
    lineend = "round"
  ) +
  geom_point(
    data = plot_tcga %>% filter(source == "local"),
    aes(x = x, y = y),
    color = "blue",
    size = 1.0,
    alpha = 0.9
  ) +
  geom_segment(
    data = median_pub,
    aes(x = xstart, xend = xend, y = med_y, yend = med_y),
    color = "grey70",
    linewidth = 0.9
  ) +
  geom_segment(
    data = median_local,
    aes(x = xstart, xend = xend, y = med_y, yend = med_y),
    color = "blue",
    linewidth = 0.9
  ) +
  geom_text(
    data = stat_df,
    aes(x = cancer_id, y = y_sig, label = signif_label),
    size = 4,
    fontface = "bold"
  ) +
  scale_x_continuous(
    breaks = seq_along(cancer_order),
    labels = cancer_order,
    expand = c(0.01, 0.01)
  ) +
  scale_y_continuous(
    breaks = 0:6,
    labels = 0:6,
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(
    x = NULL,
    y = "log10(Functional variants + 1)",
    title = "Local vs public mutation-count distribution by cancer type"
  ) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "grey40"),
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1, color = "black"),
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

ggsave(
  file.path(DIR_PUB_LOCAL, "02_public_local_distribution_by_cancer.pdf"),
  p_public_local_distribution,
  width = 8,
  height = 5
)
