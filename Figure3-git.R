# =========================================================
# 0) Packages
# =========================================================
library(tidyverse)
library(readxl)
library(patchwork)
library(scales)
library(ggalluvial)
library(ggnewscale)

# =========================================================
# 1) Input configuration and output directories
# =========================================================
input_dir <- "Inputs"
output_dir <- "Outputs"

outdir <- file.path(output_dir, "Figure3_Paired_ecDNA")
outdir_all_paired <- file.path(outdir, "AllPaired")

purrr::walk(c(outdir, outdir_all_paired), ~{
  dir.create(.x, showWarnings = FALSE, recursive = TRUE)
})

meta_xlsx <- file.path(input_dir, "metadata", "sample_metadata.xlsx")
meta_sheet <- "paired-ecDNA"

wes_records_csv <- file.path(input_dir, "GCAP", "batch_wes_run_fCNA_records_with_symbols.csv")
oncogene_list_txt <- file.path(input_dir, "database", "combined_oncogene_list.txt")
ecdna_gene_rds <- file.path(input_dir, "database", "ecDNA_genelist_all.rds")

# =========================================================
# 2) Colors
# =========================================================
cancer_cols <- c(
  "Lung cancer"          = "#A38A77",
  "Breast cancer"        = "#FFB6C1",
  "Gastric cancer"       = "#984EA3",
  "Renal cell carcinoma" = "#FF7F00",
  "Kidney cancer"        = "#FF7F00",
  "Melanoma"             = "#A65628",
  "Colorectal cancer"    = "#009FE8",
  "Cervical cancer"      = "#542788",
  "Ovarian cancer"       = "#C51B7D",
  "Esophageal cancer"    = "#008B8B"
)

dynamic_cols <- c(
  "None -> None" = "grey70",
  "Gain in Met" = "#D62728",
  "Loss in Met" = "#1F77B4",
  "Conserved"   = "#2CA02C"
)

sharing_cols <- c(
  "Shared" = "#2CA02C",
  "Primary-only" = "#1F77B4",
  "Metastasis-only" = "#D62728"
)

stage_cols <- c(
  "Pri" = "#4C78A8",
  "Met" = "#E45756"
)

# =========================================================
# 3) Helper functions
# =========================================================
norm_yes <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x %in% c("yes", "y", "true", "1")
}

parse_genes <- function(x) {
  if (length(x) == 0) return(character(0))
  if (is.na(x) || x == "/" || x == "") return(character(0))

  x <- as.character(x)
  x <- gsub("\\[|\\]|'", "", x)
  x <- gsub("\"", "", x)
  x <- trimws(x)
  if (x == "/" || x == "") return(character(0))

  out <- unlist(strsplit(x, ",\\s*"))
  out <- trimws(out)
  out <- out[out != ""]
  unique(out)
}

union_from_cols <- function(df, cols) {
  if (length(cols) == 0) {
    return(replicate(nrow(df), character(0), simplify = FALSE))
  }

  apply(df[, cols, drop = FALSE], 1, function(row) {
    genes <- character(0)
    for (v in row) genes <- union(genes, parse_genes(v))
    genes
  }) %>%
    as.list()
}

get_col <- function(df, candidates) {
  nms <- colnames(df)
  low <- tolower(nms)

  for (c in candidates) {
    idx <- which(low == tolower(c))
    if (length(idx) == 1) return(nms[idx])
  }

  NULL
}

empty_chr_list <- function(n) {
  replicate(n, character(0), simplify = FALSE)
}

normalize_listcol <- function(x) {
  purrr::map(x, function(.x) {
    if (is.null(.x)) return(character(0))
    if (length(.x) == 1 && is.na(.x)) return(character(0))

    .x <- as.character(.x)
    .x <- .x[!is.na(.x)]
    .x <- trimws(.x)
    .x[.x != ""]
  })
}

union2_safe <- function(x, y) {
  union(
    if (is.null(x)) character(0) else x,
    if (is.null(y)) character(0) else y
  )
}

as_logical_flag <- function(x) {
  x %in% c(TRUE, "TRUE", "True", "T", "1", 1)
}

flatten_list_columns <- function(df) {
  df %>%
    mutate(across(where(is.list), ~ purrr::map_chr(.x, function(z) paste(z, collapse = ", "))))
}

# =========================================================
# 4) Read paired metadata
# =========================================================
df0 <- readxl::read_excel(meta_xlsx, sheet = meta_sheet) %>%
  as.data.frame()

col_sample <- get_col(df0, c("Sample", "sample", "sample_name"))
col_pair <- get_col(df0, c("Pair.", "Pair", "pair"))
col_stage <- get_col(df0, c("Pri.Met", "Pri_Met", "PriMet", "pri_met"))
col_cancer <- get_col(df0, c("Cancer.type", "Cancer_type", "CancerType"))
col_seq <- get_col(df0, c("Seq_type", "Seq.type", "SeqType", "seq_type"))
col_used <- get_col(df0, c("used", "Used", "USE"))
col_ec_cnt <- get_col(df0, c("eCDNA_count", "ecDNA_count", "ecdna_count"))

if (any(sapply(list(col_sample, col_pair, col_stage, col_cancer, col_seq, col_used), is.null))) {
  stop("The paired metadata sheet must contain Sample, Pair, Pri/Met, Cancer type, Seq type, and used columns.")
}

onco_cols <- colnames(df0)[grepl("^ecDNA(_|\\.)?oncogene_", colnames(df0), ignore.case = TRUE)]
allgene_cols <- colnames(df0)[grepl("^All(_|\\.)?genes_", colnames(df0), ignore.case = TRUE)]

df <- df0 %>%
  mutate(.row_id = row_number()) %>%
  transmute(
    .row_id,
    Sample = trimws(.data[[col_sample]]),
    Pair = trimws(.data[[col_pair]]),
    Stage_raw = .data[[col_stage]],
    Cancer.type = trimws(.data[[col_cancer]]),
    Seq.type = toupper(trimws(.data[[col_seq]])),
    used = norm_yes(.data[[col_used]]),
    eCDNA_count = if (!is.null(col_ec_cnt)) suppressWarnings(as.numeric(.data[[col_ec_cnt]])) else NA_real_
  ) %>%
  mutate(
    Stage = case_when(
      grepl("^pri", Stage_raw, ignore.case = TRUE) ~ "Pri",
      grepl("^met", Stage_raw, ignore.case = TRUE) ~ "Met",
      TRUE ~ NA_character_
    ),
    Cancer.type = str_trim(Cancer.type),
    Cancer.type = sub("^(.)", "\\U\\1", Cancer.type, perl = TRUE)
  ) %>%
  filter(!is.na(Pair), Pair != "", !is.na(Stage))

if (length(onco_cols) > 0) {
  df[onco_cols] <- df0[df$.row_id, onco_cols, drop = FALSE]
}
if (length(allgene_cols) > 0) {
  df[allgene_cols] <- df0[df$.row_id, allgene_cols, drop = FALSE]
}

df <- df %>%
  select(-.row_id)

# =========================================================
# 5) Define allowed samples by use flag and platform
# =========================================================
allowed_wgs_samples <- df %>%
  filter(used, Seq.type == "WGS") %>%
  pull(Sample) %>%
  unique()

allowed_wes_samples <- df %>%
  filter(used, Seq.type == "WES") %>%
  pull(Sample) %>%
  unique()

# =========================================================
# 6) Read oncogene list and WES records
# =========================================================
if (!file.exists(wes_records_csv)) stop("Missing WES records: ", wes_records_csv)
if (!file.exists(oncogene_list_txt)) stop("Missing oncogene list: ", oncogene_list_txt)

oncogenes <- read.table(oncogene_list_txt, stringsAsFactors = FALSE)$V1

genes_wes0 <- readr::read_csv(wes_records_csv, show_col_types = FALSE)

low <- tolower(colnames(genes_wes0))
if (!("sample" %in% low)) {
  idx <- which(low %in% c("sample_name", "samplename", "id", "sampleid"))
  if (length(idx) == 1) colnames(genes_wes0)[idx] <- "sample"
}

low <- tolower(colnames(genes_wes0))
if (!("gene_class" %in% low)) {
  idx <- which(low %in% c("geneclass", "class"))
  if (length(idx) == 1) colnames(genes_wes0)[idx] <- "gene_class"
}

low <- tolower(colnames(genes_wes0))
if (!("symbol" %in% low)) {
  idx <- which(low %in% c("gene", "genesymbol", "symbols"))
  if (length(idx) == 1) colnames(genes_wes0)[idx] <- "SYMBOL"
}

low <- tolower(colnames(genes_wes0))
if (!("total_cn" %in% low)) {
  idx <- which(low %in% c("totalcn", "cn", "copy_number", "copynumber"))
  if (length(idx) == 1) colnames(genes_wes0)[idx] <- "total_cn"
}

need <- c("sample", "gene_class", "SYMBOL", "total_cn")
if (!all(need %in% colnames(genes_wes0))) {
  stop(
    "WES records must include: ", paste(need, collapse = ", "),
    "\nCurrent: ", paste(colnames(genes_wes0), collapse = ", ")
  )
}

genes_wes <- genes_wes0 %>%
  mutate(
    sample = trimws(sample),
    gene_class = tolower(trimws(gene_class)),
    SYMBOL = trimws(SYMBOL),
    total_cn = suppressWarnings(as.numeric(total_cn))
  ) %>%
  filter(sample %in% allowed_wes_samples)

wes_circular_samples <- genes_wes %>%
  filter(gene_class == "circular") %>%
  pull(sample) %>%
  unique()

wes_ec_summary <- genes_wes %>%
  filter(gene_class == "circular", SYMBOL != "") %>%
  group_by(sample) %>%
  summarise(
    wes_genes = list(unique(SYMBOL)),
    wes_gene_n = n_distinct(SYMBOL),
    wes_oncogenes = list(intersect(unique(SYMBOL), oncogenes)),
    wes_onco_n = length(intersect(unique(SYMBOL), oncogenes)),
    .groups = "drop"
  )

# =========================================================
# 7) Add WGS and WES evidence
# =========================================================
df <- df %>%
  mutate(
    wgs_valid = Sample %in% allowed_wgs_samples,
    wgs_ecdna_pos = wgs_valid & !is.na(eCDNA_count) & eCDNA_count > 0
  )

if (length(onco_cols) > 0) {
  tmp_onco <- union_from_cols(df, onco_cols)
  df$wgs_oncogenes <- purrr::map2(df$wgs_valid, tmp_onco, ~ if (.x) .y else character(0))
  df$wgs_onco_n <- purrr::map_int(df$wgs_oncogenes, length)
} else {
  df$wgs_oncogenes <- empty_chr_list(nrow(df))
  df$wgs_onco_n <- 0
}

if (length(allgene_cols) > 0) {
  tmp_all <- union_from_cols(df, allgene_cols)
  df$wgs_allgenes <- purrr::map2(df$wgs_valid, tmp_all, ~ if (.x) .y else character(0))
  df$wgs_gene_n <- purrr::map_int(df$wgs_allgenes, length)
} else {
  df$wgs_allgenes <- empty_chr_list(nrow(df))
  df$wgs_gene_n <- 0
}

df <- df %>%
  left_join(wes_ec_summary, by = c("Sample" = "sample")) %>%
  mutate(
    wes_valid = Sample %in% allowed_wes_samples,
    wes_ecdna_pos = wes_valid & (Sample %in% wes_circular_samples),
    wes_gene_n = ifelse(wes_valid & !is.na(wes_gene_n), wes_gene_n, 0),
    wes_onco_n = ifelse(wes_valid & !is.na(wes_onco_n), wes_onco_n, 0)
  )

df$wes_genes <- normalize_listcol(df$wes_genes)
df$wes_oncogenes <- normalize_listcol(df$wes_oncogenes)

df <- df %>%
  mutate(
    wes_genes = purrr::map2(wes_valid, wes_genes, ~ if (.x) .y else character(0)),
    wes_oncogenes = purrr::map2(wes_valid, wes_oncogenes, ~ if (.x) .y else character(0)),
    ecDNA_pos = wgs_ecdna_pos | wes_ecdna_pos,
    onco_pos = (wgs_onco_n > 0) | (wes_onco_n > 0)
  )

# =========================================================
# 8) Resolve duplicates: one row per pair and stage
# =========================================================
df_resolved <- df %>%
  mutate(
    seq_rank = case_when(
      Seq.type == "WGS" ~ 2,
      Seq.type == "WES" ~ 1,
      TRUE ~ 0
    ),
    evidence_rank = (as.integer(ecDNA_pos) * 1000) +
      (wgs_onco_n + wes_onco_n) * 10 +
      (wgs_gene_n + wes_gene_n)
  ) %>%
  group_by(Pair, Stage) %>%
  arrange(desc(used), desc(seq_rank), desc(evidence_rank), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

pair_has_both <- df_resolved %>%
  distinct(Pair, Stage) %>%
  count(Pair) %>%
  filter(n == 2) %>%
  pull(Pair)

df_pair <- df_resolved %>%
  filter(Pair %in% pair_has_both)

write.csv(flatten_list_columns(df_resolved), file.path(outdir, "paired_resolved_table.csv"), row.names = FALSE)
write.csv(flatten_list_columns(df_pair), file.path(outdir, "paired_complete_pairs.csv"), row.names = FALSE)

# =========================================================
# 9) Panel A: paired cohort by cancer type
# =========================================================
df_A <- df_pair %>%
  distinct(Pair, Cancer.type) %>%
  mutate(Cancer.type = as.character(Cancer.type))

pA <- ggplot(df_A, aes(x = forcats::fct_infreq(Cancer.type))) +
  geom_bar(aes(fill = forcats::fct_infreq(Cancer.type))) +
  scale_fill_manual(values = cancer_cols) +
  coord_flip() +
  theme_classic(base_size = 14) +
  labs(
    x = NULL,
    y = "Number of paired cases",
    title = "Paired cohort by cancer type",
    fill = "Cancer type"
  ) +
  theme(legend.position = "bottom")

ggsave(file.path(outdir, "FigureA_pairs_by_cancer.pdf"), pA, width = 6, height = 4)

# =========================================================
# 10) Panel B: ecDNA presence dynamics
# =========================================================
df_B <- df_pair %>%
  select(Pair, Stage, ecDNA_pos, Cancer.type) %>%
  mutate(ecDNA_pos = as.integer(ecDNA_pos)) %>%
  pivot_wider(names_from = Stage, values_from = ecDNA_pos) %>%
  mutate(
    category = case_when(
      Pri == 0 & Met == 0 ~ "None -> None",
      Pri == 0 & Met == 1 ~ "Gain in Met",
      Pri == 1 & Met == 0 ~ "Loss in Met",
      Pri == 1 & Met == 1 ~ "Conserved"
    ),
    category = factor(category, levels = c("None -> None", "Gain in Met", "Loss in Met", "Conserved"))
  )

df_B_count <- df_B %>%
  count(category, name = "n") %>%
  mutate(prop = n / sum(n))

pB_count <- ggplot(df_B_count, aes(x = category, y = n)) +
  geom_col(fill = "#B2182B", width = 0.8) +
  theme_classic(base_size = 14) +
  labs(x = NULL, y = "Number of pairs", title = "Dynamics of ecDNA presence (Pri -> Met)")

pB_percent <- ggplot(df_B_count, aes(x = category, y = prop)) +
  geom_col(fill = "#B2182B", width = 0.8) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = c(0, 0)) +
  theme_classic(base_size = 14) +
  labs(x = NULL, y = "Percentage of pairs", title = "Dynamics of ecDNA presence (Pri -> Met)")

df_B_long <- df_B %>%
  pivot_longer(cols = c(Pri, Met), names_to = "Stage", values_to = "ecDNA_pos") %>%
  mutate(Stage = factor(Stage, levels = c("Pri", "Met")))

pB_line <- ggplot(df_B_long, aes(x = Stage, y = ecDNA_pos, group = Pair, color = category)) +
  geom_line(alpha = 0.6, linewidth = 0.6) +
  geom_point(size = 2) +
  scale_y_continuous(breaks = c(0, 1), labels = c("ecDNA-", "ecDNA+")) +
  scale_color_manual(values = dynamic_cols) +
  theme_classic(base_size = 14) +
  labs(x = NULL, y = "ecDNA status", color = NULL, title = "Paired ecDNA status trajectories")

ggsave(file.path(outdir, "FigureB_ecDNA_dynamics_count.pdf"), pB_count, width = 5, height = 4)
ggsave(file.path(outdir, "FigureB_ecDNA_dynamics_percent.pdf"), pB_percent, width = 6, height = 4)
ggsave(file.path(outdir, "FigureB_ecDNA_dynamics_lines.pdf"), pB_line, width = 6, height = 4)

# =========================================================
# 11) Panel B supplementary: Sankey plot
# =========================================================
df_B2 <- df_B %>%
  select(Pair, Cancer.type, category)

df_sankey <- df_B2 %>%
  mutate(
    Cancer.type = str_trim(Cancer.type),
    Cancer.type = sub("^(.)", "\\U\\1", Cancer.type, perl = TRUE),
    Pri_cancer = Cancer.type,
    Met_cancer = Cancer.type
  )

cancer_levels <- df_sankey %>%
  count(Cancer.type, name = "n") %>%
  arrange(desc(n)) %>%
  pull(Cancer.type)

df_sankey <- df_sankey %>%
  mutate(
    Pri_cancer = factor(Pri_cancer, levels = cancer_levels),
    Met_cancer = factor(Met_cancer, levels = cancer_levels)
  )

plot_df <- df_sankey %>%
  count(Pri_cancer, category, Met_cancer, name = "Freq")

lab_cancer <- df_sankey %>%
  count(Pri_cancer, name = "n") %>%
  mutate(
    pct = n / sum(n),
    label = paste0(as.character(Pri_cancer), "\n", n, " (", percent(pct, accuracy = 0.1), ")")
  ) %>%
  select(name = Pri_cancer, label)

lab_category <- df_sankey %>%
  count(category, name = "n") %>%
  mutate(
    pct = n / sum(n),
    label = paste0(as.character(category), "\n", n, " (", percent(pct, accuracy = 0.1), ")")
  ) %>%
  select(name = category, label)

label_map <- bind_rows(lab_cancer, lab_category)

pB_sankey <- ggplot(
  plot_df,
  aes(axis1 = Pri_cancer, axis2 = category, axis3 = Met_cancer, y = Freq)
) +
  ggalluvial::geom_alluvium(aes(fill = category), width = 1 / 3, alpha = 0.85, color = "white") +
  scale_fill_manual(values = dynamic_cols, drop = FALSE) +
  ggnewscale::new_scale_fill() +
  ggalluvial::geom_stratum(aes(fill = after_stat(stratum)), width = 1 / 2, color = "grey40") +
  scale_fill_manual(values = c(cancer_cols, dynamic_cols), breaks = c(names(dynamic_cols), names(cancer_cols))) +
  geom_text(
    stat = "stratum",
    aes(label = label_map$label[match(after_stat(stratum), label_map$name)]),
    size = 2.5,
    lineheight = 0.95
  ) +
  scale_x_discrete(
    limits = c("Primary", "ecDNA dynamic", "Metastasis"),
    expand = c(0.06, 0.06)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

ggsave(file.path(outdir, "FigureB_ecDNA_sankey.pdf"), pB_sankey, width = 9, height = 4.5)

# =========================================================
# 12) Panel C: shared, primary-only, and metastasis-only oncogenes
# =========================================================
df_onco_long <- df_pair %>%
  transmute(
    Pair,
    Stage,
    oncos = purrr::map2(wgs_oncogenes, wes_oncogenes, union2_safe)
  ) %>%
  tidyr::unnest(oncos) %>%
  filter(!is.na(oncos), oncos != "") %>%
  distinct(Pair, Stage, oncos)

df_C <- df_onco_long %>%
  group_by(Pair, oncos) %>%
  summarise(
    Pri = any(Stage == "Pri"),
    Met = any(Stage == "Met"),
    .groups = "drop"
  ) %>%
  mutate(
    status = case_when(
      Pri & Met ~ "Shared",
      Pri & !Met ~ "Primary-only",
      !Pri & Met ~ "Metastasis-only"
    ),
    status = factor(status, levels = c("Shared", "Primary-only", "Metastasis-only"))
  )

pC <- ggplot(df_C, aes(x = status, fill = status)) +
  geom_bar() +
  scale_fill_manual(values = sharing_cols) +
  theme_classic(base_size = 14) +
  labs(x = NULL, y = "Number of ecDNA oncogenes", title = "Oncogene sharing across paired samples")

ggsave(file.path(outdir, "FigureC_oncogene_shared_bar.pdf"), pC, width = 6, height = 4)

df_C_pairdelta <- df_C %>%
  count(Pair, status, name = "n") %>%
  tidyr::pivot_wider(names_from = status, values_from = n, values_fill = 0)

write.csv(df_C_pairdelta, file.path(outdir, "oncogene_pairwise_directionality.csv"), row.names = FALSE)

# =========================================================
# 13) Panel D: gene burden in primary and metastatic samples
# =========================================================
df_D <- df_pair %>%
  transmute(
    Pair,
    Stage = factor(Stage, levels = c("Pri", "Met")),
    gene_set = purrr::map2(wgs_allgenes, wes_genes, union2_safe),
    gene_n = purrr::map_int(gene_set, length)
  )

df_D_wide <- df_D %>%
  select(Pair, Stage, gene_n) %>%
  pivot_wider(names_from = Stage, values_from = gene_n) %>%
  filter(!is.na(Pri) & !is.na(Met)) %>%
  filter(Pri > 0 | Met > 0)

test_D <- wilcox.test(df_D_wide$Met, df_D_wide$Pri, alternative = "greater", paired = TRUE)
p_label_D <- paste0("paired Wilcoxon p = ", signif(test_D$p.value, 3))
y_max_D <- max(df_D$gene_n, na.rm = TRUE)

pD <- ggplot(df_D, aes(x = Stage, y = gene_n, group = Pair)) +
  geom_line(alpha = 0.3, linewidth = 0.6) +
  geom_point(size = 2) +
  annotate(
    "text",
    x = 1.15,
    y = y_max_D + max(3, y_max_D * 0.06),
    label = p_label_D,
    hjust = 0,
    vjust = 1,
    size = 4
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = NULL,
    y = "Number of genes on ecDNA",
    title = "ecDNA gene burden (Pri vs Met)"
  )

ggsave(file.path(outdir, "FigureD_gene_burden_lines.pdf"), pD, width = 4.5, height = 4)

# =========================================================
# 14) Panel E: top oncogenes
# =========================================================
topN <- 15

df_E <- df_C %>%
  count(oncos, status, name = "n") %>%
  group_by(oncos) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = topN) %>%
  select(oncos) %>%
  left_join(df_C %>% count(oncos, status, name = "n"), by = "oncos") %>%
  mutate(
    oncos = factor(oncos, levels = rev(unique(oncos))),
    status = factor(status, levels = c("Shared", "Primary-only", "Metastasis-only"))
  )

pE <- ggplot(df_E, aes(x = oncos, y = n, fill = status)) +
  geom_col(width = 0.75, color = "black", linewidth = 0.2) +
  coord_flip() +
  scale_fill_manual(values = sharing_cols) +
  theme_classic(base_size = 14) +
  labs(
    x = NULL,
    y = "Count (pair x oncogene)",
    fill = NULL,
    title = paste0("Top ", topN, " ecDNA oncogenes (shared/gain/loss)")
  )

ggsave(file.path(outdir, "FigureE_top_oncogenes.pdf"), pE, width = 7, height = 5)

# =========================================================
# 15) Supplementary oncogene-level copy-number table
# =========================================================
sample_meta <- df %>%
  select(Sample, Pair, Stage, Cancer.type, Seq.type, used) %>%
  distinct()

cn_long0 <- genes_wes %>%
  filter(!is.na(SYMBOL), SYMBOL != "", !is.na(total_cn)) %>%
  mutate(
    Sample = sample,
    SYMBOL = trimws(SYMBOL),
    is_oncogene = SYMBOL %in% oncogenes,
    cna_group = case_when(
      gene_class == "circular" ~ "ecDNA",
      TRUE ~ "Other focal CNA"
    )
  ) %>%
  left_join(sample_meta, by = "Sample")

cn_onco <- cn_long0 %>%
  filter(is_oncogene)

write.csv(cn_long0, file.path(outdir, "supp_gene_level_cn_all.csv"), row.names = FALSE)
write.csv(cn_onco, file.path(outdir, "supp_gene_level_cn_oncogene.csv"), row.names = FALSE)

# =========================================================
# 16) Assemble paired ecDNA figures
# =========================================================
p_main <- (pA | pB_line) / (pC | pD) / pE +
  plot_annotation(title = "Paired ecDNA analysis")

p_supp <- (pB_count | pB_percent) / pB_sankey +
  plot_annotation(title = "Supplementary paired ecDNA dynamics")

ggsave(file.path(outdir, "Paired_ecDNA_main.pdf"), p_main, width = 12, height = 12)
ggsave(file.path(outdir, "Paired_ecDNA_supp.pdf"), p_supp, width = 12, height = 10)

# =========================================================
# 17) Paired gene copy-number comparison across all paired samples
# =========================================================
dat <- readRDS(ecdna_gene_rds)

prep_paired_gene_data <- function(dat_input) {
  df0 <- dat_input %>%
    mutate(is_ecDNA = as_logical_flag(is_ecDNA)) %>%
    filter(
      is_ecDNA,
      !is.na(Pair), Pair != "",
      Pri_Met %in% c("Pri", "Met")
    )

  valid_pairs <- df0 %>%
    distinct(Pair, Pri_Met, Cancer_type) %>%
    count(Pair) %>%
    filter(n == 2) %>%
    pull(Pair)

  df0 <- df0 %>%
    filter(Pair %in% valid_pairs)

  df_gene <- df0 %>%
    group_by(Pair, Pri_Met, Cancer_type, gene) %>%
    summarise(gene_cn = max(gene_cn, na.rm = TRUE), .groups = "drop")

  df_gene %>%
    pivot_wider(names_from = Pri_Met, values_from = gene_cn) %>%
    filter(!is.na(Pri) & !is.na(Met))
}

run_overall_analysis <- function(dat_input,
                                 out_prefix,
                                 plot_title,
                                 ylab_text = "Gene copy number (CN)") {
  df_wide <- prep_paired_gene_data(dat_input)

  tt <- t.test(
    df_wide$Met,
    df_wide$Pri,
    paired = TRUE,
    alternative = "greater"
  )

  df_long <- df_wide %>%
    pivot_longer(
      cols = c(Pri, Met),
      names_to = "Pri_Met",
      values_to = "gene_cn"
    ) %>%
    mutate(Pri_Met = factor(Pri_Met, levels = c("Pri", "Met")))

  ymax <- max(df_long$gene_cn, na.rm = TRUE)
  label_y <- ymax * 1.08

  p <- ggplot(df_long, aes(x = Pri_Met, y = gene_cn)) +
    geom_boxplot(
      aes(fill = Pri_Met),
      width = 0.6,
      outlier.shape = NA,
      color = "black",
      linewidth = 0.6,
      staplewidth = 0.5
    ) +
    geom_line(
      aes(group = interaction(Pair, gene)),
      alpha = 0.22,
      linewidth = 0.5,
      color = "grey55"
    ) +
    geom_point(
      aes(fill = Pri_Met),
      shape = 21,
      size = 2.1,
      color = "black",
      stroke = 0.2
    ) +
    annotate(
      "text",
      x = 1.5,
      y = label_y,
      label = paste0("paired t-test, p = ", signif(tt$p.value, 3)),
      size = 4
    ) +
    scale_fill_manual(values = stage_cols) +
    theme_classic(base_size = 14) +
    labs(
      x = NULL,
      y = ylab_text,
      title = plot_title
    ) +
    theme(legend.position = "none") +
    coord_cartesian(ylim = c(min(df_long$gene_cn, na.rm = TRUE), label_y * 1.05))

  ggsave(file.path(outdir_all_paired, paste0(out_prefix, ".pdf")), p, width = 5, height = 4.5)

  invisible(list(
    df_wide = df_wide,
    df_long = df_long,
    ttest = tt,
    plot = p
  ))
}

res_all_overall <- run_overall_analysis(
  dat_input = dat,
  out_prefix = "Figure3_ecDNA_all_gene_CN_Pri_vs_Met_paired",
  plot_title = "Paired comparison of ecDNA gene copy number",
  ylab_text = "Gene copy number (CN)"
)

res_onco_overall <- run_overall_analysis(
  dat_input = dat %>% filter(as_logical_flag(is_canonical_oncogene)),
  out_prefix = "Figure3_ecDNA_oncogene_CN_Pri_vs_Met_paired",
  plot_title = "Paired comparison of ecDNA oncogene copy number",
  ylab_text = "Oncogene copy number (CN)"
)

ggsave(
  file.path(outdir_all_paired, "Figure3_ecDNA_gene_CN_Pri_vs_Met_paired.pdf"),
  res_all_overall$plot + res_onco_overall$plot,
  width = 7,
  height = 5
)
