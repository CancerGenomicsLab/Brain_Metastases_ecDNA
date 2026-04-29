options(stringsAsFactors = FALSE)

# =========================================================
# 0) Packages
# =========================================================
library(ggplot2)
library(dplyr)
library(readxl)
library(stringr)
library(ggnewscale)
library(readr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(scales)
library(grid)

# =========================================================
# 1) Input configuration and output directory
# =========================================================
input_dir <- "Inputs"
output_dir <- "Outputs"
outdir <- file.path(output_dir, "Figure1_Atlas")
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

meta_xlsx <- file.path(input_dir, "sample_metadata.xlsx")
meta_sheet <- "ecDNA WGS Sample Sequencing Statistics-refined"
wgs_genelist_dir <- file.path(input_dir, "WGS_genelist")
wes_records_csv <- file.path(input_dir, "GCAP", "batch_wes_run_fCNA_records_with_symbols.csv")
oncogene_list_txt <- file.path(input_dir, "database", "combined_oncogene_list.txt")

# =========================================================
# 2) Colors and plotting theme
# =========================================================
cancer_cols <- c(
  "Lung cancer"          = "#A38A77",
  "Breast cancer"        = "#FFB6C1",
  "Gastric cancer"       = "#984EA3",
  "Renal cell carcinoma" = "#FF7F00",
  "Melanoma"             = "#A65628",
  "Colorectal cancer"    = "#009FE8",
  "Cervical cancer"      = "#542788",
  "Ovarian cancer"       = "#C51B7D",
  "Esophageal cancer"    = "#008B8B"
)

status_cols <- c(
  "Amplicon (Linear)" = "#DC0000",
  "ecDNA+"            = "#8B0000",
  "Negative"          = "#E0E0E0"
)

oncogene_status_cols <- c(
  "ecDNA+ (Oncogene-)" = "#F39B7F",
  "ecDNA+ (Oncogene+)" = "#DC0000"
)

ecdna_binary_cols <- c(
  "ecDNA-" = "#E0E0E0",
  "ecDNA+" = "#D53E4F"
)

theme_nature <- theme_classic() +
  theme(
    text = element_text(family = "sans", color = "black"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10, color = "black"),
    legend.position = "right"
  )


# =========================================================
# 3) Metadata table
# =========================================================
# Denominator: samples marked as used == yes, with one row per unique sample.
df_meta_raw <- readxl::read_excel(meta_xlsx, sheet = meta_sheet) %>% as.data.frame()

df_meta <- df_meta_raw %>%
  mutate(
    used = tolower(str_trim(used)),
    Sample = str_trim(Sample),
    Cancer_type = sub("^(.)", "\\U\\1", str_trim(Cancer_type), perl = TRUE),
    Seq_type = toupper(str_trim(Seq_type))
  ) %>%
  filter(used == "yes") %>%
  distinct(Sample, .keep_all = TRUE)

wgs_samples <- df_meta[df_meta$Seq_type == "WGS", ]$Sample
wes_samples <- df_meta[df_meta$Seq_type == "WES", ]$Sample

# cancer order for barplot columns
cancer_order_global <- df_meta %>% count(Cancer_type, sort = TRUE) %>% pull(Cancer_type)

# =========================================================
# 4) WGS genelist
# =========================================================
files <- list.files(wgs_genelist_dir, pattern = "\\.tsv$", full.names = TRUE)
if(length(files) == 0) stop("WGS genelist dir has no .tsv: ", wgs_genelist_dir)

df_gene_wgs <- do.call(rbind, lapply(files, function(x){
  read.table(x, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
})) %>%
  mutate(
    sample_name = str_trim(sample_name),
    feature = as.character(feature),
    gene = as.character(gene),
    gene_cn = suppressWarnings(as.numeric(gene_cn)),
    is_canonical_oncogene = as.character(is_canonical_oncogene)
  ) %>%
  filter(!is.na(gene_cn), gene_cn > 4)

df_gene_wgs_den <- df_gene_wgs %>% filter(sample_name %in% df_meta$Sample)

wgs_ecDNA_samples <- df_gene_wgs_den %>%
  filter(grepl("ecDNA", feature, ignore.case = TRUE)) %>%
  pull(sample_name) %>% unique()

wgs_amp_samples <- df_gene_wgs_den %>% pull(sample_name) %>% unique()

# =========================================================
# 5) WES/GCAP records
# =========================================================
genes_wes0 <- readr::read_csv(wes_records_csv, show_col_types = FALSE)

# Standardize column names when equivalent alternatives are present.
cn <- colnames(genes_wes0)
cn_low <- tolower(cn)

if (!("sample" %in% cn_low)) {
  idx <- which(cn_low %in% c("sample_name", "samplename", "id", "sampleid"))
  if (length(idx) == 1) colnames(genes_wes0)[idx] <- "sample"
}
if (!("gene_class" %in% cn_low)) {
  idx <- which(cn_low %in% c("geneclass", "class"))
  if (length(idx) == 1) colnames(genes_wes0)[idx] <- "gene_class"
}
if (!("symbol" %in% cn_low)) {
  idx <- which(cn_low %in% c("gene", "genesymbol", "symbols"))
  if (length(idx) == 1) colnames(genes_wes0)[idx] <- "SYMBOL"
}
if (!("total_cn" %in% cn_low)) {
  idx <- which(cn_low %in% c("totalcn", "cn", "copy_number", "copynumber"))
  if (length(idx) == 1) colnames(genes_wes0)[idx] <- "total_cn"
}

need_cols <- c("sample", "gene_class", "SYMBOL", "total_cn")
if (!all(need_cols %in% colnames(genes_wes0))) {
  stop(
    "WES records CSV must contain columns: ", paste(need_cols, collapse = ", "),
    "\nCurrent columns: ", paste(colnames(genes_wes0), collapse = ", ")
  )
}

genes_wes <- genes_wes0 %>%
  mutate(
    sample = str_trim(sample),
    gene_class = str_to_lower(str_trim(gene_class)),
    SYMBOL = str_trim(SYMBOL),
    total_cn = suppressWarnings(as.numeric(total_cn))
  ) %>%
  filter(sample %in% df_meta$Sample)

wes_class <- genes_wes %>% select(sample, gene_class) %>% distinct() %>% filter(sample %in% wes_samples)
wes_circular_samples <- wes_class %>% filter(gene_class == "circular") %>% pull(sample) %>% unique()
wes_noncircular_samples <- wes_class %>% filter(gene_class == "noncircular") %>% pull(sample) %>% unique()

# WES ecDNA gene summary (circular)
ec <- genes_wes %>%
  filter(gene_class == "circular", SYMBOL != "") %>%
  select(sample, SYMBOL) %>%
  distinct()

ec_summary <- ec %>%
  group_by(sample) %>%
  summarise(ecDNA_genes = paste(unique(SYMBOL), collapse = ", "), .groups = "drop")

write.csv(ec_summary, file.path(outdir, "ecDNA_WES_gene_summary.csv"), row.names = FALSE)

oncogenes <- read.table(oncogene_list_txt, stringsAsFactors = FALSE)$V1

ec_summary$oncogenes_in_ecDNA <- sapply(ec_summary$ecDNA_genes, function(genes_str) {
  gene_list <- unlist(strsplit(genes_str, ",\\s*"))
  hit <- intersect(gene_list, oncogenes)
  if (length(hit) > 0) paste(hit, collapse = ", ") else NA
})

write.csv(ec_summary, file.path(outdir, "ecDNA_WES_oncogene_summary.csv"), row.names = FALSE)

# =========================================================
# 6) Cross-platform sample status
# =========================================================
clean_data <- df_meta %>%
  transmute(
    Sample,
    Cancer_type,
    Status = case_when(
      Sample %in% wgs_ecDNA_samples | Sample %in% wes_circular_samples ~ "ecDNA+",
      Sample %in% wgs_amp_samples | Sample %in% wes_noncircular_samples ~ "Amplicon (Linear)",
      TRUE ~ "Negative"
    )
  )

write.csv(clean_data, file.path(outdir, "Fig1_sample_level_status_B.csv"), row.names = FALSE)

# Check WES samples with only noncircular records.
wes_should_be_amp <- setdiff(wes_noncircular_samples, wes_circular_samples)
bad_amp <- clean_data %>%
  filter(Sample %in% wes_should_be_amp, Status != "Amplicon (Linear)")

# =========================================================
# 7) Figure 1a: cancer type and genomic status donut plot
# =========================================================
cancer_order <- clean_data %>% count(Cancer_type, sort = TRUE) %>% pull(Cancer_type)
plot_data_a <- clean_data %>% group_by(Cancer_type, Status) %>% summarise(Count = n(), .groups = "drop")
plot_data_a$Cancer_type <- factor(plot_data_a$Cancer_type, levels = cancer_order)
plot_data_a$Status <- factor(plot_data_a$Status, levels = c("ecDNA+", "Amplicon (Linear)", "Negative"))

inner_data <- plot_data_a %>%
  group_by(Cancer_type) %>%
  summarise(Total = sum(Count), .groups = "drop") %>%
  mutate(
    fraction = Total / sum(Total),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n = -1)),
    label_pos = (ymin + ymax) / 2
  )

outer_data <- plot_data_a %>%
  arrange(Cancer_type, Status) %>%
  mutate(
    fraction = Count / sum(Count),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n = -1)),
    label_pos = (ymin + ymax) / 2
  )

p1a <- ggplot() +
  geom_rect(data = inner_data,
            aes(ymin=ymin, ymax=ymax, xmax=3, xmin=2, fill=Cancer_type),
            color="white", size=0.3) +
  scale_fill_manual(values=cancer_cols, name="Cancer Type", guide=guide_legend(order=1)) +
  ggnewscale::new_scale_fill() +
  geom_rect(data = outer_data,
            aes(ymin=ymin, ymax=ymax, xmax=3.5, xmin=3, fill=Status),
            color="white", size=0.2) +
  scale_fill_manual(values=status_cols, name="Genomic Status", guide=guide_legend(order=2)) +
  geom_text(data=inner_data,
            aes(x=2.5, y=label_pos,
                label=ifelse(fraction>0.03, paste0(Cancer_type,"\n(n=",Total,")"), "")),
            size=3, fontface="bold", color="white", lineheight=0.8) +
  geom_text(data=outer_data %>% filter(Status=="ecDNA+"),
            aes(x=3.25, y=label_pos, label=Count),
            size=3.5, fontface="bold", color="white") +
  coord_polar(theta="y") + xlim(c(1, 3.5)) + theme_void() +
  annotate("text", x=0, y=0, label=paste0("Total\nN=", nrow(clean_data)),
           size=3, fontface="bold") +
  theme(legend.position="right", legend.box="horizontal",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))

ggsave(file.path(outdir, "Figure1a_Donut_B.pdf"), p1a, width = 10, height = 7)

# =========================================================
# 8) Figure 1b: oncogene-positive and oncogene-negative ecDNA proportion
# =========================================================
ids_ecDNA_pos <- clean_data %>% filter(Status=="ecDNA+") %>% pull(Sample)

wgs_onco_pos_samples <- df_gene_wgs_den %>%
  filter(grepl("ecDNA", feature, ignore.case = TRUE)) %>%
  filter(grepl("^true$", is_canonical_oncogene, ignore.case = TRUE)) %>%
  pull(sample_name) %>%
  unique()

wes_onco_pos_samples <- ec_summary %>%
  filter(!is.na(oncogenes_in_ecDNA) & oncogenes_in_ecDNA != "") %>%
  pull(sample) %>%
  unique()

ids_onco_pos <- unique(c(wgs_onco_pos_samples, wes_onco_pos_samples))

stats_panel_b <- clean_data %>%
  group_by(Cancer_type) %>%
  summarise(
    Total = n(),
    N_Onco_Pos = sum(Sample %in% ids_ecDNA_pos & Sample %in% ids_onco_pos),
    N_Onco_Neg = sum(Sample %in% ids_ecDNA_pos & !Sample %in% ids_onco_pos),
    N_ecDNA_Total = sum(Sample %in% ids_ecDNA_pos),
    .groups = "drop"
  ) %>%
  mutate(
    Prop_Onco_Pos = N_Onco_Pos / Total,
    Prop_Onco_Neg = N_Onco_Neg / Total,
    Prop_Total = N_ecDNA_Total / Total
  ) %>%
  pivot_longer(
    cols = c("Prop_Onco_Pos", "Prop_Onco_Neg"),
    names_to = "Subtype",
    values_to = "Proportion"
  ) %>%
  arrange(desc(Prop_Total))

stats_panel_b$Cancer_type <- factor(stats_panel_b$Cancer_type, levels = unique(stats_panel_b$Cancer_type))
stats_panel_b$Subtype <- factor(
  stats_panel_b$Subtype,
  levels = c("Prop_Onco_Neg", "Prop_Onco_Pos"),
  labels = c("ecDNA+ (Oncogene-)", "ecDNA+ (Oncogene+)")
)

p1b <- ggplot(stats_panel_b, aes(x = Cancer_type, y = Proportion, fill = Subtype)) +
  geom_bar(stat = "identity", width = 0.7, color = "black") +
  geom_text(
    aes(y = Prop_Total, label = paste0(N_ecDNA_Total, "/", Total)),
    data = subset(stats_panel_b, Subtype == "ecDNA+ (Oncogene+)"),
    vjust = -0.5,
    size = 3
  ) +
  scale_fill_manual(values = oncogene_status_cols) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = c(0, 0)) +
  labs(x = "", y = "Proportion", fill = "ecDNA status") +
  theme_nature +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )

ggsave(file.path(outdir, "Figure1b_Prevalence_B.pdf"), p1b, width = 6, height = 4)

# =========================================================
# 9) Figure 1b simple: ecDNA-positive and ecDNA-negative counts
# =========================================================
stats_panel_b_simple <- clean_data %>%
  group_by(Cancer_type) %>%
  summarise(
    Total = n(),
    Count_Pos = sum(Status == "ecDNA+"),
    Count_Neg = sum(Status != "ecDNA+"),
    .groups = "drop"
  ) %>%
  arrange(desc(Total)) %>%
  pivot_longer(
    cols = c("Count_Neg", "Count_Pos"),
    names_to = "Status2",
    values_to = "Count"
  )

stats_panel_b_simple$Cancer_type <- factor(
  stats_panel_b_simple$Cancer_type,
  levels = unique(stats_panel_b_simple$Cancer_type)
)
stats_panel_b_simple$Status2 <- factor(
  stats_panel_b_simple$Status2,
  levels = c("Count_Neg", "Count_Pos"),
  labels = c("ecDNA-", "ecDNA+")
)

p1b_simple <- ggplot(stats_panel_b_simple, aes(x = Cancer_type, y = Count, fill = Status2)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", size = 0.3) +
  geom_text(
    data = subset(stats_panel_b_simple, Status2 == "ecDNA+"),
    aes(y = Total, label = paste0(Count, "/", Total)),
    vjust = -0.5,
    size = 3
  ) +
  scale_fill_manual(values = ecdna_binary_cols) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1)),
    breaks = scales::pretty_breaks()
  ) +
  labs(x = "", y = "Number of Samples", fill = "ecDNA Status") +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(face = "bold")
  )

ggsave(file.path(outdir, "Figure1b_Prevalence_Simple.pdf"), p1b_simple, width = 6, height = 4.5)
saveRDS(stats_panel_b_simple, file = file.path(outdir, "Figure1b_Prevalence_Simple.rds"))

# =========================================================
# 10) Figure 1b simple: metastatic samples only
# =========================================================
clean_data_met <- df_meta %>%
  filter(Pri_Met == "Met") %>%
  transmute(
    Sample,
    Cancer_type,
    Status = case_when(
      Sample %in% wgs_ecDNA_samples | Sample %in% wes_circular_samples ~ "ecDNA+",
      Sample %in% wgs_amp_samples | Sample %in% wes_noncircular_samples ~ "Amplicon (Linear)",
      TRUE ~ "Negative"
    )
  )

stats_panel_b_simple_met <- clean_data_met %>%
  group_by(Cancer_type) %>%
  summarise(
    Total = n(),
    Count_Pos = sum(Status == "ecDNA+"),
    Count_Neg = sum(Status != "ecDNA+"),
    .groups = "drop"
  ) %>%
  arrange(desc(Total)) %>%
  pivot_longer(
    cols = c("Count_Neg", "Count_Pos"),
    names_to = "Status2",
    values_to = "Count"
  )

stats_panel_b_simple_met$Cancer_type <- factor(
  stats_panel_b_simple_met$Cancer_type,
  levels = unique(stats_panel_b_simple_met$Cancer_type)
)
stats_panel_b_simple_met$Status2 <- factor(
  stats_panel_b_simple_met$Status2,
  levels = c("Count_Neg", "Count_Pos"),
  labels = c("ecDNA-", "ecDNA+")
)

p1b_simple_met <- ggplot(stats_panel_b_simple_met, aes(x = Cancer_type, y = Count, fill = Status2)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", size = 0.3) +
  geom_text(
    data = subset(stats_panel_b_simple_met, Status2 == "ecDNA+"),
    aes(y = Total, label = paste0(Count, "/", Total)),
    vjust = -0.5,
    size = 3
  ) +
  scale_fill_manual(values = ecdna_binary_cols) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1)),
    breaks = scales::pretty_breaks()
  ) +
  labs(x = "", y = "Number of Samples", fill = "ecDNA Status") +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(face = "bold")
  )

ggsave(file.path(outdir, "Figure1b_Prevalence_Simple_count_Met.pdf"), p1b_simple_met, width = 6, height = 4.5)
saveRDS(stats_panel_b_simple_met, file = file.path(outdir, "Figure1b_Prevalence_Simple_count_Met.rds"))

# =========================================================
# 11) Figure 1b: metastatic samples only, oncogene-positive/-negative proportion
# =========================================================
stats_panel_b_met <- clean_data_met %>%
  group_by(Cancer_type) %>%
  summarise(
    Total = n(),
    N_Onco_Pos = sum(Sample %in% ids_ecDNA_pos & Sample %in% ids_onco_pos),
    N_Onco_Neg = sum(Sample %in% ids_ecDNA_pos & !Sample %in% ids_onco_pos),
    N_ecDNA_Total = sum(Sample %in% ids_ecDNA_pos),
    .groups = "drop"
  ) %>%
  mutate(
    Prop_Onco_Pos = N_Onco_Pos / Total,
    Prop_Onco_Neg = N_Onco_Neg / Total,
    Prop_Total = N_ecDNA_Total / Total
  ) %>%
  pivot_longer(
    cols = c("Prop_Onco_Pos", "Prop_Onco_Neg"),
    names_to = "Subtype",
    values_to = "Proportion"
  ) %>%
  arrange(desc(Prop_Total))

stats_panel_b_met$Cancer_type <- factor(
  stats_panel_b_met$Cancer_type,
  levels = unique(stats_panel_b_met$Cancer_type)
)
stats_panel_b_met$Subtype <- factor(
  stats_panel_b_met$Subtype,
  levels = c("Prop_Onco_Neg", "Prop_Onco_Pos"),
  labels = c("ecDNA+ (Oncogene-)", "ecDNA+ (Oncogene+)")
)

p1b_met <- ggplot(stats_panel_b_met, aes(x = Cancer_type, y = Proportion, fill = Subtype)) +
  geom_bar(stat = "identity", width = 0.7, color = "black") +
  geom_text(
    aes(y = Prop_Total, label = paste0(N_ecDNA_Total, "/", Total)),
    data = subset(stats_panel_b_met, Subtype == "ecDNA+ (Oncogene+)"),
    vjust = -0.5,
    size = 3
  ) +
  scale_fill_manual(values = oncogene_status_cols) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = c(0, 0)) +
  labs(x = "", y = "Proportion", fill = "ecDNA status") +
  theme_nature +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )

ggsave(file.path(outdir, "Figure1b_Prevalence_B_ratio_Met.pdf"), p1b_met, width = 6, height = 4)
saveRDS(stats_panel_b_met, file = file.path(outdir, "Figure1b_Prevalence_B_ratio_Met.rds"))

# =========================================================
# 12) Figure 1c: WGS-only ecDNA species heterogeneity
# =========================================================
wgs_ecDNA_data <- df_gene_wgs_den %>%
  filter(grepl("ecDNA", feature, ignore.case = TRUE))

het_stats <- wgs_ecDNA_data %>%
  group_by(sample_name) %>%
  summarise(
    n_species = n_distinct(feature),
    has_onco = any(grepl("^true$", is_canonical_oncogene, ignore.case = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    category = case_when(
      n_species == 1 ~ "1 Species",
      n_species == 2 ~ "2 Species",
      TRUE ~ "≥3 Species"
    ),
    status = ifelse(has_onco, "ecDNA+ (Oncogene+)", "ecDNA+ (Oncogene-)")
  )

if (nrow(het_stats) > 0) {
  plot_data_c <- het_stats %>%
    count(category, status) %>%
    group_by(category) %>%
    mutate(prop_overall = n / nrow(het_stats),
           total_prop = sum(n) / nrow(het_stats))
  
  plot_data_c$category <- factor(plot_data_c$category, levels = c("1 Species", "2 Species", "≥3 Species"))
  plot_data_c$status <- factor(plot_data_c$status, levels = c("ecDNA+ (Oncogene-)", "ecDNA+ (Oncogene+)"))
  
  p1c <- ggplot(plot_data_c, aes(x = category, y = prop_overall, fill = status)) +
    geom_bar(stat = "identity", width = 0.6, color = "black") +
    geom_text(
      aes(y = total_prop, label = scales::percent(total_prop, accuracy = 0.1)),
      data = subset(plot_data_c, status == "ecDNA+ (Oncogene+)"),
      vjust = -0.5,
      size = 3.5
    ) +
    scale_fill_manual(values = oncogene_status_cols) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 0.8), expand = c(0, 0)) +
    labs(x = "", y = "Proportion of WGS ecDNA+ Samples", fill = "ecDNA status") +
    theme_nature +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(outdir, "Figure1c_Heterogeneity_WGSonly.pdf"), p1c, width = 4.5, height = 3.5)
} else {
  warning("No WGS ecDNA+ samples for Fig1c. Skipped.")
}

clintmp <- df_meta %>%
  filter(Seq_type == "WGS") %>%
  filter(used == "yes") %>% 
  select(Sample, Cancer_type, Pri_Met) %>%
  rename(sample_name = Sample)
clintmp$Cancer_type <- stringr::str_to_title(clintmp$Cancer_type)

plot_data_c_pri_met <- het_stats %>%
  left_join(clintmp, by = "sample_name") %>%
  count(category, Pri_Met, status) %>%
  group_by(Pri_Met) %>%
  mutate(
    prop_overall = n / sum(n),
    total_prop = scales::percent(prop_overall)
  ) %>%
  ungroup()

plot_data_c_pri_met$category <- factor(plot_data_c_pri_met$category, levels = c("1 Species", "2 Species", "≥3 Species"))
plot_data_c_pri_met$status <- factor(plot_data_c_pri_met$status, levels = c("ecDNA+ (Oncogene-)", "ecDNA+ (Oncogene+)"))

dat_tmp <- plot_data_c_pri_met[plot_data_c_pri_met$Pri_Met == "Pri", ]
p1c_Pri <- ggplot(dat_tmp, aes(x = category, y = prop_overall, fill = status)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  scale_fill_manual(values = oncogene_status_cols) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.8), expand = c(0, 0)) +
  labs(x = "", y = "Proportion of WGS ecDNA+ Samples", fill = "ecDNA status") +
  theme_nature +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dat_tmp <- plot_data_c_pri_met[plot_data_c_pri_met$Pri_Met == "Met", ]
p1c_Met <- ggplot(dat_tmp, aes(x = category, y = prop_overall, fill = status)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  scale_fill_manual(values = oncogene_status_cols) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.8), expand = c(0, 0)) +
  labs(x = "", y = "Proportion of WGS ecDNA+ Samples", fill = "ecDNA status") +
  theme_nature +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(outdir, "Figure1c_Heterogeneity_WGSonly_pri_met.pdf"), p1c_Pri+p1c_Met, width=9, height=3.5)




# =========================================================
# 13) Figure 1d: oncogene copy-number heatmaps
# =========================================================
# Build heatmap inputs with STRICT ordering + complete bar table.
build_fig1d <- function(mode=c("WGS","ALL")){
  mode <- match.arg(mode)
  
  ann_base <- df_meta %>% select(Sample, Cancer_type, Seq_type) %>% distinct()
  ann_base$Cancer_type <- sub("^(.)", "\\U\\1", ann_base$Cancer_type, perl = TRUE)
  
  # ecDNA+ samples (status defined cross-platform)
  ec_samples <- clean_data %>% filter(Status=="ecDNA+") %>% pull(Sample)
  
  # IMPORTANT: for ALL heatmap we use ONLY the platform chosen in df_meta (used==yes)
  used_wgs_samples <- ann_base %>% filter(Seq_type=="WGS", Sample %in% ec_samples) %>% pull(Sample)
  used_wes_samples <- ann_base %>% filter(Seq_type=="WES", Sample %in% ec_samples) %>% pull(Sample)
  
  # WGS oncogene CN
  wgs_long <- df_gene_wgs_den %>%
    filter(sample_name %in% used_wgs_samples) %>%
    filter(grepl("ecDNA", feature, ignore.case=TRUE)) %>%
    filter(grepl("^true$", is_canonical_oncogene, ignore.case=TRUE)) %>%
    transmute(Sample=sample_name, gene=gene, CN=as.numeric(gene_cn))
  
  # WES oncogene CN (circular + total_cn)
  wes_long <- genes_wes %>%
    filter(sample %in% used_wes_samples) %>%
    filter(gene_class=="circular", SYMBOL!="") %>%
    filter(!is.na(total_cn)) %>%
    filter(SYMBOL %in% oncogenes) %>%
    transmute(Sample=sample, gene=SYMBOL, CN=total_cn)
  
  if(mode=="WGS"){
    all_long <- wgs_long
    ann_use <- ann_base %>%
      filter(Sample %in% used_wgs_samples) %>%
      transmute(Sample, Cancer_type, Platform = Seq_type)
  } else {
    all_long <- bind_rows(wgs_long, wes_long)
    ann_use <- ann_base %>%
      filter(Sample %in% c(used_wgs_samples, used_wes_samples)) %>%
      transmute(Sample, Cancer_type, Platform = Seq_type)
  }
  
  if(nrow(all_long)==0) return(NULL)
  
  # strict gene frequency: across the SAME sample set used in this heatmap
  gene_stats <- all_long %>% distinct(Sample, gene) %>% count(gene, name="freq") %>%
    arrange(desc(freq), gene)
  sorted_genes <- gene_stats$gene
  
  # matrix rows=samples, cols=genes (CN)
  mat_df <- all_long %>%
    group_by(Sample, gene) %>%
    summarise(CN=max(CN, na.rm=TRUE), .groups="drop") %>%
    tidyr::pivot_wider(names_from=gene, values_from=CN, values_fill=0)
  
  mat <- as.matrix(mat_df[,-1])
  rownames(mat) <- mat_df$Sample
  mat <- mat[, sorted_genes, drop=FALSE]
  
  # annotation aligned STRICTLY to matrix rows
  ann_df <- ann_use %>%
    distinct(Sample, Cancer_type, Platform) %>%
    filter(Sample %in% rownames(mat)) %>%
    arrange(match(Sample, rownames(mat)))
  stopifnot(nrow(ann_df) == nrow(mat))
  
  # hard check: annotation length == matrix nrow
  stopifnot(nrow(ann_df) == nrow(mat))
  
  # bar table: complete gene × cancer, fill 0, strict order
  bar_long <- all_long %>%
    distinct(Sample, gene) %>%
    left_join(ann_use, by="Sample") %>%
    count(gene, Cancer_type, name="n") %>%
    tidyr::complete(
      gene = sorted_genes,
      Cancer_type = cancer_order_global,
      fill = list(n=0)
    ) %>%
    mutate(
      gene = factor(gene, levels=sorted_genes),
      Cancer_type = factor(Cancer_type, levels=cancer_order_global)
    ) %>%
    arrange(gene, Cancer_type)
  
  bar_df <- bar_long %>%
    tidyr::pivot_wider(names_from=Cancer_type, values_from=n, values_fill=0) %>%
    arrange(gene)
  
  bar_mat <- as.matrix(bar_df[,-1])
  rownames(bar_mat) <- as.character(bar_df$gene)
  
  # HARD CHECK: bar row order == heatmap column order
  stopifnot(identical(rownames(bar_mat), colnames(mat)))
  
  list(mat=mat, ann_df=ann_df, bar_mat=bar_mat, sorted_genes=sorted_genes)
}

plot_fig1d <- function(obj, out_pdf, title_suffix){
  if(is.null(obj)){
    warning("Fig1d skipped (no data): ", out_pdf)
    return(invisible(NULL))
  }
  mat <- obj$mat
  ann_df <- obj$ann_df
  bar_mat <- obj$bar_mat
  sorted_genes <- obj$sorted_genes
  
  # percent bars
  n_total <- nrow(mat)
  bar_mat_pct <- (bar_mat / n_total) * 100
  max_val <- max(rowSums(bar_mat_pct, na.rm=TRUE))
  my_breaks <- pretty(c(0, max_val), n=4)
  my_breaks <- my_breaks[my_breaks>=0]
  
  # align row_split + platform
  rs <- ann_df$Cancer_type[match(rownames(mat), ann_df$Sample)]
  pf <- ann_df$Platform[match(rownames(mat), ann_df$Sample)]
  if(any(is.na(rs))) stop("row_split has NA. Missing Cancer_type for some matrix rows.")
  if(any(is.na(pf))) stop("Platform has NA. Missing platform for some matrix rows.")
  
  col_fun <- circlize::colorRamp2(c(0,4,50), c("white", "#F39B7F", "#8B0000"))
  platform_cols <- c("WGS"="#4D4D4D", "WES"="#9E9E9E")
  
  left_ann <- ComplexHeatmap::rowAnnotation(
    CancerType = rs,
    Platform   = pf,
    col = list(CancerType = cancer_cols, Platform = platform_cols),
    simple_anno_size = grid::unit(0.35, "cm"),
    show_annotation_name = FALSE
  )
  
  top_ann <- ComplexHeatmap::HeatmapAnnotation(
    Frequency = ComplexHeatmap::anno_barplot(
      bar_mat_pct,
      gp = grid::gpar(fill = cancer_cols[colnames(bar_mat_pct)],
                      col="white", lwd=0.4),
      border = FALSE,
      height = grid::unit(2.2, "cm"),
      bar_width = 0.85, 
      axis_param = list(
        side = "left",
        at = my_breaks,
        labels = paste0(format(my_breaks, trim=TRUE), "%")
      )
    ),
    annotation_name_side = "left",
    annotation_name_gp = grid::gpar(fontsize=7, fontface="bold")
  )
  
  pdf(out_pdf, width=8.8, height=3.8)
  ht <- ComplexHeatmap::Heatmap(
    mat,
    name="Copy Number",
    col=col_fun,
    cluster_columns = FALSE,
    column_order = sorted_genes,
    column_names_side = "bottom",
    column_names_rot  = 90,
    column_names_gp   = grid::gpar(fontsize=8),
    cluster_rows = TRUE,
    row_split = rs,
    show_row_names = FALSE,
    row_title_rot = 0,
    row_title_gp  = grid::gpar(fontsize=9, fontface="bold"),
    rect_gp = grid::gpar(col="grey98", lwd=0.5),
    border = TRUE,
    top_annotation  = top_ann,
    left_annotation = left_ann,
    na_col="white",
    column_title = paste0("Oncogene CN Landscape ", title_suffix)
  )
  ComplexHeatmap::draw(ht, merge_legend=TRUE,
                       heatmap_legend_side="right",
                       annotation_legend_side="right")
  dev.off()
}

obj_d_wgs <- build_fig1d("WGS")
plot_fig1d(obj_d_wgs, file.path(outdir, "Figure1d_Heatmap_oncogene_WGSonly.pdf"), "(WGS-only)")

obj_d_all <- build_fig1d("ALL")
plot_fig1d(obj_d_all, file.path(outdir, "Figure1d_Heatmap_oncogene_WGSplusWES.pdf"), "(WGS+WES)")