#!/bin/bash
source ./config.sh

echo "========================================================"
echo "Starting Part 6: WES ecDNA Detection (GCAP Batch)"
echo "========================================================"
mkdir -p "$ECDNA_WES_OUT"

R_SCRIPT="$ECDNA_WES_OUT/run_gcap.R"

# Dynamically generate R script
cat << 'EOF' > "$R_SCRIPT"
library(gcap)
library(data.table)
library(org.Hs.eg.db)

base_dir <- Sys.getenv("CNV_OUT_ROOT")
outdir <- Sys.getenv("ECDNA_WES_OUT")

sample_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)

all_samples_list <- lapply(sample_dirs, function(d) {
  target_file <- list.files(d, pattern = "\\.BQSR\\.call\\.cns$", full.names = TRUE)
  if (length(target_file) == 0) return(NULL)
  
  s_name <- basename(d)
  dt <- fread(target_file[1])
  
  data.table(
    sample = s_name,
    chromosome = as.character(dt$chromosome),
    start = as.numeric(dt$start),
    end = as.numeric(dt$end),
    total_cn = as.numeric(dt$cn),
    minor_cn = NA_real_,
    purity = 1
  )
})

gcap_input_big <- rbindlist(all_samples_list)

res <- gcap.ASCNworkflow(
  data = gcap_input_big,
  genome_build = "hg38",
  model = "XGB11",
  tightness = 1L,
  gap_cn = 3L,
  overlap = 1,
  only_oncogenes = FALSE,
  outdir = outdir,
  result_file_prefix = "batch_wes_run"
)

# Gene Name Conversion
file_path <- file.path(outdir, "batch_wes_run_fCNA_records.csv")
if (file.exists(file_path)) {
    dt <- fread(file_path)
    ids_to_convert <- as.character(unique(dt$gene_id))
    gene_map <- select(org.Hs.eg.db, keys = ids_to_convert, columns = c("SYMBOL"), keytype = "ENSEMBL") 
    gene_map <- as.data.table(gene_map)
    setnames(gene_map, "ENSEMBL", "gene_id")
    
    dt_final <- merge(dt, gene_map, by = "gene_id", all.x = TRUE)
    
    out_csv <- file.path(outdir, "batch_wes_run_fCNA_records_with_symbols.csv")
    fwrite(dt_final, out_csv)
    cat("GCAP complete. Output saved to:", out_csv, "\n")
}
EOF

# Execute R script
echo "Executing GCAP R script..."
Rscript "$R_SCRIPT"