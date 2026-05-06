#!/bin/bash

# ==============================================================================
# Global Configuration
# ==============================================================================

# 1. Base Directories
export BASE_DIR="."
export INPUT_DIR="$BASE_DIR/raw_data"
export DIR="$BASE_DIR/result/genomic"

# 2. Output Directories
export WES_OUT_ROOT="$DIR/vcf_output_WES"
export WGS_OUT_ROOT="$DIR/vcf_output_WGS"
export CNV_OUT_ROOT="$DIR/CNVkit_output"
export ECDNA_WGS_OUT="$DIR/ecDNA_output_WGS"
export ECDNA_WES_OUT="$DIR/ecDNA_output_WES"

# 3. Tools and Environment
export GATK="gatk"
export CNVKIT="cnvkit.py"
export AMPSUITE="$HOME/miniconda3/envs/ampsuite/bin/AmpliconSuite-pipeline.py"
export THREADS=40             # Threads for BWA and single-sample tools
export JAVA_OPTS="-Xmx80g"    # Memory for GATK
export CNV_PARALLEL_JOBS=25   # Number of concurrent samples to process in CNVkit

# 4. Reference files and databases
export REF="$BASE_DIR/source/GRCh38.p14.genome.fa"
export PoN="$BASE_DIR/source/pon.vcf.gz"
export Gnomad="$BASE_DIR/source/gnomad.vcf.gz"
export ComVar="$BASE_DIR/source/common_variants.vcf.gz"
export FUNC_DATA="$BASE_DIR/source/funcotator_dataSources"
export KNOWN_SITES="$BASE_DIR/source/gatk_hg38/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz"

# 5. WES & CNV Intervals
export INTERVALS="$BASE_DIR/source/intervals.bed"
export CNV_REF_BED="$BASE_DIR/source/10X/exon.bed"
export CNV_ACCESS_BED="$BASE_DIR/source/cnvkit/access-5k-mappable.hg38.bed"

# 6. Software Paths
export TRIM_GALORE="$HOME/anaconda3/envs/genome_pipeline/bin/trim_galore"
export BWA_MEM2="$BASE_DIR/bin/bwa-mem2/bwa-mem2"
export SAMTOOLS="$BASE_DIR/bin/samtools"

# 7. Sample List
export SAMPLES=("sample1" "sample2")