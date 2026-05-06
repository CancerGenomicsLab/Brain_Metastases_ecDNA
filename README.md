# Mapping the dynamics of ecDNA evolution during metastasis to the brain


## Project overview
Extrachromosomal DNA (ecDNA) amplifications are key drivers of tumor evolution and associated with poor outcomes, yet their role in human brain metastases (BrM)-the most common and incurable metastatic brain tumors—remains poorly understood.  
Here, we present the first comprehensive brain metastasis ecDNA map (eMAP), revealing the landscape and clinical impact of ecDNA in BrM. 

## Study design
![Study design and cohort composition](assets/picture/Figure%201.jpg)

(A) Study design and cohort composition of the brain metastasis cohort, including cancer type distribution, paired longitudinal samples, in-house WGS, and the ecDNA analysis pipeline.  
(B) Cohort-level summary of cancer types and genomic status, including ecDNA-positive, linear amplicon, and negative samples.

## Key Findings

#### Figure 1: ecDNA is recurrently detected in brain metastases across multiple primary tumor types

#### Figure 2: ecDNA-positive brain metastases exhibit distinct genomic features and stronger oncogenic amplification

#### Figure 3: Paired longitudinal analysis reveals evolutionary dynamics of ecDNA in brain metastasis

#### Figure 4: A representative paired gastric cancer case illustrates retention and remodeling of ecDNA-associated oncogenic content during brain metastasis

#### Figure 5: Clinical association of ecDNA with brain metastasis timing and survival


## Summary
This study provides a systematic characterization of ecDNA in brain metastases across multiple tumor types. The results show that ecDNA is frequently present, associated with increased genomic instability and oncogene amplification, exhibits variable patterns between primary and metastatic tumors, and correlates with clinical outcomes.

---

## Repository structure
This repository contains R scripts used to analyze ecDNA-related genomic features, paired primary-metastatic dynamics, and survival associations in brain metastasis samples.


```text
.
├── Pipeline/
│   ├── config.sh
│   ├── 01_bam_processing.sh
│   ├── 02_wes_mutation.sh
│   ├── 03_wgs_mutation.sh
│   ├── 04_cnv_calling.sh
│   ├── 05_ecdna_wgs_ampsuite.sh
│   └── 06_ecdna_wes_gcap.sh
├── Figures/
│   ├── Figure1.R
│   ├── Figure2.R
│   ├── Figure3.R
│   └── Figure5.R
└── assets/
```


## Pipeline preprocessing

The `Pipeline/` directory contains shell scripts for upstream sequencing-data processing. Global paths, reference files, software paths, thread numbers, and sample IDs are defined in `Pipeline/config.sh`. The pipeline starts from paired FASTQ files, generates BQSR-processed BAM files, performs WES/WGS somatic mutation calling, runs CNVkit copy-number calling, and detects ecDNA/focal CNA events using AmpliconSuite for WGS and GCAP for WES.

The pipeline scripts are organized as follows:

```text
Pipeline/config.sh                   Global paths, references, tools, and sample list
Pipeline/01_bam_processing.sh        FASTQ trimming, BWA-MEM2 alignment, duplicate marking, and BQSR
Pipeline/02_wes_mutation.sh          WES Mutect2 somatic mutation calling and Funcotator annotation
Pipeline/03_wgs_mutation.sh          WGS Mutect2 somatic mutation calling and Funcotator annotation
Pipeline/04_cnv_calling.sh           CNVkit copy-number calling, segmentation, and visualization
Pipeline/05_ecdna_wgs_ampsuite.sh    WGS ecDNA detection using AmpliconSuite
Pipeline/06_ecdna_wes_gcap.sh        WES focal CNA / ecDNA inference using GCAP
```

Run the upstream pipeline from the repository root after editing `Pipeline/config.sh`:

```bash
bash Pipeline/01_bam_processing.sh
bash Pipeline/02_wes_mutation.sh
bash Pipeline/03_wgs_mutation.sh
bash Pipeline/04_cnv_calling.sh
bash Pipeline/05_ecdna_wgs_ampsuite.sh
bash Pipeline/06_ecdna_wes_gcap.sh
```

- `01_bam_processing.sh` trims paired FASTQ files with Trim Galore, aligns reads with BWA-MEM2, sorts BAM files with SAMtools, marks duplicates with GATK, and applies BQSR to generate final `.BQSR.bam` files.
- `02_wes_mutation.sh` and `03_wgs_mutation.sh` run GATK Mutect2, contamination estimation, read-orientation filtering, PASS variant selection, and Funcotator MAF annotation for WES and WGS samples, respectively.
- `04_cnv_calling.sh` runs CNVkit batch processing, segmentation, copy-number calling, mitochondrial chromosome removal, and CNV heatmap generation.
- `05_ecdna_wgs_ampsuite.sh` runs AmpliconSuite with AmpliconArchitect and AmpliconClassifier on trimmed WGS FASTQ files.
- `06_ecdna_wes_gcap.sh` generates a GCAP input table from CNVkit `.call.cns` files, runs `gcap.ASCNworkflow`, and converts Ensembl gene IDs to gene symbols for downstream figure scripts.

## Requirements
The analyses were performed in R. Required packages include:

```r
library(tidyverse)
library(tidyverse)
library(data.table)
library(readxl)
library(readr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggbeeswarm)
library(patchwork)
library(scales)
library(ggrepel)
library(ggalluvial)
library(ggnewscale)
library(survival)
library(survminer)
library(broom)
library(circlize)
library(ComplexHeatmap)
library(CINmetrics)
```
The upstream preprocessing pipeline also requires command-line tools including Trim Galore, BWA-MEM2, SAMtools, GATK, CNVkit, AmpliconSuite, GCAP, and GNU parallel. Tool paths and reference files should be configured in `Pipeline/config.sh`.

## Running the figure analysis

After upstream pipeline outputs are prepared, run the figure scripts from the repository root:

```bash
Rscript Figures/Figure1.R
Rscript Figures/Figure2.R
Rscript Figures/Figure3.R
Rscript Figures/Figure5.R
```
