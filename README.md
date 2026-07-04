# Mapping the dynamics of ecDNA evolution during metastasis to the brain


## Project overview
Extrachromosomal DNA (ecDNA) amplifications are key drivers of tumor evolution and associated with poor outcomes, yet their role in human brain metastases (BrM)-the most common and incurable metastatic brain tumors—remains poorly understood.  
Here, we present the first comprehensive brain metastasis ecDNA map (eMAP), revealing the landscape and clinical impact of ecDNA in BrM. 

## Study design
![Study design and cohort composition](assets/picture/Figure%201.jpg)

(A) Study design and cohort composition of the brain metastasis cohort, including cancer type distribution, paired longitudinal samples, in-house WGS/WES, and the ecDNA analysis pipeline.  
(B) Cohort-level summary of cancer types and genomic status, including ecDNA-positive, linear amplicon, and negative samples.

## Key Findings

#### Figure 1: ecDNA is recurrently detected in brain metastases across multiple primary tumor types

#### Figure 2: ecDNA-positive brain metastases exhibit distinct genomic features and significant oncogenic amplification

#### Figure 3: Paired longitudinal analysis reveals evolutionary dynamics of ecDNA in brain metastasis

#### Figure 4: A representative paired gastric cancer case illustrates retention and remodeling of ecDNA-associated oncogenic content during brain metastasis

#### Figure 5: Clinical association of ecDNA with brain metastasis timing and survival


## Summary
This study systematically characterizes ecDNA in brain metastases across multiple tumor types. EcDNA was recurrently detected in brain metastasis specimens, associated with increased genomic instability, oncogene amplification, and positive selection, and showed conservation, acquisition, and remodeling during primary-to-brain metastatic progression. Clinically, ecDNA positivity was associated with shorter latency to brain metastasis and worse overall survival, supporting ecDNA-associated genome remodeling as a candidate mechanism of intracranial metastatic adaptation.

---

## Repository purpose

This repository provides code associated with the manuscript. The scripts were used for sequencing-data preprocessing, ecDNA/focal amplification analysis, genomic feature analysis, paired primary–brain metastasis comparisons, and survival analysis.

The repository is intended to document the analytical workflow and support transparency of the study. It is not designed as a one-command reproducible pipeline. Running the scripts requires access to the corresponding sequencing data, processed intermediate files, reference genomes, sample metadata, and locally configured software environments.


## Repository structure
This repository contains scripts used to analyze ecDNA-related genomic features, paired primary-metastatic dynamics, and survival associations in brain metastasis samples.

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

The `Pipeline/` directory contains shell scripts used for upstream sequencing-data processing in this study. These scripts document the preprocessing strategy, including read alignment, somatic mutation calling, copy-number analysis, and ecDNA/focal CNA detection.

Global paths, reference files, software paths, thread numbers, and sample identifiers are defined in Pipeline/config.sh. The scripts were developed for the computing environment used in this study and may require adaptation before reuse in other environments.

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

The following commands illustrate the order in which the upstream scripts were used in this study after configuring `Pipeline/config.sh`:

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
The upstream preprocessing scripts also require command-line tools including Trim Galore, BWA-MEM2, SAMtools, GATK, CNVkit, AmpliconSuite, GCAP, and GNU parallel. Tool paths and reference files should be configured in `Pipeline/config.sh`.

## Running the figure analysis

The `Figures/` directory contains R scripts used to generate manuscript-related analyses and figures from processed study-specific data tables. These scripts depend on intermediate files and metadata generated during upstream analysis and are provided as reference code for the analyses reported in the manuscript.

The following commands illustrate the figure scripts associated with the main analyses:

```bash
Rscript Figures/Figure1.R
Rscript Figures/Figure2.R
Rscript Figures/Figure3.R
Rscript Figures/Figure5.R
```
