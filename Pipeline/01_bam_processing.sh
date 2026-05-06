#!/bin/bash
source ./config.sh

echo "========================================================"
echo "Starting Part 1: BAM Processing Pipeline"
echo "========================================================"

for spl in "${SAMPLES[@]}"; do
    echo "----------------------------------------"
    echo "[BAM Processing] Sample: $spl"
    
    OUT_DIR="$DIR/$spl"
    mkdir -p "$OUT_DIR"

    FQ1="$INPUT_DIR/${spl}_R1.fastq.gz" 
    FQ2="$INPUT_DIR/${spl}_R2.fastq.gz"

    echo "[1/5] Running Trim Galore..."
    $TRIM_GALORE --paired --quality 30 --length 36 --fastqc -o "$OUT_DIR" --cores 8 --basename "$spl" "$FQ1" "$FQ2"
    TRIM_FQ1="$OUT_DIR/${spl}_val_1.fq.gz"
    TRIM_FQ2="$OUT_DIR/${spl}_val_2.fq.gz"

    echo "[2/5] Running BWA-MEM2..."
    $BWA_MEM2 mem -t $THREADS -R "@RG\tID:$spl\tSM:$spl\tPL:Illumina" "$REF" "$TRIM_FQ1" "$TRIM_FQ2" | \
    $SAMTOOLS sort -@ $THREADS -o "$OUT_DIR/${spl}.sorted.bam" -

    echo "[3/5] Marking duplicates..."
    TMP_DIR="$OUT_DIR/tmp"
    mkdir -p "$TMP_DIR"
    $GATK MarkDuplicatesSpark \
        --java-options "-Djava.io.tmpdir=$TMP_DIR -Xmx2g" \
        -I "$OUT_DIR/${spl}.sorted.bam" \
        -O "$OUT_DIR/${spl}.markdup.bam" \
        -M "$OUT_DIR/${spl}.metrics.txt" \
        --spark-master "local[$THREADS]"
    rm -rf "$TMP_DIR"

    echo "[4/5] Calculating BQSR recalibration table..."
    $GATK BaseRecalibratorSpark \
        -I "$OUT_DIR/${spl}.markdup.bam" \
        -R "$REF" \
        --known-sites "$KNOWN_SITES" \
        -O "$OUT_DIR/${spl}.recal.table" \
        --disable-sequence-dictionary-validation \
        --spark-master "local[$THREADS]"

    echo "[5/5] Applying BQSR..."
    $GATK ApplyBQSRSpark \
        -R "$REF" \
        -I "$OUT_DIR/${spl}.markdup.bam" \
        -bqsr "$OUT_DIR/${spl}.recal.table" \
        -O "$OUT_DIR/${spl}.BQSR.bam" \
        --disable-sequence-dictionary-validation \
        --spark-master "local[$THREADS]"

    # Cleanup BAM intermediate files
    rm -f "$OUT_DIR/${spl}.sorted.bam" "$OUT_DIR/${spl}.markdup.bam"
    echo "Sample $spl BAM processing complete."
done