#!/bin/bash
source ./config.sh

echo "========================================================"
echo "Starting Part 2: WES Mutation Calling Pipeline"
echo "========================================================"
mkdir -p "$WES_OUT_ROOT"

for spl in "${SAMPLES[@]}"; do
    echo "----------------------------------------"
    echo "[WES Mutect2] Sample: $spl"
    bamfile="$DIR/$spl/${spl}.BQSR.bam"
    outdir="$WES_OUT_ROOT/$spl"
    mkdir -p "$outdir"

    if [ -f "$outdir/${spl}_somatic_PASS_ONLY.maf" ]; then
        echo "Sample completed, skipping..."
        continue
    fi

    echo "[1/6] Mutect2 Calling (with Intervals)..."
    $GATK --java-options "$JAVA_OPTS" Mutect2 -R "$REF" -I "$bamfile" \
        --panel-of-normals "$PoN" --germline-resource "$Gnomad" \
        -O "$outdir/unfilters.vcf.gz" --f1r2-tar-gz "$outdir/f1r2.tar.gz" \
        --native-pair-hmm-threads "$THREADS" -L "$INTERVALS"

    echo "[2/6] Calculating Contamination..."
    $GATK --java-options "$JAVA_OPTS" GetPileupSummaries -I "$bamfile" -V "$ComVar" -L "$ComVar" -O "$outdir/getpileupsummaries.table"
    $GATK --java-options "$JAVA_OPTS" CalculateContamination -I "$outdir/getpileupsummaries.table" -O "$outdir/contamination.table"
    
    echo "[3/6] Learn Read Orientation Model..."
    $GATK --java-options "$JAVA_OPTS" LearnReadOrientationModel -I "$outdir/f1r2.tar.gz" -O "$outdir/read-orientation-model.tar.gz"
    
    echo "[4/6] Filtering..."
    $GATK --java-options "$JAVA_OPTS" FilterMutectCalls -R "$REF" -V "$outdir/unfilters.vcf.gz" \
        --contamination-table "$outdir/contamination.table" --ob-priors "$outdir/read-orientation-model.tar.gz" \
        -O "$outdir/somatic_filtered.vcf.gz"

    echo "[5/6] Selecting PASS variants..."
    $GATK SelectVariants -R "$REF" -V "$outdir/somatic_filtered.vcf.gz" -O "$outdir/somatic_PASS_ONLY.vcf.gz" --exclude-filtered

    echo "[6/6] Annotating with Funcotator..."
    $GATK Funcotator -R "$REF" -V "$outdir/somatic_PASS_ONLY.vcf.gz" -O "$outdir/${spl}_somatic_PASS_ONLY.maf" \
        --output-file-format MAF --data-sources-path "$FUNC_DATA" --ref-version hg38

    if [ -f "$outdir/${spl}_somatic_PASS_ONLY.maf" ]; then
        echo "Cleaning up intermediate files..."
        rm -f "$outdir/unfilters.vcf.gz"* "$outdir/f1r2.tar.gz" "$outdir/getpileupsummaries.table"
        rm -f "$outdir/read-orientation-model.tar.gz" "$outdir/somatic_filtered.vcf.gz"*
    fi
done