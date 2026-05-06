#!/bin/bash
source ./config.sh

echo "========================================================"
echo "Starting Part 4: CNV Calling Pipeline"
echo "========================================================"
mkdir -p "$CNV_OUT_ROOT"

process_cnv_sample() {
    sample_name="$1"
    bam_path="$DIR/$sample_name/${sample_name}.BQSR.bam"
    
    if [ ! -f "$bam_path" ]; then 
        echo "[CNVkit] Error: BAM for $sample_name not found."
        return 1
    fi

    sample_dir="$CNV_OUT_ROOT/$sample_name"
    mkdir -p "$sample_dir"
    cd "$sample_dir" || exit

    echo "[CNVkit] Processing $sample_name"
    $CNVKIT batch "$bam_path" -n -t "$CNV_REF_BED" -f "$REF" --access "$CNV_ACCESS_BED" --output-reference ref.cnn -d .
    $CNVKIT segment -t 0.001 "${sample_name}.BQSR.cnr" -o "${sample_name}.BQSR.cns"
    $CNVKIT call -m none "${sample_name}.BQSR.cns" --drop-low-coverage --center mode -o "${sample_name}.BQSR.call.cns"
    sed -i '/chrM/d' "${sample_name}.BQSR.call.cns"
    $CNVKIT heatmap "${sample_name}.BQSR.call.cns" -d -o "${sample_name}_cnv.pdf"
    echo "[CNVkit] Finished $sample_name"
}

# Export function for GNU parallel
export -f process_cnv_sample

# Run in parallel
printf "%s\n" "${SAMPLES[@]}" | parallel -j "$CNV_PARALLEL_JOBS" process_cnv_sample