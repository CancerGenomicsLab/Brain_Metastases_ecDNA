#!/bin/bash
source ./config.sh

echo "========================================================"
echo "Starting Part 5: WGS ecDNA Detection (AmpliconSuite)"
echo "========================================================"
mkdir -p "$ECDNA_WGS_OUT"

for spl in "${SAMPLES[@]}"; do
    echo "----------------------------------------"
    echo "[AmpliconSuite] Processing Sample: $spl"
    
    outdir="$ECDNA_WGS_OUT/$spl"
    mkdir -p "$outdir"
    cd "$outdir" || exit

    FQ1="$DIR/$spl/${spl}_val_1.fq.gz"
    FQ2="$DIR/$spl/${spl}_val_2.fq.gz"

    if [ ! -f "$FQ1" ] || [ ! -f "$FQ2" ]; then
        echo "Error: Trimmed FastQ files for $spl not found. Skipping."
        cd "$BASE_DIR" || exit
        continue
    fi

    $AMPSUITE -s "$spl" -t $THREADS --fastqs "$FQ1" "$FQ2" --ref hg38 --run_AA --run_AC
    
    cd "$BASE_DIR" || exit
done