#!/bin/bash
set -euo pipefail

POPS=("LWK" "FIN" "CHS")

for POP in "${POPS[@]}"; do
    echo "=== Processing Population: ${POP} ==="
    echo "  Converting PLINK bed/bim/fam to ped/map..."
    plink \
        --bfile "plink/chr10_${POP}" \
        --recode \
        --out "germline/chr10_${POP}"
    
    set +e
    germline \
        -input "germline/chr10_${POP}.ped" "germline/chr10_${POP}.map" \
        -output "germline/chr10_${POP}_germline"
    GERMLINE_STATUS=$?
    set -e
    
    if [ ${GERMLINE_STATUS} -ne 1 ]; then
        echo "  WARNING: GERMLINE exited with status ${GERMLINE_STATUS} for ${POP}"
    fi
    
    echo "  Done with ${POP}"
    echo "    GERMLINE segments: germline/chr10_${POP}_germline.match"
done

echo "All GERMLINE commands processed."