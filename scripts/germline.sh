#!/bin/bash
set -euo pipefail

POPS=("LWK" "FIN" "CHS" "YRI")
BITS=(8 16)
MIN_M=(0.01 1)

while [[ $# -gt 0 ]]; do
    case "$1" in
        --bits)
            shift
            BITS=()
            while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do
                BITS+=("$1")
                shift
            done
            ;;
        --min_m)
            shift
            MIN_M=()
            while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do
                MIN_M+=("$1")
                shift
            done
            ;;
        *)
            echo "Unknown argument: $1"
            exit 1
            ;;
    esac
done

for POP in "${POPS[@]}"; do
    echo "=== Processing Population: ${POP} ==="
    
    echo "  Converting PLINK bed/bim/fam to ped/map..."
    plink \
        --bfile "../plink/chr10_${POP}" \
        --recode \
        --out "../germline/chr10_${POP}"
        
done

for POP in "${POPS[@]}"; do
    for B in "${BITS[@]}"; do
        for M in "${MIN_M[@]}"; do
            
            echo "  Running GERMLINE (bits=${B}, min_m=${M})..."
            
            OUT_PREFIX="../germline/chr10_${POP}_germline_bits${B}_minm${M}"
            
            set +e
            germline \
                -input "../germline/chr10_${POP}.ped" "germline/chr10_${POP}.map" \
                -bits ${B} \
                -min_m ${M} \
                -output "${OUT_PREFIX}"
            GERMLINE_STATUS=$?
            set -e
            
            if [ ${GERMLINE_STATUS} -ne 1 ]; then
                echo "  WARNING: GERMLINE exited with status ${GERMLINE_STATUS} for ${POP} (bits=${B}, min_m=${M})"
            fi
            
            echo "    Output segments: ${OUT_PREFIX}.match"
        done
    done
done

echo "All GERMLINE commands processed."