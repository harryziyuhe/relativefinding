#!/bin/bash
set -euo pipefail

POPS=("LWK" "FIN" "CHS" "YRI")
MAF=0.05

for POP in "${POPS[@]}"; do
    echo "=== Processing Population: ${POP} ==="
    
    VCF="../data/chr10_${POP}_biallelic_snps.vcf.gz"
    
    if [ ! -f "${VCF}" ]; then
        echo "  ERROR: VCF not found: ${VCF}"
        continue
    fi
    
    echo "   Converting VCF to PLINK bed/bim/fam..."
    plink \
        --vcf "${VCF}" \
        --maf ${MAF} \
        --double-id \
        --keep-allele-order \
        --make-bed \
        --out "../plink/chr10_${POP}"
        
    plink \
        --bfile "../plink/chr10_${POP}" \
        --genome \
        --out "../plink/chr10_${POP}_plink"
    
    echo "  Done with ${POP}"
    echo "    PLINK relatedness: plink/chr10_${POP}_plink.genome"
done

echo "All PLINK commands processed."