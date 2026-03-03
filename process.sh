#!/bin/bash
set -euo pipefail

INPUT_VCF="1kGP_high_coverage_Illumina.chr10.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
CHR="10"
PANEL_TSV="igsr-1000 genomes 30x on grch38-samples.tsv"
POPS=("LWK" "FIN" "CHS")
OUT="chr${CHR}_biallelic_snps"

###################################
# PROCESS VCF                     #
#    - subset to chromosome 10    #
#    - keep SNPs only             #
#    - keep biallelic sites only  #
###################################

#echo "=== Processing VCF ==="

#bcftools view \
#    -v snps \
#    -m2 -M2 \
#    -Oz \
#    -o "${OUT}.vcf.gz" \
#    "${INPUT_VCF}"

#tabix -p vcf "${OUT}.vcf.gz"

#echo "Created processed VCF: ${OUT}.vcf.gz"

###################################
# POPULATION SUBSETS              #
###################################

echo "=== Creating Population Subset ==="

for POP in "${POPS[@]}"; do
    echo ">>> Population: ${POP}"
    
    SAMPLE_LIST="samples_${POP}.txt"
    
    awk -F'\t' -v pop="${POP}" 'NR>1 && $4==pop {print $1}' "${PANEL_TSV}" > "${SAMPLE_LIST}"
    
    N_SAMPLES=$(wc -l < "${SAMPLE_LIST}")
    if [ "${N_SAMPLES}" -eq 0 ]; then
        echo " Warning: no samples found for ${POP}; skipping."
        rm -f "${SAMPLE_LIST}"
        continue
    fi
    echo " Found ${N_SAMPLES} records."
    
    POP_OUT="chr${CHR}_${POP}_biallelic_snps"
    
    bcftools view \
        -S "${SAMPLE_LIST}" \
        -Oz \
        -o "${POP_OUT}.vcf.gz" \
        "${OUT}.vcf.gz"
    
    tabix -p vcf "${POP_OUT}.vcf.gz"
    
    echo " Wrote ${POP_OUT}.vcf.gz"

done
