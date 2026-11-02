#!/bin/bash
#SBATCH --job-name=kraken_bracken
#SBATCH --cpus-per-task=16
#SBATCH --mem=400G
#SBATCH --time=24:00:00
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o /home/rbrennan/rw_edna_wgs/logs/%x_%j.out
#SBATCH -e /home/rbrennan/rw_edna_wgs/logs/%x_%j.err
#SBATCH --partition himem

# classify reads with kraken2 using core_nt, re-estimate abundance with Bracken,
# and pull out a domain-level summary (Bacteria / Archaea / Eukaryota / Viruses / unclassified).

# activate conda env before 'set -u' (conda's activate scripts trip on unset variables)
eval "$(conda shell.bash hook)"
conda activate kraken    

set -euo pipefail

# --- settings ---
DB="${HOME}/kraken_databases/core_nt"
QC_DIR="${HOME}/rw_edna_wgs/kraken/02_qc"
OUT_DIR="${HOME}/rw_edna_wgs/kraken/03_kraken"
CONFIDENCE=0.05        # reduces low-abundance false positives
READ_LEN=100      

SAMPLES=(
  Bric_eDNA07_S45_L006
  Bric_eDNA08_S46_L006
  Bric_eDNA09_S47_L006
  Bric_tag01_S52_L006
  Bric045_S55_L006
  Bric026_S57_L006
)

mkdir -p "${OUT_DIR}"

for SAMPLE in "${SAMPLES[@]}"; do

  echo "${SAMPLE}: Kraken start"

  kraken2 \
    --db "${DB}" \
    --threads ${SLURM_CPUS_PER_TASK} \
    --gzip-compressed \
    --confidence ${CONFIDENCE} \
    --report "${OUT_DIR}/${SAMPLE}.k2report" \
    --output "${OUT_DIR}/${SAMPLE}.k2out" \
    "${QC_DIR}/${SAMPLE}_R1.qc.fastq.gz" \
    "${QC_DIR}/${SAMPLE}_R2.qc.fastq.gz"

  echo "=== ${SAMPLE}: Bracken start $(date) ==="

  # species and genus level abundance estimation
  for LEVEL in S G; do
    bracken \
      -d "${DB}" \
      -i "${OUT_DIR}/${SAMPLE}.k2report" \
      -o "${OUT_DIR}/${SAMPLE}.bracken_${LEVEL}.txt" \
      -w "${OUT_DIR}/${SAMPLE}.bracken_${LEVEL}.k2report" \
      -r ${READ_LEN} \
      -l ${LEVEL}
  done

  # domain-level summary from the Kraken report- keeps unclassified reads
  # columns: percent of reads, read count, rank, name
  awk -F'\t' '$4=="U" || $4=="D" || $4=="R1" {gsub(/^ +/, "", $6); print $1"\t"$2"\t"$4"\t"$6}' \
    "${OUT_DIR}/${SAMPLE}.k2report" > "${OUT_DIR}/${SAMPLE}.domain_summary.tsv"

  echo "${SAMPLE} done"

done
