#!/bin/bash
#SBATCH --job-name=kraken_mapped
#SBATCH --cpus-per-task=16
#SBATCH --mem=400G
#SBATCH --time=16:00:00
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o /home/rbrennan/rw_edna_wgs/logs/%x_%j.out
#SBATCH -e /home/rbrennan/rw_edna_wgs/logs/%x_%j.err
#SBATCH --partition himem

# Classify reads that aligned well (MAPQ >= 30) to the Rice's whale genome, genome-wide.
# Samples with lots of mapped reads (e.g., the biopsy) are randomly subsampled to
# roughly TARGET_READS before filtering, to keep files and run time manageable.
#
# Reads are classified individually (single-end mode): R1, R2 and singletons all go
# into one FASTQ per sample. Pairing doesn't matter for asking "what is this read?".

module load bio/samtools/1.23
module load bio/seqtk/1.3
eval "$(conda shell.bash hook)"
conda activate kraken

set -euo pipefail

# --- settings ---
BASE="${HOME}/rw_edna_wgs"
BAM_DIR="${BASE}/PopGLen/results/mapping/bams"
BAM_SUFFIX="b_ricei.rmdup.realn.clip.bam"
DB="${HOME}/kraken_databases/core_nt"
OUT_DIR="${BASE}/kraken/05_mapped_reads"
MAPQ=30                   
TARGET_READS=5000000       # subsample samples with more mapped reads than this
CONFIDENCE=0.05

SAMPLES=(
  Bric026
  Bric045
  Bric-tag01
  Bric-eDNA07
  Bric-eDNA08
  Bric-eDNA09
)

echo -e "sample\tmapped_total\tsubsample_frac\treads_classified" > "${OUT_DIR}/read_counts.tsv"

echo -e "sample\treads_classified" > "${OUT_DIR}/read_counts.tsv"

for SAMPLE in "${SAMPLES[@]}"; do
  BAM="${BAM_DIR}/${SAMPLE}.${BAM_SUFFIX}"
  PREFIX="${OUT_DIR}/${SAMPLE}.q${MAPQ}"

  echo "starting ${SAMPLE} $(date)"

  # -F 0xF04: drop unmapped, secondary, QC-fail, duplicate, supplementary
  # seqtk keeps a random TARGET_READS reads (all of them if there are fewer)
  samtools view -u -@ 8 -F 0xF04 -q ${MAPQ} "${BAM}" |
    samtools fastq - |
    seqtk sample -s 42 - ${TARGET_READS} |
    gzip > "${PREFIX}.fq.gz"

  N=$(( $(zcat "${PREFIX}.fq.gz" | wc -l) / 4 ))
  echo -e "${SAMPLE}\t${N}" >> "${OUT_DIR}/read_counts.tsv"
  echo "    ${N} reads to classify"

  echo " "
  echo "starting kraken for ${SAMPLE}"
  kraken2 \
    --db "${DB}" \
    --threads ${SLURM_CPUS_PER_TASK} \
    --confidence ${CONFIDENCE} \
    --gzip-compressed \
    --report "${PREFIX}.k2report" \
    --output "${PREFIX}.k2out" \
    "${PREFIX}.fq.gz"
done
