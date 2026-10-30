#!/bin/bash
#SBATCH --job-name=subsample_qc
#SBATCH --array=0-5
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o /home/rbrennan/rw_edna_wgs/logs/%x_%j.out
#SBATCH -e /home/rbrennan/rw_edna_wgs/logs/%x_%j.err
#SBATCH --partition standard
 

module load bio/seqtk/1.3
module load bio/fastp/0.23.2

set -euo pipefail


# --- settings ---
RAW_DIR="${HOME}/rw_edna_wgs/data/raw"
OUT_DIR="${HOME}/rw_edna_wgs/kraken"
N_READS=10000000 # reads to keep per sample
MIN_LEN=40 

SAMPLES=(
  Bric_eDNA07_S45_L006
  Bric_eDNA08_S46_L006
  Bric_eDNA09_S47_L006
  Bric_tag01_S52_L006
  Bric045_S55_L006
  Bric026_S57_L006  
)
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

SUB_DIR="${OUT_DIR}/01_subsampled"
QC_DIR="${OUT_DIR}/02_qc"
mkdir -p "${SUB_DIR}" "${QC_DIR}"

echo "Sample: ${SAMPLE}"

N_LINES=$(( N_READS * 4 ))

# head closes the pipe early, which makes zcat exit with SIGPIPE;
# '|| true' stops set -o pipefail from killing the job
{ zcat "${RAW_DIR}/${SAMPLE}_R1_001.fastq.gz" | head -n ${N_LINES} || true; } | gzip > "${SUB_DIR}/${SAMPLE}_R1.sub.fastq.gz" &
{ zcat "${RAW_DIR}/${SAMPLE}_R2_001.fastq.gz" | head -n ${N_LINES} || true; } | gzip > "${SUB_DIR}/${SAMPLE}_R2.sub.fastq.gz" &
wait

# --- QC: adapter trimming, NovaSeq poly-G trimming, length filter ---
fastp \
  -i "${SUB_DIR}/${SAMPLE}_R1.sub.fastq.gz" \
  -I "${SUB_DIR}/${SAMPLE}_R2.sub.fastq.gz" \
  -o "${QC_DIR}/${SAMPLE}_R1.qc.fastq.gz" \
  -O "${QC_DIR}/${SAMPLE}_R2.qc.fastq.gz" \
  --detect_adapter_for_pe \
  --trim_poly_g \
  --length_required ${MIN_LEN} \
  --thread ${SLURM_CPUS_PER_TASK} \
  --json "${QC_DIR}/${SAMPLE}.fastp.json" \
  --html "${QC_DIR}/${SAMPLE}.fastp.html"

echo "Done: ${SAMPLE}"

