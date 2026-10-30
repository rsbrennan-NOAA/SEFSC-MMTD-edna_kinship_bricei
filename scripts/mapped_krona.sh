#!/bin/bash
#SBATCH --job-name=krona_mapped
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o /home/rbrennan/rw_edna_wgs/logs/%x_%j.out
#SBATCH -e /home/rbrennan/rw_edna_wgs/logs/%x_%j.err
#SBATCH --partition standard

# Krona plot of Kraken2 results for reads that mapped to the Rice's whale genome (MAPQ >= 30).
# One tab per sample.

eval "$(conda shell.bash hook)"
conda activate kraken

set -euo pipefail

IN_DIR="${HOME}/rw_edna_wgs/kraken/05_mapped_reads"
KRONA_DIR="${HOME}/rw_edna_wgs/kraken/06_mapped_krona"

SAMPLES=(
  Bric026
  Bric045
  Bric-tag01
  Bric-eDNA07
  Bric-eDNA08
  Bric-eDNA09
)

ARGS=()
for s in "${SAMPLES[@]}"; do
  kreport2krona.py -r "${IN_DIR}/${s}.q30.k2report" -o "${KRONA_DIR}/${s}.mapped_q30.krona.txt"
  ARGS+=("${KRONA_DIR}/${s}.mapped_q30.krona.txt,${s}")
done

ktImportText "${ARGS[@]}" -o "${KRONA_DIR}/mapped_q30_krona.html"

echo "Done: ${KRONA_DIR}/mapped_q30_krona.html"
