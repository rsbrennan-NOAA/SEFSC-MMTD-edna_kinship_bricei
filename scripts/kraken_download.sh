#!/bin/bash
#SBATCH --job-name=download_core_nt
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o /home/rbrennan/rw_edna_wgs/logs/%x_%j.out
#SBATCH -e /home/rbrennan/rw_edna_wgs/logs/%x_%j.err
#SBATCH --partition standard

# Stage 2: download the prebuilt Kraken2 + Bracken core_nt database (Oct 2025 build)
# from the Langmead lab index page. Streams straight into tar so the 246 GB archive
# is never stored on disk; you only need room for the ~316 GB unpacked index.
# Needs a node with internet access.

set -euo pipefail

DB_DIR="${HOME}/kraken_databases/core_nt"   # needs ~320 GB free
URL="https://genome-idx.s3.amazonaws.com/kraken/k2_core_nt_20251015.tar.gz"
MD5_URL="https://genome-idx.s3.amazonaws.com/kraken/core_nt_20251015/core_nt.md5"

mkdir -p "${DB_DIR}"
cd "${DB_DIR}"

echo "Start: $(date)"

# download and unpack in one pass
wget -q -O - "${URL}" | tar -xzf - -C "${DB_DIR}"
echo "Download + unpack done: $(date)"

# checksum the unpacked files against the published list
wget -q -O core_nt.md5 "${MD5_URL}"
md5sum -c core_nt.md5

echo "Done: $(date)"
ls -lh "${DB_DIR}"
