#!/bin/bash
#SBATCH --job-name=krona
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o /home/rbrennan/rw_edna_wgs/logs/%x_%j.out
#SBATCH -e /home/rbrennan/rw_edna_wgs/logs/%x_%j.err
#SBATCH --partition standard


eval "$(conda shell.bash hook)"
conda activate kraken

cd ~/rw_edna_wgs/kraken/03_kraken/

for s in Bric_eDNA07_S45_L006 Bric_eDNA08_S46_L006 Bric_eDNA09_S47_L006 Bric_tag01_S52_L006 Bric045_S55_L006 Bric026_S57_L006; do
  kreport2krona.py -r ${s}.k2report -o ../04_krona/${s}.krona.txt
done

ktImportText \
  ../04_krona/Bric_eDNA07_S45_L006.krona.txt,eDNA07 \
  ../04_krona/Bric_eDNA08_S46_L006.krona.txt,eDNA08 \
  ../04_krona/Bric_eDNA09_S47_L006.krona.txt,eDNA09 \
  ../04_krona/Bric_tag01_S52_L006.krona.txt,tag01 \
  ../04_krona/Bric045_S55_L006.krona.txt,Bric045 \
  ../04_krona/Bric026_S57_L006.krona.txt,Bric026 \
  -o ../04_krona/Bric_eDNA__krona.html
