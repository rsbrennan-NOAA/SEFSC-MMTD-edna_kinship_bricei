#!/bin/bash
#SBATCH --job-name=clip_redo
#SBATCH --time=1-00:00:00
#SBATCH --mem=24G
#SBATCH -c 4
#SBATCH -n 1
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o /home/rbrennan/rw_edna_wgs/logs/%x_%j.out
#SBATCH -e /home/rbrennan/rw_edna_wgs/logs/%x_%j.err
#SBATCH --partition standard
 
 
source ~/.bashrc
conda activate popglen

module load bio/bamutil/1.0.5

cd ~/rw_edna_wgs/PopGLen

bam clipOverlap \
  --in results/mapping/bams/Bric-tag01.b_ricei.rmdup.realn.bam \
  --out results/mapping/bams/Bric-tag01.b_ricei.rmdup.realn.clip.bam \
  --stats \
  2> logs/mapping/bamutil/clipoverlap/Bric-tag01.b_ricei.rmdup.realn.log

echo "Exit code was: $?"
