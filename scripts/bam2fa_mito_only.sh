#!/bin/bash

#SBATCH --job-name=bam2mt    ## job name
#SBATCH -o /home/rbrennan/rw_edna_wgs/logs/%x_%j.out
#SBATCH -e /home/rbrennan/rw_edna_wgs/logs/%x_%j.err
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -c 1   
#SBATCH --mem=20G
#SBATCH -t 7:00:00   ## walltime in mins or mins:secs or hrs:mins:secs. 
#SBATCH --ntasks=1 

##########################################################################################
module load bio/angsd/0.940
module load bio/samtools/1.11
module load bio/picard/2.23.9
module load bio/bbmap/39.06 
##############################################################################

READ_DIR=/home/rbrennan/rw_edna_wgs/PopGLen_mito/results/mapping/bams
RUN="mito_only"

OUTDIR=~/rw_edna_wgs/analysis/mito/$RUN
#OUTDIR=~/rw_edna_wgs/analysis/mito/mito_only

TEMP_DIR=/home/rbrennan/rw_edna_wgs/tmp/

Min_cov=5 # minimum coverage required for consensus mtDNA (in Angsd)

##############################################################################

cd ${READ_DIR}

echo ${READ_DIR}

echo "indiv,mito_reads" > ${OUTDIR}/mito_read_counts_${RUN}.csv

for file in *b_ricei.rmdup.realn.clip.bam 
do

echo $file

shortID=`echo $file | cut -f1 -d "."`

echo $shortID

# pull out mito reads from bam

samtools index ${READ_DIR}/${shortID}.b_ricei.rmdup.realn.clip.bam

samtools view -F4 -bh ${READ_DIR}/${shortID}.b_ricei.rmdup.realn.clip.bam \
  mitochondrion_40bp_pad -o ${OUTDIR}/${shortID}.mito.bam

samtools index ${OUTDIR}/${shortID}.mito.bam

# convert nuclear DNA from bam to fastq (2 files per bam file if paired reads)
samtools bam2fq ${OUTDIR}/${shortID}.mito.bam > ${OUTDIR}/${shortID}.${RUN}.fq
cat ${OUTDIR}/${shortID}.${RUN}.fq | grep '^@.*/1$' -A 3 --no-group-separator > ${OUTDIR}/${shortID}.${RUN}_R1.fq
cat ${OUTDIR}/${shortID}.${RUN}.fq | grep '^@.*/2$' -A 3 --no-group-separator > ${OUTDIR}/${shortID}.${RUN}_R2.fq

repair.sh in=${OUTDIR}/${shortID}.${RUN}_R1.fq in2=${OUTDIR}/${shortID}.${RUN}_R2.fq out=${OUTDIR}/${shortID}.${RUN}_repaired_R1.fq out2=${OUTDIR}/${shortID}.${RUN}_repaired_R2.fq outs=${OUTDIR}/${shortID}.${RUN}_orphan.fq


# output consensus
angsd -i ${OUTDIR}/${shortID}.mito.bam -out ${OUTDIR}/${shortID}.${RUN}.${Min_cov}x_consensus -doCounts 1 -doFasta 2 -setMinDepth $Min_cov -r mitochondrion_40bp_pad:

zcat ${OUTDIR}/${shortID}.${RUN}.${Min_cov}x_consensus.fa.gz > ${OUTDIR}/${shortID}.${RUN}.${Min_cov}x_consensus.fa

# count number of reads for each indiv mapping to mt

NUM_READS=`samtools view -c ${OUTDIR}/${shortID}.mito.bam`
echo ${shortID},${NUM_READS} >> ${OUTDIR}/mito_read_counts_${RUN}.csv

# -dofasta
# 	1: use a random (non N) base (needs -doCounts 1)
# 	2: use the most common (non N) base (needs -doCounts 1)
#	 In the case of ties a random base is chosen among the bases with the same maximum counts. N's or filtered based are ignored. The "-doCounts 1" options for allele counts is needed in order to determine the most common base. If multiple individuals are used the four bases are counted across individuals.
##############################################

sleep 2 

done

cd ${OUTDIR}

# change mitogenome sequence names to fasta file name.
for i in *.fa; do sed -i "1s/.*/>${i%.fa}/" $i; done
