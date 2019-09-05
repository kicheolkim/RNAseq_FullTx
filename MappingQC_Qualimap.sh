#!/bin/bash
#
#$ -S /bin/bash
#$ -o /netapp/home/userID/script/log_qualimap
#$ -e /netapp/home/userID/script/log_qualimap
#$ -r y
#$ -j y
#$ -l mem_free=5G
#$ -l scratch=200G
#$ -l h_rt=120:00:00
#$ -R yes
#$ -pe smp 8
## 'path' for output, 'task1' and 'task2' list must change to actual path and list
## 'labserver', 'userID', and 'projectName' must change to actual name

tasks=(0	sample1_S18_L001_R1_001.fastq.gz	sample2_S17_L001_R1_001.fastq.gz)
input="${tasks[$SGE_TASK_ID]}"

date
hostname
sample=${input%_*_*_*_*}
SCRATCH_JOB=/scratch/$USER/jobs/$JOB_ID/${sample}
mkdir -p /scratch/$USER/jobs/$JOB_ID
mkdir -p $SCRATCH_JOB

echo "COPY algined bam FILE"
cp /labserver/home/userID/projectName/bam/${sample}/${sample}_Aligned.sortedByCoord.out.bam $SCRATCH_JOB/

mkdir -p $SCRATCH_JOB/${sample}
mkdir -p $SCRATCH_JOB/tmp
export JAVA_OPTS="-Djava.io.tmpdir=$SCRATCH_JOB/tmp"
qualimap rnaseq -bam ${sample}_Aligned.sortedByCoord.out.bam -gtf /netapp/home/userID/genome/gencode.v28/gencode.v28.primary_assembly.annotation.gtf -outdir $SCRATCH_JOB/${sample} -p strand-specific-reverse -pe --java-mem-size=38G

echo "Move output files back to Leo" 
mv -f $SCRATCH_JOB/${sample} /labserver/home/userID/projectName/qc/qualimap/

echo "Remove job-specific scratch folder"
cd /scratch/
rm -rf $SCRATCH_JOB/
rmdir /scratch/$USER/jobs/$JOB_ID
rmdir /scratch/$USER/jobs
echo "========== FINISHED JOBS for ${sample} =========="
date
