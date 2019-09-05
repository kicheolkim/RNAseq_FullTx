#!/bin/bash
#
#$ -S /bin/bash
#$ -o /netapp/home/userID/script/log_map
#$ -e /netapp/home/userID/script/log_map
#$ -r y
#$ -j y
#$ -l mem_free=5G
#$ -l scratch=250G
#$ -l h_rt=240:00:00
#$ -pe smp 12
#$ -R yes
#$ -t 1-580
## 'path' for output, 'task1' and 'task2' list must change to actual path and list
## 'labserver', 'userID', and 'projectName' must change to actual name

tasks1=(0	sample1_S18_L001_R1_001.fastq.gz	sample2_S17_L001_R1_001.fastq.gz)
tasks2=(0	sample1_S18_L001_R2_001.fastq.gz	sample2_S17_L001_R2_001.fastq.gz)
input1="${tasks1[$SGE_TASK_ID]}"
input2="${tasks2[$SGE_TASK_ID]}"

date
hostname
## 0. Create job-specific scratch folder
echo R1=${input1} R2=${input2}
sample=${input1%_*_*_*_*}
echo ${sample}

SCRATCH_JOB=/scratch/$USER/jobs/$JOB_ID/${sample}
mkdir -p /scratch/$USER/jobs/$JOB_ID
mkdir -p $SCRATCH_JOB

date
echo "COPY fastq read-1 FILES"
cp /labserver/data/datasets/RNAseqData/NEB_fullTx/${input1} $SCRATCH_JOB/
date
echo "COPY fastq read-2 FILES"
cp /labserver/data/datasets/RNAseqData/NEB_fullTx/${input2} $SCRATCH_JOB/

echo "Current working directory"
cd $SCRATCH_JOB/
pwd

date
echo "START bbduk trimming..." 
mkdir -p $SCRATCH_JOB/trimmed
/netapp/home/userID/bin/bbmap/bbduk.sh in=$SCRATCH_JOB/${input1} in2=$SCRATCH_JOB/${input2} out=$SCRATCH_JOB/trimmed/${input1} out2=$SCRATCH_JOB/trimmed/${input2} ref=/netapp/home/userID/bin/bbmap/resources/adapters_neb.fa t=$NSLOTS ktrim=r k=25 mink=8 hdist=1 ftm=5 qtrim=rl trimq=10 tpe tbo
echo "END: bbduk trimming..." 

date
rm -rf /labserver/home/userID/projectName/bam/${sample}
echo "START: STAR alignment ..." 
mkdir -p $SCRATCH_JOB/${sample}_star
STAR --runThreadN $NSLOTS --genomeDir /netapp/home/userID/genome/neb/starIndex/gencodeR28 --readFilesCommand gunzip -c --readFilesIn $SCRATCH_JOB/trimmed/${input1} $SCRATCH_JOB/trimmed/${input2} --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic --varVCFfile /netapp/home/userID/genome/vcf_b151/00-All.vcf.gz --waspOutputMode SAMtag --outFileNamePrefix $SCRATCH_JOB/${sample}_star/${sample}_ --outReadsUnmapped Fastx
echo "END: STAR alignment..." 

date
echo -----
rm -rf /labserver/home/userID/projectName/counts/rsem/${sample}
echo "START: RSEM..." 
mkdir -p $SCRATCH_JOB/${sample}_rsem
rsem-calculate-expression -p $NSLOTS --bam --no-bam-output --strandedness reverse --estimate-rspd --append-names --sort-bam-by-coordinate --alignments --paired-end $SCRATCH_JOB/${sample}_star/${sample}_Aligned.toTranscriptome.out.bam /netapp/home/userID/genome/neb/rsemIndex/GenCodeR28/hg38r28 $SCRATCH_JOB/${sample}_rsem/${sample}
echo "END: RSEM..." 

date
echo "START fastqc-1..."
mkdir -p $SCRATCH_JOB/fastqc
fastqc -t $NSLOTS -o $SCRATCH_JOB/fastqc --extract $SCRATCH_JOB/trimmed/${input1}
echo "END fastqc-1..."
echo -----
echo "START fastqc-2..."
fastqc -t $NSLOTS -o $SCRATCH_JOB/fastqc --extract $SCRATCH_JOB/trimmed/${input1}
echo "END fastqc-2..."

date
echo -----
echo "Move output files back to Leo" 
mv -f $SCRATCH_JOB/${sample}_star /labserver/home/userID/projectName/bam/${sample}
mv -f $SCRATCH_JOB/${sample}_rsem /labserver/home/userID/projectName/counts/rsem/${sample}
mv -f $SCRATCH_JOB/st_gtf/*.gtf /labserver/home/userID/projectName/counts/stringtie_gtf
mv -f $SCRATCH_JOB/fastqc/* /labserver/home/userID/projectName/qc/fastqc/trimmed.fastq/

echo "Remove job-specific scratch folder"
cd /scratch/
rm -rf $SCRATCH_JOB/
rmdir /scratch/$USER/jobs/$JOB_ID
rmdir /scratch/$USER/jobs
echo "========== REMOVED AND FINISHED JOBS ========="
date
