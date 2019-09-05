#!/bin/bash
#
#$ -S /bin/bash
#$ -o /wynton/home/labName/userID/script/variant_call/log
#$ -e /wynton/home/labName/userID/script/variant_call/log
#$ -r y
#$ -j y
#$ -l mem_free=8G
#$ -l scratch=200G
#$ -l h_rt=300:00:00
#$ -R yes
#$ -pe smp 8
#$ -t 1-194
## tasks are text file, list of bam files from same subject

tasks=(0	13311d.list	22012b.list)
input="${tasks[$SGE_TASK_ID]}"

date
hostname
sample=${input%.*}
SCRATCH_JOB=/scratch/$USER/jobs/$JOB_ID/${sample}
mkdir -p /scratch/$USER/jobs/$JOB_ID
mkdir -p $SCRATCH_JOB
mkdir -p $SCRATCH_JOB/tmp
mkdir -p /wynton/scratch/userID/tmp/${sample}
export _JAVA_OPTIONS="-Djava.io.tmpdir=$SCRATCH_JOB/tmp"
export TMP_DIR=$SCRATCH_JOB/tmp

ref_fasta=/wynton/home/baranzinilab/userID/genome/gencode.v28/GRCh38.primary_assembly.genome.fa
dbSNP_vcf=/wynton/home/baranzinilab/userID/genome/dbsnp_b151/00-All_fixed.vcf.gz

cd $SCRATCH_JOB/
pwd

echo "AddOrReplaceReadGroups"
java -jar /netapp/home/userID/bin/picard.jar AddOrReplaceReadGroups I=/wynton/scratch/userID/merged_bam/${sample}_merged.bam O=$SCRATCH_JOB/${sample}.bam SO=coordinate RGLB=L001 RGPL=illumina RGPU=novaseq RGSM=${sample} TMP_DIR=$SCRATCH_JOB/tmp
rm $SCRATCH_JOB/tmp/*.*

echo "FixMateInformation"
java -jar /netapp/home/userID/bin/picard.jar FixMateInformation I=$SCRATCH_JOB/${sample}.bam O=$SCRATCH_JOB/${sample}_fixed.bam TMP_DIR=$SCRATCH_JOB/tmp
rm -f $SCRATCH_JOB/${sample}.bam
rm $SCRATCH_JOB/tmp/*.*

echo "MarkDuplicates"
java -jar /netapp/home/userID/bin/picard.jar MarkDuplicates I=$SCRATCH_JOB/${sample}_fixed.bam O=$SCRATCH_JOB/${sample}_dedupped.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=$SCRATCH_JOB/${sample}_dedupped.metrics TMP_DIR=$SCRATCH_JOB/tmp
rm -f $SCRATCH_JOB/${sample}_fixed.bam
rm $SCRATCH_JOB/tmp/*.*

echo "SplitNCigarReads"
java -Djava.io.tmpdir=$SCRATCH_JOB/tmp -jar /netapp/home/userID/bin/GenomeAnalysisTK-3.8-1-0/GenomeAnalysisTK.jar -T SplitNCigarReads -R ${ref_fasta} -I $SCRATCH_JOB/${sample}_dedupped.bam -o $SCRATCH_JOB/${sample}_split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
rm -f $SCRATCH_JOB/${sample}_dedupped.ba*
rm $SCRATCH_JOB/tmp/*.*

echo "BaseRecalibrator"
gatk --java-options "-Djava.io.tmpdir=$SCRATCH_JOB/tmp" BaseRecalibrator -R ${ref_fasta} -I $SCRATCH_JOB/${sample}_split.bam --use-original-qualities -O recal.table -known-sites ${dbSNP_vcf}
rm $SCRATCH_JOB/tmp/*.*

echo "ApplyBQSR"
gatk --java-options "-Djava.io.tmpdir=$SCRATCH_JOB/tmp" ApplyBQSR --add-output-sam-program-record -R ${ref_fasta} -I $SCRATCH_JOB/${sample}_split.bam --use-original-qualities -O /wynton/scratch/kkim/bam_bqsr/${sample}_BQSR.bam --bqsr-recal-file recal.table
rm -f $SCRATCH_JOB/${sample}_split.bam
rm $SCRATCH_JOB/tmp/*.*

echo "HaplotypeCaller"
gatk --java-options "-Djava.io.tmpdir=$SCRATCH_JOB/tmp" HaplotypeCaller -R ${ref_fasta} -I /wynton/scratch/kkim/bam_bqsr/${sample}_BQSR.bam -O $SCRATCH_JOB/${sample}.vcf.gz -dont-use-soft-clipped-bases --standard-min-confidence-threshold-for-calling 20 --dbsnp ${dbSNP_vcf}
rm $SCRATCH_JOB/tmp/*.*

echo "VariantFiltration"
gatk --java-options "-Djava.io.tmpdir=$SCRATCH_JOB/tmp" VariantFiltration --R ${ref_fasta} --V $SCRATCH_JOB/${sample}.vcf.gz --window 35 --cluster 3 --filter-name "FS" --filter "FS > 30.0" --filter-name "QD" --filter "QD < 2.0" -O /wynton/scratch/kkim/vcf_filtered/${sample}.vcf.gz
rm -f $SCRATCH_JOB/${sample}.vcf.gz*

echo "Move output files back to Leo" 
cp -f /wynton/scratch/userID/vcf_filtered/${sample}.vcf.* /labserver/home/userID/epic.neb/variantCall/indiv_vcf/

echo "Remove job-specific scratch folder"
cd /scratch/
rm -rf $SCRATCH_JOB/
rm -rf /wynton/scratch/userID/tmp/${sample}
rmdir /scratch/$USER/jobs/$JOB_ID
rmdir /scratch/$USER/jobs
echo "========== FINISHED JOBS for ${sample} =========="
date
