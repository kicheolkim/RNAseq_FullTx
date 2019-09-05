#!/bin/bash
#
#$ -S /bin/bash
#$ -o /wynton/home/baranzinilab/kkim1/script/variant_call/log
#$ -e /wynton/home/baranzinilab/kkim1/script/variant_call/log
#$ -r y
#$ -j y
#$ -l mem_free=8G
#$ -l scratch=200G
#$ -l h_rt=300:00:00
#$ -R yes
#$ -pe smp 8
#$ -t 1-194

tasks=(0	13311d.list	22012b.list	38613b.list	43113b.list	43113c.list	43113d.list	43213b.list	43213c.list	46913b.list	47313a.list	47413a.list	47513a.list	47513b.list	47713b.list	48413b.list	48413c.list	48613a.list	48613b.list	48813b.list	48813c.list	48913b.list	48913c.list	49013a.list	49013b.list	49013c.list	49113a.list	49313a.list	49713a.list	49713b.list	50913b.list	51013a.list	51113a.list	51314a.list	51314b.list	51514a.list	51514b.list	51814a.list	51814b.list	51914a.list	52114a.list	52214a.list	52314a.list	52614a.list	52814a.list	52914a.list	53114b.list	55914a.list	56014a.list	56114a.list	56314a.list	56414a.list	56514a.list	56514b.list	56914a.list	56914b.list	57114a.list	57314a.list	57414a.list	57714a.list	57814a.list	57915a.list	58015a.list	58415a.list	58515a.list	58715a.list	58915a.list	59015a.list	59215a.list	59815a.list	59915a.list	60115a.list	60615a.list	60615b.list	60815a.list	61215a.list	61415a.list	61515a.list	61915a.list	62115a.list	62315a.list	62415a.list	63315a.list	63815a.list	64215a.list	64615a.list	65015a.list	65315a.list	65415a.list	65615a.list	65715a.list	65815a.list	66215a.list	66415a.list	66915a.list	67415a.list	67515a.list	67615a.list	67815a.list	68116a.list	68216a.list	69016a.list	69316a.list	69516a.list	69716a.list	70216a.list	70416a.list	70616a.list	70916a.list	71016a.list	71416a.list	71616a.list	71916a.list	72116a.list	72516a.list	73216a.list	73416a.list	73816a.list	73916a.list	74116a.list	74216a.list	75016a.list	75216a.list	75416a.list	75616a.list	75816a.list	76016a.list	76116a.list	76416a.list	77016a.list	77616a.list	78417a.list	78417b.list	78617a.list	78817a.list	78917a.list	79017a.list	79117a.list	79217a.list	79317a.list	79617a.list	79717a.list	80017a.list	80217a.list	80317a.list	80417a.list	80517a.list	80617a.list	80717a.list	80917a.list	81017a.list	81517a.list	81917a.list	82117a.list	82317a.list	82417a.list	82517a.list	82817a.list	82917a.list	83117a.list	83217a.list	83317a.list	83417a.list	83517a.list	84017a.list	84117a.list	84217a.list	84317a.list	84417a.list	84517a.list	84617a.list	84717a.list	84817a.list	84917a.list	85017a.list	85117a.list	85217a.list	85317a.list	85417a.list	85517a.list	85617a.list	85917a.list	86117a.list	86217a.list	86417a.list	86517a.list	86617a.list	86717a.list	87017a.list	87117a.list	87317a.list	87417a.list	87617a.list	87717a.list	87817a.list)
input="${tasks[$SGE_TASK_ID]}"

date
hostname
sample=${input%.*}
SCRATCH_JOB=/scratch/$USER/jobs/$JOB_ID/${sample}
mkdir -p /scratch/$USER/jobs/$JOB_ID
mkdir -p $SCRATCH_JOB
mkdir -p $SCRATCH_JOB/tmp
mkdir -p /wynton/scratch/kkim/tmp/${sample}
export _JAVA_OPTIONS="-Djava.io.tmpdir=$SCRATCH_JOB/tmp"
export TMP_DIR=$SCRATCH_JOB/tmp

ref_fasta=/wynton/home/baranzinilab/kkim1/genome/gencode.v28/GRCh38.primary_assembly.genome.fa
dbSNP_vcf=/wynton/home/baranzinilab/kkim1/genome/dbsnp_b151/00-All_fixed.vcf.gz

cd $SCRATCH_JOB/
pwd

echo "AddOrReplaceReadGroups"
java -jar /netapp/home/kkim1/bin/picard.jar AddOrReplaceReadGroups I=/wynton/scratch/kkim/merged_bam/${sample}_merged.bam O=$SCRATCH_JOB/${sample}.bam SO=coordinate RGLB=L001 RGPL=illumina RGPU=novaseq RGSM=${sample} TMP_DIR=$SCRATCH_JOB/tmp
rm $SCRATCH_JOB/tmp/*.*

echo "FixMateInformation"
java -jar /netapp/home/kkim1/bin/picard.jar FixMateInformation I=$SCRATCH_JOB/${sample}.bam O=$SCRATCH_JOB/${sample}_fixed.bam TMP_DIR=$SCRATCH_JOB/tmp
rm -f $SCRATCH_JOB/${sample}.bam
rm $SCRATCH_JOB/tmp/*.*

echo "MarkDuplicates"
java -jar /netapp/home/kkim1/bin/picard.jar MarkDuplicates I=$SCRATCH_JOB/${sample}_fixed.bam O=$SCRATCH_JOB/${sample}_dedupped.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=$SCRATCH_JOB/${sample}_dedupped.metrics TMP_DIR=$SCRATCH_JOB/tmp
rm -f $SCRATCH_JOB/${sample}_fixed.bam
rm $SCRATCH_JOB/tmp/*.*

echo "SplitNCigarReads"
java -Djava.io.tmpdir=$SCRATCH_JOB/tmp -jar /netapp/home/kkim1/bin/GenomeAnalysisTK-3.8-1-0/GenomeAnalysisTK.jar -T SplitNCigarReads -R ${ref_fasta} -I $SCRATCH_JOB/${sample}_dedupped.bam -o $SCRATCH_JOB/${sample}_split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
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
cp -f /wynton/scratch/kkim/vcf_filtered/${sample}.vcf.* /baranzini/leo/home/kicheol/epic.neb/variantCall/indiv_vcf/

echo "Remove job-specific scratch folder"
cd /scratch/
rm -rf $SCRATCH_JOB/
rm -rf /wynton/scratch/kkim/tmp/${sample}
rmdir /scratch/$USER/jobs/$JOB_ID
rmdir /scratch/$USER/jobs
echo "========== FINISHED JOBS for ${sample} =========="
date
