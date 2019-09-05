#!/bin/bash
#
#$ -S /bin/bash
#$ -o /wynton/scratch/userID/qtl/log_fastqtl
#$ -e /wynton/scratch/userID/qtl/log_fastqtl
#$ -r y
#$ -j y
#$ -l mem_free=6G
#$ -l h_rt=120:00:00
#$ -pe smp 8
#$ -R yes
#$ -t 1-9
##### test eQTL for MS-all patients with GRCh38 - variants from merged bam file

tasks1=(0	CD4_expr_MS-untreat.bed.gz	CD8_expr_MS-untreat.bed.gz	CD14_expr_MS-untreat.bed.gz)
tasks2=(0	covar_CD4_MS-untreat.txt	covar_CD8_MS-untreat.txt	covar_CD14_MS-untreat.txt)
tasks3=(0	CD4_TrtNaive_filtered_final.recode.vcf.gz	CD8_TrtNaive_filtered_final.recode.vcf.gz	CD14_TrtNaive_filtered_final.recode.vcf.gz)
input1="${tasks1[$SGE_TASK_ID]}"
input2="${tasks2[$SGE_TASK_ID]}"
input3="${tasks3[$SGE_TASK_ID]}"

date
hostname
## 0. Create job-specific scratch folder
echo ${input1}
name=${input1%%.*}
cell=${input1%*_*_*}

date
cd /wynton/scratch/userID/qtl/results
mkdir -p /wynton/scratch/userID/qtl/results/
mkdir -p /wynton/scratch/userID/qtl/results/tmp/${name}
mkdir -p /wynton/scratch/userID/qtl/results/tmp/${name}_nominal

date
echo "START FastQTL-permutation test run..."
for c in $(seq 1 100); do
	/netapp/home/userID1/bin/FastQTL/bin/fastQTL.static --vcf /wynton/scratch/userID/qtl/variant_rnaseq/${input3} --bed /wynton/scratch/userID/qtl/pheno_expr/${input1} --cov /wynton/scratch/userID/qtl/pheno_expr/${input2} --permute 10000 100000 --out /wynton/scratch/userID/qtl/results/tmp/${name}/${name}_chunk$c\.results --log /wynton/scratch/userID/qtl/results/tmp/${name}/${name}_chunk$c\.log --chunk $c 100
	done
cat /wynton/scratch/userID/qtl/results/tmp/${name}/${name}_chunk*.results > /wynton/scratch/userID/qtl/results/${name}_permute.results
cat /wynton/scratch/userID/qtl/results/tmp/${name}/${name}_chunk*.log > /wynton/scratch/userID/qtl/results/${name}_permute.log

echo "START FastQTL-nominal test run..."
for c in $(seq 1 100); do
	/netapp/home/userID1/bin/FastQTL/bin/fastQTL.static --vcf /wynton/scratch/userID/qtl/variant_rnaseq/${input3} --bed /wynton/scratch/userID/qtl/pheno_expr/${input1} --cov /wynton/scratch/userID/qtl/pheno_expr/${input2} --out /wynton/scratch/userID/qtl/results/tmp/${name}_nominal/${name}_chunk$c\.results --log /wynton/scratch/userID/qtl/results/tmp/${name}_nominal/${name}_chunk$c\.log --chunk $c 100
	done
cat /wynton/scratch/userID/qtl/results/tmp/${name}_nominal/${name}_chunk*.results > /wynton/scratch/userID/qtl/results/${name}_nominal.results
cat /wynton/scratch/userID/qtl/results/tmp/${name}_nominal/${name}_chunk*.log > /wynton/scratch/userID/qtl/results/${name}_nominal.log

echo "END FastQTL run..." 

echo "Remove job-specific scratch folder"
cd /scratch/
rmdir /scratch/$USER/jobs/$JOB_ID
rmdir /scratch/$USER/jobs
echo "========== REMOVED AND FINISHED JOBS ========="
date
