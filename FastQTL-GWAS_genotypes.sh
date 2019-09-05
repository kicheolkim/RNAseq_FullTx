#!/bin/bash
#
#$ -S /bin/bash
#$ -o /scrapp/userID/qtl/log_fastqtl
#$ -e /scrapp/userID/qtl/log_fastqtl
#$ -r y
#$ -j y
#$ -l mem_free=5G
#$ -l h_rt=60:00:00
#$ -pe smp 8
#$ -R yes
#$ -t 1-3
##### test eQTL for healthy with GRCh38
tasks1=(0	gene_expr_filt_CD4.txt.gz	gene_expr_filt_CD8.txt.gz	gene_expr_filt_CD14.txt.gz)
tasks2=(0	gene_covar_CD4.txt	gene_covar_CD8.txt	gene_covar_CD14.txt)
input1="${tasks1[$SGE_TASK_ID]}"
input2="${tasks2[$SGE_TASK_ID]}"

date
hostname
## 0. Create job-specific scratch folder
echo ${input1}
tmp=${input1#*_*_*_*}
cell=${tmp%%.*}

date
cd /scrapp/userID/qtl/results_expr3
mkdir -p /scrapp/userID/qtl/results_expr3/tmp
mkdir -p /scrapp/userID/qtl/results_expr3/tmp/${cell}_nominal

date
echo "START FastQTL-permutation test run..."
for c in $(seq 1 100); do
	/netapp/home/userID/bin/FastQTL/bin/fastQTL.static --vcf /scrapp/userID/qtl/GWAS_mschip/UCSF-MSchip_29overlap_remapGRCh38_rehead_sorted_filt.recode.vcf.gz --bed /scrapp/userID/qtl/pheno_expr/with_GWAS_MSchip/${input1} --cov /scrapp/userID/qtl/pheno_expr/with_GWAS_MSchip/${input2} --permute 1000 100000 --out /scrapp/userID/qtl/results_expr3/tmp/${cell}_chunk$c\.results --log /scrapp/userID/qtl/results_expr3/tmp/${cell}_chunk$c\.log --chunk $c 100
	done
cat /scrapp/userID/qtl/results_expr3/tmp/${cell}_chunk*.results > /scrapp/userID/qtl/results_expr3/${cell}_expr_permute.results
cat /scrapp/userID/qtl/results_expr3/tmp/${cell}_chunk*.log > /scrapp/userID/qtl/results_expr3/${cell}_expr_permute.log

echo "START FastQTL-nominal test run..."
for c in $(seq 1 100); do
	/netapp/home/userID/bin/FastQTL/bin/fastQTL.static --vcf /scrapp/userID/qtl/GWAS_mschip/UCSF-MSchip_29overlap_remapGRCh38_rehead_sorted_filt.recode.vcf.gz --bed /scrapp/userID/qtl/pheno_expr/with_GWAS_MSchip/${input1} --cov /scrapp/userID/qtl/pheno_expr/with_GWAS_MSchip/${input2} --out /scrapp/userID/qtl/results_expr3/tmp/${cell}_nominal/${cell}_chunk$c\.results --log /scrapp/userID/qtl/results_expr3/tmp/${cell}_nominal/${cell}_chunk$c\.log --chunk $c 100
	done
cat /scrapp/userID/qtl/results_expr3/tmp/${cell}_nominal/${cell}_chunk*.results > /scrapp/userID/qtl/results_expr3/${cell}_expr_nominal.results
cat /scrapp/userID/qtl/results_expr3/tmp/${cell}_nominal/${cell}_chunk*.log > /scrapp/userID/qtl/results_expr3/${cell}_expr_nominal.log

echo "END FastQTL run..." 

echo "Remove job-specific scratch folder"
cd /scratch/
rmdir /scratch/$USER/jobs/$JOB_ID
rmdir /scratch/$USER/jobs
echo "========== REMOVED AND FINISHED JOBS ========="
date
