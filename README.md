# RNA-Seq data analysis
Code for full transcript RNA-Seq data analysis of the paper titled "Cell-Type-Specific Transcriptomics Identifies Neddylation as a Novel Therapeutic Target in Multiple Sclerosis"

### RNA-Seq library
Sequencing library was prepared using NEBNext Ultra II Directional RNA Library Prep Kit for Illumina and NEBNext® rRNA Depletion Kit (Human/Mouse/Rat). Therefore, all options software used here are adapted for stranded RNA-seq reads using dUTP method.

## Analysis workflow
BASH scripts are prepared for running in the UCSF Wynton cluster (https://ucsf-hpc.github.io/wynton/).


#### 1. Adaptor trimming and low quality sequence: BBDuk of BBTools v38.05 [Wynton-cluster/rawQC-trimming-mapping-counting.sh]
#### 2. QC of fastq file: FastQC v0.11.7 [Wynton-cluster/rawQC-trimming-mapping-counting.sh]
#### 3. Mapping to reference genome: STAR aligner v2.6.0c [Wynton-cluster/rawQC-trimming-mapping-counting.sh]
- reference genome: GRCh38.p12 with Gencode annotation (release 28)
#### 4. Gene and transcript counting: RSEM v1.3.1 [Wynton-cluster/rawQC-trimming-mapping-counting.sh]
#### 5. Statistical analysis using R and Bioconductor [R]
- R v3.5.1 and Bioconductor v3.7
#### 6. Weighted gene co-expression network analysis (WGCNA) [R/WGCNA_MS-HC.R]
- WGCNA v1.64.1 (R package)
#### 7. Variant calling: GATK4, vcftools, bcftools
- BAM files of different cell subset from same subject were merged [Wynton-cluster/VarCall_2_mergeBAM.sh]
- Variant calling procedure is adapted from GATK Best Practices for variant calling on RNAseq [Wynton-cluster/VarCall_1_filter-calling.sh]
  (https://software.broadinstitute.org/gatk/best-practices/workflow?id=11164)
#### 8. Cis-eQTL analysis: FastQTL v2.0 [Wynton-cluster/ and R/]


-----
By Kicheol Kim (August 2019 updated)
