# RNA-Seq data analysis
Code for full transcript RNA-Seq data analysis of the paper titled "Cell-Type-Specific Transcriptomics Identifies Neddylation as a Novel Therapeutic Target in Multiple Sclerosis"

### RNA-Seq library
Sequencing library was prepared using NEBNext Ultra II Directional RNA Library Prep Kit for Illumina and NEBNext® rRNA Depletion Kit (Human/Mouse/Rat). Therefore, all options software used here are adapted for stranded RNA-seq reads using dUTP method.

## Analysis workflow
BASH scripts are prepared for running in the UCSF Wynton clster (https://ucsf-hpc.github.io/wynton/).


#### 1. Adaptor trimming and low quality sequence: BBDuk of BBTools v38.05 [rawQC-trimming-mapping-counting.sh]
#### 2. QC of fastq file: FastQC v0.11.7 [rawQC-trimming-mapping-counting.sh]
#### 3. Mapping to reference genome: STAR aligner v2.6.0c [rawQC-trimming-mapping-counting.sh]
- reference genome: GRCh38.p12 with Gencode annotation (release 28)
#### 4. Gene and transcript counting: RSEM v1.3.1 [rawQC-trimming-mapping-counting.sh]
#### 5. Statistical analysis using R and Bioconductor 
- R v3.5.1 and Bioconductor v3.7
#### 6. Weighted gene co-expression network analysis (WGCNA)
- WGCNA v1.64.1 (R package)
#### 7. Variant calling: GATK4, vcftools, bcftools
#### 8. Cis-eQTL analysis: FastQTL v2.0

### 


-----
By Kicheol Kim (August 2019 updated)
