# RNA-Seq data analysis
Code for full transcript RNA-Seq data analysis of the paper titled "Cell-Type-Specific Transcriptomics Identifies Neddylation as a Novel Therapeutic Target in Multiple Sclerosis"

### RNA-Seq library
Sequencing library was prepared using NEBNext Ultra II Directional RNA Library Prep Kit for Illumina and NEBNext® rRNA Depletion Kit (Human/Mouse/Rat). Therefore, all options software used here are adapted for stranded RNA-seq reads using dUTP method. 

## Analysis workflow
#### 1. Adaptor trimming and low quality sequence: BBDuk of BBTools v38.05 (https://jgi.doe.gov/data-and-tools/bbtools/)
#### 2. QC of fastq file: FastQC v0.11.7
#### 3. Mapping to reference genome: STAR aligner v2.6.0c
- reference genome: GRCh38.p12 with Gencode annotation (release 28)
#### 4. Gene and transcript counting: RSEM v1.3.1
#### 5. Statistical analysis using R and Bioconductor 
- R version 3.5.1 and Bioconductor version 3.7

### Pre-processing


### 


-----
By Kicheol Kim (August 2019 updated)
