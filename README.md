   ![tests](https://github.com/AidanDavis47/airway-rnaseq-pipeline/actions/workflows/tests.yaml/badge.svg)

airway-rnaseq-pipeline

project for resume and to learn snakemake and get exprience with it, used claude for help on the project and learning snakemake.

A reproducible Snakemake pipeline that takes raw RNA-seq reads from the public airway dataset (GEO GSE52778) through quality control, transcript quantification and differential expression analysis, to find which genes respond to the corticosteroid dexamethasone in human airway smooth muscle cells.



How to run:
1. download code and reference
      git clone https://github.com/AidanDavis47/airway-rnaseq-pipeline.git
      cd airway-rnaseq-pipeline
      mkdir -p ref
      wget -p ref https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.transcripts.fa.gz

2. run with either conda or docker
   Conda:
      mamba env create -f environment.yml
      conda activate airway
      snakemake -n 
      snakemake --cores 4
   Docker:
      docker build -t airway-rnaseq .
      docker run --rm --user $(id -u):$(id -g) -e HOME=/tmp -v "$PWD":/pipeline airway-rnaseq snakemake --cores 4

3. Testing:
      pytest -v tests/


Folder and File Structure:
Snakefile (the pipeline rules)
config.yaml (settings)
samples.csv (sample sheet, has donors and treatments)
environment.yml (conda environment)
Dockerfile (for docker)
scripts/ (holds scripts for DESeq2, plots and creating test data for the test)
tests/ (holds the pytests and a small test data set)
figures/ (figures for the readme)


Goal of project:
What genes do Glucoccorticoids switch off and on in airway muscle cells?

Data is from GEO GSE52778, airway smooth muscle cells from 4 donors, each untreated and treated with dexamethasone (8 paired-end samples). The study's albuterol samples are not used.

Results:

Variation is caused by if the samples are treated or untreated.

PC1, which has 29% variance, seperates the treated and untreated samples cleanly. PC2, which has 19% variance, shows differences between donors. (see pca plot)

Dexamethason is found to mainly turn genes on. This is consistant with dexamethasone's mechanism of activating genes via the glucocorticoid receptor. (see volcano plot)


Refernce for data and comparing results:
Himes BE, et al. (2014). RNA-Seq transcriptome profiling identifies CRISPLD2 as a glucocorticoid responsive gene that modulates cytokine function in airway smooth muscle cells. PLoS ONE 9(6): e99625. https://doi.org/10.1371/journal.pone.0099625








-----NOTES MADE WHILE WORKING OF PROJECT-----


Key result

The top differentially expressed gene is CRISPLD2, which is upregulated about 5.9-fold after dexamethasone treatment (log2 fold change 2.56, adjusted p = 3.5 × 10⁻²⁰). This reproduces the main finding of the original study:

Himes BE et al. (2014). RNA-Seq Transcriptome Profiling Identifies CRISPLD2 as a Glucocorticoid Responsive Gene that Modulates Cytokine Function in Airway Smooth Muscle Cells. PLoS ONE. doi:10.1371/journal.pone.0099625

Pipeline overview
SRA reads ──► fastq-dump ──► FastQC ──────────────┐
                    │                             ├──► MultiQC report
                    └──► Salmon quant ────────────┘
                              ▲        │
GENCODE v48 transcripts ──► Salmon index
                    │                  ▼
                    └──► tx2gene ──► tximport + DESeq2 ──► results table



Rule	Tool	What it does
download	sra-tools (fastq-dump)	Downloads paired-end reads for each sample (first 1M reads by default)
fastqc	FastQC	Per-file read quality checks
salmon_index	Salmon	Builds a transcriptome index from GENCODE v48
salmon_quant	Salmon	Quantifies transcript expression for each sample
multiqc	MultiQC	Combines FastQC and Salmon QC into one HTML report
tx2gene	awk	Builds a transcript-to-gene map from the GENCODE FASTA headers
deseq2	tximport, DESeq2 (R)	Summarizes to gene level and tests dexamethasone vs. untreated
Experimental design

Eight samples from four donors. Each donor has one untreated and one dexamethasone-treated sample (see samples.csv).

The DESeq2 model is ~ donor + treatment. Including donor accounts for the paired design, so differences between donors don't mask the treatment effect. Genes are kept if they have at least 10 counts in at least 4 samples.

Quick start
1. Set up the environment
bash
conda env create -f environment.yml
conda activate airway
2. Download the reference transcriptome

The ref/ folder is not tracked in git, so download GENCODE v48 human transcripts first:

bash
mkdir -p ref
wget -P ref https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.transcripts.fa.gz
3. Run the pipeline
bash
snakemake -n            # dry run: show what will be executed
snakemake --cores 4     # run the pipeline
Configuration

Settings live in config.yaml:

Setting	Default	Description
samples	samples.csv	Sample sheet with sample, donor and treatment columns
reads_per_sample	1000000	Reads downloaded per sample (lower it for a faster test run)
transcripts	ref/gencode.v48.transcripts.fa.gz	Reference transcriptome
salmon_index	ref/salmon_index	Where the Salmon index is built
threads	4	Threads for Salmon
Outputs
Path	Contents
results/multiqc/multiqc_report.html	Combined QC report
results/fastqc/	Per-file FastQC reports
results/salmon/{sample}/quant.sf	Transcript-level quantification
results/deseq2/deseq2_results.csv	Gene-level results with log2 fold change, p-values and gene symbols, sorted by adjusted p-value
results/deseq2/dds.rds	Saved DESeq2 object for further analysis in R
logs/	Log file for each step



Himes BE et al. (2014), GEO accession GSE52778. Reference transcriptome from GENCODE release 48.



