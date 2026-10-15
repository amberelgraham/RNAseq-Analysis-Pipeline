# RNA-seq Analysis Pipeline | *Trypanosoma congolense*

## Project Overview

This project develops a command-line RNA sequencing (RNA-seq) analysis pipeline using Bash and AWK to investigate gene expression changes following RNA interference (RNAi) in *Trypanosoma congolense*.

The pipeline processes paired-end sequencing data, performs quality control, aligns reads to a reference genome, generates gene-level counts, calculates mean expression across biological replicates, and produces descriptive fold-change comparisons.

The workflow was developed as part of bioinformatics programming training at the University of Edinburgh and is presented here as a demonstration of computational biology, data processing, automation and reproducible analysis skills.

## Biological Background

*Trypanosoma congolense* is a protozoan parasite responsible for animal African trypanosomiasis (nagana), a disease affecting livestock and agricultural productivity in sub-Saharan Africa.

The experiment investigates the effects of knocking down a gene hypothesised to be involved in energy metabolism.

Three sample types were investigated:

- Wild-type (WT) parasites
- Clone1 containing an inducible RNAi construct
- Clone2 containing the same RNAi construct

Samples were collected at 0, 24 and 48 hours, with induced and uninduced experimental conditions.

The objective is to identify potential gene expression changes associated with RNAi induction and investigate possible compensatory biological responses.

## Analysis Workflow

**1. Quality control — FastQC**

Assess raw paired-end FASTQ sequencing data and generate quality reports.

**2. Sequencing quality assessment**

Examine FastQC outputs to assess sequencing quality and read characteristics.

**3. Reference genome alignment — Bowtie2 and SAMtools**

Build a reference genome index, align paired-end reads and convert alignment outputs into sorted, indexed BAM files.

**4. Gene-level counting — BEDTools**

Calculate the number of alignments overlapping annotated genomic regions using the *T. congolense* IL3000 gene annotation BED file.

**5. Mean gene expression — AWK**

Use sample metadata to group biological replicates by sample type, time point and treatment, then calculate mean gene counts for each experimental group.

**6. Fold-change analysis — AWK and Unix utilities**

Compare induced and uninduced experimental groups at matched time points, calculating descriptive expression ratios and log₂ fold changes.

Results are ranked by absolute log₂ fold-change magnitude and include gene identifiers and functional descriptions.

## Technologies and Skills

| Technology | Application |
|---|---|
| Bash | Pipeline automation, loops, variables and file management |
| AWK | Data extraction, grouping and numerical calculations |
| FastQC | Sequencing quality control |
| Bowtie2 | Reference genome alignment |
| SAMtools | BAM conversion, sorting and indexing |
| BEDTools | Gene-level alignment counting |
| Git and GitHub | Version control and project management |
| Linux / Unix | Command-line data processing |

## Dataset

The supplied teaching dataset contains 48 paired-end RNA-seq samples, each subsampled to 100,000 read pairs.

Reference genome: *Trypanosoma congolense* IL3000.

Gene annotations: TriTrypDB 2019 BED annotation, containing 10,312 annotated genomic regions.

Raw sequencing data and large generated alignment files are not included in this repository.

## Running the Pipeline

The workflow was developed for the University of Edinburgh MSc bioinformatics Linux server environment.

Required software includes FastQC, Bowtie2, SAMtools, BEDTools, Bash and AWK.

The main pipeline script is:

`pipeline.sh`

The downstream analysis script is:

`run_analysis.sh`

The scripts expect access to the supplied teaching dataset and its sample metadata file. They are not currently configured as a standalone workflow for arbitrary external datasets.

## Outputs

The workflow generates:

- Sequencing quality-control reports
- Sorted and indexed BAM alignment files
- Gene-level count tables
- Mean gene-count tables for experimental groups
- Ranked descriptive fold-change tables

For the supplied dataset, the analysis produced 48 gene-count files, 15 group-mean files and four matched induced-versus-uninduced fold-change comparisons.

## Limitations

This pipeline is intended as an educational and exploratory analysis workflow rather than a validated differential gene expression analysis.

- Gene counts are not normalised for sequencing depth or other technical factors.
- Fold changes are descriptive and do not establish statistical significance.
- A pseudocount of 1 is used when calculating expression ratios to avoid division by zero.
- BEDTools overlap counting does not necessarily represent unique paired-end fragments.
- The analysis assumes genes are continuous genomic regions without introns, as specified in the teaching exercise.
- Further analysis using tools such as DESeq2 would be required for statistically supported differential expression conclusions.

## Future Improvements

Potential extensions include sequencing-depth normalisation, statistical differential expression analysis, functional enrichment analysis, improved error handling and more flexible experimental comparisons.

## Project Purpose

This project demonstrates the application of programming and computational analysis to a biological research question, combining experimental design, genomic data processing, automation, reproducibility and interpretation of gene expression data

