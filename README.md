# RNA-seq Analysis Pipeline | *Trypanosoma congolense*

## Project Overview

This project develops a command-line RNA sequencing (RNA-seq) analysis pipeline using Bash and AWK to investigate gene expression changes following RNA interference (RNAi) in *Trypanosoma congolense*.


The pipeline processes paired-end sequencing data, performs quality control, aligns reads to a reference genome, generates gene-level counts, calculates mean expression across biological replicates, and produces fold-change comparisons.

The workflow was developed as part of bioinformatics (msc) at the University of Edinburgh.

## Biological Background

*Trypanosoma congolense* is a protozoan parasite responsible for animal African trypanosomiasis (nagana), a disease affecting livestock and agricultural productivity in sub-Saharan Africa.

T. congolense is a single-celled parasite transmitted by infected tsetse flies. It causes animal African trypanosomiasis (nagana), which can lead to weight loss, weakness and death

The experiment investigates the effects of knocking down a gene hypothesised to be involved in energy metabolism.

Three sample types were investigated:

Wild-type (WT) parasites
Clone1 containing an inducible RNAi construct
Clone2 containing the same RNAi construct

Samples were collected at 0, 24 and 48 hours, with induced and uninduced experimental conditions.

The objective was to identify potential gene expression changes associated with RNAi induction and investigate possible compensatory biological responses.

## Analysis Workflow

**1. Quality control FastQC**

Assess raw paired-end FASTQ sequencing data and generate quality reports.

**2. Sequencing quality assessment**

Examine FastQC outputs to assess sequencing quality and read characteristics.

**3. Reference genome alignment Bowtie2 and SAMtools**

Build a reference genome index, align paired-end reads and convert alignment outputs into sorted, indexed BAM files.

**4. Gene-level counting BEDTools**

Calculate the number of alignments overlapping annotated genomic regions using the *T. congolense* IL3000 gene annotation BED file.

**5. Mean gene expression AWK**

Use sample metadata to group biological replicates by sample type, time point and treatment, then calculate mean gene counts for each experimental group.

**6. Fold-change analysis AWK and Unix utilities**

Compare induced and uninduced experimental groups at matched time points, calculating expression ratios and log₂ fold changes.

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

https://docs.github.com/en/get-started/writing-on-github/working-with-advanced-formatting/organizing-information-with-tables
## Dataset

The supplied dataset contains 48 paired-end RNA-seq samples, each subsampled to 100,000 read pairs.

Reference genome: *Trypanosoma congolense* IL3000.

Gene annotations: TriTrypDB 2019 BED annotation, containing 10,312 annotated genomic regions.

Raw sequencing data and large generated alignment files are not included in this repository.

## Running the Pipeline

The workflow was developed in University of Edinburgh Linux server environment.

Software includes FastQC, Bowtie2, SAMtools, BEDTools, Bash and AWK.

The main pipeline script is:

`pipeline.sh`

The downstream analysis script is:

`run_analysis.sh`


## Outputs

The workflow generates:
FastQC reports showing the quality of the sequencing data
Sorted and indexed BAM alignment files
Tables showing the number of sequencing alignments overlapping each annotated gene
Mean gene count tables for experimental groups
Ranked fold-change tables

For the supplied dataset, the analysis produced 48 gene-count files, 15 group-mean files and four matched induced-versus-uninduced fold-change comparisons.

## Limitations
Unnormalised counts: Gene counts have not been adjusted for differences in sequencing depth between samples.
Fold changes show differences between experimental groups but do not determine whether changes are statistically significant.
Low-count genes. Adding 1 to the counts prevents division by zero but can affect fold changes when counts are low.
Read counting, BEDTools counts alignments overlap gene regions, which may not represent unique RNA fragments.
The analysis treats genes as continuous regions without introns, following the teaching exercise.
Statistical tools such as DESeq2 could be used to assess which genes show stronger evidence of differential expression.

## Future Improvements

Future improvements could include adjusting gene counts to allow fairer comparisons between samples, using statistical methods to identify changes in gene expression, and investigating which biological pathways may be affected.

The pipeline could also be improved to detect errors more effectively and automatically handle additional experimental groups.


