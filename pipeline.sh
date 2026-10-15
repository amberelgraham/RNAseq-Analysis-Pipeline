#!/bin/bash
#telling linux to use bash


#MyFirstPipeline
#RNA-seq analysis pipeline for trypanosoma congoLense


#1.Setting directories 
#Location of assignment
DATA_DIR="/localdisk/data/BPSM/MyFirstPipeline"
#compressed FASTQ sequencing files
FASTQ_DIR="$DATA_DIR/fastq"
#where pipeline is stored $HOME  directory
OUTPUT_DIR="$HOME/MyFirstPipeline"

#DATA_DIR this creates variable so we dont have to keep repeating the /localdisk/data/BPSM/MyFirstPipeline ie giving it a shortter name
 


#2.Preparing the reference genome 
# Create directory for reference genome
#-p prevents error if the directory already exists
mkdir -p "$OUTPUT_DIR/genome"



# Only prepare the reference genome if the Bowtie2 index does not already exist
#if is from BASH  [ ! -f FILE] "if file does not exist,  found from help if    
if [ ! -f "$OUTPUT_DIR/genome/Tcongo_index.1.bt2" ]
then

        # Copy compressed genome into project 
        cp "$DATA_DIR/Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz" "$OUTPUT_DIR/genome/"

        # Decompress reference genome while keeping compressed copy
		#-d decmpress
		#-k keeps original can find via gzip --help
        gzip -dk "$OUTPUT_DIR/genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz"

        # Build Bowtie2 index from reference genome
		#INPUT the reference .fasta
		#OUTPUT is Tcongo_index files 
		#found via bowtie2-build --help, AL suggetsion on assignment  good manual;  https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer
        bowtie2-build "$OUTPUT_DIR/genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta" \
                      "$OUTPUT_DIR/genome/Tcongo_index"

fi


#3.FAST Q - Qualiry check RAW sequencing reads 

#Creating the directory for the fastq results, the -p used again
 
mkdir -p "$OUTPUT_DIR/fastq_results"

#running fastqc on all compressed FASTQ files 
#*.fq.gz selects every file ending in .fq.gz in the FASTQ directory such asTco-230_fq.gz etc, cretaes a loop(do,done) 
# The for loop takes each file one at a time and stores its filename in the variable FILE.
# -o tells FastQC where to save the output.
for FILE in "$FASTQ_DIR"/*.fq.gz
do 
	fastqc "$FILE" -o "$OUTPUT_DIR/fastq_results"
done 





#4. Assesing FASTQ results 
#FASTQC produces a zip file for each FASTQ file 
#summary.txt inside containing PASS/warn/FAIL results for each
# unzip -p = print a file from inside the zip without extracting it
# grep -E 'FAIL|WARN' = keep FAIL or WARN lines
# cut -f1,2 = keep fields 1 and 2
# sort = put matching lines together
# uniq -c = count identical lines
#note read puts text into variables 

#I tested manually before adding to piepline and the results were 

# 92 FAIL  Per base sequence content
# 34 WARN  Overrepresented sequences
# 2 WARN   Per base sequence content
# 2 WARN   Per sequence GC content

for ZIP in "$OUTPUT_DIR"/fastq_results/*_fastqc.zip
do
        unzip -p "$ZIP" '*/summary.txt'
done | grep -E 'FAIL|WARN' | cut -f1,2 | sort | uniq -c

# Check the number of sequences in the FASTQ files.
# fastqc_data.txt contains "Total Sequences".
# grep finds that line in each FastQC report.
# sort + uniq -c count how many files had each result.
# Manual result: 96 Total Sequences 100000
# Therefore all 96 FASTQ files contained 100000 sequences.

for ZIP in "$OUTPUT_DIR"/fastq_results/*_fastqc.zip
do
        unzip -p "$ZIP" '*/fastqc_data.txt' | grep "Total Sequences"
done | sort | uniq -c

#5. prep for allignment
# Create directory for alignment results
mkdir -p "$OUTPUT_DIR/alignment_results"


#6.Reading Samples 

# Tco2.fqfiles contains the sample information and the names of the two FASTQ files belonging to each sample.
#
# tail -n +2 starts at line 2,  skips the column headings.
# output from tail becomes input for the while loop.
# read splits each line into variables: SAMPLE TYPE REP TIME TREATMENT END1 END2

# Found this through  help read


#For example 
#$SAMPLE = Tco230
#$END1   = Tco-230_1.fq.gz
#$END2   = Tco-230_2.fq.gz


# Read each sample pair from the metadata file,  END1 and END2 are the paired reads (R1 and R2).

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR/alignment_results"

#echo prints sample being processed
tail -n +2 "$FASTQ_DIR/Tco2.fqfiles" | while read SAMPLE TYPE REP TIME TREATMENT END1 END2
do
# Skip samples already aligned
if [ -s "$OUTPUT_DIR/alignment_results/${SAMPLE}_sorted.bam" ] &&
   [ -s "$OUTPUT_DIR/alignment_results/${SAMPLE}_sorted.bam.bai" ]; then
    echo "Skipping $SAMPLE"
    continue
fi

        echo "Processing $SAMPLE"

#7. Align paired reads with bowtie2 against teh T. congolense reference genome
        bowtie2 -x "$OUTPUT_DIR/genome/Tcongo_index" \
                -1 "$FASTQ_DIR/$END1" \
                -2 "$FASTQ_DIR/$END2" \
                -S "$OUTPUT_DIR/alignment_results/$SAMPLE.sam"

   # -x =  Bowtie2 index, use to bowtie2 --help
        # -1 = first paired FASTQ / R1
        # -2 = second paired FASTQ / R2
        # -S = output SAM file

   #manually tested this using Tco230.
        # Tco230 produced a 98.65% overall alignment rate.


#8. Convert SAM to BAM
# SAM = text alignment format
# Stop processing this sample if Bowtie2 failed
if [ $? -ne 0 ]; then
    echo "Alignment failed for $SAMPLE"
    continue
fi

        # BAM = binary alignment format
	#I had to search this up, use samtool and then samtools view --help
        samtools view -b \
                -o "$OUTPUT_DIR/alignment_results/$SAMPLE.bam" \
                "$OUTPUT_DIR/alignment_results/$SAMPLE.sam"

        # Sort BAM by genomic position 
        samtools sort \
                -o "$OUTPUT_DIR/alignment_results/${SAMPLE}_sorted.bam" \
                "$OUTPUT_DIR/alignment_results/$SAMPLE.bam"

        # Index sorted BAM, accesing genomic regions
        samtools index "$OUTPUT_DIR/alignment_results/${SAMPLE}_sorted.bam"

# Remove intermediate files after successful indexing
if [ -s "$OUTPUT_DIR/alignment_results/${SAMPLE}_sorted.bam.bai" ]; then
    rm -f "$OUTPUT_DIR/alignment_results/$SAMPLE.sam"
    rm -f "$OUTPUT_DIR/alignment_results/$SAMPLE.bam"
fi

done

# QUESTION 4 - Count alignments overlapping each annotated gene

BED_FILE="$DATA_DIR/TriTrypDB-46_TcongolenseIL3000_2019.bed"
COUNT_DIR="$OUTPUT_DIR/gene_counts"

mkdir -p "$COUNT_DIR"

for BAM in "$OUTPUT_DIR"/alignment_results/*_sorted.bam
do
    SAMPLE=$(basename "$BAM" _sorted.bam)

    echo "Counting genes for $SAMPLE"

    bedtools coverage \
        -a "$BED_FILE" \
        -b "$BAM" \
        -counts > "$COUNT_DIR/${SAMPLE}_counts.tsv"
done

# QUESTIONS 5 AND 6
# Calculate group mean expression and fold changes.
# Run the analysis script from the same directory as this pipeline.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

bash "$SCRIPT_DIR/run_analysis.sh"

# QUESTIONS 5 AND 6
# Calculate average gene expression and fold changes.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

bash "$SCRIPT_DIR/run_analysis.sh"

