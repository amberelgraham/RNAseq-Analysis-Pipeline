DATA_DIR="/localdisk/data/BPSM/MyFirstPipeline"
#compressed FASTQ sequencing files
FASTQ_DIR="$DATA_DIR/fastq"
#where pipeline is stored $HOME  directory
OUTPUT_DIR="$HOME/MyFirstPipeline"
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


