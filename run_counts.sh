#!/bin/bash

# Question 4: Count reads overlapping annotated genes

DATA_DIR="/localdisk/data/BPSM/MyFirstPipeline"
OUTPUT_DIR="$HOME/MyFirstPipeline"

BED_FILE="$DATA_DIR/TriTrypDB-46_TcongolenseIL3000_2019.bed"
COUNT_DIR="$OUTPUT_DIR/gene_counts"

mkdir -p "$COUNT_DIR"

# Process each sorted BAM file
for BAM in "$OUTPUT_DIR"/alignment_results/*_sorted.bam
do
    SAMPLE=$(basename "$BAM" _sorted.bam)

    echo "Counting genes for $SAMPLE"

    bedtools coverage \
        -a "$BED_FILE" \
        -b "$BAM" \
        -counts > "$COUNT_DIR/${SAMPLE}_counts.tsv"

done
