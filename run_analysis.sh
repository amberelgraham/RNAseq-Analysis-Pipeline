#!/bin/bash
set -euo pipefail

DATA_DIR="/localdisk/data/BPSM/MyFirstPipeline"
OUTPUT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
METADATA="$DATA_DIR/fastq/Tco2.fqfiles"
COUNT_DIR="$OUTPUT_DIR/gene_counts"
RESULT_DIR="$OUTPUT_DIR/expression_results"

mkdir -p "$RESULT_DIR"

# QUESTION 5: Calculate mean gene counts for each experimental grou# Groups are defined by SampleType, Time and Treatment.
# Replicate numbers are not part of the group name.

awk -v countdir="$COUNT_DIR" -v outdir="$RESULT_DIR" '
NR == 1 { next }

{
    sample = $1
    group = $2 "_" $4 "_" $5
    file = countdir "/" sample "_counts.tsv"

    # Read the gene counts for this sample.
    n = 0
        while ((getline line < file) > 0) {
        split(line, fields, "\t")
        gene = fields[4]
	genes[gene] = 1
        description[gene] = fields[5]
        count = fields[6] + 0

        sum[group, gene] += count
        n++
    }
    close(file)

    if (n == 0) {
        print "ERROR: Missing or empty count file: " file > "/dev/stderr"
        exit 1
    }

    groups[group] = 1
    replicates[group]++
    genes[gene] = 1

}
END {
    for (group in groups) {
        output = outdir "/" group "_means.tsv"
        print "GeneID\tDescription\tMeanCount" > output

        for (gene in genes) {
            mean = sum[group, gene] / replicates[group]
            print gene "\t" description[gene] "\t" mean >> output
        }
        close(output)
        print "Calculated means for " group " (" replicates[group] " replicates)"
    }
}
' "$METADATA"

# QUESTION 6: Compare induced and uninduced groups
# at matching timepoints within each clone.
# A pseudocount of 1 avoids division by zero.
# Signed log2 fold change represents increases and decreases.

for CLONE in Clone1 Clone2
do
    for TIME in 24 48
    do
        INDUCED="$RESULT_DIR/${CLONE}_${TIME}_Induced_means.tsv"
        CONTROL="$RESULT_DIR/${CLONE}_${TIME}_Uninduced_means.tsv"
        OUTPUT="$RESULT_DIR/${CLONE}_${TIME}_foldchange.tsv"

        if [ ! -f "$INDUCED" ] || [ ! -f "$CONTROL" ]; then
            echo "Skipping $CLONE at ${TIME}h: comparison group missing"
            continue
        fi

        awk -F '\t' '
        BEGIN { OFS = "\t" }
        FNR == NR {
            if (FNR > 1) control[$1] = $3
            next
        }
        FNR > 1 && ($1 in control) {
            ratio = ($3 + 1) / (control[$1] + 1)
            log2fc = log(ratio) / log(2)
            magnitude = (log2fc < 0 ? -log2fc : log2fc)

            print $1, $2, control[$1], $3, ratio, log2fc, magnitude
        }
        ' "$CONTROL" "$INDUCED" |
        sort -t "$(printf '\t')" -k7,7nr |
        awk 'BEGIN {
            OFS = "\t"
            print "GeneID", "Description", "UninducedMean", "InducedMean", "FoldChangeRatio", "Log2FoldChange", "AbsoluteLog2FC"
        }
        { print }' > "$OUTPUT"

        echo "Created $OUTPUT"
    done
done
