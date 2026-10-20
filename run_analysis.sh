#!/bin/bash
# Stop the script if a command fails, an undefined variable is used, or a command in the pipeline fails.
set -euo pipefail

# 1. SETTING DIRECTORIES
DATA_DIR="/localdisk/data/BPSM/MyFirstPipeline"
OUTPUT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
METADATA="$DATA_DIR/fastq/Tco2.fqfiles"
# Directory containing the gene counts generated in Question 4
COUNT_DIR="$OUTPUT_DIR/gene_counts"
#directory  where mean expression and fold-change results are
RESULT_DIR="$OUTPUT_DIR/expression_results"

mkdir -p "$RESULT_DIR"

# QUESTION 5: Calculate mean gene counts for each experimental groups, groups are defined by SampleType, Time and Treatment.
# Replicate numbers are not part of the group name.
# AWK processes text files

# -v passes Bash variables into AWK.
# countdir and outdir are the AWK names for the directories
awk -v countdir="$COUNT_DIR" -v outdir="$RESULT_DIR" '
# Skip the first line because it contains column headings.
NR == 1 { next }
 # Column 1 contains the sample name, e.g. Tco935.
    sample = $1
    # Column 2 = SampleType
    # Column 4 = Time
    # Column 5 = Treatment
	#Replicate number is excluded, replicates belonging to the same experimental group are averaged

{
    sample = $1
    group = $2 "_" $4 "_" $5
    file = countdir "/" sample "_counts.tsv"

    # Read the gene counts for this sample
	#reads lines one at a time until there is none left
	#Column 4 contains the gene identifier 
	# Column 5 contains the gene description.
	# Column 6 contains the alignment count.
    # Adding 0 ensures AWK treats it as a number.
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
# It stops if the count file was missing or empty
    if (n == 0) {
        print "ERROR: Missing or empty count file: " file > "/dev/stderr"
        exit 1
    }
# Record the experimental group
    groups[group] = 1
# Count the number of biological replicates in this group
    replicates[group]++
    genes[gene] = 1

}
# Calculate the mean for every gene

END {
    for (group in groups) {
        output = outdir "/" group "_means.tsv"
        print "GeneID\tDescription\tMeanCount" > output
#Mean = total count / number of replicate
        for (gene in genes) {
            mean = sum[group, gene] / replicates[group]
            print gene "\t" description[gene] "\t" mean >> output
        }
        close(output)

  #Display which groups has been processed
        print "Calculated means for " group " (" replicates[group] " replicates)"
    }
}
' "$METADATA"

# QUESTION 6: Compare induced and uninduced groups at matching timepoints (24 and 48 hours) within each clone which helps investigate changes associated with RNAi
# A pseudocount of 1 avoids division by zero
# log2 fold change then allows increases and decreases to be compared
#first loop through the two RNAi clones
for CLONE in Clone1 Clone2
do
#Second loop through the two time point
# ${CLONE} and ${TIME} insert the variable values
    for TIME in 24 48
    do
        INDUCED="$RESULT_DIR/${CLONE}_${TIME}_Induced_means.tsv"
		# Matching uninduced control grou
        CONTROL="$RESULT_DIR/${CLONE}_${TIME}_Uninduced_means.tsv"
		# Name of the final fold-change 
        OUTPUT="$RESULT_DIR/${CLONE}_${TIME}_foldchange.tsv"

 # Check whether both input files exist
        # -f checks if a regular file exists
        # ! means NOT
        # || means OR
		#if either file is missing, skip this compariso
        if [ ! -f "$INDUCED" ] || [ ! -f "$CONTROL" ]; then
            echo "Skipping $CLONE at ${TIME}h: comparison group missing"
	#continue moves to the next loop iteration
            continue
        fi
# AWK to compare mean gene counts from the two files
        # -F '\t' tells AWK that the input columns are separated by tabs
        awk -F '\t' '
	#BEGIN runs before reading the files.
    #OFS sets the output field separator to a tab
	#FNR == NR identifies records in the first file
        BEGIN { OFS = "\t" }
        FNR == NR {

#Skip the header (FNR > 1).
# $1 = gene identifier.
# $3 = mean gene count.
            if (FNR > 1) control[$1] = $3
            next
        }
 # Skip its header and only compare genes present
        FNR > 1 && ($1 in control) {
  # Calculate fold-change ratio:
            # (induced mean + 1) / (uninduced mean + 1)
            ratio = ($3 + 1) / (control[$1] + 1)
			# Calculate log2 fold change
			#dividing by log(2) converts it to log base 2 beacuse AWK log() calculates the natural logarithm
            log2fc = log(ratio) / log(2)
 # Calculate absolute log2 fold change.
            # If log2fc is negative, multiply by -1.
            # allows increases and decreases to be ranked according to the size of the expression change.
            magnitude = (log2fc < 0 ? -log2fc : log2fc)
 # Output the gene ID, gene description, control mean, induced mean, fold-change ratio, log2 fold change and absolute log2 fold change
            print $1, $2, control[$1], $3, ratio, log2fc, magnitude
        }
 #Read the uninduced file first, then the induced file
        ' "$CONTROL" "$INDUCED" |

 # Pipe the AWK output into sort.
        # -t sets the field delimiter to a tab
        # -k7,7 sorts using column 7 (absolute log2 fold change)
        # n = numerical order
        # r = reverse order, largest firs
        sort -t "$(printf '\t')" -k7,7nr |
 # Add column headings to the final output, the previous AWK and sort commands produce data rows
        awk 'BEGIN {
            OFS = "\t"
            print "GeneID", "Description", "UninducedMean", "InducedMean", "FoldChangeRatio", "Log2FoldChange", "AbsoluteLog2FC"
        }
        { print }' > "$OUTPUT"

        echo "Created $OUTPUT"
    done
done
