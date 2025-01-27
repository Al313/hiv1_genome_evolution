#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_fasta> <output_fasta>"
    exit 1
fi

# Assign input and output file names from command-line arguments

input="$1"
output="$2"


# Process the interleaved FASTA
awk '
    /^>/ {
        # Print the header and reset the sequence accumulator
        if (seq) { print seq >> output_file }
        print $0 >> output_file
        seq = ""
    }
    /^[^>]/ {
        # Append sequence lines
        seq = seq $0
    }
    END {
        # Print the last sequence
        if (seq) { print seq >> output_file }
    }
' output_file="$output" "$input"

echo "Conversion complete: Output saved to $output"

