#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_fasta> <output_fasta>"
    exit 1
fi

input="$1"
output="$2"


# Count the number of sequences and the length of the first sequence
seq_count=$(grep -c "^>" "$input")
seq_length=$(grep -v "^>" "$input" | head -1 | tr -d '\n' | wc -c)

# Start writing to the PHYLIP file
echo "$seq_count $seq_length" > "$output"

# Process each sequence and ensure proper formatting
awk '
    /^>/ {
        # Print the previous sequence if one exists
        if (seq) {
            printf("%-10s  %s\n", id, seq) >> output_file
        }
        # Extract the sequence ID (header without ">")
        id = substr($0, 2)
        seq = ""
    }
    /^[^>]/ {
        # Append sequence lines
        seq = seq $0
    }
    END {
        # Print the last sequence
        if (seq) {
            printf("%-10s  %s\n", id, seq) >> output_file
        }
    }
' output_file="$output" "$input"

echo "Conversion complete: PHYLIP file written to $output"

