#!/usr/bin/env python3

import os
import sys

# Check if both input files are provided
if len(sys.argv) != 3:
    print("Usage: {} <unique.txt> <count.txt>".format(sys.argv[0]))
    sys.exit(1)

# Assign input file paths to variables
UNIQUE_FILE = sys.argv[1]
COUNT_FILE = sys.argv[2]

# Make output directories
OUTDIR1 = "final_counts"
os.makedirs(OUTDIR1, exist_ok=True)
OUTDIR2 = "error_counts"
os.makedirs(OUTDIR2, exist_ok=True)

# Check if input files exist
if not (os.path.isfile(UNIQUE_FILE) and os.path.isfile(COUNT_FILE)):
    print("Error: Input file not found")
    sys.exit(1)

# Output file names
output_filename = os.path.join(OUTDIR1, COUNT_FILE.replace('_all_count.txt', '_FINAL_count.txt'))
error_filename = os.path.join(OUTDIR2, COUNT_FILE.replace('_all_count.txt', '_ERROR_count.txt'))

# Loop through each line of count.txt
with open(COUNT_FILE, 'r') as count_file:
    for line in count_file:
        sequence, count = line.strip().split('\t')
        found = False
        
        # Loop through each line of unique.txt
        with open(UNIQUE_FILE, 'r') as unique_file:
            for unique_line in unique_file:
                number, identifier, seq = unique_line.strip().split('\t')
                # Check if the sequence in count.txt exists within any of the sequences in unique.txt
                if sequence.lower() == seq.lower():
                    # Write to output file
                    with open(output_filename, 'a') as output_file:
                        output_file.write("{}\t{}\t{}\n".format(identifier, sequence, count))
                    found = True
                    break  # Sequence found, no need to continue searching
        
        # If sequence not found, write to error_counts.txt
        if not found:
            with open(error_filename, 'a') as error_file:
                error_file.write("{}\t{}\n".format(sequence, count))

print("Final counts written to", output_filename)
print("Unmatched sequences written to", error_filename)

#make script executable with chmod +x script.py
#run by using ./script.py /path/to/unique.txt /path/to/count.txt