#!/usr/bin/env python3

import os
import sys

# Check if both input files are provided
if len(sys.argv) != 2:
    print("Usage: {} <count.txt>".format(sys.argv[0]))
    sys.exit(1)

# Assign input file paths to variables
COUNT_FILE = sys.argv[1]

# Make output directories
OUTDIR1 = "lengths"
os.makedirs(OUTDIR1, exist_ok=True)

# Check if input files exist
if not os.path.isfile(COUNT_FILE):
    print("Error: Input file not found")
    sys.exit(1)

# Output file names
output_filename = os.path.join(OUTDIR1, os.path.basename(COUNT_FILE).replace('_FINAL_count.txt', '_length_count.txt'))

def count_sequences(COUNT_FILE):
    lengths = [0] * 8  # To store counts for sequences of length <= 50, 100, 150, ..., 350
    with open(COUNT_FILE, 'r') as file:
        for line in file:
            identifier, sequence, count = line.strip().split('\t')
            sequence_length = len(sequence)
            count = int(count)  # Convert count from string to integer
            # Multiply count by sequence length and update counts based on sequence length
            if sequence_length <= 50:
                lengths[0] += count
            elif 50 < sequence_length <= 100:
                lengths[1] += count
            elif 100 < sequence_length <= 150:
                lengths[2] += count
            elif 150 < sequence_length <= 200:
                lengths[3] += count
            elif 200 < sequence_length <= 250:
                lengths[4] += count
            elif 250 < sequence_length <= 300:
                lengths[5] += count
            elif sequence_length <= 350:
                lengths[6] += count
            else:
                lengths[7] += count
    return lengths

def write_lengths(output_filename, lengths):
    with open(output_filename, 'w') as file:
        file.write("Length <= 50: {}\n".format(lengths[0]))
        file.write("Length <= 100: {}\n".format(lengths[1]))
        file.write("Length <= 150: {}\n".format(lengths[2]))
        file.write("Length <= 200: {}\n".format(lengths[3]))
        file.write("Length <= 250: {}\n".format(lengths[4]))
        file.write("Length <= 300: {}\n".format(lengths[5]))
        file.write("Length <= 350: {}\n".format(lengths[6]))
        file.write("Length > 350: {}\n".format(lengths[7]))

def main():
    lengths = count_sequences(COUNT_FILE)
    write_lengths(output_filename, lengths)
    print("Sequence lengths written to", output_filename)

if __name__ == "__main__":
    main()

#make script executable with chmod +x script.py
#run by using ./script.py /path/to/count.txt
    