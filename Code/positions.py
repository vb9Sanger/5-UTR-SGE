#!/usr/bin/env python3

import os
import sys

# Check if both input files are provided
if len(sys.argv) != 2:
    print("Usage: {} <final_count.txt>".format(sys.argv[0]))
    sys.exit(1)

# Assign input file paths to variables
FINAL_COUNT_FILE = sys.argv[1]

# Make output directories
OUTDIR1 = "positions"
os.makedirs(OUTDIR1, exist_ok=True)

# Check if input files exist
if not os.path.isfile(FINAL_COUNT_FILE):
    print("Error: Input file not found")
    sys.exit(1)
# Output file names
output_filename = os.path.join(
    OUTDIR1,
    os.path.basename(FINAL_COUNT_FILE).replace(
        "_FINAL_count.txt", "_FINAL_count_positions.txt"
    ),

    
)


def process_final_counts(input_file, output_filename):
    lines = []  # List to store lines from the input file
    with open(input_file, "r") as fin:
        for line in fin:
            parts = line.strip().split("\t")
            identifier = parts[0]
            sequence = parts[1]
            count = parts[2]
            position = get_position(identifier)
            lines.append(
                (identifier, sequence, count, position)
            )  # Add line to the list

    # Sort the lines based on the fourth column (position)
    lines.sort(key=lambda x: x[3])

    # Write the sorted lines to the output file
    with open(output_filename, "w") as fout:
        for line in lines:
            fout.write("{}\t{}\t{}\t{}\n".format(line[0], line[1], line[2], line[3]))


def get_position(identifier):
    # Split the identifier by underscores
    parts = identifier.split("_")
    if len(parts) >= 4:  # Check if there is a number after the third underscore
        try:
            position = int(parts[3])  # Try to convert the fourth part to an integer
            return position
        except ValueError:
            pass  # If conversion to integer fails, continue to check for number after the second underscore
    if len(parts) >= 3:  # Check if there is a number after the second underscore
        try:
            position = int(parts[2])  # Try to convert the third part to an integer
            return position
        except ValueError:
            pass  # If conversion to integer fails, continue to return 0
    return (
        0  # Default position if no number is found after the second or third underscore
    )


def main():
    process_final_counts(FINAL_COUNT_FILE, output_filename)
    print("Positions written to", output_filename)


if __name__ == "__main__":
    main()
