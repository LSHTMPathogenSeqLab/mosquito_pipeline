## This script uses the ratio of coverage between the X and 3R chromosomes in Anopheles to determine sex of the mosquito

import subprocess
import pandas as pd

# Function to read BAM file names from a list file
def read_bam_list(filename):
    with open(filename, 'r') as file:
        bamfiles = file.readlines()
        bamfiles = [line.strip() for line in bamfiles]
    return bamfiles

# Function to calculate average coverage using samtools depth
def calculate_average_coverage(bamfile, chromosome):
    command = f"samtools depth -r {chromosome} {bamfile} | awk '{{sum+=$3; count++}} END {{if (count > 0) print sum/count; else print 0}}'"
    process = subprocess.run(command, shell=True, capture_output=True, text=True)
    average_coverage = float(process.stdout.strip()) if process.stdout.strip() else 0
    return average_coverage

# Read BAM files from bamlist.txt
bamfiles = read_bam_list('bamlist.txt')

# Initialize an empty list to store results
results = []

# Process each BAM file
for bamfile in bamfiles:
    # Calculate average coverage for X and 3R chromosomes
    anop_X_coverage = calculate_average_coverage(bamfile, 'X')
    r3_coverage = calculate_average_coverage(bamfile, '3R')

    # Calculate coverage ratio of X to 3R to prevent division by zero
    coverage_ratio = anop_X_coverage / r3_coverage if r3_coverage > 0 else 0

    # Determine sex based on the ratio (assuming criteria or provide your own)
    # NOTE: Adjust the criteria as necessary
    if coverage_ratio > 1:  # Example criterion, adjust based on actual use case
        sex = 'Female'
    else:
        sex = 'Male'

    # Append results, including X and 3R coverages, their ratio, and the sex call
    results.append([bamfile, anop_X_coverage, r3_coverage, coverage_ratio, sex])

# Convert results to a DataFrame and save to a file
df = pd.DataFrame(results, columns=['BamFile', 'XCoverage', 'R3Coverage', 'CoverageRatio', 'Sex'])
df.to_csv('coverage_and_sex_from_bam_X_3R.txt', index=False, sep='\t')

print('Finished processing. Results saved to coverage_and_sex_from_bam_X_3R.txt.')