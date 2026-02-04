
######################## SELECTION STATISTICS #########################


# This script is designed to be run like a jupyter notebook so that you can see what is happening at each stage
# So I recommend running each cell individually and insepcting the output
# The following uses my melas vcf and zarr file as an example, swap out for your vcf and adjust other parameters eg. chromosome names as needed.

# some selection tests only support biallelic variants, not multiallelic. 
# This filtered VCF should already be biallelic SNPs only.
# Note that the zarr file should have been made from a phased vcf
# You make the zarr file with allel.vcf_to_zarr('phased_vcf_file_name.vcf.gz', 'output_name.zarr', fields='*', overwrite=True)

# %%
import os
os.chdir('/mnt/storage11/sophie/bijagos_mosq_wgs/2019_melas_fq2vcf_gambiae_aligned/genomics_database_melas2019plusglobal/genomics_database_melas2019plusglobal_vcf/melas_2019_plusglobal_filtering')
os.getcwd()

# %%
import numpy as np
import allel
import zarr
import pandas as pd
import sys
import matplotlib.pyplot as plt
import seaborn as sns
import random

## convert phased, filtered, VCF file to zarr file
# %%
# allel.vcf_to_zarr('2022gambiaevcfphased.vcf.gz', '2022gambiaevcfphased.zarr', fields='*', overwrite=True)

# %%
callset = zarr.open('2019melasglobal_finalfiltered_gambiaealigned_phased.zarr', mode='r')
#callset.tree(expand=True)

# %%
## convert zarr file to genotype array (GenotypeDaskArray)
genotype_array = allel.GenotypeDaskArray(callset['calldata/GT'])
print(genotype_array.shape)

# %%
## import metadata
df_samples=pd.read_csv('metadata_melasplusglobal.csv',sep=',',usecols=['sample','year','country','species','island'])
df_samples.head()
df_samples.groupby(by=['country']).count

# %% Filter genotype array
# 1. Make allele count array so that we can filter for biallelic and segregating variants
ac_array = genotype_array.count_alleles(max_allele=8).compute() # ac_array is an AlleleCountsArray

# %% Filter for biallelic and segregating and store as boolean array
ac_array_filter = ac_array.is_biallelic_01()

# %% remove variants that are on Y_unplaced chromosome
chrom = callset['variants/CHROM'][:]
exclude_chrom = 'Y_unplaced'
ac_array_filter = ac_array_filter & (chrom != exclude_chrom)

# %% filter the allele counts array using the filter
filtered_ac_array = ac_array.compress(ac_array_filter, axis=0)

# %% also make a genotype dask array
filtered_gt = genotype_array.compress(ac_array_filter, axis = 0)

# %% partition samples by county using metadata, split by index value and store as array
bissau_samples = df_samples[df_samples['country'] == 'Guinea-Bissau'].index.values
cameroon_samples = df_samples[df_samples['country'] == 'Cameroon'].index.values
gambia_samples = df_samples[df_samples['country'] == 'The Gambia'].index.values

# %% select genotypes for variants and store as a genotype dask array
# this array contains many columns (1 for each mosquito), many million rows for each variant, and stores the genotype for each.
gt_bissau_samples = filtered_gt.take(bissau_samples, axis=1)
gt_cameroon_samples = filtered_gt.take(cameroon_samples, axis=1)
gt_gambia_samples = filtered_gt.take(gambia_samples, axis=1)

# %% convert this genotype array to haplotype array (we can do this because the original data was phased)
# the haplotype array is similar to the genotype array, but there are two columns per mosquito, one for each haplotype
h_bissau_seg = gt_bissau_samples.to_haplotypes().compute()
h_cameroon_seg = gt_cameroon_samples.to_haplotypes().compute()
h_gambia_seg = gt_gambia_samples.to_haplotypes().compute()

# %% also store chromosome of each variant as we need this for shading the plots later
chrom = callset['variants/CHROM'][:]
chrom_filtered = chrom.compress(ac_array_filter, axis=0)

# %% we need variant positions
pos = callset['variants/POS'][:]
pos_filtered = pos.compress(ac_array_filter, axis=0)

# %% some variants in 1000 genomes project have multiple variants at the same genomic position, 
# which causes problems for some selection tests in scikit-allel. 
# Let's check if there any of these.
count_multiple_variants = np.count_nonzero(np.diff(pos_filtered== 0))

if count_multiple_variants == 0:
    print("No cases where there are multiple variants at the same genomic position, script will continue")
else:
    print("There are multiple variants at the same genomic position. This causes problems with some selection tests using sci-kit allel.")
    sys.exit()  # This will stop the script. If you want the script to continue anyway, # out this line

# %% H12 was calculated using phased biallelic SNPs in 1000 bp windows along the genome
# SNP windows, using the garuds_h function in scikit-allel.
# Working with Bissau samples only.

# ######## Create Iterations ##########

# Your list of tuples - there are two haplotypes per mosquito, for bissau melas there are 
# 30 mosquitoes so 60 haplotypes in the haplotype array h_bissau_seg
tuples_list = [
    (0, 1), (2, 3), (4, 5), (6, 7), (8, 9), (10, 11), (12, 13), (14, 15),
    (16, 17), (18, 19), (20, 21), (22, 23), (24, 25), (26, 27), (28, 29),
    (30, 31), (32, 33), (34, 35), (36, 37), (38, 39), (40, 41), (42, 43),
    (44, 45), (46, 47), (48, 49), (50, 51), (52, 53), (54, 55), (56, 57),
    (58, 59),
]

# Number of iterations
n_iterations = 200

# Initialize a list to store h12 values for each window across all iterations
iterated_bissau_h12_values = []

# Loop for n_iterations
for _ in range(n_iterations):
    # Randomly select n-1 tuples
    selected_tuples = random.sample(tuples_list, 29)

    # Extract the column indices from the tuples
    column_indices = [idx for tup in selected_tuples for idx in tup]

    # Subset the haplotype array
    subset_bissau_hap_array = h_bissau_seg[:, column_indices]

    # Calculate h12 for the subsetted hap array
    _, iterated_real_bissau_h12, _, _ = allel.moving_garud_h(subset_bissau_hap_array, 1000)

    # Store the h12 value for each window in the list
    iterated_bissau_h12_values.append(iterated_real_bissau_h12)

# Convert the list to a numpy array
iterated_bissau_h12_values_array = np.array(iterated_bissau_h12_values)

# Calculate the mean of the h12 values for each window
mean_bissau_h12_per_window = np.mean(iterated_bissau_h12_values_array, axis=0)

# mean_h12_per_window now contains the mean h12 value for each window

# %% Check the number of windows for each array 
num_windows_bissau = mean_bissau_h12_per_window.shape[0]
print("Number of windows for bissau samples:", num_windows_bissau)

# %% The h12 values are calculated in windows of 1000 SNPs. Each SNP has a POS value which is in the POS array.
window_size = 1000  # Define the window size as 1000 SNPs
num_windows = len(pos_filtered) // window_size  # Calculate the number of windows in total

# %% Plot. Calculate the median genomic position for each window
# for each iteration (i) a segment of pos_res_seg is processed
median_positions = [np.median(pos_filtered[i * window_size: (i + 1) * window_size]) for i in range(num_windows)]

bissau_h12_chrom = {}
# Loop through each window
for i in range(num_windows):
    chrom = chrom_filtered[i * window_size]  # Assumes the chromosome for all SNPs in a single window is consistent 
    pos = median_positions[i]
    h12 = mean_bissau_h12_per_window[i]
    # Add this data to the corresponding chromosome in the dictionary
    if chrom not in bissau_h12_chrom:
        bissau_h12_chrom[chrom] = {'positions': [], 'h12': []}
    bissau_h12_chrom[chrom]['positions'].append(pos)
    bissau_h12_chrom[chrom]['h12'].append(h12)
# Now plot for each chromosome
for chrom in bissau_h12_chrom:
    plt.figure(figsize=(10, 6))
    plt.scatter(bissau_h12_chrom[chrom]['positions'], bissau_h12_chrom[chrom]['h12'], alpha=0.6)
    plt.xlabel('Median position of SNP windows across the chromosome')
    plt.ylabel('Mean H12 value in Bijagos Archipelago Anopheles melas samples')
    plt.title(f'Genome Wide Selection Scan of H12 Values across Chromosome {chrom}')
    filename = f'H12_value_200_iterations_{chrom}.png'
    plt.savefig(filename)
    plt.show()
# Save the Bissau samples dictionary which contains chromsome, positions and h12 value.
import csv
with open('bissau_h12_chrom.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Chromosome', 'Position', 'H12'])
    # Iterate through the dictionary and write data
    for chrom, data in bissau_h12_chrom.items():
        for position, h12 in zip(data['positions'], data['h12']):
            writer.writerow([chrom, position, h12])

#############################################################
### Visually inspect for peaks and select a cutoff point ###
### No peak identified so not using the below code, but if there were peaks could use a cut off (here 0.2) ###
#############################################################

# %% Extract the h-values above 0.2 with their corresponding SNP window's median genomic position
# Open a CSV file to write
with open('iterated_mean_h12_above_02_data.csv', 'w', newline='') as file:
    writer = csv.writer(file)

    # Write the header
    writer.writerow(['Chromosome', 'Position', 'mean_h12_Bissau'])

    # Iterate through the dictionary and write data for H12 values above 0.2
    for chrom, data in bissau_h12_chrom.items():
        for pos, h12 in zip(data['positions'], data['h12']):
            # Check if H12 value is above 0.2 for either resistant or not resistant
            if h12 > 0.2:
                writer.writerow([chrom, pos, h12, h12])

# %% Make a range for each peak that you want to look at genes beneath
# Write code to automate this section
# This section should create the iterated_chromosome_ranges.csv that is used in the next section

# Loop through each line of the h12_peaks_locations.txt file
while IFS=',' read -r chr start end; do
    echo "Processing: $chr from $start to $end" # Debugging line
    # Use awk to filter genes within the window from the GFF file
    awk -v chr="$chr" -v start="$start" -v end="$end" '
        $1 == chr && $4 <= end && $5 >= start {print $0}' \
        /mnt/storage11/sophie/reference_genomes/A_gam_P4_ensembl/Anopheles_gambiae.AgamP4.56.chr.gff3 >> iterated_h12_peaks_genes.txt
done < iterated_chromosome_ranges.csv
