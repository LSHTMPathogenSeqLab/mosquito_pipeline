# ########### Cross-population extended haplotype homozygosity (XPEHH) ###########

# This script is designed to be run like a jupyter notebook so that you can see what is happening at each stage
# So I recommend running each cell individually and insepcting the output
# The following uses an example vcf and zarr file as an example, swap out for your vcf and adjust other parameters eg. chromosome names as needed.


## Compute the unstandardized cross-population extended haplotype homozygosity score (XPEHH) for each variant.
## allel.xpehh(h1, h2, pos, map_pos=None, min_ehh=0.05, include_edges=False, gap_scale=20000, max_gap=200000, is_accessible=None, use_threads=True)
# create h1 and h2, selecting all variants instead of segregating variants only, which is what we did in iHS

# %%
import os
import numpy as np
import allel
import zarr
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import gffutils

# %% set wd
os.chdir('/mnt/storage11/sophie/Anisa/wgs_project/')
os.getcwd()

# %%
# convert phased, filtered, VCF file to zarr file
# already converted to zarr
#allel.vcf_to_zarr('2019melasglobal_finalfiltered_gambiaealigned_phased.vcf.gz', '2019melasglobal_finalfiltered_gambiaealigned_phased.zarr', fields='*', overwrite=True)

# %%
callset = zarr.open('ghana_gambiae_filtered_phased.zarr', mode='r')
#callset.tree(expand=True)

# %%
## convert zarr file to genotype array
gt = allel.GenotypeDaskArray(callset['calldata/GT'])
print(gt.shape)

# %%
## import metadata
df_samples=pd.read_csv('anisa_metadata.csv',sep=',',usecols=['sample','site','phenotype','species'])
df_samples.head()
df_samples.groupby(by=['site']).count

# %%
## VCF is phased so we can convert genotype arrays made earlier to haplotype array

sample_ids = callset['samples'][:]

# %% Create arrays needed for Cameroon samples
# Get sample identifiers for Cameroon samples from df_samples
cam_sample_ids = df_samples[df_samples['site'] == 'Cameroon']['sample'].values
# Find indices of these samples in the genotype array
cam_indices = np.array([np.where(sample_ids == id)[0][0] for id in cam_sample_ids if id in sample_ids])
# Verify the indices are within the correct range
print("Max index:", cam_indices.max(), "Sample array size:", len(sample_ids))
# Select genotypes for Cameroon samples using the indices
gt_cam_samples = gt.take(cam_indices, axis=1)

# %% Create arrays needed for Bijagos samples
# Get sample identifiers for Bijagos samples from df_samples
bij_sample_ids = df_samples[df_samples['site'] == 'Guinea-Bissau']['sample'].values
# Find indices of these samples in the genotype array
bij_indices = np.array([np.where(sample_ids == id)[0][0] for id in bij_sample_ids if id in sample_ids])
# Verify the indices are within the correct range
print("Max index:", bij_indices.max(), "Sample array size:", len(sample_ids))
# Select genotypes for Cameroon samples using the indices
gt_bij_samples = gt.take(bij_indices, axis=1)

# %% these are from a phased VCF so we can convert the genotype arrays to haplotype arrays

h_array_cameroon = gt_cam_samples.to_haplotypes().compute()
h_array_cameroon

h_array_bijagos = gt_bij_samples.to_haplotypes().compute()
h_array_bijagos

# %% we need variant positions
pos = callset['variants/POS'][:]
chrom = callset['variants/CHROM'][:]

# %% compute xpehh
# xpehh_raw = allel.xpehh(h_sus, h_res, pos, map_pos=None, min_ehh=0.05, include_edges=False, gap_scale=20000, max_gap=20000, is_accessible=None, use_threads=True)
xpehh_raw = allel.xpehh(h_array_bijagos, h_array_cameroon, pos, use_threads=True)
xpehh_raw

# %% look for where the biggest signal is
xpehh_hit_max = np.nanargmax(xpehh_raw)
xpehh_hit_max

# %% genomic position of top hit
pos[xpehh_hit_max]

# %% Plot the raw xp-ehh values

fig, ax = plt.subplots()
ax.hist(xpehh_raw[~np.isnan(xpehh_raw)], bins=20)
ax.set_xlabel('Raw XP-EHH')
ax.set_ylabel('Frequency (no. variants)');

# %% standardise xpehh-raw

allele_counts_array = gt.count_alleles(max_allele=3).compute()
xpehh_std = allel.standardize_by_allele_count(xpehh_raw, allele_counts_array[:, 1])

# %% plot standardised xp-ehh values

fig, ax = plt.subplots()
ax.hist(xpehh_std[0][~np.isnan(xpehh_std[0])], bins=20)
ax.set_xlabel('Raw XP-EHH')
ax.set_ylabel('Frequency (no. variants)');

# %% Plot standardized xp-ehh

# define chromosome lengths and colours 
chromosome_lengths = {
    '2L': 49364325,
    '2R': 61545105,
    '3L': 41963435,
    '3R': 53200684,
    'anop_mito': 15363,
    'anop_X': 24393108
}

# Calculate cumulative offsets for each chromosome
cumulative_lengths = {}
cumulative_length = 0
for chrom, length in chromosome_lengths.items():
    cumulative_lengths[chrom] = cumulative_length
    cumulative_length += length

# %% Plot XP-EHH

# Set threshold
bijagos_threshold = 5
cameroon_threshold = -5

# Set up the plot
fig, ax = plt.subplots(figsize=(10, 6))

# Ensure that pos, chrom, and xpehh_std are all numpy arrays to support advanced indexing
pos = np.array(callset['variants/POS'][:])
chrom = np.array(callset['variants/CHROM'][:])
xpehh_standardised_values = np.array(xpehh_std[0])

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# Set up the plot
fig, ax = plt.subplots(figsize=(10, 6))

# Define colors for each chromosome (for illustration)
chromosome_colours = {
    '2L': '#3d348b', '2R': '#f18701', '3L': '#f7b801', '3R': '#7678ed', 'anop_mito': '#f35b04', 'anop_X': '#119DA4'
}

# Set threshold
bijagos_threshold = 5
cameroon_threshold = -5

# Create a list to hold the legend patches
legend_patches = []

# Filtered chromosomes list, assuming cumulative_lengths are defined for these
filtered_chroms = ['2L', '2R', '3L', '3R', 'anop_X', 'anop_mito']

# Iterate through each chromosome to plot its variants
for unique_chrom in filtered_chroms:
    chrom_mask = chrom == unique_chrom
    
    chrom_positions = pos[chrom_mask]
    chrom_xpehh_values = xpehh_standardised_values[chrom_mask]
    
    non_nan_mask = ~np.isnan(chrom_xpehh_values)
    chrom_positions_no_nan = chrom_positions[non_nan_mask]
    chrom_xpehh_values_no_nan = chrom_xpehh_values[non_nan_mask]
    
    adjusted_positions = chrom_positions_no_nan + cumulative_lengths.get(unique_chrom, 0)

    # Conditions for plotting
    solid_mask = (chrom_xpehh_values_no_nan >= bijagos_threshold) | (chrom_xpehh_values_no_nan <= cameroon_threshold)
    faded_mask = ~solid_mask
    
    # Plot solid points for values above 5 or below -5
    ax.scatter(adjusted_positions[solid_mask], 
               chrom_xpehh_values_no_nan[solid_mask], 
               color=chromosome_colours[unique_chrom], alpha=1.0, s=10)
    
    # Plot faded points for other values
    ax.scatter(adjusted_positions[faded_mask], 
               chrom_xpehh_values_no_nan[faded_mask], 
               color=chromosome_colours[unique_chrom], alpha=0.1, s=10)
    
    # Add patch for the legend
    patch = mpatches.Patch(color=chromosome_colours[unique_chrom], label=unique_chrom)
    legend_patches.append(patch)

# Add significance threshold lines and legend
ax.axhline(y=bijagos_threshold, color='black', linestyle='--', label='Bijagos Threshold')
ax.axhline(y=cameroon_threshold, color='black', linestyle='--', label='Cameroon Threshold')
ax.legend(handles=legend_patches, title='Chromosome', bbox_to_anchor=(1.05, 1), loc='upper left')

# Set labels
ax.set_xlabel('Genomic Position (bp)')
ax.set_ylabel('XP-EHH')
plt.tight_layout()
plt.savefig('xpehh_plot_600ppi.png', dpi=600)  # Save at 600 PPI

# %% list all positions with xpehh value over or below a certain threshold

bijagos_threshold_mask = xpehh_standardised_values >= bijagos_threshold
cameroon_threshold_mask = xpehh_standardised_values <= cameroon_threshold

# %% Apply the mask to filter the data
bij_significant_chrom = chrom[bijagos_threshold_mask]
bij_significant_pos = pos[bijagos_threshold_mask]
bij_significant_xpehh = xpehh_standardised_values[bijagos_threshold_mask]

cam_significant_chrom = chrom[cameroon_threshold_mask]
cam_significant_pos = pos[cameroon_threshold_mask]
cam_significant_xpehh = xpehh_standardised_values[cameroon_threshold_mask]

# %% Combine the filtered data into a structured array
bij_significant_xpehh_data = np.column_stack((bij_significant_chrom, bij_significant_pos, bij_significant_xpehh))
cam_significant_xpehh_data = np.column_stack((cam_significant_chrom, cam_significant_pos, cam_significant_xpehh))

# %% Convert the structured array to a pandas DataFrame for easier handling
df_significant_bij_xpehh = pd.DataFrame(bij_significant_xpehh_data, columns=['Chromosome', 'Position', 'XPEHH'])
df_significant_cam_xpehh = pd.DataFrame(cam_significant_xpehh_data, columns = ['Chromosome', 'Position', 'XPEHH'])

# %% Save to csv
df_significant_bij_xpehh.to_csv(f'df_significant_bij_xpehh_bijagos_threshold_{bijagos_threshold}.csv', index=False)
df_significant_cam_xpehh.to_csv(f'df_significant_cam_xpehh_cameroon_threshold_{cameroon_threshold}.csv', index=False)

# %% bring in the gff file to understand where each of these variants is

print("Using GFF file to bring in annotations for these positions")

## Annotate the Bijagos XPEHH file
# Parameters
input_file_name = f"df_significant_bij_xpehh_bijagos_threshold_{bijagos_threshold}.csv"
output_file_name = f"df_significant_bij_xpehh_bijagos_threshold_{bijagos_threshold}_annotated.csv"
gff_file = '/mnt/storage11/sophie/reference_genomes/A_gam_P4_ensembl/Anopheles_gambiae.AgamP4.56.chr.gff3'

# Function to find and format the GFF line(s) that overlap a given position
def find_overlapping_gff_lines(chromosome, position, gff_file):
    overlapping_lines = []
    with open(gff_file, 'r') as gff:
        for line in gff:
            if line.startswith('#') or line.strip() == "":
                continue  # Skip header and empty lines
            parts = line.split('\t')
            if parts[0] == chromosome and int(parts[3]) <= position <= int(parts[4]):
                formatted_line = ",".join(parts).strip()
                overlapping_lines.append(formatted_line)
    return overlapping_lines

# Open the output file to write the annotated positions
with open(output_file_name, "w") as outfile:
    # Write the header line
    outfile.write("Chromosome,Position,XPEHH,Gff_Annotation\n")

    # Open the file containing significant iHS positions to read
    with open(input_file_name, "r") as infile:
        next(infile) #skip header line
        for line in infile:
            parts = line.strip().split(",")
            chromosome, position, XPEHH = parts[0], int(parts[1]), parts[2]
            
            # Find overlapping GFF lines for the position
            overlapping_gff_lines = find_overlapping_gff_lines(chromosome, position, gff_file)
            
            # Join all overlapping GFF lines into a single string
            gff_annotation = "; ".join(overlapping_gff_lines)
            
            # Write to the output file
            outfile.write(f"{chromosome},{position},{XPEHH},{gff_annotation}\n")

print(f"XP-EHH significant values identified and GFF annotations written here: {output_file_name}")
# %%
## Annotate the Cameroon XPEHH file
# Parameters
input_file_name = f"df_significant_cam_xpehh_cameroon_threshold_{cameroon_threshold}.csv"
output_file_name = f"df_significant_cam_xpehh_cameroon_threshold_{cameroon_threshold}_annotated.csv"
gff_file = '/mnt/storage11/sophie/reference_genomes/A_gam_P4_ensembl/Anopheles_gambiae.AgamP4.56.chr.gff3'

# Function to find and format the GFF line(s) that overlap a given position
def find_overlapping_gff_lines(chromosome, position, gff_file):
    overlapping_lines = []
    with open(gff_file, 'r') as gff:
        for line in gff:
            if line.startswith('#') or line.strip() == "":
                continue  # Skip header and empty lines
            parts = line.split('\t')
            if parts[0] == chromosome and int(parts[3]) <= position <= int(parts[4]):
                formatted_line = ",".join(parts).strip()
                overlapping_lines.append(formatted_line)
    return overlapping_lines

# Open the output file to write the annotated positions
with open(output_file_name, "w") as outfile:
    # Write the header line
    outfile.write("Chromosome,Position,XPEHH,Gff_Annotation\n")

    # Open the file containing significant iHS positions to read
    with open(input_file_name, "r") as infile:
        next(infile) #skip header line
        for line in infile:
            parts = line.strip().split(",")
            chromosome, position, XPEHH = parts[0], int(parts[1]), parts[2]
            
            # Find overlapping GFF lines for the position
            overlapping_gff_lines = find_overlapping_gff_lines(chromosome, position, gff_file)
            
            # Join all overlapping GFF lines into a single string
            gff_annotation = "; ".join(overlapping_gff_lines)
            
            # Write to the output file
            outfile.write(f"{chromosome},{position},{XPEHH},{gff_annotation}\n")

print(f"XP-EHH significant values identified and GFF annotations written here: {output_file_name}")
# %%
