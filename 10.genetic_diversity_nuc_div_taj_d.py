# This script is designed to be run like a jupyter notebook so that you can see what is happening at each stage
# So I recommend running each cell individually and insepcting the output
# The following uses my melas vcf and zarr file as an example, swap out for your vcf and adjust other parameters eg. chromosome names as needed.


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

# %%
callset = zarr.open('bijagos_only_melas_phased.zarr', mode='r')
#callset.tree(expand=True)

# %%
## convert zarr file to genotype array
gt = allel.GenotypeDaskArray(callset['calldata/GT'])
print(gt.shape)

# %%
ac = gt.count_alleles()
# %%
pos = callset['variants/POS'][:]

# %% sort pos and reorder ac
# get the indices that would sort the pos array
sorted_indices = np.argsort(pos)
# sort pos using those indices
pos_sorted = pos[sorted_indices]
# sort ac according to the same indices
ac_sorted = ac[sorted_indices]

# %% Compute nucleotide diversity
pi = allel.sequence_diversity(pos_sorted,ac_sorted)
pi

# %% Computing for only chromosome 3L.
chrom = callset['variants/CHROM'][:]
# 1. Filter for chromosome 3L first
chrom_mask = chrom == '3L'
pos_3L = pos[chrom_mask]
ac_3L = gt.count_alleles()

# %% 2. Then, sort pos_3L and reorder ac_3L accordingly
sorted_indices_3L = np.argsort(pos_3L)
pos_sorted_3L = pos_3L[sorted_indices_3L]
ac_sorted_3L = ac_3L[sorted_indices_3L]

# %% 3. Now, you can compute nucleotide diversity for chromosome 3L
pi_3L = allel.sequence_diversity(pos_sorted_3L, ac_sorted_3L)
pi_3L

# %% Now calculate in windows over the genome

pi_windowed_3L = allel.windowed_diversity(pos_sorted_3L,ac_sorted_3L,size=20000)
pi_windowed_3L
# this returns four arrays, the first is nucleotide diversity in each window, the second is the windows used
# rhe third is the number of accessible bases in each window, # the fourt is the number of variants in each window
# this is done over genomic position, the number of variants depends on the genomic window

# %% Plot nucleotide diversity (pi) and calculate mean and SD 

# Extracting nucleotide diversity values and window positions
pi_values = pi_windowed_3L[0]
windows = pi_windowed_3L[1]

# Calculating midpoints of windows for plotting
midpoints = np.mean(windows, axis=1)

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(midpoints, pi_values, color='black', marker='.', linestyle='-', linewidth=1.25, markersize=1)
plt.xlabel('Position on Chromosome 3L')
plt.ylabel('Nucleotide Diversity (π)')
plt.title('Nucleotide Diversity across Chromosome 3L')
plt.grid(True)
plt.show()

# Calculating mean and standard deviation of nucleotide diversity
mean_pi = np.mean(pi_values)
std_pi = np.std(pi_values)

print(f"Mean nucleotide diversity (π): {mean_pi}")
print(f"Standard deviation of nucleotide diversity (π): {std_pi}")

# %% Calculate tajimas d for the same region

D, windows, counts = allel.windowed_tajima_d(pos_sorted_3L, ac_sorted_3L, size=20000)

# this outputs three arrays, the first is Tajima's D, the second is the windows, the third is the number of variants in each window

# %% Plot Tajima's D over the chromosome

# Extracting nucleotide diversity values and window positions
# Filtering out NaN values from the Tajima's D calculations
valid_indices = ~np.isnan(D)  # ~np.isnan(D) creates a boolean mask where NaN values are False
filtered_td_values = D[valid_indices]
filtered_midpoints = np.mean(windows[valid_indices], axis=1)

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(filtered_midpoints, filtered_td_values, color='black', marker='.', linestyle='-', linewidth=1.25, markersize=1)
plt.xlabel('Position on Chromosome 3L')
plt.ylabel('Tajima\'s D')
plt.title('Tajima\'s D across Chromosome 3L')
plt.grid(True)
plt.show()

# Calculating mean and standard deviation of Tajima's D, excluding NaN values
mean_td = np.mean(filtered_td_values)
std_td = np.std(filtered_td_values)

print(f"Mean Tajima's D: {mean_td}")
print(f"Standard deviation of Tajima's D: {std_td}")


# %% Calculate mean and standard deviation of Tajima's D for star t of chromosome
# First 1,000,000 bp, so first 50 windows of 20,000

# Slice D array to get first 50 values, excluding any NaN values
D_first_50_windows = D[:50
valid_indices_50 = ~np.isnan(D_first_50_windows)
filtered_td_values_50 = D_first_50_windows[valid_indices_50]

# Calcualte mean and SD
mean_td_first_50 = np.mean(filtered_td_values_50)
std_td_first_50 = np.std(filtered_td_values_50)

mean_td_first_50, std_td_first_50

# %%
