# FILTERING OF COMBINED GENOTYPED VCF
# Below are a set of suggested hard and soft filters to filter your genomics database so that only high quality variants remain.
# I'd advise that you do this step by step and view your vcf as you go with 'zless' or 'less' and inspect your database
# The following uses the examples of my VCF database 2019_merged_melas.vcf.gz
# The hard filtering options are mostly from recommended GATK filtering parameters
# More information on the filtering parameters:
    # GATK filter recommendations: https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
    # excellent explanation of this hard filtering and additional soft filtering
    # https://evomics.org/learning/population-and-speciation-genomics/2020-population-and-speciation-genomics/first-steps-in-genomic-data-analysis/#ex3.4


########### FILTER 1: remove INDELS (insertions and deletions) with bcftools to keep just biallelic SNPs

bcftools query -l 2019_merged_melas.vcf.gz

bcftools view -M2 -m2 -v snps 2019_merged_melas.vcf.gz -Oz -o bi_snps_2019_merged_melas.vcf.gz

# tabix index the compressed VCF file, creates .vcf.gz.tbi
tabix -p vcf bi_snps_2019_merged_melas.vcf.gz

### Note that to retain biallelic SNPs and split these into multiallelic sites, you need to decompose:

bcftools norm -m -any input.vcf.gz -Oz -o output.vcf.gz

# filter out the contigs from the VCF file to keep just the chromosomes, note that to produce a properly bgzipped vcf file you need the -Oz flag
# replace the below chromosomes with the chromosome names for the species you are analysing

bcftools view bi_snps_2019_merged_melas.vcf.gz --regions 2L,2R,3L,3R,Mt,X,Y_unplaced | bcftools sort -Oz -o bi_snps_chr_2019_merged_melas.vcf.gzz

# view the unique chromosome names to check that this has worked properly
bcftools query -f '%CHROM\n' bi_snps_chr_gambiae_nov2022.2023_07_05.genotyped.vcf.gz | sort | uniq > unique_chromosomes_filtered.txt
tabix -p vcf bi_snps_chr_gambiae_nov2022.2023_07_05.genotyped.vcf.gz

########### FILTER 2: Filter samples to keep those with 40% of genome with > 10x coverage.
# and min-ac=1 so that all variants that remain are still variants after sample removal
# This is another opportunity to remove any samples that you think are too low quality if they are in your VCF
# For example, using the basic_statistics_analysis_wgs you can identify samples with 40% of the genome with over 10x coverage.
# Then create the file samples_40_10.txt containing the sample names that you want to keep.

bcftools view -S samples_40_10.txt bi_snps_2019_merged_melas.vcf.gz --min-ac=1 -Oz -o miss40_mac_bi_snps_2019_merged_melas.vcf.gz

tabix -p vcf miss40_mac_bi_snps_2019_merged_melas.vcf.gz

## count variants
preVARIANTS=$(bcftools view -H miss40_mac_bi_snps_2019_merged_melas.vcf.gz | wc -l)
echo "The number of variants in the vcf before gatk hard filtering is $preVARIANTS"

########### FILTER 3: 
# Use gatk VariantFiltration to tag SNPs which are poor quality and need removing
# Replace the reference -R, your vcf -V, and the output gatk_tagged VCF name.

gatk VariantFiltration \
-R VectorBase-62_AmelasCM1001059_A_Genome.fasta \
-V miss40_mac_bi_snps_2019_merged_melas.vcf.gz \
-filter "QD < 5.0" --filter-name "QD5" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
-O gatk_tagged_miss40_mac_bi_snps_2019_merged_melas.vcf.gz

## note there is a repeated warning here JEXL engine undefined variable, but it's not a problem
## because this is coming from where there are positions with no coverage
## https://gatk.broadinstitute.org/hc/en-us/community/posts/4408733963803-GATK-Variant-Filtration-undefined-variable

## remove the SNPs that have been tagged with the above filters
bcftools view -f 'PASS' gatk_tagged_miss40_mac_bi_snps_2019_merged_melas.vcf.gz -Oz -o gatk_filtered_miss40_mac_bi_snps_2019_merged_melas.vcf.gz

## count variants
VARIANTS=$(bcftools view -H gatk_filtered_miss40_mac_bi_snps_2019_merged_melas.vcf.gz | wc -l)
echo "The number of variants in the vcf after gatk hard filtering is $VARIANTS" 

## ADDITIONAL SOFT FILTERING

## After hard quality filtering, we have a data set containing only variant sites that we trust with high confidence. 
## However, so far we have only removed variants that had skewed values across all samples – 
## we haven’t yet looked at the individual genotypes at each variant site. 
## For example, we may have kept a variant that has a high quality across all of the individual samples, but there 
## may be some individuals where the read depth (FMT/DP) or the genotype quality (GQ) for the individual genotype is very low, 
## so we don’t trust the calls for these individuals

########### FILTER 4: Filter reads that have a read depth below 5 OR a genotype quality below 20
## The difference is, that instead of the whole variant site, we are now considering single genotypes 
## (the information in the “FORMAT” fields”) using -e 'FMT/DP<3 | FMT/GQ<20' and we will not remove the whole variant site, 
## but only set the respective genotypes to missing (./.) by using the bcftools filter -S . command.
## Here we use the logical operator “|” instead of “||”, since using “||” would mean that every genotype at a variant site 
## is set to missing, even if only one genotype doesn’t fulfill the threshold
bcftools filter -S . -e 'FMT/DP<5 | FMT/GQ<20' -O z -o DP5_GQ20_gatk_filtered_miss40_mac_bi_snps_2019_merged_melas.vcf.gz gatk_filtered_miss40_mac_bi_snps_2019_merged_melas.vcf.gz

## check this, the genotype here should have been set to missing ./. as we are filtering for depth below 5
bcftools query -i 'FMT/DP<5' -f '[GT=%GT\tDP=%DP\tGQ=%GQ\t]\n' DP5_GQ20_gatk_filtered_miss40_mac_bi_snps_2019_merged_melas.vcf.gz | less -S

## happy? index vcf
tabix -p vcf DP5_GQ20_gatk_filtered_miss40_mac_bi_snps_2019_merged_melas.vcf.gz

## count variants. Should be the same number as above all are filters above just set to missing, so genotypes removed but site still in VCF.
countVARIANTS=$(bcftools view -H DP5_GQ20_gatk_filtered_miss40_mac_bi_snps_2019_merged_melas.vcf.gz | wc -l)
echo "The number of variants after filtering FMT/DP<5 | FMT/GQ<20 is $countVARIANTS"

########### FILTER 5: exclude -e all sites at which no alternative alleles are called for any of the samples
bcftools filter -e 'AC==0' DP5_GQ20_gatk_filtered_miss40_mac_bi_snps_2019_merged_melas.vcf.gz -O z -o AC0_DP5_GQ20_gatk_filtered_miss40_mac_bi_snps_2019_merged_melas.vcf.gz

tabix -p vcf AC0_DP5_GQ20_gatk_filtered_miss40_mac_bi_snps_2019_merged_melas.vcf.gz

## count variants

countVARIANTSagain=$(bcftools view -H AC0_DP5_GQ20_gatk_filtered_miss40_mac_bi_snps_2019_merged_melas.vcf.gz | wc -l)
echo "The number of variants after filtering to exclude any non-alternate sites is $countVARIANTSagain"

########### FILTER 6
## Remove variants with a high amount of missing genotypes and filter on minor allele frequency
## The data set overall will now have lots of missing data, because we have replaced calls with depth <5 or quality below 20 with ./.  
## Therefore we will now remove all variants that have more than 20% missing genotypes or MAF < 0.01
## Filtering for MAF 0.01 means that remaining sites in the VCF need to have a minimum minor allele frequency of 1%
## So if a dataset contained 150 individuals (thus maximally 300 alleles per site because mosquitoes are diploid). If we filter for MAF 0.01, 
## a variant is excluded if its minor allele is found in less than 3 genotypes (1% of 300). 
## this would remove alleles where the minor allele occurs 3 times or less
## Note that depending on your sample size, you might want to adjust your MAF filter. 
# I have had the comment from a reviewer before that it needs to remove singletons (where there is only 1 occurrence of a particular SNP, but this depends on your research question)


bcftools filter -e 'F_MISSING > 0.2 || MAF <= 0.01' -O z -o F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_miss40_mac_bi_snps_2019_merged_melas.vcf.gz AC0_DP5_GQ20_gatk_filtered_miss40_mac_bi_snps_2019_merged_melas.vcf.gz

# check this has worked. The minor allele count should be 1 or above in this dataset.
bcftools query F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_miss40_mac_bi_snps_2019_merged_melas.vcf.gz -f'%AC\n' | sort -g | head

# count variants
variantcount=$(bcftools view -H F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_miss40_mac_bi_snps_2019_merged_melas.vcf.gz | wc -l)
echo "the number of variants after filtering for F_MISSING and MAF is $variantcount"

tabix -p vcf F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_miss40_mac_bi_snps_2019_merged_melas.vcf.gz

# Filtering complete!


########## PHASE VCF FILE ###########

# Phase the filtered vcf using beagle, https://faculty.washington.edu/browning/beagle/beagle_5.2_13Oct21.pdf
# This is necessary for some selection statistical analyses

# Recommend to use beagle in a separate conda environment. Conda install mamba. Mamba install beagle.
conda activate beagle
beagle -Xmx500g gt=final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.vcf.gz out=2019melasglobal_finalfiltered_gambiaealigned_phased

tabix -p vcf 2019melasglobal_finalfiltered_gambiaealigned_phased.vcf.gz

bcftools query -l 2019melasglobal_finalfiltered_gambiaealigned_phased.vcf.gz
bcftools index --stats 2019melasglobal_finalfiltered_gambiaealigned_phased.vcf.gz | cut -f1 | uniq

## Check the number of SNPs in the phased and unphased VCF files for one chromosome, these should be the same
# unphased melas SNPs in chr 2L 1725502
bcftools query -f '%CHROM\t%POS\n' final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.vcf.gz | awk '$1=="2L"' | wc -l
# phased melas SNPs in chr 2L 1725502
bcftools query -f '%CHROM\t%POS\n' 2019melasglobal_finalfiltered_gambiaealigned_phased.vcf.gz | awk '$1=="2L"' | wc -l

# This final vcf contains these chromosome names, with X and Y been renamed.
2L
2R
3L
3R
anop_mito
anop_X
Y_unplaced

###### snpEFF annotation of filtered or the filtered & phased VCF #####

### Here I used the Anopheles_gambiae database from snpEff. Change depending on species. Eg. Anopheles_darlingi_2 is the database for An. darlingi.

snpEff Anopheles_gambiae 2019melasglobal_finalfiltered_gambiaealigned_phased.vcf.gz > 2019melasglobal_finalfiltered_gambiaealigned_phased.ann.vcf.gz
tabix -p vcf 2019melasglobal_finalfiltered_gambiaealigned_phased.ann.vcf.gz

