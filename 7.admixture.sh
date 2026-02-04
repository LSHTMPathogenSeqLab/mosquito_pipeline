## Admixture

## used this vcf

cp AC_MAF_PBO_exposed_F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_filt_coluzzii.2025_02_28.genotyped.ann.vcf.gz admixture.vcf.gz
# Identify which samples are in your vcf:
bcftools query -l yourfile.vcf.gz

# STEP 1 EXTRA FILTERING for MAF if you need to do it

# STEP 2 CONVERT CHR NAMES IN VCF TO INTEGERS

zgrep -v "^#" admixture.vcf.gz | cut -f1 | sort | uniq

zcat AC_MAF_PBO_exposed_F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_filt_coluzzii.2025_02_28.genotyped.ann.vcf.gz \
| awk 'BEGIN{OFS=FS="\t"}
  /^##contig/ {
    gsub(/NC_064604.1/, "1");
    gsub(/NC_064669.1/, "2");
    gsub(/NC_064670.1/, "3");
    gsub(/NC_064671.1/, "4");
    print;
    next
  }
  /^#/ {print; next}
  {
    gsub(/^NC_064604.1$/, "1", $1);
    gsub(/^NC_064669.1$/, "2", $1);
    gsub(/^NC_064670.1$/, "3", $1);
    gsub(/^NC_064671.1$/, "4", $1);
    print
  }' | bgzip > admixture_renamed_chrs.vcf.gz

tabix -p vcf admixture_renamed_chrs.vcf.gz

zgrep -v "^#" admixture_renamed_chrs.vcf.gz | cut -f1 | sort | uniq

bcftools query -l admixture_renamed_chrs.vcf.gz

# STEP 3 MAKE BED AND BIM FILES

plink --vcf admixture_renamed_chrs.vcf.gz \
      --set-missing-var-ids @:# \
      --keep-allele-order \
      --const-fid \
      --allow-extra-chr \
      --make-bed \
      --out guinea_coluzzii_pbo_exposed

# Basic QC for Admixture

plink --bfile guinea_coluzzii_pbo_exposed \
      --geno 0.05 \
      --maf 0.01 \
      --make-bed \
      --out guinea_coluzzii_pbo_exposed_qc

# STEP 4 RUN ADMIXTURE
# Make a file called K_runs.txt with values 1 to 10 and then run the following command
# Run this in screen as it takes a while (fastq2matrix conda environment)

cat K_runs.txt | xargs -I {} sh -c 'admixture --cv=10 -j20 -s 74926 ghana_gambiae_2025_admixture_pruned.bed {} | tee log{}.cv10.seed74926.out'

# Inspect the CV values for the inflection point to check what the best K value is for your dataset
# Inspect using grep -h CV *out

# Plot admixture using R script
# Change the directory and metadata below as needed

conda create -n radmix r-essentials r-base
install.packages(c("unikn", "countrycode", "optparse"))
setwd("/mnt/storage11/sophie/admixture")

Rscript /mnt/storage11/sophie/gitrepos/mosquito_pipeline/generate_admix_barplot_colours.R \
-d /mnt/storage11/sophie/Anisa/wgs_project/admixture \
--prefix ghana_gambiae_2025_admixture_pruned \
--kval 4 \
-m metadata.tsv \
--filter_N 1 \
--label_id sampleID \
--label_region region \
--label_country country \
--label_site site \
--country_code TRUE


### ____ Troubleshooting ___ ###
bcftools query -f '[%CHROM\t%POS\t%GT\n]' F_MISS_MAF_admixture_modified.vcf.gz | \
awk '{print $1, $2, $3}' | sort | uniq -c | grep -w '1' > invariant_sites.txt

bcftools query -f '[%CHROM\t%POS\t%GT\n]' PBO_exposed_F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_filt_coluzzii.2025_02_28.genotyped.ann.vcf.gz | \
awk '{if ($3 == "./.") print $0}' | sort | uniq -c | \ awk '{if ($1 == <number_of_samples>) print $2, $3}' > all_missing_sites.txt

