
## Identifying Structural Variants Using DELLY

## Using delly call ALL

# read the sample identifiers from all_samples.txt, process the .all.bcf files, in parallel using -P 10,
# then convert this to a vcf and then convert it to a text file without the header lines for easy viewing.
# Produces individual text files for each sample

ls *.mkdup.bam | sed 's/.bam//' > delly_all_samples.txt

for f in *.mkdup.bam; do delly call -t ALL -g Anopheles_gambiae.AgamP4.dna.toplevel.fa "$f" -o "${f%.*}.all.bcf"; done

delly merge -o structural_variants.bcf *.all.bcf

bcftools view structural_variants.bcf > structural_variants.vcf

bgzip structural_variants.vcf

## Have a look at the structural_variants.vcf.gz
zgrep -v ^"##" structural_variants.vcf.gz | less

## genotype the structural variants 
cat delly_all_samples.txt | parallel -j 15 --bar "delly call -g Anopheles_gambiae.AgamP4.dna.toplevel.fa -v structural_variants.bcf -o {}.all.genotyped.bcf  {}.bam"

# %% merge all of the genotyped bcf files together
bcftools merge -m id -O b -o merged.bcf *.all.genotyped.bcf

# list samples in merged.bcf
bcftools query -l merged.bcf

# convert to vcf to see all structural variants in all samples (unfiltered)
bcftools convert -Oz -o merged.vcf.gz merged.bcf

# check sample names
bcftools view -h merged.vcf.gz | grep '^#CHROM' | cut -f10- > merged_sample_names.txt

## Filter - keep samples that were used in downstream analysis - edit delly_good_samples.txt for this (removed samples with over 20% missing data)

bcftools view -S delly_good_samples.txt merged.bcf -Oz -o merged_genotyped_structural_variants_sample_filt.vcf.gz
tabix -p vcf merged_genotyped_structural_variants_sample_filt.vcf.gz

## Filter again to remove SVs with >20% missing data
bcftools view merged_genotyped_structural_variants_sample_filt.vcf.gz | bcftools view -i 'F_PASS(GT!="mis")>0.8' -Oz -o merged_genotyped_structural_variants_sample_filt_miss20.vcf.gz
tabix -p vcf merged_genotyped_structural_variants_sample_filt_miss20.vcf.gz

# Query number of variants
bcftools view merged_genotyped_structural_variants_sample_filt_miss20.vcf.gz | grep -v "^#" | wc -l
# 74043 structural variants in total across entire genome for filtered samples and missingness

# Filter the VCF based on genes of interest
bcftools view -R genes_of_interest_+-1kb.bed merged_genotyped_structural_variants_sample_filt_miss20.vcf.gz -Oz -o genes_merged_genotyped_structural_variants_sample_filt_miss20.vcf.gz
tabix -p vcf genes_merged_genotyped_structural_variants_sample_filt_miss20.vcf.gz
bcftools view genes_merged_genotyped_structural_variants_sample_filt_miss20.vcf.gz | grep -v "^#" | wc -l
# 129 variants in genes of interest

# Extract sample names
bcftools view -h genes_merged_genotyped_structural_variants_sample_filt_miss20.vcf.gz | grep '^#CHROM' | cut -f10- > filtered_SV_vcf_sample_names.txt
# Extract information from vcf
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/SVTYPE\t[%GT\t]\n' genes_merged_genotyped_structural_variants_sample_filt_miss20.vcf.gz > structural_variants.txt

# Annotate with snpeff

snpEff Anopheles_gambiae genes_merged_genotyped_structural_variants_sample_filt_miss20.vcf.gz > snpeff_genes_merged_genotyped_structural_variants_sample_filt_miss20.vcf
bgzip snpeff_genes_merged_genotyped_structural_variants_sample_filt_miss20.vcf
tabix -p vcf snpeff_genes_merged_genotyped_structural_variants_sample_filt_miss20.vcf.gz

# Extract information from vcf with snpeff
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/SVTYPE\t%INFO/ANN\t[%GT\t]\n' snpeff_genes_merged_genotyped_structural_variants_sample_filt_miss20.vcf.gz > snpeff_structural_variants.txt
























# Query what the unique Structural Variant Types are in the dataset
# zgrep -E "SVTYPE" genes_merged_genotyped_structural_variants_sample_filt.site_filt.call_filt.FMISSING_MAF.vcf.gz | cut -f8 | grep -o "SVTYPE=[^;]*" | sort | uniq


#####  Separate by variant type

#zgrep -E '^#|SVTYPE=DUP' genes_merged_genotyped_structural_variants_sample_filt.site_filt.call_filt.vcf.gz | bgzip > duplications.vcf.gz
#bcftools view duplications.vcf.gz | grep -v "^#" | wc -l
#snpEff Anopheles_gambiae duplications.vcf.gz > duplications.ann.vcf
#bgzip duplications.ann.vcf
#tabix -p vcf duplications.ann.vcf.gz
#

### Analysis Type 1: between Bijagos and Elsewhere Fst ###

# Will now look to see if there are any in the Bijagos melas (pop1) vs non-Bijagos melas (pop2)
# Made two populations.txt files pop1.txt (Bijagos samples) and pop2.txt (cameroon/gambia samples)
vcftools --gzvcf genes_merged_genotyped_structural_variants_sample_filt.site_filt.call_filt.vcf.gz --weir-fst-pop pop1.txt --weir-fst-pop pop2.txt --out Bijagos_vs_Elsewhere



## creat mappability map

dicey chop Anopheles_gambiae.AgamP4.dna.toplevel.fa
bwa index Anopheles_gambiae.AgamP4.dna.toplevel.fa
bwa mem Anopheles_gambiae.AgamP4.dna.toplevel.fa car1002_Combined_1.fastq.gz car1002_Combined_2.fastq.gz | samtools sort -@ 8 -o srt.bam -
samtools index srt.bam 
dicey mappability2 srt.bam 
gunzip map.fa.gz && bgzip map.fa && samtools faidx map.fa.gz 



### possible snpeff addition

#
#bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/SVTYPE\t[%GT\t]\n' genes_merged_genotyped_structural_variants_sample_filt.site_filt.call_filt.ann.vcf.gz > structural_variants_snpeff.txt
#

## Looking for chromosomal inversions using compkaryo, which looks at SNPs

# use the 'test' conda environment

# works for '2La', '2Rj', '2Rb', '2Rc_col', '2Rc_gam', '2Rd', '2Ru'
python /mnt/storage11/sophie/gitrepos/compkaryo/compkaryo/compkaryo.py 2L_only_2019melasglobal_finalfiltered_gambiaealigned_phased.vcf.gz 2La -s compkaryo_samples.txt -o compkaryo_out_2La.txt --total

python /mnt/storage11/sophie/gitrepos/compkaryo/compkaryo/compkaryo.py 2R_only_2019melasglobal_finalfiltered_gambiaealigned_phased.vcf.gz 2Rb -s compkaryo_samples.txt -o compkaryo_out_2Rb.txt --total

python /mnt/storage11/sophie/gitrepos/compkaryo/compkaryo/compkaryo.py 2R_only_2019melasglobal_finalfiltered_gambiaealigned_phased.vcf.gz 2Rj -s compkaryo_samples.txt -o compkaryo_out_2Rj.txt --total
python /mnt/storage11/sophie/gitrepos/compkaryo/compkaryo/compkaryo.py 2R_only_2019melasglobal_finalfiltered_gambiaealigned_phased.vcf.gz 2Rc_gam -s compkaryo_samples.txt -o compkaryo_out_2Rc_gam.txt --total
python /mnt/storage11/sophie/gitrepos/compkaryo/compkaryo/compkaryo.py 2R_only_2019melasglobal_finalfiltered_gambiaealigned_phased.vcf.gz 2Rd -s compkaryo_samples.txt -o compkaryo_out_2Rd.txt --total
python /mnt/storage11/sophie/gitrepos/compkaryo/compkaryo/compkaryo.py 2R_only_2019melasglobal_finalfiltered_gambiaealigned_phased.vcf.gz 2Ru -s compkaryo_samples.txt -o compkaryo_out_2Ru.txt --total
