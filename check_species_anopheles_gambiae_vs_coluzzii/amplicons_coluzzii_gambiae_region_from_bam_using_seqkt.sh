#!/bin/bash
# Check species (gambiae s.s. vs coluzzii) using SINE200 insertion in SINE amplicons
# We are looking for the SINE200 X6.1 insertion which is in coluzzii but not in gambiae
# https://link.springer.com/article/10.1186/1475-2875-7-163
# We look inside the SINE amplicons that we have generated (their bam files)
set -euo pipefail

TARGET="ACATCCTAAAATAATGAATTAAATGCAGTTCCTTATTATTCTCAGCTAACCGCTATTGTTCAGCGTGAACAATGCATTAAATTCTGAGTGGAGTAAAATGTATC"
OUT="m2_an_colgamb_region_reads_seqkitfuzzy.txt"
> "$OUT"

for bam in *.bam; do
  sample=$(basename "$bam" .bam)

  # samtools fasta: filter out secondary (0x100) and supplementary (0x800) with -F 0x900
  # send samtools/seqkit messages to stderr (or silence them with 2>/dev/null if you prefer)
  count=$(
    samtools fasta -F 0x900 "$bam" 2>/dev/null \
    | seqkit grep -s -i -m 2 -p "$TARGET" -n 2>/dev/null \
    | sort -u \
    | wc -l
  )

  echo -e "${sample}\t${count}" >> "$OUT"
done
