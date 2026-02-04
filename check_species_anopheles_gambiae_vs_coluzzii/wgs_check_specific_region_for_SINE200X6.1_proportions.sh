#!/bin/bash

# Check species (An. coluzzii vs An. gambiae) in WGS data using SINE200 X6.1 insertion on chromosome X 
# https://link.springer.com/article/10.1186/1475-2875-7-163

TARGET="ACATCCTAAAATAATGAATTAAATGCAGTTCCTTATTATTCTCAGCTAACCGCTATTGTTCAGCGTGAACAATGCATTAAATTCTGAGTGGAGTAAAATGTATC"
RC_TARGET=$(echo "$TARGET" | tr 'ACGTacgt' 'TGCAtgca' | rev)

CHR="X"
INS_START=22951690
INS_END=22951793
INS_INT="${CHR}:${INS_START}-${INS_END}"   # insertion interval only (104 bp)

# Motif search window = insertion +/- buffer
BUF=500
INS_BUF_INT="${CHR}:$((INS_START-BUF))-$((INS_END+BUF))"

# Flanks (500bp each, 50bp gap)
FLANK_WIN=500
FLANK_GAP=50
LEFT_INT="${CHR}:$((INS_START-FLANK_GAP-FLANK_WIN))-$((INS_START-FLANK_GAP-1))"
RIGHT_INT="${CHR}:$((INS_END+FLANK_GAP+1))-$((INS_END+FLANK_GAP+FLANK_WIN))"

OUT="m2_coluzzii_insertion.tsv"

JOBS=6
THREADS=2
MM=2

mean_depth() {
  bam="$1"
  interval="$2"
  # mean depth across all positions in interval (includes zero-coverage positions via -a)
  samtools depth -a -r "$interval" "$bam" \
    | awk '{sum+=$3; n++} END{ if(n==0) print "NA"; else printf "%.6f", sum/n }'
}

process_one() {
  bam="$1"
  sample=$(basename "$bam" .bam)

  # Ensure BAM index exists
  if [[ ! -f "${bam}.bai" && ! -f "${bam%.bam}.bai" ]]; then
    samtools index -@ "$THREADS" "$bam" 2>/dev/null || samtools index "$bam"
  fi

  # Baseline flank coverage (all reads)
  left_mean=$(mean_depth "$bam" "$LEFT_INT")
  right_mean=$(mean_depth "$bam" "$RIGHT_INT")
  flank_mean=$(awk -v l="$left_mean" -v r="$right_mean" 'BEGIN{ if(l=="NA"||r=="NA") print "NA"; else printf "%.6f",(l+r)/2 }')

  # ---- Identify motif-containing readnames within INS+buffer ----
  tmp_names=$(mktemp -p . "${sample}.motif_names.XXXXXX")
  tmp_bam=$(mktemp -p . "${sample}.motif_only.XXXXXX.bam")

  # IMPORTANT: samtools fastq here needs BAM input on stdin -> use samtools view -bh
  samtools view -@ "$THREADS" -F 0x900 -bh "$bam" "$INS_BUF_INT" \
    | samtools fastq -n -@ "$THREADS" - \
    | seqkit grep -s -i -m "$MM" -p "$TARGET" -p "$RC_TARGET" -n \
    | awk '{name=$1; sub(/^@/,"",name); sub(/\/[12]$/,"",name); print name}' \
    | LC_ALL=C sort -u > "$tmp_names"

  motif_readnames=$(wc -l < "$tmp_names" | tr -d '[:space:]')

  # ---- Build motif-only BAM (alignments in INS+buffer whose QNAME is in tmp_names) ----
  if [[ "$motif_readnames" -eq 0 ]]; then
    motif_alignments_in_region=0
    ins_mean="NA"
  else
    samtools view -h -@ "$THREADS" -F 0x900 "$bam" "$INS_BUF_INT" \
      | awk -v names="$tmp_names" '
          BEGIN{ while((getline < names) > 0) keep[$1]=1 }
          /^@/ {print; next}
          ($1 in keep) {print}
        ' \
      | samtools view -b - > "$tmp_bam"

    motif_alignments_in_region=$(samtools view -c "$tmp_bam")

    # Index motif-only BAM and compute mean depth over INSERTION INTERVAL ONLY
    samtools index -@ "$THREADS" "$tmp_bam" 2>/dev/null || samtools index "$tmp_bam"
    ins_mean=$(mean_depth "$tmp_bam" "$INS_INT")
  fi

  rm -f "$tmp_names" "$tmp_bam" "${tmp_bam}.bai" 2>/dev/null || true

  # Percent of flank baseline
  ins_depth_as_pct_of_flanks=$(awk -v i="$ins_mean" -v f="$flank_mean" '
    BEGIN{ if(i=="NA"||f=="NA"||f==0) print "NA"; else printf "%.2f",(100*i)/f }')

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$sample" "$motif_alignments_in_region" "$motif_readnames" "$left_mean" "$right_mean" "$flank_mean" "$ins_mean" "$ins_depth_as_pct_of_flanks"
}

export -f process_one mean_depth
export TARGET RC_TARGET CHR INS_START INS_END INS_INT BUF INS_BUF_INT FLANK_WIN FLANK_GAP LEFT_INT RIGHT_INT THREADS MM

{
  printf "sample\tmotif_alignments_in_region\tmotif_readnames\tleft_mean\tright_mean\tflank_mean\tins_mean\tins_depth_as_pct_of_flanks\n"
  if command -v parallel >/dev/null 2>&1; then
    parallel -j "$JOBS" process_one ::: *.bam
  else
    ls -1 *.bam | xargs -n1 -P "$JOBS" -I{} bash -c 'process_one "$@"' _ {}
  fi
} > "$OUT"