#!/usr/bin/env bash
set -euo pipefail

# SNP CSV file. Expected columns:
# chrom,pos,id,ref,alt
SNPS="${1:-/home/kg522/data/TF_motif_analysis/variants_pass.csv}" # your prioritized variants

# Genome reference
GENOME="${2:-/home/kg522/data/ensembl-vep/vep_plugin_data/fasta/hg38.fa}" # reference genome FASTA
MOTIFS="${3:-$HOME/.conda/envs/HOMER/share/homer-4.10-0/motifs/human_motifs/human_only.motifs}" # selected human motifs

REF_FASTA="ref.fa"
ALT_FASTA="alt.fa"
TMP_FASTA="$(mktemp)"
trap 'rm -f "$TMP_FASTA"' EXIT

WINDOW=20 # window around SNP
CENTER=$WINDOW

: > "$REF_FASTA"
: > "$ALT_FASTA"

tail -n +2 "$SNPS" | while IFS=',' read -r chr pos id ref alt extra
do
    chr="$(printf '%s' "$chr" | tr -d '[:space:]')"
    pos="$(printf '%s' "$pos" | tr -d '[:space:]')"
    id="$(printf '%s' "$id" | sed 's/^[[:space:]]*//; s/[[:space:]]*$//')"
    ref="$(printf '%s' "$ref" | tr '[:lower:]' '[:upper:]' | tr -d '[:space:]')"
    alt="$(printf '%s' "$alt" | tr '[:lower:]' '[:upper:]' | tr -d '[:space:]')"

    [[ -z "$chr" ]] && continue

    if [[ -n "${extra:-}" ]]; then
        echo "ERROR: too many CSV columns for $id" >&2
        exit 1
    fi

    if [[ ! "$pos" =~ ^[0-9]+$ || -z "$id" || -z "$ref" || -z "$alt" ]]; then
        echo "ERROR: malformed CSV row: chr='$chr' pos='$pos' id='$id' ref='$ref' alt='$alt'" >&2
        exit 1
    fi

    if [[ ${#ref} -ne 1 || ${#alt} -ne 1 ]]; then
        echo "ERROR: $id is not a single-nucleotide substitution: ref='$ref' alt='$alt'" >&2
        exit 1
    fi

    start=$((pos - WINDOW))
    end=$((pos + WINDOW))

    if (( start < 1 )); then
        echo "ERROR: $id window starts before coordinate 1: ${chr}:${start}-${end}" >&2
        exit 1
    fi

    samtools faidx "$GENOME" "${chr}:${start}-${end}" > "$TMP_FASTA"

    seq="$(tail -n +2 "$TMP_FASTA" | tr -d '[:space:]' | tr '[:lower:]' '[:upper:]')"
    ref_base="${seq:$CENTER:1}"

    if [[ ${#seq} -ne $((WINDOW * 2 + 1)) ]]; then
        echo "ERROR: $id fetched ${#seq} bp, expected $((WINDOW * 2 + 1)) bp for ${chr}:${start}-${end}" >&2
        exit 1
    fi

    if [[ "$ref_base" != "$ref" ]]; then
        echo "ERROR: $id CSV ref '$ref' does not match genome base '$ref_base' at ${chr}:${pos}" >&2
        exit 1
    fi

    refseq="$seq"
    altseq="${seq:0:$CENTER}${alt}${seq:$((CENTER + 1))}"

    printf '>%s\n%s\n' "$id" "$refseq" >> "$REF_FASTA"

    printf '>%s\n%s\n' "$id" "$altseq" >> "$ALT_FASTA"

done

# Scan for human motifs using HOMER (NOT ENRICHMENT, just find matches)
findMotifs.pl "$REF_FASTA" fasta ref_scan -find "$MOTIFS" > ref_hits.txt
findMotifs.pl "$ALT_FASTA" fasta alt_scan -find "$MOTIFS" > alt_hits.txt
