"""
SNP → Contribution Scores Pipeline from CRESted
===============================================
For each SNP:
  1. Extract ±250bp reference and alternative sequences from FASTA
  2. One-hot encode both
  3. Compute per-nucleotide contribution scores (once per sequence)
  4. Save per-SNP directory with one-hot arrays and contribution scores

Output per SNP:
    {output}/{snp_id}/
        ref_onehot.npy       (500, 4)
        alt_onehot.npy       (500, 4)
        contrib_ref.npy      (500, 4)
        contrib_alt.npy      (500, 4)
        contrib_delta.npy    (500, 4)   alt - ref contribution
        meta.json            coordinates + alleles

Usage:
    python snp_contrib_pipeline.py \\
        --snp_csv  snps.csv \\
        --fasta    hg38.fa \\
        --output   results/
"""

import os
import json
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
import pysam
import crested

# ── Column names — edit if your CSV differs ──────────────────
COL_SNP_ID = "id"
COL_CHROM  = "chrom"
COL_POS    = "pos"     # 1-based
COL_REF    = "ref"
COL_ALT    = "alt"
# ─────────────────────────────────────────────────────────────


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--snp_csv",        required=True)
    p.add_argument("--fasta",          required=True)
    p.add_argument("--output",         required=True)
    p.add_argument("--model",          default="DeepHumanCortex1")
    p.add_argument("--window",         type=int, default=250,
                   help="bp each side of SNP → 500bp total (default: 250)")
    p.add_argument("--contrib_method", default="expected_integrated_grad",
                   choices=["expected_integrated_grad", "smooth_grad", "grad_times_input"])
    return p.parse_args()


def load_snps(csv_path):
    df = pd.read_csv(csv_path)
    df[COL_CHROM] = df[COL_CHROM].astype(str).apply(
        lambda c: c if c.startswith("chr") else f"chr{c}"
    )
    df[COL_REF] = df[COL_REF].str.upper()
    df[COL_ALT] = df[COL_ALT].str.upper()
    df[COL_POS] = df[COL_POS].astype(int)
    print(f"[INFO] Loaded {len(df)} SNPs")
    return df


def build_windows(df, window):
    df = df.copy()
    center           = df[COL_POS] - 1          # 1-based → 0-based
    df["win_start"]  = (center - window).clip(lower=0)
    df["win_end"]    = center + window
    df["snp_offset"] = window
    return df


def extract_sequences(df, fasta_path, seq_len):
    """
    Extract sequences and return:
      - ref_seqs, alt_seqs : lists of strings (only valid SNPs)
      - valid_idx          : list of integer positions in df that passed validation

    SNPs are skipped (not crashed) if:
      - their window is shorter than seq_len (edge-of-chromosome)
      - their chromosome is not found in the FASTA
    """
    fasta        = pysam.FastaFile(fasta_path)
    available    = set(fasta.references)
    ref_seqs     = []
    alt_seqs     = []
    valid_idx    = []
    skipped      = []

    for i, (_, row) in enumerate(tqdm(df.iterrows(), total=len(df), desc="Extracting sequences")):
        snp_id = row[COL_SNP_ID]
        chrom  = row[COL_CHROM]

        # Guard 1: chromosome not in FASTA
        if chrom not in available:
            print(f"[SKIP] {snp_id}: chromosome '{chrom}' not in FASTA")
            skipped.append((snp_id, "chromosome not found"))
            continue

        seq = fasta.fetch(chrom, int(row["win_start"]), int(row["win_end"])).upper()

        # Guard 2: sequence too short (edge-of-chromosome clipping)
        if len(seq) != seq_len:
            print(f"[SKIP] {snp_id}: got {len(seq)}bp, expected {seq_len}bp "
                  f"— likely near chromosome edge")
            skipped.append((snp_id, f"wrong length {len(seq)}bp"))
            continue

        offset = int(row["snp_offset"])

        # Guard 3: offset out of bounds (shouldn't happen, but be safe)
        if offset >= len(seq):
            print(f"[SKIP] {snp_id}: snp_offset {offset} out of range for seq len {len(seq)}")
            skipped.append((snp_id, "offset out of range"))
            continue

        # Guard 4: ref allele mismatch (warn but keep)
        if seq[offset] != row[COL_REF]:
            print(f"[WARN] {snp_id}: FASTA '{seq[offset]}' ≠ declared ref '{row[COL_REF]}' "
                  f"— keeping SNP, using FASTA base as reference")

        ref_seqs.append(seq)
        alt_seqs.append(seq[:offset] + row[COL_ALT] + seq[offset + 1:])
        valid_idx.append(i)

    fasta.close()

    if skipped:
        print(f"\n[INFO] Skipped {len(skipped)} SNPs:")
        for snp_id, reason in skipped:
            print(f"         {snp_id}: {reason}")

    # Guard 5: after filtering, ensure all sequences are the same length
    # (protects np.stack in compute_contributions from the ValueError)
    lengths = set(len(s) for s in ref_seqs)
    if len(lengths) > 1:
        raise ValueError(
            f"Sequences have inconsistent lengths after extraction: {lengths}. "
            f"This should not happen — check your BED regions."
        )

    print(f"[INFO] {len(valid_idx)}/{len(df)} SNPs passed extraction")
    return ref_seqs, alt_seqs, valid_idx


def one_hot_batch(seqs):
    """List of DNA strings → float32 array (N, L, 4)  [A, C, G, T]"""
    base_map = {"A": 0, "C": 1, "G": 2, "T": 3}
    arr = np.zeros((len(seqs), len(seqs[0]), 4), dtype=np.float32)
    for i, seq in enumerate(seqs):
        for j, base in enumerate(seq):
            idx = base_map.get(base)
            if idx is not None:
                arr[i, j, idx] = 1.0
    return arr


def compute_contributions(batch, model, output_names, method):
    """
    Compute contribution scores for a batch of sequences, separately per cell type.

    Returns a dict keyed by cell type name:
        { "AST": np.ndarray (N, 500, 4),
          "EXC_L2_3_IT": np.ndarray (N, 500, 4),
          ... }
    """
    scores_by_class = {}

    for ct_idx, ct_name in enumerate(tqdm(output_names, desc="  Contrib classes", leave=False)):
        result = crested.tl.contribution_scores(
            input=batch,
            model=model,
            all_class_names=output_names,
            target_idx=ct_idx,
            method=method,
        )
        # target_idx=int  →  result has length 1, scores at result[0]
        # squeeze removes the spurious size-1 class dimension: (N,1,500,4) → (N,500,4)
        scores = np.squeeze(result[0])

        if scores.shape != batch.shape:
            raise ValueError(
                f"Unexpected shape for class {ct_name!r}: {scores.shape} "
                f"(after squeeze). Expected {batch.shape}."
            )

        scores_by_class[ct_name] = scores   # (N, 500, 4)

    return scores_by_class   # dict: cell_type_name -> (N, 500, 4)
def main():
    args    = parse_args()
    seq_len = 2 * args.window
    os.makedirs(args.output, exist_ok=True)

    # ── Load & prepare ───────────────────────────────────────
    df = load_snps(args.snp_csv)
    df = build_windows(df, args.window)

    # extract_sequences returns valid_idx so df stays in sync with batches
    ref_seqs, alt_seqs, valid_idx = extract_sequences(df, args.fasta, seq_len)
    df = df.iloc[valid_idx].reset_index(drop=True)   # drop skipped SNPs

    if len(df) == 0:
        print("[ERROR] No valid SNPs remaining after extraction. Check your FASTA and coordinates.")
        return

    ref_batch = one_hot_batch(ref_seqs)   # (N, 500, 4)
    alt_batch = one_hot_batch(alt_seqs)   # (N, 500, 4)

    # ── Load model ───────────────────────────────────────────
    print(f"[INFO] Loading model: {args.model}")
    model_path, output_names = crested.get_model(args.model)
    model = crested.utils.load_model(model_path)

    # ── Contribution scores — one call per cell type, ref and alt ───
    print("[INFO] Computing contribution scores (ref) ...")
    contrib_ref = compute_contributions(ref_batch, model, output_names, args.contrib_method)
    # contrib_ref: dict { cell_type_name -> (N, 500, 4) }

    print("[INFO] Computing contribution scores (alt) ...")
    contrib_alt = compute_contributions(alt_batch, model, output_names, args.contrib_method)

    # ── Save per SNP ─────────────────────────────────────────
    print("[INFO] Saving outputs ...")
    for i, (_, row) in enumerate(tqdm(df.iterrows(), total=len(df), desc="Saving")):
        snp_id = str(row[COL_SNP_ID])
        outdir = os.path.join(args.output, snp_id)
        os.makedirs(outdir, exist_ok=True)

        # One-hot sequences
        np.save(f"{outdir}/ref_onehot.npy", ref_batch[i])   # (500, 4)
        np.save(f"{outdir}/alt_onehot.npy", alt_batch[i])   # (500, 4)

        # Contribution scores — one file per cell type
        for ct_name in output_names:
            ct_safe = ct_name.replace(" ", "_").replace("/", "-")
            ref_ct  = contrib_ref[ct_name][i]                      # (500, 4)
            alt_ct  = contrib_alt[ct_name][i]                      # (500, 4)
            np.save(f"{outdir}/contrib_ref_{ct_safe}.npy",   ref_ct)
            np.save(f"{outdir}/contrib_alt_{ct_safe}.npy",   alt_ct)
            np.save(f"{outdir}/contrib_delta_{ct_safe}.npy", alt_ct - ref_ct)

        with open(f"{outdir}/meta.json", "w") as f:
            json.dump({
                "snp_id":     snp_id,
                "chrom":      row[COL_CHROM],
                "pos_1based": int(row[COL_POS]),
                "win_start":  int(row["win_start"]),
                "win_end":    int(row["win_end"]),
                "snp_offset": int(row["snp_offset"]),
                "ref_allele": row[COL_REF],
                "alt_allele": row[COL_ALT],
                "cell_types": output_names,
            }, f, indent=2)

    print(f"\n[DONE] {len(df)} SNP directories written to: {args.output}")
    print(f"  Files per SNP: ref_onehot.npy, alt_onehot.npy")
    print(f"  Per cell type: contrib_ref_<CT>.npy, contrib_alt_<CT>.npy, contrib_delta_<CT>.npy")


if __name__ == "__main__":
    main()
