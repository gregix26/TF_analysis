import os
import re
import ast
import numpy as np
import pandas as pd
from tqdm import tqdm
import tfmindi as tm

# ---------------------
# CONFIG
# ---------------------
BASE_DIR = "/home/kg522/data/TF_motif_analysis/GPNMB_output_Menon_per_class"
WINDOW = 20
TOP_K = 3


def discover_cell_types(base_dir):
    pattern = re.compile(r"^contrib_ref_(.+)\.npy$")
    cell_types = set()

    for snp_id in os.listdir(base_dir):
        snp_path = os.path.join(base_dir, snp_id)
        if not os.path.isdir(snp_path):
            continue

        for filename in os.listdir(snp_path):
            match = pattern.match(filename)
            if match:
                cell_types.add(match.group(1))

    return sorted(cell_types)


def available_cell_types_for_snp(snp_path, cell_types):
    available = []

    for cell_type in cell_types:
        ref_path = os.path.join(snp_path, f"contrib_ref_{cell_type}.npy")
        alt_path = os.path.join(snp_path, f"contrib_alt_{cell_type}.npy")
        if os.path.exists(ref_path) and os.path.exists(alt_path):
            available.append(cell_type)

    return available


def normalize_motif_key(motif_name):
    if isinstance(motif_name, tuple):
        return "|".join(str(x) for x in motif_name)
    return str(motif_name)


def motif_lookup_candidates(motif_name):
    candidates = []

    def add_candidate(value):
        if value is None:
            return
        if isinstance(value, float) and pd.isna(value):
            return
        if value not in candidates:
            candidates.append(value)

    add_candidate(motif_name)
    add_candidate(normalize_motif_key(motif_name))

    if isinstance(motif_name, tuple):
        for item in reversed(motif_name):
            add_candidate(item)
            add_candidate(str(item))

    elif isinstance(motif_name, str):
        text = motif_name.strip()
        add_candidate(text)

        if text.startswith("(") and text.endswith(")"):
            try:
                parsed = ast.literal_eval(text)
                if isinstance(parsed, tuple):
                    for item in reversed(parsed):
                        add_candidate(item)
                        add_candidate(str(item))
            except (ValueError, SyntaxError):
                pass

    return candidates


def build_tf_family_lookup(motif_annotations, motif_to_dbd):
    lookup = {}

    if isinstance(motif_to_dbd, dict):
        for key, value in motif_to_dbd.items():
            for candidate in motif_lookup_candidates(key):
                lookup[candidate] = value

    if isinstance(motif_annotations, pd.DataFrame):
        candidate_name_cols = [
            "motif",
            "motif_name",
            "name",
            "motif_id",
            "id",
        ]
        candidate_family_cols = [
            "dbd",
            "dbd_family",
            "tf_family",
            "family",
            "tfclass",
        ]

        name_col = next((col for col in candidate_name_cols if col in motif_annotations.columns), None)
        family_col = next((col for col in candidate_family_cols if col in motif_annotations.columns), None)

        if name_col and family_col:
            for _, row in motif_annotations[[name_col, family_col]].dropna().iterrows():
                for candidate in motif_lookup_candidates(row[name_col]):
                    lookup[candidate] = row[family_col]

    return lookup


def resolve_tf_family(motif_name, family_lookup):
    for candidate in motif_lookup_candidates(motif_name):
        if candidate in family_lookup:
            return family_lookup[candidate]

    return "Unassigned"


def calculate_motif_scores(ref, alt, motif_collection):
    center = ref.shape[0] // 2
    start = center - WINDOW
    end = center + WINDOW

    ref_crop = ref[start:end].T
    alt_crop = alt[start:end].T

    sim_ref = tm.pp.calculate_motif_similarity([ref_crop], motif_collection)
    sim_alt = tm.pp.calculate_motif_similarity([alt_crop], motif_collection)

    motif_names = list(motif_collection.keys())
    sim_ref = pd.Series(sim_ref.toarray().flatten(), index=motif_names)
    sim_alt = pd.Series(sim_alt.toarray().flatten(), index=motif_names)
    delta_tf = sim_alt - sim_ref

    return sim_ref, sim_alt, delta_tf


def build_top_hit_rows(snp_id, cell_type, sim_ref, sim_alt, delta_tf, family_lookup):
    top_gain = delta_tf.sort_values(ascending=False).head(TOP_K)
    top_loss = delta_tf.sort_values().head(TOP_K)
    selected = pd.concat([top_gain, top_loss])

    sim_ref_map = sim_ref.to_dict()
    sim_alt_map = sim_alt.to_dict()
    rows = []

    for motif_name, delta_score in selected.items():
        rows.append({
            "SNP": snp_id,
            "Cell_type": cell_type,
            "Direction": "Gain" if float(delta_score) >= 0 else "Loss",
            "Motif": motif_name,
            "TF_family": resolve_tf_family(motif_name, family_lookup),
            "Ref_similarity": round(float(sim_ref_map[motif_name]), 6),
            "Alt_similarity": round(float(sim_alt_map[motif_name]), 6),
            "Delta_score": round(float(delta_score), 6),
            "Abs_delta_score": round(abs(float(delta_score)), 6),
        })

    return rows


print("[INFO] Loading motif resources...")
motif_collection = tm.load_motif_collection(tm.fetch_motif_collection())
motif_annotations = tm.load_motif_annotations(tm.fetch_motif_annotations())
motif_to_dbd = tm.load_motif_to_dbd(motif_annotations)
family_lookup = build_tf_family_lookup(motif_annotations, motif_to_dbd)

all_cell_types = discover_cell_types(BASE_DIR)
if not all_cell_types:
    raise RuntimeError(f"No contrib_ref/contrib_alt files found in {BASE_DIR}")

print(f"[INFO] Discovered cell types: {', '.join(all_cell_types)}")

motif_detail_rows = []
summary_rows = []

print("[INFO] Processing SNPs...")
for snp_id in tqdm(sorted(os.listdir(BASE_DIR))):
    snp_path = os.path.join(BASE_DIR, snp_id)
    if not os.path.isdir(snp_path):
        continue

    snp_cell_types = available_cell_types_for_snp(snp_path, all_cell_types)
    if not snp_cell_types:
        continue

    for cell_type in snp_cell_types:
        try:
            ref = np.load(os.path.join(snp_path, f"contrib_ref_{cell_type}.npy"))
            alt = np.load(os.path.join(snp_path, f"contrib_alt_{cell_type}.npy"))
            sim_ref, sim_alt, delta_tf = calculate_motif_scores(ref, alt, motif_collection)
            top_rows = build_top_hit_rows(
                snp_id,
                cell_type,
                sim_ref,
                sim_alt,
                delta_tf,
                family_lookup,
            )
            motif_detail_rows.extend(top_rows)

            top_df = pd.DataFrame(top_rows)
            gain_families = top_df.loc[top_df["Direction"] == "Gain", "TF_family"].astype(str)
            loss_families = top_df.loc[top_df["Direction"] == "Loss", "TF_family"].astype(str)

            summary_rows.append({
                "SNP": snp_id,
                "Cell_type": cell_type,
                "Top_gain_motifs": ";".join(top_df.loc[top_df["Direction"] == "Gain", "Motif"].astype(str)),
                "Top_loss_motifs": ";".join(top_df.loc[top_df["Direction"] == "Loss", "Motif"].astype(str)),
                "Top_gain_families": ";".join(gain_families),
                "Top_loss_families": ";".join(loss_families),
            })
        except Exception as e:
            print(f"[WARN] Skipping {snp_id} / {cell_type}: {e}")

motif_detail_df = pd.DataFrame(motif_detail_rows)
summary_df = pd.DataFrame(summary_rows)

if not motif_detail_df.empty:
    family_summary_df = (
        motif_detail_df
        .groupby(["Cell_type", "Direction", "TF_family"], dropna=False)
        .agg(
            Hit_count=("TF_family", "size"),
            Unique_SNPs=("SNP", "nunique"),
            Mean_abs_delta=("Abs_delta_score", "mean"),
            Max_abs_delta=("Abs_delta_score", "max"),
        )
        .reset_index()
        .sort_values(
            ["Cell_type", "Direction", "Unique_SNPs", "Mean_abs_delta"],
            ascending=[True, True, False, False],
        )
    )
else:
    family_summary_df = pd.DataFrame(
        columns=[
            "Cell_type",
            "Direction",
            "TF_family",
            "Hit_count",
            "Unique_SNPs",
            "Mean_abs_delta",
            "Max_abs_delta",
        ]
    )

motif_out = os.path.join(BASE_DIR, "top_motif_hits_with_tf_families.csv")
summary_out = os.path.join(BASE_DIR, "top_tf_families_per_snp_cell_type.csv")
family_out = os.path.join(BASE_DIR, "tf_family_summary_by_cell_type.csv")

motif_detail_df.to_csv(motif_out, index=False)
summary_df.to_csv(summary_out, index=False)
family_summary_df.to_csv(family_out, index=False)

print("\n=== Top motif hits with TF families ===\n")
print(motif_detail_df.to_string(index=False))

print("\n=== SNP/cell type family summary ===\n")
print(summary_df.to_string(index=False))

print("\n=== Aggregate TF family summary by cell type ===\n")
print(family_summary_df.to_string(index=False))

print(f"\n[DONE] Saved motif-level output to: {motif_out}")
print(f"[DONE] Saved SNP summary to: {summary_out}")
print(f"[DONE] Saved family summary to: {family_out}")
