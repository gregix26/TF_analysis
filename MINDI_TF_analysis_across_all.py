# Compare motif similarity between ALT and REF across all available cell types
# and summarize the strongest gain/loss effects across the combined cell-type set.

import os
import re
import numpy as np
import pandas as pd
from tqdm import tqdm
import tfmindi as tm

# ---------------------
# CONFIG
# ---------------------
BASE_DIR = "/home/kg522/data/TF_motif_analysis/"
WINDOW = 20  # number of bases on each side of the SNP to consider for motif similarity
TOP_K = 3  # top motifs across all cell types
MIN_DELTA = 0.05  # ignore very small motif-similarity changes


def empirical_percentile(values, target_value):
    values = np.asarray(values, dtype=float)
    return 100.0 * np.mean(values <= float(target_value))


def interpret_delta_strength(delta_score):
    abs_delta = abs(float(delta_score))

    if abs_delta > 0.3:
        return "very strong"
    if abs_delta > 0.2:
        return "strong"
    if abs_delta >= 0.1:
        return "moderate"
    if abs_delta >= 0.05:
        return "weak"
    return "noise"


def discover_cell_types(base_dir):
    cell_types = set()
    pattern = re.compile(r"^contrib_(?:ref|alt)_(.+)\.npy$")

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


def compute_cell_type_metrics(ref, alt, motif_collection):
    delta_full = alt - ref
    center = ref.shape[0] // 2
    s, e = center - WINDOW, center + WINDOW

    ref_crop = ref[s:e].T
    alt_crop = alt[s:e].T

    sim_ref = tm.pp.calculate_motif_similarity([ref_crop], motif_collection)
    sim_alt = tm.pp.calculate_motif_similarity([alt_crop], motif_collection)

    sim_ref = sim_ref.toarray().flatten()
    sim_alt = sim_alt.toarray().flatten()

    motif_names = list(motif_collection.keys())
    sim_ref = pd.Series(sim_ref, index=motif_names)
    sim_alt = pd.Series(sim_alt, index=motif_names)
    delta_similarity = sim_alt - sim_ref

    importance = np.abs(delta_full).sum(axis=1)
    snp_score = float(importance[center])
    snp_percentile = empirical_percentile(importance, snp_score)
    z_score = (snp_score - importance.mean()) / (importance.std() + 1e-6)

    delta_percentiles = pd.Series(
        [
            empirical_percentile(np.abs(delta_similarity.values), abs(delta_score))
            for delta_score in delta_similarity.values
        ],
        index=delta_similarity.index,
    )

    return {
        "snp_score": snp_score,
        "snp_percentile": float(snp_percentile),
        "z_score": float(z_score),
        "sim_ref": sim_ref,
        "sim_alt": sim_alt,
        "delta_similarity": delta_similarity,
        "delta_percentiles": delta_percentiles,
    }


def motif_effect(delta_similarity):
    delta_similarity = float(delta_similarity)
    if delta_similarity >= MIN_DELTA:
        return "Gain"
    if delta_similarity <= -MIN_DELTA:
        return "Loss"
    return "No change"


def build_motif_rows(snp_id, cell_type, metrics):
    sim_ref_map = metrics["sim_ref"].to_dict()
    sim_alt_map = metrics["sim_alt"].to_dict()
    delta_percentiles_map = metrics["delta_percentiles"].to_dict()

    rows = []
    for motif_name, delta_score in metrics["delta_similarity"].items():
        ref_similarity = float(sim_ref_map[motif_name])
        alt_similarity = float(sim_alt_map[motif_name])
        delta_similarity = float(delta_score)
        direction = motif_effect(delta_similarity)
        rows.append({
            "SNP": snp_id,
            "Cell_type": cell_type,
            "Direction": direction,
            "Motif": motif_name,
            "SNP_contribution": round(metrics["snp_score"], 3),
            "SNP_importance_percentile": round(metrics["snp_percentile"], 2),
            "Ref_similarity": round(ref_similarity, 6),
            "Alt_similarity": round(alt_similarity, 6),
            "Delta_similarity": round(delta_similarity, 6),
            "Delta_strength": interpret_delta_strength(delta_similarity),
            "Delta_abs_percentile": round(float(delta_percentiles_map[motif_name]), 2),
        })

    return rows


def summarize_top_hits(motif_rows_df, top_k):
    gain_df = motif_rows_df[motif_rows_df["Delta_similarity"] >= MIN_DELTA]
    loss_df = motif_rows_df[motif_rows_df["Delta_similarity"] <= -MIN_DELTA]

    top_gain_df = gain_df.nlargest(top_k, "Delta_similarity")
    top_loss_df = loss_df.nsmallest(top_k, "Delta_similarity")
    top_hits_df = pd.concat([top_gain_df, top_loss_df], ignore_index=True)

    gain_str = ";".join(
        f"{row.Motif}@{row.Cell_type}"
        for row in top_gain_df.itertuples(index=False)
    )
    loss_str = ";".join(
        f"{row.Motif}@{row.Cell_type}"
        for row in top_loss_df.itertuples(index=False)
    )

    tf_effect = float(top_hits_df["Delta_similarity"].abs().mean()) if not top_hits_df.empty else 0.0
    return top_gain_df, top_loss_df, gain_str, loss_str, tf_effect


# ---------------------
# 1. Load motif database
# ---------------------
print("[INFO] Loading motif database...")
motif_collection = tm.load_motif_collection(tm.fetch_motif_collection())

all_cell_types = discover_cell_types(BASE_DIR)
if not all_cell_types:
    raise RuntimeError(f"No contrib_ref/contrib_alt files found in {BASE_DIR}")

print(f"[INFO] Discovered cell types: {', '.join(all_cell_types)}")

results = []
motif_detail_results = []

print("[INFO] Processing SNPs across all cell types...")

for snp_id in tqdm(os.listdir(BASE_DIR)):
    snp_path = os.path.join(BASE_DIR, snp_id)
    if not os.path.isdir(snp_path):
        continue

    try:
        snp_cell_types = available_cell_types_for_snp(snp_path, all_cell_types)
        if not snp_cell_types:
            print(f"[WARN] Skipping {snp_id}: no matching ref/alt files for any cell type")
            continue

        per_cell_type_metrics = []
        snp_motif_rows = []

        for cell_type in snp_cell_types:
            ref = np.load(os.path.join(snp_path, f"contrib_ref_{cell_type}.npy"))
            alt = np.load(os.path.join(snp_path, f"contrib_alt_{cell_type}.npy"))

            metrics = compute_cell_type_metrics(ref, alt, motif_collection)
            per_cell_type_metrics.append({
                "Cell_type": cell_type,
                "SNP_importance": round(metrics["snp_score"], 3),
                "SNP_importance_percentile": round(metrics["snp_percentile"], 2),
                "Z_score": round(metrics["z_score"], 3),
            })
            snp_motif_rows.extend(build_motif_rows(snp_id, cell_type, metrics))

        if not snp_motif_rows:
            print(f"[WARN] Skipping {snp_id}: no motif rows produced")
            continue

        snp_motif_df = pd.DataFrame(snp_motif_rows)
        top_gain_df, top_loss_df, gain_str, loss_str, tf_effect = summarize_top_hits(
            snp_motif_df,
            TOP_K,
        )
        top_hits_df = pd.concat([top_gain_df, top_loss_df], ignore_index=True)

        per_cell_type_df = pd.DataFrame(per_cell_type_metrics)
        mean_importance = float(per_cell_type_df["SNP_importance"].mean())
        mean_importance_percentile = float(per_cell_type_df["SNP_importance_percentile"].mean())
        mean_z_score = float(per_cell_type_df["Z_score"].mean())
        max_importance = float(per_cell_type_df["SNP_importance"].max())
        max_z_score = float(per_cell_type_df["Z_score"].max())
        functional_score = mean_z_score + tf_effect

        results.append({
            "SNP": snp_id,
            "Cell_types_tested": len(snp_cell_types),
            "Cell_types": ";".join(snp_cell_types),
            "Lost_TFs": loss_str,
            "Gained_TFs": gain_str,
            "Mean_SNP_importance": round(mean_importance, 3),
            "Mean_SNP_importance_percentile": round(mean_importance_percentile, 2),
            "Max_SNP_importance": round(max_importance, 3),
            "Mean_Z_score": round(mean_z_score, 3),
            "Max_Z_score": round(max_z_score, 3),
            "TF_effect": round(tf_effect, 3),
            "Functional_score": round(functional_score, 3),
        })

        motif_detail_results.extend(top_hits_df.to_dict(orient="records"))

    except Exception as e:
        print(f"[WARN] Skipping {snp_id}: {e}")

# ---------------------
# 2. Build tables
# ---------------------
df = pd.DataFrame(results)
motif_detail_df = pd.DataFrame(motif_detail_results)

if not df.empty:
    df = df.sort_values("Functional_score", ascending=False)
    df["Rank"] = range(1, len(df) + 1)

if not motif_detail_df.empty:
    motif_detail_df = motif_detail_df.sort_values(
        ["SNP", "Direction", "Delta_similarity"],
        ascending=[True, True, False],
    )

# ---------------------
# 3. Save
# ---------------------
summary_out_csv = os.path.join(BASE_DIR, "ranked_functional_snps_all_cell_types.csv")
detail_out_csv = os.path.join(BASE_DIR, "ranked_functional_snps_all_cell_types_motif_details.csv")

df.to_csv(summary_out_csv, index=False)
motif_detail_df.to_csv(detail_out_csv, index=False)

# ---------------------
# 4. Display
# ---------------------
print("\n=== Ranked SNPs by Functional Impact Across All Cell Types ===\n")
print(df.to_string(index=False))

print("\n=== Top Gain/Loss Motif Details Across All Cell Types ===\n")
print(motif_detail_df.to_string(index=False))

print(f"\n[DONE] Saved summary to: {summary_out_csv}")
print(f"[DONE] Saved motif details to: {detail_out_csv}")
