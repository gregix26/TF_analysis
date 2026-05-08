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
BASE_DIR = "/home/kg522/data/TF_motif_analysis/GPNMB_output_Menon_per_class"
WINDOW = 20  # number of bases on each side of the SNP to consider for motif similarity
TOP_K = 3  # top motifs across all cell types


def empirical_percentile(values, target_value):
    values = np.asarray(values, dtype=float)
    return 100.0 * np.mean(values <= float(target_value))


def empirical_two_sided_pvalues(values):
    if isinstance(values, pd.Series):
        index = values.index
        numeric_values = values.to_numpy(dtype=float)
    else:
        index = None
        numeric_values = np.asarray(values, dtype=float)

    abs_values = np.abs(numeric_values)
    n = len(abs_values)
    return pd.Series(
        [((abs_values >= value).sum() + 1) / (n + 1) for value in abs_values],
        index=index,
    )


def benjamini_hochberg(pvalues):
    pvalues = np.asarray(pvalues, dtype=float)
    n = len(pvalues)

    if n == 0:
        return np.array([])

    order = np.argsort(pvalues)
    ranked = pvalues[order]
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)

    qvalues = np.empty(n, dtype=float)
    qvalues[order] = adjusted
    return qvalues


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
    delta_tf = sim_alt - sim_ref
    delta_abs_tf = sim_alt.abs() - sim_ref.abs()

    importance = np.abs(delta_full).sum(axis=1)
    snp_score = float(importance[center])
    snp_percentile = empirical_percentile(importance, snp_score)
    z_score = (snp_score - importance.mean()) / (importance.std() + 1e-6)

    motif_pvalues = empirical_two_sided_pvalues(delta_abs_tf)
    motif_qvalues = pd.Series(
        benjamini_hochberg(motif_pvalues.values),
        index=delta_tf.index,
    )
    motif_abs_percentiles = pd.Series(
        [
            empirical_percentile(np.abs(delta_abs_tf.values), abs(delta_score))
            for delta_score in delta_abs_tf.values
        ],
        index=delta_tf.index,
    )

    return {
        "snp_score": snp_score,
        "snp_percentile": float(snp_percentile),
        "z_score": float(z_score),
        "sim_ref": sim_ref,
        "sim_alt": sim_alt,
        "delta_tf": delta_tf,
        "delta_abs_tf": delta_abs_tf,
        "motif_pvalues": motif_pvalues,
        "motif_qvalues": motif_qvalues,
        "motif_abs_percentiles": motif_abs_percentiles,
    }


def contribution_sign(similarity_score):
    similarity_score = float(similarity_score)
    if similarity_score > 0:
        return "Positive"
    if similarity_score < 0:
        return "Negative"
    return "Neutral"


def interpret_motif_effect(ref_similarity, alt_similarity):
    ref_similarity = float(ref_similarity)
    alt_similarity = float(alt_similarity)
    strength_delta = abs(alt_similarity) - abs(ref_similarity)
    alt_sign = contribution_sign(alt_similarity)
    ref_sign = contribution_sign(ref_similarity)

    if strength_delta > 0:
        strength_direction = "Gain"
    elif strength_delta < 0:
        strength_direction = "Loss"
    else:
        strength_direction = "No change"

    if ref_sign != alt_sign and ref_sign != "Neutral" and alt_sign != "Neutral":
        return f"Sign flip to {alt_sign.lower()}"
    if strength_direction == "No change":
        return f"No strength change ({alt_sign.lower()})"
    return f"{alt_sign} motif {strength_direction.lower()}"


def motif_strength_direction(motif_strength_delta):
    motif_strength_delta = float(motif_strength_delta)
    if motif_strength_delta > 0:
        return "Gain"
    if motif_strength_delta < 0:
        return "Loss"
    return "No change"


def build_motif_rows(snp_id, cell_type, metrics):
    sim_ref_map = metrics["sim_ref"].to_dict()
    sim_alt_map = metrics["sim_alt"].to_dict()
    delta_abs_tf_map = metrics["delta_abs_tf"].to_dict()
    motif_pvalues_map = metrics["motif_pvalues"].to_dict()
    motif_qvalues_map = metrics["motif_qvalues"].to_dict()
    motif_abs_percentiles_map = metrics["motif_abs_percentiles"].to_dict()

    rows = []
    for motif_name, delta_score in metrics["delta_tf"].items():
        ref_similarity = float(sim_ref_map[motif_name])
        alt_similarity = float(sim_alt_map[motif_name])
        motif_strength_delta = float(delta_abs_tf_map[motif_name])
        direction = motif_strength_direction(motif_strength_delta)
        rows.append({
            "SNP": snp_id,
            "Cell_type": cell_type,
            "Direction": direction,
            "Motif": motif_name,
            "SNP_contribution": round(metrics["snp_score"], 3),
            "SNP_importance_percentile": round(metrics["snp_percentile"], 2),
            "Ref_similarity": round(ref_similarity, 6),
            "Alt_similarity": round(alt_similarity, 6),
            "Signed_delta_score": round(float(delta_score), 6),
            "Motif_strength_delta": round(motif_strength_delta, 6),
            "Alt_contribution_sign": contribution_sign(alt_similarity),
            "Interpretable_effect": interpret_motif_effect(ref_similarity, alt_similarity),
            "Delta_strength": interpret_delta_strength(motif_strength_delta),
            "Delta_abs_percentile": round(float(motif_abs_percentiles_map[motif_name]), 2),
            "Empirical_pvalue": round(float(motif_pvalues_map[motif_name]), 6),
            "BH_FDR": round(float(motif_qvalues_map[motif_name]), 6),
        })

    return rows


def summarize_top_hits(motif_rows_df, top_k):
    gain_df = motif_rows_df[motif_rows_df["Motif_strength_delta"] > 0]
    loss_df = motif_rows_df[motif_rows_df["Motif_strength_delta"] < 0]

    top_gain_df = gain_df.nlargest(top_k, "Motif_strength_delta")
    top_loss_df = loss_df.nsmallest(top_k, "Motif_strength_delta")
    top_hits_df = pd.concat([top_gain_df, top_loss_df], ignore_index=True)

    gain_str = ";".join(
        f"{row.Motif}@{row.Cell_type}({row.Alt_contribution_sign})"
        for row in top_gain_df.itertuples(index=False)
    )
    loss_str = ";".join(
        f"{row.Motif}@{row.Cell_type}({row.Alt_contribution_sign})"
        for row in top_loss_df.itertuples(index=False)
    )

    tf_effect = float(top_hits_df["Motif_strength_delta"].abs().mean()) if not top_hits_df.empty else 0.0
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
        ["SNP", "Direction", "Motif_strength_delta"],
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
