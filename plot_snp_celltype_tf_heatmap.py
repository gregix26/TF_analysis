import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ---------------------
# CONFIG
# ---------------------
BASE_DIR = "/home/kg522/data/TF_motif_analysis/GPNMB_output_Menon_per_class"
INPUT_CSV = os.path.join(BASE_DIR, "top_motif_hits_with_tf_families.csv")
OUTPUT_PNG = os.path.join(BASE_DIR, "snp_celltype_tf_family_heatmap.png")
OUTPUT_CSV = os.path.join(BASE_DIR, "snp_celltype_top_tf_family_heatmap_values.csv")


def load_top_hits(input_csv):
    df = pd.read_csv(input_csv)
    if df.empty:
        raise RuntimeError("No motif hits found in input CSV.")

    df = df.copy()
    df["Abs_delta_score"] = pd.to_numeric(df["Abs_delta_score"], errors="coerce")
    df["Delta_score"] = pd.to_numeric(df["Delta_score"], errors="coerce")
    df["TF_family"] = df["TF_family"].fillna("Unassigned").astype(str)

    top_hits = (
        df.sort_values(
            ["SNP", "Cell_type", "Abs_delta_score", "Delta_score"],
            ascending=[True, True, False, False],
        )
        .groupby(["SNP", "Cell_type"], as_index=False)
        .first()
    )
    return top_hits


def main():
    if not os.path.exists(INPUT_CSV):
        raise FileNotFoundError(f"Input CSV not found: {INPUT_CSV}")

    top_hits = load_top_hits(INPUT_CSV)
    top_hits.to_csv(OUTPUT_CSV, index=False)

    cell_types = sorted(top_hits["Cell_type"].drop_duplicates().tolist())
    snps = sorted(top_hits["SNP"].drop_duplicates().tolist())

    value_matrix = (
        top_hits.pivot(index="Cell_type", columns="SNP", values="Delta_score")
        .reindex(index=cell_types, columns=snps)
    )
    family_matrix = (
        top_hits.pivot(index="Cell_type", columns="SNP", values="TF_family")
        .reindex(index=cell_types, columns=snps)
    )

    values = value_matrix.to_numpy(dtype=float)
    vmax = np.nanmax(np.abs(values))
    if not np.isfinite(vmax) or vmax == 0:
        vmax = 1.0

    fig_width = max(10, 0.45 * len(snps) + 4)
    fig_height = max(5, 0.45 * len(cell_types) + 2.5)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height), constrained_layout=True)

    im = ax.imshow(values, cmap="coolwarm", vmin=-vmax, vmax=vmax, aspect="auto")

    ax.set_xticks(np.arange(len(snps)))
    ax.set_xticklabels(snps, rotation=90, fontsize=8)
    ax.set_yticks(np.arange(len(cell_types)))
    ax.set_yticklabels(cell_types, fontsize=9)
    ax.set_xlabel("SNP id")
    ax.set_ylabel("Cell types")
    ax.set_title("Top TF Family per SNP and Cell Type", fontsize=15, weight="bold")

    for i, cell_type in enumerate(cell_types):
        for j, snp_id in enumerate(snps):
            value = value_matrix.loc[cell_type, snp_id]
            family = family_matrix.loc[cell_type, snp_id]
            if pd.isna(value):
                continue

            label = str(family)
            if len(label) > 18:
                label = label[:15] + "..."
            text_color = "white" if abs(value) > 0.5 * vmax else "black"
            ax.text(j, i, label, ha="center", va="center", fontsize=7, color=text_color)

    cbar = fig.colorbar(im, ax=ax, shrink=0.9)
    cbar.set_label("Delta score", rotation=90)

    fig.savefig(OUTPUT_PNG, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"[DONE] Saved heatmap values to: {OUTPUT_CSV}")
    print(f"[DONE] Saved heatmap to: {OUTPUT_PNG}")


if __name__ == "__main__":
    main()
