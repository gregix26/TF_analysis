import os
import pandas as pd
import matplotlib.pyplot as plt

# ---------------------
# CONFIG
# ---------------------
BASE_DIR = "/home/kg522/data/TF_motif_analysis/GPNMB_output_Menon_per_class"
INPUT_CSV = os.path.join(BASE_DIR, "top_motif_hits_with_tf_families.csv")
TARGET_SNP = "rs199357"


def load_top_hits(input_csv):
    df = pd.read_csv(input_csv)
    if df.empty:
        raise RuntimeError("No motif hits found in input CSV.")

    df = df.copy()
    df["Abs_delta_score"] = pd.to_numeric(df["Abs_delta_score"], errors="coerce")
    df["Delta_score"] = pd.to_numeric(df["Delta_score"], errors="coerce")
    df["TF_family"] = df["TF_family"].fillna("Unassigned").astype(str)
    df["Direction"] = df["Direction"].fillna("NA").astype(str)

    top_hits = (
        df.sort_values(
            ["SNP", "Cell_type", "Abs_delta_score", "Delta_score"],
            ascending=[True, True, False, False],
        )
        .groupby(["SNP", "Cell_type"], as_index=False)
        .first()
    )

    top_hits = top_hits[["SNP", "Cell_type", "TF_family", "Direction", "Delta_score", "Abs_delta_score"]]
    return top_hits


def render_table_panel(ax, panel_df, title):
    ax.axis("off")
    col_labels = ["Cell type", "Top TF family", "Gain/Loss", "Delta"]
    cell_text = [
        [
            row["Cell_type"],
            row["TF_family"],
            row["Direction"],
            f"{row['Delta_score']:.3f}",
        ]
        for _, row in panel_df.iterrows()
    ]

    table = ax.table(
        cellText=cell_text,
        colLabels=col_labels,
        cellLoc="left",
        rowLoc="center",
        loc="center",
        colWidths=[0.23, 0.37, 0.16, 0.12],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.4)

    for (row, col), cell in table.get_celld().items():
        if row == 0:
            cell.set_facecolor("#dbe9f6")
            cell.set_text_props(weight="bold")
        else:
            direction = panel_df.iloc[row - 1]["Direction"]
            if direction == "Gain":
                cell.set_facecolor("#dbe9f6")
            elif direction == "Loss":
                cell.set_facecolor("#f4cccc")

    ax.set_title(title, fontsize=12, weight="bold", pad=12)


def main():
    if not os.path.exists(INPUT_CSV):
        raise FileNotFoundError(f"Input CSV not found: {INPUT_CSV}")

    top_hits = load_top_hits(INPUT_CSV)
    snp_df = top_hits[top_hits["SNP"] == TARGET_SNP].copy()
    if snp_df.empty:
        raise RuntimeError(f"No rows found for TARGET_SNP={TARGET_SNP}")

    snp_df = snp_df.sort_values(
        ["Abs_delta_score", "Delta_score", "Cell_type"],
        ascending=[False, False, True],
    ).reset_index(drop=True)

    output_csv = os.path.join(BASE_DIR, f"{TARGET_SNP}_top_tf_family_table.csv")
    output_png = os.path.join(BASE_DIR, f"{TARGET_SNP}_top_tf_family_table.png")
    snp_df[["Cell_type", "TF_family", "Direction", "Delta_score"]].to_csv(output_csv, index=False)

    fig, ax = plt.subplots(
        nrows=1,
        ncols=1,
        figsize=(12, max(3.5, 0.5 * len(snp_df) + 1.5)),
        constrained_layout=True,
    )
    render_table_panel(ax, snp_df, f"{TARGET_SNP}: top TF family by cell type")
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"[DONE] Saved table CSV to: {output_csv}")
    print(f"[DONE] Saved table plot to: {output_png}")


if __name__ == "__main__":
    main()
