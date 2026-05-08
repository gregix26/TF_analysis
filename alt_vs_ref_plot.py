# For plotting the CRESted contribution scores per SNP
# Produces two PNGs: ALT and REF with matching intervals and a visual vline where the SNP sits
# If SNP sits somewhere else, change vline position!

import numpy as np
import matplotlib.pyplot as plt
import crested

snp_dir = "/home/kg522/data/TF_motif_analysis/"

# --- load data ---
ref_scores = np.load(f"{snp_dir}/contrib_ref_VASC.npy")
alt_scores = np.load(f"{snp_dir}/contrib_alt_VASC.npy")

ref_seq = np.load(f"{snp_dir}/ref_onehot.npy")
alt_seq = np.load(f"{snp_dir}/alt_onehot.npy")

# --- add batch dim ---
ref_scores = ref_scores[np.newaxis, ...]
alt_scores = alt_scores[np.newaxis, ...]

ref_seq = ref_seq[np.newaxis, ...]
alt_seq = alt_seq[np.newaxis, ...]

# SNP position
snp_pos = 50
zoom_n_bases = 100

# Use the same y-axis scale and tick intervals for REF and ALT.
combined_scores = np.concatenate([ref_scores.ravel(), alt_scores.ravel()])
y_min = combined_scores.min()
y_max = combined_scores.max()
y_pad = (y_max - y_min) * 0.05

if y_pad == 0:
    y_pad = 0.01

y_limits = (y_min - y_pad, y_max + y_pad)
x_limits = (0, zoom_n_bases)
x_ticks = np.arange(0, zoom_n_bases + 1, 10)


def apply_matching_intervals(ax):
    ax.set_xlim(x_limits)
    ax.set_ylim(y_limits)
    ax.set_xticks(x_ticks)

# --- REF ---
crested.pl.explain.contribution_scores(
    ref_scores,
    ref_seq,
    sequence_labels=["REF"],
    class_labels=["VASC"],
    zoom_n_bases=zoom_n_bases,
    suptitle=None,
    highlight_positions=None
)
ax = plt.gca()
apply_matching_intervals(ax)

# draw vertical line
ax.axvline(
    x=snp_pos,
    color="red",
    linestyle="--",
    linewidth=2,
    alpha=0.8
)
ax.set_title("Reference")
plt.savefig(f"{snp_dir}/ref_VASC.png", dpi=300, bbox_inches="tight")
plt.close()

# --- ALT ---
crested.pl.explain.contribution_scores(
    alt_scores,
    alt_seq,
    sequence_labels=["ALT"],
    class_labels=["VASC"],
    zoom_n_bases=zoom_n_bases,
    suptitle=None,
    highlight_positions=None
)
ax = plt.gca()
apply_matching_intervals(ax)

# draw vertical line
ax.axvline(
    x=snp_pos,
    color="red",
    linestyle="--",
    linewidth=2,
    alpha=0.8
)

ax.set_title("Alternative")
plt.savefig(f"{snp_dir}/alt_VASC.png", dpi=300, bbox_inches="tight")
plt.close()
