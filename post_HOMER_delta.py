import pandas as pd

# HOMER output format:
# SequenceID   Offset   Sequence   MotifName   Strand   Score

ref = pd.read_csv("ref_hits.txt", sep="\t", header=None)
alt = pd.read_csv("alt_hits.txt", sep="\t", header=None)

ref.columns = ["id","offset","seq","motif","strand","score"]
alt.columns = ["id","offset","seq","motif","strand","score"]

ref["score"] = pd.to_numeric(ref["score"], errors="coerce")
alt["score"] = pd.to_numeric(alt["score"], errors="coerce")

ref["offset"] = pd.to_numeric(ref["offset"], errors="coerce")
alt["offset"] = pd.to_numeric(alt["offset"], errors="coerce")

# Motif lenght is needed to determine if SNP is within motif
ref["motif_len"] = ref["seq"].str.len()
alt["motif_len"] = alt["seq"].str.len()

SNP_POS = 0

ref_overlap = ref[
    (ref["offset"] <= SNP_POS) &
    (ref["offset"] + ref["motif_len"] > SNP_POS)
]

alt_overlap = alt[
    (alt["offset"] <= SNP_POS) &
    (alt["offset"] + alt["motif_len"] > SNP_POS)
]

ref_best = ref_overlap.groupby(["id","motif"])["score"].max().reset_index()
alt_best = alt_overlap.groupby(["id","motif"])["score"].max().reset_index()

df = pd.merge(
    ref_best,
    alt_best,
    on=["id","motif"],
    how="outer",
    suffixes=("_ref","_alt")
)

df["score_ref"] = df["score_ref"].fillna(0)
df["score_alt"] = df["score_alt"].fillna(0)

# Calculate score difference
df["delta"] = df["score_alt"] - df["score_ref"]

#Classify effect based on score changes
def classify(row):
    if row["score_ref"] == 0 and row["score_alt"] > 0:
        return "gained"
    elif row["score_ref"] > 0 and row["score_alt"] == 0:
        return "lost"
    elif row["delta"] > 1:
        return "strengthened"
    elif row["delta"] < -1:
        return "weakened"
    else:
        return "neutral"

df["effect"] = df.apply(classify, axis=1)

df.to_csv("motif_delta.tsv", sep="\t", index=False)