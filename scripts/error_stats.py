import pandas as pd
import numpy as np


def load_tsv(path: str) -> pd.DataFrame:
    df = pd.read_csv(
        path,
        sep=r"\s+",
        header=None,
        names=["CHROM", "POS", "REF", "ALT", "FREQ"],
        dtype=str,
    )

    df["POS"] = df["POS"].astype(int)

    freq_clean = (
        df["FREQ"]
        .astype(str)
        .str.replace("%", "", regex=False)
        .str.replace(",", ".", regex=False)
    )
    df["FREQ"] = freq_clean.astype(float) / 100.0

    return df


ctrl1 = load_tsv("results/controls/control1_rare_simple.tsv")
ctrl2 = load_tsv("results/controls/control2_rare_simple.tsv")
ctrl3 = load_tsv("results/controls/control3_rare_simple.tsv")
roommate = load_tsv("results/roommate/roommate_rare_simple.tsv")


stats = []

for name, df in [("control1", ctrl1), ("control2", ctrl2), ("control3", ctrl3)]:
    mean = df["FREQ"].mean()
    sd = df["FREQ"].std()
    stats.append((name, mean, sd))

stats_df = pd.DataFrame(stats, columns=["sample", "mean_error", "sd_error"])

print("Per-control error statistics (fraction units):")
print(stats_df.to_string(index=False))

thresholds = stats_df["mean_error"].values + 3 * stats_df["sd_error"].values
print("\nThresholds (mean + 3*sd):")
for name, thr in zip(stats_df["sample"], thresholds):
    print(f"{name}: {thr:.6f}")


t1, t2, t3 = thresholds
roommate["HIGH_CONF"] = (
    (roommate["FREQ"] > t1)
    & (roommate["FREQ"] > t2)
    & (roommate["FREQ"] > t3)
)

hc = roommate[roommate["HIGH_CONF"]]

print("\nHigh-confidence variants in roommate sample (> mean+3sd of all controls):")
if hc.empty:
    print("No variants above all three thresholds.")
else:
    print(hc[["CHROM", "POS", "REF", "ALT", "FREQ"]].to_string(index=False))
