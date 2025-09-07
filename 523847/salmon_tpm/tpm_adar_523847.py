# boxplot_tpm_FIXED.py
# ---------------------------------------------
# Analysis pipeline for ADAR gene expression (TPM):
# - log1p transform
# - Type-III ANOVA: logTPM ~ treatment * phase + C(zt)
# - Planned contrasts: Exercise vs Sedentary within each Phase (BH-FDR)
# - Hedges g effect sizes
# - New figure names (won't overwrite your originals)
# ---------------------------------------------

import os
import math
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm
from statsmodels.stats.multitest import multipletests

# ---------------------------
# Config
# ---------------------------
INPUT_CSV = "salmonTPMwideGenesDataTable.csv"  # adjust if needed
OUTPUT_DIR = "results_tpm_FIXED"
FIG_DPI = 220
RANDOM_SEED = 42

NEW_FIG_TPL = "NEW_TPM_boxplot_by_ZT_{gene}.png"
ANOVA_CSV_TPL = "ANOVA_TypeIII_{gene}.csv"
POSTHOC_CSV = "PostHoc_Exercise_vs_Sedentary_byPhase_ALLGENES.csv"
CELLCOUNT_CSV = "CellCounts_treatment_phase_zt.csv"

# Optional: your preferred order for legend colors (treatment + phase)
GROUP_PALETTE = {
    "Exercise | Rest phase (ZT3)": "#f4b6b6",     # soft red
    "Sedentary | Rest phase (ZT3)": "#bfc6ff",   # soft blue
    "Exercise | Active phase (ZT15)": "#b51d1d", # deep red
    "Sedentary | Active phase (ZT15)": "#1022c0" # deep blue
}

np.random.seed(RANDOM_SEED)
warnings.filterwarnings("ignore", category=FutureWarning)

# ---------------------------
# Helpers
# ---------------------------

def standardize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize column names and gene names."""
    df = df.copy()
    df.columns = df.columns.str.strip()

    rename = {}
    for c in ["Treatment", "treatment"]:
        if c in df.columns: rename[c] = "treatment"
    for c in ["Phase", "phase"]:
        if c in df.columns: rename[c] = "phase"
    for c in ["ZT", "zt", "time", "Time"]:
        if c in df.columns: rename[c] = "zt"
    for c in ["Sample", "sample", "Mouse", "mouse", "Mouse_ID", "mouse_id", "animal", "Animal", "ID", "Id"]:
        if c in df.columns: rename[c] = "Sample"

    df = df.rename(columns=rename)

    # Try to standardize ADAR gene columns
    gene_map = {}
    for col in df.columns:
        low = col.lower()
        if low == "adar" and col != "Adar":
            gene_map[col] = "Adar"
        elif low == "adarb1" and col != "Adarb1":
            gene_map[col] = "Adarb1"
        elif low == "adarb2" and col != "Adarb2":
            gene_map[col] = "Adarb2"
    df = df.rename(columns=gene_map)

    # Coerce categorical fields
    if "treatment" in df.columns:
        df["treatment"] = df["treatment"].astype("category")
    if "phase" in df.columns:
        df["phase"] = df["phase"].astype("category")
    if "zt" in df.columns:
        df["zt"] = pd.Categorical(df["zt"], ordered=False)

    return df

def check_required_columns(df, required):
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

def melt_to_long(df: pd.DataFrame, gene_cols) -> pd.DataFrame:
    long_df = df.melt(
        id_vars=[c for c in ["Sample", "treatment", "phase", "zt"] if c in df.columns],
        value_vars=gene_cols,
        var_name="gene",
        value_name="TPM"
    ).dropna(subset=["TPM"]).reset_index(drop=True)
    return long_df

def hedges_g(x, y):
    """Hedges g (unbiased Cohen's d) for two independent samples."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    nx, ny = len(x), len(y)
    if nx < 2 or ny < 2:
        return np.nan
    vx = np.var(x, ddof=1)
    vy = np.var(y, ddof=1)
    # pooled SD
    sp = math.sqrt(((nx - 1)*vx + (ny - 1)*vy) / (nx + ny - 2))
    if sp == 0:
        return 0.0
    d = (np.mean(x) - np.mean(y)) / sp
    # small-sample correction
    J = 1.0 - (3.0 / (4.0*(nx + ny) - 9.0))
    return J * d

def ensure_outdir(path):
    os.makedirs(path, exist_ok=True)
    return path

# ---------------------------
# Main
# ---------------------------

def main():
    ensure_outdir(OUTPUT_DIR)

    # Load & standardize
    df = pd.read_csv(INPUT_CSV, encoding="ISO-8859-1")
    df = standardize_columns(df)
    check_required_columns(df, ["treatment", "phase", "zt"])

    # Which genes exist
    gene_cols = [c for c in ["Adar", "Adarb1", "Adarb2"] if c in df.columns]
    if not gene_cols:
        raise ValueError("No ADAR gene columns found (Adar, Adarb1, Adarb2).")

    # Long format for plotting
    long_df = melt_to_long(df, gene_cols)

    # Add group label for legend
    long_df["Group"] = (long_df["treatment"].astype(str) + " | " + long_df["phase"].astype(str))

    # Quick cell-count table (treatment x phase x zt)
    counts = (long_df.groupby(["treatment", "phase", "zt"])
              .size().reset_index(name="n"))
    counts.to_csv(os.path.join(OUTPUT_DIR, CELLCOUNT_CSV), index=False)

    # ------ Figures (NEW names) ------
    sns.set_theme(style="whitegrid", context="talk")

    for gene in gene_cols:
        g = long_df[long_df["gene"] == gene].copy()

        plt.figure(figsize=(14, 5))
        # Keep raw TPM on Y for interpretability in figures
        ax = sns.boxplot(
            data=g, x="zt", y="TPM", hue="Group",
            linewidth=1, fliersize=2,
            palette=GROUP_PALETTE if set(g["Group"]).issubset(GROUP_PALETTE) else "Set2"
        )
        sns.stripplot(
            data=g, x="zt", y="TPM", hue=None, dodge=True, alpha=0.55, size=3
        )
        ax.set_title(f"TPM of {gene} by Time and Group (Treatment + Phase)")
        ax.set_xlabel("ZT (Hours)")
        ax.set_ylabel("TPM Expression")
        plt.legend(title="Group", bbox_to_anchor=(1.02, 1), loc="upper left")
        plt.tight_layout()

        fig_path = os.path.join(OUTPUT_DIR, NEW_FIG_TPL.format(gene=gene))
        plt.savefig(fig_path, dpi=FIG_DPI)
        plt.close()

    # ------ Statistics ------
    # Add log1p for modeling
    long_df["logTPM"] = np.log1p(long_df["TPM"])

    # Store ANOVA tables per gene
    all_posthoc = []

    for gene in gene_cols:
        g = long_df[long_df["gene"] == gene].copy()

        # Type-III ANOVA with treatment*phase + C(zt)
        # Sum-to-zero contrasts to get proper Type-III tests
        model = ols("logTPM ~ C(treatment, Sum) * C(phase, Sum) + C(zt)", data=g).fit()
        aov = anova_lm(model, typ=3)

        aov_path = os.path.join(OUTPUT_DIR, ANOVA_CSV_TPL.format(gene=gene))
        aov.to_csv(aov_path)

        # Planned contrasts: Exercise vs Sedentary within each Phase (controlling ZT)
        for ph in g["phase"].cat.categories if hasattr(g["phase"], "cat") else g["phase"].unique():
            sub = g[g["phase"] == ph].copy()
            # Safety: ensure both groups exist
            if sub["treatment"].nunique() < 2:
                continue

            # Small OLS within phase, controlling for ZT
            m = ols("logTPM ~ C(treatment, Sum) + C(zt)", data=sub).fit()
            aov_ph = anova_lm(m, typ=3)
            p_raw = aov_ph.loc["C(treatment, Sum)", "PR(>F)"]

            # Effect size (Hedges g) on logTPM between Exercise vs Sedentary
            try:
                ex = sub.loc[sub["treatment"].astype(str).str.lower() == "exercise", "logTPM"]
                se = sub.loc[sub["treatment"].astype(str).str.lower() == "sedentary", "logTPM"]
                g_eff = hedges_g(ex.values, se.values)
            except Exception:
                g_eff = np.nan

            all_posthoc.append({
                "Gene": gene,
                "Phase": ph,
                "Contrast": "Exercise - Sedentary (within Phase, adj. for ZT)",
                "p_raw": p_raw,
                "Hedges_g_logTPM": g_eff,
                "n_Exercise": int((sub["treatment"].astype(str).str.lower() == "exercise").sum()),
                "n_Sedentary": int((sub["treatment"].astype(str).str.lower() == "sedentary").sum())
            })

    # BH-FDR across all gene×phase contrasts
    posthoc_df = pd.DataFrame(all_posthoc)
    if not posthoc_df.empty:
        posthoc_df["p_adj_BH"] = multipletests(posthoc_df["p_raw"], method="fdr_bh")[1]
        posthoc_df["significant_0.05"] = posthoc_df["p_adj_BH"] < 0.05
        posthoc_path = os.path.join(OUTPUT_DIR, POSTHOC_CSV)
        posthoc_df.to_csv(posthoc_path, index=False)

    # Console summary
    print("\n=== Saved outputs ===")
    print(f"- Cell counts table: {os.path.join(OUTPUT_DIR, CELLCOUNT_CSV)}")
    for gene in gene_cols:
        print(f"- ANOVA Type-III ({gene}): {os.path.join(OUTPUT_DIR, ANOVA_CSV_TPL.format(gene=gene))}")
        print(f"- Figure ({gene}): {os.path.join(OUTPUT_DIR, NEW_FIG_TPL.format(gene=gene))}")
    if not posthoc_df.empty:
        print(f"- Posthoc contrasts (BH-FDR): {os.path.join(OUTPUT_DIR, POSTHOC_CSV)}")
    else:
        print("- Posthoc contrasts: not available (no valid contrasts found)")

if __name__ == "__main__":
    main()
