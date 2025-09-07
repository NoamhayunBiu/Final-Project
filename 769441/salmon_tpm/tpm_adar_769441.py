import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.formula.api import ols
import statsmodels.api as sm
from scipy.stats import kruskal

# Load and clean data
df = pd.read_csv("source_file.csv")
df.columns = df.columns.str.strip()

# Rename columns
df.rename(columns={
    'treatment (lna)': 'treatment_ina',
    'treatment (eps)': 'treatment_eps'
}, inplace=True)

# Melt dataframe
long_df = pd.melt(
    df,
    id_vars=['Sample Name', 'treatment_ina', 'treatment_eps'],
    value_vars=['Adar', 'Adarb1', 'Adarb2'],
    var_name='gene',
    value_name='TPM'
)

# Strip whitespaces and define category orders
long_df['treatment_eps'] = long_df['treatment_eps'].str.strip()
long_df['treatment_ina'] = long_df['treatment_ina'].str.strip()

eps_order = ['Rest', 'EPS 0h', 'EPS 3h']
long_df['treatment_eps'] = pd.Categorical(long_df['treatment_eps'], categories=eps_order, ordered=True)

# ✅ Filter only 'Negative control' (no Tug1 KD)
long_df = long_df[long_df['treatment_ina'] == 'Negative control']

# Plot
sns.set(style='whitegrid')

for gene in ['Adar', 'Adarb1', 'Adarb2']:
    subset = long_df[long_df['gene'] == gene]

    plt.figure(figsize=(8, 6))
    ax = sns.boxplot(
        data=subset,
        x='treatment_eps',
        y='TPM',
        palette='Set2'
    )

    plt.title(f'{gene} TPM Expression', fontsize=14)
    plt.xlabel("EPS Treatment")
    plt.ylabel("TPM Expression")

    # --- סטטיסטיקה ---
    groups = [group['TPM'].values for _, group in subset.groupby('treatment_eps')]

    # בדיקת תנאים ל-ANOVA
    unique_counts = subset['treatment_eps'].value_counts()
    if all(unique_counts >= 3):
        try:
            model = ols('TPM ~ C(treatment_eps)', data=subset).fit()
            anova_table = sm.stats.anova_lm(model, typ=2)
            p_val = anova_table['PR(>F)'][0]
            test_name = 'One-way ANOVA'
        except:
            stat, p_val = kruskal(*groups)
            test_name = 'Kruskal-Wallis'
    else:
        stat, p_val = kruskal(*groups)
        test_name = 'Kruskal-Wallis'

    # --- תיבת p-value ---
    ax.text(1.05, 0.6, f'{test_name}\np = {p_val:.4f}', transform=ax.transAxes,
            fontsize=10, fontweight='bold', verticalalignment='top', ha='center',
            bbox=dict(boxstyle='round', facecolor='lightgrey', alpha=0.8))

    plt.tight_layout()
    plt.savefig(f"{gene}_TPM_NoKD_by_eps.png", dpi=300)
    plt.show()
