import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import ttest_ind

# Load "Alu_EditingIndex_with_groups_A2G_only" or "3UTR_EditingIndex_with_groups_A2G_only"
csv_filename = "Alu_EditingIndex_with_groups_A2G_only.csv"
file_label = csv_filename.split("_")[0]
df = pd.read_csv(csv_filename)
df = df.rename(columns={"A2GEditingIndex": "EditingValue"})

# ----- Setup -----
palette = {"Exercise": "blue", "Sedentary": "orange"}
phases = df['Phase'].unique()

fig, axs = plt.subplots(1, len(phases), figsize=(14, 5), sharey=True)

for i, phase in enumerate(sorted(phases)):
    ax = axs[i]
    sub_df = df[df['Phase'] == phase]

    grouped = sub_df.groupby(['zt', 'treatment'])['EditingValue'].agg(['mean', 'std', 'count']).reset_index()
    grouped['ci95'] = 1.96 * grouped['std'] / np.sqrt(grouped['count'])

    for treatment in grouped['treatment'].unique():
        data = grouped[grouped['treatment'] == treatment]
        ax.errorbar(
            data['zt'], data['mean'], yerr=data['ci95'],
            label=treatment, fmt='o-', capsize=5, color=palette[treatment]
        )

    ax.set_title(f"Phase: {phase}")
    ax.set_xlabel("Zeitgeber Time (ZT)")
    if i == 0:
        ax.set_ylabel("Mean A2G Editing Index")
    ax.grid(True)
    ax.legend(title="Treatment")

plt.suptitle("A2G Editing Index with 95% CI – Split by Phase", fontsize=14)
plt.figtext(0.5, 0.92, f"Dataset: {file_label}", ha='center', fontsize=10, color='gray')


plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig(f"{csv_filename.rsplit('.', 1)[0]}_lineplot_CI_by_phase.png", dpi=300)
plt.show()

# ----- T-test per ZT within each Phase -----
all_results = []

for phase in sorted(phases):
    sub_df = df[df['Phase'] == phase]
    zts = sorted(sub_df['zt'].unique())

    for zt in zts:
        group_data = sub_df[sub_df['zt'] == zt]
        ex = group_data[group_data['treatment'] == 'Exercise']['EditingValue']
        se = group_data[group_data['treatment'] == 'Sedentary']['EditingValue']

        if len(ex) > 1 and len(se) > 1:  # Ensure enough data for t-test
            t_stat, p_val = ttest_ind(ex, se, equal_var=False)
        else:
            p_val = np.nan  # Not enough data

        all_results.append({'Phase': phase, 'ZT': zt, 'p-value': p_val})

# Create a DataFrame with the results
ttest_results_df = pd.DataFrame(all_results)

# Highlight significant results (p < 0.05)
ttest_results_df['significant'] = ttest_results_df['p-value'].apply(
    lambda p: "*" if pd.notnull(p) and p < 0.05 else ""
)

# Display the table
print("\nT-test results by Phase and ZT (with significance marker):")
print(ttest_results_df.to_string(index=False))
