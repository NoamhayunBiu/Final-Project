import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import kruskal

def create_editing_plot(csv_path, region_label, output_filename):
    # Read the data
    df = pd.read_csv(csv_path)
    df.columns = df.columns.str.strip()  # Clean whitespace

    # Filter relevant groups
    df_filtered = df[df["treatment (eps)"].isin(["Rest", "EPS 0h", "EPS 3h"])].copy()
    df_filtered["group"] = "Control – " + df_filtered["treatment (eps)"]

    # Statistics
    groups = df_filtered["group"].unique()
    data_groups = [df_filtered[df_filtered["group"] == g]["A2GEditingIndex"] for g in groups]
    stat, p = kruskal(*data_groups)
    p_text = f"Kruskal-Wallis test:\nstatistic = {stat:.3f}\np-value = {p:.4f}"

    # Figure
    fig, (ax_plot, ax_text) = plt.subplots(ncols=2, figsize=(12, 6), gridspec_kw={'width_ratios': [4, 1]})

    # Boxplot
    sns.boxplot(data=df_filtered, x="group", y="A2GEditingIndex", palette="Set2", ax=ax_plot)
    sns.stripplot(data=df_filtered, x="group", y="A2GEditingIndex", color='black', size=4, jitter=True, alpha=0.7, ax=ax_plot)

    # Titles
    ax_plot.set_title("A2G Editing Index by Group", fontsize=14, weight='bold', pad=30)
    ax_plot.text(0.5, 1.05, region_label, transform=ax_plot.transAxes, fontsize=12, ha='center')
    ax_plot.set_ylabel("A2G Editing Index")
    ax_plot.set_xlabel("Group")
    ax_plot.tick_params(axis='x', rotation=15)
    ax_plot.grid(axis='y', linestyle='--', alpha=0.5)

    # Statistical text on the side
    ax_text.axis('off')
    ax_text.text(0, 0.9, p_text, fontsize=12, verticalalignment='top',
                 bbox=dict(facecolor='white', edgecolor='gray', boxstyle='round,pad=0.5'))

    # Save
    plt.tight_layout()
    plt.savefig(output_filename, dpi=300)
    plt.close()

    print(f"{region_label}: Plot saved as {output_filename}")
    print(f"Statistic = {stat:.3f}, p-value = {p:.4f}\n")


# Run the function on both files
create_editing_plot("3utr_averaged_notKD_only.csv", "3′UTR", "editing_index_3utr.png")
create_editing_plot("alu_averaged_notKD_only.csv", "ALU", "editing_index_alu.png")
