import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

csv_filename = "3UTR_EditingIndex.csv"  # 
file_label = csv_filename.split("_")[0]
df = pd.read_csv(csv_filename)

# Extract base name from the file (without .csv extension)
base_filename = os.path.splitext(os.path.basename(csv_filename))[0]

# Convert the table to long format ("long")
df_long = df.melt(id_vars=["Sample"],
                  value_vars=["A2CEditingIndex", "A2GEditingIndex", "A2TEditingIndex",
                              "C2AEditingIndex", "C2GEditingIndex", "C2TEditingIndex"],
                  var_name="EditingType", value_name="EditingValue")

# Plot: comparison of editing levels between types
plt.figure(figsize=(10, 6))
sns.boxplot(x="EditingType", y="EditingValue", data=df_long)
sns.stripplot(x="EditingType", y="EditingValue", data=df_long, color='black', alpha=0.3, jitter=True)

plt.suptitle("Comparison of Potential RNA Editing Types", fontsize=14)
plt.figtext(0.5, 0.92, f"Dataset: {file_label}", ha='center', fontsize=10, color='gray')
plt.ylabel("Editing Index Value")
plt.xlabel("Editing Type")
plt.xticks(rotation=45)
plt.tight_layout()

# Save with the file name
plt.savefig(f"{base_filename}_editing_comparison_plot.png", dpi=300)

plt.show()


csv_filename = "Alu_EditingIndex.csv"  # ← Change to your file name here
file_label = csv_filename.split("_")[0]
df = pd.read_csv(csv_filename)

# Extract base name from the file (without .csv extension)
base_filename = os.path.splitext(os.path.basename(csv_filename))[0]

# Convert the table to long format ("long")
df_long = df.melt(id_vars=["Sample"],
                  value_vars=["A2CEditingIndex", "A2GEditingIndex", "A2TEditingIndex",
                              "C2AEditingIndex", "C2GEditingIndex", "C2TEditingIndex"],
                  var_name="EditingType", value_name="EditingValue")

# Plot: comparison of editing levels between types
plt.figure(figsize=(10, 6))
sns.boxplot(x="EditingType", y="EditingValue", data=df_long)
sns.stripplot(x="EditingType", y="EditingValue", data=df_long, color='black', alpha=0.3, jitter=True)

plt.suptitle("Comparison of Potential RNA Editing Types", fontsize=14)
plt.figtext(0.5, 0.92, f"Dataset: {file_label}", ha='center', fontsize=10, color='gray')
plt.ylabel("Editing Index Value")
plt.xlabel("Editing Type")
plt.xticks(rotation=45)
plt.tight_layout()

# Save with the file name
plt.savefig(f"{base_filename}_editing_comparison_plot.png", dpi=300)

plt.show()
