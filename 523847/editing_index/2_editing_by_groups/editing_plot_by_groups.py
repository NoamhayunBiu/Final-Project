import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

# Load "Alu_EditingIndex_with_groups" or "3UTR_EditingIndex_with_groups" 
csv_filename = "3UTR_EditingIndex_with_groups.csv"  # Change if needed
file_label = csv_filename.split("_")[0]
df = pd.read_csv(csv_filename)

# Extract base name
base_filename = os.path.splitext(os.path.basename(csv_filename))[0]

# Convert table to long format
df_long = df.melt(
    id_vars=["Sample", "Phase", "treatment", "zt"],
    value_vars=[
        "A2CEditingIndex", "A2GEditingIndex", "A2TEditingIndex",
        "C2AEditingIndex", "C2GEditingIndex", "C2TEditingIndex"
    ],
    var_name="EditingType",
    value_name="EditingValue"
)

# ========= Combined plot: compare Phase by treatment =========
plt.figure(figsize=(12, 6))
sns.boxplot(data=df_long, x="Phase", y="EditingValue", hue="treatment")
sns.stripplot(data=df_long, x="Phase", y="EditingValue", hue="treatment", 
              dodge=True, color='black', alpha=0.2)

plt.suptitle("RNA Editing by Phase and Treatment", fontsize=14)
plt.figtext(0.5, 0.92, f"Dataset: {file_label}", ha='center', fontsize=10, color='gray')
             
plt.ylabel("Editing Index Value")
plt.xlabel("Phase")
plt.tight_layout()
plt.legend(title="Treatment", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.savefig(f"{base_filename}_editing_by_phase_and_treatment.png", dpi=300)
plt.show()

# ========= Plot 3: Exercise + Phase by ZT =========
exercise_only = df_long[df_long["treatment"] == "Exercise"]
g = sns.catplot(
    data=exercise_only,
    x="EditingType", y="EditingValue",
    hue="Phase", col="zt",
    kind="box", col_wrap=4, height=4, aspect=1.2
)
g.fig.subplots_adjust(top=0.85)
plt.suptitle("RNA Editing in Exercise Mice by Phase and ZT", fontsize=14)
plt.figtext(0.5, 0.92, f"Dataset: {file_label}", ha='center', fontsize=10, color='gray')
g.set_xticklabels(rotation=45)
g.savefig(f"{base_filename}_exercise_phase_by_zt.png", dpi=300)
plt.show()


csv_filename = "Alu_EditingIndex_with_groups.csv"  # Change if needed
file_label = csv_filename.split("_")[0]
df = pd.read_csv(csv_filename)

# Extract base name
base_filename = os.path.splitext(os.path.basename(csv_filename))[0]

# Convert table to long format
df_long = df.melt(
    id_vars=["Sample", "Phase", "treatment", "zt"],
    value_vars=[
        "A2CEditingIndex", "A2GEditingIndex", "A2TEditingIndex",
        "C2AEditingIndex", "C2GEditingIndex", "C2TEditingIndex"
    ],
    var_name="EditingType",
    value_name="EditingValue"
)

# ========= Combined plot: compare Phase by treatment =========
plt.figure(figsize=(12, 6))
sns.boxplot(data=df_long, x="Phase", y="EditingValue", hue="treatment")
sns.stripplot(data=df_long, x="Phase", y="EditingValue", hue="treatment", 
              dodge=True, color='black', alpha=0.2)

plt.suptitle("RNA Editing by Phase and Treatment", fontsize=14)
plt.figtext(0.5, 0.92, f"Dataset: {file_label}", ha='center', fontsize=10, color='gray')
             
plt.ylabel("Editing Index Value")
plt.xlabel("Phase")
plt.tight_layout()
plt.legend(title="Treatment", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.savefig(f"{base_filename}_editing_by_phase_and_treatment.png", dpi=300)
plt.show()

# ========= Plot 3: Exercise + Phase by ZT =========
exercise_only = df_long[df_long["treatment"] == "Exercise"]
g = sns.catplot(
    data=exercise_only,
    x="EditingType", y="EditingValue",
    hue="Phase", col="zt",
    kind="box", col_wrap=4, height=4, aspect=1.2
)
g.fig.subplots_adjust(top=0.85)
plt.suptitle("RNA Editing in Exercise Mice by Phase and ZT", fontsize=14)
plt.figtext(0.5, 0.92, f"Dataset: {file_label}", ha='center', fontsize=10, color='gray')
g.set_xticklabels(rotation=45)
g.savefig(f"{base_filename}_exercise_phase_by_zt.png", dpi=300)
plt.show()
