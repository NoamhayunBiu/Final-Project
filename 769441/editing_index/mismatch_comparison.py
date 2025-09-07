import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

# Load the file
input_file = "alu_source_file.csv"  # Change to your file name if needed
df = pd.read_csv(input_file)

# Extract base name (without extension)
base_name = os.path.splitext(os.path.basename(input_file))[0]

# The columns we want to compare
editing_indices = [
    'A2CEditingIndex',
    'A2GEditingIndex',
    'A2TEditingIndex',
    'C2AEditingIndex',
    'C2GEditingIndex',
    'C2TEditingIndex'
]

# Convert the table to long format – suitable for boxplot
df_long = df[editing_indices].melt(var_name='MismatchType', value_name='EditingPercentage')

# Boxplot
plt.figure(figsize=(10, 6))
sns.boxplot(data=df_long, x='MismatchType', y='EditingPercentage')
plt.title(f'Editing Percentage by Mismatch Type ({base_name})')
plt.ylabel('Editing Percentage')
plt.xlabel('Mismatch Type')
plt.tight_layout()

# Create output file name matching the input file
output_file = f"{base_name}_editing_index_distribution.png"

# Save the plot
plt.savefig(output_file, dpi=300)
print(f"Saved plot to {output_file}")

plt.show()
