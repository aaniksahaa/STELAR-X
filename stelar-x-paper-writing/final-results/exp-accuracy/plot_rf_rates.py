#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# Set style for better-looking plots with white background
plt.style.use('default')  # Use default matplotlib style for clean white background
sns.set_palette("husl")

# Read the CSV file
df = pd.read_csv('stat-standard.csv')

# Filter for 100-taxon and 200-taxon data only
filtered_df = df[df['folder'].isin(['100-taxon', '200-taxon'])].copy()

print("Data overview:")
print(f"Total rows: {len(filtered_df)}")
print(f"Unique algorithms: {filtered_df['alg'].unique()}")
print(f"Unique taxon numbers: {filtered_df['folder'].unique()}")
print("\nData summary by algorithm and taxon:")
print(filtered_df.groupby(['folder', 'alg']).agg({
    'rf-rate': ['count', 'mean', 'std']
}).round(4))

# Calculate average RF rates for each algorithm-taxon combination
avg_rf = filtered_df.groupby(['folder', 'alg'])['rf-rate'].mean().reset_index()
std_rf = filtered_df.groupby(['folder', 'alg'])['rf-rate'].std().reset_index()

# Merge average and std
plot_data = avg_rf.merge(std_rf, on=['folder', 'alg'], suffixes=('_mean', '_std'))
plot_data['rf-rate_std'] = plot_data['rf-rate_std'].fillna(0)  # Handle cases with single data point

print("\nAverage RF rates:")
print(plot_data)

# Create subplot figure with white background
fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=True, facecolor='white')
fig.patch.set_facecolor('white')

# Define colors for each algorithm
algorithms = ['stelar', 'wqfmtree', 'astral', 'treeqmc']
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
color_map = dict(zip(algorithms, colors))

# Plot data for each taxon number
taxon_folders = ['100-taxon', '200-taxon']
taxon_labels = ['100 Taxa', '200 Taxa']

for i, (folder, label) in enumerate(zip(taxon_folders, taxon_labels)):
    ax = axes[i]
    
    # Get data for this taxon number
    taxon_data = plot_data[plot_data['folder'] == folder].copy()
    
    # Ensure all algorithms are present (fill with 0 if missing)
    for alg in algorithms:
        if alg not in taxon_data['alg'].values:
            new_row = pd.DataFrame({
                'folder': [folder],
                'alg': [alg],
                'rf-rate_mean': [0],
                'rf-rate_std': [0]
            })
            taxon_data = pd.concat([taxon_data, new_row], ignore_index=True)
    
    # Sort by algorithm order
    taxon_data['alg_order'] = taxon_data['alg'].map({alg: i for i, alg in enumerate(algorithms)})
    taxon_data = taxon_data.sort_values('alg_order')
    
    # Create bar plot
    x_pos = np.arange(len(algorithms))
    bars = ax.bar(x_pos, taxon_data['rf-rate_mean'], 
                  yerr=taxon_data['rf-rate_std'],
                  capsize=5,
                  color=[color_map[alg] for alg in taxon_data['alg']],
                  alpha=0.8,
                  edgecolor='black',
                  linewidth=0.5)
    
    # Customize subplot
    ax.set_title(f'{label}', fontsize=14, fontweight='bold')
    ax.set_xlabel('Algorithm', fontsize=12)
    if i == 0:
        ax.set_ylabel('Average RF Rate', fontsize=12)
    
    # Set x-axis labels
    ax.set_xticks(x_pos)
    ax.set_xticklabels([alg.upper() if alg != 'wqfmtree' else 'WQFM' for alg in algorithms], 
                       rotation=45, ha='right')
    
    # Add value labels on top of bars
    for j, (bar, val, std_val) in enumerate(zip(bars, taxon_data['rf-rate_mean'], taxon_data['rf-rate_std'])):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + std_val + 0.001,
                f'{val:.4f}',
                ha='center', va='bottom', fontsize=10)
    
    # Customize grid and background
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_axisbelow(True)
    ax.set_facecolor('white')

# Overall plot formatting
plt.suptitle('Average RF Rates by Algorithm and Taxon Number', fontsize=16, fontweight='bold', y=0.98)
plt.tight_layout()

# Save the plot
plt.savefig('rf_rates_comparison.png', dpi=600, bbox_inches='tight', facecolor='white')
plt.savefig('rf_rates_comparison.pdf', bbox_inches='tight', facecolor='white')

print("\nPlot saved as 'rf_rates_comparison.png' and 'rf_rates_comparison.pdf'")

# Show the plot
plt.show()
