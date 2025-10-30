#!/usr/bin/env python3
"""
Generate LaTeX tables from experiment data for ILS variation study.
"""

import csv
from collections import defaultdict

# Configuration: Define the taxon-gene-tree pairs and ILS level mappings
CONFIGURATIONS = [
    (10000, 1000),
    (15000, 1000),
    (20000, 1000),
    (25000, 1000),
]

# ILS level mapping based on spmax values
ILS_LEVELS = {
    'ILS-L1': 150000,
    'ILS-L2': 200000,
    'ILS-L3': 250000,
    'ILS-L4': 300000,
}

# Read CSV data
data = []
with open('vary-ils.csv', 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        data.append(row)

# Group data by (num_taxa, gene_trees, ils_level)
grouped_data = defaultdict(lambda: defaultdict(list))

for row in data:
    num_taxa = int(row['num-taxa'])
    gene_trees = int(row['gene-trees'])
    spmax = int(row['spmax'])
    
    # Find the ILS level
    ils_level = None
    for level, spmax_value in ILS_LEVELS.items():
        if spmax == spmax_value:
            ils_level = level
            break
    
    if ils_level and (num_taxa, gene_trees) in CONFIGURATIONS:
        key = (num_taxa, gene_trees, ils_level)
        
        # Collect relevant metrics
        try:
            grouped_data[key]['gt-gt'].append(float(row['gt-gt']))
            grouped_data[key]['gt-st'].append(float(row['gt-st']))
            grouped_data[key]['running-time'].append(float(row['running-time-s']))
            grouped_data[key]['cpu-ram'].append(float(row['max-cpu-mb']))
            grouped_data[key]['gpu-vram'].append(float(row['max-gpu-mb']))
        except ValueError as e:
            print(f"Warning: Could not parse data for {key}: {e}")
            continue

# Compute averages
averaged_data = {}
for key, metrics in grouped_data.items():
    averaged_data[key] = {
        'gt-gt': sum(metrics['gt-gt']) / len(metrics['gt-gt']) if metrics['gt-gt'] else None,
        'gt-st': sum(metrics['gt-st']) / len(metrics['gt-st']) if metrics['gt-st'] else None,
        'running-time': sum(metrics['running-time']) / len(metrics['running-time']) if metrics['running-time'] else None,
        'cpu-ram': sum(metrics['cpu-ram']) / len(metrics['cpu-ram']) if metrics['cpu-ram'] else None,
        'gpu-vram': sum(metrics['gpu-vram']) / len(metrics['gpu-vram']) if metrics['gpu-vram'] else None,
    }

# Check for missing data
missing = []
for config in CONFIGURATIONS:
    for ils_level in ILS_LEVELS.keys():
        key = (*config, ils_level)
        if key not in averaged_data or averaged_data[key]['gt-gt'] is None:
            missing.append(key)
            print(f"Warning: Missing data for {config} at {ils_level}")

# Generate LaTeX tables
latex_output = []

# Table 1: Discordance levels
latex_output.append("\\begin{table}[h]")
latex_output.append("\\caption{Discordance levels under different ILS conditions}")
latex_output.append("\\centering")
latex_output.append("\\setlength{\\tabcolsep}{8pt}")
latex_output.append("\\begin{tabular}{l c c}")
latex_output.append("\\toprule")
latex_output.append("\\textbf{ILS Level} & \\textbf{GT--GT Discordance (\\%)} & \\textbf{GT--ST Discordance (\\%)} \\\\")
latex_output.append("\\midrule")

# Average across all configurations for each ILS level
for ils_level in ['ILS-L1', 'ILS-L2', 'ILS-L3', 'ILS-L4']:
    gt_gt_values = []
    gt_st_values = []
    for config in CONFIGURATIONS:
        key = (*config, ils_level)
        if key in averaged_data and averaged_data[key]['gt-gt'] is not None:
            gt_gt_values.append(averaged_data[key]['gt-gt'])
            gt_st_values.append(averaged_data[key]['gt-st'])
    
    if gt_gt_values:
        avg_gt_gt = sum(gt_gt_values) / len(gt_gt_values)
        avg_gt_st = sum(gt_st_values) / len(gt_st_values)
        latex_output.append(f"{ils_level} & {avg_gt_gt:.1f} & {avg_gt_st:.1f} \\\\")
    else:
        latex_output.append(f"{ils_level} & --- & --- \\\\")

latex_output.append("\\bottomrule")
latex_output.append("\\end{tabular}")
latex_output.append("\\end{table}")
latex_output.append("")

# Table 2: Running time
latex_output.append("\\begin{table}[h]")
latex_output.append("\\caption{Running time of STELAR-X under different ILS conditions}")
latex_output.append("\\centering")
latex_output.append("\\setlength{\\tabcolsep}{8pt}")
latex_output.append("\\begin{tabular}{l l c c c c}")
latex_output.append("\\toprule")
latex_output.append("\\multicolumn{2}{c}{} & \\multicolumn{4}{c}{\\textbf{Running Time (s)}} \\\\")
latex_output.append("\\cmidrule(lr){3-6}")
latex_output.append("No. of Taxa & No. of Gene Trees & ILS-L1 & ILS-L2 & ILS-L3 & ILS-L4 \\\\")
latex_output.append("\\midrule")

for config in CONFIGURATIONS:
    taxa, gene_trees = config
    row = f"{taxa} & {gene_trees} & "
    values = []
    for ils_level in ['ILS-L1', 'ILS-L2', 'ILS-L3', 'ILS-L4']:
        key = (*config, ils_level)
        if key in averaged_data and averaged_data[key]['running-time'] is not None:
            values.append(f"{averaged_data[key]['running-time']:.2f}")
        else:
            values.append("---")
    row += " & ".join(values)
    latex_output.append(row + " \\\\")

latex_output.append("\\bottomrule")
latex_output.append("\\end{tabular}")
latex_output.append("\\end{table}")
latex_output.append("")

# Table 3: CPU RAM
latex_output.append("\\begin{table}[h]")
latex_output.append("\\caption{CPU RAM usage of STELAR-X under different ILS conditions}")
latex_output.append("\\centering")
latex_output.append("\\setlength{\\tabcolsep}{8pt}")
latex_output.append("\\begin{tabular}{l l c c c c}")
latex_output.append("\\toprule")
latex_output.append("\\multicolumn{2}{c}{} & \\multicolumn{4}{c}{\\textbf{CPU RAM (MB)}} \\\\")
latex_output.append("\\cmidrule(lr){3-6}")
latex_output.append("No. of Taxa & No. of Gene Trees & ILS-L1 & ILS-L2 & ILS-L3 & ILS-L4 \\\\")
latex_output.append("\\midrule")

for config in CONFIGURATIONS:
    taxa, gene_trees = config
    row = f"{taxa} & {gene_trees} & "
    values = []
    for ils_level in ['ILS-L1', 'ILS-L2', 'ILS-L3', 'ILS-L4']:
        key = (*config, ils_level)
        if key in averaged_data and averaged_data[key]['cpu-ram'] is not None:
            values.append(f"{averaged_data[key]['cpu-ram']:.2f}")
        else:
            values.append("---")
    row += " & ".join(values)
    latex_output.append(row + " \\\\")

latex_output.append("\\bottomrule")
latex_output.append("\\end{tabular}")
latex_output.append("\\end{table}")
latex_output.append("")

# Table 4: GPU VRAM
latex_output.append("\\begin{table}[h]")
latex_output.append("\\caption{GPU VRAM usage of STELAR-X under different ILS conditions}")
latex_output.append("\\centering")
latex_output.append("\\setlength{\\tabcolsep}{8pt}")
latex_output.append("\\begin{tabular}{l l c c c c}")
latex_output.append("\\toprule")
latex_output.append("\\multicolumn{2}{c}{} & \\multicolumn{4}{c}{\\textbf{GPU VRAM (MB)}} \\\\")
latex_output.append("\\cmidrule(lr){3-6}")
latex_output.append("No. of Taxa & No. of Gene Trees & ILS-L1 & ILS-L2 & ILS-L3 & ILS-L4 \\\\")
latex_output.append("\\midrule")

for config in CONFIGURATIONS:
    taxa, gene_trees = config
    row = f"{taxa} & {gene_trees} & "
    values = []
    for ils_level in ['ILS-L1', 'ILS-L2', 'ILS-L3', 'ILS-L4']:
        key = (*config, ils_level)
        if key in averaged_data and averaged_data[key]['gpu-vram'] is not None:
            values.append(f"{averaged_data[key]['gpu-vram']:.2f}")
        else:
            values.append("---")
    row += " & ".join(values)
    latex_output.append(row + " \\\\")

latex_output.append("\\bottomrule")
latex_output.append("\\end{tabular}")
latex_output.append("\\end{table}")

# Write to file
with open('ils_tables.tex', 'w') as f:
    f.write('\n'.join(latex_output))

print("LaTeX tables generated successfully in ils_tables.tex")
if missing:
    print(f"\nMissing data for {len(missing)} configurations")
