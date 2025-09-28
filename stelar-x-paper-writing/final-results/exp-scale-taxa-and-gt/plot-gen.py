import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

def detect_scaling_type(df):
    """Detect whether this is scaling taxa or scaling gene trees"""
    num_taxa_unique = len(df['num-taxa'].unique())
    gene_trees_unique = len(df['gene-trees'].unique())
    
    if num_taxa_unique > gene_trees_unique:
        return "taxa"  # Scaling number of taxa
    else:
        return "gene_trees"  # Scaling number of gene trees

def generate_scaling_plots(csv_file):
    # Read the CSV file
    df = pd.read_csv(csv_file)
    
    # Detect scaling type
    scaling_type = detect_scaling_type(df)
    
    # Extract relevant data
    num_taxa = df['num-taxa'].values
    gene_trees = df['gene-trees'].values
    running_time = df['running-time-s'].values / 3600  # Convert to hours
    cpu_ram = df['max-cpu-mb'].values / 1000  # Convert to GB
    gpu_vram = df['max-gpu-mb'].values / 1000  # Convert to GB
    
    # Determine x-axis data and labels based on scaling type
    if scaling_type == "taxa":
        x_data = num_taxa / 1000  # Convert to thousands
        x_label = 'Number of Taxa (×1000)'
        x_ticks = num_taxa / 1000
        x_tick_labels = [f'{int(x)}' for x in num_taxa/1000]
        plot_title_prefix = "Taxa Scaling"
    else:
        x_data = gene_trees / 1000  # Convert to thousands
        x_label = 'Number of Gene Trees (×1000)'
        x_ticks = gene_trees / 1000
        x_tick_labels = [f'{int(x)}' for x in gene_trees/1000]
        plot_title_prefix = "Gene Trees Scaling"
    
    # Set up the figure with academic styling
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Define colors for academic papers (distinct and colorblind-friendly)
    color_runtime = '#2E86C1'  # Blue
    color_cpu = '#E74C3C'      # Red
    color_gpu = '#28B463'      # Green
    
    # Left subplot: Running Time vs Variable
    ax1.plot(x_data, running_time, 'o-', 
             color=color_runtime, linewidth=2.5, markersize=8, 
             markerfacecolor=color_runtime, markeredgecolor='white', 
             markeredgewidth=1.5, label='Running Time')
    
    ax1.set_xlabel(x_label, fontsize=12, fontweight='bold')
    ax1.set_ylabel('Running Time (hours)', fontsize=12, fontweight='bold')
    ax1.set_title('Running Time', fontsize=14, fontweight='bold', pad=20)
    ax1.grid(True, alpha=0.3)
    ax1.tick_params(axis='both', which='major', labelsize=10)
    
    # Set x-axis ticks to match the data points
    ax1.set_xticks(x_ticks)
    ax1.set_xticklabels(x_tick_labels)
    
    # Add some padding to y-axis
    ax1.set_ylim(bottom=0)
    y_max = max(running_time) * 1.1
    ax1.set_ylim(top=y_max)
    
    # Right subplot: Memory Usage vs Variable
    ax2.plot(x_data, cpu_ram, 'o-', 
             color=color_cpu, linewidth=2.5, markersize=8,
             markerfacecolor=color_cpu, markeredgecolor='white', 
             markeredgewidth=1.5, label='CPU RAM')
    
    ax2.plot(x_data, gpu_vram, 's-', 
             color=color_gpu, linewidth=2.5, markersize=8,
             markerfacecolor=color_gpu, markeredgecolor='white', 
             markeredgewidth=1.5, label='GPU VRAM')
    
    ax2.set_xlabel(x_label, fontsize=12, fontweight='bold')
    ax2.set_ylabel('Memory Usage (GB)', fontsize=12, fontweight='bold')
    ax2.set_title('Memory Usage', fontsize=14, fontweight='bold', pad=20)
    ax2.grid(True, alpha=0.3)
    ax2.tick_params(axis='both', which='major', labelsize=10)
    ax2.legend(loc='upper left', fontsize=11, framealpha=0.9)
    
    # Set x-axis ticks to match the data points
    ax2.set_xticks(x_ticks)
    ax2.set_xticklabels(x_tick_labels)
    
    # Add some padding to y-axis
    ax2.set_ylim(bottom=0)
    y_max_mem = max(max(cpu_ram), max(gpu_vram)) * 1.1
    ax2.set_ylim(top=y_max_mem)
    
    # Adjust layout to prevent overlap
    plt.tight_layout(pad=3.0)

    # Generate appropriate filename based on CSV input
    exp_type = os.path.basename(csv_file).replace('result-scale-', '').replace('.csv', '')
    name = f"stelar-scale-{exp_type}"
    
    # Save the plot in high resolution for academic use
    plt.savefig(f'{name}.png', dpi=600, bbox_inches='tight', 
                facecolor='white', edgecolor='none')
    plt.savefig(f'{name}.pdf', bbox_inches='tight', 
                facecolor='white', edgecolor='none')
    
    # Show the plot
    plt.show()
    
    print(f"Plots saved as '{name}.png' and '{name}.pdf'")
    
    # Print some statistics based on scaling type
    print(f"\nData Summary:")
    print(f"Scaling type: {scaling_type}")
    if scaling_type == "taxa":
        print(f"Taxa range: {min(num_taxa):,} - {max(num_taxa):,}")
        print(f"Gene trees (constant): {gene_trees[0]:,}")
    else:
        print(f"Gene trees range: {min(gene_trees):,} - {max(gene_trees):,}")
        print(f"Taxa (constant): {num_taxa[0]:,}")
    print(f"Runtime range: {min(running_time):.1f} - {max(running_time):.1f} hours")
    print(f"CPU RAM range: {min(cpu_ram):.1f} - {max(cpu_ram):.1f} GB")
    print(f"GPU VRAM range: {min(gpu_vram):.1f} - {max(gpu_vram):.1f} GB")
    
    return fig, scaling_type

# def generate_alternative_plot(csv_file):
#     """Alternative version with log scale for better visualization of exponential growth"""
#     df = pd.read_csv(csv_file)
    
#     num_taxa = df['num-taxa'].values
#     running_time = df['running-time-s'].values / 3600  # Convert to hours
#     cpu_ram = df['max-cpu-mb'].values
#     gpu_vram = df['max-gpu-mb'].values
    
#     plt.style.use('seaborn-v0_8-whitegrid')
#     fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
#     # Colors
#     color_runtime = '#2E86C1'
#     color_cpu = '#E74C3C'
#     color_gpu = '#28B463'
    
#     # Left subplot with log scale
#     ax1.semilogy(num_taxa/1000, running_time, 'o-', 
#                  color=color_runtime, linewidth=2.5, markersize=8, 
#                  markerfacecolor=color_runtime, markeredgecolor='white', 
#                  markeredgewidth=1.5, label='Running Time')
    
#     ax1.set_xlabel('Number of Taxa (×1000)', fontsize=12, fontweight='bold')
#     ax1.set_ylabel('Running Time (hours, log scale)', fontsize=12, fontweight='bold')
#     ax1.set_title('Running Time (Log Scale)', fontsize=14, fontweight='bold', pad=20)
#     ax1.grid(True, alpha=0.3)
#     ax1.tick_params(axis='both', which='major', labelsize=10)
#     ax1.set_xticks(num_taxa/1000)
#     ax1.set_xticklabels([f'{int(x)}' for x in num_taxa/1000])
    
#     # Right subplot
#     ax2.plot(num_taxa/1000, cpu_ram/1000, 'o-', 
#              color=color_cpu, linewidth=2.5, markersize=8,
#              markerfacecolor=color_cpu, markeredgecolor='white', 
#              markeredgewidth=1.5, label='CPU RAM')
    
#     ax2.plot(num_taxa/1000, gpu_vram/1000, 's-', 
#              color=color_gpu, linewidth=2.5, markersize=8,
#              markerfacecolor=color_gpu, markeredgecolor='white', 
#              markeredgewidth=1.5, label='GPU VRAM')
    
#     ax2.set_xlabel('Number of Taxa (×1000)', fontsize=12, fontweight='bold')
#     ax2.set_ylabel('Memory Usage (GB)', fontsize=12, fontweight='bold')
#     ax2.set_title('Memory Usage', fontsize=14, fontweight='bold', pad=20)
#     ax2.grid(True, alpha=0.3)
#     ax2.tick_params(axis='both', which='major', labelsize=10)
#     ax2.legend(loc='upper left', fontsize=11, framealpha=0.9)
#     ax2.set_xticks(num_taxa/1000)
#     ax2.set_xticklabels([f'{int(x)}' for x in num_taxa/1000])
#     ax2.set_ylim(bottom=0)
    
#     plt.tight_layout(pad=3.0)
#     plt.savefig('stelar_scalability_logscale.png', dpi=300, bbox_inches='tight', 
#                 facecolor='white', edgecolor='none')
#     plt.savefig('stelar_scalability_logscale.pdf', bbox_inches='tight', 
#                 facecolor='white', edgecolor='none')
    
#     plt.show()
#     print("Log-scale plots saved as 'stelar_scalability_logscale.png' and 'stelar_scalability_logscale.pdf'")
    
#     return fig

def get_available_csv_files():
    """Get list of available result CSV files"""
    csv_files = []
    for file in os.listdir('.'):
        if file.startswith('result-scale-') and file.endswith('.csv'):
            csv_files.append(file)
    return sorted(csv_files)

def interactive_file_selection():
    """Interactive file selection with good UX"""
    csv_files = get_available_csv_files()
    
    if not csv_files:
        print("No result CSV files found in current directory.")
        return None
    
    if len(csv_files) == 1:
        print(f"Found CSV file: {csv_files[0]}")
        return csv_files[0]
    
    print("Available CSV files:")
    for i, file in enumerate(csv_files, 1):
        # Extract experiment type from filename for better display
        exp_type = file.replace('result-scale-', '').replace('.csv', '')
        print(f"{i}. {file} (scaling {exp_type})")
    
    while True:
        try:
            choice = input(f"\nSelect file (1-{len(csv_files)}) or press Enter for default [{csv_files[0]}]: ").strip()
            if not choice:
                return csv_files[0]
            choice_idx = int(choice) - 1
            if 0 <= choice_idx < len(csv_files):
                return csv_files[choice_idx]
            else:
                print(f"Please enter a number between 1 and {len(csv_files)}")
        except ValueError:
            print("Please enter a valid number")

if __name__ == "__main__":
    # Handle command line arguments or interactive selection
    if len(sys.argv) > 1:
        csv_file = sys.argv[1]
        if not os.path.exists(csv_file):
            print(f"Error: File '{csv_file}' not found")
            sys.exit(1)
    else:
        csv_file = interactive_file_selection()
        if csv_file is None:
            sys.exit(1)
    
    print(f"\nProcessing: {csv_file}")
    print("Generating scaling plots...")
    fig, scaling_type = generate_scaling_plots(csv_file)
