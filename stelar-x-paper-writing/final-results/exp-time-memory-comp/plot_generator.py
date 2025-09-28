#!/usr/bin/env python3
"""
Academic Performance Comparison Plot Generator
Generates publication-quality plots comparing STELAR and ASTRAL performance
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from pathlib import Path

# Set style for academic publications
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("bright")

# Configure matplotlib for publication quality
plt.rcParams.update({
    'font.size': 12,
    'font.family': 'DejaVu Sans',
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'figure.titlesize': 18,
    'lines.linewidth': 2,
    'grid.alpha': 0.3
})

def load_data():
    """Load and prepare the performance data from CSV files"""
    stelar_df = pd.read_csv('stelar.csv')
    astral_df = pd.read_csv('astral.csv')
    
    # Combine the datasets
    combined_df = pd.concat([stelar_df, astral_df], ignore_index=True)
    
    return combined_df

def create_time_comparison_plot(df):
    """Create time comparison bar chart between STELAR and ASTRAL"""
    
    # Prepare data
    stelar_data = df[df['alg'] == 'stelar'].copy()
    astral_data = df[df['alg'] == 'astral'].copy()
    
    # Sort by num-taxa for consistent ordering
    stelar_data = stelar_data.sort_values('num-taxa')
    astral_data = astral_data.sort_values('num-taxa')
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Set up bar positions
    x = np.arange(len(stelar_data))
    width = 0.35
    
    # Create bars
    bars1 = ax.bar(x - width/2, stelar_data['running-time-s'], width, 
                   label='STELAR', color='#2E86AB', alpha=0.8, edgecolor='black', linewidth=0.5)
    bars2 = ax.bar(x + width/2, astral_data['running-time-s'], width,
                   label='ASTRAL', color='#A23B72', alpha=0.8, edgecolor='black', linewidth=0.5)
    
    # Customize the plot
    ax.set_xlabel('Number of Taxa', fontweight='bold')
    ax.set_ylabel('Running Time (seconds)', fontweight='bold')
    ax.set_title('Performance Comparison: STELAR vs ASTRAL\nRunning Time Analysis', 
                 fontweight='bold', pad=20)
    ax.set_xticks(x)
    ax.set_xticklabels([f'{int(taxa)}' for taxa in stelar_data['num-taxa']])
    
    # Add value labels on bars
    def add_value_labels(bars, values):
        for bar, value in zip(bars, values):
            height = bar.get_height()
            if height > 1000:
                label = f'{height/3600:.1f}h' if height > 3600 else f'{height:.0f}s'
            else:
                label = f'{height:.1f}s'
            ax.text(bar.get_x() + bar.get_width()/2., height + height*0.01,
                   label, ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    add_value_labels(bars1, stelar_data['running-time-s'])
    add_value_labels(bars2, astral_data['running-time-s'])
    
    # Set y-axis to log scale for better visualization
    ax.set_yscale('log')
    ax.set_ylabel('Running Time (seconds, log scale)', fontweight='bold')
    
    # Customize legend
    ax.legend(loc='upper left', frameon=True, fancybox=True, shadow=True)
    
    # Add grid for better readability
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save the plot
    plt.savefig('time_comparison.png', dpi=300, bbox_inches='tight', 
                facecolor='white', edgecolor='none')
    plt.savefig('time_comparison.pdf', bbox_inches='tight', 
                facecolor='white', edgecolor='none')
    
    print("Time comparison plot saved as 'time_comparison.png' and 'time_comparison.pdf'")
    
    return fig

def create_memory_comparison_plot(df):
    """Create memory comparison with separate subplots for CPU RAM and GPU VRAM"""
    
    # Prepare data
    stelar_data = df[df['alg'] == 'stelar'].copy()
    astral_data = df[df['alg'] == 'astral'].copy()
    
    # Sort by num-taxa for consistent ordering
    stelar_data = stelar_data.sort_values('num-taxa')
    astral_data = astral_data.sort_values('num-taxa')
    
    # Convert MB to GB for better readability
    stelar_cpu = stelar_data['max-cpu-mb'] / 1024
    stelar_gpu = stelar_data['max-gpu-mb'] / 1024
    astral_cpu = astral_data['max-cpu-mb'] / 1024  
    astral_gpu = astral_data['max-gpu-mb'] / 1024
    
    # Create figure with two subplots side by side
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # Set up bar positions
    x = np.arange(len(stelar_data))
    width = 0.35
    
    # === LEFT SUBPLOT: CPU RAM ===
    bars1_cpu = ax1.bar(x - width/2, stelar_cpu, width, 
                        label='STELAR', color='#2E86AB', alpha=0.8, 
                        edgecolor='black', linewidth=0.5)
    bars2_cpu = ax1.bar(x + width/2, astral_cpu, width,
                        label='ASTRAL', color='#A23B72', alpha=0.8,
                        edgecolor='black', linewidth=0.5)
    
    # Customize CPU subplot
    ax1.set_xlabel('Number of Taxa', fontweight='bold')
    ax1.set_ylabel('CPU RAM Usage (GB)', fontweight='bold')
    ax1.set_title('CPU RAM Comparison', fontweight='bold', pad=15)
    ax1.set_xticks(x)
    ax1.set_xticklabels([f'{int(taxa)}' for taxa in stelar_data['num-taxa']])
    
    # Add value labels on CPU bars
    def add_value_labels(ax, bars, values):
        for bar, value in zip(bars, values):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + height*0.01,
                   f'{height:.1f}GB', ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    add_value_labels(ax1, bars1_cpu, stelar_cpu)
    add_value_labels(ax1, bars2_cpu, astral_cpu)
    
    ax1.legend(loc='upper left', frameon=True, fancybox=True, shadow=True)
    ax1.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    
    # === RIGHT SUBPLOT: GPU VRAM ===
    bars1_gpu = ax2.bar(x - width/2, stelar_gpu, width, 
                        label='STELAR', color='#2E86AB', alpha=0.8, 
                        edgecolor='black', linewidth=0.5)
    bars2_gpu = ax2.bar(x + width/2, astral_gpu, width,
                        label='ASTRAL', color='#A23B72', alpha=0.8,
                        edgecolor='black', linewidth=0.5)
    
    # Customize GPU subplot
    ax2.set_xlabel('Number of Taxa', fontweight='bold')
    ax2.set_ylabel('GPU VRAM Usage (GB)', fontweight='bold')
    ax2.set_title('GPU VRAM Comparison', fontweight='bold', pad=15)
    ax2.set_xticks(x)
    ax2.set_xticklabels([f'{int(taxa)}' for taxa in stelar_data['num-taxa']])
    
    # Add value labels on GPU bars
    add_value_labels(ax2, bars1_gpu, stelar_gpu)
    add_value_labels(ax2, bars2_gpu, astral_gpu)
    
    ax2.legend(loc='upper left', frameon=True, fancybox=True, shadow=True)
    ax2.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    
    # Add overall title with more spacing
    fig.suptitle('Memory Usage Comparison: STELAR vs ASTRAL', 
                 fontsize=18, fontweight='bold', y=0.98)
    
    # Adjust layout to prevent overlap
    plt.tight_layout()
    plt.subplots_adjust(top=0.85)  # Make more room for the main title
    
    # Save the plot
    plt.savefig('memory_comparison.png', dpi=300, bbox_inches='tight', 
                facecolor='white', edgecolor='none')
    plt.savefig('memory_comparison.pdf', bbox_inches='tight', 
                facecolor='white', edgecolor='none')
    
    print("Memory comparison plot saved as 'memory_comparison.png' and 'memory_comparison.pdf'")
    
    return fig

def print_summary_statistics(df):
    """Print summary statistics for the comparison"""
    print("\n" + "="*60)
    print("PERFORMANCE COMPARISON SUMMARY")
    print("="*60)
    
    stelar_data = df[df['alg'] == 'stelar']
    astral_data = df[df['alg'] == 'astral']
    
    print("\nRunning Time Analysis:")
    print("-" * 30)
    for _, stelar_row in stelar_data.iterrows():
        astral_row = astral_data[astral_data['num-taxa'] == stelar_row['num-taxa']].iloc[0]
        speedup = astral_row['running-time-s'] / stelar_row['running-time-s']
        print(f"Taxa {int(stelar_row['num-taxa'])}: STELAR {stelar_row['running-time-s']:.1f}s vs ASTRAL {astral_row['running-time-s']:.1f}s (Speedup: {speedup:.1f}x)")
    
    print("\nMemory Usage Analysis (GB):")
    print("-" * 30)
    for _, stelar_row in stelar_data.iterrows():
        astral_row = astral_data[astral_data['num-taxa'] == stelar_row['num-taxa']].iloc[0]
        stelar_total = (stelar_row['max-cpu-mb'] + stelar_row['max-gpu-mb']) / 1024
        astral_total = (astral_row['max-cpu-mb'] + astral_row['max-gpu-mb']) / 1024
        reduction = astral_total / stelar_total
        print(f"Taxa {int(stelar_row['num-taxa'])}: STELAR {stelar_total:.1f}GB vs ASTRAL {astral_total:.1f}GB (Memory Reduction: {reduction:.1f}x)")

def main():
    """Main function to generate all plots"""
    print("Loading performance data...")
    df = load_data()
    
    print("Generating time comparison plot...")
    time_fig = create_time_comparison_plot(df)
    
    print("Generating memory comparison plot...")
    memory_fig = create_memory_comparison_plot(df)
    
    # Print summary statistics
    print_summary_statistics(df)
    
    print("\nAll plots generated successfully!")
    print("Files created:")
    print("- time_comparison.png (high-resolution)")
    print("- time_comparison.pdf (vector format)")
    print("- memory_comparison.png (high-resolution)")
    print("- memory_comparison.pdf (vector format)")

if __name__ == "__main__":
    main()
