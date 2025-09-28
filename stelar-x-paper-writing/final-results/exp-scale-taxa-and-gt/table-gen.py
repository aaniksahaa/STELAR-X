import pandas as pd
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

def generate_latex_table(csv_file):
    # Read the CSV file
    df = pd.read_csv(csv_file)
    
    # Detect what type of scaling this is
    scaling_type = detect_scaling_type(df)
    
    # Extract relevant columns: num-taxa, gene-trees, running-time-s, max-cpu-mb, max-gpu-mb
    relevant_data = df[['num-taxa', 'gene-trees', 'running-time-s', 'max-cpu-mb', 'max-gpu-mb']].copy()
    
    # Generate appropriate caption based on scaling type
    if scaling_type == "taxa":
        caption = "Running time and memory usage of STELAR with varying number of taxa"
    else:
        caption = "Running time and memory usage of STELAR with varying number of gene trees"
    
    # Always use the same column order: No. of Taxa first, then No. of Gene Trees
    # This ensures consistency regardless of which variable is being scaled
    latex_code = f"""\\begin{{table}}[h]
\\caption{{{caption}}}
\\centering
\\begin{{tabular}}{{l l l l l}}
\\toprule
No. of Taxa & No. of Gene Trees & Running Time (s) & CPU RAM (MB) & GPU VRAM (MB) \\\\
\\midrule
"""
    
    # Add data rows with consistent ordering: always taxa first, gene trees second
    for _, row in relevant_data.iterrows():
        num_taxa = int(row['num-taxa'])
        gene_trees = int(row['gene-trees'])
        running_time = row['running-time-s']
        cpu_ram = row['max-cpu-mb']
        gpu_vram = row['max-gpu-mb']
        
        # Always show taxa first, then gene trees, regardless of which is being varied
        latex_code += f"{num_taxa}  & {gene_trees} & {running_time:.3f}   & {cpu_ram:.3f} & {gpu_vram:.3f} \\\\\n"
    
    # Close the table
    latex_code += """\\bottomrule
\\end{tabular}
\\end{table}
"""
    
    return latex_code, scaling_type

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
    
    # Generate the LaTeX table
    latex_table, scaling_type = generate_latex_table(csv_file)
    
    # Print to console
    print(latex_table)
    
    # Generate appropriate output filename
    exp_type = csv_file.replace('result-scale-', '').replace('.csv', '')
    output_file = f"scale-{exp_type}-table.tex"
    
    # Save to file
    with open(output_file, "w") as f:
        f.write(latex_table)
    
    print(f"\nLaTeX table code saved to '{output_file}'")
    print(f"Detected scaling type: {scaling_type}")
    
    # Show summary
    df = pd.read_csv(csv_file)
    if scaling_type == "taxa":
        variable_range = f"{df['num-taxa'].min():,} - {df['num-taxa'].max():,}"
        print(f"Taxa range: {variable_range}")
    else:
        variable_range = f"{df['gene-trees'].min():,} - {df['gene-trees'].max():,}"
        print(f"Gene trees range: {variable_range}")
