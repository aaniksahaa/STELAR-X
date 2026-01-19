import dendropy
import argparse
import sys
import matplotlib.pyplot as plt
import random
from collections import Counter

# increase recursion limit (be conservative but high enough)
sys.setrecursionlimit(10_000_000)

def inspect_tree(tree, tree_index=None):
    """Inspect a single tree and print its statistics."""
    freq_map = {}
    
    # Count nodes
    num_leaves = len(tree.leaf_nodes())
    num_internal_nodes = len(tree.internal_nodes())
    total_nodes = len(tree.nodes())
    
    # iterate over all internal nodes
    for node in tree.internal_nodes():
        num_childs = len(node.child_nodes())
        freq_map[num_childs] = freq_map.get(num_childs, 0) + 1
    
    # compute depth of the tree (max root-to-leaf path length, in edges)
    def get_depth(node):
        if not node.child_nodes():
            return 0
        return 1 + max(get_depth(child) for child in node.child_nodes())
    
    depth = get_depth(tree.seed_node)
    
    # Get furcation size at topmost node (root)
    root_furcation = len(tree.seed_node.child_nodes())
    
    # Print results
    if tree_index is not None:
        print(f"Tree {tree_index + 1}:")
    print(f"Number of leaves: {num_leaves}")
    print(f"Number of internal nodes: {num_internal_nodes}")
    print(f"Total nodes: {total_nodes}")
    print(f"Child frequency map: {freq_map}")
    print(f"Tree depth (in edges): {depth}")
    print(f"Topmost node furcation size: {root_furcation}")
    if tree_index is not None:
        print()
    
    return {
        'num_leaves': num_leaves,
        'num_internal_nodes': num_internal_nodes,
        'total_nodes': total_nodes,
        'freq_map': freq_map,
        'depth': depth,
        'root_furcation': root_furcation
    }

def plot_taxa_distribution(taxa_counts, filename, ax=None):
    """Plot histogram of taxa counts across trees."""
    if ax is None:
        plt.figure(figsize=(10, 6))
        ax = plt.gca()
    
    # Create histogram
    counts = Counter(taxa_counts)
    taxa_values = sorted(counts.keys())
    frequencies = [counts[taxa] for taxa in taxa_values]
    
    ax.bar(taxa_values, frequencies, alpha=0.7, edgecolor='black')
    ax.set_xlabel('Number of Taxa (Leaves)')
    ax.set_ylabel('Number of Trees')
    ax.set_title(f'Distribution of Taxa Counts Across {len(taxa_counts)} Trees')
    ax.grid(True, alpha=0.3)
    
    # Add text annotations on bars
    for taxa, freq in zip(taxa_values, frequencies):
        ax.text(taxa, freq + 0.1, str(freq), ha='center', va='bottom')
    
    if ax is None:
        plt.tight_layout()
        plot_filename = f"{filename}_taxa_distribution.png"
        plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
        print(f"Taxa distribution plot saved as: {plot_filename}")
        plt.show()

def plot_children_frequency(freq_map, tree_index, filename, ax=None):
    """Plot child frequency map for a specific tree."""
    if ax is None:
        plt.figure(figsize=(8, 6))
        ax = plt.gca()
    
    children_counts = sorted(freq_map.keys())
    frequencies = [freq_map[count] for count in children_counts]
    
    ax.bar(children_counts, frequencies, alpha=0.7, color='orange', edgecolor='black')
    ax.set_xlabel('Number of Children per Internal Node')
    ax.set_ylabel('Frequency (Number of Internal Nodes)')
    ax.set_title(f'Child Frequency Distribution - Tree {tree_index + 1}')
    ax.grid(True, alpha=0.3)
    
    # Add text annotations on bars
    for children, freq in zip(children_counts, frequencies):
        ax.text(children, freq + 0.1, str(freq), ha='center', va='bottom')
    
    ax.set_xticks(children_counts)
    
    if ax is None:
        plt.tight_layout()
        plot_filename = f"{filename}_tree{tree_index + 1}_children_freq.png"
        plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
        print(f"Children frequency plot saved as: {plot_filename}")
        plt.show()

def plot_combined(taxa_counts, freq_map, tree_index, filename):
    """Plot both taxa distribution and children frequency side by side."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Plot taxa distribution on left
    plot_taxa_distribution(taxa_counts, filename, ax1)
    
    # Plot children frequency on right
    plot_children_frequency(freq_map, tree_index, filename, ax2)
    
    plt.tight_layout()
    plot_filename = f"{filename}_combined_analysis.png"
    plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
    print(f"Combined analysis plot saved as: {plot_filename}")
    plt.show()

def main():
    parser = argparse.ArgumentParser(description='Inspect phylogenetic tree(s)')
    parser.add_argument('tree_file', help='Path to the tree file (.tre)')
    parser.add_argument('--all', action='store_true', 
                       help='Process all trees in the file (split by semicolons)')
    parser.add_argument('--plot-taxa', action='store_true',
                       help='Show histogram of taxa counts across all trees (requires --all)')
    parser.add_argument('--plot-children', action='store_true',
                       help='Show child frequency plot for a random tree (requires --all)')
    
    args = parser.parse_args()

    args.plot_taxa = True
    args.plot_children = True
    
    # # Validate plotting arguments
    # if (args.plot_taxa or args.plot_children) and not args.all:
    #     print("Error: Plotting options require --all flag")
    #     return
    
    # Read the tree file
    with open(args.tree_file, 'r') as f:
        content = f.read().strip()
    
    # Split by semicolons early
    tree_strings = [t.strip() for t in content.split(';') if t.strip()]
    
    # Get base filename for plot saving
    base_filename = args.tree_file.rsplit('.', 1)[0]
    
    if args.all:
        # Process all trees
        print(f"Found {len(tree_strings)} trees in the file.\n")

        mn = float('inf')
        mx = 0
        taxa_counts = []
        tree_stats = []
        
        for i, tree_string in enumerate(tree_strings):
            try:
                tree = dendropy.Tree.get(data=tree_string + ';', schema="newick")
                stats = inspect_tree(tree, i)
                tree_stats.append((i, stats))
                
                num_leaves = len(tree.leaf_nodes())
                taxa_counts.append(num_leaves)
                mx = max(mx, num_leaves)
                mn = min(mn, num_leaves)
            except Exception as e:
                print(f"Error processing tree {i + 1}: {e}")
        
        print(f"Max number of leaves: {mx}")
        print(f"Min number of leaves: {mn}")
        
        # Generate plots if requested
        if args.plot_taxa and args.plot_children and taxa_counts and tree_stats:
            # Both plots requested - show side by side
            random_tree_idx, random_stats = random.choice(tree_stats)
            print(f"\nGenerating combined analysis plot (taxa distribution + children frequency for tree {random_tree_idx + 1})...")
            plot_combined(taxa_counts, random_stats['freq_map'], random_tree_idx, base_filename)
        elif args.plot_taxa and taxa_counts:
            # Only taxa distribution
            print(f"\nGenerating taxa distribution plot...")
            plot_taxa_distribution(taxa_counts, base_filename)
        elif args.plot_children and tree_stats:
            # Only children frequency
            random_tree_idx, random_stats = random.choice(tree_stats)
            print(f"\nGenerating children frequency plot for randomly selected tree {random_tree_idx + 1}...")
            plot_children_frequency(random_stats['freq_map'], random_tree_idx, base_filename)
            
    else:
        # Process only the first tree
        if tree_strings:
            try:
                tree = dendropy.Tree.get(data=tree_strings[0] + ';', schema="newick")
                inspect_tree(tree)
            except Exception as e:
                print(f"Error processing tree: {e}")
        else:
            print("No trees found in the file.")

if __name__ == "__main__":
    main()
