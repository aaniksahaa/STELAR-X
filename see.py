import dendropy
import argparse
import os
import sys
import matplotlib.pyplot as plt
import random
from collections import Counter
from concurrent.futures import ProcessPoolExecutor

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


def iter_tree_strings_from_file_robust(path: str):
    tree_index = 0
    buffer = ""
    with open(path, "r", encoding="utf-8") as fh:
        while True:
            chunk = fh.read(8192)
            if not chunk:
                break
            buffer += chunk
            while ';' in buffer:
                semicolon_pos = buffer.index(';')
                tree_str = buffer[:semicolon_pos].strip()
                buffer = buffer[semicolon_pos + 1:]
                if not tree_str:
                    continue
                tree_index += 1
                yield (tree_index, tree_str + ';')
    remaining = buffer.strip()
    if remaining:
        tree_index += 1
        yield (tree_index, remaining + (';' if not remaining.endswith(';') else ''))


def _worker_inspect(args):
    tree_str, idx = args
    from io import StringIO
    old_out, old_err = sys.stdout, sys.stderr
    buf = StringIO()
    sys.stdout = buf
    sys.stderr = buf
    try:
        tree = dendropy.Tree.get(data=tree_str, schema="newick")
        stats = inspect_tree(tree, idx - 1)
        num_leaves = stats.get("num_leaves", len(tree.leaf_nodes()))
    except Exception as e:
        stats = None
        num_leaves = None
        print(f"Error processing tree {idx}: {e}")
    finally:
        sys.stdout = old_out
        sys.stderr = old_err
    return idx, stats, num_leaves, buf.getvalue()

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
    parser.add_argument('--num-workers', type=int, default=None,
                        help='Number of parallel worker processes when using --all (default: cpu count)')
    
    args = parser.parse_args()

    args.plot_taxa = True
    args.plot_children = True
    
    # # Validate plotting arguments
    # if (args.plot_taxa or args.plot_children) and not args.all:
    #     print("Error: Plotting options require --all flag")
    #     return
    
    # Get base filename for plot saving
    base_filename = args.tree_file.rsplit('.', 1)[0]
    
    if args.all:
        num_workers = args.num_workers if args.num_workers is not None else (os.cpu_count() or 1)
        if num_workers < 1:
            num_workers = 1
        # Process all trees
        # Count trees and stream them
        tree_iter = iter_tree_strings_from_file_robust(args.tree_file)
        tree_list = list(tree_iter)
        print(f"Found {len(tree_list)} trees in the file.\n")

        mn = float('inf')
        mx = 0
        taxa_counts = []
        tree_stats = []
        if num_workers == 1:
            for i, (idx, tree_string) in enumerate(tree_list):
                try:
                    tree = dendropy.Tree.get(data=tree_string, schema="newick")
                    stats = inspect_tree(tree, i)
                    tree_stats.append((i, stats))
                    num_leaves = len(tree.leaf_nodes())
                    taxa_counts.append(num_leaves)
                    mx = max(mx, num_leaves)
                    mn = min(mn, num_leaves)
                except Exception as e:
                    print(f"Error processing tree {i + 1}: {e}")
        else:
            work_iter = ((tree_str, idx) for idx, tree_str in tree_list)
            with ProcessPoolExecutor(max_workers=num_workers) as ex:
                for idx, stats, num_leaves, log in ex.map(_worker_inspect, work_iter, chunksize=1):
                    if log:
                        if not log.endswith("\n"):
                            log = log + "\n"
                        sys.stdout.write(log)
                        sys.stdout.flush()
                    if stats is None or num_leaves is None:
                        continue
                    tree_stats.append((idx - 1, stats))
                    taxa_counts.append(num_leaves)
                    mx = max(mx, num_leaves)
                    mn = min(mn, num_leaves)
        
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
        tree_list = list(iter_tree_strings_from_file_robust(args.tree_file))
        if tree_list:
            try:
                tree = dendropy.Tree.get(data=tree_list[0][1], schema="newick")
                inspect_tree(tree)
            except Exception as e:
                print(f"Error processing tree: {e}")
        else:
            print("No trees found in the file.")

if __name__ == "__main__":
    main()
