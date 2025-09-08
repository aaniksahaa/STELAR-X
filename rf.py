#!/usr/bin/env python3
"""
Robinson-Foulds Distance Calculator for Phylogenetic Trees

This script calculates the Robinson-Foulds (RF) distance between phylogenetic trees
in Newick format using the dendropy library. The RF distance is normalized to a 0-1 range where:
- 0 means the trees are identical (same topology)
- 1 means the trees are maximally different

Requirements:
    pip install dendropy

Usage:
    python robinson_foulds.py tree1.txt tree2.txt
    python robinson_foulds.py --compare-all *.tre *.txt
    python robinson_foulds.py true_tree_48.txt out_w.tre
"""

import sys
import argparse
from typing import List, Dict, Tuple
from itertools import combinations
import os

try:
    import dendropy
    from dendropy.calculate import treecompare
except ImportError:
    print("Error: dendropy library is required. Please install it with:")
    print("pip install dendropy")
    sys.exit(1)


def load_tree_from_file(filename: str, taxon_namespace=None) -> dendropy.Tree:
    """
    Load a phylogenetic tree from a file using dendropy.
    
    Args:
        filename: Path to the tree file in Newick format
        taxon_namespace: Shared taxon namespace for comparing trees
        
    Returns:
        dendropy.Tree: dendropy Tree object
    """
    try:
        with open(filename, 'r') as f:
            content = f.read().strip()
            # Handle multi-line files - take the first non-empty line
            lines = [line.strip() for line in content.split('\n') if line.strip()]
            if not lines:
                raise ValueError(f"Empty file: {filename}")
            newick = lines[0]
            
            # Load tree using dendropy
            tree = dendropy.Tree.get(
                data=newick,
                schema="newick",
                taxon_namespace=taxon_namespace
            )
            return tree
            
    except FileNotFoundError:
        raise FileNotFoundError(f"Tree file not found: {filename}")
    except Exception as e:
        raise ValueError(f"Error parsing tree file {filename}: {str(e)}")


def robinson_foulds_distance(tree1: dendropy.Tree, tree2: dendropy.Tree) -> float:
    """
    Calculate the normalized Robinson-Foulds distance between two trees.
    
    Args:
        tree1: First dendropy Tree object
        tree2: Second dendropy Tree object
        
    Returns:
        float: RF distance normalized to [0, 1] where:
               0 = identical trees
               1 = maximally different trees
    """
    # Ensure trees have the same taxon namespace
    if tree1.taxon_namespace is not tree2.taxon_namespace:
        raise ValueError("Trees must have the same taxon namespace for comparison")
    
    # Calculate Robinson-Foulds distance using dendropy
    rf_distance = treecompare.symmetric_difference(tree1, tree2)
    
    # Calculate maximum possible RF distance
    n_taxa = len(tree1.taxon_namespace)
    if n_taxa <= 3:
        return 0.0  # Trees with 3 or fewer taxa have no internal branches
    
    # Maximum RF distance is 2 * (n - 3) for n taxa
    max_rf = 2 * (n_taxa - 3)
    
    # Normalize to [0, 1]
    normalized_rf = rf_distance / max_rf if max_rf > 0 else 0.0
    return min(normalized_rf, 1.0)


def get_tree_info(tree: dendropy.Tree) -> Dict:
    """Get basic information about a tree."""
    leaves = [taxon.label for taxon in tree.taxon_namespace]
    internal_nodes = [node for node in tree.internal_nodes()]
    
    return {
        'n_taxa': len(leaves),
        'n_internal_nodes': len(internal_nodes),
        'taxa': sorted(leaves)
    }


def compare_trees(file1: str, file2: str, verbose: bool = False) -> float:
    """
    Compare two tree files and return RF distance.
    
    Args:
        file1: Path to first tree file
        file2: Path to second tree file
        verbose: Whether to print detailed information
        
    Returns:
        float: Normalized RF distance
    """
    # Create shared taxon namespace
    taxon_namespace = dendropy.TaxonNamespace()
    
    # Load trees with shared namespace
    tree1 = load_tree_from_file(file1, taxon_namespace)
    tree2 = load_tree_from_file(file2, taxon_namespace)
    
    if verbose:
        info1 = get_tree_info(tree1)
        info2 = get_tree_info(tree2)
        print(f"\nTree 1 ({file1}):")
        print(f"  Taxa: {info1['n_taxa']}")
        print(f"  Internal nodes: {info1['n_internal_nodes']}")
        
        print(f"\nTree 2 ({file2}):")
        print(f"  Taxa: {info2['n_taxa']}")
        print(f"  Internal nodes: {info2['n_internal_nodes']}")
        print()
    
    return robinson_foulds_distance(tree1, tree2)


def compare_all_trees(filenames: List[str], verbose: bool = False) -> Dict[Tuple[str, str], float]:
    """
    Compare all pairs of trees and return a dictionary of RF distances.
    
    Args:
        filenames: List of tree file paths
        verbose: Whether to print detailed information
        
    Returns:
        Dict: Dictionary mapping (file1, file2) tuples to RF distances
    """
    results = {}
    trees = {}
    
    print(f"Loading {len(filenames)} tree files...")
    
    # Create shared taxon namespace
    taxon_namespace = dendropy.TaxonNamespace()
    
    # Load all trees with shared namespace
    for filename in filenames:
        try:
            trees[filename] = load_tree_from_file(filename, taxon_namespace)
            info = get_tree_info(trees[filename])
            print(f"✓ Loaded {filename}: {info['n_taxa']} taxa, {info['n_internal_nodes']} internal nodes")
        except Exception as e:
            print(f"✗ Error loading {filename}: {e}")
            continue
    
    if len(trees) < 2:
        print("Error: Need at least 2 valid trees for comparison")
        return results
    
    print(f"\nComparing {len(trees)} trees ({len(list(combinations(trees.keys(), 2)))} pairwise comparisons)...")
    print("-" * 60)
    
    # Compare all pairs
    for file1, file2 in combinations(trees.keys(), 2):
        try:
            rf_dist = robinson_foulds_distance(trees[file1], trees[file2])
            results[(file1, file2)] = rf_dist
            
            # Format file names for display
            name1 = os.path.basename(file1)
            name2 = os.path.basename(file2)
            print(f"{name1:<20} vs {name2:<20}: {rf_dist:.4f}")
            
        except Exception as e:
            print(f"✗ Error comparing {file1} and {file2}: {e}")
    
    return results


def main():
    """Main function with command-line interface."""
    parser = argparse.ArgumentParser(
        description="Calculate Robinson-Foulds distance between phylogenetic trees using dendropy",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python robinson_foulds.py tree1.txt tree2.txt
    python robinson_foulds.py --compare-all *.tre *.txt
    python robinson_foulds.py true_tree_48.txt out_w.tre
    python robinson_foulds.py --compare-all --verbose true_tree_48.txt out_w.tre out_u.tre

The Robinson-Foulds distance is normalized to [0, 1]:
    0.0 = Identical tree topologies
    1.0 = Maximally different topologies
        """
    )
    
    parser.add_argument('trees', nargs='*', help='Tree files to compare')
    parser.add_argument('--compare-all', action='store_true',
                       help='Compare all possible pairs of input trees')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Print detailed information about trees')
    
    args = parser.parse_args()
    
    if not args.trees:
        print("Error: Please provide tree files to compare")
        parser.print_help()
        return 1
    
    # Check if files exist
    missing_files = [f for f in args.trees if not os.path.exists(f)]
    if missing_files:
        print(f"Error: The following files do not exist: {', '.join(missing_files)}")
        return 1
    
    try:
        if args.compare_all:
            if len(args.trees) < 2:
                print("Error: Need at least 2 trees for comparison")
                return 1
            
            results = compare_all_trees(args.trees, verbose=args.verbose)
            
            if results:
                print("\n" + "="*60)
                print("SUMMARY OF ALL COMPARISONS:")
                print("="*60)
                for (file1, file2), rf_dist in sorted(results.items(), key=lambda x: x[1]):
                    name1 = os.path.basename(file1)
                    name2 = os.path.basename(file2)
                    similarity = (1 - rf_dist) * 100
                    print(f"{name1:<20} vs {name2:<20}: RF={rf_dist:.4f} (similarity: {similarity:.1f}%)")
            
        else:
            if len(args.trees) != 2:
                print("Error: Please provide exactly 2 tree files for comparison")
                print("Or use --compare-all flag to compare multiple trees")
                return 1
            
            rf_dist = compare_trees(args.trees[0], args.trees[1], verbose=args.verbose)
            similarity = (1 - rf_dist) * 100
            
            name1 = os.path.basename(args.trees[0])
            name2 = os.path.basename(args.trees[1])
            
            print(f"Robinson-Foulds distance: {rf_dist:.4f}")
            print(f"Tree similarity: {similarity:.1f}%")
            print(f"Comparison: {name1} vs {name2}")
    
    except Exception as e:
        print(f"Error: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main()) 