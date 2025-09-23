#!/usr/bin/env python3
"""
Dataset Analysis Tool for Phylogenetic Trees

This script analyzes gene trees against a species tree by calculating:
1. Average pairwise Robinson-Foulds distance between gene trees
2. Average Robinson-Foulds distance between gene trees and species tree

Usage:
    python analyze-dataset.py gene_trees.tre species_tree.tre
    python analyze-dataset.py all_gt.tre true_st.tre --max_gt 15
    python analyze-dataset.py gene_trees.tre species_tree.tre --max_gt 5 --verbose

Requirements:
    pip install dendropy
"""

import sys
import argparse
import os
from typing import List, Tuple
from itertools import combinations

# Import functions from rf.py in the same directory
try:
    from rf import load_tree_from_file, robinson_foulds_distance, get_tree_info
    import dendropy
except ImportError as e:
    print(f"Error importing required modules: {e}")
    print("Make sure rf.py is in the same directory and dendropy is installed:")
    print("pip install dendropy")
    sys.exit(1)


def load_gene_trees(filename: str, max_trees: int = 5, taxon_namespace=None) -> List[dendropy.Tree]:
    """
    Load gene trees from a file, taking at most max_trees from the beginning.
    
    Args:
        filename: Path to the gene trees file (Newick format, one tree per line)
        max_trees: Maximum number of trees to load (default: 5)
        taxon_namespace: Shared taxon namespace for tree comparison
        
    Returns:
        List[dendropy.Tree]: List of loaded gene trees
    """
    trees = []
    
    try:
        with open(filename, 'r') as f:
            lines = [line.strip() for line in f.readlines() if line.strip()]
            
        if not lines:
            raise ValueError(f"No trees found in file: {filename}")
        
        print(f"Found {len(lines)} trees in {filename}")
        trees_to_load = min(max_trees, len(lines))
        print(f"Loading first {trees_to_load} trees...")
        
        for i, newick in enumerate(lines[:trees_to_load]):
            try:
                tree = dendropy.Tree.get(
                    data=newick,
                    schema="newick",
                    taxon_namespace=taxon_namespace
                )
                trees.append(tree)
                print(f"✓ Loaded gene tree {i+1}")
            except Exception as e:
                print(f"✗ Error loading gene tree {i+1}: {e}")
                continue
        
        if not trees:
            raise ValueError("No valid trees could be loaded")
            
        print(f"Successfully loaded {len(trees)} gene trees")
        return trees
        
    except FileNotFoundError:
        raise FileNotFoundError(f"Gene trees file not found: {filename}")
    except Exception as e:
        raise ValueError(f"Error loading gene trees from {filename}: {str(e)}")


def calculate_average_pairwise_rf(trees: List[dendropy.Tree]) -> Tuple[float, int]:
    """
    Calculate the average pairwise Robinson-Foulds distance between trees.
    
    Args:
        trees: List of dendropy Tree objects
        
    Returns:
        Tuple[float, int]: (average_rf_distance, number_of_comparisons)
    """
    if len(trees) < 2:
        return 0.0, 0
    
    distances = []
    comparisons = list(combinations(range(len(trees)), 2))
    
    print(f"\nCalculating pairwise RF distances ({len(comparisons)} comparisons)...")
    
    for i, (idx1, idx2) in enumerate(comparisons):
        try:
            rf_dist = robinson_foulds_distance(trees[idx1], trees[idx2])
            distances.append(rf_dist)
            
            if i % 10 == 0 or i == len(comparisons) - 1:  # Progress indicator
                print(f"  Progress: {i+1}/{len(comparisons)} comparisons")
                
        except Exception as e:
            print(f"✗ Error comparing trees {idx1+1} and {idx2+1}: {e}")
            continue
    
    if not distances:
        return 0.0, 0
    
    avg_distance = sum(distances) / len(distances)
    return avg_distance, len(distances)


def calculate_average_species_rf(gene_trees: List[dendropy.Tree], species_tree: dendropy.Tree) -> Tuple[float, int]:
    """
    Calculate the average Robinson-Foulds distance between gene trees and species tree.
    
    Args:
        gene_trees: List of gene tree objects
        species_tree: Species tree object
        
    Returns:
        Tuple[float, int]: (average_rf_distance, number_of_comparisons)
    """
    distances = []
    
    print(f"\nCalculating RF distances between gene trees and species tree...")
    
    for i, gene_tree in enumerate(gene_trees):
        try:
            rf_dist = robinson_foulds_distance(gene_tree, species_tree)
            distances.append(rf_dist)
            print(f"  Gene tree {i+1} vs species tree: {rf_dist:.4f}")
            
        except Exception as e:
            print(f"✗ Error comparing gene tree {i+1} with species tree: {e}")
            continue
    
    if not distances:
        return 0.0, 0
    
    avg_distance = sum(distances) / len(distances)
    return avg_distance, len(distances)


def main():
    """Main function with command-line interface."""
    parser = argparse.ArgumentParser(
        description="Analyze phylogenetic dataset by comparing gene trees with each other and with species tree",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python analyze-dataset.py all_gt.tre true_st.tre
    python analyze-dataset.py gene_trees.tre species_tree.tre --max_gt 15
    python analyze-dataset.py all_gt.tre true_st.tre --max_gt 5 --verbose

Output:
    - Average pairwise RF distance between gene trees
    - Average RF distance between gene trees and species tree
    - Number of trees analyzed and comparisons made
        """
    )
    
    parser.add_argument('gene_trees_file', help='File containing gene trees (Newick format, one per line)')
    parser.add_argument('species_tree_file', help='File containing species tree (Newick format)')
    parser.add_argument('--max_gt', type=int, default=5,
                       help='Maximum number of gene trees to analyze (default: 5)')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Print detailed information about trees')
    
    args = parser.parse_args()
    
    # Validate input files
    if not os.path.exists(args.gene_trees_file):
        print(f"Error: Gene trees file not found: {args.gene_trees_file}")
        return 1
    
    if not os.path.exists(args.species_tree_file):
        print(f"Error: Species tree file not found: {args.species_tree_file}")
        return 1
    
    if args.max_gt < 1:
        print("Error: --max_gt must be at least 1")
        return 1
    
    try:
        print("="*60)
        print("PHYLOGENETIC DATASET ANALYSIS")
        print("="*60)
        print(f"Gene trees file: {args.gene_trees_file}")
        print(f"Species tree file: {args.species_tree_file}")
        print(f"Maximum gene trees to analyze: {args.max_gt}")
        print()
        
        # Create shared taxon namespace for all trees
        taxon_namespace = dendropy.TaxonNamespace()
        
        # Load species tree
        print("Loading species tree...")
        species_tree = load_tree_from_file(args.species_tree_file, taxon_namespace)
        species_info = get_tree_info(species_tree)
        print(f"✓ Species tree loaded: {species_info['n_taxa']} taxa, {species_info['n_internal_nodes']} internal nodes")
        
        if args.verbose:
            print(f"  Taxa: {', '.join(sorted(species_info['taxa']))}")
        
        # Load gene trees
        print(f"\nLoading gene trees (max {args.max_gt})...")
        gene_trees = load_gene_trees(args.gene_trees_file, args.max_gt, taxon_namespace)
        
        if args.verbose:
            for i, tree in enumerate(gene_trees):
                info = get_tree_info(tree)
                print(f"  Gene tree {i+1}: {info['n_taxa']} taxa, {info['n_internal_nodes']} internal nodes")
        
        # Calculate average pairwise RF distance between gene trees
        avg_pairwise_rf, pairwise_comparisons = calculate_average_pairwise_rf(gene_trees)
        
        # Calculate average RF distance between gene trees and species tree
        avg_species_rf, species_comparisons = calculate_average_species_rf(gene_trees, species_tree)
        
        # Display results
        print("\n" + "="*60)
        print("RESULTS")
        print("="*60)
        print(f"Number of gene trees analyzed: {len(gene_trees)}")
        print(f"Species tree taxa: {species_info['n_taxa']}")
        print()
        
        print(f"Average pairwise RF distance (gene trees): {avg_pairwise_rf:.4f}")
        print(f"Average pairwise discordance (gene trees): {avg_pairwise_rf*100:.4f}")
        print(f"  Number of pairwise comparisons: {pairwise_comparisons}")
        if pairwise_comparisons > 0:
            avg_similarity = (1 - avg_pairwise_rf) * 100
            print(f"  Average pairwise similarity: {avg_similarity:.1f}%")
        print()
        
        print(f"Average RF distance (gene trees vs species tree): {avg_species_rf:.4f}")
        print(f"Average discordance (gene trees vs species tree): {avg_species_rf*100:.4f}")
        print(f"  Number of gene tree comparisons: {species_comparisons}")
        if species_comparisons > 0:
            avg_similarity = (1 - avg_species_rf) * 100
            print(f"  Average similarity to species tree: {avg_similarity:.1f}%")
        
        print("\n" + "="*60)
        
        # Summary interpretation
        if avg_pairwise_rf < 0.3:
            pairwise_desc = "very similar"
        elif avg_pairwise_rf < 0.6:
            pairwise_desc = "moderately similar"
        else:
            pairwise_desc = "quite different"
            
        if avg_species_rf < 0.3:
            species_desc = "very similar"
        elif avg_species_rf < 0.6:
            species_desc = "moderately similar"
        else:
            species_desc = "quite different"
        
        print(f"Gene trees are {pairwise_desc} to each other on average")
        print(f"Gene trees are {species_desc} to the species tree on average")
        
    except Exception as e:
        print(f"Error: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
