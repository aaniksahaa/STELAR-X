#!/usr/bin/env python3
"""
summary.py

Visualize a phylogenetic tree in a terminal-like tree format.
Since phylogenetic trees can be very large, this script supports a max depth parameter
to limit the visualization depth. For internal nodes beyond the max depth, a summary
of leaf taxa in the subtree is shown.

Usage:
    python summary.py tree.tre                      # Default: show trees 0-9 at depth 2
    python summary.py tree.tre -l 3                 # Show up to depth 3
    python summary.py tree.tre --index 5            # Show single tree at index 5
    python summary.py tree.tre -si 0 -ei 10         # Show trees 0-9 (first 10)
    python summary.py tree.tre -si 100 -ei 110      # Show trees 100-109
"""

import sys
sys.setrecursionlimit(10_000_000)

import argparse
import os
from concurrent.futures import ProcessPoolExecutor
from typing import List, Optional, Tuple


try:
    import dendropy
except ImportError:
    print("Error: dendropy library is required. Please install it with:")
    print("pip install dendropy")
    sys.exit(1)


# ANSI color codes for pretty output
class Colors:
    RESET = "\033[0m"
    BOLD = "\033[1m"
    GREEN = "\033[32m"
    YELLOW = "\033[33m"
    BLUE = "\033[34m"
    MAGENTA = "\033[35m"
    CYAN = "\033[36m"
    WHITE = "\033[37m"
    GRAY = "\033[90m"


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


def _worker_count_leaves(args: Tuple[int, str]) -> Tuple[int, Optional[int], str]:
    idx, tree_str = args
    from io import StringIO
    old_out, old_err = sys.stdout, sys.stderr
    buf = StringIO()
    sys.stdout = buf
    sys.stderr = buf
    try:
        tree = dendropy.Tree.get(data=tree_str, schema="newick", preserve_underscores=True)
        num_leaves = len(tree.leaf_nodes())
    except Exception as e:
        num_leaves = None
        print(f"  Tree {idx}: Error parsing - {e}")
    finally:
        sys.stdout = old_out
        sys.stderr = old_err
    return idx, num_leaves, buf.getvalue()


def load_tree_from_file(filepath: str, index: int = 0) -> Optional[dendropy.Tree]:
    """
    Load a phylogenetic tree from a file.
    
    Args:
        filepath: Path to the tree file in Newick format
        index: 0-based index of which tree to load (for files with multiple trees)
        
    Returns:
        dendropy.Tree: The loaded tree, or None if loading fails
    """
    try:
        with open(filepath, 'r') as f:
            content = f.read().strip()
        
        # Split by semicolons to handle multiple trees
        tree_strings = [t.strip() for t in content.split(';') if t.strip()]
        
        if not tree_strings:
            print(f"Error: No trees found in file '{filepath}'", file=sys.stderr)
            return None
        
        if index >= len(tree_strings):
            print(f"Error: Tree index {index} out of range. File contains {len(tree_strings)} tree(s).", file=sys.stderr)
            return None
        
        tree_str = tree_strings[index] + ';'
        tree = dendropy.Tree.get(data=tree_str, schema="newick", preserve_underscores=True)
        
        return tree
        
    except FileNotFoundError:
        print(f"Error: File not found: {filepath}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"Error loading tree: {e}", file=sys.stderr)
        return None


def get_all_leaf_labels(node: dendropy.Node) -> List[str]:
    """
    Get all leaf labels in the subtree rooted at the given node.
    
    Args:
        node: The root node of the subtree
        
    Returns:
        List of leaf taxon labels
    """
    leaves = []
    for leaf in node.leaf_iter():
        if leaf.taxon and leaf.taxon.label:
            leaves.append(leaf.taxon.label)
    return leaves


def format_leaf_summary(leaves: List[str], max_show: int = 3) -> str:
    """
    Format a summary of leaves for display.
    
    Args:
        leaves: List of leaf labels
        max_show: Maximum number of leaves to show before truncating
        
    Returns:
        Formatted string like "[A,B,C]" or "[A,B,C,...]"
    """
    n = len(leaves)
    if n == 0:
        return "[]"
    
    if n <= max_show:
        return f"[{','.join(leaves)}]"
    else:
        shown = leaves[:max_show]
        return f"[{','.join(shown)},...]"


def get_node_label(node: dendropy.Node) -> str:
    """
    Get a display label for a node.
    
    Args:
        node: The node to get label for
        
    Returns:
        Label string for the node
    """
    if node.is_leaf():
        if node.taxon and node.taxon.label:
            return node.taxon.label
        return "(unnamed leaf)"
    else:
        if node.label:
            return f"({node.label})"
        return ""


def print_tree_recursive(
    node: dendropy.Node,
    current_depth: int,
    max_depth: int,
    prefix: str,
    is_last: bool,
    is_root: bool = False
) -> List[str]:
    """
    Recursively build tree visualization lines.
    
    Args:
        node: Current node to visualize
        current_depth: Current depth in the tree (0 = root)
        max_depth: Maximum depth to visualize
        prefix: Prefix string for indentation
        is_last: Whether this is the last child of its parent
        is_root: Whether this is the root node
        
    Returns:
        List of formatted lines for the tree visualization
    """
    lines = []
    
    # Determine the branch character
    if is_root:
        branch = ""
        new_prefix = ""
    else:
        if is_last:
            branch = "└── "
            new_prefix = prefix + "    "
        else:
            branch = "├── "
            new_prefix = prefix + "│   "
    
    # Build the node representation
    if node.is_leaf():
        # Leaf node - show the taxon label
        label = get_node_label(node)
        lines.append(f"{prefix}{branch}{Colors.GREEN}{label}{Colors.RESET}")
    else:
        # Internal node
        children = node.child_nodes()
        num_leaves = len(list(node.leaf_iter()))
        node_label = get_node_label(node)
        
        if current_depth >= max_depth:
            # Beyond max depth - show summary
            leaves = get_all_leaf_labels(node)
            summary = format_leaf_summary(leaves, max_show=3)
            count_info = f"{Colors.YELLOW}({len(leaves)} leaves){Colors.RESET}"
            
            if node_label:
                lines.append(f"{prefix}{branch}{Colors.CYAN}{node_label}{Colors.RESET} {Colors.MAGENTA}{summary}{Colors.RESET} {count_info}")
            else:
                lines.append(f"{prefix}{branch}{Colors.MAGENTA}{summary}{Colors.RESET} {count_info}")
        else:
            # Within max depth - show node and recurse
            num_children = len(children)
            if node_label:
                node_info = f"{Colors.CYAN}{node_label}{Colors.RESET} "
            else:
                node_info = ""
            
            child_info = f"{Colors.GRAY}[{num_children} children, {num_leaves} leaves in subtree]{Colors.RESET}"
            
            if is_root:
                lines.append(f"{Colors.BOLD}{Colors.BLUE}ROOT{Colors.RESET} {node_info}{child_info}")
            else:
                lines.append(f"{prefix}{branch}{node_info}{child_info}")
            
            # Recurse into children
            for i, child in enumerate(children):
                is_last_child = (i == len(children) - 1)
                child_lines = print_tree_recursive(
                    child,
                    current_depth + 1,
                    max_depth,
                    new_prefix,
                    is_last_child,
                    is_root=False
                )
                lines.extend(child_lines)
    
    return lines


def visualize_tree(tree: dendropy.Tree, max_depth: int = 2) -> None:
    """
    Visualize a tree in terminal tree format.
    
    Args:
        tree: The tree to visualize
        max_depth: Maximum depth to show (default 2)
    """
    if tree is None:
        print("Error: No tree to visualize", file=sys.stderr)
        return
    
    root = tree.seed_node
    
    # Print tree statistics header
    num_leaves = len(tree.leaf_nodes())
    num_internal = len(tree.internal_nodes())
    
    # Calculate tree depth
    def get_tree_depth(node):
        if not node.child_nodes():
            return 0
        return 1 + max(get_tree_depth(child) for child in node.child_nodes())
    
    tree_depth = get_tree_depth(root)
    
    print()
    print(f"{Colors.BOLD}{'═' * 60}{Colors.RESET}")
    print(f"{Colors.BOLD}Tree Summary{Colors.RESET}")
    print(f"{'═' * 60}")
    print(f"  {Colors.CYAN}Total leaves:{Colors.RESET}         {num_leaves}")
    print(f"  {Colors.CYAN}Total internal nodes:{Colors.RESET} {num_internal}")
    print(f"  {Colors.CYAN}Tree depth:{Colors.RESET}           {tree_depth}")
    print(f"  {Colors.CYAN}Showing up to depth:{Colors.RESET}  {max_depth}")
    print(f"{'═' * 60}")
    print()
    
    # Generate and print tree visualization
    lines = print_tree_recursive(root, 0, max_depth, "", True, is_root=True)
    for line in lines:
        print(line)
    
    print()
    print(f"{Colors.GRAY}Legend: {Colors.GREEN}■ Leaf{Colors.RESET} {Colors.GRAY}| {Colors.MAGENTA}[...] Collapsed subtree{Colors.RESET} {Colors.GRAY}| {Colors.CYAN}(label) Internal node label{Colors.RESET}")
    print()


def count_trees_in_file(filepath: str) -> int:
    """Count the number of trees in a file."""
    try:
        with open(filepath, 'r') as f:
            content = f.read().strip()
        tree_strings = [t.strip() for t in content.split(';') if t.strip()]
        return len(tree_strings)
    except:
        return 0


def main():
    parser = argparse.ArgumentParser(
        description='Visualize phylogenetic tree(s) in a terminal tree format',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python summary.py tree.tre                      # Default depth of 2
    python summary.py tree.tre -l 3                 # Show up to depth 3
    python summary.py tree.tre --level 5            # Show up to depth 5
    python summary.py tree.tre --index 0            # Show first tree (0-indexed)
    python summary.py trees.tre --index 5 -l 4      # Show 6th tree up to depth 4

The tree is displayed in a vertical format similar to the 'tree' command.
For very large trees, use --level to limit the display depth.
Subtrees beyond the max depth are summarized as [leaf1,leaf2,...].
        """
    )
    
    parser.add_argument('tree_file', help='Path to the tree file (.tre)')
    parser.add_argument('-l', '--level', type=int, default=2,
                        help='Maximum depth to display (default: 2)')
    parser.add_argument('-i', '--index', type=int, default=None,
                        help='0-based index of single tree to visualize (overrides -si/-ei)')
    parser.add_argument('-si', '--start-idx', type=int, default=0,
                        help='Start index of tree range to visualize (default: 0)')
    parser.add_argument('-ei', '--end-idx', type=int, default=10,
                        help='End index (exclusive) of tree range to visualize (default: 10)')
    parser.add_argument('--no-color', action='store_true',
                        help='Disable colored output')
    parser.add_argument('--list', action='store_true',
                        help='List all trees in the file with basic info instead of visualizing')
    parser.add_argument('--num-workers', type=int, default=None,
                        help='Number of parallel worker processes for --list (default: cpu count)')
    
    args = parser.parse_args()
    
    # Disable colors if requested
    if args.no_color:
        Colors.RESET = ""
        Colors.BOLD = ""
        Colors.GREEN = ""
        Colors.YELLOW = ""
        Colors.BLUE = ""
        Colors.MAGENTA = ""
        Colors.CYAN = ""
        Colors.WHITE = ""
        Colors.GRAY = ""
    
    # List mode - show all trees in file with basic info
    if args.list:
        num_trees = count_trees_in_file(args.tree_file)
        print(f"\n{Colors.BOLD}Trees in file: {args.tree_file}{Colors.RESET}")
        print(f"{'=' * 60}")
        print(f"Total trees: {num_trees}")
        print()
        
        if num_trees > 0:
            num_workers = args.num_workers if args.num_workers is not None else (os.cpu_count() or 1)
            if num_workers < 1:
                num_workers = 1
            if num_workers == 1:
                for i, ts in iter_tree_strings_from_file_robust(args.tree_file):
                    try:
                        tree = dendropy.Tree.get(data=ts, schema="newick", preserve_underscores=True)
                        num_leaves = len(tree.leaf_nodes())
                        print(f"  Tree {i - 1}: {num_leaves} leaves")
                    except Exception as e:
                        print(f"  Tree {i - 1}: Error parsing - {e}")
            else:
                work_iter = ((idx, ts) for idx, ts in iter_tree_strings_from_file_robust(args.tree_file))
                with ProcessPoolExecutor(max_workers=num_workers) as ex:
                    for idx, num_leaves, log in ex.map(_worker_count_leaves, work_iter, chunksize=1):
                        if log:
                            if not log.endswith("\n"):
                                log = log + "\n"
                            sys.stdout.write(log)
                            sys.stdout.flush()
                        if num_leaves is None:
                            continue
                        print(f"  Tree {idx - 1}: {num_leaves} leaves")
        return
    
    # Load and visualize tree(s)
    num_trees = count_trees_in_file(args.tree_file)
    
    if num_trees == 0:
        print(f"Error: No trees found in file", file=sys.stderr)
        sys.exit(1)
    
    # Determine which trees to show
    if args.index is not None:
        # Single tree mode
        start_idx = args.index
        end_idx = args.index + 1
    else:
        # Range mode
        start_idx = args.start_idx
        end_idx = min(args.end_idx, num_trees)  # Don't exceed available trees
    
    if start_idx >= num_trees:
        print(f"Error: Start index {start_idx} out of range. File contains {num_trees} tree(s).", file=sys.stderr)
        sys.exit(1)
    
    # Process trees one by one (memory efficient)
    with open(args.tree_file, 'r') as f:
        content = f.read()
    
    tree_strings = [t.strip() for t in content.split(';') if t.strip()]
    
    trees_shown = 0
    for idx in range(start_idx, end_idx):
        if idx >= len(tree_strings):
            break
        
        tree_str = tree_strings[idx] + ';'
        try:
            tree = dendropy.Tree.get(data=tree_str, schema="newick", preserve_underscores=True)
        except Exception as e:
            print(f"Error parsing tree {idx}: {e}", file=sys.stderr)
            continue
        
        # Show which tree we're visualizing
        if end_idx - start_idx > 1 or num_trees > 1:
            print(f"\n{Colors.GRAY}(Showing tree {idx} of {num_trees} trees in file){Colors.RESET}")
        
        visualize_tree(tree, max_depth=args.level)
        trees_shown += 1
        
        # Free memory
        del tree
    
    if trees_shown == 0:
        print("No trees were displayed.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
