#!/usr/bin/env python3
"""
Script to clean Newick tree files by optionally removing branch support values.

Usage:
    python clean_tree.py input.tre output.tre                 # Remove branch support values
    python clean_tree.py input.tre output.tre --branch-support  # Keep branch support values
"""

import argparse
import re
import sys


def remove_branch_support(newick_str):
    """
    Remove branch support values from a Newick tree string.
    Branch support values appear as numbers immediately after closing parentheses.
    Example: ((A,B)95,C)100; -> ((A,B),C);
    """
    # Pattern matches numbers that appear after ) but before , or ) or ;
    # This captures branch support values
    cleaned = re.sub(r'\)(\d+)', ')', newick_str)
    return cleaned


def process_tree_file(input_path, output_path, keep_branch_support=False):
    """
    Process a Newick tree file.
    
    Args:
        input_path: Path to input tree file
        output_path: Path to output tree file
        keep_branch_support: If True, keep branch support values; if False, remove them
    """
    with open(input_path, 'r') as f:
        content = f.read()
    
    if keep_branch_support:
        # Just copy the file as-is
        output_content = content
    else:
        # Remove branch support values
        output_content = remove_branch_support(content)
    
    with open(output_path, 'w') as f:
        f.write(output_content)
    
    print(f"Processed: {input_path} -> {output_path}")
    if keep_branch_support:
        print("Branch support values: KEPT")
    else:
        print("Branch support values: REMOVED")


def main():
    parser = argparse.ArgumentParser(
        description='Clean Newick tree files by optionally removing branch support values.'
    )
    parser.add_argument('input', help='Input tree file path')
    parser.add_argument('output', help='Output tree file path')
    parser.add_argument(
        '--branch-support',
        action='store_true',
        dest='branch_support',
        help='Keep branch support values (default: remove them)'
    )
    
    args = parser.parse_args()
    
    process_tree_file(args.input, args.output, args.branch_support)


if __name__ == '__main__':
    main()
