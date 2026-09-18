#!/usr/bin/env python3
"""
Generate incomplete gene trees by randomly pruning taxa from complete Newick trees.

Usage:
  python3 test/gen_incomplete.py input.tre output.tre [options]

Options:
  --fraction F   Remove F fraction of taxa from each tree (default 0.3)
  --seed N       Random seed for reproducibility (default 42)
  --min-keep N   Minimum taxa to keep per tree (default 4)
  --stats        Print per-tree statistics to stderr

Example:
  python3 test/gen_incomplete.py all_gt.tre incomplete_gt.tre --fraction 0.3 --seed 42
"""

import argparse
import random
import sys

# ── Newick parser / pruner ────────────────────────────────────────────────────

def split_top_level(s):
    """Split string at top-level commas (not inside parentheses)."""
    parts, depth, start = [], 0, 0
    for i, c in enumerate(s):
        if   c == '(': depth += 1
        elif c == ')': depth -= 1
        elif c == ',' and depth == 0:
            parts.append(s[start:i])
            start = i + 1
    parts.append(s[start:])
    return parts

def parse_newick(s):
    """
    Parse Newick string into nested dict.
    Returns {'leaf': name} or {'children': [node, ...]}
    Branch lengths and internal node labels are stripped.
    """
    s = s.strip().rstrip(';').strip()
    if s.startswith('('):
        depth = 0
        for i, c in enumerate(s):
            if   c == '(': depth += 1
            elif c == ')':
                depth -= 1
                if depth == 0:
                    inner = s[1:i]
                    return {'children': [parse_newick(c) for c in split_top_level(inner)]}
        raise ValueError(f"Unmatched parenthesis: {s[:60]}")
    else:
        name = s.split(':')[0].strip()
        return {'leaf': name}

def get_leaves(node):
    if 'leaf' in node:
        return [node['leaf']]
    result = []
    for child in node['children']:
        result.extend(get_leaves(child))
    return result

def prune(node, to_remove):
    """
    Remove leaves whose names are in to_remove.
    Returns (pruned_node_or_None, survived).
    Unary internal nodes are collapsed (child promoted).
    """
    if 'leaf' in node:
        if node['leaf'] in to_remove:
            return None, False
        return node, True

    surviving = []
    for child in node['children']:
        c, alive = prune(child, to_remove)
        if alive:
            surviving.append(c)

    if len(surviving) == 0:
        return None, False
    if len(surviving) == 1:
        return surviving[0], True   # collapse unary node
    return {'children': surviving}, True

def to_newick(node):
    if 'leaf' in node:
        return node['leaf']
    return '(' + ','.join(to_newick(c) for c in node['children']) + ')'

# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('input',  help='Input gene tree file (one Newick per line)')
    ap.add_argument('output', help='Output incomplete gene tree file')
    ap.add_argument('--fraction', type=float, default=0.3,
                    help='Fraction of taxa to remove per tree (default 0.3)')
    ap.add_argument('--seed', type=int, default=42,
                    help='Random seed (default 42)')
    ap.add_argument('--min-keep', type=int, default=4,
                    help='Minimum taxa to keep per tree (default 4)')
    ap.add_argument('--stats', action='store_true',
                    help='Print per-tree statistics to stderr')
    args = ap.parse_args()

    rng = random.Random(args.seed)

    with open(args.input) as f:
        lines = [l.strip() for l in f if l.strip()]

    total_in, total_out = 0, 0
    with open(args.output, 'w') as fout:
        for idx, line in enumerate(lines):
            tree   = parse_newick(line)
            leaves = get_leaves(tree)
            n_orig = len(leaves)

            n_keep = max(args.min_keep, int(n_orig * (1.0 - args.fraction)))
            n_keep = min(n_keep, n_orig)   # can't keep more than we have

            keep      = set(rng.sample(leaves, n_keep))
            to_remove = set(leaves) - keep

            pruned, alive = prune(tree, to_remove)
            if not alive:
                # Degenerate: entire tree removed — write original
                pruned = tree

            newick = to_newick(pruned) + ';'
            fout.write(newick + '\n')

            total_in  += n_orig
            total_out += n_keep

            if args.stats:
                print(f"  tree {idx+1:4d}: {n_orig} → {n_keep} taxa "
                      f"(removed {n_orig - n_keep})", file=sys.stderr)

    n_trees = len(lines)
    print(f"Done: {n_trees} trees  "
          f"avg taxa {total_in/n_trees:.1f} → {total_out/n_trees:.1f}  "
          f"fraction removed ≈ {1 - total_out/total_in:.3f}",
          file=sys.stderr)
    print(f"Output: {args.output}", file=sys.stderr)

if __name__ == '__main__':
    main()
