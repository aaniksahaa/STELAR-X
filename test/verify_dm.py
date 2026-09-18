#!/usr/bin/env python3
"""
Independent Python implementation of the STELAR-X distance matrix computation.

Algorithm (same as DistanceMatrixBuilder.buildCPU):
  For each gene tree T, do a post-order traversal.
  At each internal node u with left subtree L and right subtree R:
    for each taxon a in L, taxon b in R:
      dist(a,b) += depth(a) + depth(b) - 2*depth(u)
      cooccurrence(a,b)++
  After all trees: D[a][b] = sum / cooccurrence (or inf if never co-appear).

Usage:
  python3 test/verify_dm.py <input.tre>

Output: same format as STELAR-X --verify-distance-matrix
  DISTANCE_MATRIX
  n=<count>
  taxa=name0,name1,...
  row0=d00,d01,...
  ...

To compare against STELAR-X:
  python3 test/verify_dm.py input.tre > expected.txt
  java ... --verify-distance-matrix --cpu -i input.tre > cpu.txt
  java ... --verify-distance-matrix --gpu -i input.tre > gpu.txt
  python3 test/verify_dm.py --compare expected.txt cpu.txt gpu.txt
"""

import sys
import re
import argparse
from collections import defaultdict

# ── Newick parser ─────────────────────────────────────────────────────────────

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
    Parse a Newick string into a nested dict:
      {'leaf': name}  or  {'children': [node, ...]}
    Branch lengths (:value) are stripped. Internal node labels are ignored.
    """
    s = s.strip().rstrip(';').strip()
    # Remove branch length suffix at the end: everything after the outermost token
    # For a subtree ending in ')label:length' or 'leaf:length' strip the :length part
    if s.startswith('('):
        # Find matching closing ')'
        depth = 0
        for i, c in enumerate(s):
            if c == '(': depth += 1
            elif c == ')':
                depth -= 1
                if depth == 0:
                    inner = s[1:i]
                    children = [parse_newick(c) for c in split_top_level(inner)]
                    return {'children': children}
        raise ValueError(f"Unmatched parenthesis in: {s[:80]}")
    else:
        name = s.split(':')[0].strip()
        return {'leaf': name}

def get_leaves(node):
    if 'leaf' in node:
        return [node['leaf']]
    leaves = []
    for child in node['children']:
        leaves.extend(get_leaves(child))
    return leaves

# ── Distance accumulation ─────────────────────────────────────────────────────

def compute_leaf_depths(node, depth, leaf_depths):
    """DFS: fill leaf_depths[name] = depth."""
    if 'leaf' in node:
        leaf_depths[node['leaf']] = depth
    else:
        for child in node['children']:
            compute_leaf_depths(child, depth + 1, leaf_depths)

def accumulate(node, depth, leaf_depths, dist_sum, cooccur):
    """
    Post-order traversal.
    At each internal node, accumulate distances for all (left taxon, right taxon) pairs.
    For n-ary nodes (polytomies): treat as sequence of splits using running union.
    """
    if 'leaf' in node:
        return [node['leaf']]

    subtrees = []
    for child in node['children']:
        subtrees.append(accumulate(child, depth + 1, leaf_depths, dist_sum, cooccur))

    # For each pair of subtrees (i < j), accumulate all cross-pairs
    all_left = []
    for i, st_i in enumerate(subtrees):
        for j in range(i + 1, len(subtrees)):
            st_j = subtrees[j]
            for a in st_i:
                da = leaf_depths[a]
                for b in st_j:
                    db = leaf_depths[b]
                    d = da + db - 2 * depth
                    key_ab = (a, b) if a < b else (b, a)
                    dist_sum[key_ab]  += d
                    cooccur[key_ab]   += 1
        all_left.extend(st_i)

    all_taxa = []
    for st in subtrees:
        all_taxa.extend(st)
    return all_taxa

def build_dm(trees_newick):
    """
    Returns (taxa_sorted, dist_matrix) where dist_matrix[i][j] = average distance
    between taxa_sorted[i] and taxa_sorted[j].  Returns float('inf') for never-co-appearing pairs.
    """
    # Collect all taxa
    all_taxa = set()
    parsed = []
    for nw in trees_newick:
        tree = parse_newick(nw)
        leaves = get_leaves(tree)
        all_taxa.update(leaves)
        parsed.append(tree)

    taxa = sorted(all_taxa)
    n = len(taxa)
    tid = {name: i for i, name in enumerate(taxa)}

    dist_sum = defaultdict(float)
    cooccur  = defaultdict(int)

    for tree in parsed:
        leaf_depths = {}
        compute_leaf_depths(tree, 0, leaf_depths)
        accumulate(tree, 0, leaf_depths, dist_sum, cooccur)

    # Build matrix
    mat = [[float('inf')] * n for _ in range(n)]
    for i in range(n):
        mat[i][i] = 0.0

    for (a, b), s in dist_sum.items():
        i, j = tid[a], tid[b]
        d = s / cooccur[(a, b)]
        mat[i][j] = d
        mat[j][i] = d

    return taxa, mat

# ── Output / compare ──────────────────────────────────────────────────────────

def print_dm(taxa, mat, file=sys.stdout):
    n = len(taxa)
    print("DISTANCE_MATRIX", file=file)
    print(f"n={n}", file=file)
    print("taxa=" + ",".join(taxa), file=file)
    for i in range(n):
        vals = []
        for j in range(n):
            d = mat[i][j]
            vals.append("inf" if d == float('inf') else f"{d:.6f}")
        print(f"row{i}=" + ",".join(vals), file=file)

def parse_dm_file(path):
    """Parse a DISTANCE_MATRIX file into (taxa_list, flat_matrix_dict)."""
    with open(path) as f:
        lines = [l.strip() for l in f if l.strip()]

    assert lines[0] == "DISTANCE_MATRIX", f"Bad header in {path}"
    n = int(lines[1].split('=')[1])
    taxa = lines[2].split('=')[1].split(',')

    rows = {}
    for line in lines[3:]:
        m = re.match(r'row(\d+)=(.+)', line)
        if m:
            i = int(m.group(1))
            vals = m.group(2).split(',')
            rows[i] = [float('inf') if v == 'inf' else float(v) for v in vals]

    return taxa, rows, n

def compare_dms(files, tol=1e-6):
    """
    Compare multiple DM files for agreement within tolerance.
    Order-independent: taxa are matched by name, not by position.
    Returns True if all files match the first within tol.
    """
    parsed = [parse_dm_file(f) for f in files]
    ref_taxa, ref_rows, ref_n = parsed[0]
    ref_tid = {name: i for i, name in enumerate(ref_taxa)}

    ok = True
    for path, (taxa, rows, n) in zip(files[1:], parsed[1:]):
        if n != ref_n:
            print(f"MISMATCH: n={ref_n} vs n={n}  ({files[0]} vs {path})")
            ok = False
            continue
        if set(taxa) != set(ref_taxa):
            print(f"MISMATCH: taxon sets differ  ({files[0]} vs {path})")
            ok = False
            continue

        tid = {name: i for i, name in enumerate(taxa)}
        max_diff = 0.0
        mismatches = 0

        for name_a in ref_taxa:
            for name_b in ref_taxa:
                ri = ref_tid[name_a]; rj = ref_tid[name_b]
                ci = tid[name_a];    cj = tid[name_b]
                a = ref_rows[ri][rj]
                b = rows[ci][cj]
                if a == float('inf') and b == float('inf'):
                    continue
                if a == float('inf') or b == float('inf'):
                    print(f"  MISMATCH ({name_a},{name_b}): inf vs {b if a==float('inf') else a}")
                    mismatches += 1
                    continue
                diff = abs(a - b)
                if diff > max_diff:
                    max_diff = diff
                if diff > tol:
                    print(f"  MISMATCH ({name_a},{name_b}): {a:.9f} vs {b:.9f}  diff={diff:.2e}")
                    mismatches += 1

        if mismatches > 0:
            print(f"  FAIL  {mismatches} mismatches  {files[0]} vs {path}")
            ok = False
        else:
            print(f"  OK    max_diff={max_diff:.2e}  (tol={tol:.0e})  {files[0]} vs {path}")

    return ok

# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('input', nargs='?', help='Gene tree file (Newick, one tree per line)')
    ap.add_argument('--compare', nargs='+', metavar='FILE',
                    help='Compare two or more DISTANCE_MATRIX files for agreement')
    ap.add_argument('--tol', type=float, default=1e-6,
                    help='Tolerance for floating-point comparison (default 1e-6)')
    args = ap.parse_args()

    if args.compare:
        ok = compare_dms(args.compare, tol=args.tol)
        sys.exit(0 if ok else 1)

    if not args.input:
        ap.print_help()
        sys.exit(1)

    with open(args.input) as f:
        trees = [l.strip() for l in f if l.strip()]

    taxa, mat = build_dm(trees)
    print_dm(taxa, mat)

if __name__ == '__main__':
    main()
