#!/usr/bin/env python3
"""
test_upgma.py
=============
Randomized tests for the STELAR-X UPGMA clusterer.

For each test case:
  1. Generate n random taxa and k random binary gene trees (some incomplete).
  2. Compute the similarity matrix by brute-force in Python.
  3. Run Python O(n³) UPGMA brute-force on that matrix → Python bipartitions.
  4. Run STELAR-X with --verify-upgma → STELAR-X bipartitions.
  5. Assert the two bipartition *sets* are identical.

Python UPGMA is deliberately O(n³): find-best in O(n²) each of n−1 iterations.
This is correct but slow; it serves only as a reference for n ≤ ~50.

Usage:
  python3 test/test_upgma.py [--stelarx-root PATH] [--seeds N N N ...] [-n NUM] [-v]
"""

import argparse
import os
import random
import subprocess
import sys
import tempfile

# ── Reuse tree helpers from test_similarity_matrix.py ────────────────────────

class Node:
    __slots__ = ("left", "right", "taxon", "parent", "_sub_lc")
    def __init__(self, left=None, right=None, taxon=None):
        self.left = left; self.right = right
        self.taxon = taxon; self.parent = None; self._sub_lc = None
    def is_leaf(self): return self.left is None


def random_binary_tree(taxa):
    taxa = list(taxa)
    if len(taxa) == 1:
        return Node(taxon=taxa[0])
    random.shuffle(taxa)
    k     = random.randint(1, len(taxa) - 1)
    left  = random_binary_tree(taxa[:k])
    right = random_binary_tree(taxa[k:])
    node  = Node(left=left, right=right)
    left.parent = right.parent = node
    return node


def to_newick(node):
    return node.taxon if node.is_leaf() else f"({to_newick(node.left)},{to_newick(node.right)})"


def compute_sub_lc(node):
    if node is None: return 0
    if node.is_leaf(): node._sub_lc = 1; return 1
    lc = compute_sub_lc(node.left) + compute_sub_lc(node.right)
    node._sub_lc = lc; return lc


def leaf_map(node):
    if node.is_leaf(): return {node.taxon: node}
    r = {}; r.update(leaf_map(node.left)); r.update(leaf_map(node.right)); return r


def lca(a_node, b_node):
    anc = set()
    cur = a_node
    while cur: anc.add(id(cur)); cur = cur.parent
    cur = b_node
    while cur:
        if id(cur) in anc: return cur
        cur = cur.parent
    return None


def tree_distance(a_node, b_node):
    ancestors = {}
    distance = 0
    cur = a_node
    while cur is not None:
        ancestors[id(cur)] = distance
        distance += 1
        cur = cur.parent
    distance = 0
    cur = b_node
    while cur is not None:
        if id(cur) in ancestors:
            return ancestors[id(cur)] + distance
        distance += 1
        cur = cur.parent
    raise AssertionError("disconnected tree")


def c2(x): return x * (x - 1) // 2 if x >= 2 else 0


# ── Python brute-force similarity matrix ─────────────────────────────────────

def brute_force_similarity(trees_with_taxa, all_taxa):
    """Returns n×n list-of-lists of similarity values in [0,1]."""
    n   = len(all_taxa)
    idx = {t: i for i, t in enumerate(all_taxa)}
    num = [[0.0]*n for _ in range(n)]
    den = [[0.0]*n for _ in range(n)]
    for root, present in trees_with_taxa:
        kt = len(present)
        if kt < 4: continue
        leaves = leaf_map(root)
        den_t  = c2(kt - 2)
        plist  = sorted(present)
        distances = {}
        for i, a in enumerate(plist):
            for b in plist[i + 1:]:
                distances[(a, b)] = distances[(b, a)] = tree_distance(leaves[a], leaves[b])
        for pi, a in enumerate(plist):
            for b in plist[pi+1:]:
                ia, ib = idx[a], idx[b]
                den[ia][ib] += den_t; den[ib][ia] += den_t
                others = [taxon for taxon in plist if taxon != a and taxon != b]
                num_t = 0
                for ci, c in enumerate(others):
                    for d in others[ci + 1:]:
                        ab_cd = distances[(a, b)] + distances[(c, d)]
                        ac_bd = distances[(a, c)] + distances[(b, d)]
                        ad_bc = distances[(a, d)] + distances[(b, c)]
                        if ab_cd < ac_bd and ab_cd < ad_bc:
                            num_t += 1
                if num_t == 0: continue
                num[ia][ib] += num_t; num[ib][ia] += num_t
    sim = [[0.0]*n for _ in range(n)]
    for i in range(n):
        sim[i][i] = 1.0
        for j in range(n):
            if i != j:
                sim[i][j] = num[i][j] / den[i][j] if den[i][j] > 0 else 0.0
    return sim


# ── Python O(n³) UPGMA brute-force ───────────────────────────────────────────

def upgma_bipartitions(sim, ordered_taxa):
    """
    O(n³) UPGMA reference implementation.

    sim         — n×n list-of-lists indexed by position in ordered_taxa.
    ordered_taxa — taxon names in the order used for the sim matrix rows/cols;
                   this determines tie-breaking when similarities are equal.

    Returns a set of frozensets of taxon *names* (strings), one per internal
    non-root merge (size 2..n-1).  The all-taxa cluster (size n) is excluded.
    """
    n = len(ordered_taxa)
    # Working similarity dictionary: (min_i, max_i) → float
    d = {}
    for i in range(n):
        for j in range(i+1, n):
            d[(i, j)] = sim[i][j]

    clusters = [frozenset([ordered_taxa[i]]) for i in range(n)]
    weights  = [1] * n
    active   = list(range(n))
    bips     = set()

    for _ in range(n - 1):
        # O(n²) find-best among active pairs
        best_s = -float('inf')
        I, J   = -1, -1
        for ii in range(len(active)):
            for jj in range(ii + 1, len(active)):
                i, j = active[ii], active[jj]
                key = (min(i, j), max(i, j))
                if d[key] > best_s:
                    best_s = d[key]; I, J = i, j

        wI, wJ   = weights[I], weights[J]
        merged   = clusters[I] | clusters[J]

        # Record bipartition unless it is the all-taxa cluster
        if len(merged) < n:
            bips.add(merged)

        # O(n) weighted-average update of similarities
        for k in active:
            if k == I or k == J: continue
            ki   = (min(k, I), max(k, I))
            kj   = (min(k, J), max(k, J))
            d[ki] = (d[ki] * wI + d[kj] * wJ) / (wI + wJ)

        clusters[I] = merged
        weights[I]  = wI + wJ
        active.remove(J)

    return bips


# ── Replicate TaxonRegistry first-seen ordering ───────────────────────────────

def first_seen_order(newick_lines):
    """
    Return taxa names in the order they first appear (left-to-right) across
    all Newick strings, exactly as STELAR-X's TaxonRegistry assigns IDs.
    """
    seen = {}
    for line in newick_lines:
        token = []
        for ch in line:
            if ch in '(),;':
                if token:
                    name = ''.join(token)
                    if name not in seen:
                        seen[name] = len(seen)
                    token = []
            elif ch not in ' \t\n':
                token.append(ch)
        if token:
            name = ''.join(token)
            if name not in seen:
                seen[name] = len(seen)
    return list(seen.keys())   # ordered by first appearance


# ── Parse STELAR-X --verify-upgma output ─────────────────────────────────────

def parse_stelarx_upgma(output):
    """
    Parse stdout from STELAR-X --verify-upgma.

    Returns a set of frozensets of taxon *names* (strings).
    """
    bips = set()
    for line in output.splitlines():
        if not line.startswith("bipartition="):
            continue
        names = [nm for nm in line[len("bipartition="):].split(",") if nm]
        bips.add(frozenset(names))
    return bips


# ── Random test-case generation ───────────────────────────────────────────────

def make_test_case(seed, n_taxa, k_trees, incomplete_frac=0.3):
    """
    Returns (all_taxa, gene_trees_newick, trees_with_taxa_for_python).
    gene_trees_newick is a list of Newick strings.
    """
    rng = random.Random(seed)
    all_taxa = [f"t{i}" for i in range(n_taxa)]

    newick_lines    = []
    trees_with_taxa = []

    for _ in range(k_trees):
        if rng.random() < incomplete_frac and n_taxa >= 4:
            # Drop 1..n//3 taxa
            drop = rng.randint(1, max(1, n_taxa // 3))
            taxa = rng.sample(all_taxa, n_taxa - drop)
        else:
            taxa = all_taxa[:]

        root = random_binary_tree(taxa)
        # Fix parent pointers (random_binary_tree already sets them)
        newick_lines.append(to_newick(root) + ";")
        trees_with_taxa.append((root, set(taxa)))

    return all_taxa, newick_lines, trees_with_taxa


# ── Single test run ───────────────────────────────────────────────────────────

def run_one(seed, n_taxa, k_trees, stelarx_root, verbose):
    random.seed(seed)  # ensure global random state is deterministic per seed
    all_taxa, newick_lines, trees_for_python = make_test_case(seed, n_taxa, k_trees)

    # Determine first-seen taxon order (mirrors STELAR-X TaxonRegistry)
    ordered_taxa = first_seen_order(newick_lines)

    # ── Python reference ──────────────────────────────────────────────────────
    # Build sim matrix indexed by all_taxa (0..n-1 = t0..t_{n-1})
    sim_raw = brute_force_similarity(trees_for_python, all_taxa)
    raw_idx = {t: i for i, t in enumerate(all_taxa)}

    # Reorder sim matrix to first-seen ordering so tie-breaking matches STELAR-X
    n = len(ordered_taxa)
    sim = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            ri = raw_idx[ordered_taxa[i]]
            rj = raw_idx[ordered_taxa[j]]
            sim[i][j] = sim_raw[ri][rj]

    py_bips = upgma_bipartitions(sim, ordered_taxa)

    # ── STELAR-X ──────────────────────────────────────────────────────────────
    build_dir  = os.path.join(stelarx_root, "build")
    native_dir = os.path.join(stelarx_root, "native")
    with tempfile.NamedTemporaryFile(mode='w', suffix='.tre', delete=False) as f:
        f.write('\n'.join(newick_lines) + '\n')
        tmp_tree = f.name

    try:
        cmd = ["java",
               f"-Djava.library.path={native_dir}",
               "-cp", build_dir,
               "stelarx.Main",
               "-i", tmp_tree, "--verify-upgma", "-q"]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    finally:
        os.unlink(tmp_tree)

    if result.returncode != 0:
        if verbose:
            print(f"  STELAR-X stderr:\n{result.stderr[:500]}")
        return False, f"STELAR-X exited with code {result.returncode}"

    ax_bips = parse_stelarx_upgma(result.stdout)

    # ── Compare (by taxon names, not indices) ─────────────────────────────────
    if py_bips == ax_bips:
        return True, ""

    # Build readable diff
    only_py = py_bips - ax_bips
    only_ax = ax_bips - py_bips
    msg = []
    if only_py:
        msg.append(f"  in Python only: {[sorted(s) for s in only_py]}")
    if only_ax:
        msg.append(f"  in STELAR-X only: {[sorted(s) for s in only_ax]}")
    return False, "\n".join(msg)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description="Test STELAR-X UPGMA clusterer")
    ap.add_argument("--stelarx-root", default=os.path.join(os.path.dirname(__file__), ".."),
                    help="STELAR-X project root (default: parent of test/)")
    ap.add_argument("--seeds", type=int, nargs="+",
                    help="Specific seeds to test (default: 1..20)")
    ap.add_argument("-n", "--num-tests", type=int, default=20)
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    stelarx_root = os.path.abspath(args.stelarx_root)
    seeds = args.seeds if args.seeds else list(range(1, args.num_tests + 1))

    print(f"Running {len(seeds)} UPGMA tests ...")

    passed = 0
    for seed in seeds:
        rng = random.Random(seed)
        # Vary n (5..30) and k (8..20) per seed for coverage
        n_taxa  = rng.randint(5, 30)
        k_trees = rng.randint(8, 20)

        ok, msg = run_one(seed, n_taxa, k_trees, stelarx_root, args.verbose)
        if ok:
            passed += 1
            if args.verbose:
                print(f"  seed={seed:3d}  n={n_taxa}  k={k_trees}  PASS")
        else:
            print(f"  seed={seed:3d}  n={n_taxa}  k={k_trees}  FAIL")
            print(msg)

    total = len(seeds)
    print(f"\n{passed}/{total} tests passed  {'✓' if passed == total else '✗'}")
    sys.exit(0 if passed == total else 1)


if __name__ == "__main__":
    main()
