#!/usr/bin/env python3
"""
test_similarity_matrix.py
=========================
Randomized tests for the STELAR-X similarity matrix computation.

For each test case:
  1. Generate n random taxa and k random binary gene trees (some incomplete).
  2. Compute the expected similarity matrix by brute-force in Python.
  3. Run STELAR-X with --verify-similarity-matrix and parse its output.
  4. Assert all values match within tolerance.

Brute-force definition (binary trees):
  For each tree t (with kt leaves, kt >= 4), pair (a,b), and other pair
  (c,d), use the four-point condition on independently computed tree distances
  to count whether the displayed quartet is ab|cd.  den_t = C2(kt - 2).
  sim(a,b)  = sum_t(num_t) / sum_t(den_t)  if sum_t(den_t) > 0 else 0.0
  sim(a,a)  = 1.0

Usage:
  python3 test/test_similarity_matrix.py [--stelarx-root PATH] [--mode cpu|gpu]
                                          [--seeds N N N ...] [-n NUM_TESTS] [-v]
"""

import argparse
import math
import os
import random
import subprocess
import sys
import tempfile


# ── Tree representation ───────────────────────────────────────────────────────

class Node:
    __slots__ = ("left", "right", "taxon", "parent", "_sub_lc")

    def __init__(self, left=None, right=None, taxon=None):
        self.left   = left
        self.right  = right
        self.taxon  = taxon
        self.parent = None
        self._sub_lc = None

    def is_leaf(self):
        return self.left is None


def random_binary_tree(taxa):
    """Return a random binary Node tree over the given taxa list."""
    taxa = list(taxa)
    if len(taxa) == 1:
        return Node(taxon=taxa[0])
    random.shuffle(taxa)
    k = random.randint(1, len(taxa) - 1)
    left  = random_binary_tree(taxa[:k])
    right = random_binary_tree(taxa[k:])
    node  = Node(left=left, right=right)
    left.parent  = node
    right.parent = node
    return node


def to_newick(node):
    if node.is_leaf():
        return node.taxon
    return f"({to_newick(node.left)},{to_newick(node.right)})"


# ── Sub-leaf-count ─────────────────────────────────────────────────────────

def compute_sub_lc(node):
    """Annotate every node with _sub_lc = number of leaf descendants."""
    if node is None:
        return 0
    if node.is_leaf():
        node._sub_lc = 1
        return 1
    lc = compute_sub_lc(node.left) + compute_sub_lc(node.right)
    node._sub_lc = lc
    return lc


# ── LCA ───────────────────────────────────────────────────────────────────────

def leaf_map(node):
    """Returns dict: taxon_name -> Node (leaf node)."""
    if node.is_leaf():
        return {node.taxon: node}
    result = {}
    result.update(leaf_map(node.left))
    result.update(leaf_map(node.right))
    return result


def lca(a_node, b_node):
    """Find LCA of two nodes using ancestor-set approach."""
    ancestors = set()
    cur = a_node
    while cur is not None:
        ancestors.add(id(cur))
        cur = cur.parent
    cur = b_node
    while cur is not None:
        if id(cur) in ancestors:
            return cur
        cur = cur.parent
    return None


def tree_distance(a_node, b_node):
    """Unit-edge distance, implemented independently of STELAR-X Euler data."""
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


# ── C2 helper ─────────────────────────────────────────────────────────────────

def c2(x):
    return x * (x - 1) // 2 if x >= 2 else 0


# ── Brute-force similarity matrix ─────────────────────────────────────────────

def brute_force_similarity_matrix(trees_with_taxa, all_taxa):
    """
    trees_with_taxa: list of (Node root, set of taxon names present)
    all_taxa:        ordered list of all taxa (index = taxon id)

    Returns n×n list-of-lists of similarities in [0,1].
    """
    n   = len(all_taxa)
    idx = {t: i for i, t in enumerate(all_taxa)}

    num_sum = [[0.0] * n for _ in range(n)]
    den_sum = [[0.0] * n for _ in range(n)]

    for root, present in trees_with_taxa:
        kt = len(present)
        if kt < 4:
            continue   # C2(kt-2) = 0 for kt <= 3

        leaves = leaf_map(root)
        den_t  = c2(kt - 2)
        present_list = sorted(present)
        distances = {}
        for i, a in enumerate(present_list):
            for b in present_list[i + 1:]:
                distances[(a, b)] = distances[(b, a)] = tree_distance(leaves[a], leaves[b])

        for pi, a in enumerate(present_list):
            for b in present_list[pi + 1:]:
                ia, ib  = idx[a], idx[b]
                # Every co-occurring pair receives the denominator, including
                # pairs whose LCA is the root and therefore contribute zero to
                # the numerator. Skipping those zero-numerator observations
                # incorrectly turns many ratios into 1.0.
                den_sum[ia][ib] += den_t
                den_sum[ib][ia] += den_t
                others = [taxon for taxon in present_list if taxon != a and taxon != b]
                num_t = 0
                for ci, c in enumerate(others):
                    for d in others[ci + 1:]:
                        ab_cd = distances[(a, b)] + distances[(c, d)]
                        ac_bd = distances[(a, c)] + distances[(b, d)]
                        ad_bc = distances[(a, d)] + distances[(b, c)]
                        if ab_cd < ac_bd and ab_cd < ad_bc:
                            num_t += 1
                if num_t == 0:
                    continue
                num_sum[ia][ib] += num_t
                num_sum[ib][ia] += num_t

    result = [[0.0] * n for _ in range(n)]
    for i in range(n):
        result[i][i] = 1.0
        for j in range(n):
            if i == j:
                continue
            result[i][j] = num_sum[i][j] / den_sum[i][j] if den_sum[i][j] > 0 else 0.0

    return result


# ── Run STELAR-X and parse output ─────────────────────────────────────────────

def run_stelarx_similarity_matrix(newick_lines, stelarx_root, mode):
    build_dir  = os.path.join(stelarx_root, "build")
    native_dir = os.path.join(stelarx_root, "native")

    with tempfile.NamedTemporaryFile(mode="w", suffix=".tre", delete=False) as f:
        f.write("\n".join(newick_lines) + "\n")
        input_path = f.name

    try:
        cmd = ["java",
               f"-Djava.library.path={native_dir}",
               "-cp", build_dir,
               "stelarx.Main",
               "-i", input_path,
               "--verify-similarity-matrix", "-q"]
        if mode == "gpu":
            cmd.append("--gpu-strict")
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if result.returncode != 0:
            raise RuntimeError(
                f"STELAR-X exited {result.returncode}\nstderr:\n{result.stderr[:2000]}")
        return parse_sm_output(result.stdout)
    finally:
        os.unlink(input_path)


def parse_sm_output(output):
    """
    Parse SIMILARITY_MATRIX block from stdout.
    Returns (taxa_order: list[str], matrix: list[list[float]]).
    """
    lines      = [l.strip() for l in output.splitlines() if l.strip()]
    taxa_order = []
    matrix     = []
    in_sm      = False

    for line in lines:
        if line == "SIMILARITY_MATRIX":
            in_sm = True
            continue
        if not in_sm:
            continue
        if line.startswith("n="):
            pass
        elif line.startswith("taxa="):
            taxa_order = line[5:].split(",")
        elif line.startswith("sim_row"):
            row = [float(v) for v in line.split("=", 1)[1].split(",")]
            matrix.append(row)

    if not taxa_order or not matrix:
        raise RuntimeError(
            f"Failed to parse SIMILARITY_MATRIX from output:\n{output[:1000]}")
    return taxa_order, matrix


# ── Test runner ───────────────────────────────────────────────────────────────

TOLERANCE = 1e-6


def run_test(seed, stelarx_root, mode, verbose=False):
    random.seed(seed)
    n_taxa    = random.randint(5, 40)
    k_trees   = random.randint(5, 20)
    miss_rate = random.uniform(0.0, 0.4)

    all_taxa = [f"t{i}" for i in range(n_taxa)]

    trees_with_taxa = []
    newick_lines    = []

    for _ in range(k_trees):
        n_drop  = int(n_taxa * miss_rate)
        present = all_taxa[:]
        if n_drop > 0 and len(present) - n_drop >= 2:
            present = [t for t in present if t not in random.sample(present, n_drop)]
        if len(present) < 2:
            present = all_taxa[:]
        root = random_binary_tree(present)
        trees_with_taxa.append((root, set(present)))
        newick_lines.append(to_newick(root) + ";")

    # Python brute-force answer
    expected = brute_force_similarity_matrix(trees_with_taxa, all_taxa)

    # STELAR-X answer
    taxa_order, actual_raw = run_stelarx_similarity_matrix(newick_lines, stelarx_root, mode)

    # Reorder actual to match all_taxa order
    idx_java = {t: i for i, t in enumerate(taxa_order)}
    n = len(all_taxa)
    actual = [[0.0] * n for _ in range(n)]
    for i, ta in enumerate(all_taxa):
        for j, tb in enumerate(all_taxa):
            ji, jj = idx_java[ta], idx_java[tb]
            actual[i][j] = actual_raw[ji][jj]

    # Compare
    max_err = 0.0
    worst   = (0, 0)
    for i in range(n):
        for j in range(n):
            err = abs(expected[i][j] - actual[i][j])
            if err > max_err:
                max_err = err
                worst   = (i, j)

    passed = max_err <= TOLERANCE
    if verbose or not passed:
        status = "PASS" if passed else "FAIL"
        print(f"  seed={seed:5d}  n={n_taxa:3d}  k={k_trees:3d}  "
              f"miss_rate={miss_rate:.2f}  max_err={max_err:.2e}  [{status}]")
        if not passed:
            i, j = worst
            print(f"    worst pair: ({all_taxa[i]},{all_taxa[j]})  "
                  f"expected={expected[i][j]:.8f}  actual={actual[i][j]:.8f}")
    return passed


def main():
    ap = argparse.ArgumentParser(description="Test STELAR-X similarity matrix")
    ap.add_argument("--stelarx-root", default=".", help="Path to STELAR-X root")
    ap.add_argument("--mode", choices=["cpu", "gpu"], default="cpu")
    ap.add_argument("--seeds", type=int, nargs="+",
                    help="Specific seeds to test (default: 1..20)")
    ap.add_argument("-n", "--num-tests", type=int, default=20,
                    help="Number of tests (ignored if --seeds given)")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    seeds = args.seeds if args.seeds else list(range(1, args.num_tests + 1))
    print(f"Running {len(seeds)} similarity-matrix tests (mode={args.mode}) ...")

    passed = failed = 0
    for seed in seeds:
        try:
            ok = run_test(seed, args.stelarx_root, args.mode, verbose=args.verbose)
            if ok:
                passed += 1
            else:
                failed += 1
        except Exception as exc:
            print(f"  seed={seed}  ERROR: {exc}")
            failed += 1

    total = passed + failed
    print(f"\n{passed}/{total} tests passed", end="")
    if failed:
        print(f"  ({failed} FAILED)")
        sys.exit(1)
    else:
        print("  ✓")


if __name__ == "__main__":
    main()
