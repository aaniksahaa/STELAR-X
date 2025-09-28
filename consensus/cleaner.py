#!/usr/bin/env python3
"""
resolve_clean_dendropy.py

Robustly read one or more Newick files (each file may contain multiple trees separated by semicolons),
resolve polytomies (random or deterministic), strip branch lengths and internal labels (by default),
and write cleaned Newick (multiple trees concatenated).

This version uses a conservative per-tree parser: it splits input on ';' and parses each tree string
separately with dendropy.Tree.get(..., data=...).
"""
import sys
sys.setrecursionlimit(10_000_000)

import argparse
import random
from typing import Optional, List

import dendropy

def read_trees_from_file_robust(path: str) -> List[dendropy.Tree]:
    """
    Read a Newick file that may contain multiple trees separated by semicolons.
    This function reads the file as text, splits on ';', and attempts to parse each
    non-empty chunk as a single Newick tree. Returns a list of parsed dendropy.Tree objects.
    """
    txt = open(path, "r", encoding="utf-8").read()
    # Split on semicolon. Keep everything before each semicolon as a tree chunk.
    chunks = [c.strip() for c in txt.split(';')]
    trees = []
    for i, chunk in enumerate(chunks):
        if not chunk:
            continue
        tree_str = chunk + ';'  # ensure terminating semicolon
        try:
            # parse single tree from string
            t = dendropy.Tree.get(data=tree_str, schema="newick", preserve_underscores=True)
            trees.append(t)
        except Exception as e:
            # Try one more attempt with relaxed whitespace (some weird inputs have stray BOMs or control chars)
            try:
                cleaned = ''.join(ch for ch in tree_str if ord(ch) >= 9)
                t = dendropy.Tree.get(data=cleaned, schema="newick", preserve_underscores=True)
                trees.append(t)
            except Exception as e2:
                print(f"Warning: failed to parse tree #{i+1} in file '{path}': {e2}", file=sys.stderr)
                # skip this chunk but continue
                continue
    return trees

def resolve_polytomies_dendropy(tree: dendropy.Tree, mode: str = "random", new_edge_length: Optional[float] = None, rnd: Optional[random.Random] = None):
    if mode not in ("random", "firstpair"):
        raise ValueError("mode must be 'random' or 'firstpair'")

    nodes = list(tree.postorder_node_iter())
    for node in nodes:
        while len(list(node.child_nodes())) > 2:
            children = list(node.child_nodes())
            if mode == "firstpair":
                a = children[0]
                b = children[1]
            else:
                if rnd is None:
                    i, j = random.sample(range(len(children)), 2)
                else:
                    i, j = rnd.sample(range(len(children)), 2)
                a = children[i]
                b = children[j]

            node.remove_child(a)
            node.remove_child(b)

            new_internal = dendropy.Node()
            new_internal.add_child(a)
            new_internal.add_child(b)

            if new_edge_length is not None:
                new_internal.edge.length = float(new_edge_length)

            node.add_child(new_internal)

def strip_lengths_and_labels_dendropy(tree: dendropy.Tree, strip_internal_labels: bool = True):
    for e in tree.postorder_edge_iter():
        try:
            e.length = None
        except Exception:
            pass

    for n in tree.postorder_node_iter():
        if n.is_leaf():
            continue
        if strip_internal_labels:
            n.label = None
        if hasattr(n, "annotations") and n.annotations is not None:
            try:
                n.annotations.clear()
            except Exception:
                pass

def degree_freq_and_tips(tree: dendropy.Tree):
    freq = {}
    tips = 0
    for n in tree.preorder_node_iter():
        if n.is_leaf():
            tips += 1
        else:
            d = len(list(n.child_nodes()))
            freq[d] = freq.get(d, 0) + 1
    return tips, freq

def process_trees(trees: List[dendropy.Tree], args, rnd: Optional[random.Random] = None):
    cleaned_strings = []
    for idx, tree in enumerate(trees, start=1):
        print(f"\n--- Tree #{idx} (starting) ---")
        tips_before, freq_before = degree_freq_and_tips(tree)
        print(f"Tips before: {tips_before}")
        print("Degree freq before (deg -> count):", dict(sorted(freq_before.items())))

        print("Resolving polytomies...")
        resolve_polytomies_dendropy(tree, mode=args.mode, new_edge_length=args.new_edge_length, rnd=rnd)

        tips_after, freq_after = degree_freq_and_tips(tree)
        print(f"Tips after resolve: {tips_after}")
        print("Degree freq after resolve (deg -> count):", dict(sorted(freq_after.items())))

        print("Stripping branch lengths and internal labels..." if args.strip_internal_labels else "Stripping branch lengths (keeping internal labels)...")
        strip_lengths_and_labels_dendropy(tree, strip_internal_labels=args.strip_internal_labels)

        tips_final, freq_final = degree_freq_and_tips(tree)
        print(f"Tips final: {tips_final}")
        print("Degree freq final (deg -> count):", dict(sorted(freq_final.items())))

        s = tree.as_string(schema="newick", suppress_rooting=True).strip()
        if not s.endswith(";"):
            s = s + ";"
        cleaned_strings.append(s)
        print(f"--- Tree #{idx} processed ---")
    return cleaned_strings

def main():
    p = argparse.ArgumentParser(description="Resolve polytomies and clean Newick using DendroPy.")
    p.add_argument("input", nargs="+", help="Input Newick file(s). Each file may contain multiple trees separated by semicolons.")
    p.add_argument("output", help="Output cleaned Newick file (multiple trees concatenated).")
    p.add_argument("--new-edge-length", type=float, default=None, help="Edge length for newly created internal edges (default: None)")
    p.add_argument("--no-strip-internal-labels", dest="strip_internal_labels", action="store_false", help="Do NOT strip internal node labels (by default internal labels are removed).")
    p.add_argument("--deterministic", dest="deterministic", action="store_true", help="Resolve polytomies deterministically (first-two children) instead of randomly.")
    p.add_argument("--random-seed", type=int, default=None, help="If provided and resolution is random, use this seed for reproducibility.")
    args = p.parse_args()

    args.mode = "firstpair" if args.deterministic else "random"

    rnd = None
    if args.mode == "random":
        if args.random_seed is not None:
            rnd = random.Random(args.random_seed)
            print(f"Random mode with seed {args.random_seed}")
        else:
            rnd = random.Random()
            print("Random mode with no fixed seed (non-deterministic)")

    all_cleaned = []
    total = 0
    for inp in args.input:
        print(f"\nReading file: {inp}")
        try:
            trees = read_trees_from_file_robust(inp)
        except Exception as e:
            print(f"Error reading '{inp}': {e}", file=sys.stderr)
            continue
        print(f"Found {len(trees)} tree(s) in {inp}.")
        total += len(trees)
        cleaned = process_trees(trees, args, rnd=rnd)
        all_cleaned.extend(cleaned)

    if total == 0:
        print("No trees processed. Exiting.", file=sys.stderr)
        sys.exit(1)

    print(f"\nWriting {len(all_cleaned)} cleaned tree(s) to: {args.output}")
    with open(args.output, "w", encoding="utf-8") as fh:
        for i, s in enumerate(all_cleaned, start=1):
            fh.write(s + "\n")
            # if i != len(all_cleaned):
            #     fh.write("\n")
    print("Done.")

if __name__ == "__main__":
    main()
