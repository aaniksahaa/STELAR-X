#!/usr/bin/env python3
"""
resolve_polytomy_and_clean.py

Read one or more Newick files (each file may contain multiple trees separated by semicolons),
resolve polytomies (random or deterministic), strip branch lengths and internal labels (by default),
and write cleaned Newick (multiple trees concatenated).

This script intentionally does NOT reroot trees. Rooting is handled separately by root_by_outgroups.py.
"""
import sys
sys.setrecursionlimit(10_000_000)

import argparse
import os
import random
from concurrent.futures import ProcessPoolExecutor
from typing import Optional, List, Tuple

import dendropy


class Colors:
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    RED = '\033[91m'
    YELLOW = '\033[93m'
    BOLD = '\033[1m'
    RESET = '\033[0m'


def iter_tree_strings_from_file_robust(path: str):
    """
    Generator that reads a Newick file containing multiple trees separated by semicolons.
    Yields (tree_index, tree_str) tuples one at a time to minimize memory usage.
    """
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

                tree_str = tree_str + ';'
                tree_index += 1
                yield (tree_index, tree_str)

    remaining = buffer.strip()
    if remaining:
        tree_str = remaining if remaining.endswith(';') else remaining + ';'
        tree_index += 1
        yield (tree_index, tree_str)


def parse_tree_from_string(tree_str: str, tree_idx: int, path: str) -> Optional[dendropy.Tree]:
    try:
        return dendropy.Tree.get(data=tree_str, schema="newick", preserve_underscores=True)
    except Exception:
        try:
            cleaned = ''.join(ch for ch in tree_str if ord(ch) >= 9)
            return dendropy.Tree.get(data=cleaned, schema="newick", preserve_underscores=True)
        except Exception as e2:
            print(f"Warning: failed to parse tree #{tree_idx} in file '{path}': {e2}", file=sys.stderr)
            return None


def resolve_polytomies_dendropy(
    tree: dendropy.Tree,
    mode: str = "random",
    new_edge_length: Optional[float] = None,
    rnd: Optional[random.Random] = None,
    min_children: int = 2,
):
    if mode not in ("random", "firstpair"):
        raise ValueError("mode must be 'random' or 'firstpair'")

    nodes = list(tree.postorder_node_iter())
    for node in nodes:
        while len(list(node.child_nodes())) > min_children:
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


def warn_if_polytomies_remain(freq: dict, tree_idx: int):
    remaining = {deg: cnt for deg, cnt in freq.items() if deg > 2}
    if not remaining:
        return
    total_nodes = sum(remaining.values())
    degrees = ", ".join(f"{deg}->{cnt}" for deg, cnt in sorted(remaining.items()))
    print(
        f"{Colors.YELLOW}{Colors.BOLD}⚠ Warning:{Colors.RESET} "
        f"{Colors.YELLOW}Tree #{tree_idx} still has {total_nodes} polytomy node(s) "
        f"(degree>2) after resolution: {degrees}{Colors.RESET}"
    )


def process_single_tree_from_string(
    tree_str: str,
    idx: int,
    path: str,
    mode: str,
    new_edge_length: Optional[float],
    strip_internal_labels: bool,
    base_seed: Optional[int],
) -> Tuple[Optional[str], str]:
    print(f"\n--- Tree #{idx} (starting) ---")
    tree = parse_tree_from_string(tree_str, idx, path)
    if tree is None:
        return None, ""

    tips_before, freq_before = degree_freq_and_tips(tree)
    print(f"Tips before: {tips_before}")
    print("Degree freq before (deg -> count):", dict(sorted(freq_before.items())))

    print("Resolving polytomies (>2 children)...")
    rnd = random.Random(base_seed + idx) if (mode == "random" and base_seed is not None) else None
    resolve_polytomies_dendropy(tree, mode=mode, new_edge_length=new_edge_length, rnd=rnd, min_children=2)

    tips_after, freq_after = degree_freq_and_tips(tree)
    print(f"Tips after resolve: {tips_after}")
    print("Degree freq after resolve (deg -> count):", dict(sorted(freq_after.items())))
    warn_if_polytomies_remain(freq_after, idx)

    print("Stripping branch lengths and internal labels..." if strip_internal_labels else "Stripping branch lengths (keeping internal labels)...")
    strip_lengths_and_labels_dendropy(tree, strip_internal_labels=strip_internal_labels)

    tips_final, freq_final = degree_freq_and_tips(tree)
    print(f"Tips final: {tips_final}")
    print("Degree freq final (deg -> count):", dict(sorted(freq_final.items())))

    s = tree.as_string(schema="newick", suppress_rooting=True).strip()
    if not s.endswith(";"):
        s = s + ";"

    print(f"--- Tree #{idx} processed ---")
    del tree

    return s, ""


def _worker_clean(args: Tuple[str, int, str, str, Optional[float], bool, Optional[int]]) -> Tuple[int, Optional[str], str]:
    tree_str, idx, path, mode, new_edge_length, strip_internal_labels, base_seed = args
    # Capture logs in this worker and return them as a string
    from io import StringIO
    old_out, old_err = sys.stdout, sys.stderr
    buf = StringIO()
    sys.stdout = buf
    sys.stderr = buf
    try:
        cleaned, _ = process_single_tree_from_string(
            tree_str=tree_str,
            idx=idx,
            path=path,
            mode=mode,
            new_edge_length=new_edge_length,
            strip_internal_labels=strip_internal_labels,
            base_seed=base_seed,
        )
    finally:
        sys.stdout = old_out
        sys.stderr = old_err
    return idx, cleaned, buf.getvalue()


def main():
    p = argparse.ArgumentParser(description="Resolve polytomies and clean Newick using DendroPy (no rooting).")
    p.add_argument("-i", "--input", nargs="+", required=True, help="Input Newick file(s). Each file may contain multiple trees separated by semicolons.")
    p.add_argument("-o", "--output", required=True, help="Output cleaned Newick file (multiple trees concatenated).")
    p.add_argument("--new-edge-length", type=float, default=None, help="Edge length for newly created internal edges (default: None)")
    p.add_argument("--no-strip-internal-labels", dest="strip_internal_labels", action="store_false", help="Do NOT strip internal node labels (by default internal labels are removed).")
    p.add_argument("--deterministic", dest="deterministic", action="store_true", help="Resolve polytomies deterministically (first-two children) instead of randomly.")
    p.add_argument("--random-seed", type=int, default=None, help="If provided and resolution is random, use this seed for reproducibility.")
    p.add_argument("--num-workers", type=int, default=None, help="Number of parallel worker processes (default: cpu count). Use 1 for serial.")
    args = p.parse_args()

    args.mode = "firstpair" if args.deterministic else "random"

    if args.mode == "random":
        if args.random_seed is not None:
            print(f"Random mode with seed {args.random_seed}")
        else:
            print("Random mode with no fixed seed (non-deterministic)")

    num_workers = args.num_workers if args.num_workers is not None else (os.cpu_count() or 1)
    if num_workers < 1:
        num_workers = 1
    print(f"Using {num_workers} worker(s)")

    total = 0
    with open(args.output, "w", encoding="utf-8") as out_fh:
        for inp in args.input:
            print(f"\nReading file: {inp}")
            try:
                if num_workers == 1:
                    for tree_idx, tree_str in iter_tree_strings_from_file_robust(inp):
                        cleaned_str, _ = process_single_tree_from_string(
                            tree_str=tree_str,
                            idx=tree_idx,
                            path=inp,
                            mode=args.mode,
                            new_edge_length=args.new_edge_length,
                            strip_internal_labels=args.strip_internal_labels,
                            base_seed=args.random_seed,
                        )
                        if cleaned_str is None:
                            continue
                        out_fh.write(cleaned_str + "\n")
                        out_fh.flush()
                        total += 1

                        if total % 1000 == 0:
                            print(f"\n[Progress] Processed {total} trees so far...")
                else:
                    work_iter = (
                        (tree_str, tree_idx, inp, args.mode, args.new_edge_length, args.strip_internal_labels, args.random_seed)
                        for tree_idx, tree_str in iter_tree_strings_from_file_robust(inp)
                    )
                    with ProcessPoolExecutor(max_workers=num_workers) as ex:
                        for tree_idx, cleaned_str, log in ex.map(_worker_clean, work_iter, chunksize=1):
                            if log:
                                if not log.endswith("\n"):
                                    log = log + "\n"
                                sys.stdout.write(log)
                                sys.stdout.flush()
                            if cleaned_str is None:
                                continue
                            out_fh.write(cleaned_str + "\n")
                            out_fh.flush()
                            total += 1

                            if total % 1000 == 0:
                                print(f"\n[Progress] Processed {total} trees so far...")

            except Exception as e:
                print(f"Error reading '{inp}': {e}", file=sys.stderr)
                continue

    if total == 0:
        print("No trees processed. Exiting.", file=sys.stderr)
        sys.exit(1)

    print(f"\nSuccessfully wrote {total} cleaned tree(s) to: {args.output}")
    print("Done.")


if __name__ == "__main__":
    main()
