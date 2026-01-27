#!/usr/bin/env python3
"""
root_by_outgroups.py

Root unrooted (seed-rooted) trees using a clustered outgroup, with fallbacks:
1) If outgroup taxa form a true clade, root at the edge above that clade.
2) Else, if ingroup taxa form a true clade, root at the edge above that clade.
3) Else, if a single fallback outgroup taxon (-ogs) is provided and found, root by it.
4) Else, root randomly.

Uses left-to-right leaf ordering and subtree-interval validation for fast clade checks.
Supports multiple trees per input file.
"""
import sys
sys.setrecursionlimit(10_000_000)

import argparse
import os
import random
from concurrent.futures import ProcessPoolExecutor
from typing import Dict, Iterable, List, Optional, Set, Tuple

import dendropy
try:
    from tqdm import tqdm
except Exception:  # pragma: no cover - optional dependency
    tqdm = None


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


def count_trees_in_file(path: str) -> int:
    txt = open(path, "r", encoding="utf-8").read()
    chunks = [c.strip() for c in txt.split(';')]
    return sum(1 for c in chunks if c)


def get_leaf_order(tree: dendropy.Tree) -> List[str]:
    leaves = []
    for leaf in tree.leaf_node_iter():
        if leaf.taxon:
            leaves.append(leaf.taxon.label)
    return leaves


def compute_subtree_ranges_with_nodes(tree: dendropy.Tree, leaf_order: List[str]) -> Tuple[Set[Tuple[int, int]], Dict[Tuple[int, int], dendropy.Node], Dict[str, int]]:
    label_to_idx = {label: i for i, label in enumerate(leaf_order)}
    subtree_ranges_set: Set[Tuple[int, int]] = set()
    range_to_node: Dict[Tuple[int, int], dendropy.Node] = {}

    def dfs(node: dendropy.Node) -> Tuple[int, int]:
        if node.is_leaf():
            if node.taxon and node.taxon.label in label_to_idx:
                idx = label_to_idx[node.taxon.label]
                rng = (idx, idx)
                subtree_ranges_set.add(rng)
                if rng not in range_to_node:
                    range_to_node[rng] = node
                return rng
            return (float('inf'), float('-inf'))

        left_min = float('inf')
        right_max = float('-inf')
        for child in node.child_node_iter():
            child_left, child_right = dfs(child)
            left_min = min(left_min, child_left)
            right_max = max(right_max, child_right)

        if left_min != float('inf') and right_max != float('-inf'):
            rng = (left_min, right_max)
            subtree_ranges_set.add(rng)
            if rng not in range_to_node:
                range_to_node[rng] = node
        return (left_min, right_max)

    if tree.seed_node:
        dfs(tree.seed_node)

    if leaf_order:
        rng = (0, len(leaf_order) - 1)
        subtree_ranges_set.add(rng)
        if rng not in range_to_node and tree.seed_node:
            range_to_node[rng] = tree.seed_node

    return subtree_ranges_set, range_to_node, label_to_idx


def get_consecutive_span(indices: List[int]) -> Tuple[Optional[int], Optional[int], bool]:
    if not indices:
        return None, None, False
    left = min(indices)
    right = max(indices)
    return left, right, (right - left + 1 == len(indices))


def find_group_node(
    group_taxa: Set[str],
    all_tips: Set[str],
    leaf_order: List[str],
    label_to_idx: Dict[str, int],
    subtree_ranges_set: Set[Tuple[int, int]],
    range_to_node: Dict[Tuple[int, int], dendropy.Node],
) -> Dict[str, object]:
    found_taxa = group_taxa & all_tips
    positions = [label_to_idx[t] for t in found_taxa if t in label_to_idx]
    left, right, consecutive = get_consecutive_span(positions)

    result = {
        "found_taxa": found_taxa,
        "positions": positions,
        "left": left,
        "right": right,
        "consecutive": consecutive,
        "is_subtree": False,
        "node": None,
    }

    if not found_taxa or left is None or right is None or not consecutive:
        return result

    if (left, right) in subtree_ranges_set:
        result["is_subtree"] = True
        result["node"] = range_to_node.get((left, right))
    return result


def root_at_node_edge(tree: dendropy.Tree, node: dendropy.Node) -> bool:
    if node is None or node.parent_node is None:
        return False
    try:
        edge = node.edge
        if edge is None:
            return False
        tree.reroot_at_edge(edge, length1=None, length2=None, update_bipartitions=False)
        return True
    except Exception as e:
        print(f"Warning: failed to reroot at edge: {e}", file=sys.stderr)
        return False


def root_with_single_taxon(tree: dendropy.Tree, taxon_label: str) -> bool:
    for node in tree.leaf_node_iter():
        if node.taxon and node.taxon.label == taxon_label:
            try:
                if node.edge is None:
                    return False
                tree.reroot_at_edge(node.edge, length1=None, length2=None, update_bipartitions=False)
                return True
            except Exception as e:
                print(f"Warning: failed to root with single outgroup '{taxon_label}': {e}", file=sys.stderr)
                return False
    return False


def ensure_seed_node(tree: dendropy.Tree):
    if tree.seed_node is None:
        nodes = list(tree.preorder_node_iter())
        if nodes:
            tree.seed_node = nodes[0]


def root_by_root_side(tree: dendropy.Tree, outgroup_found: Set[str], verbose: bool = True) -> bool:
    seed = tree.seed_node
    if seed is None or not outgroup_found:
        return False
    children = list(seed.child_node_iter())
    if len(children) != 2:
        return False

    counts = []
    for child in children:
        count = 0
        for leaf in child.leaf_node_iter():
            if leaf.taxon and leaf.taxon.label in outgroup_found:
                count += 1
        counts.append(count)

    total = len(outgroup_found)
    if counts[0] == total and counts[1] == 0:
        target = children[0]
    elif counts[1] == total and counts[0] == 0:
        target = children[1]
    else:
        return False

    if verbose:
        print(f"{Colors.GREEN}✓ Outgroup taxa all fall on one side of an existing root split.{Colors.RESET}")
        print("Rooting at the edge above that root side...")
    return root_at_node_edge(tree, target)


def root_tree_randomly(tree: dendropy.Tree, rnd: Optional[random.Random] = None) -> bool:
    try:
        internal_edges = [
            edge
            for edge in tree.preorder_edge_iter()
            if edge.head_node
            and edge.tail_node
            and not edge.head_node.is_leaf()
            and not edge.tail_node.is_leaf()
        ]
        candidate_edges = internal_edges
        if not candidate_edges:
            candidate_edges = [edge for edge in tree.preorder_edge_iter() if edge.head_node and edge.tail_node]

        if not candidate_edges:
            return False

        if rnd is None:
            rnd = random.Random()
        rnd.shuffle(candidate_edges)

        for edge in candidate_edges[:10]:
            try:
                tree.reroot_at_edge(edge, length1=None, length2=None, update_bipartitions=False)
                return True
            except Exception:
                continue

        # Last resort: midpoint (may fail on some malformed trees)
        tree.reroot_at_midpoint()
        return True
    except Exception as e:
        print(f"Warning: failed to randomly root tree: {e}", file=sys.stderr)
        try:
            tree.reroot_at_midpoint()
            return True
        except Exception as e2:
            print(f"Warning: failed to root at midpoint: {e2}", file=sys.stderr)
            return False


def read_taxa_file(path: str) -> Set[str]:
    taxa: Set[str] = set()
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            parts = [p.strip() for p in line.strip().split(',') if p.strip()]
            for p in parts:
                taxa.add(p)
    return taxa


def process_single_tree(
    tree: dendropy.Tree,
    idx: int,
    outgroup_taxa: Set[str],
    fallback_single: Optional[str],
    rnd: Optional[random.Random],
    verbose: bool = True,
) -> Tuple[str, Dict[str, int]]:
    stats = {
        "outgroup_clade": 0,
        "ingroup_clade": 0,
        "rooted_side": 0,
        "fallback_single": 0,
        "random": 0,
        "found_any_outgroup": 0,
        "found_all_outgroup": 0,
        "missing_outgroup": 0,
    }

    leaf_order = get_leaf_order(tree)
    all_tips = set(leaf_order)
    total_tips = len(leaf_order)

    found_outgroup = outgroup_taxa & all_tips
    missing_outgroup = outgroup_taxa - all_tips

    if verbose:
        print(f"\n{Colors.BOLD}=== Tree #{idx} ==={Colors.RESET}")
        print(f"Tips: {Colors.CYAN}{total_tips}{Colors.RESET}")
        print(f"Outgroup specified: {Colors.CYAN}{len(outgroup_taxa)}{Colors.RESET}")
        print(f"Outgroup found: {Colors.CYAN}{len(found_outgroup)}{Colors.RESET}")
        if missing_outgroup:
            print(f"{Colors.YELLOW}Missing outgroup taxa ({len(missing_outgroup)}):{Colors.RESET}")
            if len(missing_outgroup) <= 20:
                for t in sorted(missing_outgroup):
                    print(f"  - {t}")
            else:
                print("  (too many to list)")

    if found_outgroup:
        stats["found_any_outgroup"] = 1
    if outgroup_taxa and not missing_outgroup:
        stats["found_all_outgroup"] = 1
    if missing_outgroup:
        stats["missing_outgroup"] = 1

    action_taken = None
    if found_outgroup:
        ensure_seed_node(tree)
        subtree_ranges_set, range_to_node, label_to_idx = compute_subtree_ranges_with_nodes(tree, leaf_order)

        outgroup_result = find_group_node(found_outgroup, all_tips, leaf_order, label_to_idx, subtree_ranges_set, range_to_node)

        if outgroup_result["consecutive"] and outgroup_result["is_subtree"] and outgroup_result["node"] is not None:
            if verbose:
                left = outgroup_result["left"]
                right = outgroup_result["right"]
                print(f"{Colors.GREEN}✓ Outgroup forms a true clade: interval [{left}, {right}]{Colors.RESET}")
                print("Rooting at edge above outgroup clade...")
            if root_at_node_edge(tree, outgroup_result["node"]):
                stats["outgroup_clade"] = 1
                action_taken = "outgroup_clade"
            else:
                if verbose:
                    print(f"{Colors.YELLOW}⚠ Failed to reroot above outgroup clade. Will try other options.{Colors.RESET}")
        else:
            if verbose:
                if not outgroup_result["consecutive"]:
                    print(f"{Colors.RED}✗ Outgroup taxa are not consecutive in leaf order (not monophyletic).{Colors.RESET}")
                elif not outgroup_result["is_subtree"]:
                    print(f"{Colors.YELLOW}⚠ Outgroup taxa are consecutive but NOT a true subtree.{Colors.RESET}")

        if action_taken is None:
            ingroup = all_tips - found_outgroup
            if ingroup:
                ingroup_result = find_group_node(ingroup, all_tips, leaf_order, label_to_idx, subtree_ranges_set, range_to_node)
                if ingroup_result["consecutive"] and ingroup_result["is_subtree"] and ingroup_result["node"] is not None:
                    if verbose:
                        left = ingroup_result["left"]
                        right = ingroup_result["right"]
                        print(f"{Colors.GREEN}✓ Ingroup forms a true clade: interval [{left}, {right}]{Colors.RESET}")
                        print("Rooting at edge above ingroup clade...")
                    if root_at_node_edge(tree, ingroup_result["node"]):
                        stats["ingroup_clade"] = 1
                        action_taken = "ingroup_clade"
                    else:
                        if verbose:
                            print(f"{Colors.YELLOW}⚠ Failed to reroot above ingroup clade. Will try other options.{Colors.RESET}")
                else:
                    if verbose:
                        print(f"{Colors.YELLOW}⚠ Ingroup does not form a true clade either.{Colors.RESET}")
            else:
                if verbose:
                    print(f"{Colors.YELLOW}⚠ Ingroup is empty (all tips are outgroup); cannot root by split.{Colors.RESET}")

        if action_taken is None:
            if root_by_root_side(tree, found_outgroup, verbose=verbose):
                stats["rooted_side"] = 1
                action_taken = "rooted_side"
            else:
                if verbose:
                    print(f"{Colors.YELLOW}⚠ Outgroup taxa are split across root sides (or tree not rooted by an edge).{Colors.RESET}")

    if action_taken is None and fallback_single:
        ensure_seed_node(tree)
        if verbose:
            print(f"{Colors.CYAN}Trying single fallback outgroup: {fallback_single}{Colors.RESET}")
        if root_with_single_taxon(tree, fallback_single):
            stats["fallback_single"] = 1
            action_taken = "fallback_single"
            if verbose:
                print(f"{Colors.GREEN}✓ Rooted by single fallback outgroup.{Colors.RESET}")
        else:
            if verbose:
                print(f"{Colors.YELLOW}⚠ Single fallback outgroup not found; will root randomly.{Colors.RESET}")

    if action_taken is None:
        ensure_seed_node(tree)
        if verbose:
            print(f"{Colors.CYAN}Rooting randomly (fallback).{Colors.RESET}")
        if root_tree_randomly(tree, rnd=rnd):
            stats["random"] = 1
            action_taken = "random"
            if verbose:
                print(f"{Colors.GREEN}✓ Random rooting complete.{Colors.RESET}")
        else:
            if verbose:
                print(f"{Colors.RED}✗ Random rooting failed; tree left as-is.{Colors.RESET}")

    s = tree.as_string(schema="newick", suppress_rooting=True).strip()
    if not s.endswith(";"):
        s = s + ";"
    return s, stats


def process_single_tree_from_string(
    tree_str: str,
    idx: int,
    path: str,
    outgroup_taxa: Set[str],
    fallback_single: Optional[str],
    base_seed: Optional[int],
    verbose: bool,
) -> Tuple[Optional[str], Dict[str, int]]:
    tree = parse_tree_from_string(tree_str, idx, path)
    if tree is None:
        return None, {
            "outgroup_clade": 0,
            "ingroup_clade": 0,
            "rooted_side": 0,
            "fallback_single": 0,
            "random": 0,
            "found_any_outgroup": 0,
            "found_all_outgroup": 0,
            "missing_outgroup": 0,
        }
    rnd = random.Random(base_seed + idx) if base_seed is not None else None
    return process_single_tree(
        tree=tree,
        idx=idx,
        outgroup_taxa=outgroup_taxa,
        fallback_single=fallback_single,
        rnd=rnd,
        verbose=verbose,
    )


def _worker_root(args: Tuple[str, int, str, Set[str], Optional[str], Optional[int], bool]) -> Tuple[int, Optional[str], Dict[str, int], str]:
    tree_str, idx, path, outgroup_taxa, fallback_single, base_seed, verbose = args
    from io import StringIO
    old_out, old_err = sys.stdout, sys.stderr
    buf = StringIO()
    sys.stdout = buf
    sys.stderr = buf
    try:
        rooted, stats = process_single_tree_from_string(
            tree_str=tree_str,
            idx=idx,
            path=path,
            outgroup_taxa=outgroup_taxa,
            fallback_single=fallback_single,
            base_seed=base_seed,
            verbose=verbose,
        )
    finally:
        sys.stdout = old_out
        sys.stderr = old_err
    return idx, rooted, stats, buf.getvalue()


def main():
    p = argparse.ArgumentParser(description="Root trees using clustered outgroups with graceful fallbacks.")
    p.add_argument("-i", "--input", required=True, help="Input Newick file (may contain multiple trees)")
    p.add_argument("-o", "--output", required=True, help="Output Newick file (rooted trees)")
    p.add_argument("-og", "--outgroup", default=None, help="Comma-separated list of outgroup taxa")
    p.add_argument("-ogf", "--outgroup-file", default=None, help="File containing outgroup taxa (one per line or comma-separated)")
    p.add_argument("-ogs", "--outgroup-single", default=None, help="Single fallback outgroup taxon if outgroup set fails")
    p.add_argument("--random-seed", type=int, default=None, help="Seed for random rooting (fallback)")
    p.add_argument("-q", "--quiet", action="store_true", help="Quiet mode (minimal output)")
    p.add_argument("--num-workers", type=int, default=None, help="Number of parallel worker processes (default: cpu count). Use 1 for serial.")
    args = p.parse_args()

    outgroup_taxa: Set[str] = set()
    if args.outgroup:
        outgroup_taxa.update(t.strip() for t in args.outgroup.split(',') if t.strip())
    if args.outgroup_file:
        try:
            outgroup_taxa.update(read_taxa_file(args.outgroup_file))
        except Exception as e:
            print(f"Error reading outgroup file: {e}", file=sys.stderr)
            sys.exit(1)

    if not outgroup_taxa and not args.outgroup_single and not args.quiet:
        print(f"{Colors.YELLOW}⚠ No outgroup set provided; will use single fallback if given or random rooting.{Colors.RESET}")

    try:
        num_trees = count_trees_in_file(args.input)
    except Exception as e:
        print(f"Error reading tree file: {e}", file=sys.stderr)
        sys.exit(1)

    if num_trees == 0:
        print(f"Error: No valid trees found in {args.input}", file=sys.stderr)
        sys.exit(1)

    if not args.quiet:
        print(f"{Colors.BOLD}Found {num_trees} tree(s) in file{Colors.RESET}")

    if args.random_seed is not None and not args.quiet:
        print(f"{Colors.CYAN}Random seed set to {args.random_seed}{Colors.RESET}")

    num_workers = args.num_workers if args.num_workers is not None else (os.cpu_count() or 1)
    if num_workers < 1:
        num_workers = 1
    if not args.quiet:
        print(f"{Colors.CYAN}Using {num_workers} worker(s){Colors.RESET}")

    totals = {
        "trees": 0,
        "outgroup_clade": 0,
        "ingroup_clade": 0,
        "rooted_side": 0,
        "fallback_single": 0,
        "random": 0,
        "found_any_outgroup": 0,
        "found_all_outgroup": 0,
        "missing_outgroup": 0,
    }

    with open(args.output, "w", encoding="utf-8") as out_fh:
        if num_workers == 1:
            tree_iter = iter_tree_strings_from_file_robust(args.input)
            if tqdm is not None:
                tree_iter = tqdm(tree_iter, total=num_trees, desc="Rooting trees", unit="tree", leave=True)
            elif not args.quiet:
                print(f"{Colors.YELLOW}⚠ tqdm not available; progress bar disabled.{Colors.RESET}")

            for tree_idx, tree_str in tree_iter:
                rooted_str, stats = process_single_tree_from_string(
                    tree_str=tree_str,
                    idx=tree_idx,
                    path=args.input,
                    outgroup_taxa=outgroup_taxa,
                    fallback_single=args.outgroup_single,
                    base_seed=args.random_seed,
                    verbose=not args.quiet,
                )
                if rooted_str is None:
                    continue
                out_fh.write(rooted_str + "\n")
                out_fh.flush()

                totals["trees"] += 1
                for k in stats:
                    totals[k] += stats[k]

                if tqdm is None and totals["trees"] % 1000 == 0 and not args.quiet:
                    print(f"\n[Progress] Processed {totals['trees']} trees so far...")
        else:
            work_iter = (
                (tree_str, tree_idx, args.input, outgroup_taxa, args.outgroup_single, args.random_seed, not args.quiet)
                for tree_idx, tree_str in iter_tree_strings_from_file_robust(args.input)
            )
            progress = None
            if tqdm is not None:
                progress = tqdm(total=num_trees, desc="Rooting trees", unit="tree", leave=True)
            elif not args.quiet:
                print(f"{Colors.YELLOW}⚠ tqdm not available; progress bar disabled.{Colors.RESET}")

            with ProcessPoolExecutor(max_workers=num_workers) as ex:
                for tree_idx, rooted_str, stats, log in ex.map(_worker_root, work_iter, chunksize=1):
                    if log and not args.quiet:
                        if not log.endswith("\n"):
                            log = log + "\n"
                        sys.stdout.write(log)
                        sys.stdout.flush()
                    if rooted_str is None:
                        if progress is not None:
                            progress.update(1)
                        continue
                    out_fh.write(rooted_str + "\n")
                    out_fh.flush()

                    totals["trees"] += 1
                    for k in stats:
                        totals[k] += stats[k]

                    if progress is not None:
                        progress.update(1)
                    elif totals["trees"] % 1000 == 0 and not args.quiet:
                        print(f"\n[Progress] Processed {totals['trees']} trees so far...")

            if progress is not None:
                progress.close()

    if totals["trees"] == 0:
        print("No trees processed. Exiting.", file=sys.stderr)
        sys.exit(1)

    if not args.quiet:
        print(f"\n{Colors.BOLD}{'='*60}{Colors.RESET}")
        print(f"{Colors.BOLD}=== OVERALL SUMMARY ({totals['trees']} trees) ==={Colors.RESET}")
        print(f"Outgroup clade rooted:   {Colors.CYAN}{totals['outgroup_clade']}{Colors.RESET}")
        print(f"Ingroup clade rooted:    {Colors.CYAN}{totals['ingroup_clade']}{Colors.RESET}")
        print(f"Root-side (rooted tree): {Colors.CYAN}{totals['rooted_side']}{Colors.RESET}")
        print(f"Single fallback rooted:  {Colors.CYAN}{totals['fallback_single']}{Colors.RESET}")
        print(f"Random rooted:           {Colors.CYAN}{totals['random']}{Colors.RESET}")
        if outgroup_taxa:
            print(f"Trees with any outgroup found: {Colors.CYAN}{totals['found_any_outgroup']}{Colors.RESET}")
            print(f"Trees with all outgroup found: {Colors.CYAN}{totals['found_all_outgroup']}{Colors.RESET}")
            print(f"Trees missing outgroups:       {Colors.CYAN}{totals['missing_outgroup']}{Colors.RESET}")
        print(f"{Colors.BOLD}Output written to:{Colors.RESET} {args.output}")


if __name__ == "__main__":
    main()
