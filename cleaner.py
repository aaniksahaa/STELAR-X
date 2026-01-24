#!/usr/bin/env python3
"""
cleaner.py

Robustly read one or more Newick files (each file may contain multiple trees separated by semicolons),
resolve polytomies (random or deterministic), strip branch lengths and internal labels (by default),
and write cleaned Newick (multiple trees concatenated).

This version uses a conservative per-tree parser: it splits input on ';' and parses each tree string
separately with dendropy.Tree.get(..., data=...).

python cleaner.py avian-gt-363.tre avian-gt-363-cleaned.tre --outgroup "STRCAM"
python cleaner.py 1kp.tre 1kp-og.tre --outgroup "BWVJ,CKXF,Cyame_v1.0,FTRP,IEHF,IHJY,IKIZ,IKWM,JEBK,JJZR,LJPN,LLXJ,OBUY,PVGP,PWKQ,PYDB,RSOF,RTLC,SBLT,UGPM,URSB,UYFR,VNAL,VZWX,WEJN,XAXW,YSBD,ZJOJ,ZULJ"
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

def find_mrca_of_taxa(tree: dendropy.Tree, taxon_names: List[str]):
    """
    Find the Most Recent Common Ancestor (MRCA) of multiple taxa.
    
    Args:
        tree: The dendropy.Tree object
        taxon_names: List of taxon names to find MRCA for
        
    Returns:
        The MRCA node if found, None otherwise
    """
    # Find all the leaf nodes corresponding to the taxon names
    target_nodes = []
    found_taxa = []
    
    for node in tree.leaf_node_iter():
        if node.taxon and node.taxon.label in taxon_names:
            target_nodes.append(node)
            found_taxa.append(node.taxon.label)
    
    if len(target_nodes) == 0:
        print(f"Warning: None of the outgroup taxa {taxon_names} found in tree.", file=sys.stderr)
        return None
    
    if len(target_nodes) < len(taxon_names):
        missing_taxa = set(taxon_names) - set(found_taxa)
        print(f"Warning: Some outgroup taxa not found in tree: {missing_taxa}", file=sys.stderr)
        print(f"Proceeding with found taxa: {found_taxa}", file=sys.stderr)
    
    if len(target_nodes) == 1:
        # Only one taxon found, return it directly
        print(f"Only one outgroup taxon found: {found_taxa[0]}", file=sys.stderr)
        return target_nodes[0]
    
    try:
        # Use dendropy's mrca method with correct parameter name
        # Create a list of taxa objects from the found nodes
        target_taxa = [node.taxon for node in target_nodes if node.taxon]
        
        if len(target_taxa) == 0:
            print(f"Warning: No valid taxa objects found for MRCA calculation", file=sys.stderr)
            return None
        
        # Use the correct parameter name for dendropy's mrca method
        mrca = tree.mrca(taxa=target_taxa)
        
        if mrca:
            print(f"Successfully found MRCA for {len(found_taxa)} outgroup taxa", file=sys.stderr)
        
        return mrca
        
    except Exception as e:
        print(f"Warning: Failed to find MRCA of taxa {found_taxa}: {e}", file=sys.stderr)
        
        # Fallback: try to find MRCA using a different approach
        try:
            # Alternative approach: find common ancestor by traversing up from first node
            if len(target_nodes) >= 2:
                # Start with the first two nodes and find their common ancestor
                current_mrca = target_nodes[0]
                
                for i in range(1, len(target_nodes)):
                    # Find path from current_mrca to root
                    path1 = []
                    node = current_mrca
                    while node:
                        path1.append(node)
                        node = node.parent_node
                    
                    # Find path from target_nodes[i] to root
                    path2 = []
                    node = target_nodes[i]
                    while node:
                        path2.append(node)
                        node = node.parent_node
                    
                    # Find the first common node in the paths (MRCA)
                    path1_set = set(path1)
                    common_ancestor = None
                    for node in path2:
                        if node in path1_set:
                            common_ancestor = node
                            break
                    
                    if common_ancestor:
                        current_mrca = common_ancestor
                    else:
                        print(f"Warning: Could not find common ancestor", file=sys.stderr)
                        return None
                
                print(f"Successfully found MRCA using fallback method for {len(found_taxa)} taxa", file=sys.stderr)
                return current_mrca
            
        except Exception as e2:
            print(f"Warning: Fallback MRCA method also failed: {e2}", file=sys.stderr)
        
        return None

def root_tree_with_outgroup(tree: dendropy.Tree, outgroup_spec: str):
    """
    Root the tree using the specified outgroup taxon(s).
    
    Args:
        tree: The dendropy.Tree object to root
        outgroup_spec: Either a single taxon name or comma-separated list of taxon names
        
    Returns:
        True if rooting was successful, False otherwise
    """
    # Parse outgroup specification (could be single taxon or comma-separated list)
    outgroup_names = [name.strip() for name in outgroup_spec.split(',') if name.strip()]
    
    if len(outgroup_names) == 0:
        print(f"Warning: No valid outgroup taxa specified in '{outgroup_spec}'", file=sys.stderr)
        return False
    
    if len(outgroup_names) == 1:
        # Single outgroup - use original logic
        outgroup_name = outgroup_names[0]
        outgroup_node = None
        for node in tree.leaf_node_iter():
            if node.taxon and node.taxon.label == outgroup_name:
                outgroup_node = node
                break
        
        if outgroup_node is None:
            print(f"Warning: Outgroup taxon '{outgroup_name}' not found in tree. Skipping rooting.", file=sys.stderr)
            return False
        
        try:
            tree.to_outgroup_position(outgroup_node)
            return True
        except Exception as e:
            print(f"Warning: Failed to root tree with outgroup '{outgroup_name}': {e}", file=sys.stderr)
            return False
    
    else:
        # Multiple outgroups - find MRCA and root there
        print(f"Multiple outgroups specified: {outgroup_names}")
        mrca_node = find_mrca_of_taxa(tree, outgroup_names)
        
        if mrca_node is None:
            print(f"Warning: Could not find MRCA of outgroup taxa. Skipping rooting.", file=sys.stderr)
            return False
        
        try:
            # Root the tree so that the MRCA of outgroups is positioned as outgroup
            if mrca_node.is_leaf():
                # If MRCA is a leaf (only one outgroup taxon found), use single outgroup logic
                tree.to_outgroup_position(mrca_node)
            else:
                # If MRCA is internal, root at the edge leading to this node
                if mrca_node.edge:
                    tree.reroot_at_edge(mrca_node.edge, length1=None, length2=None, update_bipartitions=False)
                else:
                    # Fallback: try to position one of the outgroup taxa as outgroup
                    for node in tree.leaf_node_iter():
                        if node.taxon and node.taxon.label in outgroup_names:
                            tree.to_outgroup_position(node)
                            break
            return True
            
        except Exception as e:
            print(f"Warning: Failed to root tree with multiple outgroups {outgroup_names}: {e}", file=sys.stderr)
            return False

def root_tree_randomly(tree: dendropy.Tree, rnd: Optional[random.Random] = None):
    """
    Root the tree at a random location.
    
    Args:
        tree: The dendropy.Tree object to root
        rnd: Optional random number generator for reproducible results
        
    Returns:
        True if rooting was successful, False otherwise
    """
    try:
        # Get all internal edges (excluding leaf edges)
        internal_edges = [edge for edge in tree.preorder_edge_iter() if edge.head_node and not edge.head_node.is_leaf()]
        
        if not internal_edges:
            # If no internal edges, try to root at midpoint
            tree.reroot_at_midpoint()
            return True
        
        # Select a random internal edge to root at
        if rnd is not None:
            selected_edge = rnd.choice(internal_edges)
        else:
            selected_edge = random.choice(internal_edges)
        
        # Root the tree at the selected edge
        tree.reroot_at_edge(selected_edge, length1=None, length2=None, update_bipartitions=False)
        return True
        
    except Exception as e:
        print(f"Warning: Failed to randomly root tree: {e}", file=sys.stderr)
        # Fallback to midpoint rooting
        try:
            tree.reroot_at_midpoint()
            return True
        except Exception as e2:
            print(f"Warning: Failed to root tree at midpoint: {e2}", file=sys.stderr)
            return False

def process_trees(trees: List[dendropy.Tree], args, rnd: Optional[random.Random] = None):
    cleaned_strings = []
    for idx, tree in enumerate(trees, start=1):
        print(f"\n--- Tree #{idx} (starting) ---")
        tips_before, freq_before = degree_freq_and_tips(tree)
        print(f"Tips before: {tips_before}")
        print("Degree freq before (deg -> count):", dict(sorted(freq_before.items())))

        # Root the tree with outgroup if specified, otherwise root randomly
        rooting_success = False
        if args.outgroup:
            outgroup_names = [name.strip() for name in args.outgroup.split(',') if name.strip()]
            if len(outgroup_names) == 1:
                print(f"Rooting tree with single outgroup: {args.outgroup}")
            else:
                print(f"Rooting tree with multiple outgroups: {outgroup_names}")
                print(f"Will find MRCA of these taxa and root accordingly")
            
            rooting_success = root_tree_with_outgroup(tree, args.outgroup)
            if rooting_success:
                if len(outgroup_names) == 1:
                    print("Tree successfully rooted with single outgroup")
                else:
                    print("Tree successfully rooted using MRCA of multiple outgroups")
            else:
                print("Failed to root tree with outgroup(s), falling back to random rooting")
        
        # If no outgroup specified or outgroup rooting failed, root randomly
        if not rooting_success:
            print("Rooting tree randomly...")
            random_rooting_success = root_tree_randomly(tree, rnd=rnd)
            if random_rooting_success:
                print("Tree successfully rooted randomly")
            else:
                print("Failed to root tree randomly, proceeding with original tree structure")

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
    p.add_argument("--outgroup", type=str, default=None, help="Outgroup taxon name(s) for rooting the tree. Can be a single taxon name or comma-separated list of multiple taxa. For multiple taxa, the tree will be rooted at the MRCA of all specified outgroups.")
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
