#!/usr/bin/env python3
"""
Polytomy simulator (polytomy-design.md §8.5).

Reads a binary gene-tree file and, for each tree, randomly collapses internal
edges — contracting an internal (non-root) node into its parent so the parent
adopts the node's children directly, creating a polytomy.  Writes the resulting
(polytomous) Newick, one per line.

Usage:
  simulate_polytomous.py IN.tre OUT.tre [--prob P] [--seed S]

--prob is the per-internal-edge collapse probability (default 0.3).
The output is guaranteed parseable by STELAR-X (and verify_weights.py).
"""
import sys, re, random, argparse

def _top_commas(s):
    depth, pos = 0, []
    for i, c in enumerate(s):
        if   c == '(': depth += 1
        elif c == ')': depth -= 1
        elif c == ',' and depth == 1: pos.append(i)
    return pos

def parse(s):
    """Newick → nested structure: leaf = str; internal = list[children]."""
    s = s.strip().rstrip(';').strip()
    s = re.sub(r':[^,)]*', '', s)
    if not s.startswith('('):
        return s.strip()
    commas = _top_commas(s)
    childstrs, prev = [], 1
    for c in commas:
        childstrs.append(s[prev:c]); prev = c + 1
    childstrs.append(s[prev:-1])
    return [parse(cs) for cs in childstrs]

def collapse(node, prob, rng, is_root=True):
    """Recursively collapse internal child edges with probability `prob`."""
    if isinstance(node, str):
        return node
    # first recurse into children
    node = [collapse(c, prob, rng, False) for c in node]
    # then splice internal children into this node (collapse the edge above them)
    new_children = []
    for c in node:
        if (not isinstance(c, str)) and len(c) >= 2 and rng.random() < prob:
            new_children.extend(c)        # contract: adopt grandchildren
        else:
            new_children.append(c)
    return new_children

def serialize(node):
    if isinstance(node, str):
        return node
    return '(' + ','.join(serialize(c) for c in node) + ')'

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('inp'); ap.add_argument('outp')
    ap.add_argument('--prob', type=float, default=0.3)
    ap.add_argument('--seed', type=int, default=1)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    n_in = n_poly = 0
    with open(a.inp) as f, open(a.outp, 'w') as g:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'): continue
            n_in += 1
            t = parse(line)
            t = collapse(t, a.prob, rng)
            # count polytomies (internal nodes with >2 children, root too)
            def maxdeg(nd):
                if isinstance(nd, str): return 0
                d = len(nd);
                return max([d] + [maxdeg(c) for c in nd])
            if maxdeg(t) > 2: n_poly += 1
            g.write(serialize(t) + ';\n')
    print(f"simulate_polytomous: {n_in} trees in, {n_poly} now have a polytomy "
          f"(prob={a.prob}, seed={a.seed}) → {a.outp}")

if __name__ == '__main__':
    main()
