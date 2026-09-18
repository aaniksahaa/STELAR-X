#!/usr/bin/env python3
"""
Reference STELAR-X weight verifier.

Independently recomputes:
  1. Gene-tree tripartitions (M1|M2|M3) with frequencies.
  2. Candidate clusters (sub(u) and S\\sub(u) for every non-root u).
  3. QI weight for each (candidate_split, tripartition) pair using:
       - lgA = |A ∩ Lg_GT|         (row sum, correct for incomplete trees)
       - c2  = sz3 - a2 - b2       (column M3 constraint, correct formula)
  4. Final split scores and inferred species tree.

Matches STELAR-X exactly.  Any discrepancy indicates a bug in either.
"""

import sys, re
from collections import defaultdict
from itertools import permutations as iperms

# ─── Newick parser ────────────────────────────────────────────────────────────

def _top_commas(s):
    """Positions of ',' at parenthesis depth 1."""
    depth, pos = 0, []
    for i, c in enumerate(s):
        if   c == '(': depth += 1
        elif c == ')': depth -= 1
        elif c == ',' and depth == 1: pos.append(i)
    return pos

# Node representation:  leaf = frozenset([name]);  internal = (children_tuple, own_leaves)
def _node(children):
    lv = frozenset()
    for c in children:
        lv = lv | leaves(c)
    return (tuple(children), lv)

def parse_newick(s, is_root=True):
    s = s.strip().rstrip(';').strip()
    # strip branch lengths / internal labels at this level
    s = re.sub(r':[^,)]*', '', s)   # remove ':0.12' annotations
    if not s.startswith('('):
        return frozenset([s.strip()])
    commas = _top_commas(s)
    # split into child substrings at depth-1 commas
    childstrs, prev = [], 1
    for c in commas:
        childstrs.append(s[prev:c]); prev = c + 1
    childstrs.append(s[prev:-1])
    children = [parse_newick(cs, False) for cs in childstrs]
    nc = len(children)
    # Mirror STELAR-X TreeParser.validateAndConvert (polytomy-design.md §3.2):
    if nc == 2:
        return _node(children)                              # binary internal
    if nc == 3 and is_root:
        inner = _node([children[1], children[2]])           # arbitrary binary rooting
        return _node([children[0], inner])
    if nc >= 4 and is_root:
        inner = _node(children[1:])                         # polytomous inner; complement=child0
        return _node([children[0], inner])
    # nc >= 3 and not root → polytomous node
    return _node(children)

def leaves(node):
    if isinstance(node, frozenset): return node
    return node[1]

def subtrees(node):
    """All internal nodes in the tree, post-order."""
    if isinstance(node, frozenset): return []
    result = []
    for c in node[0]:
        result += subtrees(c)
    result.append(node)
    return result

# ─── d-Partitions ──────────────────────────────────────────────────────────────

def extract_tripartitions(parsed_trees):
    """
    For each non-root internal node u of gene tree g (leaf set Lg) with children
    c0..c_{k-1}:  the d-partition (d = k+1) is
      M_i = sub(c_i)            for i = 0..k-1
      M_{d-1} = Lg - sub(u)     (complement)
    Returns: dict{ key -> frequency }.  For binary (d=3) the key matches STELAR-X's
    PartitionHash exactly (sort {M0,M1}, complement separate); for d≥4 the key is
    order-invariant over all d parts.
    """
    triparts = defaultdict(int)
    for (tree, lg) in parsed_trees:
        for node in subtrees(tree):
            own = leaves(node)
            if own == lg:
                continue   # skip root
            comp = lg - own
            if not comp:
                continue   # shouldn't happen for non-root
            child_parts = [leaves(c) for c in node[0]]
            if len(child_parts) == 2:
                m1, m2 = child_parts
                key = tuple(sorted([m1, m2], key=lambda x: sorted(x))) + (comp,)
            else:
                key = tuple(sorted(child_parts + [comp], key=lambda x: sorted(x)))
            triparts[key] += 1
    return triparts

# ─── Clusters ─────────────────────────────────────────────────────────────────

def extract_clusters(parsed_trees, S):
    """
    For every non-root internal node u in every gene tree:
      register sub(u)   and   S\\sub(u).
    Also: singletons are always valid (leaf DP nodes).
    Returns: dict{ frozenset -> count }
    """
    clusters = defaultdict(int)
    # Singletons are always in the candidate set
    for t in S:
        clusters[frozenset([t])] += 0  # ensure present with count>=0
    for (tree, lg) in parsed_trees:
        nodes = subtrees(tree)
        root_leaves = lg
        for node in nodes:
            own = leaves(node)
            if own == root_leaves: continue
            # sub(u)  (for polytomous nodes this is the whole subtree; children's
            # own sets are registered when subtrees() recurses into them — NO combos)
            clusters[own] += 1
            # super-complement S \ sub(u)
            sc = S - own
            if sc: clusters[sc] += 1
        # also register Lg and S\Lg for incomplete trees (Type 3 DP)
        if lg != S:
            clusters[lg]       += 1
            clusters[S - lg]   += 1
    return clusters

# ─── QI computation ───────────────────────────────────────────────────────────

def two_qi(a, b, c):
    """2*QI = sum over distinct (i,j,k) of a[i]*b[j]*c[k]*(a[i]+b[j]+c[k]-3),
    the brute O(d³) definition for a d-partition (polytomy-design.md §4).
    For d=3 this equals the binary 6-permutation formula."""
    d = len(a)
    res = 0
    for i in range(d):
        for j in range(d):
            if j == i: continue
            for k in range(d):
                if k == i or k == j: continue
                ai, bj, ck = a[i], b[j], c[k]
                s = ai + bj + ck - 3
                if s > 0:
                    res += ai * bj * ck * s
    return res


def score_split(A_set, B_set, S, triparts, verbose=False):
    """
    Score the species-tree split (A_set | B_set), where C = S - A - B.

    For each d-partition (M_0|…|M_{d-1}) with frequency f, compute explicitly:
      a[i] = |A ∩ M_i|,  b[i] = |B ∩ M_i|,  c[i] = |M_i| - a[i] - b[i]
    (A,B disjoint ⇒ c[i] ≥ 0; the complement part M_{d-1} is stored explicitly so
    a[i]/b[i] are the true intersections — equivalent to STELAR-X's row-constraint
    derivation, but valid for any degree d, binary or polytomous).

    Returns: score = (1/2) * sum_P freq * 2*QI
    """
    two_score = 0
    for parts, freq in triparts.items():
        a = [len(A_set & M) for M in parts]
        b = [len(B_set & M) for M in parts]
        c = [len(M) - a[i] - b[i] for i, M in enumerate(parts)]
        if any(x < 0 for x in c):
            continue
        tqi = two_qi(a, b, c)
        if verbose:
            szs = '|'.join(str(len(M)) for M in parts)
            print(f"    PART d={len(parts)} sz={szs} a={a} b={b} c={c} 2*QI={tqi} freq={freq}")
        two_score += freq * tqi

    return two_score // 2

# ─── DP ───────────────────────────────────────────────────────────────────────

def infer(S, clusters, scores_map):
    """
    ASTRAL inference DP.
    dp[cluster] = max over all valid splits B|R of cluster:
                    score(B, R) + dp[B] + dp[R]
    A valid split of cluster P: B ⊂ P, R = P-B, both B and R in clusters ∪ {S}.
    """
    all_cl = set(clusters) | {S}

    # Build adjacency: parent → list of (B, R)
    adj = defaultdict(list)
    for P in all_cl:
        for B in clusters:
            if not (B < P): continue
            R = P - B
            if R not in all_cl: continue
            if len(B) > len(R): continue    # canonical: smaller or lex-first
            elif len(B) == len(R) and sorted(B) > sorted(R): continue
            adj[P].append((B, R))

    memo = {}
    best = {}

    def dp(P):
        if P in memo: return memo[P]
        if len(P) == 1:
            memo[P] = 0; return 0
        best_sc = -1
        for (B, R) in adj.get(P, []):
            key = (min(B,R,key=lambda x:sorted(x)), max(B,R,key=lambda x:sorted(x)))
            sc = scores_map.get(key, 0) + dp(B) + dp(R)
            if sc > best_sc:
                best_sc = sc
                best[P] = (B, R)
        memo[P] = max(best_sc, 0)
        return memo[P]

    total = dp(S)

    def newick(P):
        if len(P) == 1: return next(iter(P))
        if P not in best: return '(' + ','.join(sorted(P)) + ')'
        B, R = best[P]
        return '(' + newick(B) + ',' + newick(R) + ')'

    return total, newick(S) + ';'

# ─── Main ──────────────────────────────────────────────────────────────────────

def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('input', help='Gene-tree file (one Newick per line)')
    ap.add_argument('--dump-splits',  action='store_true')
    ap.add_argument('--dump-triparts', action='store_true')
    ap.add_argument('--dump-clusters', action='store_true')
    ap.add_argument('--verbose-score', action='store_true',
                    help='Print per-(split,partition) QI values')
    ap.add_argument('--compare', metavar='TRUE_NEWICK',
                    help='Compute RF against this tree string')
    args = ap.parse_args()

    # --- parse ---
    raw = [l.strip() for l in open(args.input) if l.strip() and not l.startswith('#')]
    parsed = []
    for line in raw:
        try:
            t = parse_newick(line)
            lg = leaves(t)
            parsed.append((t, lg))
        except Exception as e:
            print(f"[warn] skip: {e}", file=sys.stderr)
    if not parsed: sys.exit("No trees parsed")

    S = frozenset().union(*(lg for _, lg in parsed))
    n = len(S)
    print(f"=== {len(parsed)} gene trees, {n} taxa: {sorted(S)} ===")

    # --- tripartitions ---
    triparts = extract_tripartitions(parsed)
    print(f"[Tripartitions] {len(triparts)} unique")
    if args.dump_triparts:
        for parts, f in sorted(triparts.items(), key=lambda x: (-x[1], sorted(x[0][0]))):
            body = ' | '.join(str(sorted(M)) for M in parts)
            szs  = '|'.join(str(len(M)) for M in parts)
            print(f"  freq={f:3d}  d={len(parts)}  {body}  sz={szs}")

    # --- clusters ---
    clusters = extract_clusters(parsed, S)
    print(f"[Clusters] {len(clusters)} unique")
    if args.dump_clusters:
        for cl, f in sorted(clusters.items(), key=lambda x:(len(x[0]),sorted(x[0]))):
            print(f"  sz={len(cl):2d}  freq={f:3d}  {sorted(cl)}")

    # --- score all candidate splits ---
    all_cl = set(clusters) | {S}
    scores_map = {}
    n_scored = 0
    for P in all_cl:
        for B in clusters:
            if not (B < P): continue
            R = P - B
            if R not in all_cl: continue
            if len(B) > len(R): continue
            if len(B) == len(R) and sorted(B) > sorted(R): continue
            key = (min(B,R,key=lambda x:sorted(x)), max(B,R,key=lambda x:sorted(x)))
            if key in scores_map: continue
            if args.verbose_score:
                print(f"SPLIT sz={len(B)}|{len(R)}  {sorted(B)} | {sorted(R)}")
            sc = score_split(B, R, S, triparts, verbose=args.verbose_score)
            if args.verbose_score:
                print(f"  => score={sc}")
            scores_map[key] = sc
            n_scored += 1

    print(f"[Scores] {n_scored} candidate splits scored")
    if args.dump_splits:
        for (B, R), sc in sorted(scores_map.items(), key=lambda x: -x[1]):
            print(f"  score={sc:8d}  {sorted(B)} | {sorted(R)}")

    # --- DP ---
    total, tree = infer(S, clusters, scores_map)
    print(f"[Inference] quartet score = {total}")
    print(f"Species tree: {tree}")

    # --- RF ---
    if args.compare:
        try:
            import dendropy
            from dendropy.calculate import treecompare
            tns = dendropy.TaxonNamespace()
            t1 = dendropy.Tree.get(data=args.compare, schema='newick', taxon_namespace=tns)
            t2 = dendropy.Tree.get(data=tree, schema='newick', taxon_namespace=tns)
            rf = treecompare.symmetric_difference(t1, t2)
            mx = max(1, 2*(n-3))
            print(f"RF = {rf}  (norm={rf/mx:.4f}, similarity={1-rf/mx:.1%})")
        except ImportError:
            print("(dendropy not available)")

if __name__ == '__main__':
    main()
