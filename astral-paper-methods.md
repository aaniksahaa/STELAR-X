According to documents from **2014 (ASTRAL-I)**, **2015 (ASTRAL-II)**, and **2018 (ASTRAL-III)**, here is a **methods-only** description of what ASTRAL does and how each version differs.   

---

## Core ASTRAL idea (all versions)

### Inputs

* A set **G** of **k** input **gene trees** (typically unrooted; later versions allow partially resolved / polytomous gene trees). 

### Optimization objective: Maximum quartet agreement

* ASTRAL seeks the **species tree** (t) that **maximizes the total number of induced quartet topologies shared with the gene trees**:
  [
  \max_t \sum_{g\in G} |Q(g)\cap Q(t)|
  ]
  where (Q(\cdot)) is the set of induced quartets. This unconstrained problem is NP-hard, so ASTRAL solves a **constrained** version. 

### Constrained search space via a bipartition set (X)

* ASTRAL restricts the output tree’s **bipartitions** (splits) to come from a user-provided/constructed set (X) (with closure: if (A\in X), then (L\setminus A\in X)). The DP only considers tripartitions derivable from (X). 

---

## ASTRAL dynamic programming (ASTRAL-I/II/III share this structure)

### Tripartitions as DP “atoms”

* An unrooted internal node corresponds to a **tripartition** (T=(A|B|C)). ASTRAL assigns a **weight** (w(T)) measuring how many gene-tree quartets “support” that tripartition, then uses DP to assemble a maximum-scoring tree. 

### DP recursion

* Let (L) be the full taxon set, and let (V(A)) be the best score for a subtree spanning taxa (A). ASTRAL recursively splits (A) into (A') and (A\setminus A') (only when these clusters are permitted by (X)) and adds the tripartition score that connects them to the “outside” (L\setminus A):
  [
  V(A)=\max_{A'\subset A,\ (A'|A-A'|L-A)\in Y}\Big(V(A')+V(A-A')+w(A'|A-A'|L-A)\Big)
  ]
  (equivalent forms appear across versions). 

---

## How tripartition weights (w(T)) are defined (quartet-intersection scoring)

### Scoring a tripartition against a gene-tree internal node/partition

* Each gene-tree internal node defines a **partition** (M=(M_1|\dots|M_d)) (with (d=3) for binary nodes and (d>3) for polytomies). For a candidate species-tree tripartition (T=(A|B|C)), the contribution of that gene-tree node is based on
  [
  QI(T,M)
  ]
  which (by design) computes **twice** the number of quartet topologies shared between a tree containing only (T) and a tree containing only (M). 

* The total tripartition weight is then
  [
  w(T)=\sum_{g\in G}\sum_{M\in N(g)} \frac{1}{2}QI(T,M)
  ]
  where (N(g)) is the set of internal nodes/partitions of gene tree (g). 

---

## ASTRAL-I (2014): “ASTRAL: genome-scale coalescent-based species tree estimation”

### What it adds/establishes

1. **Constrained MQSST formulation + DP** to avoid enumerating quartets explicitly; uses (X) (often all bipartitions seen in gene trees). 
2. **Preprocessing**: enumerate allowed tripartitions induced by (X), compute (W(N)) (tripartition weights) by comparing each allowed tripartition (N) to each gene-tree tripartition (M) using a quartet-intersection count (QI(N,M)), then run DP/backtrack. 

---

## ASTRAL-II (2015): “ASTRAL-II: … with many hundreds of taxa and thousands of genes”

ASTRAL-II keeps the same objective + DP, but introduces three main method changes. 

### (1) Faster weight computation: from (O(n^2k)) to (O(nk)) per tripartition

* ASTRAL-I computed intersections using bitsets costing (O(n)) per gene-tree node, giving (O(n^2k)) per tripartition in the worst case.
* ASTRAL-II computes, for each gene tree, the needed intersection counts by a single **postorder traversal** using a stack, accumulating ((|X\cap \text{subtree}|,|Y\cap \text{subtree}|,|Z\cap \text{subtree}|)) and updating (w) along the way (Algorithm 1 “WEIGHT”). This yields (O(nk)) time per tripartition. 

### (2) Expand the constraint set (X) (larger search space) using heuristics

ASTRAL-II adds bipartitions to (X) beyond those directly present in gene trees, using:

* A **quartet-based similarity matrix** between taxa, computed by traversing gene trees and accumulating similarity contributions; then build a **UPGMA** tree and add its bipartitions to (X). 
* **Greedy consensus** trees at multiple thresholds; for each polytomy in a consensus, resolve it multiple ways (UPGMA-based and random-sampling/greedy-based resolutions) and add implied bipartitions to (X). (Algorithm 3 in ASTRAL-II describes this workflow.) 

### (3) Allow polytomous input gene trees

* ASTRAL-II generalizes scoring so that a gene-tree node of degree (d) is treated as representing many tripartitions; the quartet-intersection computation (QI(T,M)) is generalized accordingly (Equation (3) in the ASTRAL-II text). 

---

## ASTRAL-III (2018): “ASTRAL-III: polynomial time species tree reconstruction from partially resolved gene trees”

ASTRAL-III keeps the same DP framework, but redesigns weight computation and (X)-construction to guarantee polynomial scaling (and to handle polytomies efficiently).

### (1) Re-states the problem + DP formally (same structure)

* Uses the same recursion on clusters (A), and the same concept that tripartition scores (w(\cdot)) drive the DP. 

### (2) Makes (X) provably bounded: (|X| = O(nk))

* ASTRAL-II’s heuristic expansions can, in principle, make (|X|) grow too large (even exponentially). ASTRAL-III modifies the heuristics so (|X|) grows at most linearly in (n) and (k). 

### (3) Scores polytomies much faster: from (\Theta(d^3)) to (\Theta(d))

* ASTRAL-II used a cubic-in-degree scoring for a degree-(d) polytomy.
* ASTRAL-III derives an equivalent formulation that computes (QI(T,M)) in time (\Theta(d)) using precomputed sums of (a_i,b_i,c_i) intersection counts (Lemma 1). 

### (4) Reuses repeated gene-tree nodes via a “polytree” of unique clusters

* Instead of scoring each tripartition against each node in each gene tree, ASTRAL-III collects **unique partitions** (E) with multiplicities and defines
  [
  w(T)=\sum_{(M,c)\in E} c\cdot QI(T,M)
  ]
  reducing redundant computations. 

* It builds a DAG-like “polytree” of unique clusters and computes intersections (|W\cap A|,|W\cap B|,|W\cap C|) by memoized postorder traversal, enabling (w(T)) computation in (\Theta(D)) where (D) is the sum of degrees of unique partitions. 

### (5) Overall running time guarantee

* Using a bound (|Y|=O(|X|^{1.726})) plus the (\Theta(D)) scoring, ASTRAL-III achieves:
  [
  O\big(D(nk)^{1.726}\big)
  ]
  for both binary and multifurcating gene trees. 

### (6) DP pruning (“trimming”) to skip provably-unneeded subproblems

* ASTRAL-III introduces an upper bound (U(A)) on (V(A)) and uses it to avoid expanding splits that cannot beat the current best (plus a two-stage “(\alpha)-trimming” refinement). 

### (7) Output add-ons noted in ASTRAL-III (not the DP itself)

* ASTRAL (as a system) can compute **branch lengths in coalescent units** and a **branch support measure** called **local posterior probability** (details are elsewhere, but the capability is explicitly noted). 

---

If you want, I can rewrite the above into a **paper-style “Methods” section** (continuous prose, no bullets) with the same content and citations.


**ASTRAL-MP** (“Parallel ASTRAL”) is *not* a new optimization criterion—it still solves the same **maximum quartet agreement** species-tree problem as ASTRAL-III using the same **DP over constrained clusters** and the same **tripartition-weighting framework**. What changes is **how the expensive steps are implemented** so ASTRAL can scale to much larger (n) (taxa) and/or (k) (gene trees). 

## What ASTRAL-MP does (methods)

### 1) Same core objective + DP as ASTRAL-III

* Input: (k) unrooted gene trees on leaf set (S), (|S|=n). 
* Defines **weighted quartet score** for a candidate species tree (T) as the number of gene-tree–induced quartets that match (T) (extended to missing data / polytomies by counting only fully-resolved quartets). 
* Uses the same **tripartition weight** (w(P)) based on summed quartet intersections (QI(P,M)) over gene-tree internal nodes/partitions, and the same DP recursion constrained by a cluster set (X) (and induced valid split-pairs (Y)). 

So: **same answer in principle**, but **much faster** in practice.

---

### 2) New randomized algorithm for *cluster partitioning* (building (Y) from (X))

In the DP, for each cluster (A\in X), ASTRAL must find all ways to split it into (A') and (A\setminus A') such that **both pieces are in (X)** (i.e., find valid pairs in (Y)). ASTRAL-III did this with costly pairwise intersection checks among subsets of (A). 

ASTRAL-MP speeds this up using a **randomized hashing / Abelian-group homomorphism trick**:

* Represent clusters as 0/1 vectors (a,b,c\in{0,1}^n).
* Observe: (B|C) is a partition of (A) iff (b+c=a) (component-wise).
* Map vectors into a finite Abelian group (G) via a homomorphism (\varphi) so that
  [
  B|C \text{ partitions } A ;\Rightarrow; \varphi(c)=\varphi(a)-\varphi(b).
  ]
* Then to test whether (A\setminus B\in X), you just check whether (\varphi(a)-\varphi(b)) hits a hash-table of ({\varphi(x):x\in X}), avoiding the expensive pairwise set operations. 

They analyze collision probability as *astronomically small* under suitable group size, and explicitly **check for collisions at the end**; if one is detected, rerun with fresh randomness (they report it never happened in their tests). 

---

### 3) Parallelization strategy (CPU threads + SIMD + GPU)

ASTRAL-MP profiles ASTRAL-III into major steps (similarity matrix, greedy consensus, polytomy processing, cluster partitions, tripartition weights, branch scoring, etc.) and then parallelizes/optimizes most of them. 

**(a) Tripartition weight calculation (the bottleneck)**

* Weight calculation (w(P)) is the dominant runtime component, so ASTRAL-MP attacks it with:

  * **CPU multithreading** (compute many (w(P)) values in parallel),
  * **AVX2 vectorization** (compute multiple integer operations per instruction; implemented in C and invoked from Java), 
  * **GPU parallelization**: an OpenCL kernel computes batches of (w(P)) values (data-parallel “many tripartitions at once”), connected from Java via JOCL. 
  * On GPU they choose a compact **array representation** of gene trees rather than the polytree to fit memory constraints. 

**(b) Similarity matrix**

* Parallelized by splitting the set of gene trees across threads; each thread traverses its subset and accumulates contributions. 

**(c) Cluster partitioning**

* Parallelized so different threads handle partition candidates (e.g., by cardinality classes) for each cluster (A). 

**(d) Post-tree branch scoring**

* Branch support (localPP) and branch lengths can be computed independently across branches, so they are parallelizable after the tree is built. 

---

## One-sentence summary

**ASTRAL-MP = ASTRAL-III’s same quartet-score DP, plus (i) a randomized hash-based cluster-partitioning algorithm to speed building (Y), and (ii) heavy parallelization/acceleration (CPU threads, AVX2 SIMD, and OpenCL GPU) especially for tripartition weight computation.** 


