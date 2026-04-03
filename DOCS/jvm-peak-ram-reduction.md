# JVM Peak RAM Reduction in Multi-Phase Pipelines

## Motivation

ASTRAL-X (and similar multi-phase phylogenetic inference tools) runs as a single long-lived JVM process whose
pipeline looks roughly like:

```
Parse → Hash → Clusters → Partitions → DP → Weight Calc → Inference
```

Each phase creates large intermediate data structures.  By the time a late phase starts, several of those
structures are **logically dead** — the pipeline will never touch them again — but the JVM may not have
reclaimed their memory yet.  The result is a RAM peak that substantially exceeds what any single phase
actually needs.

Two JVM behaviours make this worse than it sounds:

### 1. Deferred GC under a large `-Xmx`

Modern collectors (G1GC, ZGC, Shenandoah) are *concurrent* and *lazy*: they collect only when heap pressure
crosses internal thresholds.  When you set `-Xmx 512g` on a machine with 64 GB of physical RAM, the JVM sees
virtually unlimited headroom and feels no urgency.  Dead objects from Phase 3 can still occupy memory when
Phase 6 starts its own large allocations, causing both generations of data to coexist — and the process RSS
to spike.

### 2. JNI array pinning

When a Java array is passed to a native (JNI) method, the garbage collector **must not move or reclaim it**
for the duration of that call.  If you later null the reference in Java, the GC still cannot collect it until
the native call returns.  This means that during a GPU kernel call, all input arrays are effectively live
regardless of whether Java code has released its references.  The window where both input and output arrays
coexist in RAM can be several times larger than either alone.

---

## The Two-Part Fix

### Part A — Explicit null at phase boundaries

As soon as a large object will not be used again, explicitly set its reference to `null` in the calling code.
Do this at the **call site** (e.g. `Main.java`), not inside the class that owns the object, so the intent is
visible in the phase-flow narrative.

```java
// Example: hasher is only used to build the prefix arrays (pref).
// Once pref is constructed, hasher's data structures (~hundreds of MB) are dead.
hasher = null;   // allow GC to reclaim before the next large allocation

// Example: pref (PrefixHashArrays, ~3 GB) is only needed through Phase 5.
// Phase 6 (weight calculation) allocates its own large working set.
pref = null;     // free ~3 GB before Phase 6 allocates ~8–12 GB
```

The reference itself costs nothing to null.  The benefit is that the next GC cycle (whenever it runs) can
actually reclaim the object rather than treating it as live.

### Part B — Explicit GC hint at expensive phase transitions

After nulling the dead references, call `System.gc()` to *hint* the JVM to collect now, before the next phase
allocates its large working set.

```java
long gcHeapBefore = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
System.gc();
long gcHeapAfter  = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
Logging.debug("Pre-Phase-6 GC hint: heap %d MB → %d MB (freed %d MB)",
    gcHeapBefore / 1_000_000, gcHeapAfter / 1_000_000,
    (gcHeapBefore - gcHeapAfter) / 1_000_000);
```

`System.gc()` is a hint, not a command.  It is suppressed by `-XX:+DisableExplicitGC` (verify this is not set
on the target machine: `java -XX:+PrintFlagsFinal -version 2>&1 | grep DisableExplicitGC`).  When it works,
it gives the JVM a chance to compact the heap *before* the next phase inflates it, reducing the peak
coexistence window.

> **When to place the hint:** just before a phase that is known to allocate a large, long-lived working set
> — i.e. a phase whose allocations will themselves suppress future GC by filling the heap.

### Part C — Null input arrays inside GPU/native call wrappers

This is a specific case of Part A that applies whenever a method passes large arrays to a native call and
then processes the return value in a loop.

```java
long[] twoScores = GPUWeightCalculator.computeWeightsGPU(
    splitsData, partsData, orderings, invIndex, ...);

// The native call has returned.  splitsData / partsData / orderings / invIndex
// are no longer referenced by any live code path.  Null them now so the GC
// can reclaim ~7 GB while we iterate over twoScores.
splitsData = null;
partsData  = null;
orderings  = null;
invIndex   = null;

for (int i = 0; i < numSplits; i++) {
    scoreArray[i] = twoScores[i] / 2L;
}
// Do NOT null twoScores here — it is still in use.
```

The key safety check: confirm that the nulled variables are **local parameters** or fields that are not read
again after the native call.  The return value (`twoScores`) must remain live until the loop completes.

---

## How to Identify Candidates in a New Codebase

Walk the phase sequence and build a table:

| Object | Created | Last read | Can be nulled after |
|--------|---------|-----------|---------------------|
| `Hasher` / hash seeds | Phase 2 | Phase 2 (to build prefix arrays) | Phase 2 completes |
| Prefix / hash arrays | Phase 2 | Phase 5 (DP) | Phase 5 completes |
| Cluster table | Phase 3 | Phase 5 (DP) | Phase 5 completes |
| DP table | Phase 5 | Phase 5 (inference) | Phase 5 completes |
| GPU input arrays | Phase 6 | Phase 6 (before JNI call returns) | Immediately after JNI returns |

Anything whose "can be nulled after" column is well before the phase that causes the peak RAM spike is a
candidate for the fix.

---

## Observed Effect in ASTRAL-X (n = 100 K taxa, k = 1 000 gene trees)

| Without fix | With fix |
|-------------|----------|
| ~45–48 GB at Phase 6 start | ~21 GB at Phase 6 start |
| Peak overlaps Phase 3–5 garbage + Phase 6 working set | Phase 3–5 garbage collected before Phase 6 inflates |

Three separate null + GC actions together freed ~24–27 GB of peak coexistence:

1. `hasher = null` (Phase 2 → 3 boundary): reclaims hash seed structures
2. `pref = null` + `System.gc()` (Phase 5 → 6 boundary): reclaims ~3.2 GB prefix arrays +
   whatever Phase 3–5 intermediates survived
3. `splitsData/partsData/orderings/invIndex = null` inside `WeightTable.computeScoresGPU()`:
   reclaims ~7–8 GB of GPU input arrays while the output loop runs

---

## Applicability to STELAR-X

STELAR-X has a similar multi-phase structure (parsing → bipartition manager → weight calculation → DP →
tree reconstruction).  The same pattern applies:

- **`MemoryEfficientBipartitionManager`** (prefix sums, XOR arrays, hash-to-bipartition map) is only needed
  through the weight calculation phase.  Nulling it before the DP phase begins frees the largest intermediate
  structures before the DP memo table inflates.

- **`rangeBipWeights`** and **`clusterHashToRangeBips`** are needed through DP but can be released once
  `reconstructTree()` returns.

- If STELAR-X is ever run with a large `-Xmx` on a memory-constrained machine (or if future versions add GPU
  acceleration with JNI input arrays), the `System.gc()` hint at the bipartition-manager → DP boundary would
  provide the same benefit seen in ASTRAL-X.

No changes are needed for correctness in either codebase; this is purely a memory-pressure management
technique that reduces peak RSS without affecting results.

---

## Checklist Before Applying

1. **Verify `DisableExplicitGC` is not set** on the target JVM.
2. **Confirm object lifetimes** by tracing every read of the candidate variable — use IDE Find Usages or `grep`.
3. **Check for early-return paths** (e.g. `--verify-*` flags that exit mid-pipeline): ensure the null
   assignment is placed *after* the last early-return that reads the object.
4. **Do not null return values in use**: inside a post-JNI processing loop, only null the input arrays,
   never the output array that the loop is iterating.
5. **Log before/after heap size** when placing a `System.gc()` hint — this confirms the mechanism is
   working and gives you a quantified win to report.
