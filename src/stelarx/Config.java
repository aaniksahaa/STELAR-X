package stelarx;

public class Config {
    public enum ComputeMode      { AUTO, CPU, GPU }
    public enum SearchMode       { LOCAL, FULL }
    /**
     * Which matrix is used to guide taxon insertion when auto-completing
     * incomplete gene trees (Phase 1b).
     *   SIMILARITY — quartet-based similarity matrix (default, matches ASTRAL-MP)
     *   DISTANCE   — topological distance matrix (legacy behaviour)
     */
    public enum CompletionMethod { SIMILARITY, DISTANCE }

    /**
     * How the GPU weight kernel counts gene-tree ↔ candidate-cluster intersections.
     * CLI: --weight-intersection-method {prefix-sum | smaller-side-traversal}
     *   PREFIX_SUM             — default. Per-tree leaf prefix sums (with its own
     *                            shared/auto-global adaptive sub-paths); each
     *                            intersection is an O(1) prefix difference. Builds
     *                            O(L) prefix working memory per active block.
     *   SMALLER_SIDE_TRAVERSAL — legacy. One thread per split; each thread walks the
     *                            smaller side of every intersection element-by-element.
     *                            No prefix sums are built at all (no per-thread state),
     *                            saving that memory.
     *   BITSET                 — low-taxa fast path. Every cluster and gene-tree part is
     *                            materialized once as a global-taxon bitset (W = ceil(n/64)
     *                            64-bit words); each core intersection is popcount(A & M)
     *                            over W words — O(1) for small n. One thread per split; no
     *                            orderings/invIndex or prefix arrays needed in the kernel.
     *                            Best when n is small and gene-tree count is large; the win
     *                            shrinks as W grows (n ≳ a few hundred). Same QI math and
     *                            bit-identical scores to the other two methods.
     *   SIMPLE_TREE_WALK       — many-candidate fast path. One thread per split walks a
     *                            resident flat postorder token stream of all gene trees
     *                            sequentially, maintaining a small O(n) per-thread stack of
     *                            (|node∩A|,|node∩B|,|node|) triples; every non-root internal
     *                            node's child partition is scored in O(1) (no
     *                            prefix arrays, no dedup, no cross-tree parallelism). Lean
     *                            kernel; wins when the candidate set is huge (full search on
     *                            large gene-tree counts), where per-split O(n·k) but tiny
     *                            constant beats the prefix path's per-split prefix rebuild.
     *                            GPU stack is capped at n; larger n falls back to CPU. Same
     *                            QI math and bit-identical scores to the other methods.
     */
    public enum WeightIntersectionMethod { PREFIX_SUM, SMALLER_SIDE_TRAVERSAL, BITSET, SIMPLE_TREE_WALK }

    /**
     * Numeric type used for weight scores when the taxon set is large enough that
     * exact 64-bit integers would overflow (see WeightTable.needsDoubleAccumulation).
     * Below that threshold scores are always exact LONG regardless of this setting.
     *   DOUBLE — 64-bit floating point (default); faster in practice despite FP64
     *            throttling on consumer GPUs, because INT128 emulation requires more
     *            instructions and 2× memory bandwidth for JNI transport.
     *   INT128 — exact 128-bit integers; use when exact scores are required
     *            (--large-n-score-type int128).
     */
    public enum LargeScoreType { INT128, DOUBLE }

    /** Set operation used by the parser-backed taxa extraction mode. */
    public enum TaxaSetMode { UNION, INTERSECTION }

    private static Config instance;

    private String inputFile;
    private String outputFile;
    private String logFile;
    private String scoreSpeciesTreeFile;
    /** Optional inference/scoring taxon allow-list (one name per non-empty line). */
    private String taxaFile;
    /** Parser-backed utility mode: write taxa from the input tree file and exit. */
    private boolean extractTaxa = false;
    private TaxaSetMode taxaSetMode = TaxaSetMode.UNION;
    /** AUTO probes the bundled CUDA backend and falls back safely to CPU. */
    private ComputeMode computeMode = ComputeMode.AUTO;
    private ComputeMode requestedComputeMode = ComputeMode.AUTO;
    private boolean gpuStrict = false;
    private String computeModeDetail = "not resolved";
    private int threadCount = Runtime.getRuntime().availableProcessors();
    private int numHashSeeds = 2;
    private long baseSeed = 0xDEADBEEFCAFEL;
    private int verbosity = 1; // 0=quiet 1=INFO 2=DEBUG 3=TRACE
    /** STELAR-X is intrinsically rooted; retained only for CLI compatibility. */
    private boolean treatAsUnrooted = false;
    /**
     * Inference-only input policy. Final triplet scoring always preserves genuine
     * input multifurcations, independently of this setting.
     */
    private boolean keepPolytomyDuringInference = false;
    private SearchMode searchMode = SearchMode.LOCAL;
    private WeightIntersectionMethod weightIntersectionMethod = WeightIntersectionMethod.PREFIX_SUM;
    private LargeScoreType largeScoreType = LargeScoreType.INT128;

    /**
     * GPU split-batching control.
     *   true  (default) — adaptive: batch size computed from free VRAM at runtime.
     *   false           — disabled: all splits sent in one kernel launch (original behaviour).
     */
    private boolean gpuBatch = true;

    /**
     * Prune candidate splits whose parent cluster is unreachable from the DP root
     * before the weight step.  The top-down inference DP only scores splits of
     * clusters it reaches from the root, so splits of never-reached clusters cost
     * weight time for nothing.  A cheap reachability BFS over the transitions graph
     * (DPTable.reachableClusters) marks the needed subset; filtering to it is
     * result-preserving.  Default on; disable with --no-prune-search-space to
     * measure the full (unfiltered) candidate count.
     */
    private boolean pruneUnreachableSplits = true;

    /** Deprecated unrooted compatibility switch; STELAR-X always rejects it. */
    private boolean anchorOutgroup = false;

    /** Global taxon id used as the outgroup anchor (default 0). */
    private int anchorTaxon = 0;

    /**
     * Manual GPU batch size override (ignored when gpuBatch=false).
     *   0 (default) — method-specific auto sizing (free-VRAM occupancy for the
     *                 legacy scorers; bounded scratch for simple-tree-walk).
     *   > 0         — use exactly this many splits per kernel launch.
     */
    private int gpuBatchSize = 0;

    /**
     * Explicit number of GPU batches (highest priority when > 0).
     * batchSize is computed as ceil(numSplits / gpuNumBatches) at runtime.
     * Overrides gpuBatchSize and the auto-VRAM logic.
     *   0 (default) — not set; fall back to gpuBatchSize or auto.
     */
    private int gpuNumBatches = 0;


    /**
     * Fraction of free VRAM (after static upload) to use for the batch buffer.
     * This is the DEFAULT auto batching mode — the native code queries free VRAM
     * via cudaMemGetInfo after uploading static data, then allocates:
     *
     *   batchSize = floor(freeVRAM × gpuVramFraction / 48 B)
     *
     * This adapts automatically to whatever GPU and dataset are in use.
     * Default 0.75 (use 75% of remaining free VRAM for splits + scores buffers,
     * leaving 25% headroom for driver, kernel stack, page tables).
     *
     * Configured via --gpu-vram-occupancy-factor.  Must be in (0, 1].
     */
    private double gpuVramFraction = 0.75;

    /**
     * VRAM control factor for GPU weight-calculation split batching.
     * Manual override — resident-relative sizing:
     *
     *   resident   = mem(orderings) + mem(invIndex) + mem(parts)
     *   mem(batch) = F × resident
     *   batchSize  = F × resident / 48 B
     *
     * Hardware-independent (same batchSize on any GPU).  Only active when
     * explicitly set via --gpu-vram-control-factor; otherwise the auto
     * free-VRAM adaptive path (gpuVramFraction) is used.
     *
     * Priority: --no-gpu-batch  >  --gpu-batches  >  --gpu-batch-size
     *         >  --gpu-vram-control-factor  >  method-specific auto sizing
     */
    private double gpuVramControlFactor    = 1.0;
    private boolean gpuVramControlFactorSet = false;

    /**
     * Maximum automatic batch scratch allocation for the GPU simple-tree-walk
     * scorer.  The optimized kernel stages both candidate sides in a transposed
     * [word][split] layout, so its scratch cost is proportional to
     * 2 * wordsPerSet * numberOfSplits.  Capping that staging area prevents a
     * large-taxon dataset from consuming most otherwise-free VRAM merely to
     * reduce the number of equivalent kernel launches.
     *
     * Manual batching flags retain priority over this automatic cap.
     * Configured via --gpu-treewalk-vram-cap-mb. Default: 512 MiB.
     */
    private int gpuTreeWalkVramCapMiB = 512;

    /**
     * GPU output buffer size for the cross-tree DP state-space construction phase
     * (Phase 5b), stored in bytes.  Each transition triple occupies 12 bytes
     * (3 × sizeof(int)), so the number of triples the buffer can hold is
     * gpuDpOutputCapBytes / 12.
     *
     * Default: 120 MB = 10 000 000 triples.
     *
     * Sub-batching normally guarantees no overflow, but if a dataset has an
     * extraordinarily large single size-bin the kernel will overflow and print
     * a CRITICAL WARNING.  Raise this value (e.g. "1g") to avoid the overflow
     * at the cost of more VRAM.
     *
     * Configured via --gpu-dp-state-space-construction-output-cap.
     * Accepts memory-unit suffixes: k/K (×10³), m/M (×10⁶), g/G (×10⁹).
     * Examples: "120m"  "1.2g"  "500k"  "1500000000"
     */
    private long gpuDpOutputCapBytes = 128_000_000L; // 128 MB default

    /**
     * Minimum seconds between DP progress bar updates.
     * Configured via --gpu-dp-state-space-progress-time-interval.
     */
    private double gpuDpProgressInterval = 1.0;

    /**
     * Maximum number of progress bar print steps for DP phase.
     * When > 0, switches to step-based mode (% advancement) and disables the
     * time-interval trigger entirely.
     * 0 = not set; use time-interval mode.
     * Default: 1000 (step mode, print every 0.1% advancement).
     * Configured via --gpu-dp-state-space-progress-max-steps.
     */
    private int gpuDpProgressMaxSteps = 1000;

    private Config() {}

    public static Config getInstance() {
        if (instance == null) instance = new Config();
        return instance;
    }
    public static void reset() { instance = null; }

    public String getInputFile()      { return inputFile; }
    public void setInputFile(String f){ this.inputFile = f; }
    public String getOutputFile()     { return outputFile; }
    public void setOutputFile(String f){ this.outputFile = f; }
    public String getLogFile()        { return logFile; }
    public void setLogFile(String f)  { this.logFile = f; }
    public String getScoreSpeciesTreeFile()      { return scoreSpeciesTreeFile; }
    public void setScoreSpeciesTreeFile(String f){ this.scoreSpeciesTreeFile = f; }
    public boolean isScoreOnly()       { return scoreSpeciesTreeFile != null; }
    public String getTaxaFile()        { return taxaFile; }
    public void setTaxaFile(String f)  { this.taxaFile = f; }
    public boolean isExtractTaxa()     { return extractTaxa; }
    public void setExtractTaxa(boolean v) { this.extractTaxa = v; }
    public TaxaSetMode getTaxaSetMode() { return taxaSetMode; }
    public void setTaxaSetMode(TaxaSetMode m) { this.taxaSetMode = m; }
    public ComputeMode getComputeMode()          { return computeMode; }
    public ComputeMode getRequestedComputeMode() { return requestedComputeMode; }
    public void setComputeMode(ComputeMode m) {
        this.computeMode = m;
        this.requestedComputeMode = m;
    }
    /** Set the resolved mode without losing whether AUTO/GPU/CPU was requested. */
    public void resolveComputeMode(ComputeMode m, String detail) {
        this.computeMode = m;
        this.computeModeDetail = detail == null ? "" : detail;
    }
    public String getComputeModeDetail()       { return computeModeDetail; }
    public boolean isGpuStrict()               { return gpuStrict; }
    public void setGpuStrict(boolean v)         { this.gpuStrict = v; }
    public int getThreadCount()               { return threadCount; }
    public void setThreadCount(int t)         { this.threadCount = Math.max(1, t); }
    public int getNumHashSeeds()              { return numHashSeeds; }
    public void setNumHashSeeds(int m)        { this.numHashSeeds = m; }
    public long getBaseSeed()                 { return baseSeed; }
    public int getVerbosity()                 { return verbosity; }
    public void setVerbosity(int v)           { this.verbosity = v; }
    public boolean getTreatAsUnrooted()       { return treatAsUnrooted; }
    public void setTreatAsUnrooted(boolean u) { this.treatAsUnrooted = u; }
    public boolean isKeepPolytomyDuringInference() {
        return keepPolytomyDuringInference;
    }
    public void setKeepPolytomyDuringInference(boolean keep) {
        this.keepPolytomyDuringInference = keep;
    }
    public SearchMode getSearchMode()          { return searchMode; }
    public void setSearchMode(SearchMode s)   { this.searchMode = s; }
    public WeightIntersectionMethod getWeightIntersectionMethod()        { return weightIntersectionMethod; }
    public void setWeightIntersectionMethod(WeightIntersectionMethod m)  { this.weightIntersectionMethod = m; }

    public LargeScoreType getLargeScoreType()        { return largeScoreType; }
    public void setLargeScoreType(LargeScoreType t)  { this.largeScoreType = t; }
    public boolean isGpuBatch()               { return gpuBatch; }
    public boolean isPruneUnreachableSplits()           { return pruneUnreachableSplits; }
    public void setPruneUnreachableSplits(boolean v)    { this.pruneUnreachableSplits = v; }
    public boolean isAnchorOutgroup()                   { return anchorOutgroup; }
    public void setAnchorOutgroup(boolean v)            { this.anchorOutgroup = v; }
    public int getAnchorTaxon()                         { return anchorTaxon; }
    public void setAnchorTaxon(int t)                   { this.anchorTaxon = t; }
    public void setGpuBatch(boolean b)        { this.gpuBatch = b; }
    public int getGpuBatchSize()              { return gpuBatchSize; }
    public void setGpuBatchSize(int s)        { this.gpuBatchSize = s; }
    public int getGpuNumBatches()             { return gpuNumBatches; }
    public void setGpuNumBatches(int n)       { this.gpuNumBatches = Math.max(1, n); }
    public double getGpuVramFraction()            { return gpuVramFraction; }
    public void setGpuVramFraction(double f)      { this.gpuVramFraction = Math.max(0.01, Math.min(1.0, f)); }

    // GPU weight-kernel progress-bar update interval (seconds).  -1 = auto (TTY ~2 s
    // overwriting line; non-TTY ~300 s newline).  Precedence: this flag > env
    // STELARX_GPU_PROGRESS_SEC > auto default.
    private double gpuProgressIntervalSec = -1.0;
    public double getGpuProgressIntervalSec()        { return gpuProgressIntervalSec; }
    public void setGpuProgressIntervalSec(double s)  { this.gpuProgressIntervalSec = s; }
    public double getGpuVramControlFactor()       { return gpuVramControlFactor; }
    public boolean isGpuVramControlFactorSet()    { return gpuVramControlFactorSet; }
    public void setGpuVramControlFactor(double f) { this.gpuVramControlFactor = Math.max(0.001, Math.min(1.0, f)); this.gpuVramControlFactorSet = true; }
    public int getGpuTreeWalkVramCapMiB()         { return gpuTreeWalkVramCapMiB; }
    public void setGpuTreeWalkVramCapMiB(int cap) { this.gpuTreeWalkVramCapMiB = Math.max(1, cap); }
    /** Raw byte count of the GPU DP output buffer. */
    public long getGpuDpOutputCapBytes()      { return gpuDpOutputCapBytes; }

    /**
     * Number of transition triples the GPU output buffer can hold
     * (= bytes / 12, clamped to at least 1).
     */
    public int getGpuDpOutputCapTriples()     { return (int) Math.max(1, gpuDpOutputCapBytes / 12); }

    /**
     * Set the GPU DP output-buffer cap from a human-readable memory string.
     * Accepts optional suffixes k/K (×1 000), m/M (×1 000 000), g/G (×1 000 000 000).
     * A bare integer is interpreted as bytes.
     * Examples: "120m", "1.2g", "500k", "1500000000"
     */
    public void setGpuDpStateSpaceConstructionOutputCap(String spec) {
        String s = spec.trim().toLowerCase();
        double value;
        long multiplier;
        if (s.endsWith("g")) {
            value = Double.parseDouble(s.substring(0, s.length() - 1));
            multiplier = 1_000_000_000L;
        } else if (s.endsWith("m")) {
            value = Double.parseDouble(s.substring(0, s.length() - 1));
            multiplier = 1_000_000L;
        } else if (s.endsWith("k")) {
            value = Double.parseDouble(s.substring(0, s.length() - 1));
            multiplier = 1_000L;
        } else {
            value = Double.parseDouble(s);
            multiplier = 1L;
        }
        long bytes = (long)(value * multiplier);
        this.gpuDpOutputCapBytes = Math.max(12L, bytes); // at least 1 triple
    }

    public double getGpuDpProgressInterval()          { return gpuDpProgressInterval; }
    public void setGpuDpProgressInterval(double v)    { this.gpuDpProgressInterval = Math.max(0.0, v); }
    public int getGpuDpProgressMaxSteps()             { return gpuDpProgressMaxSteps; }
    public void setGpuDpProgressMaxSteps(int v)       { this.gpuDpProgressMaxSteps = Math.max(1, v); }

    // Completion flags
    private boolean autoCompleteIncompleteTrees = false;
    public boolean isAutoCompleteIncompleteTrees()          { return autoCompleteIncompleteTrees; }
    public void setAutoCompleteIncompleteTrees(boolean v)   { this.autoCompleteIncompleteTrees = v; }

    /** Which matrix guides taxon insertion in tree completion. Default: SIMILARITY. */
    private CompletionMethod completionMethod = CompletionMethod.SIMILARITY;
    public CompletionMethod getCompletionMethod()             { return completionMethod; }
    public void setCompletionMethod(CompletionMethod m)       { this.completionMethod = m; }

    /**
     * Tile side-length B for the GPU distance-matrix kernel.
     * Controls GPU VRAM for the output tile: B² × 12 bytes.
     * Default 0 = auto: B = min(n, ceil(sqrt(n * k))), capped by available VRAM.
     * Configured via --gpu-dist-tile-size.
     */
    private int gpuDistTileSizeB = 0;
    public int  getGpuDistTileSizeB()          { return gpuDistTileSizeB; }
    public void setGpuDistTileSizeB(int b)     { this.gpuDistTileSizeB = Math.max(0, b); }

    /**
     * Maximum GPU memory reserved for one similarity-matrix tree-data batch.
     * The output tile is separate and tiny for the usual completion datasets.
     * A bounded default avoids consuming a fixed fraction of a large GPU merely
     * to reduce the number of otherwise equivalent tree batches.
     * Configured via --gpu-sim-vram-cap-mb. Default: 512 MiB for the established
     * dense path. The automatic large-N packed path may raise its batching
     * ceiling when the user did not explicitly set this option; native code
     * still clamps it to currently free VRAM.
     */
    private int gpuSimilarityVramCapMiB = 512;
    private boolean gpuSimilarityVramCapExplicit = false;
    public int  getGpuSimilarityVramCapMiB()          { return gpuSimilarityVramCapMiB; }
    public void setGpuSimilarityVramCapMiB(int cap)   {
        this.gpuSimilarityVramCapMiB = Math.max(1, cap);
        this.gpuSimilarityVramCapExplicit = true;
    }
    public boolean isGpuSimilarityVramCapExplicit()  { return gpuSimilarityVramCapExplicit; }

    // Standalone deployment self-check. Does not require an input dataset.
    private boolean diagnose = false;
    public boolean isDiagnose()          { return diagnose; }
    public void setDiagnose(boolean v)   { this.diagnose = v; }

    // Testing flags
    private boolean verifyParse = false;
    public boolean isVerifyParse()          { return verifyParse; }
    public void setVerifyParse(boolean v)   { this.verifyParse = v; }

    private boolean verifyHash = false;
    public boolean isVerifyHash()           { return verifyHash; }
    public void setVerifyHash(boolean v)    { this.verifyHash = v; }

    private boolean verifyClusters = false;
    public boolean isVerifyClusters()       { return verifyClusters; }
    public void setVerifyClusters(boolean v){ this.verifyClusters = v; }

    private boolean verifyPartitions = false;
    public boolean isVerifyPartitions()        { return verifyPartitions; }
    public void setVerifyPartitions(boolean v) { this.verifyPartitions = v; }

    private boolean verifyDPSpace = false;
    public boolean isVerifyDPSpace()           { return verifyDPSpace; }
    public void setVerifyDPSpace(boolean v)    { this.verifyDPSpace = v; }

    private boolean verifyWeights = false;
    public boolean isVerifyWeights()           { return verifyWeights; }
    public void setVerifyWeights(boolean v)    { this.verifyWeights = v; }

    private boolean verifyDistanceMatrix = false;
    public boolean isVerifyDistanceMatrix()          { return verifyDistanceMatrix; }
    public void setVerifyDistanceMatrix(boolean v)   { this.verifyDistanceMatrix = v; }

    private boolean verifySimilarityMatrix = false;
    public boolean isVerifySimilarityMatrix()        { return verifySimilarityMatrix; }
    public void setVerifySimilarityMatrix(boolean v) { this.verifySimilarityMatrix = v; }

    private boolean verifyUpgma = false;
    public boolean isVerifyUpgma()          { return verifyUpgma; }
    public void setVerifyUpgma(boolean v)   { this.verifyUpgma = v; }

    private boolean verifyGreedyConsensus = false;
    public boolean isVerifyGreedyConsensus()        { return verifyGreedyConsensus; }
    public void setVerifyGreedyConsensus(boolean v) { this.verifyGreedyConsensus = v; }

    // Greedy consensus build + polytomy resolution + emission to X is an
    // INCOMPLETE, EXPERIMENTAL feature.  It is OFF by default: Phase 3.5 is
    // skipped entirely (no compute, no memory) unless this flag is set.
    private boolean consensusExperimental = false;
    public boolean isConsensusExperimental()        { return consensusExperimental; }
    public void setConsensusExperimental(boolean v) { this.consensusExperimental = v; }

    // Step B (sampleAndResolve) per-gene-tree restriction route.
    //   true  (default) = O(d log d) induced/auxiliary-tree via Euler-tour LCA
    //   false           = O(n) full postorder walk (legacy reference route)
    // Both produce the identical emission set; the fast route avoids touching all
    // n leaves when only d≤31 reps matter.
    private boolean stepBFastRestriction = true;
    public boolean isStepBFastRestriction()        { return stepBFastRestriction; }
    public void setStepBFastRestriction(boolean v) { this.stepBFastRestriction = v; }

    // ── Optional extra X-enrichment around polytomies (OFF by default) ──────────
    // Both reproduce ASTRAL-MP heuristics but materially ENLARGE the candidate set X
    // (and thus DP + weight cost), so they are opt-in.

    // D1: ASTRAL-MP getQuadraticBitsets — for each polytomy arm, add the nested
    // nearest-neighbour "balls" (k-NN unions by induced rep similarity) as candidate
    // (multi-range) clusters. O(d²) extra candidates per polytomy/round.
    private boolean stepBQuadraticNnBalls = false;
    public boolean isStepBQuadraticNnBalls()        { return stepBQuadraticNnBalls; }
    public void setStepBQuadraticNnBalls(boolean v) { this.stepBQuadraticNnBalls = v; }

    // D2: ASTRAL-MP resolveLinearly leftover step — after the sampled mini-greedy
    // build, randomly pair-merge any still-unresolved multifurcation and add the
    // intermediate unions as candidate clusters.
    private boolean stepBRandomLeftoverResolution = false;
    public boolean isStepBRandomLeftoverResolution()        { return stepBRandomLeftoverResolution; }
    public void setStepBRandomLeftoverResolution(boolean v) { this.stepBRandomLeftoverResolution = v; }

    // "Lift the bar": process polytomies of ANY degree, matching ASTRAL-MP, which
    // never drops a polytomy — its polytomySizeLimit only disables the quadratic
    // NN-balls for over-limit polytomies; the linear path (Step A UPGMA + Step B
    // mini-greedy + UPGMA-on-reps) still runs for every polytomy.
    //   OFF (default): pool drops d > sizeLimit, and Step B is capped at d ≤ 31.
    //   ON: every polytomy is processed.  d ≤ 31 uses the unchanged int path;
    //       d > 31 uses a long[]-bitmap Step B path; UPGMA on g > 31 uses the
    //       exact O(d²) nearest-neighbour-chain (MiniUPGMA.buildFast).  Quadratic
    //       NN-balls remain disabled for d > 31 (mirrors ASTRAL-MP's size gate).
    // Turning this ON is purely ADDITIVE: small-polytomy emissions are identical;
    // only the large-polytomy candidates are added.  ENLARGES X (and DP/weight cost).
    private boolean stepBProcessLargePolytomies = false;
    public boolean isStepBProcessLargePolytomies()        { return stepBProcessLargePolytomies; }
    public void setStepBProcessLargePolytomies(boolean v) { this.stepBProcessLargePolytomies = v; }

    // Gene-tree polytomy X-enrichment (ASTRAL-MP "mechanism B",
    // addBipartitionsFromSignleIndTreesToX :172-227): resolve each INPUT gene-tree
    // polytomy against the UPGMA guide tree (3 samples) and add the resulting
    // arm-union (multi-range) clusters to X.  Distinct from the d-partition QI
    // SCORING (always on for polytomous inputs).  Default OFF: it enlarges X.
    private boolean resolveInputGeneTreePolytomies = false;
    public boolean isResolveInputGeneTreePolytomies()        { return resolveInputGeneTreePolytomies; }
    public void setResolveInputGeneTreePolytomies(boolean v) { this.resolveInputGeneTreePolytomies = v; }

    // ── Comparison / debug dump flags ─────────────────────────────────────────
    private String dumpClustersFile = null;
    public String getDumpClustersFile()       { return dumpClustersFile; }
    public void setDumpClustersFile(String f) { this.dumpClustersFile = f; }

    private String dumpCompletedTreesFile = null;
    public String getDumpCompletedTreesFile()       { return dumpCompletedTreesFile; }
    public void setDumpCompletedTreesFile(String f) { this.dumpCompletedTreesFile = f; }
}
