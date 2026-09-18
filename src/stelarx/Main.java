package stelarx;

import stelarx.cluster.ClusterTable;
import stelarx.completion.DistanceMatrix;
import stelarx.completion.DistanceMatrixBuilder;
import stelarx.completion.SimilarityMatrix;
import stelarx.completion.SimilarityMatrixBuilder;
import stelarx.completion.TreeCompleter;
import stelarx.completion.UPGMAClusterer;
import stelarx.gpu.GPUDistanceMatrix;
import stelarx.gpu.GPUSimilarityMatrix;
import stelarx.dp.DPTable;
import stelarx.dp.Inference;
import stelarx.greedy.EmissionBridge;
import stelarx.greedy.GreedyConsensus;
import stelarx.greedy.GreedyConsensusVerifier;
import stelarx.gpu.GPUDPBuilder;
import stelarx.gpu.GPUWeightCalculator;
import stelarx.partition.PartitionTable;
import stelarx.weight.WeightTable;
import stelarx.hash.PrefixHashArrays;
import stelarx.hash.TaxonHasher;
import stelarx.taxon.TaxonRegistry;
import stelarx.tree.Tree;
import stelarx.tree.TreeParser;
import stelarx.tree.TreeRestrictor;
import stelarx.tree.TreeTaxa;
import stelarx.util.Threading;

import stelarx.cluster.Cluster;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class Main {
    public static final String VERSION = Version.current();

    public static void main(String[] args) {
        FatalReporter.prepareCrashDirectory();
        TerminalLog terminalLog;
        try {
            terminalLog = TerminalLog.installFromArgs(args);
        } catch (TerminalLog.SetupException e) {
            System.err.println("Error: " + e.getMessage());
            System.exit(2);
            return;
        }
        try {
            run(args);
        } catch (Throwable t) {
            FatalReporter.report(t, args);
            System.exit(1);
        } finally {
            if (terminalLog != null) terminalLog.close();
        }
    }

    private static void run(String[] args) throws Exception {
        Config cfg = Config.getInstance();
        if (!parseArgs(args, cfg)) { printUsage(); System.exit(1); }
        if (cfg.getTreatAsUnrooted()) {
            throw new IllegalArgumentException(
                "STELAR-X accepts rooted gene trees only; --unrooted is not supported");
        }
        if (cfg.isAnchorOutgroup()) {
            throw new IllegalArgumentException(
                "--anchor-outgroup is an unrooted optimization and is not valid in STELAR-X");
        }
        validatePathCollisions(cfg);

        Logging.setLevel(cfg.getVerbosity());
        if (cfg.isExtractTaxa()) {
            runTaxaExtraction(cfg);
            return;
        }
        resolveComputeMode(cfg);
        Banner.print(cfg);

        if (cfg.isDiagnose()) {
            RuntimeDiagnostics.print(cfg);
            return;
        }

        Threading.start(cfg.getThreadCount());
        long t0 = System.nanoTime();
        String finalTripletScore = null;
        boolean analysisCompleted = false;
        boolean inferenceInputHasPolytomy = false;

        try {
            if (cfg.isScoreOnly() && cfg.getTaxaFile() != null) {
                finalTripletScore = runTaxonRestrictedScoreOnly(cfg);
                analysisCompleted = true;
                return;
            }

            // ── Phase 1: Parse gene trees ─────────────────────────────────────
            long t1 = PhaseLogger.begin("Phase 1  Parse gene trees", false);
            TaxonRegistry registry;
            List<Tree> trees;
            if (cfg.getTaxaFile() == null) {
                registry = new TaxonRegistry();
                boolean keepPolytomy = cfg.isScoreOnly()
                    || cfg.isKeepPolytomyDuringInference();
                TreeParser.ParsedGeneTrees parsed = TreeParser.parseGeneTreesDetailed(
                    cfg.getInputFile(), registry, keepPolytomy);
                trees = parsed.trees();
                inferenceInputHasPolytomy = parsed.detectedPolytomyNodeCount() > 0;
            } else {
                RestrictedInferenceInput restricted =
                    parseTaxonRestrictedInferenceInput(cfg);
                registry = restricted.registry();
                trees = restricted.trees();
                inferenceInputHasPolytomy = restricted.hasPolytomy();
            }
            PhaseLogger.end("Phase 1  Parse gene trees", t1, false);

            if (cfg.isVerifyParse()) {
                Phase1Verifier.dump(trees, registry, cfg.getOutputFile());
                return;
            }

            // ── Phase 1b: Auto-complete incomplete gene trees (optional) ──────
            // Entered only when --autocomplete-incomplete-gene-trees or
            // --verify-distance-matrix is explicitly requested.  The baseline
            // (complete trees, no flag) skips this block entirely — no library
            // load, no stream scan, zero overhead.
            //
            // IMPORTANT: originalTrees is saved BEFORE completion and is used for
            // rooted child-partition extraction (Phase 4) and weight calculation (Phase 6).
            // X (ClusterTable) and DP transitions are built from the completed trees
            // so that all bipartitions span the full taxon set — exactly what ASTRAL-MP
            // does.  Weight calculation must use the ORIGINAL gene trees (as ASTRAL-MP
            // does via inference.trees = originalInompleteGeneTrees) so triplet scores
            // reflect actual gene-tree signal, not the artificially inserted taxa.
            List<Tree> originalTrees = trees; // always points to pre-completion trees
            SimilarityMatrix similarityMatrix = null; // visible to Phase 3.5 (Step A)
            Tree upgmaGuideTree = null; // retained for S3 gene-tree-polytomy enrichment
            if (cfg.isAutoCompleteIncompleteTrees() || cfg.isVerifyDistanceMatrix()
                    || cfg.isVerifySimilarityMatrix() || cfg.isVerifyUpgma()) {
                boolean gpuDist = (cfg.getComputeMode() == Config.ComputeMode.GPU)
                                  && GPUDistanceMatrix.tryLoad();
                boolean gpuSim  = (cfg.getComputeMode() == Config.ComputeMode.GPU)
                                  && GPUSimilarityMatrix.tryLoad();

                if (cfg.isVerifyDistanceMatrix()) {
                    DistanceMatrix dm = gpuDist
                        ? DistanceMatrixBuilder.buildGPU(trees, registry.size())
                        : DistanceMatrixBuilder.buildCPU(trees, registry.size());
                    dumpDistanceMatrix(dm, registry);
                    return;
                }

                if (cfg.isVerifySimilarityMatrix()) {
                    SimilarityMatrix sm = gpuSim
                        ? SimilarityMatrixBuilder.buildGPU(trees, registry.size())
                        : SimilarityMatrixBuilder.buildCPU(trees, registry.size());
                    dumpSimilarityMatrix(sm, registry);
                    return;
                }

                if (cfg.isVerifyUpgma()) {
                    SimilarityMatrix sm = gpuSim
                        ? SimilarityMatrixBuilder.buildGPU(trees, registry.size())
                        : SimilarityMatrixBuilder.buildCPU(trees, registry.size());
                    int n = registry.size();
                    Tree upgmaTree = UPGMAClusterer.build(sm, trees.size());
                    dumpUpgmaBipartitions(upgmaTree, registry);
                    return;
                }

                long incompleteCount = trees.stream().filter(t -> !t.isComplete).count();
                boolean useSim = cfg.getCompletionMethod() == Config.CompletionMethod.SIMILARITY;
                boolean gpuActive = useSim ? gpuSim : gpuDist;
                long t1b = PhaseLogger.begin("Phase 1b Auto-complete gene trees ("
                    + (useSim ? "similarity" : "distance") + ") + UPGMA guide", gpuActive);

                // Always build the similarity matrix: needed for UPGMA guide tree regardless of
                // completion method or whether any trees are actually incomplete.
                SimilarityMatrix smForUpgma = gpuSim
                    ? SimilarityMatrixBuilder.buildGPU(trees, registry.size())
                    : SimilarityMatrixBuilder.buildCPU(trees, registry.size());
                similarityMatrix = smForUpgma; // retained for Phase 3.5 (greedy consensus Step A)

                if (incompleteCount > 0) {
                    // The four-point algorithm always needs the similarity matrix for the
                    // scoring formula (sim[x][a] + sim[b][c] - ...).  smForUpgma.sim is
                    // already available regardless of completionMethod.
                    // dist is used only to build sortedRows (nearest-neighbour order);
                    // for SIMILARITY mode we reuse smForUpgma.dist (= 1 - sim).
                    if (useSim) {
                        trees = TreeCompleter.completeAll(trees, smForUpgma, registry.size());
                    } else {
                        if (smForUpgma.isPacked()) {
                            throw new IllegalArgumentException("Large-N distance-guided completion "
                                + "requires a segmented distance matrix; use the exact default "
                                + "--completion-method similarity for this dataset");
                        }
                        DistanceMatrix dm = gpuDist
                            ? DistanceMatrixBuilder.buildGPU(trees, registry.size())
                            : DistanceMatrixBuilder.buildCPU(trees, registry.size());
                        trees = TreeCompleter.completeAll(trees, smForUpgma.sim, dm.dist,
                            registry.size());
                    }
                    Logging.info("Phase 1b: using original incomplete trees for weight scoring, completed trees for X");

                    if (cfg.getDumpCompletedTreesFile() != null) {
                        dumpCompletedTrees(trees, registry, cfg.getDumpCompletedTreesFile());
                    }
                } else {
                    Logging.info("Phase 1b: all gene trees already complete");
                }

                // Build UPGMA guide tree from similarity matrix and append to the completed
                // trees list so Phase 2/3 include its bipartitions in the cluster set X.
                // It is NOT added to originalTrees, so child-partition scoring (Phase 4/6)
                // is unaffected.
                int nTaxa = registry.size();
                upgmaGuideTree = UPGMAClusterer.build(smForUpgma, trees.size());
                trees = new ArrayList<>(trees);
                trees.add(upgmaGuideTree);
                Logging.info("Phase 1b: UPGMA guide tree (%d taxa) added to cluster search space", nTaxa);

                PhaseLogger.end("Phase 1b Auto-complete gene trees", t1b, gpuActive);
            }
            // After Phase 1b:
            //   trees         = completed gene trees (or original if no autocomplete / no incomplete)
            //   originalTrees = original gene trees (same reference as trees when no autocomplete)

            // ── Phase 2: Taxon hashing + prefix arrays ────────────────────────
            // pref     — built from completed trees; used for ClusterTable and DPTable
            // prefParts — built from original trees; used for rooted child-partition scoring
            // When no autocomplete (originalTrees == trees), prefParts == pref (same object).
            long t2 = PhaseLogger.begin("Phase 2  Taxon hashing", false);
            TaxonHasher hasher = new TaxonHasher(
                registry.size(), cfg.getNumHashSeeds(), cfg.getBaseSeed());
            PrefixHashArrays pref = new PrefixHashArrays(trees, hasher);
            PrefixHashArrays prefParts = (originalTrees == trees)
                ? pref
                : new PrefixHashArrays(originalTrees, hasher);
            PhaseLogger.end("Phase 2  Taxon hashing", t2, false);

            if (cfg.isVerifyHash()) {
                Phase2Verifier.dump(trees, registry, hasher, pref, cfg.getOutputFile());
                return;
            }
            // hasher is retained — Phase 3.5 (greedy consensus) needs per-taxon
            // hashes to build consensus-tree prefix arrays whose signatures
            // match those derived from the gene-tree prefix arrays (cross-source
            // signature parity, design §7.2 / verification §13.5).

            if (cfg.isScoreOnly()) {
                finalTripletScore = runScoreOnly(
                    cfg, registry, originalTrees, prefParts, hasher);
                analysisCompleted = true;
                return;
            }

            // ── Phase 3: Cluster extraction -> X (from COMPLETED trees) ──────
            // Anchor-free X stores only the orientation of each bipartition that
            // excludes the fixed anchor.  Both local and full modes use the
            // corresponding single anchored root when the option is enabled.
            boolean anchorFreeX = cfg.isAnchorOutgroup();
            long t3 = PhaseLogger.begin("Phase 3  Cluster extraction", false);
            ClusterTable clusterTable = new ClusterTable(trees, pref, registry.size(), anchorFreeX);
            PhaseLogger.end("Phase 3  Cluster extraction", t3, false);

            if (cfg.isVerifyClusters()) {
                Phase3Verifier.dump(trees, registry, pref, clusterTable, cfg.getOutputFile());
                return;
            }

            if (cfg.getDumpClustersFile() != null) {
                dumpClusters(clusterTable, trees, registry, cfg.getDumpClustersFile());
            }

            // ── Phase 3.6: Gene-tree polytomy X-enrichment (mechanism B, opt-in) ──
            // Resolve each INPUT gene-tree polytomy against the UPGMA guide tree into
            // arm-union (multi-range) clusters (ASTRAL-MP addBipartitionsFromSignleIndTreesToX).
            // Distinct from rooted-partition scoring; gated since it enlarges X.
            if (cfg.isResolveInputGeneTreePolytomies()
                    && anyPolytomous(trees, originalTrees.size())) {
                long t36 = PhaseLogger.begin("Phase 3.6 Gene-tree polytomy enrichment", false);
                int nT = registry.size();
                if (similarityMatrix == null) {
                    similarityMatrix = SimilarityMatrixBuilder.buildCPU(trees, nT);
                }
                Tree guide = upgmaGuideTree != null
                    ? upgmaGuideTree
                    : UPGMAClusterer.build(similarityMatrix, trees.size());
                stelarx.greedy.GeneTreePolytomySampler.run(
                    trees, originalTrees.size(), guide, pref, nT, pref.numSeeds(),
                    cfg.getBaseSeed() ^ 0xC0FFEEL, clusterTable);
                PhaseLogger.end("Phase 3.6 Gene-tree polytomy enrichment", t36, false);
            }

            // ── Phase 3.5: Greedy consensus + polytomy resolution (EXPERIMENTAL) ──
            // INCOMPLETE feature — emission of resolved polytomies into X is not
            // yet wired up.  Skipped entirely by default: no compute, no memory.
            // Enabled only via --consensus-experimental (build path) or the
            // explicit --verify-greedy-consensus entry point.
            //
            // Use ONLY the gene trees (no UPGMA guide tree) so the bipartition
            // frequencies match ASTRAL-MP's `addExtraBipartitionByHeuristics`
            // (which runs greedy consensus over the gene trees alone).
            //
            // `originalTrees.size()` is the gene-tree count regardless of
            // completion mode:
            //   - autocomplete OFF:  trees == originalTrees (no UPGMA appended)
            //   - autocomplete ON :  trees = completed gene trees + UPGMA
            //                        originalTrees still references the
            //                        pre-completion / pre-UPGMA list size.
            // Exemplar list the WEIGHT phase uses for cluster position lookups.
            // Defaults to the completed/gene trees; the consensus bridge may append
            // consensus snapshot trees (membership-only exemplars) for synthesized
            // multi-range clusters. The DP (Phase 5) keeps using `trees` directly,
            // so consensus trees are never mined for local transitions.
            List<Tree> weightClusterTrees = trees;
            if (cfg.isVerifyGreedyConsensus() || cfg.isConsensusExperimental()) {
                List<Tree> geneTreesForGreedy = trees.subList(0, originalTrees.size());

                if (cfg.isVerifyGreedyConsensus()) {
                    GreedyConsensusVerifier.dump(geneTreesForGreedy, registry, clusterTable,
                                                 pref, hasher, similarityMatrix, cfg.getOutputFile());
                    return;
                }
                long t35 = PhaseLogger.begin("Phase 3.5 Greedy consensus build + polytomy resolution", false);
                GreedyConsensus.Result gcResult =
                    GreedyConsensus.build(clusterTable, geneTreesForGreedy, pref, hasher,
                                           similarityMatrix, registry.size());
                PhaseLogger.end("Phase 3.5 Greedy consensus build + polytomy resolution", t35, false);

                // ── Bridge emissions into X (Tier-1 lookup / Tier-2 synthesize) ──
                List<Tree> ext = new ArrayList<>(trees);
                int[] bridged = EmissionBridge.bridge(gcResult.emissions, clusterTable,
                                                      ext, registry.size(),
                                                      anchorFreeX, cfg.getAnchorTaxon());
                // Only switch to the extended list if exemplar trees were actually
                // appended — preserves the `trees == originalTrees` identity that the
                // weight path's autocomplete detection relies on when nothing changed.
                if (ext.size() > trees.size()) weightClusterTrees = ext;
                Logging.info("Consensus emission → X: %d already in X (tier-1), "
                    + "%d synthesized multi-range (tier-2); +%d exemplar trees",
                    bridged[0], bridged[1], ext.size() - trees.size());
                if (cfg.getSearchMode() != Config.SearchMode.FULL && bridged[1] > 0) {
                    Logging.info("Note: synthesized multi-range clusters gain DP transitions only "
                        + "via Mode 2 (--search-mode full); in local mode they remain inert.");
                }
            }

            // Overflow-streamed similarity preprocessing uses large, short-lived
            // Java arrays.  Nothing after Phase 3.5 consumes the matrix, so drop
            // its final reference and request one collection before Phase 4's
            // allocation-heavy candidate scan.  This is deliberately gated to
            // the 100k-class overflow path: established one-shot runs receive no
            // extra full-GC pause. The JVM may ignore the hint when explicit GC
            // is disabled, without affecting correctness.
            boolean reclaimStreamedSimilarity = similarityMatrix != null
                && similarityMatrix.usedStreamedHostBatches();
            // The matrix has no consumers after Phase 3.5. Always end its
            // lifetime here; only overflow-streamed runs pay for an immediate GC.
            similarityMatrix = null;
            if (reclaimStreamedSimilarity) {
                long heapBefore = Runtime.getRuntime().totalMemory()
                    - Runtime.getRuntime().freeMemory();
                System.gc();
                long heapAfter = Runtime.getRuntime().totalMemory()
                    - Runtime.getRuntime().freeMemory();
                Logging.info("Released streamed similarity host buffers before Phase 4: "
                        + "heap %d MiB -> %d MiB",
                    heapBefore >> 20, heapAfter >> 20);
            }

            // ── Phase 4: Gene-tree child-partition extraction (from ORIGINAL trees) ──
            // Uses originalTrees so rooted child partitions reflect actual gene-tree signal.
            long t4 = PhaseLogger.begin("Phase 4  Rooted child-partition extraction", false);
            PartitionTable partTable = new PartitionTable(originalTrees, prefParts);
            PhaseLogger.end("Phase 4  Rooted child-partition extraction", t4, false);

            if (cfg.isVerifyPartitions()) {
                Phase4Verifier.dump(originalTrees, registry, prefParts, partTable, cfg.getOutputFile());
                return;
            }

            // ── Phase 5: DP search space (from COMPLETED trees) ───────────────
            long t5 = PhaseLogger.begin("Phase 5  DP local transitions", false);
            DPTable dpTable = new DPTable(trees, pref, clusterTable);
            PhaseLogger.end("Phase 5  DP local transitions", t5, false);

            // ── Phase 5b: Cross-tree transitions (Mode 2, optional) ───────────
            if (cfg.getSearchMode() == Config.SearchMode.FULL) {
                boolean gpuDP = (cfg.getComputeMode() == Config.ComputeMode.GPU)
                                && GPUDPBuilder.tryLoad();
                long t5b = PhaseLogger.begin("Phase 5b Cross-tree transitions", gpuDP);
                dpTable.addCrossTreeTransitions(clusterTable, gpuDP);
                PhaseLogger.end("Phase 5b Cross-tree transitions", t5b, gpuDP);
            }

            // ── Anchored-outgroup root (exact for unrooted inference) ─────────
            // Replace the all-taxa root's entire transition set with the single
            // split ({anchor} | S\{anchor}).  Every unrooted tree is representable
            // rooted on the anchor's pendant edge, so the optimum is unchanged; the
            // reachability prune then drops the now-unreachable with-anchor cluster
            // orientations. Local emission includes the leaf-induced Type-2
            // rotation needed to enter S\{anchor}; full mode additionally has all
            // cross-tree resolutions.
            if (cfg.isAnchorOutgroup()) {
                dpTable.applyAnchoredRoot(clusterTable.getAnchorHash());
            }

            if (cfg.isVerifyDPSpace()) {
                Phase5Verifier.dump(trees, registry, pref, clusterTable, dpTable, cfg.getOutputFile());
                return;
            }
            pref = null;     // no longer needed after Phase 5; free before Phase 6
            prefParts = null; // likewise (may be same object as pref — both nulled safely)

            // Hint JVM to collect Phase 3-5 intermediates before Phase 6 allocates its working set.
            // System.gc() is a hint — JVM may ignore it if -XX:+DisableExplicitGC is set.
            long gcHeapBefore = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
            System.gc();
            long gcHeapAfter = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
            Logging.debug("Pre-Phase-6 GC hint: heap %d MB → %d MB (freed %d MB)",
                gcHeapBefore / 1_000_000, gcHeapAfter / 1_000_000,
                (gcHeapBefore - gcHeapAfter) / 1_000_000);

            // ── Phase 6: Weight calculation ───────────────────────────────────
            boolean gpuWeight = (cfg.getComputeMode() == Config.ComputeMode.GPU)
                                && GPUWeightCalculator.isLoaded();
            long t6 = PhaseLogger.begin("Phase 6  Weight calculation", gpuWeight);
            // weightClusterTrees = completed trees (+ any consensus exemplar trees from
            //                      the emission bridge) for cluster exemplar position lookups
            // originalTrees      = original trees (for rooted gene-tree triplet scoring)
            WeightTable weightTable = new WeightTable(dpTable, partTable, clusterTable, weightClusterTrees, originalTrees);
            PhaseLogger.end("Phase 6  Weight calculation", t6, gpuWeight);

            if (cfg.isVerifyWeights()) {
                Phase6Verifier.dump(trees, registry, clusterTable, dpTable, weightTable, cfg.getOutputFile());
                return;
            }

            // ── Phase 7: Inference DP + tree reconstruction ───────────────────
            long t7 = PhaseLogger.begin("Phase 7  Inference", false);
            Inference inference = new Inference();
            String speciesTree = inference.run(
                dpTable, weightTable, clusterTable, weightClusterTrees,
                registry, hasher);
            finalTripletScore = inference.getLastTripletScore();
            PhaseLogger.end("Phase 7  Inference", t7, false);

            // Persist the inferred topology before the independent final scoring
            // pass.  If scoring encounters an environmental/runtime failure, the
            // expensive inference result remains available for a score-only retry.
            if (cfg.getOutputFile() != null) {
                try (java.io.PrintStream out = new java.io.PrintStream(
                        new java.io.FileOutputStream(cfg.getOutputFile()))) {
                    out.println(speciesTree);
                }
                Logging.info("Species tree written to %s", cfg.getOutputFile());
            } else {
                System.out.println(speciesTree);
            }

            boolean recomputeFinalScore = inferenceInputHasPolytomy
                && !cfg.isKeepPolytomyDuringInference();
            if (recomputeFinalScore) {
                // The inference state can be enormous.  Drop all topology/weight
                // references before re-reading unresolved gene trees so the final
                // fixed-tree pass does not require both full working sets at once.
                weightTable = null;
                partTable = null;
                dpTable = null;
                clusterTable = null;
                weightClusterTrees = null;
                trees = null;
                originalTrees = null;
                hasher = null;
                System.gc();

                long t8 = PhaseLogger.begin(
                    "Phase 8  Final triplet scoring against unresolved input", false);
                FinalScoringInput scoringInput = parseFinalScoringInput(cfg, registry);
                Tree scoredSpeciesTree = TreeParser.parseSpeciesTreeNewick(
                    speciesTree, scoringInput.registry());
                TaxonHasher scoringHasher = new TaxonHasher(
                    scoringInput.registry().size(), cfg.getNumHashSeeds(), cfg.getBaseSeed());
                PrefixHashArrays scoringPref = new PrefixHashArrays(
                    scoringInput.trees(), scoringHasher);
                finalTripletScore = calculateFixedTreeScore(
                    cfg, scoringInput.registry(), scoringInput.trees(), scoringPref,
                    scoringHasher, scoredSpeciesTree, "Final score", "Final scoring");
                PhaseLogger.end(
                    "Phase 8  Final triplet scoring against unresolved input", t8, false);
                Logging.info("Final triplet score = %s  [recomputed from native input "
                    + "polytomies]", finalTripletScore);
            } else {
                String reason = cfg.isKeepPolytomyDuringInference()
                    ? "inference preserved native input polytomies"
                    : "no unresolved input polytomies detected";
                Logging.info("Final triplet score = %s  [reused inference DP score; %s]",
                    finalTripletScore, reason);
            }

            analysisCompleted = true;

        } finally {
            Threading.shutdown();
            long ms = (System.nanoTime() - t0) / 1_000_000;
            Logging.info("Total time: %d ms", ms);
            if (analysisCompleted) {
                PhaseLogger.printRunSummary(finalTripletScore, ms,
                    cfg.getComputeMode() == Config.ComputeMode.GPU);
            }
        }
    }

    private static boolean parseArgs(String[] args, Config cfg) {
        for (int i = 0; i < args.length; i++) {
            if (args[i].startsWith("--log-file=")) {
                String path = args[i].substring("--log-file=".length());
                if (path.isBlank()) return false;
                cfg.setLogFile(path);
                continue;
            }
            switch (args[i]) {
                case "-i","--input"    -> { if (++i>=args.length) return false; cfg.setInputFile(args[i]); }
                case "-o","--output"   -> { if (++i>=args.length) return false; cfg.setOutputFile(args[i]); }
                case "--log-file"      -> { if (++i>=args.length) return false; cfg.setLogFile(args[i]); }
                case "-c", "--score", "--species-tree", "--score-species-tree" -> {
                    if (++i>=args.length) return false;
                    cfg.setScoreSpeciesTreeFile(args[i]);
                }
                case "--taxa-file", "--species-list", "--species-list-file" -> {
                    if (++i>=args.length) return false;
                    cfg.setTaxaFile(args[i]);
                }
                case "--extract-taxa" -> cfg.setExtractTaxa(true);
                case "--taxa-set", "--taxa-operation" -> {
                    if (++i>=args.length) return false;
                    if (args[i].equalsIgnoreCase("union"))
                        cfg.setTaxaSetMode(Config.TaxaSetMode.UNION);
                    else if (args[i].equalsIgnoreCase("intersection"))
                        cfg.setTaxaSetMode(Config.TaxaSetMode.INTERSECTION);
                    else {
                        System.err.println("Unknown --taxa-set: " + args[i]
                            + " (expected: union | intersection)");
                        return false;
                    }
                }
                case "-t", "-T", "--threads", "--num-threads" -> {
                    if (++i>=args.length) return false;
                    cfg.setThreadCount(Integer.parseInt(args[i]));
                }
                case "--auto"          -> cfg.setComputeMode(Config.ComputeMode.AUTO);
                case "--cpu"           -> cfg.setComputeMode(Config.ComputeMode.CPU);
                case "--gpu"           -> cfg.setComputeMode(Config.ComputeMode.GPU);
                case "--gpu-strict"    -> { cfg.setComputeMode(Config.ComputeMode.GPU); cfg.setGpuStrict(true); }
                case "--search-space"  -> {
                    if (++i >= args.length) return false;
                    try { CliPresets.applySearchSpace(args[i], cfg); }
                    catch (IllegalArgumentException e) {
                        System.err.println(e.getMessage());
                        return false;
                    }
                }
                case "--search-mode"   -> {
                    if (++i >= args.length) return false;
                    if (args[i].equalsIgnoreCase("full"))
                        cfg.setSearchMode(Config.SearchMode.FULL);
                    else if (args[i].equalsIgnoreCase("local"))
                        cfg.setSearchMode(Config.SearchMode.LOCAL);
                    else {
                        System.err.println("Unknown --search-mode: " + args[i]
                            + " (expected: local | full)");
                        return false;
                    }
                }
                case "--anchor-outgroup", "--anchor" ->
                    cfg.setAnchorOutgroup(true);
                case "--no-anchor-outgroup", "--no-anchor" ->
                    cfg.setAnchorOutgroup(false);
                case "--anchor-taxon" -> {
                    if (++i >= args.length) return false;
                    try { cfg.setAnchorTaxon(Integer.parseInt(args[i])); }
                    catch (NumberFormatException e) {
                        System.err.println("Invalid --anchor-taxon (expected an integer taxon id): " + args[i]);
                        return false;
                    }
                }
                case "--no-prune-search-space", "--no-prune-unreachable" ->
                    cfg.setPruneUnreachableSplits(false);
                case "--prune-search-space", "--prune-unreachable" ->
                    cfg.setPruneUnreachableSplits(true);
                case "--intersection-method", "--im", "--weight-intersection-method" -> {
                    if (++i >= args.length) return false;
                    try { CliPresets.applyIntersectionMethod(args[i], cfg); }
                    catch (IllegalArgumentException e) {
                        System.err.println(e.getMessage());
                        return false;
                    }
                }
                case "--large-n-score-type", "--large-score-type" -> {
                    if (++i >= args.length) return false;
                    String t = args[i].toLowerCase().replace('_', '-');
                    switch (t) {
                        case "int128", "i128", "exact" ->
                            cfg.setLargeScoreType(Config.LargeScoreType.INT128);
                        case "double", "fp64", "float" ->
                            cfg.setLargeScoreType(Config.LargeScoreType.DOUBLE);
                        default -> {
                            System.err.println("Unknown --large-n-score-type: " + args[i]
                                + "  (expected: int128 | double)");
                            return false;
                        }
                    }
                }
                case "-vv"             -> cfg.setVerbosity(Logging.DEBUG);
                case "-vvv"            -> cfg.setVerbosity(Logging.TRACE);
                case "-q","--quiet"    -> cfg.setVerbosity(Logging.QUIET);
                case "-m","--seeds"    -> { if (++i>=args.length) return false; cfg.setNumHashSeeds(Integer.parseInt(args[i])); }
                case "--rooted"        -> cfg.setTreatAsUnrooted(false);
                case "--unrooted"      -> cfg.setTreatAsUnrooted(true);
                case "--keep-polytomy-during-inference" ->
                    cfg.setKeepPolytomyDuringInference(true);
                case "--no-gpu-batch"    -> cfg.setGpuBatch(false);
                case "--gpu-batch-size"  -> { if (++i>=args.length) return false; cfg.setGpuBatchSize(Integer.parseInt(args[i])); }
                case "--gpu-batches"     -> { if (++i>=args.length) return false; cfg.setGpuNumBatches(Integer.parseInt(args[i])); }
                case "--gpu-vram-control-factor"   -> { if (++i>=args.length) return false; cfg.setGpuVramControlFactor(Double.parseDouble(args[i])); }
                case "--gpu-vram-occupancy-factor" -> { if (++i>=args.length) return false; cfg.setGpuVramFraction(Double.parseDouble(args[i])); }
                case "--gpu-treewalk-vram-cap-mb"  -> { if (++i>=args.length) return false; cfg.setGpuTreeWalkVramCapMiB(Integer.parseInt(args[i])); }
                case "--gpu-progress-interval"     -> { if (++i>=args.length) return false; cfg.setGpuProgressIntervalSec(Double.parseDouble(args[i])); }
                case "--gpu-dp-state-space-construction-output-cap" -> { if (++i>=args.length) return false; cfg.setGpuDpStateSpaceConstructionOutputCap(args[i]); }
                case "--gpu-dp-state-space-progress-time-interval"  -> { if (++i>=args.length) return false; cfg.setGpuDpProgressInterval(Double.parseDouble(args[i])); }
                case "--gpu-dp-state-space-progress-max-steps"      -> { if (++i>=args.length) return false; cfg.setGpuDpProgressMaxSteps(Integer.parseInt(args[i])); }
                case "--verify-parse"  -> cfg.setVerifyParse(true);
                case "--verify-hash"      -> cfg.setVerifyHash(true);
                case "--verify-clusters"    -> cfg.setVerifyClusters(true);
                case "--verify-partitions" -> cfg.setVerifyPartitions(true);
                case "--verify-dp"         -> cfg.setVerifyDPSpace(true);
                case "--verify-weights"    -> cfg.setVerifyWeights(true);
                case "--verify-distance-matrix"    -> cfg.setVerifyDistanceMatrix(true);
                case "--verify-similarity-matrix"  -> cfg.setVerifySimilarityMatrix(true);
                case "--verify-upgma"              -> cfg.setVerifyUpgma(true);
                case "--verify-greedy-consensus"   -> cfg.setVerifyGreedyConsensus(true);
                case "--consensus-experimental"    -> cfg.setConsensusExperimental(true);
                case "--stepb-restriction"         -> {
                    if (++i >= args.length) return false;
                    String v = args[i].toLowerCase();
                    if (v.equals("dlogd") || v.equals("fast"))      cfg.setStepBFastRestriction(true);
                    else if (v.equals("n") || v.equals("full"))     cfg.setStepBFastRestriction(false);
                    else { System.err.println("--stepb-restriction expects dlogd|n"); return false; }
                }
                case "--stepb-quadratic-nn-balls"          -> cfg.setStepBQuadraticNnBalls(true);
                case "--stepb-random-leftover-resolution"  -> cfg.setStepBRandomLeftoverResolution(true);
                case "--stepb-process-large-polytomies"    -> cfg.setStepBProcessLargePolytomies(true);
                case "--resolve-input-gene-tree-polytomies" -> cfg.setResolveInputGeneTreePolytomies(true);
                case "--autocomplete-incomplete-gene-trees" -> cfg.setAutoCompleteIncompleteTrees(true);
                case "--completion-method" -> {
                    if (++i >= args.length) return false;
                    if (args[i].equalsIgnoreCase("distance"))
                        cfg.setCompletionMethod(Config.CompletionMethod.DISTANCE);
                    else if (args[i].equalsIgnoreCase("similarity"))
                        cfg.setCompletionMethod(Config.CompletionMethod.SIMILARITY);
                    else {
                        System.err.println("Unknown --completion-method: " + args[i]
                            + " (expected: similarity | distance)");
                        return false;
                    }
                }
                case "--dump-clusters"         -> { if (++i>=args.length) return false; cfg.setDumpClustersFile(args[i]); }
                case "--dump-completed-gene-trees" -> { if (++i>=args.length) return false; cfg.setDumpCompletedTreesFile(args[i]); }
                case "--gpu-dist-tile-size" -> { if (++i>=args.length) return false; cfg.setGpuDistTileSizeB(Integer.parseInt(args[i])); }
                case "--gpu-sim-vram-cap-mb" -> { if (++i>=args.length) return false; cfg.setGpuSimilarityVramCapMiB(Integer.parseInt(args[i])); }
                case "--diagnose"       -> cfg.setDiagnose(true);
                case "-v", "--version"  -> { Banner.printVersion(); System.exit(0); }
                case "-h","--help"     -> { printUsage(); System.exit(0); }
                default -> { System.err.println("Unknown arg: " + args[i]); return false; }
            }
        }
        if (cfg.getInputFile() == null && !cfg.isDiagnose()) return false;
        if (cfg.isExtractTaxa() && cfg.isScoreOnly()) {
            System.err.println("--extract-taxa cannot be combined with --score-species-tree");
            return false;
        }
        if (cfg.isExtractTaxa() && cfg.getTaxaFile() != null) {
            System.err.println("--extract-taxa cannot be combined with --taxa-file");
            return false;
        }
        if (cfg.isDiagnose() && cfg.getTaxaFile() != null) {
            System.err.println("--taxa-file requires a gene-tree analysis input");
            return false;
        }
        return true;
    }

    /** Resolve AUTO/GPU requests using the bundled CUDA backend itself. */
    private static void resolveComputeMode(Config cfg) {
        Config.ComputeMode requested = cfg.getRequestedComputeMode();
        if (requested == Config.ComputeMode.CPU) {
            cfg.resolveComputeMode(Config.ComputeMode.CPU, "explicit --cpu");
            return;
        }

        GPUWeightCalculator.Probe probe = GPUWeightCalculator.probe();
        if (probe.cudaAvailable()) {
            String selection = requested == Config.ComputeMode.AUTO
                ? "auto-selected CUDA device " + probe.deviceName()
                : "CUDA device " + probe.deviceName();
            cfg.resolveComputeMode(Config.ComputeMode.GPU, selection);
            return;
        }

        String reason = probe.detail();
        if (requested == Config.ComputeMode.GPU && cfg.isGpuStrict()) {
            throw new IllegalStateException("--gpu-strict requested, but CUDA is unavailable: " + reason);
        }
        cfg.resolveComputeMode(Config.ComputeMode.CPU,
            (requested == Config.ComputeMode.AUTO ? "automatic CPU fallback: " : "GPU unavailable; CPU fallback: ")
            + reason);
        if (requested == Config.ComputeMode.GPU) {
            Logging.warn("GPU requested but unavailable; continuing on CPU. %s", reason);
        } else {
            Logging.info("CUDA unavailable; using CPU. %s", reason);
        }
    }

    /**
     * Print distance matrix to stdout in a machine-parseable format:
     *   DISTANCE_MATRIX
     *   n=<count>
     *   taxa=name0,name1,...
     *   row0=d00,d01,...
     *   row1=d10,d11,...
     *   ...
     */
    private static void dumpDistanceMatrix(stelarx.completion.DistanceMatrix dm,
                                            TaxonRegistry registry) {
        int n = dm.n;
        StringBuilder taxa = new StringBuilder("taxa=");
        for (int i = 0; i < n; i++) {
            if (i > 0) taxa.append(',');
            taxa.append(registry.getName(i));
        }
        System.out.println("DISTANCE_MATRIX");
        System.out.println("n=" + n);
        System.out.println(taxa);
        for (int i = 0; i < n; i++) {
            StringBuilder row = new StringBuilder("row").append(i).append('=');
            for (int j = 0; j < n; j++) {
                if (j > 0) row.append(',');
                double d = dm.dist[i * n + j];
                if (d == Double.MAX_VALUE) row.append("inf");
                else row.append(String.format("%.6f", d));
            }
            System.out.println(row);
        }
    }

    /**
     * Print similarity matrix to stdout in the same machine-parseable format as
     * dumpDistanceMatrix — just with SIMILARITY_MATRIX header and sim_rowN keys.
     */
    private static void dumpSimilarityMatrix(stelarx.completion.SimilarityMatrix sm,
                                              TaxonRegistry registry) {
        int n = sm.n;
        StringBuilder taxa = new StringBuilder("taxa=");
        for (int i = 0; i < n; i++) {
            if (i > 0) taxa.append(',');
            taxa.append(registry.getName(i));
        }
        System.out.println("SIMILARITY_MATRIX");
        System.out.println("n=" + n);
        System.out.println(taxa);
        for (int i = 0; i < n; i++) {
            StringBuilder row = new StringBuilder("sim_row").append(i).append('=');
            for (int j = 0; j < n; j++) {
                if (j > 0) row.append(',');
                row.append(String.format("%.8f", sm.getSim(i, j)));
            }
            System.out.println(row);
        }
    }

    /**
     * Dump UPGMA bipartitions to stdout.
     *
     * Format:
     *   UPGMA_BIPARTITIONS
     *   n=<count>
     *   taxa=name0,name1,...
     *   bipartition=nameA,nameB,...   (one line per internal non-root node, names sorted)
     *
     * Each bipartition is the sorted set of taxon names in that subtree.
     * Root (all-taxa) is skipped automatically (rangeSize == n).
     */
    /** True iff any of the first {@code numGeneTrees} trees contains a polytomous node. */
    private static boolean anyPolytomous(java.util.List<Tree> trees, int numGeneTrees) {
        for (int g = 0; g < numGeneTrees && g < trees.size(); g++) {
            if (hasPolytomousNode(trees.get(g).root)) return true;
        }
        return false;
    }

    private static boolean hasPolytomousNode(stelarx.tree.TreeNode node) {
        if (node == null || node.isLeaf()) return false;
        if (node.isPolytomous()) return true;
        return hasPolytomousNode(node.left) || hasPolytomousNode(node.right);
    }

    private static void dumpUpgmaBipartitions(Tree upgmaTree, TaxonRegistry registry) {
        int n = registry.size();
        StringBuilder taxaLine = new StringBuilder("taxa=");
        for (int i = 0; i < n; i++) {
            if (i > 0) taxaLine.append(',');
            taxaLine.append(registry.getName(i));
        }
        System.out.println("UPGMA_BIPARTITIONS");
        System.out.println("n=" + n);
        System.out.println(taxaLine);

        // Post-order walk; print subtree leaf set for every internal non-root node
        dumpUpgmaNode(upgmaTree.root, upgmaTree, registry, n);
    }

    private static void dumpUpgmaNode(stelarx.tree.TreeNode node, Tree tree,
                                       TaxonRegistry registry, int n) {
        if (node.isLeaf()) return;
        dumpUpgmaNode(node.left,  tree, registry, n);
        dumpUpgmaNode(node.right, tree, registry, n);
        if (node.isRoot()) return;  // skip all-taxa bipartition

        // Collect and sort taxon names in [rangeStart, rangeEnd)
        java.util.List<String> names = new java.util.ArrayList<>(node.rangeSize());
        for (int pos = node.rangeStart; pos < node.rangeEnd; pos++) {
            names.add(registry.getName(tree.postorderArray[pos]));
        }
        java.util.Collections.sort(names);
        System.out.println("bipartition=" + String.join(",", names));
    }

    /**
     * Dump completed gene trees to a file, one Newick per line.
     * Ordering matches the original input gene tree order.
     */
    static void dumpCompletedTrees(List<Tree> trees, TaxonRegistry registry,
                                    String outFile) throws IOException {
        try (PrintStream out = new PrintStream(new FileOutputStream(outFile))) {
            for (Tree t : trees) {
                out.println(t.toNewick(registry));
            }
        }
        Logging.info("Completed gene trees written to %s (%d trees)", outFile, trees.size());
    }

    /**
     * Dump all clusters in ClusterTable to a file in canonical sorted format.
     * Each line: {A,B,C} with taxon names sorted alphabetically, lines sorted lexicographically.
     * This format matches the ASTRAL-MP --dump-clusters output for head-to-head comparison.
     */
    static void dumpClusters(ClusterTable clusterTable, List<Tree> trees,
                              TaxonRegistry registry, String outFile) throws IOException {
        List<String> lines = new ArrayList<>(clusterTable.size());
        for (ClusterTable.Entry e : clusterTable.entries()) {
            Cluster ex = e.exemplar;
            Tree tree = trees.get(ex.treeIndex);
            int[] arr = tree.postorderArray;
            List<String> taxa = new ArrayList<>(e.hash.size);
            if (!ex.complement) {
                for (int i = ex.left; i < ex.right; i++)
                    taxa.add(registry.getName(arr[i]));
            } else {
                for (int i = 0; i < tree.leafCount; i++) {
                    if (i >= ex.left && i < ex.right) continue;
                    taxa.add(registry.getName(arr[i]));
                }
            }
            Collections.sort(taxa);
            lines.add("{" + String.join(",", taxa) + "}");
        }
        Collections.sort(lines);
        try (PrintStream out = new PrintStream(new FileOutputStream(outFile))) {
            for (String line : lines) out.println(line);
        }
        Logging.info("Cluster dump written to %s (%d entries)", outFile, lines.size());
    }

    private static void printUsage() {
        Banner.printTitle(System.err);
        System.err.print("""
            Usage:
              stelarx -i <rooted_gene_trees.tre> [-o <species_tree.tre>] [options]
              stelarx -i <rooted_gene_trees.tre> --score-species-tree <rooted_species_tree.tre> [options]
              stelarx -i <trees.tre> --extract-taxa [-o <taxa.txt>]
              stelarx --diagnose

            General:
              -i, --input FILE                 Input gene trees (one Newick tree per line)
              -o, --output FILE                Output species tree (stdout when omitted)
              --log-file FILE                  Save run messages to FILE (progress remains terminal-only)
              -c, --score-species-tree FILE    Score one supplied species tree and exit
              --taxa-file FILE                 Restrict inference to these gene-tree taxa,
                                                 or both inputs in score-only mode
                                                 (one name per non-empty line)
              --extract-taxa                   Write input taxa, one name per line, and exit
              --taxa-set union|intersection    Multi-tree extraction operation (default: union)
              -t, -T, --threads, --num-threads N
                                                 CPU worker threads (default: available cores)
              --auto                           Automatically use CUDA or fall back to CPU (default)
              --cpu                            Force CPU execution; do not require GPU libraries
              --gpu                            Prefer CUDA; warn and fall back to CPU if unavailable
              --gpu-strict                     Require CUDA; fail before reading input if unavailable
              --diagnose                       Print installation/hardware diagnostics and exit
              -v, --version                    Print version and exit
              -h, --help                       Show this help and exit
              -q, --quiet                      Quiet logging
              -vv | -vvv                       Debug or trace logging

            Search and scoring:
              --search-space S1..S3            Friendly search-space preset (default: S1)
              --intersection-method I1..I4     Friendly scoring method (default: I2)
              --im I1..I4                      Short form of --intersection-method
              --search-mode local|full         Legacy/advanced DP search control
              --weight-intersection-method M   Legacy alias; named values remain supported
              --large-n-score-type T            int128 (exact) | double
              --no-prune-search-space           Disable reachability pruning
              --rooted                          Rooted input treatment (required and default)
              --keep-polytomy-during-inference   Preserve input polytomies while inferring;
                                                  final scoring always preserves them
                                                 (default: deterministic binary refinement)
              -m, --seeds N                     Number of cluster-hash seeds

            Incomplete trees and enrichment:
              --autocomplete-incomplete-gene-trees
              --completion-method similarity|distance
              --consensus-experimental
              --stepb-restriction dlogd|n
              --stepb-quadratic-nn-balls
              --stepb-random-leftover-resolution
              --stepb-process-large-polytomies
              --resolve-input-gene-tree-polytomies

            GPU memory and batching:
              --no-gpu-batch
              --gpu-batch-size N
              --gpu-batches N
              --gpu-vram-occupancy-factor F
              --gpu-vram-control-factor F
              --gpu-treewalk-vram-cap-mb MiB    Default: 512
              --gpu-sim-vram-cap-mb MiB         Default: 512; large-N auto may raise within free VRAM
              --gpu-dist-tile-size N
              --gpu-dp-state-space-construction-output-cap SIZE
              --gpu-progress-interval SECONDS

            Verification/debug outputs:
              --verify-parse | --verify-hash | --verify-clusters
              --verify-partitions | --verify-dp | --verify-weights
              --verify-distance-matrix | --verify-similarity-matrix
              --verify-upgma | --verify-greedy-consensus
              --dump-clusters FILE
              --dump-completed-gene-trees FILE

            The standalone distribution includes its own Java runtime. CUDA is optional;
            a missing/incompatible NVIDIA driver always has a CPU fallback unless
            --gpu-strict is used.
            """);
    }

    private static void runTaxaExtraction(Config cfg) throws IOException {
        if (cfg.getOutputFile() != null
                && sameNormalizedPath(cfg.getInputFile(), cfg.getOutputFile())) {
            throw new IllegalArgumentException(
                "Taxa output file must differ from the input tree file");
        }
        int count = TreeTaxa.writeExtracted(
            cfg.getInputFile(), cfg.getOutputFile(), cfg.getTaxaSetMode());
        String destination = cfg.getOutputFile() == null ? "stdout" : cfg.getOutputFile();
        Logging.info("Extracted %d taxa by %s across the input trees; written to %s",
            count, cfg.getTaxaSetMode().name().toLowerCase(), destination);
    }

    /**
     * Build the inference universe from the allow-list and parse only its
     * induced gene-tree leaves. The initial token scan determines which listed
     * names actually occur; absent names are intentionally never registered.
     */
    private static RestrictedInferenceInput parseTaxonRestrictedInferenceInput(
            Config cfg) throws IOException {
        if (cfg.getOutputFile() != null
                && sameNormalizedPath(cfg.getTaxaFile(), cfg.getOutputFile())) {
            throw new IllegalArgumentException(
                "Species-tree output file must differ from --taxa-file");
        }

        TreeTaxa.TaxaList taxaList = TreeTaxa.readTaxaList(cfg.getTaxaFile());
        java.util.LinkedHashSet<String> requested = taxaList.names();
        TreeTaxa.SelectionScan coverage =
            TreeTaxa.scanSelection(cfg.getInputFile(), requested);
        java.util.LinkedHashSet<String> effective = new java.util.LinkedHashSet<>();
        for (String name : requested) {
            if (coverage.selectedUnion().contains(name)) effective.add(name);
        }

        int requestedCount = requested.size();
        int presentCount = coverage.selectedUnion().size();
        int absent = requestedCount - presentCount;

        Logging.info("Taxon filter report:");
        Logging.info("  Taxa file: %d unique name(s)%s", requestedCount,
            taxaList.duplicateLines() == 0 ? ""
                : String.format(" (%d duplicate line(s) ignored)",
                    taxaList.duplicateLines()));
        Logging.info("  Gene trees: %d tree(s); listed taxa missing per tree "
                + "mean=%.2f (%.3f%%), min=%d, max=%d",
            coverage.treeCount(), coverage.meanMissing(),
            percentage(coverage.meanMissing(), requestedCount),
            coverage.minMissing(), coverage.maxMissing());
        Logging.info("  Listed taxa absent from every gene tree: %d (%.3f%%)",
            absent, percentage(absent, requestedCount));
        Logging.info("  Ignored unlisted leaf occurrences: %d",
            coverage.ignoredLeafOccurrences());
        Logging.info("  Effective inference universe: %d taxa", effective.size());
        if (absent > 0) {
            Logging.warn("The inference universe excludes listed taxa absent from every "
                + "gene tree; no placement is invented for them");
        }
        if (effective.size() < 3) {
            throw new IllegalArgumentException("Fewer than three taxa from --taxa-file "
                + "occur in the gene-tree union");
        }

        TaxonRegistry registry = new TaxonRegistry();
        for (String name : effective) registry.register(name);
        registry.lock();
        if (cfg.isAnchorOutgroup() && cfg.getAnchorTaxon() >= registry.size()) {
            throw new IllegalArgumentException("--anchor-taxon " + cfg.getAnchorTaxon()
                + " is outside the filtered taxon range [0,"
                + (registry.size() - 1) + "]");
        }

        TreeParser.RestrictedGeneTrees parsed = TreeParser.parseRestrictedGeneTrees(
            cfg.getInputFile(), registry, cfg.isKeepPolytomyDuringInference());
        Logging.info("Taxon restriction retained %d/%d induced gene tree(s); discarded "
                + "%d with fewer than two selected taxa",
            parsed.trees().size(), parsed.sourceTreeCount(), parsed.droppedTreeCount());
        Logging.info("Mode: INFERENCE with taxa file (outside and globally absent taxa ignored)");
        return new RestrictedInferenceInput(registry, parsed.trees(),
            parsed.detectedPolytomyNodeCount() > 0);
    }

    private record RestrictedInferenceInput(TaxonRegistry registry, List<Tree> trees,
                                            boolean hasPolytomy) {}

    /**
     * Re-read exactly the scoring universe while retaining native input
     * multifurcations.  The inference registry is reused for a taxa-file run so
     * globally absent listed names remain excluded exactly as they were during
     * inference.
     */
    private static FinalScoringInput parseFinalScoringInput(
            Config cfg, TaxonRegistry inferenceRegistry) throws IOException {
        if (cfg.getTaxaFile() == null) {
            TaxonRegistry scoringRegistry = new TaxonRegistry();
            List<Tree> scoringTrees = TreeParser.parseGeneTrees(
                cfg.getInputFile(), scoringRegistry, true);
            return new FinalScoringInput(scoringRegistry, scoringTrees);
        }

        TreeParser.RestrictedGeneTrees parsed = TreeParser.parseRestrictedGeneTrees(
            cfg.getInputFile(), inferenceRegistry, true);
        Logging.info("Final scoring retained %d/%d taxa-restricted gene tree(s)",
            parsed.trees().size(), parsed.sourceTreeCount());
        return new FinalScoringInput(inferenceRegistry, parsed.trees());
    }

    private record FinalScoringInput(TaxonRegistry registry, List<Tree> trees) {}

    /**
     * Opt-in fixed-tree scoring on a named taxon subset.  This path is kept
     * separate from ordinary inference and strict score-only mode so neither
     * incurs filtering branches, extra scans, or changed taxon semantics.
     */
    private static String runTaxonRestrictedScoreOnly(Config cfg) throws IOException {
        if (cfg.getOutputFile() != null
                && sameNormalizedPath(cfg.getTaxaFile(), cfg.getOutputFile())) {
            throw new IllegalArgumentException(
                "Score output file must differ from --taxa-file");
        }
        TreeTaxa.TaxaList taxaList = TreeTaxa.readTaxaList(cfg.getTaxaFile());
        java.util.LinkedHashSet<String> requested = taxaList.names();
        TreeTaxa.Scan geneCoverage = TreeTaxa.scan(cfg.getInputFile(), requested);
        TreeTaxa.Scan speciesCoverage = TreeTaxa.scan(
            cfg.getScoreSpeciesTreeFile(), requested);
        if (speciesCoverage.treeCount() != 1) {
            throw new IllegalArgumentException(
                "Species tree file must contain exactly one Newick tree: "
                + cfg.getScoreSpeciesTreeFile());
        }

        java.util.LinkedHashSet<String> effective = new java.util.LinkedHashSet<>();
        for (String name : requested) {
            if (geneCoverage.union().contains(name)
                    && speciesCoverage.union().contains(name)) {
                effective.add(name);
            }
        }

        int requestedCount = requested.size();
        int absentGeneUnion = requestedCount
            - intersectionSize(requested, geneCoverage.union());
        int absentSpecies = requestedCount
            - intersectionSize(requested, speciesCoverage.union());
        int ignoredGeneTaxa = geneCoverage.union().size()
            - intersectionSize(geneCoverage.union(), requested);
        int ignoredSpeciesTaxa = speciesCoverage.union().size()
            - intersectionSize(speciesCoverage.union(), requested);

        Logging.info("Taxon filter report:");
        Logging.info("  Taxa file: %d unique name(s)%s", requestedCount,
            taxaList.duplicateLines() == 0 ? ""
                : String.format(" (%d duplicate line(s) ignored)", taxaList.duplicateLines()));
        Logging.info("  Gene trees: %d tree(s), %d union taxa; listed taxa missing per tree "
                + "mean=%.2f (%.3f%%), min=%d, max=%d",
            geneCoverage.treeCount(), geneCoverage.union().size(),
            geneCoverage.meanMissing(), percentage(geneCoverage.meanMissing(), requestedCount),
            geneCoverage.minMissing(), geneCoverage.maxMissing());
        Logging.info("  Listed taxa absent from every gene tree: %d (%.3f%%)",
            absentGeneUnion, percentage(absentGeneUnion, requestedCount));
        Logging.info("  Species tree: %d listed taxa missing (%.3f%%)",
            absentSpecies, percentage(absentSpecies, requestedCount));
        Logging.info("  Ignored outside taxa: gene-tree union=%d, species tree=%d",
            ignoredGeneTaxa, ignoredSpeciesTaxa);
        Logging.info("  Effective common scoring universe: %d taxa", effective.size());
        if (absentGeneUnion > 0 || absentSpecies > 0) {
            Logging.warn("The effective scoring universe excludes listed taxa absent from "
                + "the gene-tree union or species tree; no placement is invented for them");
        }
        if (effective.size() < 3) {
            throw new IllegalArgumentException("Fewer than three taxa from --taxa-file "
                + "occur in both the gene-tree union and species tree");
        }

        TaxonRegistry targetRegistry = new TaxonRegistry();
        for (String name : effective) targetRegistry.register(name);
        targetRegistry.lock();
        if (cfg.isAnchorOutgroup() && cfg.getAnchorTaxon() >= targetRegistry.size()) {
            throw new IllegalArgumentException("--anchor-taxon " + cfg.getAnchorTaxon()
                + " is outside the filtered taxon range [0,"
                + (targetRegistry.size() - 1) + "]");
        }

        long tg = PhaseLogger.begin("Score filter  Parse and restrict gene trees", false);
        TaxonRegistry sourceRegistry = new TaxonRegistry();
        List<Tree> sourceTrees = TreeParser.parseGeneTrees(
            cfg.getInputFile(), sourceRegistry, true);
        TreeRestrictor.GeneResult restrictedGenes = TreeRestrictor.restrictGeneTrees(
            sourceTrees, sourceRegistry, targetRegistry);
        List<Tree> geneTrees = restrictedGenes.trees();
        Logging.info("Taxon restriction retained %d/%d gene tree(s); dropped %d with "
                + "fewer than three selected taxa (zero rooted-triplet contribution)",
            geneTrees.size(), sourceTrees.size(), restrictedGenes.droppedTreeCount());
        PhaseLogger.end("Score filter  Parse and restrict gene trees", tg, false);

        long ts = PhaseLogger.begin("Score filter  Parse and restrict species tree", false);
        TreeParser.StandaloneTree sourceSpecies = TreeParser.parseStandaloneTree(
            cfg.getScoreSpeciesTreeFile());
        Tree speciesTree = TreeRestrictor.restrictSpeciesTree(
            sourceSpecies.tree(), sourceSpecies.registry(), targetRegistry);
        Logging.info("Restricted supplied species tree: %d leaves", speciesTree.leafCount);
        PhaseLogger.end("Score filter  Parse and restrict species tree", ts, false);

        long th = PhaseLogger.begin("Score filter  Taxon hashing", false);
        TaxonHasher hasher = new TaxonHasher(targetRegistry.size(),
            cfg.getNumHashSeeds(), cfg.getBaseSeed());
        PrefixHashArrays genePref = new PrefixHashArrays(geneTrees, hasher);
        PhaseLogger.end("Score filter  Taxon hashing", th, false);

        Logging.info("Mode: SCORE-ONLY with taxa file (score supplied induced species tree; "
            + "no species-tree inference)");
        return scorePreparedSpeciesTree(
            cfg, targetRegistry, geneTrees, genePref, hasher, speciesTree);
    }

    private static int intersectionSize(java.util.Set<String> a,
                                        java.util.Set<String> b) {
        java.util.Set<String> smaller = a.size() <= b.size() ? a : b;
        java.util.Set<String> larger = smaller == a ? b : a;
        int count = 0;
        for (String name : smaller) if (larger.contains(name)) count++;
        return count;
    }

    private static double percentage(double count, int total) {
        return total == 0 ? 0.0 : 100.0 * count / total;
    }

    private static boolean sameNormalizedPath(String first, String second) {
        return java.nio.file.Path.of(first).toAbsolutePath().normalize().equals(
            java.nio.file.Path.of(second).toAbsolutePath().normalize());
    }

    /** Reject read/write aliases before any analysis output can replace an input. */
    private static void validatePathCollisions(Config cfg) {
        String[][] reads = {
            {"input tree file", cfg.getInputFile()},
            {"species tree file", cfg.getScoreSpeciesTreeFile()},
            {"taxa file", cfg.getTaxaFile()}
        };
        String[][] writes = {
            {"output file", cfg.getOutputFile()},
            {"log file", cfg.getLogFile()},
            {"cluster dump", cfg.getDumpClustersFile()},
            {"completed-tree dump", cfg.getDumpCompletedTreesFile()}
        };
        for (String[] write : writes) {
            if (write[1] == null) continue;
            for (String[] read : reads) {
                if (read[1] != null && sameNormalizedPath(write[1], read[1])) {
                    throw new IllegalArgumentException(
                        write[0] + " must differ from the " + read[0]);
                }
            }
        }
        for (int i = 0; i < writes.length; i++) {
            if (writes[i][1] == null) continue;
            for (int j = i + 1; j < writes.length; j++) {
                if (writes[j][1] != null
                        && sameNormalizedPath(writes[i][1], writes[j][1])) {
                    throw new IllegalArgumentException(
                        writes[i][0] + " must differ from the " + writes[j][0]);
                }
            }
        }
    }

    private static String runScoreOnly(Config cfg, TaxonRegistry registry,
                                       List<Tree> geneTrees, PrefixHashArrays genePref,
                                       TaxonHasher hasher) throws IOException {
        Logging.info("Mode: SCORE-ONLY (score supplied species tree; no species-tree inference)");
        // The weight table dispatches GPU/CPU internally from the compute mode, so
        // both paths are supported here.  The GPU path falls back to CPU on its own
        // if it reports infeasible, so no forced downgrade is needed.

        long ts = PhaseLogger.begin("Score mode  Parse supplied species tree", false);
        Tree speciesTree = TreeParser.parseSpeciesTree(cfg.getScoreSpeciesTreeFile(), registry);
        PhaseLogger.end("Score mode  Parse supplied species tree", ts, false);

        return scorePreparedSpeciesTree(
            cfg, registry, geneTrees, genePref, hasher, speciesTree);
    }

    private static String scorePreparedSpeciesTree(Config cfg, TaxonRegistry registry,
                                                   List<Tree> geneTrees,
                                                   PrefixHashArrays genePref,
                                                   TaxonHasher hasher,
                                                   Tree speciesTree) throws IOException {
        String score = calculateFixedTreeScore(
            cfg, registry, geneTrees, genePref, hasher, speciesTree,
            "Score mode", "Score-only");
        String line = "TRIPLET_SCORE: " + score;
        System.out.println(line);
        if (cfg.getOutputFile() != null) {
            try (PrintStream out = new PrintStream(new FileOutputStream(cfg.getOutputFile()))) {
                out.println(line);
            }
            Logging.info("Triplet score written to %s", cfg.getOutputFile());
        }
        return score;
    }

    private static String calculateFixedTreeScore(Config cfg,
                                                   TaxonRegistry registry,
                                                   List<Tree> geneTrees,
                                                   PrefixHashArrays genePref,
                                                   TaxonHasher hasher,
                                                   Tree speciesTree,
                                                   String phasePrefix,
                                                   String scoreLogLabel) {
        List<Tree> speciesTrees = java.util.List.of(speciesTree);

        long tp = PhaseLogger.begin(phasePrefix + "  Species-tree hashing", false);
        PrefixHashArrays speciesPref = new PrefixHashArrays(speciesTrees, hasher);
        PhaseLogger.end(phasePrefix + "  Species-tree hashing", tp, false);

        long tc = PhaseLogger.begin(phasePrefix + "  Species-tree clusters", false);
        ClusterTable speciesClusters = new ClusterTable(speciesTrees, speciesPref, registry.size());
        PhaseLogger.end(phasePrefix + "  Species-tree clusters", tc, false);

        long tg = PhaseLogger.begin(phasePrefix + "  Gene-tree rooted child partitions", false);
        // Simple tree-walk consumes the original gene-tree topology directly; it
        // never reads PartitionTable.  Avoid materialising that very large table
        // in score-only mode (notably, the 9,524-taxon angiosperm data otherwise
        // needs tens of GiB merely to reach the weight kernel).
        PartitionTable genePartitions = null;
        if (cfg.getWeightIntersectionMethod()
                != Config.WeightIntersectionMethod.SIMPLE_TREE_WALK) {
            genePartitions = new PartitionTable(geneTrees, genePref);
        } else {
            Logging.info("%s tree-walk: skipped unused gene-tree PartitionTable",
                scoreLogLabel);
        }
        PhaseLogger.end(phasePrefix + "  Gene-tree rooted child partitions", tg, false);

        long td = PhaseLogger.begin(phasePrefix + "  Fixed-tree DP transitions", false);
        DPTable speciesDP = new DPTable(speciesTrees, speciesPref, speciesClusters);
        if (cfg.isAnchorOutgroup()) {
            speciesDP.applyAnchoredRoot(speciesClusters.getAnchorHash());
        }
        PhaseLogger.end(phasePrefix + "  Fixed-tree DP transitions", td, false);

        boolean gpuWeight = cfg.getComputeMode() == Config.ComputeMode.GPU
                            && GPUWeightCalculator.isLoaded();
        long tw = PhaseLogger.begin(phasePrefix + "  Weight calculation", gpuWeight);
        WeightTable weightTable = new WeightTable(speciesDP, genePartitions, speciesClusters,
                                                  speciesTrees, geneTrees);
        PhaseLogger.end(phasePrefix + "  Weight calculation", tw, gpuWeight);

        Inference scorer = new Inference();
        return scorer.scoreFixedTree(speciesDP, weightTable, scoreLogLabel);
    }
}
