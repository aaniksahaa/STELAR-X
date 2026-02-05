package core;

import java.util.*;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.atomic.AtomicLong;

import com.sun.jna.Memory;

import preprocessing.GeneTrees;
import tree.RangeBipartition;
import tree.Tree;
import tree.TreeNode;
import utils.Config;
import utils.Threading;
import core.WeightCalculator.WeightCalcLib;
import core.WeightCalculator.WeightCalcLib.CompactBipartition;

/**
 * Efficient triplet score calculator for species tree vs gene trees.
 * 
 * Uses the same range-based bipartition representation as gene trees.
 * The species tree is treated as an additional tree in the inverse index,
 * allowing O(min(|A|, |B|)) intersection calculations.
 * 
 * Supports CPU_SINGLE, CPU_PARALLEL, and GPU_PARALLEL modes.
 * No BitSets used - everything is range-based for efficiency.
 */
public class SpeciesTreeScorer {
    
    private final GeneTrees geneTrees;
    private final int numGeneTrees;
    private final int numTaxa;
    
    // Extended inverse index: includes gene trees + species tree
    // Species tree is at index = numGeneTrees
    private int[][] extendedInverseIndex;
    private int[][] extendedOrderings;
    private int[] treeSizes;
    private int[] missingOffsets;
    private int[] missingCounts;
    private int[] missingTaxaFlat;
    private List<RangeBipartition> geneTreeRangeList;
    private int[] geneTreeFrequencies;
    
    // Species tree bipartitions as RangeBipartitions
    private List<RangeBipartition> speciesBipartitions;
    
    public SpeciesTreeScorer(GeneTrees geneTrees) {
        this.geneTrees = geneTrees;
        this.numGeneTrees = geneTrees.geneTrees.size();
        this.numTaxa = geneTrees.realTaxaCount;
    }
    
    /**
     * Calculate triplet score between species tree and gene trees.
     * 
     * @param speciesTree The species tree to score
     * @return Total triplet score
     */
    public double calculateScore(Tree speciesTree) {
        System.out.println("\n=== Species Tree Score Calculation ===");
        System.out.println("Computation mode: " + Config.COMPUTATION_MODE);
        System.out.println("Score type: " + Config.SCORE_TYPE);
        System.out.println("Gene trees: " + numGeneTrees);
        System.out.println("Gene tree bipartitions: " + geneTrees.rangeBipartitions.size());
        System.out.println("Taxa: " + numTaxa);
        
        long startTime = System.currentTimeMillis();
        
        // Step 1: Build species tree ordering and add to extended inverse index
        buildExtendedInverseIndex(speciesTree);
        
        // Step 2: Extract bipartitions from species tree as RangeBipartitions
        extractSpeciesBipartitions(speciesTree);
        System.out.println("Species tree bipartitions: " + speciesBipartitions.size());

        // Step 3: Build stable gene-tree bipartition list (sorted by tree index)
        buildGeneTreeRangeList();
        
        // Step 4: Calculate total score based on computation mode
        double totalScore;
        
        switch (Config.COMPUTATION_MODE) {
            case GPU_PARALLEL:
                System.out.println("Using GPU_PARALLEL mode for score calculation");
                totalScore = calculateScoreGPU();
                break;
            case CPU_PARALLEL:
                System.out.println("Using CPU_PARALLEL mode for score calculation");
                totalScore = calculateScoreCPUParallel();
                break;
            case CPU_SINGLE:
            default:
                System.out.println("Using CPU_SINGLE mode for score calculation");
                totalScore = calculateScoreCPUSingle();
                break;
        }
        
        long endTime = System.currentTimeMillis();
        
        System.out.println("Time: " + (endTime - startTime) + " ms");
        System.out.println("===============================================\n");
        
        return totalScore;
    }
    
    /**
     * CPU single-threaded score calculation.
     */
    private double calculateScoreCPUSingle() {
        double totalScore = 0.0;
        
        for (RangeBipartition speciesBip : speciesBipartitions) {
            double bipScore = calculateBipartitionWeight(speciesBip);
            totalScore += bipScore;
        }
        
        return totalScore;
    }
    
    /**
     * CPU parallel score calculation.
     */
    private double calculateScoreCPUParallel() {
        int numBips = speciesBipartitions.size();
        int numThreads = Runtime.getRuntime().availableProcessors();
        
        Threading.startThreading(numThreads);
        
        int chunkSize = Math.max(1, (numBips + numThreads - 1) / numThreads);
        int actualThreads = Math.min(numThreads, (numBips + chunkSize - 1) / chunkSize);
        
        CountDownLatch latch = new CountDownLatch(actualThreads);
        AtomicLong[] partialScores = new AtomicLong[actualThreads];
        
        for (int i = 0; i < actualThreads; i++) {
            partialScores[i] = new AtomicLong(0);
            final int threadId = i;
            final int startIdx = i * chunkSize;
            final int endIdx = Math.min(startIdx + chunkSize, numBips);
            
            Threading.execute(() -> {
                try {
                    double localScore = 0.0;
                    for (int j = startIdx; j < endIdx; j++) {
                        localScore += calculateBipartitionWeight(speciesBipartitions.get(j));
                    }
                    partialScores[threadId].set(Double.doubleToLongBits(localScore));
                } finally {
                    latch.countDown();
                }
            });
        }
        
        try {
            latch.await();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new RuntimeException("Score calculation interrupted", e);
        } finally {
            Threading.shutdown();
        }
        
        double totalScore = 0.0;
        for (int i = 0; i < actualThreads; i++) {
            totalScore += Double.longBitsToDouble(partialScores[i].get());
        }
        
        return totalScore;
    }
    
    /**
     * GPU parallel score calculation.
     * Reuses the existing launchCompactWeightCalculation kernel.
     */
    private double calculateScoreGPU() {
        System.out.println("Preparing data for GPU kernel...");
        
        int numSpeciesBips = speciesBipartitions.size();
            int numGeneTreeBips = geneTreeRangeList.size();
        int totalTrees = numGeneTrees + 1;  // Gene trees + species tree
        
        System.out.println("GPU Parameters:");
        System.out.println("  Species bipartitions: " + numSpeciesBips);
        System.out.println("  Gene tree bipartitions: " + numGeneTreeBips);
        System.out.println("  Total trees (gene + species): " + totalTrees);
        System.out.println("  Taxa: " + numTaxa);
        
        try {
            // CRITICAL: Use Structure.toArray() for contiguous memory allocation
            // This is required by JNA for passing structure arrays to native code
            
            // Prepare species tree bipartitions using contiguous memory
            CompactBipartition speciesTemplate = new CompactBipartition();
            CompactBipartition[] speciesBipsArray = 
                (CompactBipartition[]) speciesTemplate.toArray(numSpeciesBips);
            
            for (int i = 0; i < numSpeciesBips; i++) {
                RangeBipartition rb = speciesBipartitions.get(i);
                speciesBipsArray[i].geneTreeIndex = rb.geneTreeIndex;  // = numGeneTrees for species tree
                speciesBipsArray[i].leftStart = rb.leftStart;
                speciesBipsArray[i].leftEnd = rb.leftEnd;
                speciesBipsArray[i].rightStart = rb.rightStart;
                speciesBipsArray[i].rightEnd = rb.rightEnd;
            }
            System.out.println("Species bipartitions array prepared (contiguous memory)");
            
            // Prepare gene tree bipartitions using contiguous memory
            CompactBipartition geneTreeTemplate = new CompactBipartition();
            CompactBipartition[] geneTreeBipsArray = 
                (CompactBipartition[]) geneTreeTemplate.toArray(numGeneTreeBips);
            int[] frequencies = new int[numGeneTreeBips];
            
            for (int idx = 0; idx < geneTreeRangeList.size(); idx++) {
                RangeBipartition rb = geneTreeRangeList.get(idx);
                geneTreeBipsArray[idx].geneTreeIndex = rb.geneTreeIndex;
                geneTreeBipsArray[idx].leftStart = rb.leftStart;
                geneTreeBipsArray[idx].leftEnd = rb.leftEnd;
                geneTreeBipsArray[idx].rightStart = rb.rightStart;
                geneTreeBipsArray[idx].rightEnd = rb.rightEnd;
                frequencies[idx] = geneTreeFrequencies[idx];
            }
            System.out.println("Gene tree bipartitions array prepared (contiguous memory)");
            
            // Flatten extended inverse index and orderings for GPU
            // GPU expects flat arrays: [tree0_taxa..., tree1_taxa..., ..., speciesTree_taxa...]
            long memSize = (long) totalTrees * numTaxa * 4;
            System.out.println("Allocating GPU memory: " + (memSize * 2 / (1024 * 1024)) + " MB for indices");
            
            Memory inverseIndexMem = new Memory(memSize);
            Memory orderingMem = new Memory(memSize);
            
            for (int t = 0; t < totalTrees; t++) {
                for (int taxon = 0; taxon < numTaxa; taxon++) {
                    long offset = ((long) t * numTaxa + taxon) * 4;
                    inverseIndexMem.setInt(offset, extendedInverseIndex[t][taxon]);
                }
                
                int[] ordering = extendedOrderings[t];
                for (int pos = 0; pos < numTaxa; pos++) {
                    long offset = ((long) t * numTaxa + pos) * 4;
                    if (ordering != null && pos < ordering.length) {
                        orderingMem.setInt(offset, ordering[pos]);
                    } else {
                        orderingMem.setInt(offset, -1);  // Sentinel for unused positions
                    }
                }
            }
            System.out.println("Inverse index and orderings flattened");
            Memory treeSizesMem = new Memory((long) totalTrees * 4);
            Memory missingOffsetsMem = new Memory((long) totalTrees * 4);
            Memory missingCountsMem = new Memory((long) totalTrees * 4);
            long missingBytes = Math.max(1, missingTaxaFlat.length) * 4L;
            Memory missingTaxaMem = new Memory(missingBytes);

            for (int t = 0; t < totalTrees; t++) {
                treeSizesMem.setInt((long) t * 4, treeSizes[t]);
                missingOffsetsMem.setInt((long) t * 4, missingOffsets[t]);
                missingCountsMem.setInt((long) t * 4, missingCounts[t]);
            }
            for (int i = 0; i < missingTaxaFlat.length; i++) {
                missingTaxaMem.setInt((long) i * 4, missingTaxaFlat[i]);
            }
            
            // Output array for weights
            double[] weights = new double[numSpeciesBips];
            
            System.out.println("Launching GPU kernel...");
            
            // Call the existing compact weight calculation kernel
            WeightCalcLib.INSTANCE.launchCompactWeightCalculation(
                speciesBipsArray,
                geneTreeBipsArray,
                frequencies,
                weights,
                inverseIndexMem,
                orderingMem,
                treeSizesMem,
                missingOffsetsMem,
                missingCountsMem,
                missingTaxaMem,
                missingTaxaFlat.length,
                numSpeciesBips,
                numGeneTreeBips,
                totalTrees,
                numTaxa,
                (Config.SCORE_TYPE == Config.ScoreType.QUARTET ? 1 : 0)
            );
            
            // Sum up weights to get total score
            double totalScore = 0.0;
            for (double w : weights) {
                totalScore += w;
            }
            
            System.out.println("GPU kernel completed successfully");
            return totalScore;
            
        } catch (Exception e) {
            System.err.println("GPU kernel failed: " + e.getMessage());
            e.printStackTrace();
            System.err.println("Falling back to CPU_PARALLEL mode...");
            return calculateScoreCPUParallel();
        }
    }
    
    /**
     * Build extended inverse index that includes the species tree.
     * Species tree gets index = numGeneTrees.
     */
    private void buildExtendedInverseIndex(Tree speciesTree) {
        int totalTrees = numGeneTrees + 1;
        extendedInverseIndex = new int[totalTrees][numTaxa];
        extendedOrderings = new int[totalTrees][];
        treeSizes = new int[totalTrees];
        
        // Copy gene tree inverse indices from existing InverseIndexManager
        // We'll build our own for gene trees too for simplicity
        for (int t = 0; t < numGeneTrees; t++) {
            Tree geneTree = geneTrees.geneTrees.get(t);
            List<Integer> ordering = new ArrayList<>();
            collectLeavesInOrder(geneTree.root, ordering);
            
            extendedOrderings[t] = ordering.stream().mapToInt(Integer::intValue).toArray();
            treeSizes[t] = extendedOrderings[t].length;
            
            // Initialize with -1 (sentinel for missing taxa)
            Arrays.fill(extendedInverseIndex[t], -1);
            
            // Fill in actual positions
            for (int pos = 0; pos < extendedOrderings[t].length; pos++) {
                int taxonId = extendedOrderings[t][pos];
                if (taxonId >= 0 && taxonId < numTaxa) {
                    extendedInverseIndex[t][taxonId] = pos;
                }
            }
        }
        
        // Add species tree at index numGeneTrees
        int speciesIdx = numGeneTrees;
        List<Integer> speciesOrdering = new ArrayList<>();
        collectLeavesInOrder(speciesTree.root, speciesOrdering);
        
        extendedOrderings[speciesIdx] = speciesOrdering.stream().mapToInt(Integer::intValue).toArray();
        treeSizes[speciesIdx] = extendedOrderings[speciesIdx].length;
        
        Arrays.fill(extendedInverseIndex[speciesIdx], -1);
        for (int pos = 0; pos < extendedOrderings[speciesIdx].length; pos++) {
            int taxonId = extendedOrderings[speciesIdx][pos];
            if (taxonId >= 0 && taxonId < numTaxa) {
                extendedInverseIndex[speciesIdx][taxonId] = pos;
            }
        }
        
        System.out.println("Extended inverse index built for " + totalTrees + " trees");

        buildMissingTaxaListsExtended();
    }
    
    /**
     * Collect leaves in left-to-right order (inorder traversal).
     */
    private void collectLeavesInOrder(TreeNode node, List<Integer> ordering) {
        if (node == null) return;
        
        if (node.isLeaf()) {
            if (node.taxon != null) {
                ordering.add(node.taxon.id);
            }
            return;
        }
        
        if (node.childs != null) {
            for (TreeNode child : node.childs) {
                collectLeavesInOrder(child, ordering);
            }
        }
    }

    private void buildMissingTaxaListsExtended() {
        int totalTrees = numGeneTrees + 1;
        missingCounts = new int[totalTrees];
        int totalMissing = 0;

        for (int treeIdx = 0; treeIdx < totalTrees; treeIdx++) {
            int missing = 0;
            for (int taxonId = 0; taxonId < numTaxa; taxonId++) {
                if (extendedInverseIndex[treeIdx][taxonId] == -1) {
                    missing++;
                }
            }
            missingCounts[treeIdx] = missing;
            totalMissing += missing;
        }

        missingOffsets = new int[totalTrees];
        int running = 0;
        for (int treeIdx = 0; treeIdx < totalTrees; treeIdx++) {
            missingOffsets[treeIdx] = running;
            running += missingCounts[treeIdx];
        }

        missingTaxaFlat = new int[totalMissing];
        for (int treeIdx = 0; treeIdx < totalTrees; treeIdx++) {
            int offset = missingOffsets[treeIdx];
            int idx = 0;
            for (int taxonId = 0; taxonId < numTaxa; taxonId++) {
                if (extendedInverseIndex[treeIdx][taxonId] == -1) {
                    missingTaxaFlat[offset + idx] = taxonId;
                    idx++;
                }
            }
        }
    }

    private void buildGeneTreeRangeList() {
        List<Map.Entry<RangeBipartition, Integer>> entries =
            new ArrayList<>(geneTrees.rangeBipartitions.entrySet());
        entries.sort(Comparator
            .comparingInt((Map.Entry<RangeBipartition, Integer> e) -> e.getKey().geneTreeIndex)
            .thenComparingInt(e -> e.getKey().leftStart)
            .thenComparingInt(e -> e.getKey().rightStart));

        geneTreeRangeList = new ArrayList<>(entries.size());
        geneTreeFrequencies = new int[entries.size()];
        for (int i = 0; i < entries.size(); i++) {
            geneTreeRangeList.add(entries.get(i).getKey());
            geneTreeFrequencies[i] = entries.get(i).getValue();
        }
    }
    
    /**
     * Extract bipartitions from species tree as RangeBipartitions.
     * Species tree has treeIndex = numGeneTrees.
     */
    private void extractSpeciesBipartitions(Tree speciesTree) {
        speciesBipartitions = new ArrayList<>();
        if (speciesTree.root != null) {
            extractBipartitionsRecursive(speciesTree.root, numGeneTrees);
        }
    }
    
    /**
     * Recursively extract bipartitions.
     * Returns [minPos, maxPos] in the species tree ordering.
     */
    private int[] extractBipartitionsRecursive(TreeNode node, int treeIndex) {
        if (node.isLeaf()) {
            if (node.taxon != null) {
                int pos = extendedInverseIndex[treeIndex][node.taxon.id];
                return new int[]{pos, pos};
            }
            return null;
        }
        
        if (node.childs == null || node.childs.size() < 2) {
            // Not a valid binary internal node
            int minPos = Integer.MAX_VALUE;
            int maxPos = Integer.MIN_VALUE;
            if (node.childs != null) {
                for (TreeNode child : node.childs) {
                    int[] range = extractBipartitionsRecursive(child, treeIndex);
                    if (range != null) {
                        minPos = Math.min(minPos, range[0]);
                        maxPos = Math.max(maxPos, range[1]);
                    }
                }
            }
            return (minPos == Integer.MAX_VALUE) ? null : new int[]{minPos, maxPos};
        }
        
        // Binary node - get ranges for left and right subtrees
        int[] leftRange = extractBipartitionsRecursive(node.childs.get(0), treeIndex);
        int[] rightRange = extractBipartitionsRecursive(node.childs.get(1), treeIndex);
        
        // Handle additional children (multifurcating)
        for (int i = 2; i < node.childs.size(); i++) {
            int[] extraRange = extractBipartitionsRecursive(node.childs.get(i), treeIndex);
            if (extraRange != null && rightRange != null) {
                rightRange[0] = Math.min(rightRange[0], extraRange[0]);
                rightRange[1] = Math.max(rightRange[1], extraRange[1]);
            } else if (extraRange != null) {
                rightRange = extraRange;
            }
        }
        
        if (leftRange != null && rightRange != null) {
            // Create RangeBipartition for this internal node
            // Note: ranges are [inclusive, inclusive], convert to [inclusive, exclusive)
            RangeBipartition bip = new RangeBipartition(
                treeIndex,
                leftRange[0], leftRange[1] + 1,   // left range
                rightRange[0], rightRange[1] + 1  // right range
            );
            speciesBipartitions.add(bip);
            
            // Return combined range
            return new int[]{
                Math.min(leftRange[0], rightRange[0]),
                Math.max(leftRange[1], rightRange[1])
            };
        }
        
        return (leftRange != null) ? leftRange : rightRange;
    }
    
    /**
     * Calculate weight for a species tree bipartition.
     * Sum over all gene tree bipartitions.
     */
    private double calculateBipartitionWeight(RangeBipartition speciesBip) {
        double totalScore = 0.0;

        if (Config.SCORE_TYPE == Config.ScoreType.QUARTET) {
            int lastTree = -1;
            int a = 0;
            int b = 0;
            int sg = 0;

            for (int i = 0; i < geneTreeRangeList.size(); i++) {
                RangeBipartition geneBip = geneTreeRangeList.get(i);
                int frequency = geneTreeFrequencies[i];
                int g = geneBip.geneTreeIndex;

                if (g != lastTree) {
                    a = getRangePresenceCount(
                        speciesBip.geneTreeIndex, speciesBip.leftStart, speciesBip.leftEnd, g);
                    b = getRangePresenceCount(
                        speciesBip.geneTreeIndex, speciesBip.rightStart, speciesBip.rightEnd, g);
                    sg = treeSizes[g];
                    lastTree = g;
                }

                double score = calculateRangeQuartetScore(speciesBip, geneBip, a, b, sg);
                totalScore += score * frequency;
            }
        } else {
            for (int i = 0; i < geneTreeRangeList.size(); i++) {
                RangeBipartition geneBip = geneTreeRangeList.get(i);
                int frequency = geneTreeFrequencies[i];

                double score = calculateRangeTripletScore(speciesBip, geneBip);
                totalScore += score * frequency;
            }
        }

        return totalScore;
    }
    
    /**
     * Calculate score between two RangeBipartitions using inverse index.
     * Same formula as MemoryOptimizedWeightCalculator.
     */
    private double calculateRangeTripletScore(RangeBipartition bip1, RangeBipartition bip2) {
        // Calculate four intersection sizes
        int aa = getRangeIntersectionSize(
            bip1.geneTreeIndex, bip1.leftStart, bip1.leftEnd,
            bip2.geneTreeIndex, bip2.leftStart, bip2.leftEnd);
        
        int bb = getRangeIntersectionSize(
            bip1.geneTreeIndex, bip1.rightStart, bip1.rightEnd,
            bip2.geneTreeIndex, bip2.rightStart, bip2.rightEnd);
        
        int ab = getRangeIntersectionSize(
            bip1.geneTreeIndex, bip1.leftStart, bip1.leftEnd,
            bip2.geneTreeIndex, bip2.rightStart, bip2.rightEnd);
        
        int ba = getRangeIntersectionSize(
            bip1.geneTreeIndex, bip1.rightStart, bip1.rightEnd,
            bip2.geneTreeIndex, bip2.leftStart, bip2.leftEnd);
        
        // Triplet scoring formula
        double score1 = 0;
        if (aa + bb >= 2) {
            score1 = aa * bb * (aa + bb - 2) / 2.0;
        }
        
        double score2 = 0;
        if (ab + ba >= 2) {
            score2 = ab * ba * (ab + ba - 2) / 2.0;
        }
        
        return score1 + score2;
    }

    private long quartetF(long a, long b, long c) {
        return (a + b + c - 3L) * a * b * c;
    }

    private long quartetFF(long a1, long a2, long a3,
                           long b1, long b2, long b3,
                           long c1, long c2, long c3) {
        return quartetF(a1, b2, c3) + quartetF(a1, b3, c2)
             + quartetF(a2, b1, c3) + quartetF(a2, b3, c1)
             + quartetF(a3, b1, c2) + quartetF(a3, b2, c1);
    }

    private double calculateRangeQuartetScore(RangeBipartition candidate, RangeBipartition geneTree,
                                              int a, int b, int sg) {
        int g = geneTree.geneTreeIndex;

        int ax = getRangeIntersectionSize(
            candidate.geneTreeIndex, candidate.leftStart, candidate.leftEnd,
            g, geneTree.leftStart, geneTree.leftEnd);
        int ay = getRangeIntersectionSize(
            candidate.geneTreeIndex, candidate.leftStart, candidate.leftEnd,
            g, geneTree.rightStart, geneTree.rightEnd);
        int bx = getRangeIntersectionSize(
            candidate.geneTreeIndex, candidate.rightStart, candidate.rightEnd,
            g, geneTree.leftStart, geneTree.leftEnd);
        int by = getRangeIntersectionSize(
            candidate.geneTreeIndex, candidate.rightStart, candidate.rightEnd,
            g, geneTree.rightStart, geneTree.rightEnd);

        int x = geneTree.leftSize();
        int y = geneTree.rightSize();
        int z = sg - x - y;

        int az = a - ax - ay;
        int bz = b - bx - by;

        int cx = x - ax - bx;
        int cy = y - ay - by;
        int cz = z - az - bz;

        long a1 = ax, a2 = ay, a3 = az;
        long b1 = bx, b2 = by, b3 = bz;
        long c1 = cx, c2 = cy, c3 = cz;

        return (double) quartetFF(a1, a2, a3, b1, b2, b3, c1, c2, c3);
    }
    
    /**
     * Calculate intersection size between two ranges using inverse index.
     * Complexity: O(min(|range1|, |range2|))
     */
    private int getRangeIntersectionSize(int tree1, int start1, int end1,
                                         int tree2, int start2, int end2) {
        if (start1 < 0 || end1 <= start1 || start2 < 0 || end2 <= start2) {
            return 0;
        }
        
        int size1 = end1 - start1;
        int size2 = end2 - start2;
        
        // Iterate over smaller range
        if (size1 <= size2) {
            return countIntersection(tree1, start1, end1, tree2, start2, end2);
        } else {
            return countIntersection(tree2, start2, end2, tree1, start1, end1);
        }
    }

    private int getRangePresenceCount(int sourceTree, int start, int end, int targetTree) {
        if (start < 0 || end < start) {
            return 0;
        }
        int rangeSize = end - start;
        if (rangeSize == 0) return 0;

        int missingCount = missingCounts[targetTree];
        if (missingCount == 0) {
            return rangeSize;
        }

        if (missingCount < rangeSize) {
            int missingInRange = 0;
            int offset = missingOffsets[targetTree];
            for (int i = 0; i < missingCount; i++) {
                int taxonId = missingTaxaFlat[offset + i];
                int posInSource = extendedInverseIndex[sourceTree][taxonId];
                if (posInSource >= start && posInSource < end) {
                    missingInRange++;
                }
            }
            return rangeSize - missingInRange;
        } else {
            int[] ordering = extendedOrderings[sourceTree];
            if (ordering == null) return 0;
            int maxPos = Math.min(end, ordering.length);
            int count = 0;
            for (int pos = start; pos < maxPos; pos++) {
                int taxonId = ordering[pos];
                if (taxonId < 0 || taxonId >= numTaxa) continue;
                if (extendedInverseIndex[targetTree][taxonId] != -1) {
                    count++;
                }
            }
            return count;
        }
    }
    
    /**
     * Count intersection by iterating over smaller range.
     */
    private int countIntersection(int smallTree, int smallStart, int smallEnd,
                                  int largeTree, int largeStart, int largeEnd) {
        int count = 0;
        int[] smallOrdering = extendedOrderings[smallTree];
        
        if (smallOrdering == null) return 0;
        
        int maxPos = Math.min(smallEnd, smallOrdering.length);
        
        for (int pos = smallStart; pos < maxPos; pos++) {
            int taxonId = smallOrdering[pos];
            
            if (taxonId < 0 || taxonId >= numTaxa) continue;
            
            int posInLargeTree = extendedInverseIndex[largeTree][taxonId];
            
            // Check if taxon exists in large tree and falls within range
            if (posInLargeTree >= 0 && posInLargeTree >= largeStart && posInLargeTree < largeEnd) {
                count++;
            }
        }
        
        return count;
    }
}
