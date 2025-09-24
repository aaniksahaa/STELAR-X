package preprocessing;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.atomic.AtomicInteger;

import utils.BitSet;
import utils.Threading;
import taxon.Taxon;
import tree.Tree;
import tree.TreeNode;
import tree.STBipartition;
import tree.RangeSTBipartition;
import tree.RangeHashCalculator;
import tree.GeneTreeOrderingGenerator;

/**
 * Memory-efficient GeneTrees that uses range-based representation during processing
 * to identify unique bipartitions, then converts only the unique ones to BitSets.
 * 
 * This reduces memory usage during the gene tree processing phase while keeping
 * the rest of the pipeline (weight calculation, DP inference) unchanged.
 */
public class MemoryEfficientGeneTrees extends GeneTrees {

    // Memory optimization data structures (used only during processing)
    private int[][] geneTreeOrderings;
    private RangeHashCalculator hashCalculator;
    private Map<RangeSTBipartition, Integer> rangeSTBipartitions; // Temporary range-based

    /**
     * Extracts taxon names from a single Newick tree string.
     */
    private void parseTaxa(String newickLine, Set<String> taxaSet){
        newickLine = newickLine.replaceAll("\\s", "");
    
        int n = newickLine.length();
        int i = 0;
    
        while(i < n){
            char curr = newickLine.charAt(i);
            if(curr == '('){
                i++;
            }
            else if(curr == ')'){
                i++;
                i = skipBranchInfo(newickLine, i);
            }
            else if(curr == ',' || curr == ';'){
                i++;
            }
            else{
                StringBuilder taxaName = new StringBuilder();
                while(i < n){
                    char c = newickLine.charAt(i);
                    if(c == ')' || c == ',' || c == ':' || c == ';'){
                        break;
                    }
                    taxaName.append(c);
                    i++;
                }
                
                String label = taxaName.toString();
                if (!label.isEmpty()) {
                    taxaSet.add(label);
                }
                
                i = skipBranchInfo(newickLine, i);
            }
        }
    }
    
    private int skipBranchInfo(String newickLine, int startIndex) {
        int i = startIndex;
        int n = newickLine.length();
        
        if (i < n && newickLine.charAt(i) == ':') {
            i++;
            while (i < n) {
                char c = newickLine.charAt(i);
                if (c == ',' || c == ')' || c == ';') {
                    break;
                }
                i++;
            }
        }
        
        return i;
    }

    public Map<String, Taxon> readTaxaNames() throws FileNotFoundException{
        Set<String> taxaSet = new HashSet<>();

        Scanner scanner = new Scanner(new File(this.path));
        while(scanner.hasNextLine()){
            String line = scanner.nextLine();
            if(line.trim().length() == 0) continue;
            parseTaxa(line, taxaSet);
        }
        scanner.close();

        this.taxaMap = new HashMap<>();
        for(var x : taxaSet){
            Taxon taxon = new Taxon(x);
            taxaMap.put(x, taxon);
        }

        return taxaMap;
    }

    public void readGeneTrees(double[][] distanceMatrix) throws FileNotFoundException{
        readGeneTreesParallel(distanceMatrix);
    }

    /**
     * Memory-efficient parallel gene tree reading that uses range-based processing
     * to identify unique bipartitions before converting to BitSets.
     */
    private void readGeneTreesParallel(double[][] distanceMatrix) throws FileNotFoundException{
        // First, read all lines from the file
        List<String> lines = new ArrayList<>();
        Scanner scanner = new Scanner(new File(path));
        while (scanner.hasNextLine()) {
            String line = scanner.nextLine();
            if(line.trim().length() != 0) {
                lines.add(line);
            }
        }
        scanner.close();
        
        if (lines.isEmpty()) {
            System.err.println("Warning: No valid gene trees found in input file.");
            return;
        }
        
        System.out.println("Processing " + lines.size() + " gene trees with memory-efficient approach...");
        
        // Thread-safe data structures
        List<Tree> threadSafeGeneTrees = new ArrayList<>();
        AtomicInteger totalInternalNodes = new AtomicInteger(0);
        AtomicInteger skippedTrees = new AtomicInteger(0);
        
        int numThreads = Runtime.getRuntime().availableProcessors();
        Threading.startThreading(numThreads);
        
        int chunkSize = (lines.size() + numThreads - 1) / numThreads;
        CountDownLatch latch = new CountDownLatch(numThreads);

        System.out.println("Using " + numThreads + " threads for parallel gene tree processing");
        
        // First pass: Parse all trees (same as original)
        for (int i = 0; i < numThreads; i++) {
            final int startIdx = i * chunkSize;
            final int endIdx = Math.min(startIdx + chunkSize, lines.size());
            final int threadId = i;
            
            Threading.execute(() -> {
                try {
                    List<Tree> localTrees = new ArrayList<>();
                    int localInternalNodes = 0;
                    int localSkippedTrees = 0;
                    
                    System.out.println("Thread " + threadId + " processing trees " + startIdx + " to " + (endIdx-1));
                    
                    for (int j = startIdx; j < endIdx; j++) {
                        String line = lines.get(j);
                        try {
                            Tree tree = new Tree(line, this.taxaMap);
                            localTrees.add(tree);
                            localInternalNodes += tree.nodes.size() - tree.leavesCount;
                        } catch (RuntimeException e) {
                            localSkippedTrees++;
                        }
                    }
                    
                    // Merge local results
                    synchronized (threadSafeGeneTrees) {
                        threadSafeGeneTrees.addAll(localTrees);
                    }
                    
                    totalInternalNodes.addAndGet(localInternalNodes);
                    skippedTrees.addAndGet(localSkippedTrees);
                    
                    System.out.println("Thread " + threadId + " completed: processed " + 
                                     localTrees.size() + " trees, skipped " + localSkippedTrees);
                    
                } finally {
                    latch.countDown();
                }
            });
        }
        
        try {
            latch.await();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new RuntimeException("Gene tree processing was interrupted", e);
        } finally {
            Threading.shutdown();
        }
        
        // Copy results and initialize data structures
        this.geneTrees.addAll(threadSafeGeneTrees);
        this.taxonIdToLabel = new String[this.taxaMap.size()];
        this.taxa = new Taxon[this.taxaMap.size()];
        this.realTaxaCount = this.taxaMap.size();

        for(var x : this.taxaMap.entrySet()){
            taxonIdToLabel[x.getValue().id] = x.getKey();
            taxa[x.getValue().id] = x.getValue();
        }

        if (skippedTrees.get() > 0) {
            System.err.println("Warning: Skipped " + skippedTrees.get() + " invalid or empty gene trees.");
        }

        // MEMORY OPTIMIZATION: Use range-based approach to identify unique bipartitions
        System.out.println("Phase 1: Generating gene tree orderings for memory optimization...");
        this.geneTreeOrderings = GeneTreeOrderingGenerator.generateOrderings(this.geneTrees, this.realTaxaCount);
        this.hashCalculator = new RangeHashCalculator(this.geneTrees, this.geneTreeOrderings);
        
        System.out.println("Phase 2: Identifying unique bipartitions using range-based approach...");
        this.rangeSTBipartitions = new HashMap<>();
        generateRangeSTBipartitions();
        
        System.out.println("Phase 3: Converting only unique bipartitions to BitSets...");
        this.stBipartitions = new HashMap<>();
        convertUniqueRangesToBitSets();
        
        // Clean up temporary data structures to free memory
        this.geneTreeOrderings = null;
        this.hashCalculator = null;
        this.rangeSTBipartitions = null;

        // Output statistics
        System.out.println("taxon count: " + this.taxaMap.size());
        System.out.println("Gene trees count: " + geneTrees.size());
        System.out.println("total internal nodes: " + totalInternalNodes.get());
        System.out.println("unique STBipartitions (memory optimized): " + stBipartitions.size());
    }

    /**
     * Generate range-based bipartitions from all gene trees to identify unique ones.
     */
    private void generateRangeSTBipartitions() {
        int numThreads = Runtime.getRuntime().availableProcessors();
        Threading.startThreading(numThreads);
        
        Map<RangeSTBipartition, Integer> globalRangeSTBipartitions = new ConcurrentHashMap<>();
        int treesPerThread = (geneTrees.size() + numThreads - 1) / numThreads;
        CountDownLatch latch = new CountDownLatch(numThreads);
        
        for (int threadIdx = 0; threadIdx < numThreads; threadIdx++) {
            final int startTreeIdx = threadIdx * treesPerThread;
            final int endTreeIdx = Math.min(startTreeIdx + treesPerThread, geneTrees.size());
            
            Threading.execute(() -> {
                try {
                    Map<RangeSTBipartition, Integer> localRangeSTBipartitions = new HashMap<>();
                    
                    for (int treeIdx = startTreeIdx; treeIdx < endTreeIdx; treeIdx++) {
                        Tree tree = geneTrees.get(treeIdx);
                        generateRangeSTBipartitionsForTree(tree, treeIdx, localRangeSTBipartitions);
                    }
                    
                    // Compute hashes for all local bipartitions
                    for (RangeSTBipartition rangeBip : localRangeSTBipartitions.keySet()) {
                        rangeBip.computeHashes(hashCalculator);
                    }
                    
                    // Merge with global map
                    for (Map.Entry<RangeSTBipartition, Integer> entry : localRangeSTBipartitions.entrySet()) {
                        globalRangeSTBipartitions.merge(entry.getKey(), entry.getValue(), Integer::sum);
                    }
                    
                } finally {
                    latch.countDown();
                }
            });
        }
        
        try {
            latch.await();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new RuntimeException("Range bipartition generation was interrupted", e);
        } finally {
            Threading.shutdown();
        }
        
        this.rangeSTBipartitions = globalRangeSTBipartitions;
        System.out.println("Generated " + rangeSTBipartitions.size() + " unique range-based bipartitions");
    }

    /**
     * Generate range-based bipartitions for a single tree.
     */
    private void generateRangeSTBipartitionsForTree(Tree tree, int treeIndex, 
                                                   Map<RangeSTBipartition, Integer> localRangeSTBipartitions) {
        if (tree.isRooted()) {
            generateRangeSTBipartitionsUtil(tree.root, tree, treeIndex, localRangeSTBipartitions);
        }
    }

    /**
     * Recursive utility to generate range bipartitions from tree structure.
     */
    private int[] generateRangeSTBipartitionsUtil(TreeNode node, Tree tree, int treeIndex,
                                                 Map<RangeSTBipartition, Integer> localRangeSTBipartitions) {
        if (node.isLeaf()) {
            // Find position of this taxon in the gene tree ordering
            int[] ordering = geneTreeOrderings[treeIndex];
            for (int i = 0; i < ordering.length; i++) {
                if (ordering[i] == node.taxon.id) {
                    return new int[]{i, i + 1};  // Single taxon range [i, i+1)
                }
            }
            return new int[]{-1, -1}; // Taxon not found
        }
        
        if (node.childs != null && node.childs.size() == 2 && tree.isRooted()) {
            // Get ranges for both children
            int[] leftRange = generateRangeSTBipartitionsUtil(node.childs.get(0), tree, treeIndex, localRangeSTBipartitions);
            int[] rightRange = generateRangeSTBipartitionsUtil(node.childs.get(1), tree, treeIndex, localRangeSTBipartitions);
            
            if (leftRange[0] != -1 && rightRange[0] != -1) {
                // Check if ranges are contiguous (adjacent)
                boolean leftFirst = (leftRange[1] == rightRange[0]);
                boolean rightFirst = (rightRange[1] == leftRange[0]);
                
                if (leftFirst || rightFirst) {
                    // Create range bipartition for one of the clusters
                    int minPos = Math.min(leftRange[0], rightRange[0]);
                    int maxPos = Math.max(leftRange[1], rightRange[1]);
                    
                    // Create bipartition for the left cluster
                    RangeSTBipartition rangeBip = new RangeSTBipartition(treeIndex, 
                        leftFirst ? leftRange[0] : rightRange[0], 
                        leftFirst ? leftRange[1] : rightRange[1]);
                    localRangeSTBipartitions.merge(rangeBip, 1, Integer::sum);
                    
                    return new int[]{minPos, maxPos};
                }
            }
        }
        
        // Handle other cases (multi-way split, etc.)
        if (node.childs != null && !node.childs.isEmpty()) {
            int minPos = Integer.MAX_VALUE;
            int maxPos = Integer.MIN_VALUE;
            
            for (TreeNode child : node.childs) {
                int[] childRange = generateRangeSTBipartitionsUtil(child, tree, treeIndex, localRangeSTBipartitions);
                if (childRange[0] != -1) {
                    minPos = Math.min(minPos, childRange[0]);
                    maxPos = Math.max(maxPos, childRange[1]);
                }
            }
            
            if (minPos != Integer.MAX_VALUE) {
                return new int[]{minPos, maxPos};
            }
        }
        
        return new int[]{-1, -1};
    }

    /**
     * Convert only the unique range-based bipartitions to BitSets.
     * This is where we save memory - only unique bipartitions become BitSets.
     */
    private void convertUniqueRangesToBitSets() {
        System.out.println("Converting " + rangeSTBipartitions.size() + " unique range bipartitions to BitSets...");
        
        for (Map.Entry<RangeSTBipartition, Integer> entry : rangeSTBipartitions.entrySet()) {
            RangeSTBipartition rangeBip = entry.getKey();
            int frequency = entry.getValue();
            
            // Convert range to BitSet
            BitSet cluster1 = new BitSet(realTaxaCount);
            BitSet cluster2 = new BitSet(realTaxaCount);
            
            // Get taxa for cluster1 (the range)
            int[] ordering = geneTreeOrderings[rangeBip.geneTreeIndex];
            for (int i = rangeBip.startPos; i < rangeBip.endPos; i++) {
                if (i < ordering.length) {
                    cluster1.set(ordering[i]);
                }
            }
            
            // Get taxa for cluster2 (the complement)
            for (int i = 0; i < ordering.length; i++) {
                if (i < rangeBip.startPos || i >= rangeBip.endPos) {
                    cluster2.set(ordering[i]);
                }
            }
            
            // Create STBipartition and add to final map
            STBipartition stbip = new STBipartition(cluster1, cluster2);
            stBipartitions.merge(stbip, frequency, Integer::sum);
        }
        
        System.out.println("Successfully converted to " + stBipartitions.size() + " unique BitSet-based bipartitions");
    }

    // Constructors
    public MemoryEfficientGeneTrees(String path) throws FileNotFoundException {
        super(path);
    }

    public MemoryEfficientGeneTrees(String path, Map<String, Taxon> taxaMap) throws FileNotFoundException {
        super(path, taxaMap);
    }

    /**
     * Generate candidate bipartitions (same as original GeneTrees).
     */
    public List<STBipartition> generateCandidateBipartitions() {
        return generateCandidateBipartitions(true);
    }
    
    public List<STBipartition> generateCandidateBipartitions(boolean useExpansion) {
        List<STBipartition> candidates = new ArrayList<>();
        Set<BitSet> processedClusters = new HashSet<>();
        
        System.out.println("\n ***************** Adding all bipartitions from gene trees (memory optimized) *****************\n");
        for (STBipartition stb : stBipartitions.keySet()) {
            candidates.add(stb);
            processedClusters.add(stb.cluster1);
            processedClusters.add(stb.cluster2);
        }
        
        System.out.println("\n ***************** Adding complementary bipartitions for each cluster *****************\n");
        for (STBipartition stb : stBipartitions.keySet()) {
            BitSet cluster1 = stb.cluster1;
            BitSet cluster2 = stb.cluster2;
            
            if (!processedClusters.contains(cluster1)) {
                BitSet complement = new BitSet(realTaxaCount);
                complement.set(0, realTaxaCount);
                complement.andNot(cluster1);
                if (complement.cardinality() > 0) {
                    candidates.add(new STBipartition(cluster1, complement));
                    processedClusters.add(cluster1);
                    processedClusters.add(complement);
                }
            }
            
            if (!processedClusters.contains(cluster2)) {
                BitSet complement = new BitSet(realTaxaCount);
                complement.set(0, realTaxaCount);
                complement.andNot(cluster2);
                if (complement.cardinality() > 0) {
                    candidates.add(new STBipartition(cluster2, complement));
                    processedClusters.add(cluster2);
                    processedClusters.add(complement);
                }
            }
        }
        
        // Apply expansion if enabled (same as original)
        if (useExpansion && utils.BipartitionExpansionConfig.isExpansionEnabled()) {
            System.out.println("\n ***************** Applying bipartition expansion *****************\n");
            expansion.BipartitionExpansionManager expansionManager = 
                new expansion.BipartitionExpansionManager(this);
            candidates = expansionManager.expandBipartitions(candidates);
        }
        
        return candidates;
    }
}
