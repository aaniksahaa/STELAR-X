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

import utils.Threading;
import taxon.Taxon;
import tree.Tree;
import tree.TreeNode;
import tree.RangeSTBipartition;
import tree.RangeHashCalculator;
import tree.GeneTreeOrderingGenerator;

/**
 * Memory-optimized version of GeneTrees using range-based bipartition representation.
 * 
 * Key optimizations:
 * 1. Replaces BitSet-based STBipartitions with range-based RangeSTBipartitions
 * 2. Uses (geneTreeIndex, startPos, endPos) representation instead of BitSets
 * 3. Reduces memory usage from O(n²k) to O(nk) where n=taxa, k=gene trees
 * 4. Uses hash-based equality checking for cross-tree bipartition comparison
 */
public class MemoryOptimizedGeneTrees {

    // Core data structures
    public ArrayList<Tree> geneTrees;                    // Parsed and preprocessed gene trees
    public String[] taxonIdToLabel;                      // ID to label mapping for output
    public Taxon[] taxa;                                 // Array of all real taxa (indexed by ID)
    public Map<RangeSTBipartition, Integer> rangeSTBipartitions; // Range-based bipartition frequency map
    public Map<String, Taxon> taxaMap;                   // Label to RealTaxon mapping
    public int realTaxaCount;                            // Total number of real taxa
    public String path;                                  // Input file path

    // Memory optimization data structures
    public int[][] geneTreeOrderings;                    // [geneTreeIndex][position] = taxonId
    public RangeHashCalculator hashCalculator;           // Hash calculator for range comparison
    
    // Compatibility mode for gradual migration
    public boolean useRangeRepresentation = true;        // Toggle between old/new representation

    /**
     * Extracts taxon names from a single Newick tree string.
     */
    private void parseTaxa(String newickLine, Set<String> taxaSet) {
        newickLine = newickLine.replaceAll("\\s", "");
    
        int n = newickLine.length();
        int i = 0;
    
        // Parse Newick format to extract taxon names
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

    /**
     * Reads and standardizes taxon names across all gene trees.
     */
    public Map<String, Taxon> readTaxaNames() throws FileNotFoundException {
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

    /**
     * Reads and preprocesses all gene trees for memory-optimized representation.
     */
    public void readGeneTrees(double[][] distanceMatrix) throws FileNotFoundException {
        readGeneTreesParallel(distanceMatrix);
    }

    /**
     * Parallel version of gene tree reading with memory optimization.
     */
    private void readGeneTreesParallel(double[][] distanceMatrix) throws FileNotFoundException {
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
        
        System.out.println("Processing " + lines.size() + " gene trees with memory optimization...");
        
        // Thread-safe data structures
        List<Tree> threadSafeGeneTrees = new ArrayList<>();
        Map<RangeSTBipartition, Integer> threadSafeRangeSTBipartitions = new ConcurrentHashMap<>();
        AtomicInteger totalInternalNodes = new AtomicInteger(0);
        AtomicInteger skippedTrees = new AtomicInteger(0);
        
        int numThreads = Runtime.getRuntime().availableProcessors();
        Threading.startThreading(numThreads);
        
        int chunkSize = (lines.size() + numThreads - 1) / numThreads;
        CountDownLatch latch = new CountDownLatch(numThreads);

        System.out.println("Using " + numThreads + " threads for parallel gene tree processing");
        
        // First pass: Parse all trees
        for (int i = 0; i < numThreads; i++) {
            final int startIdx = i * chunkSize;
            final int endIdx = Math.min(startIdx + chunkSize, lines.size());
            final int threadId = i;
            
            Threading.execute(() -> {
                try {
                    List<Tree> localTrees = new ArrayList<>();
                    int localInternalNodes = 0;
                    int localSkippedTrees = 0;
                    
                    System.out.println("Thread " + threadId + " parsing trees " + startIdx + " to " + (endIdx-1));
                    
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
                    
                    System.out.println("Thread " + threadId + " completed parsing: " + 
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
            throw new RuntimeException("Gene tree parsing was interrupted", e);
        } finally {
            Threading.shutdown();
        }
        
        // Copy parsed trees and initialize data structures
        this.geneTrees.addAll(threadSafeGeneTrees);
        this.taxonIdToLabel = new String[this.taxaMap.size()];
        this.taxa = new Taxon[this.taxaMap.size()];
        this.realTaxaCount = this.taxaMap.size();

        for(var x : this.taxaMap.entrySet()){
            taxonIdToLabel[x.getValue().id] = x.getKey();
            taxa[x.getValue().id] = x.getValue();
        }

        // Generate gene tree orderings for range-based representation
        System.out.println("Generating gene tree orderings for memory optimization...");
        this.geneTreeOrderings = GeneTreeOrderingGenerator.generateOrderings(this.geneTrees, this.realTaxaCount);
        
        // Initialize hash calculator
        this.hashCalculator = new RangeHashCalculator(this.geneTrees, this.geneTreeOrderings);
        
        // Debug: Print sample orderings
        GeneTreeOrderingGenerator.printOrderings(this.geneTreeOrderings, this.taxonIdToLabel);
        
        // Second pass: Generate range-based bipartitions
        System.out.println("Generating range-based bipartitions...");
        generateRangeSTBipartitions(threadSafeRangeSTBipartitions);
        
        this.rangeSTBipartitions = threadSafeRangeSTBipartitions;

        if (skippedTrees.get() > 0) {
            System.err.println("Warning: Skipped " + skippedTrees.get() + " invalid or empty gene trees.");
        }

        // Output preprocessing statistics
        System.out.println("taxon count: " + this.taxaMap.size());
        System.out.println("Gene trees count: " + geneTrees.size());
        System.out.println("total internal nodes: " + totalInternalNodes.get());
        System.out.println("unique range STBipartitions: " + rangeSTBipartitions.size());
        
        // Validate orderings
        boolean valid = GeneTreeOrderingGenerator.validateOrderings(this.geneTreeOrderings, this.realTaxaCount);
        System.out.println("Gene tree orderings validation: " + (valid ? "PASSED" : "FAILED"));
    }

    /**
     * Generates range-based bipartitions from all gene trees.
     */
    private void generateRangeSTBipartitions(Map<RangeSTBipartition, Integer> globalRangeSTBipartitions) {
        int numThreads = Runtime.getRuntime().availableProcessors();
        Threading.startThreading(numThreads);
        
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
    }

    /**
     * Generates range-based bipartitions for a single tree.
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
            return new int[]{-1, -1}; // Taxon not found (shouldn't happen)
        }
        
        if (node.childs.size() == 2 && tree.isRooted()) {
            // Get ranges for both children
            int[] leftRange = generateRangeSTBipartitionsUtil(node.childs.get(0), tree, treeIndex, localRangeSTBipartitions);
            int[] rightRange = generateRangeSTBipartitionsUtil(node.childs.get(1), tree, treeIndex, localRangeSTBipartitions);
            
            if (leftRange[0] != -1 && rightRange[0] != -1) {
                // Determine the overall range covered by this subtree
                int minPos = Math.min(leftRange[0], rightRange[0]);
                int maxPos = Math.max(leftRange[1], rightRange[1]);
                
                // Check if this forms a contiguous range (valid bipartition)
                if (areRangesContiguous(leftRange, rightRange)) {
                    // Create range bipartition for the left cluster
                    RangeSTBipartition rangeBip = new RangeSTBipartition(treeIndex, minPos, maxPos);
                    localRangeSTBipartitions.merge(rangeBip, 1, Integer::sum);
                }
                
                return new int[]{minPos, maxPos};
            }
        } else {
            // Multi-way split or unary node - merge all child ranges
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
     * Checks if two ranges can form a contiguous bipartition.
     */
    private boolean areRangesContiguous(int[] range1, int[] range2) {
        // For simplicity, we'll consider ranges contiguous if they don't overlap
        // and together form a valid bipartition. This is a simplified check.
        return (range1[1] <= range2[0]) || (range2[1] <= range1[0]);
    }

    // Constructors
    public MemoryOptimizedGeneTrees(String path) throws FileNotFoundException {
        this.geneTrees = new ArrayList<>();
        this.rangeSTBipartitions = new HashMap<>();
        this.path = path;
    }

    public MemoryOptimizedGeneTrees(String path, Map<String, Taxon> taxaMap) throws FileNotFoundException {
        this.geneTrees = new ArrayList<>();
        this.rangeSTBipartitions = new HashMap<>();
        this.path = path;
        this.taxaMap = taxaMap;
        this.readGeneTrees(null);
    }

    /**
     * Generates candidate bipartitions using the new range-based representation.
     */
    public List<RangeSTBipartition> generateRangeCandidateBipartitions() {
        List<RangeSTBipartition> candidates = new ArrayList<>();
        
        System.out.println("\n***************** Adding all range bipartitions from gene trees *****************\n");
        for (RangeSTBipartition rangeBip : rangeSTBipartitions.keySet()) {
            candidates.add(rangeBip);
        }
        
        System.out.println("Generated " + candidates.size() + " range-based candidate bipartitions");
        return candidates;
    }
    
    /**
     * Gets the hash calculator for range operations.
     */
    public RangeHashCalculator getHashCalculator() {
        return hashCalculator;
    }
    
    /**
     * Gets the gene tree orderings.
     */
    public int[][] getGeneTreeOrderings() {
        return geneTreeOrderings;
    }
}
