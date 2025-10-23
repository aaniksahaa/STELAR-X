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
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.atomic.AtomicInteger;

import utils.BitSet;
import utils.Threading;
import taxon.Taxon;
import tree.Tree;
import tree.RangeBipartition;
import tree.MemoryEfficientBipartitionManager;

/**
 * GeneTrees: Preprocessing and management of input gene trees for wQFM-TREE.
 * 
 * This class handles the critical first step of the wQFM-TREE algorithm: preprocessing
 * the input gene trees for efficient use in Algorithm 2 scoring. It performs several
 * key operations described in the paper:
 * 
 * 1. Taxa name extraction and standardization across all gene trees
 * 2. Gene tree parsing from Newick format with consistent taxon mapping
 * 3. Optional polytomy resolution for binary tree structure requirement
 * 4. Tripartition frequency calculation for consensus tree construction
 * 5. Data structure initialization for efficient Algorithm 2 operations
 * 
 * The preprocessing ensures that all gene trees use consistent taxon IDs and
 * are structured appropriately for the O(n²k log n) scoring complexity described
 * in Section 2.2. The tripartition data supports Algorithm 1 for initial
 * bipartition generation.
 */
public class GeneTrees {

    // Core data structures
    public ArrayList<Tree> geneTrees;               // Parsed and preprocessed gene trees
    public String[] taxonIdToLabel;                 // ID to label mapping for output
    public Taxon[] taxa;                        // Array of all real taxa (indexed by ID)
    // public Map<String, TreeNode> triPartitions;     // Tripartition frequency data for consensus
    public Map<RangeBipartition, Integer> rangeBipartitions; // RangeBipartition frequency map
    public Map<String, Taxon> taxaMap;          // Label to RealTaxon mapping
    public int realTaxaCount;                       // Total number of real taxa
    public String path;                             // Input file path

    /**
     * Extracts taxon names from a single Newick tree string.
     * 
     * This helper method parses a Newick string to collect all taxon names,
     * building the comprehensive set of taxa that appear across all gene trees.
     * This is essential for creating consistent taxon IDs used throughout
     * the algorithm.
     */
    private void parseTaxa(String newickLine, Set<String> taxaSet){
        newickLine = newickLine.replaceAll("\\s", "");
    
        int n = newickLine.length();
        int i = 0;
    
        // Parse Newick format to extract taxon names
        while(i < n){
            char curr = newickLine.charAt(i);
            if(curr == '('){
                // Start of internal node - skip
                i++;
            }
            else if(curr == ')'){
                // End of internal node - skip any support values or branch lengths
                i++;
                i = skipBranchInfo(newickLine, i);
            }
            else if(curr == ',' || curr == ';'){
                // Separators - skip
                i++;
            }
            else{
                // Taxon name - extract and add to set
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
                
                // Skip any branch length information
                i = skipBranchInfo(newickLine, i);
            }
        }
    }
    
    private int skipBranchInfo(String newickLine, int startIndex) {
        int i = startIndex;
        int n = newickLine.length();
        
        // Skip branch length (starts with ':')
        if (i < n && newickLine.charAt(i) == ':') {
            i++; // skip ':'
            // Skip the numeric value (including scientific notation)
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
     * 
     * This method performs the first pass through all gene trees to:
     * 1. Collect all unique taxon names that appear in any gene tree
     * 2. Assign consistent IDs to taxa for use throughout the algorithm
     * 3. Create the taxaMap used for consistent tree parsing
     * 
     * The consistent taxon ID assignment is crucial for Algorithm 2 scoring
     * efficiency, as it enables O(1) taxon lookup operations.
     */
    public Map<String, Taxon> readTaxaNames() throws FileNotFoundException{
        Set<String> taxaSet = new HashSet<>();

        // First pass: collect all taxon names
        Scanner scanner = new Scanner(new File(this.path));
        while(scanner.hasNextLine()){
            String line = scanner.nextLine();
            if(line.trim().length() == 0) continue;
            parseTaxa(line, taxaSet);
        }
        scanner.close();

        // Create RealTaxon objects with consistent IDs
        this.taxaMap = new HashMap<>();
        for(var x : taxaSet){
            Taxon taxon = new Taxon(x);
            taxaMap.put(x, taxon);
        }

        return taxaMap;
    }

    /**
     * Reads and preprocesses all gene trees for wQFM-TREE algorithm.
     * 
     * This method performs the main gene tree preprocessing described in Section 2:
     * 1. Parses each gene tree from Newick format using consistent taxon mapping
     * 2. Optionally resolves polytomies to ensure binary tree structure
     * 3. Calculates tripartition frequencies for Algorithm 1 (consensus tree construction)
     * 4. Prepares trees for efficient Algorithm 2 scoring operations
     * 
     * The preprocessing ensures that:
     * - All trees use consistent taxon IDs for O(1) lookup in scoring
     * - Tree structure supports efficient quartet evaluation
     * - Tripartition data enables consensus-based initial bipartition generation
     * 
     * @param distanceMatrix Optional distance matrix for polytomy resolution
     */
    public void readGeneTrees(double[][] distanceMatrix) throws FileNotFoundException{
        readGeneTreesParallel(distanceMatrix);
    }

    /**
     * Memory-efficient parallel version of gene tree reading.
     * Uses range-based bipartition representation to reduce memory usage from O(n²k) to O(nk).
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
        
        // Thread-safe data structures for tree parsing
        List<Tree> threadSafeGeneTrees = new ArrayList<>();
        AtomicInteger totalInternalNodes = new AtomicInteger(0);
        AtomicInteger skippedTrees = new AtomicInteger(0);
        
        int numThreads = Runtime.getRuntime().availableProcessors();
        Threading.startThreading(numThreads);
        
        // Calculate optimal number of threads to avoid invalid ranges
        int chunkSize = Math.max(1, (lines.size() + numThreads - 1) / numThreads);
        int actualThreads = Math.min(numThreads, (lines.size() + chunkSize - 1) / chunkSize);
        CountDownLatch latch = new CountDownLatch(actualThreads);

        System.out.println("Using " + actualThreads + " threads for parallel gene tree parsing");
        
        // Process gene trees in parallel chunks - PARSING ONLY
        for (int i = 0; i < actualThreads; i++) {
            final int startIdx = i * chunkSize;
            final int endIdx = Math.min(startIdx + chunkSize, lines.size());
            final int threadId = i;
            
            Threading.execute(() -> {
                try {
                    List<Tree> localTrees = new ArrayList<>();
                    int localInternalNodes = 0;
                    int localSkippedTrees = 0;
                    
                    // Validate range before processing
                    if (startIdx >= endIdx || startIdx >= lines.size()) {
                        System.out.println("Thread " + threadId + " skipped - invalid range [" + startIdx + ", " + endIdx + ")");
                        return;
                    }
                    
                    System.out.println("Thread " + threadId + " parsing trees " + startIdx + " to " + (endIdx-1));
                    
                    for (int j = startIdx; j < endIdx; j++) {
                        String line = lines.get(j);
                        try {
                            // Parse tree with consistent taxon mapping
                            Tree tree = new Tree(line, this.taxaMap);
                            localTrees.add(tree);
                            
                            // Track internal node count for complexity analysis
                            localInternalNodes += tree.nodes.size() - tree.leavesCount;
                        } catch (RuntimeException e) {
                            localSkippedTrees++;
                            System.err.println("Thread " + threadId + " skipped invalid tree at line " + (j + 1) + ": " + e.getMessage());
                            // Continue with next tree
                        }
                    }
                    
                    // Merge local results into thread-safe structures
                    synchronized (threadSafeGeneTrees) {
                        threadSafeGeneTrees.addAll(localTrees);
                    }
                    
                    totalInternalNodes.addAndGet(localInternalNodes);
                    skippedTrees.addAndGet(localSkippedTrees);
                    
                    System.out.println("Thread " + threadId + " completed: parsed " + 
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
        
        // Copy results to instance variables
        this.geneTrees.addAll(threadSafeGeneTrees);
        
        if (skippedTrees.get() > 0) {
            System.err.println("Warning: Skipped " + skippedTrees.get() + " invalid or empty gene trees.");
        }
        
        // Build efficient lookup arrays for Algorithm 2 operations
        this.taxonIdToLabel = new String[this.taxaMap.size()];
        this.taxa = new Taxon[this.taxaMap.size()];
        this.realTaxaCount = this.taxaMap.size();

        // Populate lookup arrays indexed by taxon ID
        for(var x : this.taxaMap.entrySet()){
            taxonIdToLabel[x.getValue().id] = x.getKey();
            taxa[x.getValue().id] = x.getValue();
        }

        // Now use the memory-efficient bipartition manager to process RangeBipartitions
        System.out.println("\n=== Memory-Efficient RangeBipartition Processing ===");
        MemoryEfficientBipartitionManager bipartitionManager = 
            new MemoryEfficientBipartitionManager(this.geneTrees, this.realTaxaCount);
        
        // Process gene trees using range-based representation
        this.rangeBipartitions = bipartitionManager.processGeneTreesParallel();
        
        // Output preprocessing statistics
        System.out.println("\n=== Gene Tree Processing Complete ===");
        System.out.println("Taxon count: " + this.taxaMap.size());
        System.out.println("Gene trees count: " + geneTrees.size());
        System.out.println("Total internal nodes: " + totalInternalNodes.get());
        System.out.println("Unique RangeBipartitions: " + rangeBipartitions.size());
        System.out.println("\n" + bipartitionManager.getStatistics());
    }

    /**
     * Constructor: Initialize gene trees preprocessing from file path.
     * 
     * This constructor sets up the GeneTrees object for processing gene trees
     * from the specified file. The actual preprocessing happens in subsequent
     * method calls to readTaxaNames() and readGeneTrees().
     */
    public GeneTrees(String path) throws FileNotFoundException{

        this.geneTrees = new ArrayList<>();
        // this.triPartitions = new HashMap<>();
        this.rangeBipartitions = new HashMap<>();
        this.path = path;
    }


    /**
     * Constructor: Initialize with existing taxon mapping.
     * 
     * This constructor is used when taxon mapping is already established,
     * immediately proceeding to gene tree preprocessing.
     */
    public GeneTrees(String path, Map<String, Taxon> taxaMap) throws FileNotFoundException{

        this.geneTrees = new ArrayList<>();
        // this.triPartitions = new HashMap<>();
        this.rangeBipartitions = new HashMap<>();
        this.path = path;
        this.taxaMap = taxaMap;

        this.readGeneTrees(null);
    }
    
    // Legacy STBipartition calculation methods removed - now using MemoryEfficientBipartitionManager

    /**
     * Generates candidate bipartitions for the inference algorithm.
     * This method creates a list of candidate bipartitions from the gene tree bipartitions,
     * ensuring that each candidate is a valid bipartition of the taxa set.
     * 
     * @return List of candidate bipartitions for inference
     */
    public List<RangeBipartition> generateCandidateBipartitions() {
        return generateCandidateBipartitions(true);
    }
    
    /**
     * Generates candidate bipartitions with optional expansion.
     * 
     * @param useExpansion Whether to apply bipartition expansion techniques
     * @return List of candidate bipartitions for inference
     */
    public List<RangeBipartition> generateCandidateBipartitions(boolean useExpansion) {
        List<RangeBipartition> candidates = new ArrayList<>();
        
        // Simply return all unique range bipartitions from gene trees
        System.out.println("\n ***************** Adding all bipartitions from gene trees *****************\n");
        for (RangeBipartition range : rangeBipartitions.keySet()) {
            candidates.add(range);
        }
        
        // Note: Complementary bipartition generation is more complex with ranges
        // and may not be necessary since we're working with subtree bipartitions
        // that already represent the natural structure of the trees.
        // If needed, this can be implemented by:
        // 1. Converting ranges to taxon sets
        // 2. Computing complements
        // 3. Finding ranges that represent those complements
        // But this would defeat the memory optimization purpose.
        
        // Apply bipartition expansion if enabled
        // TODO: Update BipartitionExpansionManager to work with RangeBipartition
        // if (useExpansion && utils.BipartitionExpansionConfig.isExpansionEnabled()) {
        //     System.out.println("\n ***************** Applying bipartition expansion *****************\n");
        //     expansion.BipartitionExpansionManager expansionManager = 
        //         new expansion.BipartitionExpansionManager(this);
        //     candidates = expansionManager.expandBipartitions(candidates);
        // }
        
        return candidates;
    }
    
    // Removed addValidBipartition method - no longer used with memory-efficient approach
} 


