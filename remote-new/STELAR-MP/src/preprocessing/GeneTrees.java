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

import utils.Config;
import utils.BitSet;
import taxon.Taxon;
import tree.Tree;
import tree.TreeNode;
import tree.STBipartition;

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
    public Map<String, TreeNode> triPartitions;     // Tripartition frequency data for consensus
    public Map<STBipartition, Integer> stBipartitions; // STBipartition frequency map
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
        int internalNodesCount = 0;
        int skippedTrees = 0;

        Scanner scanner = new Scanner(new File(path));

        // Process each gene tree
        while (scanner.hasNextLine()) {

            String line = scanner.nextLine();
            if(line.trim().length() == 0) continue;
            
            try {
                // Parse tree with consistent taxon mapping
                var tree = new Tree(line, this.taxaMap);
                
                // Trees are assumed to be properly rooted and binary
                // No polytomy resolution needed
                
                // Calculate tripartition frequencies for Algorithm 1
                tree.calculateFrequencies(triPartitions);
                
                // Calculate STBipartitions for rooted trees
                calculateSTBipartitions(tree);
                
                geneTrees.add(tree);
                
                // Track internal node count for complexity analysis
                internalNodesCount += tree.nodes.size() - tree.leavesCount;
            } catch (RuntimeException e) {
                skippedTrees++;
                // Continue with next tree
            }
        }
        
        scanner.close();
        
        if (skippedTrees > 0) {
            System.err.println("Warning: Skipped " + skippedTrees + " invalid or empty gene trees.");
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

        // Output preprocessing statistics
        System.out.println( "taxon count : " + this.taxaMap.size());
        System.out.println("Gene trees count : " + geneTrees.size());
        System.out.println( "total internal nodes : " + internalNodesCount);
        System.out.println( "unique partitions : " + triPartitions.size());
        System.out.println( "unique STBipartitions : " + stBipartitions.size());

        // if(internalNodesCount == 50000){
        //     System.out.println("No polytomy, skipping");
        //     System.exit(-1);
        // }

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
        this.triPartitions = new HashMap<>();
        this.stBipartitions = new HashMap<>();
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
        this.triPartitions = new HashMap<>();
        this.stBipartitions = new HashMap<>();
        this.path = path;
        this.taxaMap = taxaMap;

        this.readGeneTrees(null);
    }
    
    private BitSet calculateSTBipartitionsUtil(TreeNode node, Tree tree) {
        if (node.isLeaf()) {
            BitSet leafSet = new BitSet(realTaxaCount);
            leafSet.set(node.taxon.id);
            return leafSet;
        }
        
        if (node.childs.size() == 2 && tree.isRooted()) {
            BitSet leftCluster = calculateSTBipartitionsUtil(node.childs.get(0), tree);
            BitSet rightCluster = calculateSTBipartitionsUtil(node.childs.get(1), tree);
            
            STBipartition stb = new STBipartition(leftCluster, rightCluster);
            stBipartitions.merge(stb, 1, Integer::sum);
            
            BitSet nodeCluster = (BitSet) leftCluster.clone();
            nodeCluster.or(rightCluster);
            return nodeCluster;
        } else {
            BitSet nodeCluster = new BitSet(realTaxaCount);
            for (TreeNode child : node.childs) {
                BitSet childCluster = calculateSTBipartitionsUtil(child, tree);
                nodeCluster.or(childCluster);
            }
            return nodeCluster;
        }
    }
    
    private void calculateSTBipartitions(Tree tree) {
        if (tree.isRooted()) {
            calculateSTBipartitionsUtil(tree.root, tree);
        }
    }

    /**
     * Generates candidate bipartitions for the inference algorithm.
     * This method creates a list of candidate bipartitions from the gene tree bipartitions,
     * ensuring that each candidate is a valid bipartition of the taxa set.
     * 
     * @return List of candidate bipartitions for inference
     */
    public List<STBipartition> generateCandidateBipartitions() {
        List<STBipartition> candidates = new ArrayList<>();
        Set<BitSet> processedClusters = new HashSet<>();
        
        // First, add all bipartitions from gene trees
        for (STBipartition stb : stBipartitions.keySet()) {
            candidates.add(stb);
            processedClusters.add(stb.cluster1);
            processedClusters.add(stb.cluster2);
        }
        
        // Generate complementary bipartitions for each cluster
        for (STBipartition stb : stBipartitions.keySet()) {
            BitSet cluster1 = stb.cluster1;
            BitSet cluster2 = stb.cluster2;
            
            // Add complementary bipartition if not already processed
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
        
        // Add all possible valid bipartitions for small clusters (size <= 3)
        for (int size = 1; size <= 3; size++) {
            for (int i = 0; i < realTaxaCount; i++) {
                BitSet cluster1 = new BitSet(realTaxaCount);
                cluster1.set(i);
                
                // Add all possible combinations of size 'size'
                for (int j = i + 1; j < realTaxaCount; j++) {
                    if (size > 1) {
                        cluster1.set(j);
                        if (size > 2) {
                            for (int k = j + 1; k < realTaxaCount; k++) {
                                cluster1.set(k);
                                addValidBipartition(cluster1, candidates, processedClusters);
                                cluster1.clear(k);
                            }
                        } else {
                            addValidBipartition(cluster1, candidates, processedClusters);
                        }
                        cluster1.clear(j);
                    } else {
                        addValidBipartition(cluster1, candidates, processedClusters);
                    }
                }
            }
        }
        
        return candidates;
    }
    
    private void addValidBipartition(BitSet cluster1, List<STBipartition> candidates, Set<BitSet> processedClusters) {
        if (!processedClusters.contains(cluster1)) {
            BitSet cluster2 = new BitSet(realTaxaCount);
            cluster2.set(0, realTaxaCount);
            cluster2.andNot(cluster1);
            
            if (cluster2.cardinality() > 0) {
                candidates.add(new STBipartition(cluster1, cluster2));
                processedClusters.add(cluster1);
                processedClusters.add(cluster2);
            }
        }
    }
} 


